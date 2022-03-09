import argparse, os, errno, subprocess
import GenomeManipulator as gm

base_mutant_analyze_script = """
PURGE_BATCH

<<LOAD_MUTANTS>>

RECALC

DETAIL <<OUTPUT_FNAME>> <<DETAIL_ARGS>>
"""
mutant_detail_args = ["sequence", "viable"] # Ask for this info + tasks

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def read_avida_dat_file(path, backfill_missing_fields=False):
    content = None
    with open(path, "r") as fp:
        content = fp.read().strip().split("\n")
    legend_start = 0
    legend_end = 0
    # Where does the legend table start?
    for line_i in range(0, len(content)):
        line = content[line_i].strip()
        if line == "# Legend:":         # Handles analyze mode detail files.
            legend_start = line_i + 1
            break
        if "#  1:" in line:             # Handles time.dat, mutation_rate.dat, etc.
            legend_start = line_i
            break
    # For each line in legend table, extract field
    fields = []
    for line_i in range(legend_start, len(content)):
        line = content[line_i].strip()
        if line == "":
            legend_end = line_i
            break
        # patch 3-input logic tasks because avida file format is nonsense
        if "Logic 3" in line:
            line = line.split("(")[0]

        fields.append( line.split(":")[-1].strip().lower().replace(" ", "_") )
    data = []
    for line_i in range(legend_end, len(content)):
        line = content[line_i].strip()
        if line == "": continue
        data_line = line.split(" ")
        if len(data_line) > len(fields):
            print("found more items than there are fields!")
            print(fields)
            print(data_line)
            exit(-1)
        elif backfill_missing_fields:
            num_backfill = len(fields) - len(data_line)
            for _ in range(num_backfill): data_line.append("")
        elif len(data_line) != len(fields):
            print("data fields mismatch!")
            print(fields)
            print(data_line)
            exit(-1)
        data.append({field:value for field,value in zip(fields, data_line)})
    return data

def main():
    parser = argparse.ArgumentParser(
        description="Generate + analyze all mutants for given sequence"
    )
    parser.add_argument("--inst_set", type=str, help="Path to the instruction set .cfg file.")
    parser.add_argument("--input", type=str, help="Path to the input .dat file with the ancestral genotype to landscape.")
    parser.add_argument("--steps", type=int, default=1, choices=[1,2], help="How many mutational steps should we take?")
    parser.add_argument("--num_tasks", type=int, default=6, help="How many (avida) tasks to output performance?")
    parser.add_argument("--dump", type=str, default="./dump/", help="Where to dump analyze file generated for this?")
    parser.add_argument("--analysis_output", type=str, default="mutants.dat", help="File name for the detail file from analyzing mutants.")

    parser.add_argument("--run_avida_analysis", action="store_true", help="Should this script run the mutants through analyze mode? If not, just output the analyze cfg file.")
    parser.add_argument("--avida_args", type=str, default="", help="Only used if 'run_avida_analysis' is true. Any command line arguments that should be used when running Avida's analyze mode.")
    parser.add_argument("--run_dir", type=str, default="./", help="Only used if 'run_avida_analysis' is true. Where should we run avida analyze mode (this is where the detail file will be output)? ")

    # Parse command line arguments
    args = parser.parse_args()
    mutation_steps = args.steps
    inst_set_fpath = args.inst_set
    input_fpath = args.input
    dump_dir = args.dump
    num_tasks = args.num_tasks

    run_avida_analysis = args.run_avida_analysis
    avida_args = args.avida_args
    run_dir = args.run_dir
    analysis_output_fname = args.analysis_output

    # Verify that the given instruction set file exists
    if not os.path.exists(inst_set_fpath):
        print(f"Unable to find instruction set file: {inst_set_fpath}.")
        exit(-1)
    # Verify input file exists.
    if not os.path.exists(input_fpath):
        print(f"Unable to find input file: {input_fpath}")
        exit(-1)

    # Create the dump directory if it doesn't already exist.
    mkdir_p(dump_dir)
    mutant_detail_fpath = os.path.join(dump_dir, analysis_output_fname)

    # Create a genome manipulator to do the landscaping
    manipulator = gm.GenomeManipulator(inst_set_fpath)

    # Read input file
    data = read_avida_dat_file(input_fpath)
    # Verify contents of loaded data
    if len(data) == 0:
        print("Failed to find any genomes in given input.")
        exit(-1)
    elif len(data) > 1:
        print("Found more than one genome, landscaping only the first genome in the file.")
    # Only want one focal genome.
    data = data[0]
    sequence = data["genome_sequence"]

    mutants = []
    # Compute one-step mutants
    one_step_mutants = manipulator.generate_all_point_mutants(sequence)
    if mutation_steps == 2:
        # Compute all two-step mutants
        for mutant in one_step_mutants:
            mutants.append(mutant)
            second_step_mutants = manipulator.generate_all_point_mutants(mutant)
            mutants += second_step_mutants
    else:
        # Otherwise, one-step is all we care about.
        mutants = one_step_mutants
    mutants.append(sequence)
    print(f"Number of mutants generated: {len(mutants)}")

    # Build analyze file to run mutants through avida test cpus
    analyze_file_fpath = os.path.join(dump_dir, "analyze_mutants.cfg")
    mutant_load_str = "\n".join([f"LOAD_SEQUENCE {genome}" for genome in mutants]) + "\n"
    detail_args = mutant_detail_args + [f"task.{i}" for i in range(num_tasks)]

    mutant_analyze_content = base_mutant_analyze_script.replace("<<OUTPUT_FNAME>>", analysis_output_fname)
    mutant_analyze_content = mutant_analyze_content.replace("<<DETAIL_ARGS>>", " ".join(detail_args))
    mutant_analyze_content = mutant_analyze_content.replace("<<LOAD_MUTANTS>>", mutant_load_str)

    with open(analyze_file_fpath, "w") as fp:
        fp.write(mutant_analyze_content)

    # Run mutants through analyze mode?
    if not run_avida_analysis:
        print("Done.")
        return

    if not os.path.exists(run_dir):
        print(f"Unable to find run_dir: {run_dir}")

    avida_cmd = f"./avida {avida_args} -set VERBOSITY 0 -set ANALYZE_FILE {analyze_file_fpath} -a"
    print(f"Running mutants through Avida's analyze mode:")
    print(f" - Run directory: {run_dir}")
    print(f" - Avida run command: {avida_cmd}")
    analyze_mode_cmd = f"cd {run_dir}; {avida_cmd}"
    subprocess.run(analyze_mode_cmd, shell=True)

    print("Done.")

if __name__ == "__main__":
    main()