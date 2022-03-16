"""
This script generates (and analyzes) all knockouts for each genotype in a target .dat file
"""

import os, argparse, subprocess
import utilities as utils

base_knockout_analysis="""
PURGE_BATCH

<<LOAD_KNOCKOUTS>>

RECALC

DETAIL <<OUTPUT_FNAME>> <<DETAIL_ARGS>>
"""

ko_detail_args = ["sequence", "viable", "gest_time"] # Ask for this info + tasks


def main():
    parser = argparse.ArgumentParser(
        description="Generate + analyze all knockouts genotypes for a given "
    )
    parser.add_argument("--inst_set", type=str, help="Path to the instruction set .cfg file. Must include nop-x.")
    parser.add_argument("--input", type=str, help="Path to the input .dat file with genotypes to perform knockouts on.")
    parser.add_argument("--num_tasks", type=int, default=6, help="How many (avida) tasks to output performance?")
    parser.add_argument("--avida_args", type=str, default="", help="Only used if 'run_avida_analysis' is true. Any command line arguments that should be used when running Avida's analyze mode.")
    parser.add_argument("--run_dir", type=str, default="./", help="Only used if 'run_avida_analysis' is true. Where should we run avida analyze mode (this is where the detail file will be output)? ")

    # Parse commandline arguments
    args = parser.parse_args()
    inst_set_fpath = args.inst_set
    input_fpath = args.input
    num_tasks = args.num_tasks

    avida_args = args.avida_args
    run_dir = args.run_dir

    # Verify that the given instruction set file exists
    if not os.path.exists(inst_set_fpath):
        print(f"Unable to find instruction set file: {inst_set_fpath}.")
        exit(-1)

    # Verify input file exists.
    if not os.path.exists(input_fpath):
        print(f"Unable to find input file: {input_fpath}")
        exit(-1)

    # Verify that run directory exists
    if not os.path.exists(run_dir):
        print(f"Unable to find run directory: {run_dir}")
        exit(-1)

    ################################################################
    # Parse instruction set file
    inst_map = utils.read_avida_inst_set_file(inst_set_fpath)
    # Verify that 'nop-x' is included.
    if not "nop-X" in inst_map:
        print(f"'nop-X' instruction not included in provided instruction set file. Found instructions: {inst_map}")
        exit(-1)
    ################################################################

    # Read the input file
    data = utils.read_avida_dat_file(input_fpath)
    # Verify the contents of the loaded data
    if not len(data):
        print("Failed to find any genomes in the given input.")
        exit(-1)

    knockout_analyze_content = ""

    # For each genotype,
    # - Generate all knockouts
    # - Write knockouts to knockout analyze content
    # - Detail knockouts, include genotype ID in file name
    # Run avida -a
    # Back-label genotype ID of source genotype + position of the knockout

    # Loop over genotypes to generate knockouts and produce an avida analysis script.
    for gi in range(len(data)):
        # Use gi + tree depth to uniquely identify each genotype
        genotype_data = data[gi]
        sequence = genotype_data["genome_sequence"]
        # tree_depth = genotype_data["tree_depth"]
        knockouts = [list(sequence) for _ in range(len(sequence))]
        for pos in range(len(sequence)):
            knockouts[pos][pos] = inst_map["nop-X"]
        knockouts = ["".join(ko_seq) for ko_seq in knockouts]
        # First, add original; then add all knockouts.
        load_knockouts_str = f"LOAD_SEQUENCE {sequence}\n"
        load_knockouts_str += "\n".join([f"LOAD_SEQUENCE {seq}" for seq in knockouts])
        detail_args_str = " ".join(ko_detail_args + [f"task.{i}" for i in range(num_tasks)])
        output_fname_str = f"knockouts/knockouts_id-{gi}.dat"

        genotype_content = base_knockout_analysis
        genotype_content = genotype_content.replace("<<LOAD_KNOCKOUTS>>", load_knockouts_str)
        genotype_content = genotype_content.replace("<<OUTPUT_FNAME>>", output_fname_str)
        genotype_content = genotype_content.replace("<<DETAIL_ARGS>>", detail_args_str)

        knockout_analyze_content += genotype_content + "\n"

    # Write knockouts avida analyze file.
    ko_analyze_fname = "analyze_knockouts.cfg"
    with open(os.path.join(run_dir, ko_analyze_fname), "w") as fp:
        fp.write(knockout_analyze_content)

    # Run Avida analyze mode
    inst_set_abs_path = os.path.abspath(inst_set_fpath)
    avida_cmd = f"./avida {avida_args} -def INST_SET {inst_set_abs_path} -set ANALYZE_FILE {ko_analyze_fname} -a > knockouts.log"
    print(f"Running knockouts through Avida's analyze mode:")
    print(f" - Run directory: {run_dir}")
    print(f" - Avida run command: {avida_cmd}")
    analyze_mode_cmd = f"cd {run_dir}; {avida_cmd}"
    subprocess.run(analyze_mode_cmd, shell=True)

    # Back fill fields
    # TODO









if __name__ == "__main__":
    main()