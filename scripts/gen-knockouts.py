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
    parser.add_argument("--output", type=str, default="knockouts.csv", help="Name of the final output file (dumped in run directory).")
    parser.add_argument("--cleanup", action="store_true", help="Should intermediate files be removed?")

    # Parse commandline arguments
    args = parser.parse_args()
    inst_set_fpath = args.inst_set
    input_fpath = args.input
    num_tasks = args.num_tasks
    output_fname = args.output
    do_cleanup = args.cleanup

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
    ko_dat_files = []
    ko_dat_file_info = []
    for gi in range(len(data)):
        # Use gi + tree depth to uniquely identify each genotype
        genotype_data = data[gi]
        sequence = genotype_data["genome_sequence"]
        tree_depth = genotype_data["tree_depth"]
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
        # Save dat file location for later.
        ko_dat_files.append(os.path.join(run_dir, "data", output_fname_str))
        ko_dat_file_info.append({"gi": gi, "tree_depth": tree_depth, "sequence_len": len(sequence)})

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

    # Back fill fields, merge into one csv, cleanup .dat files
    # Fields: sequence, viability, gest_time, task.n, ko_pos, tree_depth
    merged_ko_lines = []
    merged_ko_header = None
    for ko_i in range(len(ko_dat_files)):
        ko_dat_file = ko_dat_files[ko_i]
        info = ko_dat_file_info[ko_i]
        tree_depth = info["tree_depth"]
        genotype_id = info["gi"]
        ko_data = utils.read_avida_dat_file(ko_dat_file)
        # Set ko_pos for each line in data.
        ko_data[0]["ko_pos"] = "-1"
        for i in range( info["sequence_len"] ):
            ko_data[i+1]["ko_pos"] = i
        for line in ko_data:
            line["tree_depth"] = tree_depth
            line["genotype_id"] = genotype_id

        merged_ko_fields = list(ko_data[0].keys())
        merged_ko_fields.sort()
        if merged_ko_header == None:
            merged_ko_header = merged_ko_fields
        elif merged_ko_header != merged_ko_fields:
            print("Header mismatch!")
            exit(-1)
        for line in ko_data:
            merged_ko_lines.append(",".join(str(line[field]) for field in merged_ko_header))

    # Write out consolidated knockouts file
    with open(os.path.join(run_dir, "data", output_fname), "w") as fp:
        fp.write(",".join(merged_ko_header) + "\n")
        fp.write("\n".join(merged_ko_lines))

    # Clean up intermediate files
    if do_cleanup:
        clean_up_cmd = ";".join(f"rm {file}" for file in ko_dat_files)
        subprocess.run(clean_up_cmd, shell=True)
        subprocess.run(f"rmdir {os.path.join(run_dir, 'data', 'knockouts')}", shell=True)
        subprocess.run(f"rm {os.path.join(run_dir, ko_analyze_fname)}", shell=True)
    print("Done.")

















if __name__ == "__main__":
    main()