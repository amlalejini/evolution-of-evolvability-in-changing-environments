'''
Aggregate data
'''

import argparse, os, sys, errno, subprocess, csv

run_identifier = "RUN_"

traits_primary = ["not","nand","and","ornot","or","andnot"]
traits_env_a = {"not", "and", "or"}
traits_env_b = {"nand", "ornot", "andnot"}

profile_env_a = "101010"
profile_env_b = "010101"
profile_all = "111111"

# because we want smaller file sizes, only keep fields that we want to look at
# time_data_time_series_fields = ["average_generation"]

mutation_rates_time_series_fields = [
    "average_copy_mutation_rate",
    "variance_in_copy_mutation_rate",
    "standard_deviation_in_copy_mutation_rate"
]

"""
This is functionally equivalent to the mkdir -p [fname] bash command
"""
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def extract_params_cmd_log(path):
    content = None
    with open(path, "r") as fp:
        content = fp.read().strip()
    content = content.replace("./avida", "")
    params = [param.strip() for param in content.split("-set") if param.strip() != ""]
    cfg = {param.split(" ")[0]:param.split(" ")[1] for param in params}
    return cfg

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

def read_csv(file_path):
    content = None
    with open(file_path, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    content = content[1:]
    lines = [{header[i]: l[i] for i in range(len(header))} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
    return lines

def simple_match_coeff(a, b):
    if len(a) != len(b):
        print(f"Length mismatch! {a} {b}")
        exit(-1)
    return sum(ai==bi for ai,bi in zip(a,b))

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--update", type=int, help="Update to pull data for?")
    parser.add_argument("--time_series_range", type=int, help="The range (in updates) to collect time series data?", nargs=2)

    # Parse command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump
    update = args.update
    time_series_range = args.time_series_range

    # Verify that the given data directory exits
    if not os.path.exists(data_dir):
        print("Unable to find data directory.")
        exit(-1)

    # Create the directory to dump aggregated data (if it doesn't already exist)
    mkdir_p(dump_dir)

    # Aggregate run directories.
    run_dirs = [run_dir for run_dir in os.listdir(data_dir) if run_identifier in run_dir]
    print(f"Found {len(run_dirs)} run directories.")

    # Create file to hold time series data
    time_series_content = []    # This will hold all the lines to write out for a single run; written out for each run.
    time_series_header = None   # Holds the time series file header (verified for consistency across runs)
    time_series_fpath = os.path.join(dump_dir, f"time_series_u{time_series_range[0]}-u{time_series_range[1]}.csv")
    time_series_update_set = None
    with open(time_series_fpath, "w") as fp:
        fp.write("")

    # Create data structure to hold summary information for each run (1 element per run)
    summary_header = None
    summary_content_lines = []

    # Only keep lines that fall within specified time series range.
    def keep_line(u): return u <= time_series_range[1] and u >= time_series_range[0]

    # Loop over runs, aggregating data from each.
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)

        # Skip over (but make note of) incomplete runs.
        if not os.path.exists(os.path.join(run_path, 'data', 'analysis')):
            print('Skipping: ', run_path)
            continue

        summary_info = {}                   # Hold summary information about run. Indexed by field.
        time_series_info = {}               # Hold time series information. Indexed by update.

        print(f"Processing: {run_path}")

        ############################################################
        # Extract commandline configuration settings (from cmd.log file)
        cmd_log_path = os.path.join(run_path, "cmd.log")
        cmd_params = extract_params_cmd_log(cmd_log_path)

        # Infer environment parameterization from event file name.
        env_cond = cmd_params["EVENT_FILE"].strip(".cfg").replace("events_", "").lower()

        summary_info["env_cond"] = env_cond
        summary_info["update"] = update

        for field in cmd_params:
            summary_info[field] = cmd_params[field]
        ############################################################

        ############################################################
        # Extract .spop file info. Create a lookup table that contains information about each genotype in the spop file.
        # This lets us cross reference information from analyze mode with information from .spop files.
        # - What are the unique keys shared across analyze mode output/.spop file?
        #   - tuple(sequence, update born)

        # Read spop file.
        spop_data = read_avida_dat_file(os.path.join(run_path, "data", f"detail-{update}.spop"), True)

        # Each entry in the lookup table is indexed off of update born and genome sequence.
        spop_lookup_table = {
            tuple([line["update_born"], line["genome_sequence"]]):line
            for line in spop_data
        }

        # Convenient function for accessing the spop_lookup_table.
        def spop_lookup(update_born, sequence, field):
            return spop_lookup_table[tuple([update_born, sequence])][field]

        # Clear raw data now that we're done with it.
        spop_data = None
        ############################################################

        ############################################################
        # Extract environment-specific final dominant information.
        dominant_data = read_avida_dat_file(os.path.join(run_path, "data", "analysis", "final_dominant.dat"))

        # Final dominant data files should only have one genotype in them.
        if len(dominant_data) != 1:
            print("Unexpected number of genotypes in final_dominant data file(s).")
            exit(-1)
        dominant_data = dominant_data[0]

        # Collect dominant genotype data.
        summary_info["dom_genome_length"] = dominant_data["genome_length"]
        summary_info["dom_generation_born"] = spop_lookup(
            update_born=dominant_data["update_born"],
            sequence=dominant_data["genome_sequence"],
            field="generation_born"
        )

        # Calculate task phenotype as a bitstring
        phenotype = "".join([dominant_data[trait] for trait in traits_primary])
        phenotype_task_order = ";".join(traits_primary)

        # Calculate match score for both environment a and b
        match_score_env_a = simple_match_coeff(phenotype, profile_env_a)
        match_score_env_b = simple_match_coeff(phenotype, profile_env_b)

        # Save summary information for final_dominant.dat
        summary_info["dom_phenotype"] = phenotype
        summary_info["dom_phenotype_task_order"] = phenotype_task_order
        summary_info["dom_match_score_env_a"] = match_score_env_a
        summary_info["dom_match_score_env_b"] = match_score_env_b
        ############################################################

        ############################################################
        # Extract data from dominant lineage
        lineage_data = read_avida_dat_file(os.path.join(run_path, "data", "analysis", "lineage.dat"))

        # Grab dominant lineage length (genotypes)
        summary_info["dom_lineage_length_genotypes"] = len(lineage_data)

        # Create data structure to hold phenotype information over time.
        task_profiles_ot = [{
            "profile": None,
            "match-score_env-a": None,
            "match-score_env-b": None,
            "match-score_all": None
        } for _ in range(len(lineage_data))]

        # Collect task profile/mutation information along lineage
        total_mutations = 0
        for i in range(len(lineage_data)):
            # Parent dist is a measure of mutational distance because genomes are fixed length.
            parent_dist = int(lineage_data[i]["parent_distance"])
            total_mutations += parent_dist

            # Collect ancestor phenotype
            ancestor_phenotype = "".join([lineage_data[i][trait] for trait in traits_primary])

            ancestor_match_score_env_a = simple_match_coeff(ancestor_phenotype, profile_env_a)
            ancestor_match_score_env_b = simple_match_coeff(ancestor_phenotype, profile_env_b)
            ancestor_match_score_all = simple_match_coeff(ancestor_phenotype, profile_all)

            task_profiles_ot[i]["profile"] = ancestor_phenotype
            task_profiles_ot[i]["muts_from_parent"] = parent_dist
            task_profiles_ot[i]["match-score_env-a"] = ancestor_match_score_env_a
            task_profiles_ot[i]["match-score_env-b"] = ancestor_match_score_env_b
            task_profiles_ot[i]["match-score_all"] = ancestor_match_score_all

        # Save summary info about mutation accumulation
        summary_info["dom_lineage_total_mut_cnt"] = total_mutations

        # Calculate task profile volatility
        task_profile_volatility = 0
        for i in range(len(task_profiles_ot)):
            if i:
                current_profile = task_profiles_ot[i]["profile"]
                previous_traits = task_profiles_ot[i-1]["profile"]
                task_profile_volatility += int(current_profile != previous_traits)
        summary_info["dom_lineage_task_volatility"] = task_profile_volatility

        lineage_data = None
        ############################################################

        ############################################################
        # Extract information from mutation_rates.dat
        mutation_rates_data = read_avida_dat_file(os.path.join(run_path, "data", "mutation_rates.dat"))

        # Summary information from mutation rates data
        final_mutation_rate_data = [line for line in mutation_rates_data if int(line["update"]) == update][0]

        summary_info["final_average_copy_mutation_rate"] = final_mutation_rate_data["average_copy_mutation_rate"]
        summary_info["final_variance_in_copy_mutation_rate"] = final_mutation_rate_data["variance_in_copy_mutation_rate"]
        summary_info["final_standard_deviation_in_copy_mutation_rate"] = final_mutation_rate_data["standard_deviation_in_copy_mutation_rate"]

        # Collect time series data from mutation rates
        mutation_rates_data_ts = {line["update"]: {field: line[field] for field in mutation_rates_time_series_fields} for line in mutation_rates_data if keep_line(int(line["update"]))}
        mutation_rates_updates = set(mutation_rates_data_ts.keys())

        # If first time series data, configure time series information. Otherwise, check for consistency.
        if time_series_update_set == None:
            time_series_update_set = mutation_rates_updates
        if len(time_series_info) == 0:
            for u in time_series_update_set: time_series_info[u] = {}
        if mutation_rates_updates != time_series_update_set:
            print("Resolution mismatch for mutation_rates.dat!")

        # Add mutation rates data to time series info
        for u in time_series_update_set:
            for field in mutation_rates_data_ts[u]:
                time_series_info[u]["mutation_rates_" + field] = mutation_rates_data_ts[u][field]

        # Clear out mutation rates data
        mutation_rates_data = None
        mutation_rates_data_ts = None
        ############################################################

        ############################################################
        # Extract time.dat data
        time_data = read_avida_dat_file(os.path.join(run_path, "data", "time.dat"))

        # Summery information
        final_time_data = [line for line in time_data if int(line["update"]) == update][0]
        summary_info["time_average_generation"] = final_time_data["average_generation"]

        time_data = None # release time_data
        ############################################################

        ############################################################
        # Output time series data for this run
        # Add extra fields
        for u in time_series_info:
            time_series_info[u]["update"] = u # Make sure that update is a field on every line
            time_series_info[u]["RANDOM_SEED"] = summary_info["RANDOM_SEED"]
            time_series_info[u]["MAX_HEAD_COPY_COST"] = summary_info["MAX_HEAD_COPY_COST"]
            time_series_info[u]["COPY_MUT_PROB"] = summary_info["COPY_MUT_PROB"]
            time_series_info[u]["env_cond"] = summary_info["env_cond"]

        # Compute time series header from time_series_info
        time_series_fields = list(time_series_info[str(time_series_range[0])].keys())
        time_series_fields.sort()
        # If we haven't written the header, write it.
        write_header = False
        if time_series_header == None:
            write_header = True
            time_series_header = ",".join(time_series_fields)
        elif time_series_header != ",".join(time_series_fields):
            print("Time series header mismatch!")
            exit(-1)

        # Write time series content line-by-line
        time_series_content = []
        update_order = list(map(int, time_series_info.keys()))
        update_order.sort()
        for u in update_order:
            time_series_content.append(",".join([str(time_series_info[str(u)][field]) for field in time_series_fields]))
        with open(time_series_fpath, "a") as fp:
            if write_header: fp.write(time_series_header)
            fp.write("\n")
            fp.write("\n".join(time_series_content))
        time_series_content = []
        ############################################################

        ############################################################
        # Add summary_info to aggregate content
        summary_fields = list(summary_info.keys())
        summary_fields.sort()
        if summary_header == None:
            summary_header = summary_fields
        elif summary_header != summary_fields:
            print("Header mismatch!")
            exit(-1)
        summary_line = [str(summary_info[field]) for field in summary_fields]
        summary_content_lines.append(",".join(summary_line))
        ############################################################

    # write out aggregate data
    with open(os.path.join(dump_dir, "aggregate.csv"), "w") as fp:
        out_content = ",".join(summary_header) + "\n" + "\n".join(summary_content_lines)
        fp.write(out_content)

if __name__ == "__main__":
    main()
