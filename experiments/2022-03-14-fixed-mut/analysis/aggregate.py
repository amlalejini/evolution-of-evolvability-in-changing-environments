'''
Aggregate data
'''

import argparse, os, sys, errno, csv, pathlib
sys.path.append(os.path.join(pathlib.Path(os.path.dirname(os.path.abspath(__file__))).parents[2], "scripts"))
import utilities as utils

run_identifier = "RUN_"

tasks_primary = ["not","nand","and","ornot","or","andnot"]
tasks_env_a = {"not", "and", "or"}
tasks_env_b = {"nand", "ornot", "andnot"}

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

average_data_time_series_fields = [
    "fitness",
    "true_replication_rate_(based_on_births/update,_time-averaged)"
]

pop_tasks_time_series_fields = [task for task in tasks_primary]

phylo_time_series_fields = [
    "mean_evolutionary_distinctiveness",
    "max_evolutionary_distinctiveness",
    "variance_evolutionary_distinctiveness",
    "mean_pairwise_distance",
    "max_pairwise_distance",
    "variance_pairwise_distance",
    "current_phylogenetic_diversity",
    "sum_pairwise_distance",
    "mrca_depth",
    "diversity",
    "mrca_changes"
]

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
    utils.mkdir_p(dump_dir)

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

    incomplete_runs = []

    # Only keep lines that fall within specified time series range.
    def keep_line(u): return u <= time_series_range[1] and u >= time_series_range[0]

    # Loop over runs, aggregating data from each.
    total_runs = len(run_dirs)
    cur_run_i = 0
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)

        summary_info = {}                   # Hold summary information about run. Indexed by field.
        time_series_info = {}               # Hold time series information. Indexed by update.

        cur_run_i += 1
        print(f"Processing ({cur_run_i}/{total_runs}): {run_path}")

        # Skip over (but make note of) incomplete runs.
        required_files = [
            os.path.join("data", "phylodiversity.csv"),
            os.path.join("data", "analysis", "final_dominant.dat"),
            os.path.join("data", "analysis", "lineage.dat"),
            os.path.join("data", "mutation_rates.dat"),
            os.path.join("data", "average.dat"),
            os.path.join("data", "tasks.dat"),
            os.path.join("data", f"detail-{update}.spop"),
            os.path.join("cmd.log")
        ]
        incomplete = any([not os.path.exists(os.path.join(run_path, req)) for req in required_files])
        if incomplete:
            print("  - Failed to find all required files!")
            incomplete_runs.append(run_dir)
            continue

        ############################################################
        # Extract commandline configuration settings (from cmd.log file)
        cmd_log_path = os.path.join(run_path, "cmd.log")
        cmd_params = utils.extract_params_cmd_log(cmd_log_path)

        # Infer environment parameterization from event file name.
        env_id = cmd_params["EVENT_FILE"].strip(".cfg").replace("events_", "").lower()
        env_cfg = {pair.split("-")[0]:pair.split("-")[1] for pair in env_id.split("_") if not "run-" in pair}
        env_cfg_fields = sorted(list(env_cfg.keys()))
        env_condition = "_".join([f"{field}-{env_cfg[field]}" for field in env_cfg_fields])
        env_type = env_cfg["env"]
        env_chg_rate = env_cfg["rate"] if "rate" in env_cfg else "none"

        summary_info["env_condition"] = env_condition
        summary_info["env_type"] = env_type
        summary_info["env_chg_rate"] = env_chg_rate

        for field in cmd_params:
            summary_info[field] = cmd_params[field]
        ############################################################

        ############################################################
        # Extract .spop file info. Create a lookup table that contains information about each genotype in the spop file.
        # This lets us cross reference information from analyze mode with information from .spop files.
        # - What are the unique keys shared across analyze mode output/.spop file?
        #   - tuple(sequence, update born)

        # Read spop file.
        spop_data = utils.read_avida_dat_file(os.path.join(run_path, "data", f"detail-{update}.spop"), True)

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
        dominant_data = utils.read_avida_dat_file(os.path.join(run_path, "data", "analysis", "final_dominant.dat"))

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
        phenotype = "".join([dominant_data[trait] for trait in tasks_primary])
        phenotype_task_order = ";".join(tasks_primary)

        # Calculate match score for both environment a and b
        match_score_env_a = utils.simple_match_coeff(phenotype, profile_env_a)
        match_score_env_b = utils.simple_match_coeff(phenotype, profile_env_b)

        # Save summary information for final_dominant.dat
        summary_info["dom_phenotype"] = phenotype
        summary_info["dom_phenotype_task_order"] = phenotype_task_order
        summary_info["dom_match_score_env_a"] = match_score_env_a
        summary_info["dom_match_score_env_b"] = match_score_env_b
        ############################################################

        ############################################################
        # Extract data from dominant lineage
        lineage_data = utils.read_avida_dat_file(os.path.join(run_path, "data", "analysis", "lineage.dat"))

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
            ancestor_phenotype = "".join([lineage_data[i][trait] for trait in tasks_primary])

            ancestor_match_score_env_a = utils.simple_match_coeff(ancestor_phenotype, profile_env_a)
            ancestor_match_score_env_b = utils.simple_match_coeff(ancestor_phenotype, profile_env_b)
            ancestor_match_score_all = utils.simple_match_coeff(ancestor_phenotype, profile_all)

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
        mutation_rates_data = utils.read_avida_dat_file(os.path.join(run_path, "data", "mutation_rates.dat"))

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
        # time_data = utils.read_avida_dat_file(os.path.join(run_path, "data", "time.dat"))

        # # Summery information
        # final_time_data = [line for line in time_data if int(line["update"]) == update][0]
        # summary_info["time_average_generation"] = final_time_data["average_generation"]

        # time_data = None # release time_data
        ############################################################

        ############################################################
        # Extract average.dat data
        average_data = utils.read_avida_dat_file(os.path.join(run_path, "data", "average.dat"))
        final_average_data = [line for line in average_data if int(line["update"]) == update][0]

        average_data_ts = {line["update"]: {field: line[field] for field in average_data_time_series_fields} for line in average_data if keep_line(int(line["update"]))}
        # Summary info
        summary_info["avg_fitness"] = final_average_data["fitness"]
        summary_info["avg_true_replication_rate"] = final_average_data["true_replication_rate_(based_on_births/update,_time-averaged)"]
        # Average fitness over time
        for u in time_series_update_set:
            time_series_info[u]["avg_fitness"] = average_data_ts[u]["fitness"]
            time_series_info[u]["avg_true_replication_rate"] = average_data_ts[u]["true_replication_rate_(based_on_births/update,_time-averaged)"]

        average_data = None
        average_data_ts = None
        ############################################################

        ############################################################
        # Extract tasks.dat data
        tasks_data = utils.read_avida_dat_file(os.path.join(run_path, "data", "tasks.dat"))
        final_tasks_data = [line for line in tasks_data if int(line["update"]) == update][0]
        tasks_data_ts = {line["update"]: {field: line[field] for field in pop_tasks_time_series_fields} for line in tasks_data if keep_line(int(line["update"]))}

        # Summary info
        for field in pop_tasks_time_series_fields:
            summary_info[f"pop_tasks_{field}"] = final_tasks_data[field]

        # Time series info
        for u in time_series_update_set:
            for field in pop_tasks_time_series_fields:
                time_series_info[u][f"pop_tasks_{field}"] = tasks_data_ts[u][field]

        tasks_data = None
        tasks_data_ts = None
        ############################################################

        ############################################################
        # Extract phylodiversity data
        phylodiversity_path = os.path.join(run_path, "data", "phylodiversity.csv")
        phylodiversity_data = utils.read_csv(phylodiversity_path)
        final_phylo_data = [line for line in phylodiversity_data if int(line["update"]) == update][0]
        phylo_data_ts = {line["update"]: {field: line[field] for field in phylo_time_series_fields} for line in phylodiversity_data if keep_line(int(line["update"]))}

        # Collect run summary information
        for field in final_phylo_data:
            if field == "update": continue
            summary_info[f"phylo_{field}"] = final_phylo_data[field]

        # Time series information
        for u in time_series_update_set:
            for field in phylo_time_series_fields:
                time_series_info[u][f"phylo_{field}"] = phylo_data_ts[u][field]

        # Clear our phylo data variables
        phylodiversity_data = None
        final_phylo_data = None
        phylo_data_ts = None
        ############################################################

        ############################################################
        # Output time series data for this run
        # Add extra fields
        for u in time_series_info:
            time_series_info[u]["update"] = u # Make sure that update is a field on every line
            time_series_info[u]["RANDOM_SEED"] = summary_info["RANDOM_SEED"]
            time_series_info[u]["MAX_HEAD_COPY_COST"] = summary_info["MAX_HEAD_COPY_COST"]
            time_series_info[u]["COPY_MUT_PROB"] = summary_info["COPY_MUT_PROB"]
            time_series_info[u]["env_condition"] = summary_info["env_condition"]
            time_series_info[u]["env_type"] = summary_info["env_type"]
            time_series_info[u]["env_chg_rate"] = summary_info["env_chg_rate"]

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

    # Write out incomplete runs, sort them!
    incomplete_runs.sort()
    with open(os.path.join(dump_dir, "incomplete_runs_agg.log"), "w") as fp:
        fp.write("\n".join(incomplete_runs))

if __name__ == "__main__":
    main()
