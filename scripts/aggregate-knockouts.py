'''
Process + aggregate knockout analysis information.
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

lineage_arch_long_cfg_fields = [
    "env_condition",
    "env_type",
    "env_chg_rate",
    "RANDOM_SEED",
    "COPY_MUT_PROB"
]

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--no_arch", action="store_true", help="Don't aggregate architectures.")

    # Parse command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump
    agg_arch = not args.no_arch

    # Verify that the given data directory exits
    if not os.path.exists(data_dir):
        print("Unable to find data directory.")
        exit(-1)

    # Create the directory to dump aggregated data (if it doesn't already exist)
    utils.mkdir_p(dump_dir)

    # Aggregate run directories.
    run_dirs = [run_dir for run_dir in os.listdir(data_dir) if run_identifier in run_dir]
    print(f"Found {len(run_dirs)} run directories.")

    run_summary_header = None
    run_summary_lines = []

    if agg_arch:
        lineage_arch_long_header = None
        # lineage_arch_long_lines = []
        lineage_arch_long_fpath = os.path.join(dump_dir, "lineage_arch_long.csv")
        with open(lineage_arch_long_fpath, "w") as fp:
            fp.write("")

    # Loop over runs, aggregating data from each.
    total_runs = len(run_dirs)
    cur_run_i = 0
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)
        # Skip over (but make note of) incomplete runs.
        if not os.path.exists(os.path.join(run_path, 'data', 'knockouts.csv')):
            print("Skipping: ", run_path)
            print(" - Failed to find 'knockouts.csv'")
            continue

        run_summary_info = {}
        lineage_arch_long_info = []

        cur_run_i += 1
        print(f"Processing ({cur_run_i}/{total_runs}): {run_path}")

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

        run_summary_info["env_condition"] = env_condition
        run_summary_info["env_type"] = env_type
        run_summary_info["env_chg_rate"] = env_chg_rate

        for field in cmd_params:
            run_summary_info[field] = cmd_params[field]
        ############################################################

        ############################################################
        # Read knockouts.csv
        knockouts_data = utils.read_csv(os.path.join(run_path, "data", "knockouts.csv"))

        # Collect genotype ids (i.e., how many knockout analyses are represented in this file?)
        genotype_ids = list({line["genotype_id"] for line in knockouts_data})
        # Collect the original genotypes (organized by genotype id)
        orig_genotypes = {line["genotype_id"]:line for line in knockouts_data if line["ko_pos"] == "-1"}
        # Collect all the knockouts (organized by genotype id)
        knockouts_by_genotype_id = {gid:[line for line in knockouts_data if line["genotype_id"] == gid and line["ko_pos"] != "-1"] for gid in genotype_ids}
        genotype_summary_info = {gid:{} for gid in genotype_ids}
        # For each genotype, process knockout analysis
        for genotype_id in genotype_ids:
            knockouts = knockouts_by_genotype_id[genotype_id]
            orig_genotype = orig_genotypes[genotype_id]
            orig_tasks_performed = {task for task in tasks_primary if int(orig_genotype[task]) > 0}
            orig_genotype["tasks_performed"] = orig_tasks_performed
            position_info = {}
            viability_sites = 0     # Calculate number of sites important for viability
            task_sites = 0          # Calculate number of sites important for tasks (any)
            multi_task_sites = 0    # Calculate number of multi-task sites (overlap)
            task_env_a_sites = 0    # Calculate number of sites important for env-a tasks
            task_env_b_sites = 0    # Calculate number of sites important for env-b tasks
            for ko in knockouts:
                pos = int(ko["ko_pos"]) # What position are we getting information about w/this knockout?
                info = {}               # Will store info about this position.
                # What tasks does this knockout perform?
                ko_tasks_performed = {task for task in tasks_primary if int(ko[task]) > 0}
                ko["tasks_performed"] = ko_tasks_performed
                # Is this position important for viability?
                info["viability"] = not bool(int(ko["is_viable_(0/1)"]))
                # If ko is viable, what tasks (if any) is this position important for?
                info["tasks"] = {}
                info["encodes_env_a_task"] = False
                info["encodes_env_b_task"] = False
                if ko["is_viable_(0/1)"] == "1":
                    info["tasks"] = orig_tasks_performed - ko_tasks_performed
                    info["encodes_env_a_task"] = len(tasks_env_a.intersection(info["tasks"])) > 0
                    info["encodes_env_b_task"] = len(tasks_env_b.intersection(info["tasks"])) > 0
                # Add stats about site
                viability_sites += int(info["viability"])
                task_sites += int(len(info["tasks"]) > 0)
                multi_task_sites += int(len(info["tasks"]) > 1)
                task_env_a_sites += int(info["encodes_env_a_task"])
                task_env_b_sites += int(info["encodes_env_b_task"])
                # Set info about site
                position_info[pos] = info

            # Save summary information about this genotype
            genotype_summary_info[genotype_id]["viability_sites"] = viability_sites
            genotype_summary_info[genotype_id]["task_sites"] = task_sites
            genotype_summary_info[genotype_id]["multi_task_sites"] = multi_task_sites
            genotype_summary_info[genotype_id]["task_env_a_sites"] = task_env_a_sites
            genotype_summary_info[genotype_id]["task_env_b_sites"] = task_env_b_sites
            genotype_summary_info[genotype_id]["position_info"] = position_info

        # print(genotype_summary_info)

        # Identify extant genotype
        extant_genotype_id = None
        for gid in orig_genotypes:
            if extant_genotype_id == None:
                extant_genotype_id = gid
                continue
            cur_extant_update_born = int(orig_genotypes[extant_genotype_id]["update_born"])
            update_born = int(orig_genotypes[gid]["update_born"])
            if update_born > cur_extant_update_born:
                extant_genotype_id = gid

        # Summarize run information
        # - Extant summary info
        run_summary_info["extant_num_viability_sites"] = genotype_summary_info[extant_genotype_id]["viability_sites"]
        run_summary_info["extant_num_task_sites"] = genotype_summary_info[extant_genotype_id]["task_sites"]
        run_summary_info["extant_num_multi_task_sites"] = genotype_summary_info[extant_genotype_id]["multi_task_sites"]
        run_summary_info["extant_num_task_env_a_sites"] = genotype_summary_info[extant_genotype_id]["task_env_a_sites"]
        run_summary_info["extant_num_task_env_b_sites"] = genotype_summary_info[extant_genotype_id]["task_env_b_sites"]
        run_summary_info["extant_num_tasks_performed"] = len(orig_genotypes[extant_genotype_id]["tasks_performed"])
        ############################################################

        ############################################################
        # Contribute to lineage architecture (long-form) file
        # orig_genotypes
        # genotype_summary_info
        if agg_arch:
            for genotype_id in genotype_ids:
                orig_genotype = orig_genotypes[genotype_id]
                genotype_info = genotype_summary_info[genotype_id]
                tree_depth = orig_genotype["tree_depth"]
                position_info = genotype_info["position_info"]
                positions = [int(pos) for pos in position_info.keys()]
                positions.sort()
                for pos in positions:
                    line_info = {field:run_summary_info[field] for field in lineage_arch_long_cfg_fields}
                    cyclic_category = ""
                    if position_info[pos]["viability"]:
                        cyclic_category = "viability"
                    elif position_info[pos]["encodes_env_a_task"] and position_info[pos]["encodes_env_b_task"]:
                        cyclic_category = "env_ab"
                    elif position_info[pos]["encodes_env_a_task"]:
                        cyclic_category = "env_a"
                    elif position_info[pos]["encodes_env_b_task"]:
                        cyclic_category = "env_b"
                    else:
                        cyclic_category = "neutral"
                    tasks = "[" + ";".join(task for task in tasks_primary if task in position_info[pos]["tasks"]) + "]"
                    viability = position_info[pos]["viability"]
                    line_info["genotype_id"] = genotype_id
                    line_info["tree_depth"] = tree_depth
                    line_info["site"] = pos
                    line_info["cyclic_category"] = cyclic_category
                    line_info["tasks"] = tasks
                    line_info["viability"] = viability
                    lineage_arch_long_info.append(line_info)

            lineage_arch_long_fields = list(lineage_arch_long_info[0].keys())
            lineage_arch_long_fields.sort()
            write_header = False
            if lineage_arch_long_header == None:
                write_header = True
                lineage_arch_long_header = lineage_arch_long_fields
            elif lineage_arch_long_fields != lineage_arch_long_header:
                print("Lineage architecture long format header mismatch!")
                exit(-1)
            with open(lineage_arch_long_fpath, "a") as fp:
                if write_header: fp.write(",".join(lineage_arch_long_header) + "\n")
                fp.write("\n".join([ ",".join( [str(info[field]) for field in lineage_arch_long_header] ) for info in lineage_arch_long_info ]))
            lineage_arch_long_info = []
        ############################################################


        ############################################################
        # Add lines to run summary
        run_summary_fields = sorted(list(run_summary_info.keys()))
        if run_summary_header == None:
            run_summary_header = run_summary_fields
        elif run_summary_header != run_summary_fields:
            print("Run summary header mismatch!")
            exit(-1)
        run_summary_line = [str(run_summary_info[field]) for field in run_summary_fields]
        run_summary_lines.append(",".join(run_summary_line))
        ############################################################

    ############################################################
    # Write out run summary data
    fname = "knockouts_run_summary.csv"
    with open(os.path.join(dump_dir, fname), "w") as fp:
        out_content = ",".join(run_summary_header) + "\n" + "\n".join(run_summary_lines)
        fp.write(out_content)
    run_summary_lines = None
    ############################################################


if __name__ == "__main__":
    main()
