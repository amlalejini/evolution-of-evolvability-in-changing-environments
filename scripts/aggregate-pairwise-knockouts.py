'''
Process + aggregate pairwise knockout analysis information.
'''

import argparse, os, sys, errno, csv, pathlib
from operator import xor
sys.path.append(os.path.join(pathlib.Path(os.path.dirname(os.path.abspath(__file__))).parents[2], "scripts"))
import utilities as utils

run_identifier = "RUN_"

tasks_primary = ["not","nand","and","ornot","or","andnot"]
tasks_env_a = {"not", "and", "or"}
tasks_env_b = {"nand", "ornot", "andnot"}

profile_env_a = "101010"
profile_env_b = "010101"
profile_all = "111111"

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--ko_file", type=str, default="pairwise-knockouts.csv", help="Name of file where knockout information for run is stored.")
    parser.add_argument("--output_id", type=str, default="", help="String to append to output file names.")

    # Parse command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump
    ko_fname = args.ko_file
    output_id = args.output_id

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

    # Loop over runs, aggregating data from each.
    total_runs = len(run_dirs)
    cur_run_i = 0
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)
        # Skip over (but make note of) incomplete runs.
        if not os.path.exists(os.path.join(run_path, 'data', ko_fname)):
            print("Skipping: ", run_path)
            print(f" - Failed to find '{ko_fname}'")
            continue

        run_summary_info = {}

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
        # Read knockouts data
        knockouts_data = utils.read_csv(os.path.join(run_path, "data", ko_fname))

        # Collect genotype ids (i.e., how many genotype knockout analyses are represented in the data)
        genotype_ids = list({line["genotype_id"] for line in knockouts_data})

        # Collect the original genotypes (organized by genotype id)
        orig_genotypes = {line["genotype_id"]:line for line in knockouts_data if (line["ko_pos_0"] == "-1") or (line["ko_pos_1"] == "-1")}

        # Identify the extant genotype
        extant_genotype_id = None
        for gid in orig_genotypes:
            if extant_genotype_id == None:
                extant_genotype_id = gid
                continue
            cur_extant_update_born = int(orig_genotypes[extant_genotype_id]["update_born"])
            update_born = int(orig_genotypes[gid]["update_born"])
            if update_born > cur_extant_update_born:
                extant_genotype_id = gid

        # Only going to work with the extant genotype
        orig_genotype = orig_genotypes[extant_genotype_id]

        # Collect lines associated with focal (extant) knockout analysis
        knockouts = []
        for line in knockouts_data:
            # If line is not part of extant genotype analysis, skip
            if line["genotype_id"] != extant_genotype_id:
                continue
            # If line is the original genotype, skip
            if line["ko_pos_0"] == "-1" or line["ko_pos_1"] == "-1":
                continue
            knockouts.append(line)

        orig_tasks_performed = {task for task in tasks_primary if int(orig_genotype[task]) > 0}
        # orig_genotype["tasks_performed"] = orig_tasks_performed

        def p(a, b=None):
            if b==None:
                return frozenset({a})
            return frozenset({a,b})

        site_info = {}
        sites = set()
        site_pairs = set()
        # Analyze each knockout analysis (which may represent a single site or two sites knocked out).
        for ko in knockouts:
            ko_pos_0 = int(ko["ko_pos_0"])
            ko_pos_1 = int(ko["ko_pos_1"])
            pos_idx = p(ko_pos_0, ko_pos_1)
            # Add this position/pair to sites/site_pairs. Skip over if already processed.
            if ko_pos_0 == ko_pos_1:
                if pos_idx in sites: continue
                sites.add(pos_idx)
            else:
                if pos_idx in site_pairs: continue
                site_pairs.add(pos_idx)

            info = {} # Store information about this ko
            # What tasks does this knockout perform?
            ko_tasks_performed =  {task for task in tasks_primary if int(ko[task]) > 0}
            ko["tasks_performed"] = ko_tasks_performed
            # Are these position(s) important for viability?
            info["viability"] = not bool(int(ko["is_viable_(0/1)"]))
            # If ko is viable, what tasks (if any) are knocked out?
            info["ko_tasks"] = {}
            if ko["is_viable_(0/1)"] == "1":
                info["ko_tasks"] = orig_tasks_performed - ko_tasks_performed
            # Store information about this task globally
            site_info[pos_idx] = info

        # `task_sites` stores sites important for each task. Includes all sites marked as redundant.
        task_sites = {task:set() for task in tasks_primary}
        # Stores sites marked as redundant
        task_redundancy_sites = {task:set() for task in tasks_primary}
        # Stores sites that can 'recover' a lost task as a result of knocking out a different site
        task_recovery_sites = {task:set() for task in tasks_primary}

        # Stores sites important for viability
        viability_sites = set()
        viability_redundancy_sites = set()
        viability_recovery_sites = set()

        # Analyze single-site knockouts
        for site_idx in sites:
            assert len(site_idx) == 1, "Single-site index should be length-1."
            info = site_info[site_idx]
            site = list(site_idx)[0]
            # Is this site important for viability?
            if info["viability"]:
                viability_sites.add(site)
            # Which tasks is this site important for?
            # (note: ko_tasks won't have anything in it if this site is critical for viability)
            for task in info["ko_tasks"]:
                task_sites[task].add(site)

        # Analyse paired knockouts
        for site_idx in site_pairs:
            assert len(site_idx) == 2, "Paired site indices should be length-2."
            site_0, site_1 = list(site_idx)
            assert site_0 != site_1, "Paired sites should be different."
            info = site_info[site_idx]
            # Already identified important single-sites. Need to identify:
            # - Redundant sites
            # - Recovery sites

            # Check viability first. If recovers viability, does not count as recovering tasks.
            s0_viability = site_0 in viability_sites
            s1_viability = site_1 in viability_sites
            s0s1_viability = info["viability"]

            # Viability redundancy
            viability_redundancy = (not s0_viability) and (not s1_viability) and s0s1_viability
            if viability_redundancy:
                viability_redundancy_sites.add(site_0)
                viability_redundancy_sites.add(site_1)

            # Viability recovery
            recovers_s0_viability = (s0_viability) and (not s1_viability) and (not s0s1_viability)
            recovers_s1_viability = (not s0_viability) and (s1_viability) and (not s0s1_viability)
            recovers_s1s0_viability = (s0_viability) and (s1_viability) and (not s0s1_viability)
            recovers_viability = recovers_s0_viability or recovers_s1_viability or recovers_s1s0_viability
            if recovers_s0_viability: viability_recovery_sites.add(site_1)
            if recovers_s1_viability: viability_recovery_sites.add(site_0)
            if recovers_s1s0_viability:
                viability_recovery_sites.add(site_0)
                viability_recovery_sites.add(site_1)

            # If this pair recovers viability, ignore what it does for task phenotype
            if recovers_viability: continue

            # Check for task redundancy and task recovery
            for task in tasks_primary:
                s0_coding = site_0 in task_sites[task]
                s1_coding = site_1 in task_sites[task]
                s0s1_coding = task in info["ko_tasks"]
                # Check for redundancy
                task_redundancy = (not s0_coding) and (not s1_coding) and (s0s1_coding)
                if task_redundancy:
                    task_redundancy_sites[task].add(site_0)
                    task_redundancy_sites[task].add(site_1)
                # Check for task recovery
                recovers_s0_task = (s0_coding) and (not s1_coding) and (not s0s1_coding)
                recovers_s1_task = (not s0_coding) and (s1_coding) and (not s0s1_coding)
                recovers_s1s0_task = (s0_coding) and (s1_coding) and (not s0s1_coding)
                # recovers_task = recovers_s0_task or recovers_s1_task or recovers_s1s0_task
                if recovers_s0_task: task_recovery_sites[task].add(site_1)
                if recovers_s1_task: task_recovery_sites[task].add(site_0)
                if recovers_s1s0_task:
                    task_recovery_sites[task].add(site_0)
                    task_recovery_sites[task].add(site_1)

        # For each site, add:
        # - nonredundant_coded_tasks (which tasks does this site code for)
        # - all_coded_tasks (which tasks does this site code for)
        # - recovered_tasks
        # sites = sorted([ for site_idx in sites])
        info_by_site = {list(site_idx)[0]:site_info[site_idx] for site_idx in site_info}
        num_multi_task_sites = 0
        for site in info_by_site:
            all_coded_tasks = {task for task in tasks_primary if (site in task_sites[task]) or (site in task_redundancy_sites[task])}
            recovered_tasks = {task for task in tasks_primary if (site in task_recovery_sites[task])}
            num_multi_task_sites += int(len(all_coded_tasks) > 1)
            info_by_site[site]["all_coded_tasks"] = all_coded_tasks
            info_by_site[site]["recovered_tasks"] = recovered_tasks

        # todo - Add redundant sites to coding sites (for both tasks and viability)
        all_nonredundant_task_coding_sites = {site for task in task_sites for site in task_sites[task]}
        all_redundant_task_coding_sites = {site for task in task_redundancy_sites for site in task_redundancy_sites[task]}
        all_task_coding_sites = all_nonredundant_task_coding_sites.union(all_redundant_task_coding_sites)
        all_task_recovery_sites = {site for task in task_recovery_sites for site in task_recovery_sites[task]}

        all_viability_sites = viability_sites.union(viability_redundancy_sites)

        # -- info about original genotype --
        run_summary_info["num_tasks_performed"] = len(orig_tasks_performed)
        run_summary_info["tasks_performed"] = ";".join(sorted(orig_tasks_performed))
        # -- viability coding site info --
        run_summary_info["num_nonredundant_viability_sites"] = len(viability_sites)
        run_summary_info["num_redundant_viability_sites"] = len(viability_redundancy_sites)
        run_summary_info["num_viability_sites"] = len(all_viability_sites)
        run_summary_info["num_viability_recovery_sites"] = len(viability_recovery_sites)
        # -- task coding site info ---
        # num_task_nonredundant_sites (all sites important when knocked out alone)
        run_summary_info["num_nonredundant_task_sites"] = len(all_nonredundant_task_coding_sites)
        # num_task_redundancy_sites (all redundancy sites)
        run_summary_info["num_redundant_task_sites"] = len(all_redundant_task_coding_sites)
        # num_task_sites (all task sites + redundancy sites)
        run_summary_info["num_task_coding_sites"] = len(all_task_coding_sites)
        # num_task_recovery_sites
        run_summary_info["num_task_recovery_sites"] = len(all_task_recovery_sites)
        # num multi task coding sites
        run_summary_info["num_multi_task_sites"] = num_multi_task_sites

        # -- per-task coding site info --
        for task in tasks_primary:
            run_summary_info[f"num_{task}_task_coding_sites"] = len(task_sites[task].union(task_redundancy_sites[task]))
            run_summary_info[f"num_{task}_task_recovery_sites"] = len(task_recovery_sites[task])

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
    fname = "pairwise_knockouts_run_summary.csv" if output_id == "" else f"pairwise_knockouts_run_summary_{output_id}.csv"
    with open(os.path.join(dump_dir, fname), "w") as fp:
        out_content = ",".join(run_summary_header) + "\n" + "\n".join(run_summary_lines)
        fp.write(out_content)
    run_summary_lines = None
    ############################################################


if __name__ == "__main__":
    main()