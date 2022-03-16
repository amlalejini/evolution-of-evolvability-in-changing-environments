'''
Aggregate mutational landscaping data
'''

import argparse, os, sys, pathlib, itertools, math, pandas
from typing import FrozenSet
from scipy.stats import entropy
sys.path.append(os.path.join(pathlib.Path(os.path.dirname(os.path.abspath(__file__))).parents[2], "scripts"))
import utilities as utils

run_identifier = "RUN_"

tasks_primary = ["not","nand","and","ornot","or","andnot"]
tasks_env_a = {"not", "and", "or"}
tasks_env_b = {"nand", "ornot", "andnot"}
tasks_all_pairs = [frozenset(pair) for pair in itertools.combinations(tasks_primary, 2)]
# print([pair for pair in itertools.product(tasks_primary, tasks_primary)])
task_id_map = {tasks_primary[i]:i for i in range(len(tasks_primary))}

profile_env_a = "101010"
profile_env_b = "010101"
profile_all = "111111"

treatment_identifiers = [
    "env_condition",
    "env_type",
    "env_chg_rate",
    "COPY_MUT_PROB"
]

task_cooccur_identifiers = [
    "env_condition",
    "env_type",
    "env_chg_rate",
    "COPY_MUT_PROB",
    "RANDOM_SEED"
]

run_command_excludes = {}

def main():
    parser = argparse.ArgumentParser(description="Aggregate mutation landscape data.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--mutant_file", type=str, default="mutants.dat", help="Name of the dat file that contains the mutants to analyze.")
    parser.add_argument("--output_id", type=str, default="", help="This identifier will be appended to the end of output files.")

    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump

    mutant_file = args.mutant_file
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

    # Create data structures to store lines to be saved to file.
    run_summary_header = None
    run_summary_content_lines = []

    treatment_cooccurrence_summary_info = {}      # {treatment_id: [{replicate_info}, {replicate_info}, ...]}
    treatment_id_cfg_map = {}

    task_cooccurrence_header = None
    task_cooccurrence_content_lines = []

    # Loop over runs, aggregating data from each.
    total_runs = len(run_dirs)
    cur_run_i = 0
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)

        # Skip over (but make note of) incomplete runs.
        if not os.path.exists(os.path.join(run_path, 'data', 'analysis')):
            print('Skipping: ', run_path)
            continue

        run_summary_info = {}                   # Hold summary information about run. Indexed by field.
        task_cooccurrence_info = []             # List of dictionaries

        cur_run_i += 1
        print(f"Processing ({cur_run_i}/{total_runs}): {run_path}")

        ############################################################
        # Extract commandline configuration settings (from cmd.log file)
        cmd_log_path = os.path.join(run_path, "cmd.log")
        cmd_params = utils.extract_params_cmd_log(cmd_log_path)

        for field in cmd_params:
            if field in run_command_excludes: continue
            run_summary_info[field] = cmd_params[field]

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

        # Figure out which treatment this run should contribute information to
        treatment_id = "__".join([f"{field}-{run_summary_info[field]}" for field in treatment_identifiers])
        if not treatment_id in treatment_cooccurrence_summary_info:
            treatment_cooccurrence_summary_info[treatment_id] = {pair:[] for pair in itertools.product(tasks_primary, tasks_primary)}
            treatment_id_cfg_map[treatment_id] = {field:run_summary_info[field] for field in treatment_identifiers}
        ############################################################

        ############################################################
        # Load 1-step mutants
        mutant_data = utils.read_avida_dat_file(os.path.join(run_path, "data", mutant_file))
        total_mutants = len(mutant_data)
        num_viable = 0
        viable_phenotypes = []
        task_occurrence_counts = {task: 0 for task in tasks_primary}
        task_cooccurrence_counts = {frozenset(pair): 0 for pair in itertools.combinations(tasks_primary, 2)}
        for mutant in mutant_data:
            viable = int(mutant["is_viable_(0/1)"])
            # If this mutant isn't viable, skip it.
            if viable == 0: continue
            # Otherwise, increment number of viable mutants.
            num_viable += 1
            # Count occurrence and co-occurrence
            for task_i in range(0, len(tasks_primary)):
                task_i_name = tasks_primary[task_i]
                task_i_occurs = int(mutant[task_i_name]) > 0
                task_occurrence_counts[task_i_name] += int(task_i_occurs)
                for task_j in range(task_i+1, len(tasks_primary)):
                    if task_i == task_j: continue
                    task_j_name = tasks_primary[task_j]
                    task_j_occurs = int(mutant[task_j_name]) > 0
                    task_cooccurrence_counts[frozenset({task_i_name, task_j_name})] += int( task_i_occurs and task_j_occurs )
            phenotype = "".join([str(int(mutant[task]) > 0) for task in tasks_primary])
            viable_phenotypes.append(phenotype)

        mutant_data = None

        # Calculate probabilities of each task appearing in viable mutants.
        task_occurrence_probs = {task: (task_occurrence_counts[task] / num_viable) if num_viable > 0 else 0 for task in tasks_primary}

        # Backfill counts of a task appearing with itself (always).
        for task in tasks_primary:
            task_cooccurrence_counts[frozenset({task})] = task_occurrence_counts[task]

        # Calculate joint probabilities of each task pair appearing in viable
        task_cooccurrence_probs = {pair:0 for pair in task_cooccurrence_counts}
        for pair in task_cooccurrence_counts:
            task_cooccurrence_probs[pair] = task_cooccurrence_counts[pair] / num_viable if num_viable > 0 else 0

        # Calculate pointwise mutual information: log2( (p(t1, t2)) / ( p(t1)*p(t2) ) )
        task_pmi = {pair:0 for pair in task_cooccurrence_probs}
        for pair in task_pmi:
            task_i, task_j = list(pair) if len(pair) > 1 else (list(pair)[0], list(pair)[0])
            # If task_i == task_j, don't multiply probabilities together
            denom = task_occurrence_probs[task_i] * task_occurrence_probs[task_j] if len(pair) > 1 else task_occurrence_probs[task_i]
            task_pmi[pair] = task_cooccurrence_probs[pair] / denom if denom != 0 else 0
            if task_pmi[pair] != 0:
                task_pmi[pair] = math.log(task_pmi[pair], 2)

        # Calculate join self-information: -log2( p(t1, t2) )
        task_jsi = {pair: 0 for pair in task_cooccurrence_probs}
        for pair in task_jsi:
            cooccur_prob = task_cooccurrence_probs[pair]
            jsi = -1 * math.log(cooccur_prob, 2) if cooccur_prob != 0 else 0
            task_jsi[pair] = jsi

        # Calculate normalized pmi
        task_npmi = {pair: 0 for pair in task_cooccurrence_probs}
        for pair in task_npmi:
            pmi = task_pmi[pair]
            jsi = task_jsi[pair]
            task_npmi[pair] = pmi / jsi if jsi != 0 else 0

        # Expected pmi & unexpected pmi
        expected_pairs = [frozenset(pair) for pair in itertools.combinations(tasks_env_a, 2)] + [frozenset(pair) for pair in itertools.combinations(tasks_env_b, 2)]
        # print(f"expected pairs: {expected_pairs}")
        unexpected_pairs = [frozenset(pair) for pair in itertools.product(tasks_env_a, tasks_env_b)]
        # print(f"unexpected pairs: {unexpected_pairs}")
        expected_npmi = sum([task_npmi[pair] for pair in expected_pairs])
        unexpected_npmi = sum([-1 * task_npmi[pair] for pair in unexpected_pairs])
        total_npmi = expected_npmi + unexpected_npmi

        # Calculate number of unique viable phenotypes
        viable_phenotypes_set = set(viable_phenotypes)
        num_unique_viable_phenotypes = len(viable_phenotypes_set)

        # Calculate entropy of viable phenotypes
        viable_phenotypes_series = pandas.Series(viable_phenotypes)
        viable_phenotypes_counts = viable_phenotypes_series.value_counts()
        viable_phenotypes_entropy = entropy(viable_phenotypes_counts, base=2)

        # print(f"  Viability: {num_viable}/{total_mutants}")
        # print(f"  Num phenotypes: {num_unique_viable_phenotypes}")
        # print(f"  Entropy of phenotypes: {viable_phenotypes_entropy}")
        # print(f"  Task occurrences: {task_occurrence_counts}")
        # print(f"  Task probabilities: {task_occurrence_probs}")
        # print(f"  Task co-occurrences: {task_cooccurrence_counts}")
        # print(f"  Task co-occurrences: {task_cooccurrence_probs}")
        # print(f"  Task pmi: {task_pmi}")
        # print(f"  Task jsi: {task_jsi}")
        # print(f"  Task jsi: {task_npmi}")

        # Save run summary info
        run_summary_info["total_mutants"] = total_mutants
        run_summary_info["num_viable"] = num_viable
        run_summary_info["prop_viable"] = num_viable / total_mutants
        run_summary_info["num_unique_viable_phenotypes"] = num_unique_viable_phenotypes
        run_summary_info["entropy_viable_phenotypes"] = viable_phenotypes_entropy
        run_summary_info["expected_npmi"] = expected_npmi
        run_summary_info["unexpected_npmi"] = unexpected_npmi
        run_summary_info["total_npmi"] = total_npmi

        # Save task cooccurrence (not including the full matrix currently)
        for pair in task_pmi:
            # for each pair, add a line to task_cooccurrence_info
            task_1, task_2 = list(pair) if len(pair) > 1 else (list(pair)[0], list(pair)[0])
            pair_info = {field:run_summary_info[field] for field in task_cooccur_identifiers}
            pair_info["task_1_id"] = task_id_map[task_1]
            pair_info["task_2_id"] = task_id_map[task_2]
            pair_info["task_1_name"] = task_1
            pair_info["task_2_name"] = task_2
            pair_info["pmi"] = task_pmi[pair]
            pair_info["npmi"] = task_npmi[pair]
            pair_info["joint_prob"] = task_cooccurrence_probs[pair]
            pair_info["joint_count"] = task_cooccurrence_counts[pair]
            pair_info["task_1_prob"] = task_occurrence_probs[task_1]
            pair_info["task_2_prob"] = task_occurrence_probs[task_2]
            pair_info["task_1_count"] = task_occurrence_counts[task_1]
            pair_info["task_2_count"] = task_occurrence_counts[task_2]
            task_cooccurrence_info.append(pair_info)

        # Save task co-occurrence for treatment summary (include full pairwise matrix)
        for task_ij in itertools.product(tasks_primary, tasks_primary):
            task_1, task_2 = task_ij
            pair = frozenset(task_ij)
            pair_info = {}
            pair_info["task_1_id"] = task_id_map[task_1]
            pair_info["task_2_id"] = task_id_map[task_2]
            pair_info["task_1_name"] = task_1
            pair_info["task_2_name"] = task_2
            pair_info["task_1_prob"] = task_occurrence_probs[task_1]
            pair_info["task_2_prob"] = task_occurrence_probs[task_2]
            pair_info["task_1_count"] = task_occurrence_counts[task_1]
            pair_info["task_2_count"] = task_occurrence_counts[task_2]
            pair_info["pmi"] = task_pmi[pair]
            pair_info["npmi"] = task_npmi[pair]
            pair_info["joint_prob"] = task_cooccurrence_probs[pair]
            pair_info["joint_count"] = task_cooccurrence_counts[pair]
            treatment_cooccurrence_summary_info[treatment_id][task_ij].append(pair_info)
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
        run_summary_content_lines.append(",".join(run_summary_line))
        ############################################################

        ############################################################
        # Add lines to task co-occurrence
        task_cooccurrence_fields = sorted(list(task_cooccurrence_info[0]))
        if task_cooccurrence_header == None:
            task_cooccurrence_header = task_cooccurrence_fields
        elif task_cooccurrence_header != task_cooccurrence_fields:
            print("Task co-occurrence header mismatch!")
            exit(-1)
        for info in task_cooccurrence_info:
            line = [str(info[field]) for field in task_cooccurrence_fields]
            task_cooccurrence_content_lines.append(",".join(line))
        ############################################################

    ############################################################
    # Write out run summary data
    fname = f"landscape_run_summary_{output_id}.csv" if output_id != "" else "landscape_run_summary.csv"
    with open(os.path.join(dump_dir, fname), "w") as fp:
        out_content = ",".join(run_summary_header) + "\n" + "\n".join(run_summary_content_lines)
        fp.write(out_content)
    run_summary_content_lines = None
    ############################################################

    ############################################################
    # Write out task co-occurrence data
    fname = f"landscape_task_cooccurrence_{output_id}.csv" if output_id != "" else "landscape_task_cooccurrence.csv"
    with open(os.path.join(dump_dir, fname), "w") as fp:
        out_content = ",".join(task_cooccurrence_header) + "\n" + "\n".join(task_cooccurrence_content_lines)
        fp.write(out_content)
    task_cooccurrence_content_lines = None
    ############################################################

    ############################################################
    # Summarize treatments
    # NOTE: The way that I organized this data ended up being a little awkward. Not sure I would implement the treatment summary in the same way if I went back and did it again.
    treatment_summary_lines = []
    treatment_summary_header = None
    num_treatments = len(treatment_cooccurrence_summary_info)
    print(f"Found {num_treatments} treatment(s) from {len(run_dirs)} runs.")
    task_pairs = tasks_all_pairs + [frozenset({task}) for task in tasks_primary]
    treatments = [treatment_id for treatment_id in treatment_cooccurrence_summary_info]
    for treatment_id in treatments:
        for pair in itertools.product(tasks_primary, tasks_primary):
            # One line for this pair for this treatment
            task_1, task_2 = pair
            info = {field:treatment_id_cfg_map[treatment_id][field] for field in treatment_id_cfg_map[treatment_id]}
            info["task_1_id"] = task_id_map[task_1]
            info["task_2_id"] = task_id_map[task_2]
            info["task_1_name"] = task_1
            info["task_2_name"] = task_2
            avg_fields = [
                "task_1_prob",
                "task_2_prob",
                "task_1_count",
                "task_2_count",
                "pmi",
                "npmi",
                "joint_prob",
                "joint_count"
            ]
            for field in avg_fields: info[field+"_avg"] = []
            for rep in treatment_cooccurrence_summary_info[treatment_id][pair]:
                for field in avg_fields:
                    info[field+"_avg"].append(float(rep[field]))
            for field in avg_fields:
                info[field+"_avg"] = float(sum(info[field+"_avg"])) / len(info[field+"_avg"])
            fields = sorted(list(info.keys()))
            if treatment_summary_header == None:
                treatment_summary_header = fields
            elif treatment_summary_header != fields:
                print("Treatment co-occurrence summary header mismatch!")
                exit(-1)
            treatment_summary_lines.append(",".join([str(info[field]) for field in treatment_summary_header]))

    fname = f"landscape_task_cooccurrence_treatment_summary_{output_id}.csv" if output_id != "" else "landscape_task_cooccurrence_treatment_summary.csv"
    with open(os.path.join(dump_dir, fname), "w") as fp:
        out_content = ",".join(treatment_summary_header) + "\n" + "\n".join(treatment_summary_lines)
        fp.write(out_content)
    treatment_summary_lines = None
    ############################################################

if __name__ == "__main__":
    main()