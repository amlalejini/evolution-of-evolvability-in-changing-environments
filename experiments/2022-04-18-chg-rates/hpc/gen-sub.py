'''
Generate slurm job submission scripts - one per condition
'''

import argparse, os, sys, errno, subprocess, csv, copy
from asyncio import base_events
from email.policy import default
from pyvarco import CombinationCollector

import markov_events

default_seed_offset = 1000
default_account = "devolab"
default_num_replicates = 30
default_job_time_request = "120:00:00"
default_job_mem_request = "4G"
default_total_updates = 300000

job_name = "04-18"
executable = "avida"

base_script_filename = "./base_script.txt"
base_events_filename = "./base_events_file.cfg"
# base_run_logic="""
# if [[ ${SLURM_ARRAY_TASK_ID} -eq <<RUN_ARRAY_ID>> ]] ; then

# fi
# """

# Create combo object to collect all conditions we'll run
combos = CombinationCollector()

fixed_parameters = {
    "DIVIDE_INS_PROB": "0.0",
    "DIVIDE_DEL_PROB": "0.0",
    "OFFSPRING_SIZE_RANGE": "1.0",
    "STERILIZE_UNSTABLE": "1.0",
    "COSTLY_HEAD_COPY": "0",
    "MAX_HEAD_COPY_COST": "0.0",
    "BIRTH_METHOD": "4",
    "WORLD_X": "150",
    "WORLD_Y": "150",
    "META_COPY_MUT": "1.0",
    "META_STD_DEV": "0.01",
    "MUT_RATE_SOURCE": "2",
    # "PHYLOGENY_SNAPSHOT_RES": "300000",
    # "SYSTEMATICS_RES": "100",
    # "FORCE_MRCA_COMP": "0",
    "ANALYZE_FILE": "analyze.cfg",
}

special_decorators = ["__DYNAMIC", "__COPY_OVER"]

combos.register_var("environment__DYNAMIC")
combos.register_var("COPY_MUT_PROB")

combos.add_val(
    "COPY_MUT_PROB",
    [
        "0.0001",
        "0.001",
        "0.01"
    ]
)

combos.add_val(
    "environment__DYNAMIC",
    [
        "env-cycling_rate-50",
        "env-cycling_rate-100",
        "env-cycling_rate-200",
        "env-cycling_rate-300",
        "env-cycling_rate-500",
        "env-cycling_rate-1000",
        "env-cycling_rate-2000",
        "env-cycling_rate-3000",
        "env-cycling_rate-5000",

        "env-random_rate-50",
        "env-random_rate-100",
        "env-random_rate-200",
        "env-random_rate-300",
        "env-random_rate-500",
        "env-random_rate-1000",
        "env-random_rate-2000",
        "env-random_rate-3000",
        "env-random_rate-5000",

        "env-constA",
        "env-constB"
    ]
)


'''
This is functionally equivalent to the mkdir -p [fname] bash command
'''
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the output directory for phase one of each run?")
    parser.add_argument("--config_dir", type=str, help="Where is the configuration directory for experiment?")
    parser.add_argument("--repo_dir", type=str, help="Where is the repository for this experiment?")
    parser.add_argument("--job_dir", type=str, default=None, help="Where to output these job files? If none, put in 'jobs' directory inside of the data_dir")
    parser.add_argument("--replicates", type=int, default=default_num_replicates, help="How many replicates should we run of each condition?")
    parser.add_argument("--seed_offset", type=int, default=default_seed_offset, help="Value to offset random number seeds by")
    parser.add_argument("--force_events_gen", action="store_true", help="Should we force events file generation even if files already exist?")
    parser.add_argument("--account", type=str, default=default_account, help="Value to use for the slurm ACCOUNT")
    parser.add_argument("--time_request", type=str, default=default_job_time_request, help="How long to request for each job on hpc?")
    parser.add_argument("--mem", type=str, default=default_job_mem_request, help="How much memory to request for each job?")
    parser.add_argument("--updates", type=int, default=default_total_updates, help="How many updates to run avida?")

    parser.add_argument("--patch", action="store_true", help="Patch unfinished runs. Only generate submission scripts for incomplete runs.")

    # Load in command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    config_dir = args.config_dir
    job_dir = args.job_dir
    repo_dir = args.repo_dir
    num_replicates = args.replicates
    hpc_account = args.account
    seed_offset = args.seed_offset
    force_events_gen = args.force_events_gen
    job_time_request = args.time_request
    job_memory_request = args.mem
    patch_mode = args.patch
    total_updates = args.updates

    # Load in the base slurm file
    base_sub_script = ""
    with open(base_script_filename, 'r') as fp:
        base_sub_script = fp.read()

    # Load in the base events file
    base_events_content = ""
    with open(base_events_filename, "r") as fp:
        base_events_content = fp.read()

    # Get list of all combinations to run
    combo_list = combos.get_combos()

    # Calculate how many jobs we have, and what the last id will be
    num_jobs = num_replicates * len(combo_list)
    print(f'Generating {num_jobs} across {len(combo_list)} files!')
    print(f' - Data directory: {data_dir}')
    print(f' - Config directory: {config_dir}')
    print(f' - Repository directory: {repo_dir}')
    print(f' - Replicates: {num_replicates}')
    print(f' - Account: {hpc_account}')
    print(f' - Time Request: {job_time_request}')
    print(f' - Seed offset: {seed_offset}')
    print(f' - Force events generation: {force_events_gen}')
    print(f' - Total updates? {total_updates}')
    print(f' - Patch mode? {patch_mode}')

    # If no job_dir provided, default to data_dir/jobs
    if job_dir == None:
        job_dir = os.path.join(data_dir, "jobs")

    # Create job file for each condition
    cur_job_id = 0
    cond_i = 0
    generated_files = set()
    for condition_dict in combo_list:
        cur_seed = seed_offset + (cur_job_id * num_replicates)
        filename_prefix = f'RUN_C{cond_i}'
        file_str = base_sub_script
        file_str = file_str.replace("<<TIME_REQUEST>>", job_time_request)
        file_str = file_str.replace("<<MEMORY_REQUEST>>", job_memory_request)
        file_str = file_str.replace("<<JOB_NAME>>", job_name)
        file_str = file_str.replace("<<CONFIG_DIR>>", config_dir)
        file_str = file_str.replace("<<REPO_DIR>>", repo_dir)
        file_str = file_str.replace("<<EXEC>>", executable)
        file_str = file_str.replace("<<JOB_SEED_OFFSET>>", str(cur_seed))
        file_str = file_str.replace("<<ACCOUNT_NAME>>", hpc_account)

        ###################################################################
        # Configure the run
        ###################################################################
        file_str = file_str.replace("<<RUN_DIR>>", \
            os.path.join(data_dir, f'{filename_prefix}_'+'${SEED}'))

        # Format commandline arguments for the run
        run_param_info = {key:condition_dict[key] for key in condition_dict if not any([dec in key for dec in special_decorators])}
        # Add fixed paramters
        for param in fixed_parameters:
            if param in run_param_info: continue
            run_param_info[param] = fixed_parameters[param]
        # Set random number seed
        run_param_info["RANDOM_SEED"] = '${SEED}'

        # Handle dynamic params (will vary experiment by experiment)
        environment = condition_dict["environment__DYNAMIC"]
        env_cfg = {cfg.split("-")[0]:cfg.split("-")[1] for cfg in environment.split("_")}

        ###################################################################
        # -- Build environment file (for each replicate). --
        ###################################################################
        events_desc_path = os.path.join(config_dir, "environment_states")
        events_path = os.path.join(config_dir, "events")
        events_file_arg = ""
        events_desc_arg = ""
        for ri in range(num_replicates):
            replicate_seed = cur_seed + ri

            events_fname = ""
            events_desc_fname = ""

            if "const" in env_cfg["env"]:
                # Constant environments can share a single events file.
                events_fname = f"events_{environment}.cfg"
                events_desc_fname = f"{environment}.csv"
                events_file_arg = events_fname
                events_desc_arg = events_desc_fname
            else:
                # Changing environments need to be run-specific.
                events_fname = f"events_{environment}_run-{replicate_seed}.cfg"
                events_desc_fname = f"{environment}_run-{replicate_seed}.csv"
                events_file_arg = "events_"+environment+"_run-${SEED}.cfg"
                events_desc_arg = environment+"_run-${SEED}.csv"

            events_fpath = os.path.join(events_path, events_fname)
            events_desc_fpath = os.path.join(events_desc_path, events_desc_fname)
            if events_fpath in generated_files: continue

            # If the events file doesn't exist or we're forcing events file generation,
            # generate a new events file.
            if ( not os.path.exists(events_fpath) ) or ( not os.path.exists(events_desc_fpath) ) or ( force_events_gen ):
                mkdir_p(events_desc_path)
                mkdir_p(events_path)

                # Start with the base events content.
                events_file_content = base_events_content
                events_file_content = events_file_content.replace("<<TOTAL_UPDATES>>", str(total_updates))

                events_gen_results = None
                if env_cfg["env"] == "constA":
                    events_gen_results=markov_events.gen_const_env_events(total_updates, "A")
                elif env_cfg["env"] == "constB":
                    events_gen_results=markov_events.gen_const_env_events(total_updates, "B")
                elif env_cfg["env"] == "cycling":
                    events_gen_results=markov_events.gen_cyclic_env_events(
                        total_updates=total_updates,
                        cycle_period=int(env_cfg["rate"]),
                        seed=replicate_seed
                    )
                elif env_cfg["env"] == "random":
                    events_gen_results=markov_events.gen_random_env_events(
                        total_updates=total_updates,
                        chg_period=int(env_cfg["rate"]),
                        seed=replicate_seed
                    )
                events_file_content += "\n"
                events_file_content += events_gen_results["event_content"]

                events_desc_fields = list(events_gen_results["event_desc"][0].keys())
                events_desc_fields.sort()
                events_desc_header = ",".join(events_desc_fields)
                events_desc_content = f"{events_desc_header}\n" + "\n".join([",".join([str(desc[field]) for field in events_desc_fields]) for desc in events_gen_results["event_desc"]])
                with open(events_desc_fpath, "w") as fp:
                    fp.write(events_desc_content)
                with open(events_fpath, "w") as fp:
                    fp.write(events_file_content)
                generated_files.add(events_fpath)

        # Set the event file accordingly
        run_param_info["EVENT_FILE"] = events_file_arg
        ###################################################################

        ###################################################################
        # Build avida commandline parameters string
        ###################################################################
        fields = list(run_param_info.keys())
        fields.sort()
        set_params = [f"-set {field} {run_param_info[field]}" for field in fields]
        copy_params = [condition_dict[key] for key in condition_dict if "__COPY_OVER" in key]
        run_params = " ".join(set_params + copy_params)
        ###################################################################

        # Add run commands to run the experiment
        cfg_run_commands = ''
        # Copy events description and events file into run directory
        cfg_run_commands += "cp ${CONFIG_DIR}/events/" + events_file_arg + " ${RUN_DIR}/\n"
        cfg_run_commands += "cp ${CONFIG_DIR}/environment_states/" + events_desc_arg + " ${RUN_DIR}/environment_states.csv\n"
        # Set the run
        cfg_run_commands += f'RUN_PARAMS="{run_params}"\n'

        # By default, add all commands to submission file.
        array_id_run_info = {
            array_id: {
                "avida": True,
                "avida_analyze_mode": True,
                "landscape_step-1": False,
                "landscape_step-2": False,
                "knockouts": False,
                "pairwise-knockouts": True
            }
            for array_id in range(1, num_replicates+1)
        }
        array_id_to_seed = {array_id:(cur_seed + (array_id - 1)) for array_id in array_id_run_info}
        if patch_mode:
            for array_id in array_id_run_info:
                replicate_seed = array_id_to_seed[array_id]
                run_path = os.path.join(data_dir, f"RUN_C{cond_i}_{replicate_seed}")
                array_id_run_info[array_id]["avida"] = not os.path.exists(os.path.join(run_path, "data", f"detail-{total_updates}.spop"))
                array_id_run_info[array_id]["avida_analyze_mode"] = not os.path.exists(os.path.join(run_path, "data", "analysis", "lineage.dat"))
                array_id_run_info[array_id]["landscape_step-1"] = not os.path.exists(os.path.join(run_path, "data", "mutants_step-1.dat"))
                array_id_run_info[array_id]["landscape_step-2"] = not os.path.exists(os.path.join(run_path, "data", "mutants_step-2.dat"))
                array_id_run_info[array_id]["knockouts"] = not os.path.exists(os.path.join(run_path, "data", "knockouts.csv"))
                array_id_run_info[array_id]["pairwise-knockouts"] = not os.path.exists(os.path.join(run_path, "data", "pairwise-knockouts.csv"))

        # Track which array ids need to be included. If none, don't need to output this file.
        active_array_ids = []
        inactive_array_ids = []
        run_sub_logic = ""
        for array_id in range(1, num_replicates+1):
            # If this run is totally done, make note and continue.
            if not any([array_id_run_info[array_id][field] for field in array_id_run_info[array_id]]):
                inactive_array_ids.append(array_id)
                continue
            # This run is not done already. Make note.
            active_array_ids.append(array_id)

            run_logic = "if [[ ${SLURM_ARRAY_TASK_ID} -eq "+str(array_id)+" ]] ; then\n"

            # (1) Run Avida
            run_commands = ''
            if array_id_run_info[array_id]["avida"]:
                run_commands += 'echo "./${EXEC} ${RUN_PARAMS}" > cmd.log\n'
                run_commands += './${EXEC} ${RUN_PARAMS} > run.log\n'
                run_commands += 'mv ./*.csv ./data/ \n'

            # (2) Run Avida in analyze mode
            analysis_commands = ''
            if array_id_run_info[array_id]["avida_analyze_mode"]:
                analysis_commands += './${EXEC} ${RUN_PARAMS} -a\n'

            # (3) Run mutational landscaping
            landscape_commands = ""
            landscape_commands += "cd ${REPO_DIR}\n"
            # -- step 1 mutants --
            if array_id_run_info[array_id]["landscape_step-1"]:
                gen_1step_mutants_cmd = 'python scripts/gen-mutants.py --steps 1 --analysis_output mutants_step-1.dat --inst_set ${CONFIG_DIR}/instset-heads.cfg --input ${RUN_DIR}/data/analysis/final_dominant.dat --dump ${RUN_DIR} --avida_args "${RUN_PARAMS}" --num_tasks 6 --run_avida_analysis --run_dir ${RUN_DIR}'
                landscape_commands += gen_1step_mutants_cmd + "\n"
            # -- step 2 mutants --
            if array_id_run_info[array_id]["landscape_step-2"]:
                gen_2step_mutants_cmd = 'python scripts/gen-mutants.py --steps 2 --analysis_output mutants_step-2.dat --inst_set ${CONFIG_DIR}/instset-heads.cfg --input ${RUN_DIR}/data/analysis/final_dominant.dat --dump ${RUN_DIR} --avida_args "${RUN_PARAMS}" --num_tasks 6 --run_avida_analysis --run_dir ${RUN_DIR}'
                landscape_commands += gen_2step_mutants_cmd + "\n"
            landscape_commands += "cd ${RUN_DIR}\n"

            # (4) Generate knockouts
            knockout_commands = ""
            knockout_commands += "cd ${REPO_DIR}\n"
            if array_id_run_info[array_id]["knockouts"]:
                knockout_commands += 'python scripts/gen-knockouts.py --inst_set ${CONFIG_DIR}/instset-heads-ko.cfg --input ${RUN_DIR}/data/analysis/lineage.dat --num_tasks 6 --avida_args "${RUN_PARAMS}" --run_dir ${RUN_DIR} --output knockouts.csv --cleanup'
                knockout_commands += "\n"
            if array_id_run_info[array_id]["pairwise-knockouts"]:
                knockout_commands += 'python scripts/gen-pairwise-knockouts.py --inst_set ${CONFIG_DIR}/instset-heads-ko.cfg --input ${RUN_DIR}/data/analysis/final_dominant.dat --num_tasks 6 --avida_args "${RUN_PARAMS}" --run_dir ${RUN_DIR} --output pairwise-knockouts.csv --cleanup'
                knockout_commands += "\n"
            knockout_commands += "cd ${RUN_DIR}\n"

            run_logic += run_commands
            run_logic += analysis_commands
            run_logic += landscape_commands
            run_logic += knockout_commands
            run_logic += "fi\n\n"
            run_sub_logic += run_logic

        # -- Set the SLURM array id range parameter --
        array_id_range_param = ""
        if len(active_array_ids) == num_replicates:
            array_id_range_param = f"1-{num_replicates}"
        else:
            array_id_range_param = ",".join([str(array_id) for array_id in active_array_ids])

        # -- add run commands to file str --
        file_str = file_str.replace("<<ARRAY_ID_RANGE>>", array_id_range_param)
        file_str = file_str.replace("<<CFG_RUN_COMMANDS>>", cfg_run_commands)
        file_str = file_str.replace("<<RUN_COMMANDS>>", run_sub_logic)

        ###################################################################
        # Write job submission file (if any of the array ids are active)
        ###################################################################
        # Report active/inactive
        print(f"RUN_C{cond_i}:")
        print(f" - Active: " + ", ".join([f"RUN_C{cond_i}_{array_id_to_seed[array_id]}" for array_id in active_array_ids]))
        print(f" - Inactive: " + ", ".join([f"RUN_C{cond_i}_{array_id_to_seed[array_id]}" for array_id in inactive_array_ids]))

        if len(active_array_ids):
            mkdir_p(job_dir)
            with open(os.path.join(job_dir, f'{filename_prefix}.sb'), 'w') as fp:
                fp.write(file_str)

        # Update condition id and current job id
        cur_job_id += 1
        cond_i += 1

if __name__ == "__main__":
    main()
