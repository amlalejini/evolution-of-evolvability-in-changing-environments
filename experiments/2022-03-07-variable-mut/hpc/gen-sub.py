'''
Generate slurm job submission scripts - one per condition
'''

import argparse, os, sys, errno, subprocess, csv, copy
from asyncio import base_events
from email.policy import default
from pyvarco import CombinationCollector

import markov_events

default_seed_offset = 2000
default_account = "devolab"
default_num_replicates = 30
default_job_time_request = "120:00:00"
default_job_mem_request = "4G"

total_updates = 300000
job_name = "03-07"
executable = "avida"

base_script_filename = './base_script.txt'
base_events_filename = "./base_events_file.cfg"

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
    "ANALYZE_FILE": "analyze.cfg"
}

special_decorators = ["__DYNAMIC", "__COPY_OVER"]

combos.register_var("environment__DYNAMIC")
combos.register_var("COPY_MUT_PROB")


combos.add_val(
    "COPY_MUT_PROB",
    [
        "0.0001",
        "0.001",
        "0.0316"
    ]
)

combos.add_val(
    "environment__DYNAMIC",
    [
        "env-cycling_rate-30",
        "env-cycling_rate-300",
        "env-cycling_rate-3000",
        "env-random_rate-30",
        "env-random_rate-300",
        "env-random_rate-3000",
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
        file_str = file_str.replace("<<ARRAY_ID_RANGE>>", f"1-{num_replicates}")
        file_str = file_str.replace("<<MEMORY_REQUEST>>", job_memory_request)
        file_str = file_str.replace("<<JOB_NAME>>", job_name)
        file_str = file_str.replace("<<CONFIG_DIR>>", config_dir)
        file_str = file_str.replace("<<REPO_DIR>>", repo_dir)
        file_str = file_str.replace("<<EXEC>>", executable)
        file_str = file_str.replace("<<JOB_SEED_OFFSET>>", str(cur_seed))
        file_str = file_str.replace("<<ACCOUNT_NAME>>", hpc_account)

        # ===================================================
        # Configure the run
        # ===================================================
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
        ###################################################################

        run_param_info["EVENT_FILE"] = events_file_arg

        fields = list(run_param_info.keys())
        fields.sort()
        set_params = [f"-set {field} {run_param_info[field]}" for field in fields]
        copy_params = [condition_dict[key] for key in condition_dict if "__COPY_OVER" in key]

        run_params = " ".join(set_params + copy_params)

        # Add run commands to run the experiment
        run_commands = ''
        # Copy events description and events file into run directory
        run_commands += "cp ${CONFIG_DIR}/events/" + events_file_arg + " ${RUN_DIR}/\n"
        run_commands += "cp ${CONFIG_DIR}/environment_states/" + events_desc_arg + " ${RUN_DIR}/environment_states.csv\n"
        # Run avida
        run_commands += f'RUN_PARAMS="{run_params}"\n'

        #####################################################################
        # -- Commenting out the bits that actually run Avida
        run_commands += '# echo "./${EXEC} ${RUN_PARAMS}" > cmd.log\n'
        run_commands += '# ./${EXEC} ${RUN_PARAMS} > run.log\n'
        #####################################################################

        file_str = file_str.replace("<<RUN_COMMANDS>>", run_commands)

        # ===================================================
        # Configure analyze mode
        # ===================================================
        # Add commands to run avida analyze mode
        analysis_commands = ''
        # analysis_commands += f'RUN_PARAMS="{run_params}"\n' # Only include if analyze mode parameters don't match run paramters.

        #####################################################################
        # -- Commenting out the bits that actually run Avida
        analysis_commands += '# ./${EXEC} ${RUN_PARAMS} -a\n'
        #####################################################################

        gen_1step_mutants_cmd = '# python scripts/gen-mutants.py --steps 1 --analysis_output mutants_step-1.dat --inst_set ${CONFIG_DIR}/instset-heads.cfg --input ${RUN_DIR}/data/analysis/final_dominant.dat --dump ${RUN_DIR} --avida_args "${RUN_PARAMS}" --num_tasks 6 --run_avida_analysis --run_dir ${RUN_DIR}'
        gen_2step_mutants_cmd = '# python scripts/gen-mutants.py --steps 2 --analysis_output mutants_step-2.dat --inst_set ${CONFIG_DIR}/instset-heads.cfg --input ${RUN_DIR}/data/analysis/final_dominant.dat --dump ${RUN_DIR} --avida_args "${RUN_PARAMS}" --num_tasks 6 --run_avida_analysis --run_dir ${RUN_DIR}'
        knockout_cmd = 'python scripts/gen-knockouts.py --inst_set ${CONFIG_DIR}/instset-heads-ko.cfg --input ${RUN_DIR}/data/analysis/lineage.dat --num_tasks 6 --avida_args "${RUN_PARAMS}" --run_dir ${RUN_DIR} --output knockouts.csv --cleanup'

        analysis_commands += "cd ${REPO_DIR}\n"
        analysis_commands += gen_1step_mutants_cmd + "\n"
        analysis_commands += gen_2step_mutants_cmd + "\n"
        analysis_commands += knockout_cmd + "\n"
        analysis_commands += "cd ${RUN_DIR}\n"

        file_str = file_str.replace("<<ANALYSIS_COMMANDS>>", analysis_commands)

        # ===================================================
        # Write job submission file
        # ===================================================
        mkdir_p(job_dir)
        with open(os.path.join(job_dir, f'{filename_prefix}.sb'), 'w') as fp:
            fp.write(file_str)

        cur_job_id += 1
        cond_i += 1

if __name__ == "__main__":
    main()
