'''
Generate slurm job submission scripts - one per condition
'''

import argparse, os, sys, errno, subprocess, csv
from email.policy import default
from pyvarco import CombinationCollector

default_seed_offset = 1000
default_account = "devolab"
default_num_replicates = 30

job_time_request = "120:00:00"
job_memory_request = "4G"
job_name = "03-02"
executable = "avida"

base_script_filename = './base_script.txt'

# Create combo object to collect all conditions we'll run
combos = CombinationCollector()

combos.register_var("EVENTS_FILE__COPY_OVER")
combos.register_var("COPY_MUT_PROB")
combos.register_var("COSTLY_COPY__COPY_OVER")
combos.register_var("WORLD_SIZE__COPY_OVER")
combos.register_var("FIXED_LENGTH__COPY_OVER")
combos.register_var("EVO_MUT_RATES__COPY_OVER")

combos.add_val(
    "FIXED_LENGTH__COPY_OVER",
    [
        "-set DIVIDE_INS_PROB 0.0 -set DIVIDE_DEL_PROB 0.0 -set OFFSPRING_SIZE_RANGE 1.0 -set STERILIZE_UNSTABLE 1.0"
    ]
)

combos.add_val(
    "WORLD_SIZE__COPY_OVER",
    [
        "-set WORLD_X 150 -set WORLD_Y 150"
    ]
)

combos.add_val(
    "EVO_MUT_RATES__COPY_OVER",
    [
        "-set META_COPY_MUT 1.0 -set META_STD_DEV 0.01 -set MUT_RATE_SOURCE 2"
    ]
)

combos.add_val(
    "COSTLY_COPY__COPY_OVER",
    [
        "-set COSTLY_HEAD_COPY 1 -set MAX_HEAD_COPY_COST 0.3",
        "-set COSTLY_HEAD_COPY 1 -set MAX_HEAD_COPY_COST 0.1",
        "-set COSTLY_HEAD_COPY 1 -set MAX_HEAD_COPY_COST 0.03",
        "-set COSTLY_HEAD_COPY 1 -set MAX_HEAD_COPY_COST 0.01",
        "-set COSTLY_HEAD_COPY 1 -set MAX_HEAD_COPY_COST 0.003",
        "-set COSTLY_HEAD_COPY 0 -set MAX_HEAD_COPY_COST 0.0"
    ]
)

combos.add_val(
    "COPY_MUT_PROB",
    [
        "0.0001",
        "0.001",
        "0.0316"
    ]
)

combos.add_val(
    "EVENTS_FILE__COPY_OVER",
    [
        "-set EVENT_FILE events_const-a.cfg -set ANALYZE_FILE analyze_env-a.cfg",
        "-set EVENT_FILE events_const-b.cfg -set ANALYZE_FILE analyze_env-b.cfg"
    ]
)



# Load in the base slurm file
with open(base_script_filename, 'r') as fp:
    base_sub_script = fp.read()

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
    parser.add_argument("--job_dir", type=str, default=None, help="Where to output these job files? If none, put in 'jobs' directory inside of the data_dir")
    parser.add_argument("--replicates", type=int, default=default_num_replicates, help="How many replicates should we run of each condition?")
    parser.add_argument("--account", type=str, default=default_account, help="Value to use for the slurm ACCOUNT")
    parser.add_argument("--seed_offset", type=int, default=default_seed_offset, help="Value to offset random number seeds by")

    # Load in command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    config_dir = args.config_dir
    job_dir = args.job_dir
    num_replicates = args.replicates
    hpc_account = args.account
    seed_offset = args.seed_offset

    # Get list of all combinations to run
    combo_list = combos.get_combos()
    # Calculate how many jobs we have, and what the last id will be
    num_jobs = num_replicates * len(combo_list)
    print(f'Generating {num_jobs} across {len(combo_list)} files!')
    print(f' - Data directory: {data_dir}')
    print(f' - Config directory: {config_dir}')
    print(f' - Replicates: {num_replicates}')
    print(f' - Account: {hpc_account}')
    print(f' - Seed offset: {seed_offset}')

    # If no job_dir provided, default to data_dir/jobs
    if job_dir == None:
        job_dir = os.path.join(data_dir, "jobs")

    # Create job file for each condition
    cur_job_id = 0
    cond_i = 0
    for condition_dict in combo_list:
        cur_seed = seed_offset + (cur_job_id * num_replicates)
        filename_prefix = f'RUN_C{cond_i}'
        file_str = base_sub_script
        file_str = file_str.replace("<<TIME_REQUEST>>", job_time_request)
        file_str = file_str.replace("<<ARRAY_ID_RANGE>>", f"1-{num_replicates}")
        file_str = file_str.replace("<<MEMORY_REQUEST>>", job_memory_request)
        file_str = file_str.replace("<<JOB_NAME>>", job_name)
        file_str = file_str.replace("<<CONFIG_DIR>>", config_dir)
        file_str = file_str.replace("<<EXEC>>", executable)
        file_str = file_str.replace("<<JOB_SEED_OFFSET>>", str(cur_seed))
        file_str = file_str.replace("<<ACCOUNT_NAME>>", hpc_account)

        # ===================================================
        # Configure the run
        # ===================================================
        file_str = file_str.replace("<<RUN_DIR>>", \
            os.path.join(data_dir, f'{filename_prefix}_'+'${SEED}'))

        # Format commandline arguments for the run
        run_param_info = {}

        run_param_info = {key:condition_dict[key] for key in condition_dict if not "__COPY_OVER" in key}
        run_param_info["RANDOM_SEED"] = '${SEED}'

        fields = list(run_param_info.keys())
        fields.sort()

        set_params = [f"-set {field} {run_param_info[field]}" for field in fields]
        copy_params = [condition_dict[key] for key in condition_dict if "__COPY_OVER" in key]

        run_params = " ".join(set_params + copy_params)

        # Add run commands to run the experiment
        run_commands = ''
        run_commands += f'RUN_PARAMS="{run_params}"\n'
        run_commands += 'echo "./${EXEC} ${RUN_PARAMS}" > cmd.log\n'
        run_commands += './${EXEC} ${RUN_PARAMS} > run.log\n'

        file_str = file_str.replace("<<RUN_COMMANDS>>", run_commands)

        # ===================================================
        # Configure analyze mode
        # ===================================================
        # Add commands to run avida analyze mode
        analysis_commands = ''
        analysis_commands += f'RUN_PARAMS="{run_params}"\n'
        analysis_commands += './${EXEC} ${RUN_PARAMS} -a\n'
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
