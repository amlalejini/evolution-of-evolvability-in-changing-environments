#!/usr/bin/env bash

REPLICATES=10
EXP_SLUG=2022-03-07-variable-mut
ACCOUNT=devolab
SEED_OFFSET=2000
JOB_TIME=00:05:00

SCRATCH_EXP_DIR=/mnt/scratch/lalejini/data/evo-evolvability
REPO_DIR=/mnt/home/lalejini/devo_ws/evolution-of-evolvability-in-changing-environments
HOME_EXP_DIR=${HOME_EXP_DIR}/experiments

DATA_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}
JOB_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}/jobs
CONFIG_DIR=${HOME_EXP_DIR}/${EXP_SLUG}/hpc/config

python3 gen-sub.py --time_request ${JOB_TIME} --data_dir ${DATA_DIR} --config_dir ${CONFIG_DIR} --repo_dir ${REPO_DIR} --replicates ${REPLICATES} --job_dir ${JOB_DIR} --account ${ACCOUNT} --seed_offset ${SEED_OFFSET}