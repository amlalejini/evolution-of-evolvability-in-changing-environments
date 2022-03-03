#!/usr/bin/env bash

REPLICATES=10
EXP_SLUG=2022-03-costly-copy
ACCOUNT=devolab
SEED_OFFSET=1000

SCRATCH_EXP_DIR=/mnt/scratch/lalejini/data/evo-evolvability
HOME_EXP_DIR=/mnt/home/lalejini/devo_ws/evolution-of-evolvability-in-changing-environments/experiments

DATA_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}
JOB_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}/jobs
CONFIG_DIR=${HOME_EXP_DIR}/${EXP_SLUG}/hpc/config

python3 gen-sub.py --data_dir ${DATA_DIR} --config_dir ${CONFIG_DIR} --replicates ${REPLICATES} --job_dir ${JOB_DIR} --account ${ACCOUNT} --seed_offset ${SEED_OFFSET}