#!/bin/bash
#SBATCH -p long
#SBATCH -n 128
#SBATCH -t 24:00:00
#SBATCH -J SimulateSplits
#SBATCH --array=0-7%8

module load neuron/develop
source activate neuron-develop

mpirun -n 128 python sim_network_split.py $SLURM_ARRAY_TASK_ID 5 $SLURM_ARRAY_JOB_ID
