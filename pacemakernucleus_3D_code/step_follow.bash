#!/bin/bash
#SBATCH -p long
#SBATCH -n 300
#SBATCH -t 25:00:00
#SBATCH --mem=200Gb
#SBATCH -J New_Follow_Step

module load neuron/develop
source activate neuron-develop

mpirun -n 300 python run_follow.py
