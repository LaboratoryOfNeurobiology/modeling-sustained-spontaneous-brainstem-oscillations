#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mem=10Gb
#SBATCH -J Split_Iteration

module load neuron/develop
source activate neuron-develop

python Split_iteration.py 850980 5 8
