#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=200Gb
#SBATCH -J New_Iteration

module load neuron/develop
source activate neuron-develop

python generate_new_it.py 898944 5
