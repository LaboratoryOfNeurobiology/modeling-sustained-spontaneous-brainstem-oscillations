#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=200Gb
#SBATCH -J New_Follow_Step

module load neuron/develop
source activate neuron-develop

python gen_follow_step.py 174619 5 
