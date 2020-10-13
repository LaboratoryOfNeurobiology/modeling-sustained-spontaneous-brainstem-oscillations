#!/bin/bash
#SBATCH -p short
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J Initialize_grid

module load neuron/develop
source activate neuron-develop

python Initialize_grid.py
