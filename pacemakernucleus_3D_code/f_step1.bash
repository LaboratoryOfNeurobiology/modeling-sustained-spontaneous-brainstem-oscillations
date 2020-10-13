#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH --mem=20Gb
#SBATCH -t 01:00:00
#SBATCH -J Initialize_grid

module load neuron/develop
source activate neuron-develop

python Initialize_Follow_Grid.py 898944 787454

