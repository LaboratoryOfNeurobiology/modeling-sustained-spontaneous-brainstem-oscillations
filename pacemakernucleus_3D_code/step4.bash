#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH --mem=15Gb
#SBATCH -t 01:00:00
#SBATCH -J Aggregate_Results

module load neuron/develop
source activate neuron-develop

python aggregate_results.py 098595 5 8
