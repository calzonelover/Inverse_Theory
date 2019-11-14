#!/bin/bash

#SBATCH -J InvLM # Job name 
#SBATCH -o log_InvLM.out
#SBATCH --partition long
#SBATCH -n 1 # Total number tasks per node
#SBATCH --cpus-per-task=4
#SBATCH -t 70:00:00 # Run time (hh:mm:ss) - 1.0 hours

module load prun/1.3
module load anaconda/3.7

python -W ignore main.py

