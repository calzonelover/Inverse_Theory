#!/bin/bash

#SBATCH -J Regularization # Job name 
#SBATCH --partition long
#SBATCH -o log_Regularization.out # Name of stdout output file (%j becomes %jobId)
#SBATCH -N 1 # Total number of nodes requested
#SBATCH -t 72:00:00 # Run time (hh:mm:ss) - 1.0 hours

module load anaconda/3.7

python main.py

