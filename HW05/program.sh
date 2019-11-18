#!/bin/bash

#SBATCH -J InvCG # Job name 
#SBATCH -o log_InvCG.out
#SBATCH -N 1 # nodes
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --nodelist=compute-c3
#SBATCH --mem=60G
#SBATCH -t 20:00:00 # Run time (hh:mm:ss) - 1.0 hours

module load prun/1.3 ohpc gnu8/8.3.0 autotools
module load anaconda/3.7

export QT_QPA_PLATFORM='offscreen'

python -W ignore main.py

