#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

#   This SLURM batch script merges multiple POD5 files in the current 
#   working directory into a single merged POD5 file using the 
#   `pod5 merge` command-line tool

# Load the required environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pod5_env                 # Activate your environment

# Conversion command
pod5 merge *.pod5 -o merged.pod5

echo "POD5 merge completed."

