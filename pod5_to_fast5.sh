#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

# Load the required environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pod5_env                 # Activate your environment

#this will produce a warning, but can still run regardless

# Conversion command
pod5 convert fast5 ./fast5_fail --output pod5_fail

echo "POD5 to FAST5 conversion completed."

