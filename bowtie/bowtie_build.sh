#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

# Load the required environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bowtie                 # Activate your environment

# Conversion command
bowtie2-build GCA_000001405.15_GRCh38_genomic.fna GRCh38

echo "indexing reference genome"

