#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate python             # Activate your environment

bwa index combined_hg38-tRNAs.fa
