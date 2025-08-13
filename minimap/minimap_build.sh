#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate minimap               # Activate your environment
#this will produce a warning, but can still run regardless

minimap2 -d GRCh38_trna.mmi hg38-tRNAs.fa
