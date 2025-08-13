#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

# Load the required environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bowtie 

# Define paths
SAM_DIR="/net/bmc-lab8/data/lab/calo/users/ddang/240725Cal/R24-5134/20240830_1459_2B_PAY10856_11fa3af5/bowtie_alignment"      
BAM_DIR="/net/bmc-lab8/data/lab/calo/users/ddang/240725Cal/R24-5134/20240830_1459_2B_PAY10856_11fa3af5/bam_alignment"

for sam in "$SAM_DIR"/*.sam; do
    base=$(basename "$sam" .sam)
    samtools view -bS "$sam" | samtools sort -o "$BAM_DIR/${base}_sorted.bam"
done

echo "BAM Complete!"
