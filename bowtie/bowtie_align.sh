#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

# ===== Required files and resources =====
# 1. FASTQ input files
#    - Location: set in $FASTQ_DIR
#    - Format: gzip-compressed FASTQ (.fastq.gz) files containing sequencing reads
#
# 2. Bowtie2 genome index files
#    - Location: set in $INDEX
#    - Must be built beforehand using `bowtie2-build`
#
# 3. Output directory
#    - Location: set in $OUTPUT_DIR
#    - Will contain SAM, BAM, and BAM index files

# Load the required environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bowtie 

# Define paths
FASTQ_DIR="/net/bmc-lab8/data/lab/calo/users/ddang/240725Cal/R24-5134/20240830_1459_2B_PAY10856_11fa3af5/fastq_fail"      
OUTPUT_DIR="/net/bmc-lab8/data/lab/calo/users/ddang/240725Cal/R24-5134/20240830_1459_2B_PAY10856_11fa3af5/bowtie_alignment"
INDEX="GRCh38_index"  

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Align, convert, sort, and index for each file
for fq in "$FASTQ_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz) 
    
    echo "Processing $base..."
    
    # Alignment (SAM output)
    bowtie2 -x "$INDEX" -U "$fq" -S "$OUTPUT_DIR/${base}.sam"
    
    # Convert SAM to BAM
    samtools view -Sb "$OUTPUT_DIR/${base}.sam" > "$OUTPUT_DIR/${base}.bam"
    
    # Sort BAM
    samtools sort "$OUTPUT_DIR/${base}.bam" -o "$OUTPUT_DIR/${base}_sorted.bam"
    
    # Index BAM
    samtools index "$OUTPUT_DIR/${base}_sorted.bam"
    
    # (Optional) Remove intermediate files to save space
    rm "$OUTPUT_DIR/${base}.sam" "$OUTPUT_DIR/${base}.bam"
    
    echo "$base processing complete!"
done

echo "All alignments and BAM indexing complete!"
