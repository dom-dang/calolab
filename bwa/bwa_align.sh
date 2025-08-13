#!/bin/bash
#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                    
#SBATCH --mail-type=END          
#SBATCH --mail-user=ddang@mit.edu       
#############################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate BWA             # Activate your environment
#this will produce a warning, but can still run regardless

# Define input files
REFERENCE="GRCh38"
READS="r24-5135_dorado_v2.fastq"
OUTPUT_PREFIX="R24-5134_bwa_v2"

# Step 1: Run BWA MEM for Nanopore reads
echo "Running BWA MEM..."
# this line takes too much time 
#bwa mem -W13 -k6 -T20 -xont2d $REFERENCE $READS > ${OUTPUT_PREFIX}.sam

#v1 BWA: 
#bwa mem -W13 -xont2d $REFERENCE $READS > ${OUTPUT_PREFIX}.sam

#v2 BWA:
bwa mem -W13 -k6 -xont2d $REFERENCE $READS > ${OUTPUT_PREFIX}.sam

# -W13: Sets the band width for banded alignment to 13. This can affect alignment speed and accuracy (smaller = faster but potentially less accurate).
# -k6: Sets the minimum seed length to 6. Seeds are short exact matches that BWA uses to start the alignment process. Shorter seeds can make it more sensitive, but potentially more prone to false positives.
#-xont2d: Tells BWA to use preset options for Oxford Nanopore 2D reads, which are long and error-prone. This adjusts internal scoring parameters for better performance with that data type.
#-T20: Sets the minimum score to output an alignment to 20. This filters out weak alignments.


# Step 2: Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view -Sb ${OUTPUT_PREFIX}.sam > ${OUTPUT_PREFIX}.bam

# Step 3: Sort BAM file
echo "Sorting BAM file..."
samtools sort ${OUTPUT_PREFIX}.bam -o ${OUTPUT_PREFIX}_sorted.bam

# Step 4: Index sorted BAM file
echo "Indexing BAM file..."
samtools index ${OUTPUT_PREFIX}_sorted.bam


echo "Alignment completed! Output: ${OUTPUT_PREFIX}_sorted.bam"
