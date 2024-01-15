#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -M rackham
#SBATCH -J bowtie2

# modules
module load bioinfo-tools
module load samtools/1.17
module load bowtie2
module load deepTools

# go to dir of fastqs
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/

# add bowtie2 index and fastq files
INDEX=/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10_bowtie_index/mm10
DATA_DIR=/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/fastq/FASTQ
FASTQ1="Bulk_TIP_12H_2_S31_L001_R1_001.fastq.gz"
FASTQ2="Bulk_TIP_12H_2_S31_L001_R2_001.fastq.gz"

# workflow
# 1. alignment 2. bam converting 3. bam coverage
base=$(basename $FASTQ1 _L001_R1_001.fastq.gz)
echo $base
echo $base >> bowtie_log.txt
bowtie2 -x $INDEX -p 8 --very-sensitive-local -1 $DATA_DIR/$FASTQ1 -2 $DATA_DIR/$FASTQ2 2>>bowtie_log.txt | samtools view -S -b -q 10 -o $base.bam - 
echo "sorting $base.bam"
samtools sort $base.bam -T $base -o $base.bam
echo "indexing"
samtools index $base.bam
echo "create RPGC bigwig"
bamCoverage --bam $base.bam -o ${base}_RPGC.bigwig \
  --binSize 10 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 2652783500

