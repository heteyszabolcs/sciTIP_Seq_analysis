#!/bin/bash -l

#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -M rackham
#SBATCH -J cellranger

# modules
module load bioinfo-tools
module load cellranger/6.1.2

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/fastq/deindexed_fastq/row_A/

# call cellranger count
cellranger count --fastqs=/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/fastq/deindexed_fastq/row_A/ \
--id=Row_1A \
--sample=Row_1A \
--transcriptome=/sw/data/Chromium/cellranger-data/2020-A/refdata-gex-mm10-2020-A \
--chemistry ARC-v1
