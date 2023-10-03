#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -M rackham
#SBATCH -J sciTIP_Seq_pipeline

module load bioinfo-tools
module load samtools/1.17
module load BEDTools/2.29.2
module load bowtie2
module load R_packages
module load ucsc-utilities

#working directory
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

Rscript all.R -i $1 -w $2
