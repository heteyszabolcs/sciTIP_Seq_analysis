#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -M rackham
#SBATCH -J remove_T7_dups

module load bioinfo-tools
module load samtools/1.17
module load R_packages

#working directory
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

Rscript removeT7dups.R
