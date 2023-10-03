#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -M rackham
#SBATCH -J pseudotime_analysis

module load bioinfo-tools
module load R_packages

# work dir
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

Rscript pseudotime_analysis.R
