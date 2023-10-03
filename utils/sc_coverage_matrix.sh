#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -M rackham
#SBATCH -J coverage

module load bioinfo-tools
module load BEDTools/2.29.2
module load R_packages
module load ucsc-utilities

#working directory
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

Rscript sc_coverage_matrix.R
