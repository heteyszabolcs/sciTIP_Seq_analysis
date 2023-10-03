#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -M rackham
#SBATCH -J SNPsplit

module load bioinfo-tools
module load R_packages

#working directory
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

#i: path to fastq files
#w: working directory, where the script will create the folders
#o: output for alignment files

Rscript snp_split.R -i $1 -w $2 -o $3
