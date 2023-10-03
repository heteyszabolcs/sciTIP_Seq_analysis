#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -M rackham
#SBATCH -J allele_spec

module load bioinfo-tools
module load samtools/1.17
module load bowtie2
module load R_packages
module load ucsc-utilities

#working directory
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

#i: path to fastq files
#w: working directory, where the script will create the folders
#o: output for alignment files
#b: path to bowtie2 index

Rscript bowtie2_allele_spec.R -i $1 -w $2 -o $3 -b $4


