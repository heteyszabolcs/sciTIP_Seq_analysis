#!/usr/bin/env Rscript
system("module load bioinfo-tools")
system("module load R_packages")

library("tidyverse")
library("glue")

# folders
original_dir = "/proj/snic2020-6-3/LINXUAN/sciTIP_sc_analysis/fastqs_single_cell_row_H/"
new_fastq_dir = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/fastq/deindexed_fastq/row_H/"

# copy into the working fastq directory
system(glue("mkdir {new_fastq_dir}"))
system(glue("cp {original_dir}* {new_fastq_dir}"))

# fastqs of one row
fastqs = list.files(new_fastq_dir)

for(i in fastqs) {
  parts = strsplit(i, "_")[[1]]
  newname = paste("Row", parts[2], parts[4], parts[5], parts[6], parts[7], sep = "_")
  system(glue("mv {new_fastq_dir}{i} {new_fastq_dir}{newname}"))
}


