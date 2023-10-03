#!/bin/bash -l

#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -M rackham
#SBATCH -J cellranger

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

# modules
module load bioinfo-tools
module load cellranger-ATAC/2.0.0

# call cellranger-atac count
cellranger-atac count --id=GSM6070562_ --fastqs=/crex/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/scATAC-Seq/fastq --reference=/sw/data/Chromium/cellranger-ATAC-data/2.0.0/rackham/refdata-cellranger-arc-mm10-2020-A-2.0.0
