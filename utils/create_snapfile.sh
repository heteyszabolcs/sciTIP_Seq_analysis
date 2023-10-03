#!/bin/bash -l 
#SBATCH -A naiss2023-22-84 
#SBATCH -p core 
#SBATCH -n 4 
#SBATCH -t 12:00:00 
#SBATCH -M rackham
#SBATCH -J snaptools

# activate virtenv
source /home/szabolcs/.pyenv/versions/3.8.1/envs/snaptools/bin/activate
pyv="$(python -V)"
echo "$pyv"

# change to workdir
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/scATAC-Seq

# run snaptools snap-pre
#snaptools snap-pre \
#	--input-file=GSM6070562_mESC_scATAC_Seq_fragments.bed.gz \
#	--output-snap=GSM6070562_mESC_scATAC_Seq.snap \
#	--genome-name=mm10 \
#	--genome-size=mm10.chrom.sizes.txt \
#	--min-mapq=30 \
#	--min-flen=50 \
#	--max-flen=1000 \
#	--keep-chrm=TRUE \
#	--keep-single=FALSE \
#	--keep-secondary=False \
#	--overwrite=True \
#	--max-num=20000 \
#	--min-cov=500 \
#	--verbose=True

# create cell-by-bin matrices
snaptools snap-add-bmat  \
	--snap-file=GSM6070562_mESC_scATAC_Seq.snap  \
	--bin-size-list 5000 10000 15000 20000  \
	--verbose=True

