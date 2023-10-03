#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -M rackham
#SBATCH -J snaptools

# activate snaptools pyenv
source /home/szabolcs/.pyenv/versions/3.8.1/envs/snaptools/bin/activate

# work dir
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

# with bowtie2
#snaptools align-paired-end --input-reference /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10_bowtie_index/mm10 \
#	--input-fastq1 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_fastqs/cluster_5/cluster_5_R1.fastq.gz \
#	--input-fastq2 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_fastqs/cluster_5/cluster_5_R2.fastq.gz \
#	--output-bam /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_fastqs/cluster_5/cluster_5.bam \
#	--aligner=bowtie2 \
#	--path-to-aligner=/home/szabolcs/aligners/bowtie2-2.4.2 \
#	--read-fastq-command=zcat \
#	--min-cov=0 \
#	--num-threads=5 \
#	--if-sort=True \
#	--tmp-folder=./ \
#	--overwrite=TRUE

# create snap
#snaptools snap-pre \
#        --input-file=/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster5/cluster5.name_sorted.bam \
#        --output-snap=/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster5/cluster5.snap \
#        --genome-name=mm10 \
#        --genome-size=/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt \
#        --min-mapq=30 \
#        --min-flen=0 \
#        --max-flen=1000 \
#        --keep-chrm=TRUE \
#        --keep-single=TRUE \
#        --keep-secondary=False \
#        --overwrite=True \
#        --max-num=1000000 \
#        --min-cov=100 \
#        --verbose=True \

# bin snap file
snaptools snap-add-bmat  \
	--snap-file=/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster0/cluster0.snap  \
	--bin-size-list 80000 100000  \
	--verbose=True
