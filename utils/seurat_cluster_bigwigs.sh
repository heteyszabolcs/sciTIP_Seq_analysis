#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -M rackham
#SBATCH -J pseudobulking

module load bioinfo-tools
module load samtools/1.17
module load deepTools

# mkdir 
mkdir -p /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.2_res0.5/$1
dir="/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.2_res0.5/$1"

# go to bam directory
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230316_H3.2/mapped_mm10

# copy all bam files 
while read p; do 
    cp $p*fixmate.bam $dir
done < $2

cd $dir

# samtools
samtools merge $1.bam *.bam
samtools sort $1.bam -o $1.sorted.bam
samtools index $1.sorted.bam
bamCoverage --bam $1.sorted.bam -o ${1}_RPGC.bigwig \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 265278350

rm $dir/*fixmate.bam
