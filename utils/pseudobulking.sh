#!/bin/bash -l
#SBATCH -A naiss2023-22-84 #
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -M rackham
#SBATCH -J pseudobulking

module load bioinfo-tools
module load R_packages
module load samtools/1.17
module load deepTools

#working directory
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

# group sorted bam files
#Rscript pseudobulking.R

# merge and convert to bigwigs
cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_6h
#samtools merge EpiLC_6h_pseudobulk.bam *.bam
samtools sort EpiLC_6h_pseudobulk.bam -o EpiLC_6h_pseudobulk_sorted.bam
samtools index EpiLC_6h_pseudobulk_sorted.bam
bamCoverage --bam EpiLC_6h_pseudobulk_sorted.bam -o EpiLC_6h_pseudobulk_RPGC.bigwig \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 265278350

# create bedgraph format for SEACR peak caller
samtools sort -n EpiLC_6h_pseudobulk_sorted.bam -o EpiLC_6h_pseudobulk_name_sorted.bam

bedtools bamtobed -bedpe -i EpiLC_6h_pseudobulk_name_sorted.bam > EpiLC_6h_pseudobulk_sorted.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' EpiLC_6h_pseudobulk_sorted.bed > EpiLC_6h_pseudobulk_clean.bed
cut -f 1,2,6 EpiLC_6h_pseudobulk_clean.bed | sort -k1,1 -k2,2n -k3,3n > EpiLC_6h_pseudobulk_fragments.bed
bedtools genomecov -bg -i EpiLC_6h_pseudobulk_fragments.bed -g /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt > EpiLC_6h_pseudobulk.bedgraph

rm EpiLC_6h_pseudobulk_sorted.bed
rm EpiLC_6h_pseudobulk_clean.bed
rm EpiLC_6h_pseudobulk_fragments.bed

cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_12h
#samtools merge EpiLC_12h_pseudobulk.bam *.bam
samtools sort EpiLC_12h_pseudobulk.bam -o EpiLC_12h_pseudobulk_sorted.bam
samtools index EpiLC_12h_pseudobulk_sorted.bam
bamCoverage --bam EpiLC_12h_pseudobulk_sorted.bam -o EpiLC_12h_pseudobulk_RPGC.bigwig \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 265278350

# create bedgraph format for SEACR peak caller
samtools sort -n EpiLC_12h_pseudobulk_sorted.bam -o EpiLC_12h_pseudobulk_name_sorted.bam

bedtools bamtobed -bedpe -i EpiLC_12h_pseudobulk_name_sorted.bam > EpiLC_12h_pseudobulk_sorted.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' EpiLC_12h_pseudobulk_sorted.bed > EpiLC_12h_pseudobulk_clean.bed
cut -f 1,2,6 EpiLC_12h_pseudobulk_clean.bed | sort -k1,1 -k2,2n -k3,3n > EpiLC_12h_pseudobulk_fragments.bed
bedtools genomecov -bg -i EpiLC_12h_pseudobulk_fragments.bed -g /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt > EpiLC_12h_pseudobulk.bedgraph

rm EpiLC_12h_pseudobulk_sorted.bed
rm EpiLC_12h_pseudobulk_clean.bed
rm EpiLC_12h_pseudobulk_fragments.bed

cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_24h
#samtools merge EpiLC_24h_pseudobulk.bam *.bam
samtools sort EpiLC_24h_pseudobulk.bam -o EpiLC_24h_pseudobulk_sorted.bam
samtools index EpiLC_24h_pseudobulk_sorted.bam
bamCoverage --bam EpiLC_24h_pseudobulk_sorted.bam -o EpiLC_24h_pseudobulk_RPGC.bigwig \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 265278350

# create bedgraph format for SEACR peak caller
samtools sort -n EpiLC_24h_pseudobulk_sorted.bam -o EpiLC_24h_pseudobulk_name_sorted.bam

bedtools bamtobed -bedpe -i EpiLC_24h_pseudobulk_name_sorted.bam > EpiLC_24h_pseudobulk_sorted.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' EpiLC_24h_pseudobulk_sorted.bed > EpiLC_24h_pseudobulk_clean.bed
cut -f 1,2,6 EpiLC_24h_pseudobulk_clean.bed | sort -k1,1 -k2,2n -k3,3n > EpiLC_24h_pseudobulk_fragments.bed
bedtools genomecov -bg -i EpiLC_24h_pseudobulk_fragments.bed -g /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt > EpiLC_24h_pseudobulk.bedgraph

rm EpiLC_24h_pseudobulk_sorted.bed
rm EpiLC_24h_pseudobulk_clean.bed
rm EpiLC_24h_pseudobulk_fragments.bed

cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_48h
#samtools merge EpiLC_48h_pseudobulk.bam *.bam
samtools sort EpiLC_48h_pseudobulk.bam -o EpiLC_48h_pseudobulk_sorted.bam
samtools index EpiLC_48h_pseudobulk_sorted.bam
bamCoverage --bam EpiLC_48h_pseudobulk_sorted.bam -o EpiLC_48h_pseudobulk_RPGC.bigwig \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 265278350

# create bedgraph format for SEACR peak caller
samtools sort -n EpiLC_48h_pseudobulk_sorted.bam -o EpiLC_48h_pseudobulk_name_sorted.bam

bedtools bamtobed -bedpe -i EpiLC_48h_pseudobulk_name_sorted.bam > EpiLC_48h_pseudobulk_sorted.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' EpiLC_48h_pseudobulk_sorted.bed > EpiLC_48h_pseudobulk_clean.bed
cut -f 1,2,6 EpiLC_48h_pseudobulk_clean.bed | sort -k1,1 -k2,2n -k3,3n > EpiLC_48h_pseudobulk_fragments.bed
bedtools genomecov -bg -i EpiLC_48h_pseudobulk_fragments.bed -g /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/mm10.chrom.sizes.txt > EpiLC_48h_pseudobulk.bedgraph

rm EpiLC_48h_pseudobulk_sorted.bed
rm EpiLC_48h_pseudobulk_clean.bed
rm EpiLC_48h_pseudobulk_fragments.bed
