#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools

cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

# computeMatrix reference-point -o \
 # ../results/deeptools/matrix_SMC1_GRHL2.gz \
 # -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Chen_et_al_EpiLC_GSE93147/fastq/EpiLC_GRHL2_rep1_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Chen_et_al_EpiLC_GSE93147/fastq/EpiLC_GRHL2_rep2_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Chen_et_al_EpiLC_GSE93147/fastq/EpiLC_SMC1_rep1_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Chen_et_al_EpiLC_GSE93147/fastq/EpiLC_SMC1_rep2_RPGC.bigwig \
 # -R ../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster1.bed \
 # ../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster2.bed \
 # -b 3000 -a 3000 \
 # --referencePoint center \
 # --sortRegions keep \
 # --samplesLabel "GRHL2, 1" "GRHL2, 2", "SMC1, 1", "SMC1, 2" \
 # --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_SMC1_GRHL2.gz \
 -out "../results/deeptools/Yang_EpiLC_K27ac_kmeans-clusters-SMC1_GRHL2.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 0 0 0 0 0 --yMax 20 20 20 20 \
 --zMin 0 0 0 0 0 0 --zMax 20 20 20 20 \
 -z "cluster 1" "cluster 2" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel "" \