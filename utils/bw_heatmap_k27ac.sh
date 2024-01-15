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
 # ../results/deeptools/matrix_mESC_K27ac.gz \
 # -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Creyghton_et_al_mESC_GSE24165/Creyghton_mESC_K27ac_RPGC.bigwig \
 # -R ../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster1.bed \
 # ../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster2.bed \
 # -b 3000 -a 3000 \
 # --referencePoint center \
 # --sortRegions keep \
 # --samplesLabel "mESC K27ac Creyghton et al" \
 # --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_mESC_K27ac.gz \
 -out "../results/deeptools/mESC_K27ac_over_epilc_enh-kmeans.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 --yMax 100 \
 --zMin 0 --zMax 100 \
 --sortRegions keep \
 -z "cluster_1" "cluster_2" \
 --yAxisLabel "" \
 --xAxisLabel "" \