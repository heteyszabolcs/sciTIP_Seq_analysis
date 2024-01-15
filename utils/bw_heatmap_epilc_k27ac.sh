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

computeMatrix reference-point -o \
 ../results/deeptools/matrix_epilc_K27ac_epilc_enh2.gz \
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_0h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_1h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_6h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_12h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_24h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_36h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_48h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_72h_RPGC.bigwig \
 -R ../data/bed/GSE56138_EpiLC_only_enhancers-p300_K27ac.bed \
 -b 3000 -a 3000 \
 --referencePoint center \
 --sortRegions keep \
 --samplesLabel "EpiLC K27ac 0h" "EpiLC K27ac 1h" "EpiLC K27ac 6h" "EpiLC K27ac 12h" "EpiLC K27ac 24h" "EpiLC K27ac 36h" "EpiLC K27ac 48h" "EpiLC K27ac 72h" --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_epilc_K27ac_epilc_enh2.gz \
 -out "../results/deeptools/Yang_EpiLC_K27ac_over_epilc_enh2-kmeans.pdf" \
 --outFileSortedRegions "../results/deeptools/Yang_EpiLC_K27ac_over_epilc_enh2-kmeans.tab" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 0 0 0 0 0 0 0 --yMax 50 50 50 50 50 50 50 50 \
 --zMin 0 0 0 0 0 0 0 0 --zMax 50 50 50 50 50 50 50 50 \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --kmeans 2