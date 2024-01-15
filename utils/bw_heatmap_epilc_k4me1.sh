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
 ../results/deeptools/matrix_epilc_K4me1_epilc_enh.gz \
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_0h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_1h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_6h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_12h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_24h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_36h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_48h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_72h_RPGC.bigwig \
 -R ../data/bed/GSE56138_EpiLC_only_enhancers-Oct4_p300_K27ac.bed \
 -b 3000 -a 3000 \
 --referencePoint center \
 --sortRegions keep \
 --samplesLabel "EpiLC K4me1 0h" "EpiLC K4me1 1h" "EpiLC K4me1 6h" "EpiLC K4me1 12h" "EpiLC K4me1 24h" "EpiLC K4me1 36h" "EpiLC K4me1 48h" "EpiLC K4me1 72h" --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_epilc_K4me1_epilc_enh.gz \
 -out "../results/deeptools/Yang_EpiLC_K4me1_over_epilc_enh.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 0 0 0 0 0 0 0 --yMax 40 40 40 40 40 40 40 40 \
 --zMin 0 0 0 0 0 0 0 0 --zMax 40 40 40 40 40 40 40 40 \
 -z "EpiLC enhancers" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel ""