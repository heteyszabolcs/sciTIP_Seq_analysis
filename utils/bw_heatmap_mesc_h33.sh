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
 ../results/deeptools/matrix_mesc_h33_epilc_cm_enh.gz \
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Navarro_et_al_mESC_GSE149080/mESC_H33_RPGC.bigwig \
 -R ../data/bed/ESC_Enhancer_CruzMolina_only_enhancers.bed \
 ../data/bed/GSE56138_EpiLC_only_enhancers-Oct4_p300_K27ac.bed \
 -b 3000 -a 3000 \
 --referencePoint center \
 --sortRegions keep \
 --samplesLabel "mESC H3.3" --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_mesc_h33_epilc_cm_enh.gz \
 -out "../results/deeptools/mESC_H3.3_over_epilc_cm_enh.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 --yMax 20 \
 --zMin 0 --zMax 20 \
 -z "EpiLC enhancers" "Cruz-Molina mESC enhancers" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel ""