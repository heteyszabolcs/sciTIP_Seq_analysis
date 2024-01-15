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
 # ../results/deeptools/matrix_ATAC_cm.gz \
 # -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bigwig/ATAC-H33WT.mm10.bw \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_ATAC-Seq/EpiLC/ATAC_EpiLC_0H_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_ATAC-Seq/EpiLC/ATAC_EpiLC_6H_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_ATAC-Seq/EpiLC/ATAC_EpiLC_12H_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_ATAC-Seq/EpiLC/ATAC_EpiLC_24H_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_ATAC-Seq/EpiLC/ATAC_EpiLC_48H_RPGC.bigwig \
 # -R ../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed \
 # -b 3000 -a 3000 \
 # --referencePoint center \
 # --samplesLabel "mESC" "EpiLC 0h" "EpiLC 6h" "EpiLC 12h" "EpiLC 24h" "EpiLC 48h" \
 # --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_ATAC_cm.gz \
 -out "../results/deeptools/mESC_cruz-molina_enh-clusters-ATAC.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 0 0 0 0 0 --yMax 200 200 200 200 200 200 \
 --zMin 0 0 0 0 0 0 --zMax 200 200 200 200 200 200 \
 -z "Cruz-Molina active mESC enhancers" \
 --yAxisLabel "" \
 --xAxisLabel "" \