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
 ../results/deeptools/matrix_cm.gz \
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_EpiLC_GSE56138/H3K27ac/H3K27ac_L001_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_EpiLC_GSE56138/p300/p300_L001_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/EpiLC_6h_pseudobulk_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/EpiLC_12h_pseudobulk_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/EpiLC_24h_pseudobulk_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/EpiLC_48h_pseudobulk_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/mESC_GSE59189/H3.3_wt/mESC_H3.3_L001_R1_001_RPGC.bigwig \
 -R ../data/bed/ESC_Enhancer_CruzMolina_only_enhancers-ordered.bed \
 -b 3000 -a 3000 \
 --referencePoint center \
 --sortRegions keep \
 --samplesLabel "EpiLC H3K27ac" "EpiLC p300" "6h" "12h" "24h" "48h" "mESC H3.3" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_cm.gz \
 -out "../results/deeptools/CruzMolina_enh.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 0 0 0 0 0 0 --yMax 200 100 10 10 10 10 10 \
 --zMin 0 0 0 0 0 0 0 --zMax 200 100 10 10 10 10 10 \
 -z "Cruz-Molina mESC enhancers" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel "" \