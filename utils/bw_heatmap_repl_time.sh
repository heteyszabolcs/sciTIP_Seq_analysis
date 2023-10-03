#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools
module load R_packages

cd /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/utils

computeMatrix reference-point -o \
 ../results/deeptools/matrix_repl_time_hm.mat.gz \
 -S ../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster0/cluster0_RPGC.bigwig \
 ../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster1/cluster1_RPGC.bigwig \
 ../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster2/cluster2_RPGC.bigwig \
 ../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster3/cluster3_RPGC.bigwig \
 ../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster4/cluster4_RPGC.bigwig \
 ../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster5/cluster5_RPGC.bigwig \
 -R ../data/bed/RepliTime_early_mm10.bed \
 ../data/bed/RepliTime_mid_mm10.bed \
 ../data/bed/RepliTime_late_mm10.bed \
 -b 3000 -a 3000 \
 --referencePoint center \
 --samplesLabel "cluster 0" "cluster 1" "cluster 2" "cluster 3" "cluster 4" "cluster 5" --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_repl_time_hm.mat.gz \
 -out "../results/deeptools/H3.2_repl_time_hm.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --yMin 0 \
  -z "early" "mid" "late" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 # --yMax 100 \
 # --zMax 100 \

 plotHeatmap -m ../results/deeptools/matrix_repl_time_hm.mat.gz \
 -out "../results/deeptools/H3.2_repl_time_hm.png" \
 --refPointLabel "H3.2" \
 --heatmapHeight 14 \
 --colorMap "origin" \
 --yMin 0 \
  -z "early" "mid" "late" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 # --yMax 100 \
 # --zMax 100 \