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

computeMatrix scale-regions -o ../results/deeptools/matrix_epilc_k27ac_IAPez.gz \
 -b 2500 -a 2500 \
 --regionBodyLength 5000 \
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_0h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_1h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_6h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_12h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_24h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_36h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_48h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_72h_RPGC.bigwig \
 -R ../data/bed/Navarro_2020_IAPEz_consensus_mm10.bed \
 -b 3000 -a 3000 \
 --sortRegions keep \
 --samplesLabel "K27ac EpiLC 0h" "K27ac EpiLC 1h" "K27ac EpiLC 6h" "K27ac EpiLC 12h" "K27ac EpiLC 24h" "K27ac EpiLC 36h" "K27ac EpiLC 48h" "K27ac EpiLC 72h" --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_epilc_k27ac_IAPez.gz \
 -out "../results/deeptools/EpiLC_H3.3K27ac_over_IAPez-Navarro2020.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --refPointLabel "IAPez" \
 --yMin 0 0 0 0 0 --yMax 10 10 10 10 10 10 10 10 \
 --zMin 0 0 0 0 0 --zMax 10 10 10 10 10 10 10 10 \
 -z "IAPez mm10 (Navarro, 2020)" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel ""