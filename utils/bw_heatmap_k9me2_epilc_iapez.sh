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

computeMatrix scale-regions -o ../results/deeptools/matrix_epilc_k9me2_IAPez.gz \
 -b 2500 -a 2500 \
 --regionBodyLength 5000 \
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_0h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_1h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_6h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_12h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_24h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_36h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_48h_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K9me2/EpiLC_H3K9me2_72h_RPGC.bigwig \
 -R ../data/bed/Navarro_2020_IAPEz_consensus_mm10.bed \
 -b 3000 -a 3000 \
 --sortRegions keep \
 --samplesLabel "K9me2 EpiLC 0h" "K9me2 EpiLC 1h" "K9me2 EpiLC 6h" "K9me2 EpiLC 12h" "K9me2 EpiLC 24h" "K9me2 EpiLC 36h" "K9me2 EpiLC 48h" "K9me2 EpiLC 72h" --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_epilc_k9me2_IAPez.gz \
 -out "../results/deeptools/EpiLC_H3.3K9me2_over_IAPez-Navarro2020.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --refPointLabel "IAPez" \
 --yMin 0 0 0 0 0 --yMax 2 2 2 2 2 2 2 2 \
 --zMin 0 0 0 0 0 --zMax 2 2 2 2 2 2 2 2 \
 -z "IAPez mm10 (Navarro, 2020)" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel ""