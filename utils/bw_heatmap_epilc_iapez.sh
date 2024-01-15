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

computeMatrix scale-regions -o ../results/deeptools/matrix_epilc_IAPez.gz \
 -b 2500 -a 2500 \
 --regionBodyLength 5000 \
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Navarro_et_al_mESC_GSE149080/mESC_H33_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_0H_1_S26_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_6H_1_S28_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_12H_1_S30_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_24H_1_S32_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_48H_1_S34_RPGC.bigwig \
 -R ../data/bed/UCSC_RepeatMasker_IAPEz_mm10.bed \
 -b 3000 -a 3000 \
 --sortRegions keep \
 --samplesLabel "mESC" "EpiLC 0h" "EpiLC 6h" "EpiLC 12h" "EpiLC 24h" "EpiLC 48h" --skipZeros --missingDataAsZero
 
 plotHeatmap -m ../results/deeptools/matrix_epilc_IAPez.gz \
 -out "../results/deeptools/bulk_EpiLC_H3.3_over_IAPez-UCSC_IAPEz.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --refPointLabel "IAPez" \
 --yMin 0 0 0 0 0 --yMax 10 2 2 2 2 2 \
 --zMin 0 0 0 0 0 --zMax 10 2 2 2 2 2 \
 -z "IAPez (UCSC RepeatMasker)" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel ""