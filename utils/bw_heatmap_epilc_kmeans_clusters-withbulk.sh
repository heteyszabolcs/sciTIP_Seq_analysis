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
 -S /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/Oct4/Oct4_L001_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/Oct4/Oct4_L002_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/p300/p300_L001_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/p300/p300_L002_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/H3K27ac/H3K27ac_L001_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/H3K27ac/H3K27ac_L002_R1_001_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Nocente_et_al_mESC_GSE210780/Flag-Brg1/MNase_ChIPSeq_Flag-Brg1_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Navarro_et_al_mESC_GSE149080/mESC_H33_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_0H_1_S26_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_6H_1_S28_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_12H_1_S30_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_24H_1_S32_RPGC.bigwig \
 /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/bulk_TIP-Seq/Bulk_TIP_48H_1_S34_RPGC.bigwig \
 -R ../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster1.bed \
 ../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster2.bed \
 -b 3000 -a 3000 \
 --referencePoint center \
 --sortRegions keep \
 --samplesLabel "EpiLC Oct4, 1" "EpiLC Oct4, 2" "EpiLC p300, 1" "EpiLC p300, 2" "EpiLC K27ac, 1" "EpiLC K27ac, 2" "mESC Brg1" "mESC" "EpiLC 0h" "EpiLC 6h" "EpiLC 12h" "EpiLC 24h" "EpiLC 48h" \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_cm.gz \
 -out "../results/deeptools/Yang_EpiLC_K27ac_kmeans-clusters-withbulk_tip.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 0 0 0 0 0 --yMax 100 100 100 100 250 250 10 40 40 40 40 40 40 \
 --zMin 0 0 0 0 0 0 --zMax 100 100 100 100 250 250 10 40 40 40 40 40 40 \
 -z "cluster 1" "cluster 2" \
 --sortRegions keep \
 --yAxisLabel "" \
 --xAxisLabel "" \