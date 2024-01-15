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



# /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_0h_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_1h_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_6h_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_12h_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_24h_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_36h_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_48h_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K4me1/EpiLC_H3K4me1_72h_RPGC.bigwig \

# computeMatrix reference-point -o \
 # ../results/deeptools/matrix_cm_summary.gz \
 # -S  /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/Oct4/Oct4_L001_R1_001_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/Oct4/Oct4_L002_R1_001_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/p300/p300_L001_R1_001_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/p300/p300_L002_R1_001_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/H3K27ac/H3K27ac_L001_R1_001_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/H3K27ac/H3K27ac_L002_R1_001_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Nocente_et_al_mESC_GSE210780/Flag-Brg1/MNase_ChIPSeq_Flag-Brg1_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/ChIP-Seq/Navarro_et_al_mESC_GSE149080/mESC_H33_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_6h/EpiLC_6h_pseudobulk_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_12h/EpiLC_12h_pseudobulk_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_24h/EpiLC_24h_pseudobulk_RPGC.bigwig \
 # /proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_48h/EpiLC_48h_pseudobulk_RPGC.bigwig \
 # -R ../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed \
 # -b 3000 -a 3000 \
 # --referencePoint center \
 # --samplesLabel "EpiLC Oct4, 1" "EpiLC Oct4, 2" "EpiLC p300, 1" "EpiLC p300, 2" "EpiLC K27ac, 1" "EpiLC K27ac, 2" "mESC Brg1" "mESC H3.3" "EpiLC 6h" "EpiLC 12h" "EpiLC 24h" "EpiLC 48h" \
 # --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_cm_summary.gz \
 -out "../results/deeptools/cruz-molina_enhancers.pdf" \
 --refPointLabel "origin" \
 --heatmapHeight 14 \
 --colorMap "Reds" \
 --refPointLabel "enh" \
 --yMin 0 0 0 0 0 0 0 0 --yMax 120 120 120 120 250 250 10 20 4 4 4 4 \
 --zMin 0 0 0 0 0 0 0 0 --zMax 120 120 120 120 250 250 10 20 4 4 4 4 \
 -z "Cruz-Molina active mESC enhancers" \
 --yAxisLabel "" \
 --xAxisLabel "" \