library("tidyverse")


bam_dir = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/mapped_mm10_20230510_EpiLC"
setwd(bam_dir)
bams = list.files(bam_dir,
                  full.names = FALSE, pattern = ".fixmkdup.nsort.bam")

EpiLC_6h = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_6h"
system(paste0("mkdir -p ", EpiLC_6h))
EpiLC_12h = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_12h"
system(paste0("mkdir -p ", EpiLC_12h))
EpiLC_24h = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_24h"
system(paste0("mkdir -p ", EpiLC_24h))
EpiLC_48h = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/20230510_EpiLC/EpiLC_48h"
system(paste0("mkdir -p ", EpiLC_48h))

for(bam in bams) {
  sample_id = as.numeric(strsplit(bam, split = "_")[[1]][3])
  if(sample_id >= 1 & sample_id <= 96) {
    system(paste0("cp ", bam, " ", EpiLC_6h))
  } else if(sample_id > 96 & sample_id <= 192) {
    system(paste0("cp ", bam, " ", EpiLC_12h))
  } else if(sample_id > 192 & sample_id <= 288) {
    system(paste0("cp ", bam, " ", EpiLC_24h))
  } else if(sample_id > 288 & sample_id <= 384) {
    system(paste0("cp ", bam, " ", EpiLC_48h))
  }
}

