if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "GenomicRanges",
               "wigglescout") 

order_by_signal = function(bw = "../data/ChIP-Seq/Yang_et_al_EpiLC_GSE117896/H3K27ac/EpiLC_H3K27ac_0h_RPGC.bigwig",
                           bed = "../data/bed/ESC_Enhancer_CruzMolina_only_enhancers.bed") {
  
  gr = fread(bed)
  gr = GRanges(
    seqnames = gr$V1,
    ranges = IRanges(
      start = gr$V2,
      end = gr$V3
    )
  )
  
  
  df = bw_loci(bwfiles = bw, loci = gr)
  df = as.data.frame(df)
  df = df[order(df[,ncol(df)], decreasing = TRUE),]
  
  df = df %>% dplyr::select(seqnames, start, end)
  write_tsv(df, bed, 
            col_names = FALSE)
  print("ready")
  return(head(df))
}

bed = order_by_signal(bw = "../results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/EpiLC_6h_pseudobulk_RPGC.bigwig",
                bed = "../data/bed/UCSC_RepeatMasker_IAPEz_mm10.bed")
