# packages
suppressPackageStartupMessages({
  library("rtracklayer")
  library("glue")
  library("data.table")
  library("tidyverse")
})

# liftover function
mm9_to_mm10_liftover = function(bw = "../data/bw/GSM5625015_NPC_G4_CnT_Rep2_Batch2.bw",
                                chrom_sizes = "../data/mm10.chrom.sizes.txt",
                                chain = "../data/bw/mm9ToMm10.over.chain",
                                export = "../data/bw/GSM5625015_NPC_G4_CnT_Rep2_Batch2_mm10.bw") {
  mm10_chrom = fread(chrom_sizes, header = FALSE)
  bigwig = import(bw)
  chain = import.chain(chain)
  lo = liftOver(bigwig, chain)
  lo = unlist(lo)
  genome(lo) = "mm10"
  sl = tibble(chrom = names(seqlengths(lo)))
  sl = sl %>% left_join(., mm10_chrom, by = c("chrom" = "V1"))
  seqlengths(lo) = sl$V2
  hits = findOverlaps(lo, drop.self = TRUE)
  lo = lo[-queryHits(hits)]
  export(object = lo, con = export)

}

mm9_to_mm10_liftover(bw = "../data/bigwig/Yang_2019/GSM3314701_H3K27ac-0h.bw",
                     chrom_sizes = "../data/mm10.chrom.sizes.txt",
                     chain = "../data/bigwig/mm9ToMm10.over.chain",
                     export = "../data/bigwig/Yang_2019/GSM3314701_H3K27ac-0h_mm10.bw")



bws = list.files("../data/bigwig/Yang_2019/", full.names = TRUE, pattern = "*h.bw")

output_names = unname(sapply(bws, function(x) {
  strsplit(x, split = ".bw")[[1]][1]
}))
output_names = unname(sapply(output_names, function(x) {
  paste0(x, "_mm10.bw")
}))

for(i in seq(1, length(bws))) {
  mm9_to_mm10_liftover(
    bw = bws[i],
    chrom_sizes = "../data/mm10.chrom.sizes.txt",
    chain = "../data/bigwig/mm9ToMm10.over.chain",
    export = output_names[i]
  )
  
}



