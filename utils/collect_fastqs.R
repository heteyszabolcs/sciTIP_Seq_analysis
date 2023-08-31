# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
})

clusters = c("cluster_0", "cluster_1", "cluster_2", "cluster_3", "cluster_4","cluster_5")
cluster_ids = list.files("../results/Seurat/", pattern = "H3.2_res0.5_cluster*", full)
fastqs = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/data/fastq/deindexed_fastq_20230316_H3.2/"
output = "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/Seurat/cluster_fastqs/"

collect_fastq = function(cluster_ids, cluster) {
  cluster_ids = fread(cluster_ids, header = FALSE)
  
  for (id in cluster_ids$V1) {
    id = str_split(id, pattern = "1_001")[[1]][1]
    system(paste0("mkdir ", output, cluster))
    system(paste0("cp ", fastqs, id, "*", " ", output, cluster))
    
  }

  print(paste0(cluster, " is done!"))
}

for(i in seq(length(clusters))) {
  collect_fastq(cluster_ids = cluster_ids[i], cluster = clusters[i])
}




