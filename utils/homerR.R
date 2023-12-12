# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
})

# Zhifen et al.
zhifen = fread("../data/Zhifen_et_al-pluripotency_gene_families_Fig1D.txt", header = FALSE)


# HOMER results
epilc_k27ac_cl1 = fread("../results/HOMER/Yang_k27ac-cluster1/knownResults.txt")
epilc_k27ac_cl1 = epilc_k27ac_cl1 %>% separate(`Motif Name`, into = "gene", sep = "\\(")
epilc_k27ac_cl1_filt = epilc_k27ac_cl1 %>% dplyr::filter(gene %in% zhifen$V1)

# HOMER results
epilc_k27ac_cl2 = fread("../results/HOMER/Yang_k27ac-cluster2/knownResults.txt")
epilc_k27ac_cl2 = epilc_k27ac_cl2 %>% separate(`Motif Name`, into = "gene", sep = "\\(")
epilc_k27ac_cl2_filt = epilc_k27ac_cl2 %>% dplyr::filter(gene %in% zhifen$V1)


joined = epilc_k27ac_cl1_filt %>% inner_join(., epilc_k27ac_cl2_filt, by = "gene", suffix = c(".cl1", ".cl2")) %>% 
  dplyr::select(gene, p_value_cluster1 = "P-value.cl1", p_value_cluster2 = "P-value.cl2") %>% 
 mutate(p_value_cluster1 = -log10(p_value_cluster1),  p_value_cluster2 = -log10(p_value_cluster2)) 

hm_input = joined %>% pivot_longer(cols = "p_value_cluster1":"p_value_cluster2", names_to = "cluster", values_to = "-log10pvalue")
