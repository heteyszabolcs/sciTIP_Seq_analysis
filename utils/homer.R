# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
})

result_folder = "../results/HOMER/"

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
  dplyr::select(gene, "H3K27ac cluster 1" = "P-value.cl1",  "H3K27ac cluster 2" = "P-value.cl2") %>% 
 mutate(`H3K27ac cluster 1` = -log10(`H3K27ac cluster 1`),  `H3K27ac cluster 2` = -log10(`H3K27ac cluster 2`)) 

hm_input = joined %>% pivot_longer(cols = "H3K27ac cluster 1":"H3K27ac cluster 2", names_to = "cluster", values_to = "-log10pvalue")
y_order = hm_input %>% group_by(gene) %>% summarize(mean = mean(`-log10pvalue`)) %>% arrange(mean) %>% pull(gene) %>% unique 
y_order = factor(hm_input$gene, levels = y_order)

hm = ggplot(hm_input, aes(x = cluster, y = y_order, fill = `-log10pvalue`)) +
  geom_tile(color = "white",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#075AFF",
    mid = "#FFFFCC",
    high = "#FF0000",
    midpoint = 30,
    limits = c(0, 60)
  ) +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "-log10(p-value)", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 12,
      angle = 90,
      hjust = 1,
      vjust = 0.6
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()
print(hm)

ggsave(
  glue("{result_folder}Zhifen_et_al-pluripotency_genes_motif_analysis-OSTN.pdf"),
  plot = hm,
  width = 5,
  height = 7
)

hm_input = hm_input %>% dplyr::filter(!gene == "OCT4-SOX2-TCF-NANOG")
y_order = hm_input %>% group_by(gene) %>% summarize(mean = mean(`-log10pvalue`)) %>% arrange(mean) %>% pull(gene) %>% unique 
y_order = factor(hm_input$gene, levels = y_order)

hm2 = ggplot(hm_input, aes(x = cluster, y = y_order, fill = `-log10pvalue`)) +
  geom_tile(color = "white",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#075AFF",
    mid = "#FFFFCC",
    high = "#FF0000",
    midpoint = 20,
    limits = c(0, 40)
  ) +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "-log10(p-value)", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 12,
      angle = 90,
      hjust = 1,
      vjust = 0.6
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()
print(hm2)

ggsave(
  glue("{result_folder}Zhifen_et_al-pluripotency_genes_motif_analysis.pdf"),
  plot = hm2,
  width = 5,
  height = 7
)





