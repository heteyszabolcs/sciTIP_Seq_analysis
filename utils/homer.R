# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("GenomicRanges")
  library("ggpubr")
  library("topGO")
  library("EnsDb.Mmusculus.v79")
})

# export
result_folder = "../results/HOMER/"

# Zhifen et al.
zhifen = fread("../data/Zhifen_et_al-pluripotency_gene_families_Fig1D.txt",
               header = FALSE)

# Yang et al. RNA-Seq
epilc_rna = fread("../results/bulk_RNA-Seq/Yang_et_al-EpiLC_RNA-Seq_TPMs.tsv")

# EpiLC H3.3 peaks
h33_peaks = "../results/SEACR/"
h33_peaks = list.files(h33_peaks, pattern = ".bed", full.names = TRUE)
create_h33_gr = function(peak_set, name) {
  peak_set = fread(peak_set)
  length = peak_set %>% mutate(diff = V3 - V2) %>% pull(diff) %>% mean
  sd = peak_set %>% mutate(diff = V3 - V2) %>% pull(diff) %>% sd
  print(
    glue(
      "Mean length is {as.character(round(length, 1))} 
    with SD {as.character(round(sd, 1))}"
    )
  )
  peak_set$V7 = name
  peak_set = GRanges(
    seqnames = peak_set$V1,
    ranges = IRanges(
      start = peak_set$V2,
      end = peak_set$V3,
      names = peak_set$V7
    )
  )
  return(peak_set)
}
h33_peaks_12h = create_h33_gr(peak_set = h33_peaks[1], name = "H3.3_12h")
h33_peaks_24h = create_h33_gr(peak_set = h33_peaks[2], name = "H3.3_24h")
h33_peaks_48h = create_h33_gr(peak_set = h33_peaks[3], name = "H3.3_48h")
h33_peaks_6h = create_h33_gr(peak_set = h33_peaks[4], name = "H3.3_6h")

# overlaps
h33_peaks = Reduce(subsetByOverlaps,
             list(h33_peaks_6h, h33_peaks_12h, h33_peaks_24h, h33_peaks_48h))

## HOMER results
# motif analysis
epilc_k27ac_cl1 = fread("../results/HOMER/Yang_k27ac-cluster1/homerResults.txt")
epilc_k27ac_cl2 = fread("../results/HOMER/Yang_k27ac-cluster2/homerResults.txt")
# annotation
cl1_annot = fread("../results/HOMER/k27ac_kmeans-cluster1-annotated.tsv")
cl1_annot = cl1_annot %>% dplyr::select(-starts_with("Peak"))
cl2_annot = fread("../results/HOMER/k27ac_kmeans-cluster2-annotated.tsv")
cl2_annot = cl2_annot %>% dplyr::select(-starts_with("Peak"))

joined = epilc_k27ac_cl2 %>% full_join(.,
                                       epilc_k27ac_cl1,
                                       by = "best_match",
                                       suffix = c(".cl2", ".cl1")) %>%
  dplyr::select(best_match,
                "H3K27ac cluster 1" = "P-value.cl1",
                "H3K27ac cluster 2" = "P-value.cl2") %>%
  mutate(
    `H3K27ac cluster 1` = -log10(`H3K27ac cluster 1`),
    `H3K27ac cluster 2` = -log10(`H3K27ac cluster 2`)
  )

hm_input = joined %>% pivot_longer(
  cols = "H3K27ac cluster 1":"H3K27ac cluster 2",
  names_to = "cluster",
  values_to = "-log10pvalue"
)
y_order = hm_input %>% group_by(best_match) %>% summarize(mean = mean(`-log10pvalue`)) %>% arrange(mean) %>% pull(best_match) %>% unique
y_order = factor(hm_input$best_match, levels = y_order)

hm = ggplot(hm_input, aes(x = cluster, y = y_order, fill = `-log10pvalue`)) +
  geom_tile(color = "black",
            lwd = 0.2,
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
  labs(fill = "-log10(p-value)", title = "Motif analysis (HOMER)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.6),
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
  glue(
    "{result_folder}Yang_et_al_K27_cluster_motif_analysis.pdf"
  ),
  plot = hm,
  width = 7,
  height = 7
)

# HOMER annotations
# distances to TSS
cl1_annot = cl1_annot %>% mutate(cluster = "cluster 1")
cl2_annot = cl2_annot %>% mutate(cluster = "cluster 2")
annots = rbind(cl1_annot, cl2_annot)
annots %>% dplyr::select(`Distance to TSS`, cluster) %>%
  ggplot(., aes(x = cluster, y = `Distance to TSS`, fill = cluster)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("#a6bddb", "#a6bddb")) +
  #ylim(0, 100) +
  labs(
    title = "Distance to TSS (HOMER annotation)",
    y = "distance to TSS (bp)",
    x = "",
    fill = "H3K27ac"
  ) +
  guides(fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black")
  ) +
  stat_compare_means(label.y = 1e+06, label.x = 1.2, size = 5) 

ggsave(
  glue(
    "{result_folder}H3K27ac_clusters-Distance_to_TSS-bp.pdf"
  ),
  plot = last_plot(),
  width = 6,
  height = 5
)

# create genomic ranges for EpiLC K27ac clusters
length = cl1_annot %>% mutate(diff = End - Start) %>% pull(diff) %>% mean
sd = cl1_annot %>% mutate(diff = End - Start) %>% pull(diff) %>% sd
print(
  glue(
    "Mean length is {as.character(round(length, 1))} 
    with SD {as.character(round(sd, 1))}"
  )
)

cl1_annot_gr = GRanges(
  seqnames = cl1_annot$Chr,
  ranges = IRanges(
    start = cl1_annot$Start,
    end = cl1_annot$End,
    names = cl1_annot$name
  )
)
cl2_annot_gr = GRanges(
  seqnames = cl2_annot$Chr,
  ranges = IRanges(
    start = cl2_annot$Start,
    end = cl2_annot$End,
    names = cl2_annot$name
  )
)

# EpiLC enhancer K27ac clusters overlaps with H3.3 peaks
ol = findOverlaps(cl2_annot_gr, h33_peaks, type = "any")
ol_cl2 = cl2_annot[queryHits(ol)]
ol_cl2_genes = ol_cl2 %>% dplyr::filter(abs(`Distance to TSS`) < 10000) %>% pull(`Gene Name`)
ol_cl2_genes_rna = epilc_rna %>% dplyr::filter(gene_name %in% ol_cl2_genes)

hm_input = ol_cl2_genes_rna %>% pivot_longer(
  cols = "TPM-EpiLC_1h":"TPM-EpiLC_72h",
  names_to = "timepoint",
  values_to = "TPM"
)
y_order = hm_input %>% group_by(gene_name) %>% 
  summarize(sd = sd(TPM)) %>% arrange(sd) %>% pull(gene_name) %>% unique
y_order = factor(hm_input$gene_name, levels = y_order)
x_order = factor(hm_input$timepoint, levels = unique(hm_input$timepoint))
cluster2_hm = ggplot(hm_input, aes(x = x_order, y = y_order, fill = log2(TPM+1))) +
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#075AFF",
    mid = "#FFFFCC",
    high = "#FF0000",
    midpoint = log2(max(ceiling(hm_input$TPM)))/2,
    limits = c(min(ceiling(hm_input$TPM)), log2(max(ceiling(hm_input$TPM)))+1)
  ) +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "log2(TPM)", title = "EpiLC enhancer cluster 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.6),
    axis.text.x = element_text(
      color = "black",
      size = 6,
      angle = 90,
      hjust = 1,
      vjust = 0.6
    ),
    axis.text.y = element_text(color = "black", size = 5),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()
print(cluster2_hm)

ggsave(
  glue(
    "{result_folder}H3K27ac_cluster2-RNASeq_TPMs_of_closestgenes.pdf"
  ),
  plot = last_plot(),
  width = 6,
  height = 7
)

ol = findOverlaps(cl1_annot_gr, h33_peaks, type = "any")
ol_cl1 = cl1_annot[queryHits(ol)]
ol_cl1_genes = ol_cl1 %>% dplyr::filter(abs(`Distance to TSS`) < 10000) %>% pull(`Gene Name`)
ol_cl1_genes_rna = epilc_rna %>% dplyr::filter(gene_name %in% ol_cl1_genes)

hm_input = ol_cl1_genes_rna %>% pivot_longer(
  cols = "TPM-EpiLC_1h":"TPM-EpiLC_72h",
  names_to = "timepoint",
  values_to = "TPM"
)
y_order = hm_input %>% group_by(gene_name) %>% 
  summarize(sd = sd(TPM)) %>% arrange(sd) %>% pull(gene_name) %>% unique
y_order = factor(hm_input$gene_name, levels = y_order)
x_order = factor(hm_input$timepoint, levels = unique(hm_input$timepoint))

cluster1_hm = ggplot(hm_input, aes(x = x_order, y = y_order, fill = log2(TPM+1))) +
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#075AFF",
    mid = "#FFFFCC",
    high = "#FF0000",
    midpoint = log2(max(ceiling(hm_input$TPM)))/2,
    limits = c(min(ceiling(hm_input$TPM)), log2(max(ceiling(hm_input$TPM)))+1)
  ) +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "log2(TPM)", title = "EpiLC enhancer cluster 1") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.6),
    axis.text.x = element_text(
      color = "black",
      size = 6,
      angle = 90,
      hjust = 1,
      vjust = 0.6
    ),
    axis.text.y = element_text(color = "black", size = 7),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()
print(cluster1_hm)

ggsave(
  glue(
    "{result_folder}H3K27ac_cluster1-RNASeq_TPMs_of_closestgenes.pdf"
  ),
  plot = last_plot(),
  width = 6,
  height = 7
)

# input: output of HOMER
annots_gr = GRanges(
  seqnames = annots$Chr,
  ranges = IRanges(
    start = annots$Start,
    end = annots$End,
    names = annots$cluster
  )
)
ol = findOverlaps(annots_gr, h33_peaks, type = "any")
epilc_enh_w_h33peaks = annots[queryHits(ol)]
epilc_enh_w_h33peaks_df = as.data.frame(epilc_enh_w_h33peaks)

create_input = function(homer_annot) {
  genes = homer_annot %>% dplyr::filter(abs(`Distance to TSS`) < 10000) %>%
    pull(`Gene Name`)
  
  input = rep(0, dim(epilc_enh_w_h33peaks)[1])
  names(input) = epilc_enh_w_h33peaks_df$`Gene Name`
  input[names(input) %in% genes] = 1
  
  return(input)
}

create_go_matrix = function(genes, colname) {
  # find biological process ontology
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    # use biological process ontology
    allGenes = genes,
    geneSelectionFun = function(x)
      (x == 1),
    annot = annFUN.org,
    mapping = "org.Mm.eg.db",
    ID = "symbol"
  )
  
  # Fisher test
  resultFisher <-
    runTest(GOdata, algorithm = "elim", statistic = "fisher")
  out <-
    GenTable(GOdata,
             Fisher = resultFisher,
             topNodes = 10,
             numChar = 60)
  
  out = out %>% dplyr::select(Term, Fisher)
  colnames(out) = c("Term", colname)
  
  return(out)
  
}

topgo_cl1 = create_go_matrix(
  create_input(
    homer_annot = cl1_annot
  ),
  colname = "K27ac cluster 1"
)
topgo_cl2 = create_go_matrix(
  create_input(
    homer_annot = cl2_annot
  ),
  colname = "K27ac cluster 2"
)

go_outputs = list(topgo_cl1, topgo_cl2)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
full = full %>% distinct(Term, .keep_all = TRUE) %>% 
  mutate(`K27ac cluster 1` = str_replace(`K27ac cluster 1`, "< ", ""), 
         `K27ac cluster 2` = str_replace(`K27ac cluster 2`, "< ", ""))

go_hm_input = full %>% pivot_longer(
  cols = `K27ac cluster 1`:`K27ac cluster 2`,
  names_to = "cluster",
  values_to = "p-value"
) %>% mutate(`p-value` = as.numeric(`p-value`)) %>% 
  dplyr::filter(Term != "biological_process")

y_order = go_hm_input %>% group_by(Term) %>% 
  mutate(Diff = `p-value` - lag(`p-value`)) %>% 
  drop_na() %>% arrange(Diff) %>% 
   pull(Term)
y_order = factor(go_hm_input$Term, levels = y_order)
x_order = factor(go_hm_input$cluster, levels = c("K27ac cluster 1", "K27ac cluster 2"))

go_hm = ggplot(go_hm_input, aes(x = x_order, y = y_order, fill = `p-value`)) +
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#de2d26",
    mid = "#fdbb84",
    high = "#fee8c8",
    midpoint = mean(c(min(go_hm_input$`p-value`), 0.010)),
    limits = c(min(go_hm_input$`p-value`), 0.010)
  ) +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "p-value", title = "topGO analysis", subtitle = "EpiLC enhancers with H3.3") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1),
    plot.subtitle = element_text(hjust = 0.8),
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
print(go_hm)

ggsave(
  glue(
    "{result_folder}H3K27ac_clusters-topGO_analysis.pdf"
  ),
  plot = last_plot(),
  width = 6,
  height = 5
)
