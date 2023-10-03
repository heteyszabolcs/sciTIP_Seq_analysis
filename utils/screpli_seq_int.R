# packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "glue",
  "tidyverse",
  "wigglescout",
  "GenomicFeatures",
  "GenomicRanges",
  "data.table",
  "purrr",
  "RColorBrewer",
  "ggpubr",
  "pROC",
  "ComplexHeatmap",
  "circlize",
  "caret",
  "caTools",
  "mlbench"
)

# output folder
result_folder = "../results/scRepli-Seq/"

# Seurat cluster bigwigs
cl0 = "../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster0_RPGC.bigwig"
cl1 = "../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster1_RPGC.bigwig"
cl2 = "../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster2_RPGC.bigwig"
cl3 = "../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster3_RPGC.bigwig"
cl4 = "../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster4_RPGC.bigwig"
cl5 = "../results/Seurat/cluster_bigwigs/H3.2_res0.5/cluster5_RPGC.bigwig"

bed = fread("../data/GSE108556_scRepli-Seq/early_25percentile_S_phase/GSM3521507_P361_39_1_80k_100S_MAP_2HMM2_eps0.01.binary.bedGraph")
bed = bed %>% dplyr::select(V1, V2, V3)
write_tsv(bed,
          "../results/scRepli-Seq/screpli_seq_80k_mm9_bins.bed",
          col_names = FALSE)

# H3.2 heatmap at 80 kb resolution
heatmap_h3.2 = function(chrom = "chr16") {
  bin_80 = bw_loci(c(cl0, cl1, cl2, cl3, cl4, cl5),
                   loci = "../results/scRepli-Seq/screpli_seq_80k_mm9_bins.bed")
  bin_80 = as.data.frame(bin_80)
  patterns = paste(c("_random", "X", "Y", "M"), collapse = "|")
  bin_80 = bin_80[grep(
    x = as.character(bin_80$seqnames),
    pattern = patterns,
    value = FALSE,
    invert = TRUE
  ), ]
  bin_80 = bin_80 %>% dplyr::select(-width,-strand) %>%
    mutate(
      seqnames = as.character(seqnames),
      start = as.character(start),
      end = as.character(end)
    )
  
  mat = bin_80 %>% dplyr::filter(seqnames == chrom)
  mat = mat %>% dplyr::select(starts_with("cluster"),-seqnames,-start,-end) %>%
    as.matrix %>% t
  rownames(mat) = c("cluster 0",
                      "cluster 1",
                      "cluster 2",
                      "cluster 3",
                      "cluster 4",
                      "cluster 5")
  col_fun = colorRamp2(c(0, 0.5, 1), c("white", "#fc9272", "red"))
  pdf(
    file = glue("{result_folder}H3.2_RPGC_{chrom}_80kb_hm.pdf"),
    width = 6,
    height = 3
  )
  hm = Heatmap(
    mat,
    column_title = glue("{chrom} - H3.2 80 kb bins"),
    row_title = "",
    name = "H3.2 RPGC",
    col = col_fun,
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    #heatmap_width = unit(3, "cm"),
    heatmap_height = unit(3, "cm"),
    show_row_names = TRUE,
    show_column_names = FALSE,
    show_heatmap_legend = TRUE,
    use_raster = TRUE
  )
  print(hm)
  dev.off()
}

# binary heatmap about scRepli-Seq
# stage: type of scRepli-Seq sample (early, mid, late S phase)
# chrom: chromosome of interest of the mm9
heatmap_screpliseq = function(stage = "early", chrom = "chr16") {
  # scRepli-seq integration
  screpliseq_mids = list.files(
    "../data/GSE108556_scRepli-Seq/midS_phase/",
    patter = "binary.bedGraph",
    full.names = TRUE
  )
  screpliseq_earlys = list.files(
    "../data/GSE108556_scRepli-Seq/early_25percentile_S_phase/",
    patter = "binary.bedGraph",
    full.names = TRUE
  )
  screpliseq_lates = list.files(
    "../data/GSE108556_scRepli-Seq/late_75percentile_S_phase/",
    patter = "binary.bedGraph",
    full.names = TRUE
  )
  
  list_of_inputs = function(stage = stage) {
    if (stage == "early") {
      data = screpliseq_earlys
      return(data)
    }
    if (stage == "mid") {
      data = screpliseq_mids
      return(data)
    }
    if (stage == "late") {
      data = screpliseq_lates
      return(data)
    }
  }
  
  screpli_seqs = lapply(list_of_inputs(stage = stage), fread)
  screpli_seq = screpli_seqs %>%
    purrr::reduce(left_join, by = c("V1", "V2", "V3"))
  colnames(screpli_seq) = c("chr", "start", "end",
                            sapply(seq(1, length(screpli_seqs)), function(x) {
                              paste0(stage, "_", x)
                            }))
  
  mat = screpli_seq %>% dplyr::filter(chr == chrom) %>%
    dplyr::select(-chr,-start,-end) %>%
    as.matrix %>% t
  rownames(mat) = c(
    glue("{stage} 1"),
    glue("{stage} 2"),
    glue("{stage} 3"),
    glue("{stage} 4"),
    glue("{stage} 5")
  )
  col_fun = colorRamp2(c(-1, 0, 1), c("black", "grey", "#fec44f"))
  lgd3 = Legend(
    labels = c("repl.", "no data", "non-repl."),
    legend_gp = gpar(
      fill = 7:9,
      color = c("black", "grey", "#fec44f")
    ),
    title = "state"
  )
  
  pdf(
    file = glue("{result_folder}scRepliSeq_{stage}_{chrom}_bin_hm.pdf"),
    width = 6,
    height = 3
  )
  hm = Heatmap(
    mat,
    column_title = glue("{chrom} - scRepli-Seq, {stage} S phase"),
    row_title = "",
    name = "scRepli-Seq",
    col = col_fun,
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    #heatmap_width = unit(3, "cm"),
    heatmap_height = unit(3, "cm"),
    show_row_names = TRUE,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE,
    use_raster = TRUE
  )
  draw(hm)
  draw(lgd3, x = unit(0.85, "npc"), y = unit(0.15, "npc"))
  dev.off()
  
}

logreg_model = function(stage = "early",
                        chrom = "chr1",
                        only.coefs = FALSE,
                        bal.acc. = TRUE,
                        sensitivity = TRUE
) {
  # scRepli-seq integration
  screpliseq_mids = list.files(
    "../data/GSE108556_scRepli-Seq/midS_phase/",
    patter = "binary.bedGraph",
    full.names = TRUE
  )
  screpliseq_earlys = list.files(
    "../data/GSE108556_scRepli-Seq/early_25percentile_S_phase/",
    patter = "binary.bedGraph",
    full.names = TRUE
  )
  screpliseq_lates = list.files(
    "../data/GSE108556_scRepli-Seq/late_75percentile_S_phase/",
    patter = "binary.bedGraph",
    full.names = TRUE
  )
  
  list_of_inputs = function(stage = stage) {
    if (stage == "early") {
      data = screpliseq_earlys
      return(data)
    }
    if (stage == "mid") {
      data = screpliseq_mids
      return(data)
    }
    if (stage == "late") {
      data = screpliseq_lates
      return(data)
    }
  }
  
  bed = fread(list_of_inputs(stage = stage)[1])
  bed = bed %>% dplyr::select(V1, V2, V3)
  write_tsv(bed,
            "../results/scRepli-Seq/temp_bed.bed",
            col_names = FALSE)
  
  bin_80 = bw_loci(c(cl0, cl1, cl2, cl3, cl4, cl5),
                   loci = "../results/scRepli-Seq/temp_bed.bed")
  bin_80 = as.data.frame(bin_80)
  patterns = paste(c("_random", "X", "Y", "M"), collapse = "|")
  bin_80 = bin_80[grep(
    x = as.character(bin_80$seqnames),
    pattern = patterns,
    value = FALSE,
    invert = TRUE
  ), ]
  bin_80 = bin_80 %>% dplyr::select(-width,-strand) %>%
    mutate(
      seqnames = as.character(seqnames),
      start = as.character(start),
      end = as.character(end)
    )
  
  screpli_seqs = lapply(list_of_inputs(stage = stage), fread)
  screpli_seq = screpli_seqs %>%
    purrr::reduce(inner_join, by = c("V1", "V2", "V3"))
  colnames(screpli_seq) = c("chr", "start", "end",
                            sapply(seq(1, length(screpli_seqs)), function(x) {
                              paste0(stage, "_", x)
                            }))
  screpli_seq = screpli_seq[grep(
    x = as.character(screpli_seq$chr),
    pattern = patterns,
    value = FALSE,
    invert = TRUE
  ), ]
  
  # find those regions that replicate in all cells (these are class 1!)
  conservative = screpli_seq %>%
    mutate(sum = rowSums(across(starts_with(glue("{stage}_"))))) %>% 
    mutate(conservative = ifelse(sum >= 2, 1, 0)) %>% pull(conservative)
  
  bin_80 = bin_80 %>% mutate(target = conservative)
  mat = bin_80 %>% dplyr::filter(seqnames == chrom) %>% 
    unite(seqnames, seqnames:end, sep = "_")
  
  ## make logistic regression model (glm)
  # splitting the data set into ratio 0.80:0.20
  set.seed(42)
  split = sample.split(mat$target, SplitRatio = 0.80)
  # creating training dataset and test dataset (by caTools)
  mat = mat %>% dplyr::select(-seqnames) %>% 
    mutate(target = as.factor(target))
  training = subset(mat, split == TRUE) %>% na.omit()
  test = subset(mat, split == FALSE) %>% na.omit()
  
  # logreg model with cross validations by caret 
  glm = train(
    form = target ~ .,
    data = training,
    trControl = trainControl(method = "cv", number = 5),
    method = "glm",
    family = "binomial"
  )
  pred = predict(glm, test)
  
  if(only.coefs == TRUE) {
    coefs = coefs %>% mutate(chrom = chrom, stage = stage)
    return(coefs)
  }
  
  # performance metrics
  if (bal.acc. == TRUE) {
    balanced_acc =
      round(confusionMatrix(data = pred, reference = test$target)$byClass[["Balanced Accuracy"]],
            2)
    print(balanced_acc)
    balanced_acc = tibble(chrom = chrom,
                          stage = stage,
                          balanced_acc = balanced_acc)
    return(balanced_acc)
  }
  
  if (sensitivity == TRUE) {
    sensitivity =
      round(confusionMatrix(data = pred, reference = test$target)$byClass[["Sensitivity"]],
            2)
    print(sensitivity)
    sensitivity = tibble(chrom = chrom,
                          stage = stage,
                          sensitivity = sensitivity)
    return(sensitivity)
  }
  
  # logreg model by glm function (base R)
  model_glm = glm(data = training, formula = target ~ ., family = "binomial")
  print(model_glm)
  coefs = tibble(coef = coef(model_glm), cluster = names(coef(model_glm)))
  coefs = coefs %>% dplyr::filter(cluster != "(Intercept)")
  
  # return coefficients
  if(only.coefs == TRUE) {
    coefs = coefs %>% mutate(chrom = chrom, stage = stage)
    return(coefs)
  }
  
  coef_bars = ggplot(data = coefs, aes(x = cluster, y = coef)) +
    geom_bar(
      stat = "identity",
      color = "steelblue",
      fill = "#deebf7",
      linewidth = 1.2
    ) +
    ylim(-50, 50) +
    theme_minimal() +
    labs(title = glue("coefficients (logreg model)"),
         subtitle = glue("{chrom}, S phase: {stage}"),
         x = "",
         y = "coefficient") +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(
        size = 13,
        color = "black",
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )
  print(coef_bars)
  ggsave(
    plot = last_plot(),
    glue("{result_folder}logreg_model_{stage}_{chrom}_coef_bars.pdf"),
    width = 5,
    height = 5,
    dpi = 500,
  )
  
  
  # predict and plot ROC curve
  test_prob = predict(model_glm, newdata = test, type = "response")
  pdf(
    file = glue("{result_folder}logreg_model_{stage}_{chrom}_ROCcurve.pdf"),
    width = 4,
    height = 4
  )
  roc = roc(test$target ~ test_prob, plot = TRUE, print.auc = TRUE)
  dev.off()
  
  return(roc)
  
}

# run functions on representative chromosomes
heatmap_h3.2(chrom = "chr1")
heatmap_screpliseq(stage = "early", chrom = "chr1")
heatmap_screpliseq(stage = "mid", chrom = "chr1")
heatmap_screpliseq(stage = "late", chrom = "chr1")

heatmap_h3.2(chrom = "chr16")
heatmap_screpliseq(stage = "early", chrom = "chr16")
heatmap_screpliseq(stage = "mid", chrom = "chr16")
heatmap_screpliseq(stage = "late", chrom = "chr16")

logreg_model(stage = "early", chrom = "chr1")
logreg_model(stage = "mid", chrom = "chr1")
logreg_model(stage = "late", chrom = "chr1")

logreg_model(stage = "early", chrom = "chr16")
logreg_model(stage = "mid", chrom = "chr16")
logreg_model(stage = "late", chrom = "chr16")

# collect all coefficients for each chromosome and stages
chroms = unname(sapply(seq(1, 19, 1), function(x) {paste0("chr", as.character(x))}))
stages = c("early", "mid", "late")

balanced_accs = lapply(chroms, function(x) {
  lapply(stages, function(y) {
    print(paste0(x, " - ", y))
    output = logreg_model(chrom = x,
                          stage = y,
                          only.coefs = FALSE,
                          bal.acc. = TRUE)
    return(output)
  })
})

sensitivities = lapply(chroms, function(x) {
  lapply(stages, function(y) {
    print(paste0(x, " - ", y))
    output = logreg_model(chrom = x,
                          stage = y,
                          only.coefs = FALSE,
                          bal.acc. = FALSE,
                          sensitivity = TRUE)
    return(output)
  })
})

coefs = lapply(chroms, function(x) {
  lapply(stages, function(y) {
    print(paste0(x, " - ", y))
    output = logreg_model(chrom = x,
                 stage = y,
                 only.coefs = TRUE)
    return(output)
  })
})

# visualize distribution of coefficients
vis_input = bind_rows(coefs)

order = factor(vis_input$stage, levels = c("early", "mid", "late"))
ggplot(vis_input, aes(x = cluster, y = coef, fill = order)) +
  geom_boxplot(alpha = 1) +
  scale_color_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  ylim(-50, 50) +
  labs(
    title = "logistic regression coefficients",
    x = "",
    y = "coefficient",
    fill = "S phase"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 0.9, vjust = 1),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) 

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_coefficients_bp.pdf"),
  width = 5,
  height = 5,
  dpi = 500,
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_coefficients_bp.png"),
  width = 5,
  height = 5
)

order = factor(vis_input$stage, levels = c("early", "mid", "late"))
facet = ggplot(vis_input, aes(x = order, y = coef, color = stage)) +
  geom_boxplot() +
  facet_wrap( ~ cluster) +
  scale_color_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  ylim(-50, 50) +
  labs(
    title = "logistic regression coefficients",
    subtitle = "(over all chromosomes)",
    x = "",
    y = "coefficient",
    color = "S phase"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  guides(color = "none") +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(
      size = 15,
      color = "black",
      angle = 45,
      hjust = 0.9,
      vjust = 1
    ),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) + stat_compare_means(label = "p.signif", label.y = 45, label.x = 1.9)
facet

ggsave(
  plot = facet,
  glue("{result_folder}logreg_model_coefficients_facet.pdf"),
  width = 7,
  height = 7,
  dpi = 500,
)

ggsave(
  plot = facet,
  glue("{result_folder}logreg_model_coefficients_facet.png"),
  width = 7,
  height = 7
)

# visualize distribution of balanced accuracy
vis_input_bacc = bind_rows(balanced_accs)

order = factor(vis_input_bacc$stage, levels = c("early", "mid", "late"))
comps = list( c("early", "late"), c("early", "mid"), c("late", "mid") )
ggplot(vis_input_bacc, aes(x = order, y = balanced_acc, fill = order)) +
  geom_boxplot(alpha = 1) +
  scale_color_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  ylim(0, 1.1) +
  labs(
    title = "distribution of balanced accuracy",
    x = "",
    y = "balanced accuracy",
    fill = "S phase"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    axis.text.x = element_text(size = 18, color = "black", angle = 0, hjust = 0.9, vjust = 1),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) +
  stat_compare_means(comparisons = comps, label.y = c(0.92, 0.99, 1.07)) +
  stat_compare_means(label.y = 1.2)

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_balanced_accs_bp.pdf"),
  width = 5,
  height = 5,
  dpi = 500,
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_balanced_accs_bp.png"),
  width = 5,
  height = 5
)

order = factor(vis_input_bacc$chrom, levels =
                 sapply(seq(1, 19, 1), function(x) {
                   paste0("chr", as.character(x))
                 }))
ggplot(vis_input_bacc, aes(x = order, y = balanced_acc)) +
  geom_boxplot(alpha = 1, fill = "#9ecae1") +
  ylim(0, 1.0) +
  labs(
    title = "distribution of balanced accuracy",
    x = "",
    y = "balanced accuracy",
    fill = "S phase"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 0.9, vjust = 1),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) 

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_balanced_accs_chromlevel_bp.pdf"),
  width = 10,
  height = 5,
  dpi = 500,
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_balanced_accs_chromlevel_bp.png"),
  width = 10,
  height = 5
)

# visualize distribution of balanced accuracy
vis_input_sens = bind_rows(sensitivities)
order = factor(vis_input_sens$stage, levels = c("early", "mid", "late"))
comps = list( c("early", "late"), c("early", "mid"), c("late", "mid") )
ggplot(vis_input_sens, aes(x = order, y = sensitivity, fill = order)) +
  geom_boxplot(alpha = 1) +
  scale_color_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  ylim(0, 1.1) +
  labs(
    title = "distribution of true positive rates",
    x = "",
    y = "sensitivity (true positive rate)",
    fill = "S phase"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    axis.text.x = element_text(size = 18, color = "black", angle = 0, hjust = 0.9, vjust = 1),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) +
  stat_compare_means(comparisons = comps, label.y = c(1.10, 1.04, 0.96)) +
  stat_compare_means(label.y = 1.20)
  
ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_sens_bp.pdf"),
  width = 5,
  height = 5,
  dpi = 500,
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_sens_bp.png"),
  width = 5,
  height = 5
)

order = factor(vis_input_sens$chrom, levels =
                 sapply(seq(1, 19, 1), function(x) {
                   paste0("chr", as.character(x))
                 }))
ggplot(vis_input_sens, aes(x = order, y = sensitivity)) +
  geom_boxplot(alpha = 1, fill = "#9ecae1") +
  ylim(0, 1.0) +
  labs(
    title = "distribution of true positive rates",
    x = "",
    y = "sensitivity (true positive rate)",
    fill = "S phase"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 0.9, vjust = 1),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) 

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_sens_chromlevel_bp.pdf"),
  width = 10,
  height = 5,
  dpi = 500,
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}logreg_model_sens_chromlevel_bp.png"),
  width = 10,
  height = 5
)
