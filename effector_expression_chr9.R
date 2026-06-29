library(data.table)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

setwd("~/Desktop/BLR_data/New_518_560/")

# Inputs

effector_list_file <- "560A_chr9_effector_list.txt"
effector_bed_file  <- "effector_genes_chr9A_by_suffix.bed"
expr_file          <- "Ph560A_DESeq2_normalized_counts.tsv"
te_file            <- "560A_TE_all_chr13Aupdated.txt"

out_prefix <- "Chr9A_TE_border_effector_expression"

old_col  <- "#4D6FA3"
new_col  <- "#C7646A"
te_col   <- "#8C72B1"
grey_col <- "#D9D9D9"


eff_list <- fread(effector_list_file)
setnames(eff_list, 1, "gene_id")

eff_bed <- fread(effector_bed_file, header = FALSE)
# BED col names as chr, start, end, gene_id
eff_bed <- eff_bed[, .(
  chr = V1,
  start = as.integer(V2),
  end = as.integer(V3),
  gene_id = V4
)]

te <- fread(te_file, header = FALSE)

te <- te[, .(
  chr = V1,
  start = as.integer(V2),
  end = as.integer(V3),
  TE_class = V4
)]

expr <- fread(expr_file)

setnames(expr, 1, "gene_id")


# Keep Chr9A effectors

eff <- eff_bed %>%
  filter(gene_id %in% eff_list$gene_id)

# Calculate nearest TE distance

eff_gr <- GRanges(
  seqnames = eff$chr,
  ranges = IRanges(eff$start, eff$end),
  gene_id = eff$gene_id
)

te <- te[chr == "Chr9A"]
te <- te[!is.na(start) & !is.na(end) & start < end]

te_gr <- GRanges(
  seqnames = te$chr,
  ranges = IRanges(start = te$start, end = te$end),
  TE_class = te$TE_class
)

nearest_te <- distanceToNearest(eff_gr, te_gr, ignore.strand = TRUE)

te_dist <- data.table(
  gene_id = mcols(eff_gr)$gene_id[queryHits(nearest_te)],
  nearest_TE_class = mcols(te_gr)$TE_class[subjectHits(nearest_te)],
  dist_to_TE = mcols(nearest_te)$distance
)

eff_anno <- as.data.table(eff) %>%
  left_join(te_dist, by = "gene_id") %>%
  mutate(
    TE_border = case_when(
      dist_to_TE == 0 ~ "overlapping TE",
      dist_to_TE <= 2000 ~ "<=2 kb from TE",
      dist_to_TE <= 5000 ~ "2–5 kb from TE",
      TRUE ~ ">5 kb from TE"
    )
  )


# Expression matrix

expr_eff <- expr %>%
  filter(gene_id %in% eff_anno$gene_id) %>%
  left_join(eff_anno, by = "gene_id")

expr_cols <- grep("^Ph[0-9]+_Rep[0-9]+$", colnames(expr_eff), value = TRUE)

expr_long <- expr_eff %>%
  pivot_longer(
    cols = all_of(expr_cols),
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  mutate(
    isolate_group = case_when(
      grepl("^Ph491|^Ph518", sample) ~ "older isolates",
      grepl("^Ph612|^Ph685", sample) ~ "CNV-gained isolates",
      TRUE ~ "other"
    ),
    log2_expr = log2(norm_count + 1)
  ) %>%
  filter(isolate_group %in% c("older isolates", "CNV-gained isolates"))


# Summarise expression per gene, old vs new isolayes

gene_summary <- expr_long %>%
  group_by(gene_id, TE_border, dist_to_TE, nearest_TE_class, isolate_group) %>%
  summarise(mean_log2_expr = mean(log2_expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = isolate_group,
    values_from = mean_log2_expr
  ) %>%
  mutate(
    delta_new_minus_old = `CNV-gained isolates` - `older isolates`
  )


# Gene-level summary old to new isolates shift

gene_summary <- expr_long %>%
  group_by(
    gene_id, chr, start, end,
    TE_border, dist_to_TE, nearest_TE_class,
    isolate_group
  ) %>%
  summarise(
    mean_log2_expr = mean(log2_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = isolate_group,
    values_from = mean_log2_expr
  ) %>%
  mutate(
    # IMPORTANT: old -> new direction
    delta_old_to_new = `older isolates` - `CNV-gained isolates`,
    mid = (start + end) / 2,
    mean_expr_all = rowMeans(
      cbind(`older isolates`, `CNV-gained isolates`),
      na.rm = TRUE
    )
  )


# plot genomic expression shift, older -> new

blue <- "#4D6FA3"   
red  <- "#C7646A"  

pA <- ggplot(gene_summary, aes(x = mid / 1e6, y = delta_old_to_new)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.45) +
  geom_point(
    aes(
      size = mean_expr_all,
      colour = delta_old_to_new
    ),
    alpha = 0.85
  ) +
  scale_colour_gradient2(
    low = blue,
    mid = "grey85",
    high = red,
    midpoint = 0,
    name = "Expression shift\nolder - new"
  ) +
  scale_size_continuous(
    name = "Mean expression\nlog2(count + 1)",
    range = c(1.8, 5.8)
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Chr9A position (Mb)",
    y = "Mean log2 expression difference\nolder - new",
    title = "A  TE-proximal Chr9A effector expression shifts"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )



# expression in each isolate, average replicates first

expr_gene_isolate <- expr_long %>%
  mutate(
    isolate = sub("_Rep[0-9]+$", "", sample),
    isolate = factor(
      isolate,
      levels = c("Ph491", "Ph518", "Ph612", "Ph685")
    ),
    isolate_label = case_when(
      isolate == "Ph491" ~ "Ph491\n(older)",
      isolate == "Ph518" ~ "Ph518\n(older)",
      isolate == "Ph612" ~ "Ph612\n(new)",
      isolate == "Ph685" ~ "Ph685\n(new)"
    ),
    isolate_label = factor(
      isolate_label,
      levels = c(
        "Ph491\n(older)",
        "Ph518\n(older)",
        "Ph612\n(new)",
        "Ph685\n(new)"
      )
    )
  ) %>%
  group_by(gene_id, isolate, isolate_label, isolate_group) %>%
  summarise(
    mean_log2_expr = mean(log2_expr, na.rm = TRUE),
    .groups = "drop"
  )

# Mixed model statistics oldvs new comparison

library(lme4)
library(lmerTest)
library(emmeans)

expr_gene_isolate <- expr_long %>%
  mutate(
    isolate = sub("_Rep[0-9]+$", "", sample),
    isolate = factor(isolate, levels = c("Ph491", "Ph518", "Ph612", "Ph685")),
    isolate_group = factor(
      isolate_group,
      levels = c("older isolates", "CNV-gained isolates")
    ),
    isolate_label = case_when(
      isolate == "Ph491" ~ "Ph491\n(older)",
      isolate == "Ph518" ~ "Ph518\n(older)",
      isolate == "Ph612" ~ "Ph612\n(new)",
      isolate == "Ph685" ~ "Ph685\n(new)"
    ),
    isolate_label = factor(
      isolate_label,
      levels = c("Ph491\n(older)", "Ph518\n(older)", "Ph612\n(new)", "Ph685\n(new)")
    )
  ) %>%
  group_by(gene_id, isolate, isolate_label, isolate_group) %>%
  summarise(
    mean_log2_expr = mean(log2_expr, na.rm = TRUE),
    .groups = "drop"
  )

# mixed model old/new effect
model_group <- lmer(
  mean_log2_expr ~ isolate_group + (1 | gene_id),
  data = expr_gene_isolate
)

anova_group <- anova(model_group)
group_p <- anova_group["isolate_group", "Pr(>F)"]

emm_group <- emmeans(model_group, ~ isolate_group)
group_contrast <- pairs(emm_group)
print(emm_group)
print(group_contrast)

# group letters a or b if significant
group_letters <- data.frame(
  isolate_group = c("older isolates", "CNV-gained isolates"),
  group_letter = if (group_p < 0.05) c("a", "b") else c("a", "a")
)

label_pos <- expr_gene_isolate %>%
  group_by(isolate_label, isolate_group) %>%
  summarise(
    y = max(mean_log2_expr, na.rm = TRUE) + 0.8,
    .groups = "drop"
  ) %>%
  left_join(group_letters, by = "isolate_group")


pB <- ggplot(
  expr_gene_isolate,
  aes(x = isolate_label, y = mean_log2_expr, fill = isolate_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.58,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.16,
    size = 1.1,
    alpha = 0.4,
    colour = "black"
  ) +
  geom_text(
    data = label_pos,
    aes(x = isolate_label, y = y, label = group_letter),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2.5,
    y = max(label_pos$y, na.rm = TRUE) + 0.9,
    label = paste0(
      "Mixed model old vs new p = ",
      formatC(group_p, format = "e", digits = 2)
    ),
    size = 3.8
  ) +
  scale_fill_manual(values = c(
    "older isolates" = blue,
    "CNV-gained isolates" = red
  )) +
  coord_cartesian(
    ylim = c(
      min(expr_gene_isolate$mean_log2_expr, na.rm = TRUE),
      max(label_pos$y, na.rm = TRUE) + 1.6
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Mean DESeq2-normalized expression\nlog2(count + 1)",
    title = "B  Expression of Chr9A effectors in each isolate"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10)
  )


# heatmap of expressed/most shifted effectors

library(pheatmap)

all_effector_genes <- gene_summary %>%
  arrange(start) %>%
  pull(gene_id)

heat_dt <- expr_long %>%
  filter(gene_id %in% all_effector_genes) %>%
  group_by(gene_id, isolate) %>%
  summarise(
    mean_expr = mean(log2_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = isolate,
    values_from = mean_expr
  )

heat_dt <- heat_dt %>%
  mutate(gene_id = factor(gene_id, levels = all_effector_genes)) %>%
  arrange(gene_id)

heat_mat <- as.matrix(heat_dt[, -1])
rownames(heat_mat) <- heat_dt$gene_id

heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z[is.na(heat_mat_z)] <- 0

ann_col <- data.frame(
  Group = ifelse(
    colnames(heat_mat_z) %in% c("Ph491", "Ph518"),
    "older isolates",
    "CNV-gained isolates"
  )
)
rownames(ann_col) <- colnames(heat_mat_z)

ann_colors <- list(
  Group = c(
    "older isolates" = "#4D6FA3",
    "CNV-gained isolates" = "#C7646A"
  )
)

pdf("Chr9A_all_TE_proximal_effector_expression_heatmap.pdf",
    width = 7,
    height = 13)

my_colors <- colorRampPalette(c("#7B8DBF", "white", "#57B893"))(100)

pheatmap(
  heat_mat_z,
  color = my_colors,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 5.2,
  fontsize_col = 10,
  border_color = NA,
  main = "TE-proximal Chr9A effector expression"
)

dev.off()

#############
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(scales)


# reuse te_chr9, all_genes, eff_anno

te_chr9 <- te_chr9[chr == "Chr9A" & !is.na(start) & !is.na(end) & start < end]
all_genes <- all_genes[chr == "Chr9A" & !is.na(start) & !is.na(end) & start < end]

# distance to nearest left/right TE

get_flanking_te_dist <- function(gene_dt, te_dt, label) {
  gene_dt <- as.data.table(gene_dt)
  te_dt <- as.data.table(te_dt)
  
  out <- gene_dt[, {
    g_start <- start
    g_end <- end
    
    left_te <- te_dt[end < g_start]
    right_te <- te_dt[start > g_end]
    
    left_dist <- if (nrow(left_te) > 0) min(g_start - left_te$end) else NA_real_
    right_dist <- if (nrow(right_te) > 0) min(right_te$start - g_end) else NA_real_
    
    .(
      gene_id = gene_id,
      left_TE_dist_bp = left_dist,
      right_TE_dist_bp = right_dist
    )
  }, by = .(chr, start, end)]
  
  out[, group := label]
  out[, left_TE_dist_kb := left_TE_dist_bp / 1000]
  out[, right_TE_dist_kb := right_TE_dist_bp / 1000]
  
  out
}

all_gene_te <- get_flanking_te_dist(all_genes, te_chr9, "All genes")
eff_te <- get_flanking_te_dist(eff_anno, te_chr9, "Candidate effectors")

plot_dt <- rbindlist(list(all_gene_te, eff_te), fill = TRUE)
plot_dt <- plot_dt[
  !is.na(left_TE_dist_kb) &
    !is.na(right_TE_dist_kb) &
    left_TE_dist_kb > 0 &
    right_TE_dist_kb > 0
]


plot_density_panel <- function(dt, title_text, max_count = NULL) {
  p <- ggplot(dt, aes(x = left_TE_dist_kb, y = right_TE_dist_kb)) +
    geom_hex(bins = 35) +
    scale_x_log10(
      limits = c(0.01, 1000),
      breaks = c(0.01, 0.1, 1, 10, 100, 1000),
      labels = c("0.01", "0.1", "1", "10", "100", "1000")
    ) +
    scale_y_log10(
      limits = c(0.01, 1000),
      breaks = c(0.01, 0.1, 1, 10, 100, 1000),
      labels = c("0.01", "0.1", "1", "10", "100", "1000")
    ) +
    scale_fill_gradient(
      low = "#7B8DBF",
      high = "#57B893",
      name = "Gene\ncount"
    ) +
    theme_bw(base_size = 13) +
    labs(
      title = title_text,
      x = "Distance to nearest left-side TE (kb)",
      y = "Distance to nearest right-side TE (kb)"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  p
}

p_all <- plot_density_panel(
  plot_dt[group == "All genes"],
  "All genes"
)

p_eff <- plot_density_panel(
  plot_dt[group == "Candidate effectors"],
  "Candidate effectors"
)

p_TE_density <- p_all + p_eff + plot_annotation(tag_levels = "A")


# Chr9A TE-distance analysis for non-effector genes vs candidate effectors

library(data.table)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/BLR_data/New_518_560/")

te_file        <- "560A_TE_all_chr13Aupdated.txt"
all_gene_file  <- "Chr9A_all_genes.bed"
effector_file  <- "effector_genes_chr9A_by_suffix.bed"

out_prefix <- "Chr9A_effector_TE_distance_overlap_aware"

gene_col <- "#4D6FA3"
eff_col  <- "#C7646A"

plot_floor_kb <- 0.001

te <- fread(te_file, header = FALSE, fill = TRUE)

te <- te[, .(
  chr = V1,
  start = suppressWarnings(as.numeric(V2)),
  end = suppressWarnings(as.numeric(V3)),
  TE_class = V4
)]

te_chr9 <- te[
  chr == "Chr9A" &
    !is.na(start) &
    !is.na(end) &
    start < end
]

cat("Chr9A TE features:", nrow(te_chr9), "\n")

all_genes <- fread(all_gene_file, header = FALSE, fill = TRUE)

all_genes <- all_genes[, .(
  chr = V1,
  start = suppressWarnings(as.numeric(V2)),
  end = suppressWarnings(as.numeric(V3)),
  gene_id = V4
)]

all_genes <- all_genes[
  chr == "Chr9A" &
    !is.na(start) &
    !is.na(end) &
    start < end
]

cat("Chr9A all genes:", nrow(all_genes), "\n")

# candidate effectors

effectors <- fread(effector_file, header = FALSE, fill = TRUE)

effectors <- effectors[, .(
  chr = V1,
  start = suppressWarnings(as.numeric(V2)),
  end = suppressWarnings(as.numeric(V3)),
  gene_id = V4
)]

effectors <- effectors[
  chr == "Chr9A" &
    !is.na(start) &
    !is.na(end) &
    start < end
]

cat("Chr9A candidate effectors:", nrow(effectors), "\n")

get_flanking_te_dist <- function(gene_dt, te_dt, group_name) {
  
  gene_dt <- as.data.table(copy(gene_dt))
  te_dt <- as.data.table(copy(te_dt))
  
  out <- gene_dt[, {
    
    g_start <- start
    g_end <- end
    
    left_te    <- te_dt[end < g_start]
    right_te   <- te_dt[start > g_end]
    overlap_te <- te_dt[start <= g_end & end >= g_start]
    
    left_dist <- if (nrow(left_te) > 0) {
      min(g_start - left_te$end)
    } else {
      NA_real_
    }
    
    right_dist <- if (nrow(right_te) > 0) {
      min(right_te$start - g_end)
    } else {
      NA_real_
    }
    
    overlaps_TE <- nrow(overlap_te) > 0
    
    nearest_dist <- if (overlaps_TE) {
      0
    } else {
      min(c(left_dist, right_dist), na.rm = TRUE)
    }
    
    if (!is.finite(nearest_dist)) nearest_dist <- NA_real_
    
    .(
      gene_id = gene_id,
      left_TE_dist_bp = left_dist,
      right_TE_dist_bp = right_dist,
      overlaps_TE = overlaps_TE,
      nearest_TE_bp = nearest_dist
    )
    
  }, by = .(chr, start, end)]
  
  out[, group := group_name]
  
  out[, left_TE_dist_kb := left_TE_dist_bp / 1000]
  out[, right_TE_dist_kb := right_TE_dist_bp / 1000]
  out[, nearest_TE_kb := nearest_TE_bp / 1000]
  
  out[, left_TE_dist_kb_plot := fifelse(
    is.na(left_TE_dist_kb) | left_TE_dist_kb <= 0,
    NA_real_,
    left_TE_dist_kb
  )]
  
  out[, right_TE_dist_kb_plot := fifelse(
    is.na(right_TE_dist_kb) | right_TE_dist_kb <= 0,
    NA_real_,
    right_TE_dist_kb
  )]
  
  out[, nearest_TE_kb_plot := fifelse(
    is.na(nearest_TE_kb),
    NA_real_,
    pmax(nearest_TE_kb, plot_floor_kb)
  )]
  
  out
}

# Calculate distances

all_gene_te <- get_flanking_te_dist(
  gene_dt = all_genes,
  te_dt = te_chr9,
  group_name = "All genes"
)

eff_te <- get_flanking_te_dist(
  gene_dt = effectors,
  te_dt = te_chr9,
  group_name = "Candidate effectors"
)

non_eff_te <- all_gene_te[!gene_id %in% eff_te$gene_id]
non_eff_te[, group := "Non-effector genes"]

cat("Non-effector background genes:", nrow(non_eff_te), "\n")

plot_dt <- rbindlist(list(non_eff_te, eff_te), fill = TRUE)

plot_dt[, group := factor(
  group,
  levels = c("Non-effector genes", "Candidate effectors")
)]


# Statistics of TE overlap enrichment


overlap_table <- table(plot_dt$group, plot_dt$overlaps_TE)
print(overlap_table)

fisher_overlap <- fisher.test(overlap_table)

overlap_summary <- plot_dt[, .(
  n_genes = .N,
  n_overlap_TE = sum(overlaps_TE, na.rm = TRUE),
  n_non_overlap_TE = sum(!overlaps_TE, na.rm = TRUE),
  prop_overlap_TE = mean(overlaps_TE, na.rm = TRUE)
), by = group]

print(overlap_summary)


# Statistics of nearest TE distance including overlaps

eff_dist_all <- eff_te[!is.na(nearest_TE_kb), nearest_TE_kb]
non_eff_dist_all <- non_eff_te[!is.na(nearest_TE_kb), nearest_TE_kb]

wilcox_all_two <- wilcox.test(
  eff_dist_all,
  non_eff_dist_all,
  alternative = "two.sided",
  exact = FALSE
)

wilcox_all_closer <- wilcox.test(
  eff_dist_all,
  non_eff_dist_all,
  alternative = "less",
  exact = FALSE
)

wilcox_all_farther <- wilcox.test(
  eff_dist_all,
  non_eff_dist_all,
  alternative = "greater",
  exact = FALSE
)

median_eff_all <- median(eff_dist_all, na.rm = TRUE)
median_non_eff_all <- median(non_eff_dist_all, na.rm = TRUE)


# Statistics of nearest TE distance excluding overlaps

eff_dist_no_overlap <- eff_te[
  overlaps_TE == FALSE & !is.na(nearest_TE_kb),
  nearest_TE_kb
]

non_eff_dist_no_overlap <- non_eff_te[
  overlaps_TE == FALSE & !is.na(nearest_TE_kb),
  nearest_TE_kb
]

wilcox_no_overlap_two <- wilcox.test(
  eff_dist_no_overlap,
  non_eff_dist_no_overlap,
  alternative = "two.sided",
  exact = FALSE
)

wilcox_no_overlap_closer <- wilcox.test(
  eff_dist_no_overlap,
  non_eff_dist_no_overlap,
  alternative = "less",
  exact = FALSE
)

wilcox_no_overlap_farther <- wilcox.test(
  eff_dist_no_overlap,
  non_eff_dist_no_overlap,
  alternative = "greater",
  exact = FALSE
)

median_eff_no_overlap <- median(eff_dist_no_overlap, na.rm = TRUE)
median_non_eff_no_overlap <- median(non_eff_dist_no_overlap, na.rm = TRUE)


summary_stats <- data.table(
  comparison = "Candidate effectors vs non-effector genes",
  
  n_effectors_all = length(eff_dist_all),
  n_non_effectors_all = length(non_eff_dist_all),
  
  n_effectors_non_overlap = length(eff_dist_no_overlap),
  n_non_effectors_non_overlap = length(non_eff_dist_no_overlap),
  
  effector_overlap_fraction = mean(eff_te$overlaps_TE, na.rm = TRUE),
  non_effector_overlap_fraction = mean(non_eff_te$overlaps_TE, na.rm = TRUE),
  fisher_overlap_p_value = fisher_overlap$p.value,
  fisher_overlap_odds_ratio = unname(fisher_overlap$estimate),
  
  median_effector_distance_all_kb = median_eff_all,
  median_non_effector_distance_all_kb = median_non_eff_all,
  wilcox_all_two_sided_p = wilcox_all_two$p.value,
  wilcox_all_effectors_closer_p = wilcox_all_closer$p.value,
  wilcox_all_effectors_farther_p = wilcox_all_farther$p.value,
  
  median_effector_distance_non_overlap_kb = median_eff_no_overlap,
  median_non_effector_distance_non_overlap_kb = median_non_eff_no_overlap,
  wilcox_non_overlap_two_sided_p = wilcox_no_overlap_two$p.value,
  wilcox_non_overlap_effectors_closer_p = wilcox_no_overlap_closer$p.value,
  wilcox_non_overlap_effectors_farther_p = wilcox_no_overlap_farther$p.value
)

print(summary_stats)


# left vs right TE distance

plot_2d <- plot_dt[
  !is.na(left_TE_dist_kb_plot) &
    !is.na(right_TE_dist_kb_plot)
]

plot_point_panel <- function(dt, title_text, point_col) {
  
  ggplot(dt, aes(x = left_TE_dist_kb_plot, y = right_TE_dist_kb_plot)) +
    geom_point(
      colour = point_col,
      size = 1.5,
      alpha = 0.55
    ) +
    scale_x_log10(
      limits = c(0.001, 1000),
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
      labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
    ) +
    scale_y_log10(
      limits = c(0.001, 1000),
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
      labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = title_text,
      x = "Distance to nearest left-side TE (kb)",
      y = "Distance to nearest right-side TE (kb)"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

p_all <- plot_point_panel(
  plot_2d[group == "Non-effector genes"],
  "Non-effector genes",
  gene_col
)

p_eff <- plot_point_panel(
  plot_2d[group == "Candidate effectors"],
  "Candidate effectors",
  eff_col
)

p_AB <- p_all + p_eff + plot_annotation(tag_levels = "A")

# TE overlap proportion

p_overlap <- ggplot(
  overlap_summary,
  aes(x = group, y = prop_overlap_TE, fill = group)
) +
  geom_col(width = 0.6, alpha = 0.8, colour = "black", linewidth = 0.3) +
  geom_text(
    aes(
      label = paste0(
        n_overlap_TE,
        "/",
        n_genes,
        "\n",
        round(prop_overlap_TE * 100, 1),
        "%"
      )
    ),
    vjust = -0.35,
    size = 3.8
  ) +
  scale_fill_manual(values = c(
    "Non-effector genes" = gene_col,
    "Candidate effectors" = eff_col
  )) +
  scale_y_continuous(
    limits = c(0, max(overlap_summary$prop_overlap_TE, na.rm = TRUE) * 1.25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Genes overlapping annotated TEs",
    subtitle = paste0(
      "Fisher exact p = ",
      formatC(fisher_overlap$p.value, format = "e", digits = 2)
    ),
    x = NULL,
    y = "Proportion overlapping TE"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

# nearest TE distance, all genes included

box_all_dt <- plot_dt[!is.na(nearest_TE_kb_plot)]

label_all <- paste0(
  "All genes included",
  "\nMedian candidate effector = ",
  round(median_eff_all, 3),
  " kb",
  "\nMedian non-effector = ",
  round(median_non_eff_all, 3),
  " kb",
  "\nTwo-sided Wilcoxon p = ",
  formatC(wilcox_all_two$p.value, format = "e", digits = 2)
)

p_box_all <- ggplot(
  box_all_dt,
  aes(x = group, y = nearest_TE_kb_plot, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    alpha = 0.75,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(shape = overlaps_TE),
    width = 0.18,
    size = 1.1,
    alpha = 0.35,
    colour = "black"
  ) +
  scale_y_log10(
    limits = c(plot_floor_kb, 1000),
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c("0/0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  scale_fill_manual(values = c(
    "Non-effector genes" = gene_col,
    "Candidate effectors" = eff_col
  )) +
  scale_shape_manual(values = c(
    "FALSE" = 16,
    "TRUE" = 1
  )) +
  annotate(
    "text",
    x = 1.5,
    y = 80,
    label = label_all,
    size = 3.5
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Nearest TE distance including TE-overlapping genes",
    x = NULL,
    y = "Distance to nearest TE (kb, log10)",
    fill = "Group",
    shape = "Overlaps TE"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

# nearest TE distance, non-overlapping genes only

box_no_overlap_dt <- plot_dt[
  overlaps_TE == FALSE &
    !is.na(nearest_TE_kb)
]

label_no_overlap <- paste0(
  "TE-overlapping genes excluded",
  "\nMedian candidate effector = ",
  round(median_eff_no_overlap, 3),
  " kb",
  "\nMedian non-effector = ",
  round(median_non_eff_no_overlap, 3),
  " kb",
  "\nTwo-sided Wilcoxon p = ",
  formatC(wilcox_no_overlap_two$p.value, format = "e", digits = 2)
)

p_box_no_overlap <- ggplot(
  box_no_overlap_dt,
  aes(x = group, y = nearest_TE_kb, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    alpha = 0.75,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.18,
    size = 1.1,
    alpha = 0.35,
    colour = "black"
  ) +
  scale_y_log10(
    limits = c(plot_floor_kb, 1000),
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  scale_fill_manual(values = c(
    "Non-effector genes" = gene_col,
    "Candidate effectors" = eff_col
  )) +
  annotate(
    "text",
    x = 1.5,
    y = 80,
    label = label_no_overlap,
    size = 3.5
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Nearest TE distance among non-overlapping genes",
    x = NULL,
    y = "Distance to nearest TE (kb, log10)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

combined <- (p_all + p_eff) / (p_overlap + p_box_no_overlap) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))


##############Chr9A effector expression analysis

  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)


setwd("~/Desktop/BLR_data/New_518_560/")

effector_list_file <- "560A_chr9_effector_list.txt"
effector_bed_file  <- "effector_genes_chr9A_by_suffix.bed"
expr_file          <- "Ph560A_DESeq2_normalized_counts.tsv"
te_file            <- "560A_TE_all_chr13Aupdated.txt"

out_prefix <- "Chr9A_effector_expression_isolate_level"

old_col <- "#4D6FA3"
new_col <- "#C7646A"

group_cols <- c(
  "older isolates" = old_col,
  "CNV-gained isolates" = new_col
)

eff_list <- fread(effector_list_file)
setnames(eff_list, 1, "gene_id")
eff_list <- unique(eff_list[, .(gene_id)])

eff_bed <- fread(effector_bed_file, header = FALSE, fill = TRUE)

eff_bed <- eff_bed[, .(
  chr = V1,
  start = suppressWarnings(as.integer(V2)),
  end = suppressWarnings(as.integer(V3)),
  gene_id = V4
)]

eff_bed <- eff_bed[
  chr == "Chr9A" &
    !is.na(start) &
    !is.na(end) &
    start < end
]

te <- fread(te_file, header = FALSE, fill = TRUE)

te <- te[, .(
  chr = V1,
  start = suppressWarnings(as.integer(V2)),
  end = suppressWarnings(as.integer(V3)),
  TE_class = V4
)]

te <- te[
  chr == "Chr9A" &
    !is.na(start) &
    !is.na(end) &
    start < end
]

expr <- fread(expr_file)
setnames(expr, 1, "gene_id")

cat("Effector IDs in list:", nrow(eff_list), "\n")
cat("Chr9A effector BED rows:", nrow(eff_bed), "\n")
cat("Chr9A TE features:", nrow(te), "\n")

# Chr9A effectors only

eff <- eff_bed %>%
  filter(gene_id %in% eff_list$gene_id) %>%
  distinct(gene_id, .keep_all = TRUE)

cat("Unique Chr9A candidate effectors:", nrow(eff), "\n")

# Nearest TE distance

eff_gr <- GRanges(
  seqnames = eff$chr,
  ranges = IRanges(start = eff$start, end = eff$end),
  gene_id = eff$gene_id
)

te_gr <- GRanges(
  seqnames = te$chr,
  ranges = IRanges(start = te$start, end = te$end),
  TE_class = te$TE_class
)

nearest_te <- distanceToNearest(
  eff_gr,
  te_gr,
  ignore.strand = TRUE
)

te_dist_raw <- data.table(
  gene_id = mcols(eff_gr)$gene_id[queryHits(nearest_te)],
  nearest_TE_class = mcols(te_gr)$TE_class[subjectHits(nearest_te)],
  dist_to_TE = mcols(nearest_te)$distance
)

te_dist <- te_dist_raw[
  order(gene_id, dist_to_TE)
][
  , .SD[1], by = gene_id
]

eff_anno <- as.data.table(eff) %>%
  left_join(te_dist, by = "gene_id") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(
    TE_border = case_when(
      is.na(dist_to_TE) ~ "no TE found",
      dist_to_TE == 0 ~ "overlapping TE",
      dist_to_TE <= 2000 ~ "<=2 kb from TE",
      dist_to_TE <= 5000 ~ "2-5 kb from TE",
      TRUE ~ ">5 kb from TE"
    ),
    TE_border = factor(
      TE_border,
      levels = c(
        "overlapping TE",
        "<=2 kb from TE",
        "2-5 kb from TE",
        ">5 kb from TE",
        "no TE found"
      )
    ),
    mid = (start + end) / 2
  )

cat("Annotated effectors:", nrow(eff_anno), "\n")

fwrite(
  eff_anno,
  paste0(out_prefix, "_effector_TE_annotation.tsv"),
  sep = "\t"
)


# Expression table

expr_eff <- expr %>%
  filter(gene_id %in% eff_anno$gene_id) %>%
  left_join(eff_anno, by = "gene_id")

expr_cols <- grep("^Ph[0-9]+_Rep[0-9]+$", colnames(expr_eff), value = TRUE)

cat("Expression columns detected:\n")
print(expr_cols)

expr_long <- expr_eff %>%
  pivot_longer(
    cols = all_of(expr_cols),
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  mutate(
    isolate = sub("_Rep[0-9]+$", "", sample),
    replicate = sub("^Ph[0-9]+_", "", sample),
    isolate_group = case_when(
      isolate %in% c("Ph491", "Ph518") ~ "older isolates",
      isolate %in% c("Ph612", "Ph685") ~ "CNV-gained isolates",
      TRUE ~ NA_character_
    ),
    log2_expr = log2(norm_count + 1)
  ) %>%
  filter(!is.na(isolate_group))

cat("Unique effectors in expression table:", length(unique(expr_long$gene_id)), "\n")
cat("Replicate-level expression rows:", nrow(expr_long), "\n")

# Average replicates

expr_gene_isolate <- expr_long %>%
  group_by(
    gene_id, chr, start, end, mid,
    TE_border, dist_to_TE, nearest_TE_class,
    isolate, isolate_group
  ) %>%
  summarise(
    mean_log2_expr = mean(log2_expr, na.rm = TRUE),
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  mutate(
    isolate = factor(
      isolate,
      levels = c("Ph491", "Ph518", "Ph612", "Ph685")
    ),
    isolate_group = factor(
      isolate_group,
      levels = c("older isolates", "CNV-gained isolates")
    ),
    isolate_label = case_when(
      isolate == "Ph491" ~ "Ph491\n(older)",
      isolate == "Ph518" ~ "Ph518\n(older)",
      isolate == "Ph612" ~ "Ph612\n(new)",
      isolate == "Ph685" ~ "Ph685\n(new)"
    ),
    isolate_label = factor(
      isolate_label,
      levels = c(
        "Ph491\n(older)",
        "Ph518\n(older)",
        "Ph612\n(new)",
        "Ph685\n(new)"
      )
    )
  )

cat("Gene-isolate rows:", nrow(expr_gene_isolate), "\n")
cat("Expected maximum rows:", length(unique(expr_gene_isolate$gene_id)) * 4, "\n")


# Statistics accounting for gene repeated measures

# Gene fixed-effect model:
# tests whether expression differs between older and CNV-gained isolates
# while accounting for baseline expression differences among genes.

lm_group <- lm(
  mean_log2_expr ~ gene_id + isolate_group,
  data = expr_gene_isolate
)

lm_group_summary <- summary(lm_group)

group_coef_name <- grep(
  "^isolate_group",
  rownames(lm_group_summary$coefficients),
  value = TRUE
)

group_effect <- lm_group_summary$coefficients[group_coef_name, "Estimate"]
group_p <- lm_group_summary$coefficients[group_coef_name, "Pr(>|t|)"]

# Isolate-specific effect:
# tests whether the four isolates differ after accounting for gene identity.

lm_isolate <- lm(
  mean_log2_expr ~ gene_id + isolate,
  data = expr_gene_isolate
)

anova_isolate <- anova(lm_isolate)
isolate_p <- anova_isolate["isolate", "Pr(>F)"]

# TE interaction:
# tests whether old/new group difference depends on TE-distance class.

lm_no_interaction <- lm(
  mean_log2_expr ~ gene_id + isolate_group,
  data = expr_gene_isolate
)

lm_te_interaction <- lm(
  mean_log2_expr ~ gene_id + isolate_group + isolate_group:TE_border,
  data = expr_gene_isolate
)

anova_te_interaction <- anova(lm_no_interaction, lm_te_interaction)
te_interaction_p <- anova_te_interaction$`Pr(>F)`[2]

# Descriptive statistics

expr_group_stats <- expr_gene_isolate %>%
  group_by(isolate, isolate_group) %>%
  summarise(
    n_genes = n(),
    median_log2_expr = median(mean_log2_expr, na.rm = TRUE),
    mean_log2_expr = mean(mean_log2_expr, na.rm = TRUE),
    .groups = "drop"
  )

te_border_stats <- expr_gene_isolate %>%
  group_by(TE_border, isolate_group) %>%
  summarise(
    n_values = n(),
    n_genes = n_distinct(gene_id),
    median_log2_expr = median(mean_log2_expr, na.rm = TRUE),
    mean_log2_expr = mean(mean_log2_expr, na.rm = TRUE),
    .groups = "drop"
  )

model_stats <- data.table(
  n_genes = length(unique(expr_gene_isolate$gene_id)),
  n_gene_isolate_values = nrow(expr_gene_isolate),
  group_effect_CNV_gained_vs_older = group_effect,
  gene_fixed_effect_group_p = group_p,
  gene_fixed_effect_isolate_p = isolate_p,
  TE_border_by_group_interaction_p = te_interaction_p
)


# Heatmap of isolate-level expression

heatmap_dt <- expr_gene_isolate %>%
  select(gene_id, isolate, mean_log2_expr) %>%
  pivot_wider(
    names_from = isolate,
    values_from = mean_log2_expr
  )

heatmap_mat <- as.matrix(heatmap_dt[, c("Ph491", "Ph518", "Ph612", "Ph685")])
rownames(heatmap_mat) <- heatmap_dt$gene_id

heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

annotation_col <- data.frame(
  Group = factor(
    c("older isolates", "older isolates", "CNV-gained isolates", "CNV-gained isolates"),
    levels = c("older isolates", "CNV-gained isolates")
  )
)

rownames(annotation_col) <- c("Ph491", "Ph518", "Ph612", "Ph685")

ann_colors <- list(
  Group = group_cols
)

dev.off()

g = pheatmap(
  heatmap_mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 5,
  fontsize_col = 11,
  main = "TE-proximal Chr9A effector expression",
  filename = paste0(out_prefix, "_heatmap_scaled_expression.pdf"),
  width = 8,
  height = 10
)

g2 = pheatmap(
  heatmap_mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 5,
  fontsize_col = 11,
  main = "TE-proximal Chr9A effector expression",
  filename = paste0(out_prefix, "_heatmap_scaled_expression.png"),
  width = 8,
  height = 10
)


# top  old-vs-new contrasted effectors

heatmap_dt <- expr_gene_isolate %>%
  select(gene_id, isolate, mean_log2_expr) %>%
  pivot_wider(
    names_from = isolate,
    values_from = mean_log2_expr
  )

heatmap_dt <- heatmap_dt %>%
  filter(
    !is.na(Ph491),
    !is.na(Ph518),
    !is.na(Ph612),
    !is.na(Ph685)
  ) %>%
  mutate(
    older_mean = rowMeans(cbind(Ph491, Ph518), na.rm = TRUE),
    new_mean = rowMeans(cbind(Ph612, Ph685), na.rm = TRUE),
    old_new_delta = older_mean - new_mean,
    abs_old_new_delta = abs(old_new_delta)
  ) %>%
  arrange(desc(abs_old_new_delta)) %>%
  slice_head(n = 35)

cat("Top contrasted genes in heatmap:", nrow(heatmap_dt), "\n")

heatmap_mat <- as.matrix(
  heatmap_dt[, c("Ph491", "Ph518", "Ph612", "Ph685")]
)

rownames(heatmap_mat) <- heatmap_dt$gene_id

# Row-scale expression
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

annotation_col <- data.frame(
  Group = factor(
    c(
      "older isolates",
      "older isolates",
      "CNV-gained isolates",
      "CNV-gained isolates"
    ),
    levels = c("older isolates", "CNV-gained isolates")
  )
)

rownames(annotation_col) <- c("Ph491", "Ph518", "Ph612", "Ph685")

ann_colors <- list(
  Group = group_cols
)

dev.off ()

my_colors <- colorRampPalette(c("#7B8DBF", "white", "#57B893"))(100)

g3=pheatmap(
  heatmap_mat_scaled,
  color = my_colors,  
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  fontsize_col = 11,
  main = "Top old-new contrasted Chr9A effector genes",
  filename = paste0(out_prefix, "_heatmap_old_new_contrasted_scaled.pdf"),
  width = 7,
  height = 8
)


# isolate-level expression across Chr9A
# no collapsed old/new delta

label_A <- paste0(
  "Gene-adjusted isolate effect p = ",
  formatC(isolate_p, format = "e", digits = 2)
)

pA <- ggplot(
  expr_gene_isolate,
  aes(x = mid / 1e6, y = mean_log2_expr)
) +
  geom_point(
    aes(colour = isolate_group),
    size = 1.8,
    alpha = 0.75
  ) +
  facet_wrap(
    ~ isolate_label,
    nrow = 1
  ) +
  scale_colour_manual(values = group_cols) +
  annotate(
    "text",
    x = min(expr_gene_isolate$mid / 1e6, na.rm = TRUE),
    y = max(expr_gene_isolate$mean_log2_expr, na.rm = TRUE),
    label = label_A,
    hjust = 0,
    vjust = 1,
    size = 3.4
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Chr9A position (Mb)",
    y = "Mean expression\nlog2(count + 1)",
    title = "A  Chr9A effector expression varies among isolates",
    colour = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

# expression by isolate
# one point = one gene averaged across replicates in one isolate

label_B <- paste0(
  "Gene-adjusted group p = ",
  formatC(group_p, format = "e", digits = 2)
)

pB <- ggplot(
  expr_gene_isolate,
  aes(x = isolate_label, y = mean_log2_expr, fill = isolate_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.58,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.16,
    size = 1.2,
    alpha = 0.45,
    colour = "black"
  ) +
  annotate(
    "text",
    x = 2.5,
    y = max(expr_gene_isolate$mean_log2_expr, na.rm = TRUE),
    label = label_B,
    hjust = 0.5,
    vjust = 1,
    size = 3.6
  ) +
  scale_fill_manual(values = group_cols) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Mean DESeq2-normalized expression\nlog2(count + 1)",
    title = "Expression of Chr9A effectors in each isolate"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10)
  )

# expression by TE-distance class and group
# no collapsed old-new delta

plot_te_dt <- expr_gene_isolate %>%
  filter(!is.na(TE_border)) %>%
  droplevels()

label_C <- paste0(
  "Gene-adjusted TE class x group p = ",
  formatC(te_interaction_p, format = "e", digits = 2)
)

pC <- ggplot(
  plot_te_dt,
  aes(x = TE_border, y = mean_log2_expr, fill = isolate_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.65,
    position = position_dodge(width = 0.75),
    linewidth = 0.5
  ) +
  geom_point(
    aes(colour = isolate_group),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1.1,
    alpha = 0.45
  ) +
  annotate(
    "text",
    x = 1,
    y = max(plot_te_dt$mean_log2_expr, na.rm = TRUE),
    label = label_C,
    hjust = 0,
    vjust = 1,
    size = 3.6
  ) +
  scale_fill_manual(values = group_cols) +
  scale_colour_manual(values = group_cols) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Mean expression\nlog2(count + 1)",
    title = "Effector expression by TE-distance class",
    fill = "Group",
    colour = "Group"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

combined_ABC <- pA / (pB + pC) +
  plot_layout(heights = c(0.8, 1.1))


