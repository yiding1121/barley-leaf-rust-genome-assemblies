library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)

setwd("~/Desktop/BLR_data/New_518_560/")
# ==== Load chromosome sizes ====

# ==== Load and count genes ====
gff <- read_tsv("Ph560A_braker.EVM.gff3", comment = "#", col_names = FALSE)
gff <- read_tsv("Ph560B_braker.EVM.gff3", comment = "#", col_names = FALSE)

colnames(gff)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
genes <- gff %>% filter(type == "gene")

gene_counts <- genes %>%
  group_by(seqid) %>%
  summarise(gene_count = n(), .groups = "drop")

# ==== Load and summarize TE ====

# ==== Load and fix chromosome sizes ====

# ==== Load TE data correctly ====
te <- read.table("560A_TE_all_chr13Aupdated.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
te <- read.table("560B_TE_all.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(te) <- c("chrom", "start", "end", "type")
te$start <- as.numeric(te$start)
te$end <- as.numeric(te$end)
te <- te[!is.na(te$start) & !is.na(te$end), ]

# ==== Create GRanges and merge overlapping regions ====
te <- te[te$end >= te$start - 1, ]
te_gr <- GRanges(seqnames = te$chrom, ranges = IRanges(start = te$start, end = te$end))
te_gr_merged <- reduce(te_gr)  # ✅ this step removes overlapping bases

# ==== Calculate total TE length per chromosome ====
te_chr_length <- as.data.frame(te_gr_merged) %>%
  group_by(seqnames) %>%
  summarise(total_te_bp = sum(width)) %>%
  rename(chrom = seqnames)

# ==== Calculate % TE density ====

chrom_sizes <- read.table("blr560A_new.chrom.sizes", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chrom_sizes <- read.table("blr560B.chrom.sizes", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

te_chr_df <- te_chr_length %>%
  mutate(chrom = as.character(chrom)) %>%
  left_join(chrom_sizes %>% rename(length = length), by = "chrom") %>%
  mutate(te_density = (total_te_bp / length) * 100)

te_chr_df <- te_chr_length %>%
  mutate(chrom = as.character(chrom)) %>%
  left_join(chrom_sizes %>% rename(length = max_pos), by = "chrom") %>%
  mutate(te_density = (total_te_bp / length) * 100)

unique(chrom_sizes$chrom)

# Only keep chrom and te_density from te_chr_df
te_chr_df_clean <- te_chr_df %>% select(chrom, te_density)

# ==== Load and compute GC content ====
# ==== Load GC content table (100 kb binned) ====
gc_bins <- read.table("blr560A_gc_cov_100kb.tsv", header = TRUE, sep = "\t")  # update filename if needed
gc_bins <- read.table("blr560B_gc_cov_100kb.tsv", header = TRUE, sep = "\t")  # update filename if needed

# Compute average GC content per chromosome
gc_content <- gc_bins %>%
  group_by(chrom) %>%
  summarise(gc = mean(GC, na.rm = TRUE), .groups = "drop")

# ==== Load E1 values ====
# Load E1 file and ensure E1 is numeric
e1 <- read_tsv("./HiC/560A_cis_eigs_100kb.tsv", col_names = TRUE)
e1 <- read_tsv("./HiC/560B_cis_eigs_100kb.tsv", col_names = TRUE)

# Convert E1 to numeric (this handles any accidental character entries)
e1$E1 <- as.numeric(e1$E1)

e1_avg <- e1 %>%
  group_by(chrom) %>%
  summarise(E1_mean = mean(E1, na.rm = TRUE), .groups = "drop")

# ==== Combine all ====
# ==== Load chromosome sizes correctly once ====
#chrom_sizes <- read.table("blr560A_new.chrom.sizes", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("chrom", "length")
chrom_sizes$length <- as.numeric(chrom_sizes$length)
chrom_sizes$length_Mb <- chrom_sizes$length / 1e6

head(chrom_sizes)

# ==== Combine all features into summary ===

# Clearly join tables now (no conflict):
summary_df <- chrom_sizes %>%
  left_join(gene_counts, by = c("chrom" = "seqid")) %>%
  left_join(te_chr_df_clean, by = "chrom") %>%
  left_join(gc_content, by = "chrom") %>%
  left_join(e1_avg, by = "chrom") %>%
  mutate(
    gene_density = gene_count / length_Mb
  )
# Pivot to long format
summary_long <- summary_df %>%
  select(chrom, gene_density, te_density, gc, E1_mean) %>%
  pivot_longer(-chrom, names_to = "feature", values_to = "value")

# Faceted plot
ggplot(summary_long, aes(x = chrom, y = value, fill = feature)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~feature, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Genome Feature Distributions per Chromosome (HapA)",
    x = "Chromosome",
    y = "Feature Value"
  ) 

# === Reorder chromosomes: put Chr9A last or highlight in the middle ===
library(forcats)
# Reorder chromosomes and highlight Chr9A
summary_long <- summary_df %>%
  select(chrom, gene_density, te_density, gc, E1_mean) %>%
  pivot_longer(-chrom, names_to = "feature", values_to = "value") %>%
  mutate(
    chrom = fct_relevel(chrom, c(setdiff(sort(unique(summary_df$chrom)), "Chr9A"), "Chr9A")),
    highlight = ifelse(chrom == "Chr9A", "highlight", "normal")
  )

# === Define custom fill colors ===
feature_colors <- c(
  E1_mean = "#F8766D",      # coral/red
  gc = "#7CAE00",           # green
  gene_density = "#00BFC4", # teal
  te_density = "#C77CFF"    # purple
)

# === Plot ===
ggplot(summary_long, aes(x = chrom, y = value, fill = highlight)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(value, 1)), vjust = -0.3, size = 3) +
  facet_wrap(~feature, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("highlight" = "#7B8DBE", "normal" = "#57B899")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 14)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds top space
  labs(
    title = "Genome Feature Distributions per Chromosome (HapA)",
    x = "Chromosome",
    y = "Feature Value"
  )
#########################circos plot
#######################################
library(rtracklayer)

genes_gr <- import("Ph560A_braker.EVM.chr13Aupdated.gff3")
genes_gr <- import("Ph560B_braker.EVM.gff3")

genes_gr <- genes_gr[genes_gr$type == "gene"]

table(strand(genes_gr))  # Confirm you now get '+' and '-' values
# Subset genes by strand
gene_plus  <- genes_gr[strand(genes_gr) == "+"]
gene_minus <- genes_gr[strand(genes_gr) == "-"]

# 5' ends: start for +, end for -
gene_5prime <- c(
  resize(gene_plus, width = 1, fix = "start"),
  resize(gene_minus, width = 1, fix = "end")
)

# 3' ends: end for +, start for -
gene_3prime <- c(
  resize(gene_plus, width = 1, fix = "end"),
  resize(gene_minus, width = 1, fix = "start")
)

# Ensure order matches original genes_gr
gene_5prime <- gene_5prime[order(c(which(strand(genes_gr) == "+"), which(strand(genes_gr) == "-")))]
gene_3prime <- gene_3prime[order(c(which(strand(genes_gr) == "+"), which(strand(genes_gr) == "-")))]

dist_5prime <- distanceToNearest(gene_5prime, te_gr_merged)
dist_3prime <- distanceToNearest(gene_3prime, te_gr_merged)

min_dist <- pmin(mcols(dist_5prime)$distance, mcols(dist_3prime)$distance)
######################################

# Compute distances to nearest TE for 5' and 3' ends
dist_5prime <- distanceToNearest(gene_5prime, te_gr_merged)
dist_3prime <- distanceToNearest(gene_3prime, te_gr_merged)

distances_5prime_kb <- mcols(dist_5prime)$distance / 1000
distances_3prime_kb <- mcols(dist_3prime)$distance / 1000

# Gene metadata
genes_gr <- genes_gr[!grepl("^unscaffold", as.character(seqnames(genes_gr)))]
genes_gr <- dropSeqlevels(genes_gr, setdiff(seqlevels(genes_gr), seqlevels(bins_gr)), pruning.mode = "coarse")
genes_meta <- as.data.frame(genes_gr)


te_gr_merged
# Add distances in kb
genes_meta$dist_5prime_kb <- distances_5prime_kb
genes_meta$dist_3prime_kb <- distances_3prime_kb
genes_meta$chr_group <- ifelse(genes_meta$seqnames == "Chr9A", "Chr9A", "Other")

# Tidy format
library(tidyr)

genes_long_kb <- genes_meta %>%
  select(seqnames, chr_group, dist_5prime_kb, dist_3prime_kb) %>%
  pivot_longer(cols = c(dist_5prime_kb, dist_3prime_kb),
               names_to = "end", values_to = "distance_kb") %>%
  mutate(end = ifelse(end == "dist_5prime_kb", "5′ end", "3′ end"))


ggplot(genes_long_kb, aes(x = end, y = distance_kb, fill = chr_group)) +
  geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Chr9B" = "#FF6F61", "Other" = "gray70")) +
  facet_wrap(~chr_group, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distance from 5′ and 3′ Gene Ends to Nearest TE (in kb)",
    subtitle = "Chr9AB vs Other Chromosomes",
    x = "Gene End",
    y = "Distance to Nearest TE (kb)"
  )


wilcox.test(distance_kb ~ chr_group, data = genes_long_kb %>% filter(end == "5′ end"))
wilcox.test(distance_kb ~ chr_group, data = genes_long_kb %>% filter(end == "3′ end"))

# Compute summary stats: mean per group
summary_stats <- genes_long_kb %>%
  group_by(chr_group, end) %>%
  summarise(mean_kb = mean(distance_kb, na.rm = TRUE),
            median_kb = median(distance_kb, na.rm = TRUE),
            .groups = "drop")

# Plot
ggplot(genes_long_kb, aes(x = end, y = distance_kb, fill = chr_group)) +
  geom_violin(trim = FALSE, alpha = 0.6, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, alpha = 0.5, position = position_dodge(width = 0.9), outlier.shape = NA) +
  # Add mean points
  geom_point(data = summary_stats,
             aes(x = end, y = mean_kb, group = chr_group),
             position = position_dodge(width = 0.9),
             shape = 18, size = 3, color = "black") +
  # Optionally add median lines
  geom_segment(data = summary_stats,
               aes(x = as.numeric(as.factor(end)) - 0.25 + (chr_group == "Chr9B") * 0.5,
                   xend = as.numeric(as.factor(end)) - 0.25 + (chr_group == "Chr9B") * 0.5,
                   y = median_kb, yend = median_kb),
               color = "black", size = 0.6) +
  scale_fill_manual(values = c("Chr9B" = "#FF6F61", "Other" = "gray70")) +
  facet_wrap(~chr_group, nrow = 1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  labs(
    title = "Distance from 5′ and 3′ Gene Ends to Nearest TE (kb)",
    subtitle = "Means (dots) and medians (lines) per group — highlighting Chr9A",
    x = "Gene End",
    y = "Distance to Nearest TE (kb)"
  )

################################
ggplot(genes_df, aes(x = chr_group, y = log10(min_end_TE_distance + 1), fill = chr_group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("Chr9B" = "#FF6F61", "Other" = "gray70")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distance from Gene Ends to Nearest TE",
    subtitle = "Chr9B genes tend to be closer to TEs at 5'/3' ends",
    x = "Chromosome Group",
    y = "log10(Min(5',3') TE Distance + 1)"
  )

##################
# Filter to ≤ 3 kb
genes_filtered <- genes_long_kb %>%
  filter(distance_kb <= 3)

# Compute mean distances for annotations
summary_stats <- genes_filtered %>%
  group_by(chr_group, end) %>%
  summarise(
    mean = mean(distance_kb, na.rm = TRUE),
    .groups = "drop"
  )

# Determine y position for text (relative to max visible range)
summary_stats$y_pos <- 3 * 0.95

# Plot
ggplot(genes_long_kb, aes(x = end, y = distance_kb, fill = chr_group)) +
  geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, alpha = 0.4, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "black",
               position = position_dodge(width = 0.9)) +
  geom_text(data = summary_stats,
            aes(x = end, y = y_pos, label = paste0("Mean: ", sprintf("%.2f", mean))),
            position = position_dodge(width = 0.9),
            size = 3.5,
            inherit.aes = FALSE) +
  scale_fill_manual(values = c("Chr9A" = "#7B8DBE", "Other" = "#57B899")) +
  facet_wrap(~ chr_group, nrow = 1) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distance from 5′ and 3′ Gene Ends to Nearest TE",
    subtitle = "Chr9A vs Other Chromosomes (Distances ≤ 3 kb)",
    x = "Gene End",
    y = "Distance to Nearest TE (kb)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )
##################
# Calculate means for annotation
summary_stats <- genes_long_kb %>%
  group_by(chr_group, end) %>%
  summarise(
    mean = mean(distance_kb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_pos = max(genes_long_kb$distance_kb, na.rm = TRUE) * 0.9)

# Plot all distances
ggplot(genes_long_kb, aes(x = end, y = distance_kb, fill = chr_group)) +
  geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, alpha = 0.4, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "black",
               position = position_dodge(width = 0.9)) +
  geom_text(data = summary_stats,
            aes(x = end, y = y_pos, label = paste0("Mean: ", sprintf("%.2f", mean))),
            size = 3.5, vjust = -0.7, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Chr9A" = "#7B8DBE", "Other" = "#57B899")) +
  facet_wrap(~ chr_group, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distance from 5′ and 3′ Gene Ends to Nearest TE",
    subtitle = "Chr9A vs Other Chromosomes",
    x = "Gene End",
    y = "Distance to Nearest TE (kb)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )



# === Plot TE % coverage in intergenic regions ===
ggplot(intergenic_df, aes(x=chr_group, y=te_overlap_pct, fill=chr_group)) +
  geom_violin(trim=FALSE, alpha=0.7) +
  geom_boxplot(width=0.1, alpha=0.5) +
  theme_minimal(base_size=14) +
  scale_fill_manual(values=c("Chr9A"="#FF6F61", "Other"="gray70")) +
  labs(
    title="% TE Coverage in Intergenic Regions",
    subtitle="Intergenic regions on Chr9B have higher TE coverage",
    x="Chromosome Group",
    y="TE Coverage (%)"
  )

wilcox.test(width ~ chr_group, data = intergenic_df)


# Median intergenic length per group
intergenic_df %>%
  group_by(chr_group) %>%
  summarise(median_intergenic = median(intergenic_length))

# Effect size as difference in medians
diff_median <- intergenic_df %>%
  group_by(chr_group) %>%
  summarise(median_intergenic = median(intergenic_length)) %>%
  summarise(effect_size = diff(median_intergenic))

diff_median


# Compute intergenic regions as gaps between sorted genes
genes_sorted <- sort(genes_gr)
intergenic_gr <- gaps(genes_sorted)

# Filter to usable intergenic regions (exclude ends and unscaffolded chromosomes)
intergenic_df <- as.data.frame(intergenic_gr) %>%
  filter(as.character(seqnames) %in% chrom_sizes$chrom) %>%
  mutate(chr_group = ifelse(seqnames == "Chr9A", "Chr9A", "Other"))


# Find overlaps between intergenic regions and TEs
overlaps <- findOverlaps(intergenic_gr, te_gr_merged)
intersect_ranges <- pintersect(intergenic_gr[queryHits(overlaps)], te_gr_merged[subjectHits(overlaps)])

# Total TE coverage per intergenic region
te_overlap_bp <- tapply(width(intersect_ranges), queryHits(overlaps), sum)

# Assign TE basepair overlap to each region
intergenic_df$te_overlap_bp <- 0
intergenic_df$te_overlap_bp[as.numeric(names(te_overlap_bp))] <- te_overlap_bp

intergenic_df <- intergenic_df %>%
  rename(intergenic_length = width)

intergenic_df$te_overlap_pct <- (intergenic_df$te_overlap_bp / intergenic_df$intergenic_length) * 100


ggplot(intergenic_df, aes(x = chr_group, y = intergenic_length / 1000, fill = chr_group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("Chr9A" = "#FF6F61", "Other" = "gray70")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Intergenic Region Lengths",
    x = "Chromosome Group",
    y = "Intergenic Region Length (kb)"
  )

ggplot(intergenic_df, aes(x=chr_group, y=te_overlap_pct, fill=chr_group)) +
  geom_violin(trim=FALSE, alpha=0.7) +
  geom_boxplot(width=0.1, alpha=0.5) +
  theme_minimal(base_size=14) +
  scale_fill_manual(values=c("Chr9A"="#FF6F61", "Other"="gray70")) +
  labs(
    title="% TE Coverage in Intergenic Regions",
    x="Chromosome Group",
    y="TE Coverage (%)"
  )

wilcox.test(te_overlap_pct ~ chr_group, data = intergenic_df)


intergenic_df %>%
  group_by(chr_group) %>%
  summarise(
    n = n(),
    mean_kb = mean(intergenic_length, na.rm = TRUE) / 1000,
    median_kb = median(intergenic_length, na.rm = TRUE) / 1000,
    sd_kb = sd(intergenic_length, na.rm = TRUE) / 1000,
    .groups = "drop"
  )

intergenic_df %>%
  group_by(chr_group) %>%
  summarise(
    n = n(),
    mean_TE_pct = mean(te_overlap_pct, na.rm = TRUE),
    median_TE_pct = median(te_overlap_pct, na.rm = TRUE),
    sd_TE_pct = sd(te_overlap_pct, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(intergenic_df, aes(x = chr_group, y = intergenic_length / 1000, fill = chr_group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("Chr9B" = "#FF6F61", "Other" = "gray70")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Intergenic Region Length Comparison",
    subtitle = "Chr9A has longer intergenic regions",
    x = "Chromosome Group",
    y = "Intergenic Length (kb)"
  )

ggplot(intergenic_df, aes(x = intergenic_length / 1000, y = te_overlap_pct, color = chr_group)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("Chr9B" = "#FF6F61", "Other" = "gray70")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TE Coverage vs. Intergenic Region Size",
    x = "Intergenic Region Length (kb)",
    y = "TE Coverage (%)"
  )

# Total intergenic basepairs & total TE-covered basepairs
intergenic_df %>%
  group_by(chr_group) %>%
  summarise(
    total_intergenic_bp = sum(intergenic_length),
    total_te_bp = sum(te_overlap_bp, na.rm = TRUE),
    percent_te_occupied = total_te_bp / total_intergenic_bp * 100,
    .groups = "drop"
  )

####################compartments and TE
# Use existing GC or E1 bins as your reference bin set
bins_df <- read.table("blr560A_gc_cov_100kb.tsv", header = TRUE, sep = "\t")
bins_gr <- GRanges(seqnames = bins_df$chrom, ranges = IRanges(start = bins_df$start + 1, end = bins_df$end))

# TE coverage in each bin
overlaps <- findOverlaps(bins_gr, te_gr_merged)
te_intersect <- pintersect(bins_gr[queryHits(overlaps)], te_gr_merged[subjectHits(overlaps)])
te_bp_per_bin <- tapply(width(te_intersect), queryHits(overlaps), sum)

# Build full bin-level TE data frame
bins_df$te_bp <- 0
bins_df$te_bp[as.integer(names(te_bp_per_bin))] <- te_bp_per_bin
bins_df$te_pct <- bins_df$te_bp / 1e5 * 100  # % coverage in each 100kb bin

e1_df <- read.table("HiC/560A_cis_eigs_100kb.tsv", header = TRUE, sep = "\t")
e1_df$E1 <- as.numeric(e1_df$E1)

# Merge by bin position
bin_te_e1 <- bins_df %>%
  inner_join(e1_df, by = c("chrom", "start", "end"))

cor_test <- cor.test(bin_te_e1$E1, bin_te_e1$te_pct, method = "spearman")
print(cor_test)

library(ggplot2)
ggplot(bin_te_e1, aes(x = E1, y = te_pct)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", se = FALSE, color = "darkred") +
  theme_minimal(base_size = 14) +
  labs(
    title = "TE Coverage vs. Hi-C Compartment (E1)",
    x = "E1 Compartment Value",
    y = "TE Coverage per 100kb Bin (%)"
  )

#bin_te_e1$E1_quantile <- ntile(bin_te_e1$E1, 4)

#ggplot(bin_te_e1, aes(x = factor(E1_quantile), y = te_pct)) +
  geom_violin(fill = "lightblue") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TE Coverage Across E1 Quartiles",
    x = "E1 Quartile",
    y = "TE Coverage (%)"
  )


  library(GenomicRanges)
  
  # Step 1: Define your bins as GRanges
  bins_gr <- GRanges(seqnames = bin_te_e1$chrom,
                     ranges = IRanges(start = bin_te_e1$start + 1, end = bin_te_e1$end))  # 1-based
  
  # Step 2: Count number of genes overlapping each bin
  gene_counts <- countOverlaps(bins_gr, genes_gr)  # genes_gr should be all genes as GRanges
  
  # Step 3: Assign to your bin-level data
  bin_te_e1$gene_count <- gene_counts
  
  # Step 4: (Optional) Rename to make plotting label clear
  bin_te_e1$gene_density <- gene_counts  # This is now genes per 100kb
  
  # Fit a linear model
  model <- lm(te_pct ~ E1 + GC + gene_density, data = bin_te_e1)
  summary(model)
  
  # Remove rows with NAs in required predictors before fitting the model
  clean_data <- bin_te_e1 %>%
    filter(!is.na(E1), !is.na(GC), !is.na(gene_density))
  
  # Fit the model on cleaned data
  model <- lm(te_pct ~ E1 + GC + gene_density, data = clean_data)
  
  # Add predictions to the same cleaned data
  clean_data$predicted <- predict(model)
  
  # Optional: plot predicted vs observed
  ggplot(clean_data, aes(x = predicted, y = te_pct)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_minimal(base_size = 14) +
    labs(title = "Predicted vs Observed TE Coverage",
         x = "Predicted TE %", y = "Observed TE %")
  
  bin_te_e1$E1_quartile <- ntile(bin_te_e1$E1, 4)
  bin_te_e1$E1_quartile <- factor(bin_te_e1$E1_quartile, labels = c("Q1 (lowest E1)", "Q2", "Q3", "Q4 (highest E1)"))
  
  ggplot(bin_te_e1, aes(x = E1_quartile, y = te_pct)) +
    geom_violin(fill = "lightgray", alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red") +
    theme_minimal(base_size = 14) +
    labs(
      title = "TE Coverage by E1 Quartile",
      x = "E1 Quartile (Hi-C Compartment Signal)",
      y = "TE Coverage (%)"
    )
 
  
  # === Step 1: Bin E1 into quartiles ===
  bin_te_e1 <- bin_te_e1 %>%
    mutate(
      E1_q = ntile(E1, 4),
      E1_q = factor(E1_q, labels = c("Q1 (lowest E1)", "Q2", "Q3", "Q4 (highest E1)"))
    )
  
  # === Step 2: Prepare long-format for multiple feature plots ===
  bin_te_e1_long <- bin_te_e1 %>%
    pivot_longer(cols = c(te_pct, GC, gene_density),
                 names_to = "feature",
                 values_to = "value")   
  
  bin_te_e1_long <- bin_te_e1_long %>%
    filter(!is.na(E1_q))
  # === Step 3: Boxplot per E1 quartile and feature ===
  ggplot(bin_te_e1_long, aes(x = E1_q, y = value, fill = E1_q)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    facet_wrap(~feature, scales = "free_y", ncol = 1,
               labeller = as_labeller(c(
                 te_pct = "TE Coverage (%)",
                 GC = "GC Content",
                 gene_density = "Gene Density (genes/100kb)"
               ))) +
    scale_fill_brewer(palette = "RdBu", direction = -1) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    labs(
      title = "Variation of Genomic Features Across Hi-C E1 Compartments",
      x = "E1 Quartile",
      y = "Feature Value"
    )
  
  
  summary_stats <- bin_te_e1_long %>%
    group_by(E1_q, feature) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_stats)
  
  
  kruskal_results <- bin_te_e1_long %>%
    group_by(feature) %>%
    summarise(
      p_value = kruskal.test(value ~ E1_q)$p.value,
      .groups = "drop"
    )
  
  kruskal_results
  
  
  # Prepare medians for annotation
  medians <- bin_te_e1_long %>%
    group_by(E1_q, feature) %>%
    summarise(med = median(value, na.rm = TRUE), .groups = "drop")

  # 1. Remove NA quartile bins
  bin_te_e1_long_clean <- bin_te_e1_long %>% filter(!is.na(E1_q))
  
  # 2. Compute p-values per feature
  pvals <- bin_te_e1_long_clean %>%
    group_by(feature) %>%
    summarise(p = formatC(kruskal.test(value ~ E1_q)$p.value, format = "e", digits = 2))
  
  # 3. Get maximum y for each facet to position labels above boxplots
  y_pos <- bin_te_e1_long_clean %>%
    group_by(feature) %>%
    summarise(y = max(value, na.rm = TRUE) * 1.1)
  
  # 4. Combine y positions and p-values into annotation data
  annotations <- left_join(pvals, y_pos, by = "feature")
  
  # 5. Plot with per-facet p-value labels using geom_text()
  ggplot(bin_te_e1_long_clean, aes(x = E1_q, y = value, fill = E1_q)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    facet_wrap(~feature, scales = "free_y", ncol = 1) +
    geom_text(data = annotations, aes(x = 1.5, y = y, label = paste0("p = ", p)), inherit.aes = FALSE, size = 5) +
    scale_fill_brewer(palette = "RdBu", direction = -1) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = "Genomic Feature Distributions Across E1 Quartiles",
      x = "E1 Quartile",
      y = "Feature Value"
    )
  
  
  
  
  # Filter out NA quartile bins
  bin_te_e1_long_clean <- bin_te_e1_long %>% 
    filter(E1_q %in% c("Q1 (lowest E1)", "Q4 (highest E1)"))  # Focus on Q1 vs Q4
  
  # Set factor levels to order Q1 < Q4
  bin_te_e1_long_clean$E1_q <- factor(bin_te_e1_long_clean$E1_q, 
                                      levels = c("Q1 (lowest E1)", "Q4 (highest E1)"))
  
  # Pairwise Wilcoxon test between Q1 and Q4, per feature
  pairwise_results <- bin_te_e1_long_clean %>%
    group_by(feature) %>%
    summarise(
      p = wilcox.test(value ~ E1_q)$p.value,
      label = paste0("p = ", formatC(p, format = "e", digits = 2)),
      .groups = "drop"
    )
  
  # Determine Y position for placing p-values
  y_pos <- bin_te_e1_long_clean %>%
    group_by(feature) %>%
    summarise(y = max(value, na.rm = TRUE) * 1.05)
  
  # Combine positions and labels
  annotations <- left_join(pairwise_results, y_pos, by = "feature")
  
  # Final plot
  ggplot(bin_te_e1_long_clean, aes(x = E1_q, y = value, fill = E1_q)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
    geom_text(data = annotations, aes(x = 1.5, y = y, label = label),
              inherit.aes = FALSE, size = 5) +
    facet_wrap(~feature, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Q1 (lowest E1)" = "#4575b4", 
                                 "Q4 (highest E1)" = "#d73027")) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    ) +
    labs(
      title = "Q1 vs Q4 Comparison Across Genomic Features",
      subtitle = "Pairwise Wilcoxon p-values for E1 extremes",
      x = "E1 Quartile",
      y = "Feature Value"
    )
  
  
  
  # Ensure E1_q is properly ordered factor
  bin_te_e1_long$E1_q <- factor(bin_te_e1_long$E1_q,
                                levels = c("Q1 (lowest E1)", "Q2", "Q3", "Q4 (highest E1)"))
  
  # Optional: filter out NA bins
  bin_te_e1_long_clean <- bin_te_e1_long %>% 
    filter(!is.na(E1_q))
  
  # Generate all pairwise comparisons
  quartile_levels <- levels(bin_te_e1_long_clean$E1_q)
  comparisons <- combn(quartile_levels, 2, simplify = FALSE)
  
  # Plot with all pairwise Wilcoxon test annotations
  ggplot(bin_te_e1_long_clean, aes(x = E1_q, y = value, fill = E1_q)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.signif", label.y.npc = "top", 
                       size = 4, tip.length = 0.01) +
    facet_wrap(~feature, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "RdBu", direction = -1) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    ) +
    labs(
      title = "Genomic Feature Distributions Across E1 Quartiles",
      subtitle = "All pairwise comparisons (Wilcoxon, p.signif)",
      x = "E1 Quartile",
      y = "Feature Value"
    )
  
  # Run and store results
  stat_table <- compare_means(value ~ E1_q, data = bin_te_e1_long_clean,
                              method = "wilcox.test", group.by = "feature", 
                              p.adjust.method = "BH", pairwise = TRUE)
  
  # View the result
  print(stat_table)  
  
#########################################
  # === Filter for Chr9A bins only
  chr9_bins <- bin_te_e1_long_clean %>%
    filter(chrom == "Chr9A") %>%
    mutate(
      bin_size_kb = (end - start + 1) / 1000,
      te_pct = 100 * te_bp / (end - start + 1),
      gene_density = gene_count / (bin_size_kb),  # genes per 100kb if bin_size_kb = 100
      E1_quartile = ntile(E1, 4)
    ) %>%
    mutate(
      E1_quartile = factor(E1_quartile,
                           levels = 1:4,
                           labels = c("Q1 (lowest E1)", "Q2", "Q3", "Q4 (highest E1)"))
    )
  
  # === Convert to long format for plotting
  chr9_long_clean <- chr9_bins %>%
    filter(!is.na(E1_quartile)) %>%
    select(E1_quartile, feature, value)
  # === Plot without NA
  ggplot(chr9_long_clean, aes(x = E1_quartile, y = value)) +
    geom_boxplot(aes(fill = E1_quartile), alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
    facet_wrap(~feature, scales = "free_y", ncol = 1, labeller = as_labeller(
      c(GC = "GC Content", gene_density = "Gene Density (genes/100kb)", te_pct = "TE Coverage (%)")
    )) +
    scale_fill_brewer(palette = "RdBu") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Chr9B: Genomic Feature Distributions Across E1 Quartiles",
      x = "E1 Quartile",
      y = "Feature Value"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "none"
    )
  
  
  comparisons <- combn(quartile_levels, 2, simplify = FALSE)
  
  label_y_df <- data.frame(
    group1 = rep(sapply(comparisons, `[[`, 1), 3),
    group2 = rep(sapply(comparisons, `[[`, 2), 3),
    feature = rep(c("GC", "gene_density", "te_pct"), each = length(comparisons)),
    y.position = c(
      rep(0.49, length(comparisons)),      # for GC
      rep(0.22, length(comparisons)),      # for gene_density (genes/100kb scale)
      rep(105, length(comparisons))        # for te_pct
    )
  )
  
  # Faceted boxplot for GC, gene_density, te_pct across E1 quartiles in Chr9A
  ggplot(chr9_long_clean, aes(x = E1_quartile, y = value, fill = E1_quartile)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
    facet_wrap(~feature, scales = "free_y", ncol = 1,
               labeller = as_labeller(c(
                 GC = "GC Content",
                 gene_density = "Gene Density (genes/100kb)",
                 te_pct = "TE Coverage (%)"
               ))) +
    stat_compare_means(
      comparisons = comparisons,
      method = "t.test",
      label = "p.signif",
      group.by = "feature",  # ✅ ensures pairwise t-tests per facet
      size = 4
    ) +
    scale_fill_manual(values = c(
      "Q1 (lowest E1)" = "#B2182B", 
      "Q2" = "#EF8A62", 
      "Q3" = "#67A9CF", 
      "Q4 (highest E1)" = "#2166AC"
    )) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Chr9B: Genomic Feature Distributions Across E1 Quartiles",
      subtitle = "Pairwise t-tests by feature",
      x = "E1 Quartile",
      y = "Feature Value"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  table(bin_te_e1$gene_count)
  summary(bin_te_e1$gene_count)
  hist(bin_te_e1$gene_count, breaks = 50)