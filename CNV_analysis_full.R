library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)

# Read data
depth <- read_tsv("../BLR_data/Avr_GWAS/hap1_depth_matrix.tsv", col_names = FALSE)
depth2 <- read_tsv("../BLR_data/Avr_GWAS/hap2_depth_matrix.tsv", col_names = FALSE)

depth

# Add column names: CHROM, POS, and your actual sample names
# Assign real column names: CHROM, POS, sample names
colnames(depth) <- c("CHROM", "POS",
                     "BLR478_hap1", "BLR479_hap1", "BLR482_hap1", "BLR483_hap1", "BLR484_hap1",
                     "BLR485_hap1", "BLR486_hap1", "BLR487_hap1", "BLR490_hap1", "BLR491_hap1",
                     "BLR506_hap1", "BLR507_hap1", "BLR520_hap1", "BLR531_hap1", "BLR537_hap1",
                     "BLR542_hap1", "BLR545_hap1", "BLR560_hap1", "BLR561_hap1", "BLR570_hap1",
                     "BLR577_hap1", "BLR584_hap1", "BLR608_hap1", "BLR612_hap1", "BLR623_hap1",
                     "BLR626_hap1", "BLR627_hap1", "BLR637_hap1", "BLR639_hap1", "BLR640_hap1",
                     "BLR660_hap1", "BLR671_hap1", "BLR672_hap1", "BLR675_hap1", "BLR685_hap1",
                     "BLR690_hap1", "BLR691_hap1", "P7VM253_hap1", "BLR488_hap1", "BLR489_hap1"
)

colnames(depth2) <- c("CHROM", "POS",
                     "BLR478_hap2", "BLR479_hap2", "BLR482_hap2", "BLR483_hap2", "BLR484_hap2",
                     "BLR485_hap2", "BLR486_hap2", "BLR487_hap2", "BLR490_hap2", "BLR491_hap2",
                     "BLR506_hap2", "BLR507_hap2", "BLR520_hap2", "BLR531_hap2", "BLR537_hap2",
                     "BLR542_hap2", "BLR545_hap2", "BLR560_hap2", "BLR561_hap2", "BLR570_hap2",
                     "BLR577_hap2", "BLR584_hap2", "BLR608_hap2", "BLR612_hap2", "BLR623_hap2",
                     "BLR626_hap2", "BLR627_hap2", "BLR637_hap2", "BLR639_hap2", "BLR640_hap2",
                     "BLR660_hap2", "BLR671_hap2", "BLR672_hap2", "BLR675_hap2", "BLR685_hap2",
                     "BLR690_hap2", "BLR691_hap2", "P7VM253_hap2", "UnnamedSample","BLR488_hap2", "BLR489_hap2"
)

depth <- depth[, !is.na(colnames(depth))]
#depth2 <- depth2[, !is.na(colnames(depth2))]


colnames(depth)

library(dplyr)
depth <- depth %>% select(-UnnamedSample)
depth %>% dplyr::select(-UnnamedSample)

depth_long <- depth %>%
  pivot_longer(cols = -c(CHROM, POS), names_to = "Sample", values_to = "Depth")

depth_long

depth_long$Depth <- as.numeric(depth_long$Depth)

depth_long <- depth_long %>%
  group_by(Sample) %>%
  mutate(NormDepth = Depth / median(Depth, na.rm = TRUE)) %>%
  ungroup()

ggplot(depth_long %>% filter(CHROM == "chr_9A"), aes(x = POS, y = NormDepth, color = Sample)) +
  geom_line(alpha = 0.6) +
  labs(title = "Copy Number Proxy for chr_1B", x = "Genomic Position", y = "Normalized Depth") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if too crowded


depth_long_binned <- depth_long %>%
  mutate(bin = floor(POS / 1000) * 1000) %>%
  group_by(CHROM, bin, Sample) %>%
  summarise(NormDepth = median(NormDepth, na.rm = TRUE), .groups = "drop")

library(ggplot2)
library(forcats)

# Reorder chromosomes (optional: adjust if your chromosomes have a specific order)
depth_long_binned$CHROM <- factor(
  depth_long_binned$CHROM,
  levels = paste0("chr_", c(1:18, "1A", "2A", "3A", "4A", "5A", "6A", "7A", "8A", "9A", "10A", "11A", "12A", "13A", "14A", "15A", "16A", "17A", "18A"))
)

depth_long_binned$CHROM <- factor(
  depth_long_binned$CHROM,
  levels = paste0("chr_", c(1:18, "1B", "2B", "3B", "4B", "5B", "6B", "7B", "8B", "9B", "10B", "11B", "12B", "13B", "14B", "15B", "16B", "17B", "18B"))
)


# Plot with clearer labels
library(scales)  # for label_number

ggplot(depth_long_binned, aes(x = bin / 1e6, y = NormDepth, color = Sample)) +  # convert bin to Mb
  geom_line(alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~ CHROM, scales = "free_x", ncol = 6) +
  labs(
    title = "Normalized Read Depth by Chromosome",
    x = "Genomic Position (Mb)",
    y = "Normalized Depth",
    color = "Sample"
  ) +
  scale_x_continuous(labels = label_number(suffix = " Mb", scale = 1)) +  # Mb labels
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 8)
  )

# 1. Get unique sample names
samples <- unique(depth_long_binned$Sample)

# 2. Assign visually distinct colors to each sample
sample_colors <- setNames(
  rainbow(length(samples)),  # you can substitute with other palettes too
  samples
)

# 3. Plot using manual color scale
ggplot(depth_long_binned, aes(x = bin / 1e6, y = NormDepth, color = Sample)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ CHROM, scales = "free_x") +
  scale_color_manual(values = sample_colors) +
  labs(
    title = "Normalized Read Depth by Chromosome",
    x = "Genomic Position (Mb)",
    y = "Normalized Depth"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#####Circos plot????

#Heatmap now
library(dplyr)

chrom_depth <- depth_long %>%
  group_by(Sample, CHROM) %>%
  summarise(MeanNormDepth = mean(NormDepth, na.rm = TRUE)) %>%
  ungroup()

depth_long <- depth_long %>%
  group_by(Sample) %>%
  mutate(NormDepth = scale(Depth)) %>%  # subtract mean, divide by SD
  ungroup()
#####normalised method change
depth_long <- depth_long %>%
  group_by(Sample) %>%
  mutate(NormDepth = Depth / mean(Depth, na.rm = TRUE)) %>%
  ungroup()

###depth_long <- depth_long %>%
  #group_by(Sample) %>%
  #mutate(NormDepth = (Depth - median(Depth, na.rm = TRUE)) /
           #mad(Depth, constant = 1, na.rm = TRUE)) %>%
  #ungroup()

depth_wide <- chrom_depth %>%
  pivot_wider(names_from = CHROM, values_from = MeanNormDepth)

# Turn Sample into rownames
depth_matrix <- as.data.frame(depth_wide)
rownames(depth_matrix) <- depth_matrix$Sample
depth_matrix$Sample <- NULL
depth_matrix

library(pheatmap)

pheatmap(as.matrix(depth_matrix),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Normalized Copy Number Heatmap (per Chromosome)",
         fontsize_row = 8,
         fontsize_col = 10)

dev.off()

library(DNAcopy)

#combined_segments$estimated_CN <- round(2 * 2^combined_segments$seg.mean, 2)

# Collapse segments into average per chromosome


################################

bin_size <- 10000  # 100kb

depth_binned <- depth_long %>%
  mutate(Bin = floor(POS / bin_size) * bin_size) %>%
  group_by(Sample, CHROM, Bin) %>%
  summarise(MeanDepth = mean(Depth, na.rm = TRUE), .groups = "drop")

depth_binned <- depth_binned %>%
  group_by(Sample) %>%
  mutate(NormDepth = MeanDepth / median(MeanDepth, na.rm = TRUE),
         log2depth = log2(pmax(NormDepth, 1e-3))) %>%
  ungroup()

segment_sample <- function(df, sample_name) {
  df_sample <- df %>%
    filter(Sample == sample_name) %>%
    arrange(CHROM, Bin) %>%
    distinct(CHROM, Bin, .keep_all = TRUE)
  
  cna_obj <- CNA(genomdat = df_sample$log2depth,
                 chrom = as.character(df_sample$CHROM),
                 maploc = df_sample$Bin,
                 data.type = "logratio",
                 sampleid = sample_name)
  
  smoothed <- smooth.CNA(cna_obj)
  segmented <- segment(smoothed, verbose = 0)
  segmented$output
}

samples <- unique(depth_binned$Sample)
samples
all_segments <- list()

summary(depth_binned$NormDepth)


for (sample in samples) {
  seg_df <- segment_sample(depth_binned, sample)
  all_segments[[sample]] <- seg_df
  
  write.table(seg_df, file = paste0("../BLR_data/Avr_GWAS/cnv_segments", sample, "_segments.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

combined_segments <- bind_rows(all_segments, .id = "Sample")

ggplot(combined_segments, aes(x = seg.mean)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Distribution of Segment log2 Ratios",
       x = "Segment Mean (log2 CN ratio)", y = "Count")

combined_segments <- combined_segments %>%
  mutate(CNV_type = case_when(
    seg.mean > 0.6 ~ "Duplication",
    seg.mean < -0.6 ~ "Deletion",
    TRUE ~ "Normal"
  ))

combined_segments

# Per sample summary
table(combined_segments$CNV_type)

# Per-sample barplot
ggplot(combined_segments, aes(x = Sample, fill = CNV_type)) +
  geom_bar() +
  theme_minimal() +
  scale_fill_manual(values = c("Deletion" = "blue", "Duplication" = "red", "Normal" = "gray")) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Number of CNV Segments per Sample", y = "Count")

############

bed_export <- combined_segments %>%
  mutate(start = loc.start,
         end = loc.end,
         chrom = paste0("chr", chrom),
         name = paste0(Sample, "_", round(seg.mean, 3))) %>%
  select(chrom, start, end, name, seg.mean)

write.table(bed_export, "../BLR_data/Avr_GWAS/all_hap1_cnv_segments.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Step 1: Get average log2 CN ratio per chromosome per sample
cnv_chr_matrix <- combined_segments %>%
  group_by(Sample, chrom) %>%
  summarise(avg_log2CN = mean(seg.mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = chrom, values_from = avg_log2CN)

# Convert to matrix for heatmap
cnv_mat <- as.data.frame(cnv_chr_matrix)
rownames(cnv_mat) <- cnv_mat$Sample
cnv_mat <- cnv_mat[, -1]  # drop 'Sample' column

# Step 3: Plot heatmap
pheatmap::pheatmap(cnv_mat,
                   scale = "column",
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   color = colorRampPalette(c("blue", "white", "red"))(100),
                   main = "Average CNV per Chromosome (10kb bins)",
                   fontsize_row = 6,
                   fontsize_col = 10) 

dev.off()

library(readr)
library(dplyr)

# Read the .fam or population assignment file
group_df <- read_table("../BLR_data/BLR560_assembly/evalmix/pruned_hap2_out.fam.txt", col_names = FALSE)
group_df <- read.table("../BLR_data/BLR560_assembly/Hap1_based_pop_info.fam", header=FALSE)

# Check structure — often .fam files look like:
# FID IID PID MID Sex Phenotype (we just care about FID or IID)
head(group_df) 

# Ensure column names
colnames(group_df)[1:2] <- c("Sample", "Group")  # Or adjust if needed

# Ensure it's a plain data.frame
group_df <- as.data.frame(group_df)

group_df

# Clean whitespace (just in case)
group_df$Sample <- trimws(group_df$Sample)
rownames(cnv_mat) <- trimws(rownames(cnv_mat))

# Build annotation_row (matching samples only)
annotation_row <- group_df %>%
  filter(Sample %in% rownames(cnv_mat)) %>%
  distinct(Sample, Group)

# Reformat: use sample names as rownames
annotation_row <- as.data.frame(annotation_row)
rownames(annotation_row) <- annotation_row$Sample
annotation_row <- annotation_row[, "Group", drop = FALSE]

# Final check
str(annotation_row)
head(annotation_row) 

# Assign colors to populations
group_levels <- unique(annotation_row$Group)

custom_colors <- c(
  "pop1" = "#57B893",  # blue-green
  "pop2" = "#F87850",  # orange-red
  "pop3" = "#7B8DBF",  # soft blue
  "pop4" = "#D771B6",  # pink-purple
  "pop5" = "#B2DF8A"   # light green
)

group_colors <- list(Group = custom_colors)


#group_colors <- list(Group = setNames(
  #brewer.pal(max(3, length(group_levels)), "Set1")[1:length(group_levels)],
  #group_levels
#))

# Plot!
pheatmap(
  cnv_mat,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward",
  annotation_row = annotation_row,
  annotation_colors = group_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Average CNV per Chromosome (10kb bins)",
  fontsize_row = 6,
  fontsize_col = 10
)

library(GenomicRanges)
library(rtracklayer)

cnv_df <- combined_segments

# Example: load CNVs as GRanges and overlap with genes
cnv_gr <- GRanges(seqnames = cnv_df$chrom,
                  ranges = IRanges(start = cnv_df$loc.start, end = cnv_df$loc.end),
                  seg.mean = cnv_df$seg.mean)

gene_gr <- import("../BLR_data/BLR560_assembly/BLR560_hap1_gffread.gff3")

# Find overlaps
overlaps <- findOverlaps(cnv_gr, gene_gr)
cnv_genes <- gene_gr[subjectHits(overlaps)]

cnv_genes

# Keep only transcript or gene features
gene_features <- cnv_genes[mcols(cnv_genes)$type %in% c("transcript", "gene")]



# Extract unique gene IDs
gene_ids <- unique(mcols(gene_features)$geneID)

# Preview or save
head(gene_ids)

####map gene to cnv type 
# Step 1: Keep only gene or transcript features
gene_features <- gene_gr[mcols(gene_gr)$type %in% c("gene", "transcript")]

# Step 2: Find overlaps with CNVs
overlaps <- findOverlaps(cnv_gr, gene_features)

# Step 3: Build mapping: which CNV overlaps which gene
cnv_hits <- cnv_df[queryHits(overlaps), ]
gene_hits <- gene_features[subjectHits(overlaps)]

# Step 4: Combine info
gene_cnv_df <- data.frame(
  Sample = cnv_hits$Sample,
  CNV_type = cnv_hits$CNV_type,
  seg.mean = cnv_hits$seg.mean,
  chr = as.character(seqnames(gene_hits)),
  start = start(gene_hits),
  end = end(gene_hits),
  geneID = mcols(gene_hits)$geneID,
  transcriptID = mcols(gene_hits)$ID
)

# Optional: drop duplicates if needed
gene_cnv_df <- gene_cnv_df[!duplicated(gene_cnv_df[c("Sample", "geneID", "CNV_type")]), ]
gene_cnv_df
gene_cnv_filtered <- gene_cnv_df %>%
  filter(CNV_type != "Normal")

gene_cnv_filtered_dup <- gene_cnv_filtered %>%
  filter(CNV_type != "Deletion")

gene_cnv_filtered_del <- gene_cnv_filtered %>%
  filter(CNV_type != "Duplication")

head(gene_cnv_filtered)

write.table(
  gene_cnv_filtered_dup,
  file = "../BLR_data/BLR560_assembly/hap1_gene_cnv_dup_only.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_cnv_filtered_del,
  file = "../BLR_data/BLR560_assembly/hap1_gene_cnv_del_only.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

library(dplyr)

# Step 1: Remove transcript redundancy (distinct Sample–gene pairs)
sample_gene_list <- gene_cnv_filtered_dup %>%
  distinct(Sample, geneID)

# Step 2: Count how many samples each gene appears in
gene_sample_counts <- sample_gene_list %>%
  count(geneID, name = "sample_count") %>%
  arrange(desc(sample_count))

# Step 3: Filter by a threshold (e.g. genes in ≥3 samples)
shared_genes <- gene_sample_counts %>%
  filter(sample_count >= 3)

# View top shared genes
shared_genes

library(ggtext)

# Join group info
gene_cnv_annotated <- gene_cnv_filtered %>%
  left_join(annotation_row %>% rownames_to_column("Sample"), by = "Sample")

# Create colored x-axis labels using your group_colors$Group
gene_cnv_annotated <- gene_cnv_annotated %>%
  mutate(Sample_colored = paste0(
    "<span style='color:", group_colors$Group[Group], "'>", Sample, "</span>"
  ))

ggplot(gene_cnv_annotated, aes(x = Sample_colored, fill = CNV_type)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("Deletion" = "#7B8DBF", "Duplication" = "#57B893")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "CNV-Affected Genes in hap1 per Sample",
    y = "Number of Genes", x = "Sample"
  ) +
  theme(
    axis.text.x = element_markdown(angle = 90, hjust = 1, size = 8),  # <- Group-colored x-axis labels
    legend.position = "bottom"
  )


ggplot(gene_cnv_filtered, aes(x = Sample, fill = CNV_type)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("Deletion" = "#7B8DBF", "Duplication" = "#57B893")) +
  theme_minimal() +
  labs(title = "CNV-Affected Genes in hap1 per Sample",
       y = "Number of Genes", x = "Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

########Heatmap: Gene × Sample matrix (binary: duplicated, deleted)
# Step 1: Binary matrix of gene × sample
gene_matrix <- gene_cnv_filtered %>%
  distinct(geneID, Sample) %>%     # remove duplicates
  mutate(flag = 1) %>%
  pivot_wider(names_from = Sample, values_from = flag, values_fill = 0) %>%
  column_to_rownames("geneID")


write.table(
  gene_matrix,
  file = "../BLR_data/BLR560_assembly/hap1_cyp51_cnv.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Step 2: Heatmap
pheatmap::pheatmap(as.matrix(gene_matrix),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   clustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   color = colorRampPalette(c("white", "red"))(100),
                   fontsize_row = 6,
                   fontsize_col = 6,
                   main = "Presence of CNV-Affected Genes per Sample")


# Step 1: Subset your CNV table
highlight_genes <- c("g8344","g8347")

# Subset from the full CNV-gene annotation table
all_gene_states <- gene_cnv_df %>%
  filter(geneID %in% highlight_genes) %>%
  select(Sample, geneID, seg.mean) %>%
  pivot_wider(names_from = Sample, values_from = seg.mean)

# Format for heatmap
all_gene_states <- as.data.frame(all_gene_states)
rownames(all_gene_states) <- all_gene_states$geneID
all_gene_states$geneID <- NULL

all_gene_states

#gene_cnv_filtered %>% filter(transcriptID == "g8344") %>% select(Sample, seg.mean)

# Plot
pheatmap::pheatmap(as.matrix(all_gene_states),
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   clustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   color = colorRampPalette(c("#7B8DBF", "white", "red"))(99),
                   breaks = seq(-1.5, 1.5, length.out = 100),
                   fontsize_row = 10,
                   fontsize_col = 8,
                   main = "CNV log2 ratios for g8347 and g8344 across all samples")


pheatmap::pheatmap(as.matrix(all_gene_states),
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   clustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   color = colorRampPalette(c("#7B8DBF", "white", "#57B893"))(99),
                   breaks = seq(-0.1, 0.1, length.out = 100),
                   fontsize_row = 10,
                   fontsize_col = 8,
                   main = "CNV log2 ratios for g8089 across all samples")

dev.off()

pheatmap::pheatmap(
  as.matrix(all_gene_states),
  cluster_cols = TRUE,
  cluster_rows = FALSE,
  clustering_distance_cols = "euclidean",
  clustering_method = "median",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-0.1, 0.45, length.out = 100),  # tighter range
  main = "CNV log2 ratios for g8089 (rescaled)"
)

zscore_matrix <- t(scale(t(as.matrix(all_gene_states))))

pheatmap::pheatmap(
  zscore_matrix,
  cluster_cols = TRUE,
  cluster_rows = FALSE,
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "CNV log2 ratios for g8089 (z-score normalized)"
)



write.table(
  all_gene_states,
  file = "../BLR_data/BLR560_assembly/hap1_cyp51_cnv.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 1. Filter to deletions only
deleted_genes <- gene_cnv_filtered %>%
  filter(CNV_type == "Deletion") %>%
  distinct(geneID, Sample) %>%   # remove any duplicate hits
  mutate(present = 1) %>%
  pivot_wider(names_from = Sample, values_from = present, values_fill = 0)

# 2. Install and load ComplexUpset if needed
# install.packages("ComplexUpset")
library(ComplexUpset)

# 3. Convert to data.frame
deleted_genes_df <- as.data.frame(deleted_genes)

# 4. Plot UpSet (excluding genes present in <2 samples, optional)
upset(deleted_genes_df, 
      intersect = colnames(deleted_genes_df)[-1],  # all sample columns
      min_size = 2,                                # only shared genes
      name = "Deleted Genes",
      width_ratio = 0.2,
      sort_sets = FALSE)


# Step 1: Filter only meaningful CNVs (optional)
cnv_pca_input <- combined_segments %>%
  filter(CNV_type %in% c("Deletion", "Duplication"))  # skip Normal if needed

# Step 2: Create unique region IDs
cnv_pca_input <- cnv_pca_input %>%
  mutate(region_id = paste(chrom, loc.start, loc.end, sep = "_"))

# Step 3: Pivot to wide format (samples as rows, regions as columns)
cnv_matrix <- cnv_pca_input %>%
  select(Sample, region_id, seg.mean) %>%
  pivot_wider(names_from = region_id, values_from = seg.mean, values_fill = 0) %>%
  column_to_rownames("Sample")


cnv_matrix_scaled <- scale(cnv_matrix, center = TRUE, scale = TRUE)
pca_result <- prcomp(cnv_matrix_scaled)

cnv_matrix_no_outliers <- cnv_matrix[!rownames(cnv_matrix) %in% c("BLR542_hap1", "BLR545_hap1"), ]

# Remove constant columns (zero variance)
cnv_matrix_filtered <- cnv_matrix_no_outliers[, apply(cnv_matrix_no_outliers, 2, function(x) sd(x) > 0)]
pca_result <- prcomp(cnv_matrix_filtered , scale. = TRUE)

# Step 6: Prepare PCA results
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(pca_df)

pca_df$Group <- annotation_row[rownames(pca_df), "Group"]


# Calculate group centroids
group_centroids <- pca_df %>%
  group_by(Group) %>%
  summarize(centroid_PC1 = mean(PC1), centroid_PC2 = mean(PC2), .groups = "drop")

# Join centroids back to pca_df
pca_with_centroids <- left_join(pca_df, group_centroids, by = "Group")

ggplot(pca_with_centroids, aes(x = PC1, y = PC2, color = Group)) +
  # Cross axes with ticks
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  
  # Dashed lines from samples to centroids
  geom_segment(aes(xend = centroid_PC1, yend = centroid_PC2), 
               color = "gray60", linewidth = 0.5, linetype = "dashed") +
  
  # Points for samples
  geom_point(size = 3, alpha = 0.9) +
  
  # Sample labels
  ggrepel::geom_text_repel(aes(label = Sample), 
                           size = 3, box.padding = 0.5, 
                           segment.color = "gray70", max.overlaps = 100) +
  
  # Manual group colors
  scale_color_manual(values = custom_colors) +
  
  # Title and axis labels with variance %
  labs(
    title = "PCA of CNV Profiles (Deletions + Duplications)",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), " %)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), " %)")
  ) +
  
  # Theme tweaks for axis visibility
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none"
  )

ggplot(pca_with_centroids, aes(x = PC1, y = PC2, color = Group)) +
  # Crosshairs at origin
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  
  # Dotted lines from points to centroids
  geom_segment(aes(xend = centroid_PC1, yend = centroid_PC2), 
               color = "gray60", linewidth = 0.4, linetype = "dashed") +
  
  # Sample points
  geom_point(size = 3, alpha = 0.9) +
  
  # Labels for each sample
  geom_text_repel(aes(label = Sample), size = 3, box.padding = 0.3, 
                  segment.color = "gray60", max.overlaps = 100) +
  
  # Colors for groups
  scale_color_manual(values = custom_colors) +
  
  # Axis labels with variance explained
  labs(
    title = "PCA of CNV Profiles (Deletions + Duplications)",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), " %)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), " %)")
  ) +
  
  # Make background white with black axis lines and border
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # Border box!
    axis.line = element_blank(),  # Already using hline/vline
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none"
  )

# Filter Deletions and Duplications separately
cnv_deletions <- combined_segments %>%
  filter(CNV_type == "Deletion")

cnv_duplications <- combined_segments %>%
  filter(CNV_type == "Duplication")

# Convert to wide format (genes/regions as features)
deletion_matrix <- cnv_deletions %>%
  select(Sample, chrom, loc.start, loc.end, seg.mean) %>%
  mutate(region = paste0(chrom, ":", loc.start, "-", loc.end)) %>%
  select(Sample, region, seg.mean) %>%
  pivot_wider(names_from = region, values_from = seg.mean, values_fill = 0)

duplication_matrix <- cnv_duplications %>%
  select(Sample, chrom, loc.start, loc.end, seg.mean) %>%
  mutate(region = paste0(chrom, ":", loc.start, "-", loc.end)) %>%
  select(Sample, region, seg.mean) %>%
  pivot_wider(names_from = region, values_from = seg.mean, values_fill = 0)

# Set rownames for PCA
deletion_mat <- as.data.frame(deletion_matrix)
duplication_mat <- as.data.frame(duplication_matrix)

rownames(deletion_mat) <- deletion_mat$Sample
rownames(duplication_mat) <- duplication_mat$Sample

deletion_mat <- deletion_mat[ , -1]
duplication_mat <- duplication_mat[ , -1]

# Remove zero-variance columns
deletion_mat <- deletion_mat[, apply(deletion_mat, 2, sd) > 0]
duplication_mat <- duplication_mat[, apply(duplication_mat, 2, sd) > 0]

pca_del <- prcomp(deletion_mat, scale. = TRUE)
pca_dup <- prcomp(duplication_mat, scale. = TRUE)

pca_df_del <- as.data.frame(pca_del$x)
pca_df_del$Sample <- rownames(pca_df_del)
pca_df_del <- left_join(pca_df_del, group_df, by = "Sample")

pca_df_dup <- as.data.frame(pca_dup$x)
pca_df_dup$Sample <- rownames(pca_df_dup)
pca_df_dup <- left_join(pca_df_dup, group_df, by = "Sample")

library(ggrepel)

ggplot(pca_df_del, aes(x = PC1, y = PC2, color = Group)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = Sample), size = 3, box.padding = 0.3) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "PCA of CNV Deletions",
    x = paste0("PC1 (", round(summary(pca_del)$importance[2,1]*100, 1), " %)"),
    y = paste0("PC2 (", round(summary(pca_del)$importance[2,2]*100, 1), " %)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )

# Step 1: Filter to deletions and duplications
cnv_filtered <- combined_segments %>%
  filter(CNV_type %in% c("Deletion", "Duplication"))

# Step 2: Summarize CNV segments per sample into matrix
cnv_matrix <- cnv_filtered %>%
  select(Sample, chrom, loc.start, loc.end, seg.mean, CNV_type) %>%
  mutate(region = paste(chrom, loc.start, loc.end, sep = "_")) %>%
  select(Sample, region, seg.mean, CNV_type) %>%
  pivot_wider(names_from = region, values_from = seg.mean, values_fill = 0)

# Step 3: Separate into deletion and duplication matrices
deletion_mat <- cnv_matrix %>% filter(CNV_type == "Deletion") %>% select(-CNV_type)
duplication_mat <- cnv_matrix %>% filter(CNV_type == "Duplication") %>% select(-CNV_type)

# Step 4: Convert to matrix and remove zero-variance columns
mat_filter <- function(df) {
  mat <- as.data.frame(df)
  rownames(mat) <- mat$Sample
  mat <- mat[, -1]
  mat <- mat[, apply(mat, 2, sd) > 0]
  return(mat)
}

del_mat <- mat_filter(deletion_mat)
dup_mat <- mat_filter(duplication_mat)

# Step 5: Remove outliers (optional)
remove_outliers_pca <- function(mat) {
  pca <- prcomp(mat, scale. = TRUE)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$Sample <- rownames(mat)
  scores <- scores %>%
    filter(abs(PC1) < 3 * sd(PC1), abs(PC2) < 3 * sd(PC2))
  mat[scores$Sample, ]
}

del_mat <- remove_outliers_pca(del_mat)
dup_mat <- remove_outliers_pca(dup_mat)

# Remove zero-variance columns
del_mat <- del_mat[, apply(del_mat, 2, sd) > 0]
dup_mat <- dup_mat[, apply(dup_mat, 2, sd) > 0]

# Step 6: Run PCA again
pca_del <- prcomp(del_mat, scale. = TRUE)
pca_dup <- prcomp(dup_mat, scale. = TRUE)

# Step 7: Prepare plotting data
pca_df <- function(pca, type) {
  df <- as.data.frame(pca$x[, 1:2])
  df$Sample <- rownames(df)
  df$Type <- type
  df
}

del_df <- pca_df(pca_del, "Deletion")
dup_df <- pca_df(pca_dup, "Duplication")

pca_all <- bind_rows(del_df, dup_df)

# Step 8: Merge with group info
group_df <- read.table("../BLR_data/BLR560_assembly/Hap1_based_pop_info.fam", header = FALSE)
colnames(group_df)[1:2] <- c("Sample", "Group")

pca_all <- left_join(pca_all, group_df, by = "Sample")

# Step 9: Plot
plot_pca <- function(df, title, color_map) {
  ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    scale_color_manual(values = color_map) +
    theme_minimal(base_size = 13) +
    labs(title = title,
         x = paste0("PC1 (", round(summary(prcomp(df[, c("PC1", "PC2")]))$importance[2, 1] * 100, 1), " %)"),
         y = paste0("PC2 (", round(summary(prcomp(df[, c("PC1", "PC2")]))$importance[2, 2] * 100, 1), " %)")) +
    theme(legend.position = "bottom")
}

# Define color palette
custom_colors <- c(
  "pop1" = "#57B893",
  "pop2" = "#F87850",
  "pop3" = "#7B8DBF",
  "pop4" = "#D771B6",
  "pop5" = "#B2DF8A"
)

# Create plots
plot_del <- plot_pca(filter(pca_all, Type == "Deletion"), "PCA: CNV Deletions", custom_colors)
plot_dup <- plot_pca(filter(pca_all, Type == "Duplication"), "PCA: CNV Duplications", custom_colors)




library(ggplot2)
library(dplyr)
library(ggrepel)

# --- Step 1: Helper to calculate group centroids ---
add_centroids <- function(df) {
  df %>%
    group_by(Group) %>%
    summarize(centroid_PC1 = mean(PC1), centroid_PC2 = mean(PC2), .groups = "drop") %>%
    left_join(df, by = "Group")
}

# --- Step 2: Updated plotting function ---
plot_pca_with_centroids <- function(df, title, color_map) {
  df <- add_centroids(df)
  
  ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
    # Crosshairs
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    
    # Dotted line to centroid
    geom_segment(aes(xend = centroid_PC1, yend = centroid_PC2), 
                 linetype = "dashed", linewidth = 0.3, color = "gray60") +
    
    # Points
    geom_point(size = 3, alpha = 0.9) +
    
    # Sample labels
    geom_text_repel(aes(label = Sample), size = 2, box.padding = 0.4,
                    segment.color = "gray60", max.overlaps = 100) +
    
    # Colors and axes
    scale_color_manual(values = color_map) +
    labs(
      title = title,
      x = paste0("PC1 (", round(summary(prcomp(df[, c("PC1", "PC2")]))$importance[2, 1] * 100, 1), " %)"),
      y = paste0("PC2 (", round(summary(prcomp(df[, c("PC1", "PC2")]))$importance[2, 2] * 100, 1), " %)")
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14, face = "bold"),
      legend.position = "none"
    )
}

plot_del <- plot_pca_with_centroids(filter(pca_all, Type == "Deletion"), 
                                    "PCA of CNV Deletion", 
                                    custom_colors)

plot_dup <- plot_pca_with_centroids(filter(pca_all, Type == "Duplication"), 
                                    "PCA of CNV Duplication", 
                                    custom_colors)

library(patchwork)
plot_del + plot_dup + plot_layout(ncol = 2)
