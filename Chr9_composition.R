library(ggplot2)
library(dplyr)
library(readr)

extract_genes_from_gff <- function(file, sample_label) {
  gff <- read_tsv(file, comment = "#", col_names = FALSE)
  colnames(gff)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  gff %>%
    filter(type == "gene") %>%
    mutate(sample = sample_label) %>%
    select(sample, seqid, start, end)
}

gff_560A <- extract_genes_from_gff("../Ph560A_chr9A.gff3", "Ph560A")
gff_560B <- extract_genes_from_gff("../Ph560B_chr9B.gff3", "Ph560B")
gff_518A <- extract_genes_from_gff("../Ph518A_chr9A.gff3", "Ph518A")
gff_518B <- extract_genes_from_gff("../Ph518B_chr9B.gff3", "Ph518B")

gene_data <- bind_rows(gff_560A, gff_560B, gff_518A, gff_518B)

ggplot(gene_data) +
  geom_segment(aes(x = start, xend = end, y = 0, yend = 0), color = "steelblue", linewidth = 0.7) +
  facet_wrap(~ sample, ncol = 1, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(x = "Chromosome Position (bp)", y = "Genes")


library(tidyverse)

# Read all GFFs and label them
extract_genes <- function(file, sample_label) {
  gff <- read_tsv(file, comment = "#", col_names = FALSE)
  colnames(gff)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  gff %>%
    filter(type == "gene") %>%
    mutate(sample = sample_label) %>%
    select(sample, seqid, start, end)
}

df_560A <- extract_genes("../Ph560A_chr9A.gff3", "Ph560A")
df_560B <- extract_genes("../Ph560B_chr9B.gff3", "Ph560B")
df_518A <- extract_genes("../Ph518A_chr9A.gff3", "Ph518A")
df_518B <- extract_genes("../Ph518B_chr9B.gff3", "Ph518B")

gene_data <- bind_rows(df_560A, df_560B, df_518A, df_518B)


gene_bed <- gene_data %>%
  mutate(chrom = seqid,
         start = start,
         end = end,
         name = paste0(sample, "_", row_number())) %>%
  select(chrom, start, end, name)

gene_bed <- gene_data %>%
  mutate(chr = paste0(sample, "_", seqid))

# Ensure correct column order and types
gene_bed_fixed <- gene_data %>%
  mutate(chr = paste0(sample, "_", seqid),
         start = as.integer(start),
         end = as.integer(end)) %>%
  select(chr, start, end)

# Save as valid tab-separated BED file
write.table(gene_bed_fixed, "../genes_for_karyoplot.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

library(karyoploteR)

# Start clean from the original full gene_data
gene_bed <- gene_data %>%
  mutate(chr = paste0(sample, "_", seqid),
         start = as.numeric(start),
         end = as.numeric(end)) %>%
  select(chr, start, end, sample)

# Create custom genome GRanges safely
custom.genome <- gene_bed %>%
  group_by(chr) %>%
  summarise(start = 1, end = max(end), .groups = "drop") %>%
  as.data.frame() %>%
  toGRanges()

# Gene GRanges
gene_gr <- toGRanges("../genes_for_karyoplot.bed")

# 1. Setup karyotype layout using your custom genome
kp <- plotKaryotype(genome = custom.genome, plot.type = 1)  # one per line

# 2. Draw chromosome ideograms (gray bar w/ borders)
kpAddBaseNumbers(kp, tick.dist = 2e6, add.units = TRUE)

# 3. Overlay gene features ON the bar (r0 = 0, r1 = 1 means full bar)
kpPlotRegions(kp, data = gene_gr,
              col = "#1f77b4", border = NA,
              r0 = 0, r1 = 1)  # full ideogram height

# Density (top half)
kpPlotDensity(kp, data = gene_gr, col = "#1f77b4", window.size = 50000, r0 = 0.5, r1 = 1)

# Gene ticks (bottom half)
kpPlotRegions(kp, data = gene_gr, col = "#1f77b4", r0 = 0, r1 = 0.4)




###############################
library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)

# Load chromosome size
fai <- read.table("../560A_chr9.fa.fai", header = FALSE)
chr9_len <- fai[fai$V1 == "Chr9A", 2]
custom_genome <- toGRanges(data.frame(chr="Chr9A", start=1, end=chr9_len))

# Load data
genes <- read.table("../Ph560A_chr9A_genes.bed", header = FALSE)
colnames(genes) <- c("chr", "start", "end", "gene_id")
gr_genes <- makeGRangesFromDataFrame(genes)

effectors <- read.table("../effector_genes_chr9A_by_suffix.bed", header = FALSE)
colnames(effectors) <- c("chr", "start", "end", "gene_id")
gr_effectors <- makeGRangesFromDataFrame(effectors)

breakpoints <- read.table("../560A_approx_breakpoints.bed", header = FALSE)
colnames(breakpoints) <- c("chr", "start", "end")
gr_bps <- makeGRangesFromDataFrame(breakpoints)

# Plot setup
kp <- plotKaryotype(genome = custom_genome, chromosomes = "Chr9A", plot.type = 1)

# Draw a horizontal chromosome bar
kpRect(kp, data=custom_genome, y0=0.45, y1=0.55, col="#D3D3D3", border="black")

# Draw genes as black ticks inside the bar
kpRect(kp, data=gr_genes, y0=0.45, y1=0.55, col="black", border=NA)

# Draw effectors in green over the same y-position (for overlay)
kpRect(kp, data=gr_effectors, y0=0.45, y1=0.55, col="forestgreen", border=NA)

# Draw breakpoints as red dashed vertical lines through the bar
kpSegments(kp, data=gr_bps, y0=0.40, y1=0.60, col="red", lwd=1.5, lty=2)

# Optional title
title("Chr9A: Genes (black), Effectors (green), Breakpoints (red)")

##############################


# Read in gene/effector/breakpoint data
genes <- read.table("../Ph560A_chr9A_genes.bed", header = FALSE, col.names = c("chr", "start", "end", "gene"))
effectors <- read.table("../effector_genes_chr9A_by_suffix.bed", header = FALSE, col.names = c("chr", "start", "end", "gene"))
bps <- read.table("../560A_approx_breakpoints.bed", header = FALSE, col.names = c("chr", "start", "end"))
# Load chromosome index
fai <- read.table("../560A_chr9.fa.fai", header = FALSE, stringsAsFactors = FALSE)
colnames(fai)[1:2] <- c("chr", "length")

# Create genome info as GRanges object
custom_genome <- GRanges(seqnames = fai$chr, ranges = IRanges(start = 1, end = fai$length))
autoplot(custom_genome)

# Label types
genes$type <- "Gene"
effectors$type <- "Effector"
features <- bind_rows(genes, effectors)



library(ggbio)
library(GenomicRanges)

dev.off()
# Convert to GRanges
gr_genes <- makeGRangesFromDataFrame(genes)
gr_effectors <- makeGRangesFromDataFrame(effectors)
gr_bps <- GRanges(seqnames = bps$chr, ranges = IRanges(bps$start, bps$end))

# Plot with tracks
tracks(
  Genes = autoplot(gr_genes, aes(color = "steelblue")),
  Effectors = autoplot(gr_effectors, aes(color = "green")),
  Breakpoints = autoplot(gr_bps, aes(color = "red"))
) +
  ggtitle("Chr9A: Genes, Effectors, and Breakpoints")



library(ggbio)
library(GenomicRanges)

# Convert to GRanges
gr_genes <- makeGRangesFromDataFrame(genes)
gr_effectors <- makeGRangesFromDataFrame(effectors)
gr_bps <- GRanges(seqnames = bps$chr, ranges = IRanges(start = bps$start, end = bps$end))

# Add metadata column for fixed y position (e.g., all at y = 0.5)
gr_genes$y <- 0.5
gr_effectors$y <- 0.5
gr_bps$y <- 0.5

# Genes track
p_genes <- autoplot(gr_genes, geom = "segment", aes(y = y, yend = y), color = "black", size = 10)

# Effectors track
p_effectors <- autoplot(gr_effectors, geom = "segment", aes(y = y, yend = y), color = "forestgreen", size = 10)

# Breakpoints track
p_bps <- autoplot(gr_bps, geom = "segment", aes(y = y, yend = y), color = "red", size = 10)

# Combine tracks
tracks(
  "Genes" = p_genes,
  "Effectors" = p_effectors,
  "Breakpoints" = p_bps,
  heights = c(1, 1, 1)
) + ggtitle("Chr9A: Genes, Effectors, and Breakpoints")


# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# Load gene BED
genes_df <- read_tsv("../560B_genes_at_breakpoints_2.5kb.bed", col_names = c("chr", "start", "end", "gene_id")) %>%
  mutate(y = "Genes", type = "Gene")

# Load effector BED
effectors_df <- read_tsv("../effector_genes_chr9B_by_suffix.bed", col_names = c("chr", "start", "end", "gene_id")) %>%
  mutate(y = "Effectors", type = "Effector")

# Load breakpoints BED (ensure non-zero length)
breakpoints_df <- read_tsv("../560B_approx_breakpoints.bed", col_names = c("chr", "start", "end")) %>%
  mutate(end = ifelse(end == start, start + 1, end),  # ensure width
         y = "Breakpoints",
         type = "Breakpoint")

# Combine all tracks
plot_df <- bind_rows(genes_df, effectors_df, breakpoints_df)

# Plot
ggplot(plot_df) +
  geom_segment(aes(x = start, xend = end, y = y, yend = y, color = type), size = 20) +
  scale_color_manual(
    name = "Feature",
    values = c("Gene" = "black", "Effector" = "forestgreen", "Breakpoint" = "red")
  ) +
  labs(
    title = "Chr9A: Genes, Effectors, Breakpoints",
    x = "Genomic Position (bp)",
    y = "Feature"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb"))


# Artificially widen features for better x-axis visibility
widen_factor <- 10000  # you can tune this value

plot_df_wide <- plot_df
plot_df_wide$start <- plot_df_wide$start - widen_factor / 2
plot_df_wide$end <- plot_df_wide$end + widen_factor / 2

# Plot with widened features (thicker on x-axis)
ggplot(plot_df_wide) +
  geom_segment(aes(x = start, xend = end, y = y, yend = y, color = type), linewidth = 10) +
  scale_color_manual(
    name = "Feature",
    values = c("Gene" = "#57B899", "Effector" = "#D771B6", "Breakpoint" = "#7B8DBE")
  ) +
  labs(
    title = "Chr9A: Genes, Effectors, Breakpoints",
    x = "Genomic Position (bp)",
    y = "Feature"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb"))
