library(ggplot2)
library(dplyr)
library(readr)

# Read in gene/effector/breakpoint data
genes <- read.table("../Ph560A_chr9A_genes.bed", header = FALSE, col.names = c("chr", "start", "end", "gene"))
effectors <- read.table("../effector_genes_chr9A_by_suffix.bed", header = FALSE, col.names = c("chr", "start", "end", "gene"))
bps <- read.table("../560A_approx_breakpoints.bed", header = FALSE, col.names = c("chr", "start", "end"))
# Load chrom index
fai <- read.table("../560A_chr9.fa.fai", header = FALSE, stringsAsFactors = FALSE)
colnames(fai)[1:2] <- c("chr", "length")

# Create genome info as GRanges object
custom_genome <- GRanges(seqnames = fai$chr, ranges = IRanges(start = 1, end = fai$length))
autoplot(custom_genome)

# Label
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

# Plottracks
tracks(
  Genes = autoplot(gr_genes, aes(color = "steelblue")),
  Effectors = autoplot(gr_effectors, aes(color = "green")),
  Breakpoints = autoplot(gr_bps, aes(color = "red"))
) +
  ggtitle("Chr9A genes, effectors and breakpoints")


library(ggbio)
library(GenomicRanges)

# Convert to GRanges
gr_genes <- makeGRangesFromDataFrame(genes)
gr_effectors <- makeGRangesFromDataFrame(effectors)
gr_bps <- GRanges(seqnames = bps$chr, ranges = IRanges(start = bps$start, end = bps$end))

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
) + ggtitle("Chr9A genes, effectors, and breakpoints")


# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# gene BED
genes_df <- read_tsv("../560B_genes_at_breakpoints_2.5kb.bed", col_names = c("chr", "start", "end", "gene_id")) %>%
  mutate(y = "Genes", type = "Gene")

# effector BED
effectors_df <- read_tsv("../effector_genes_chr9B_by_suffix.bed", col_names = c("chr", "start", "end", "gene_id")) %>%
  mutate(y = "Effectors", type = "Effector")

# breakpoints BED
breakpoints_df <- read_tsv("../560B_approx_breakpoints.bed", col_names = c("chr", "start", "end")) %>%
  mutate(end = ifelse(end == start, start + 1, end),  # ensure width
         y = "Breakpoints",
         type = "Breakpoint")

# Combine all tracks
plot_df <- bind_rows(genes_df, effectors_df, breakpoints_df)

ggplot(plot_df) +
  geom_segment(aes(x = start, xend = end, y = y, yend = y, color = type), size = 20) +
  scale_color_manual(
    name = "Feature",
    values = c("Gene" = "black", "Effector" = "forestgreen", "Breakpoint" = "red")
  ) +
  labs(
    title = "Chr9A genes, effectors, breakpoints",
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
