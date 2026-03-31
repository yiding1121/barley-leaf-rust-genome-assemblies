library(dplyr)
library(ggplot2)
library(tidyr)
library(GenomicRanges)

setwd("~/BLR_data/")

# === 1. Load chromosome sizes ===
chrom_sizes <- read.table("blr560A.chrom.sizes", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("chrom", "length")
chrom_sizes$length_Mb <- chrom_sizes$length / 1e6

# === 2. Load and count genes ===
gff <- read.table("Ph560A_braker.EVM.chr13Aupdated.gff3", header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(gff)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
genes <- gff %>% filter(type == "gene")

gene_counts <- genes %>%
  group_by(seqid) %>%
  summarise(gene_count = n(), .groups = "drop")

# === 3. Load and summarize non-overlapping TEs ===
te <- read.table("560A_TE_all.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(te) <- c("chrom", "start", "end", "type")
te <- te %>% mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start) & !is.na(end))

# Merge overlapping TE regions
te_gr <- reduce(GRanges(seqnames = te$chrom, IRanges(te$start, te$end)))

te_chr_length <- as.data.frame(te_gr) %>%
  group_by(seqnames) %>%
  summarise(total_te_bp = sum(width)) %>%
  rename(chrom = seqnames)

# Calculate % TE density
te_chr_df <- te_chr_length %>%
  left_join(chrom_sizes, by = "chrom") %>%
  mutate(te_density = (total_te_bp / length) * 100) %>%
  select(chrom, te_density) # Avoid conflicting columns

# === 4. Load GC content ===
gc_bins <- read.table("blr560a_gc_cov_100kb.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gc_content <- gc_bins %>%
  group_by(chrom) %>%
  summarise(gc = mean(GC, na.rm = TRUE), .groups = "drop")

# === 5. Load E1 compartment values ===
e1 <- read.table("./HiC/560A_cis_eigs_100kb.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
e1$E1 <- as.numeric(e1$E1)

e1_avg <- e1 %>%
  group_by(chrom) %>%
  summarise(E1_mean = mean(E1, na.rm = TRUE), .groups = "drop")

# === 6. Combine clearly into final dataframe ===
summary_df <- chrom_sizes %>%
  left_join(gene_counts, by = c("chrom" = "seqid")) %>%
  left_join(te_chr_df, by = "chrom") %>%
  left_join(gc_content, by = "chrom") %>%
  left_join(e1_avg, by = "chrom") %>%
  mutate(gene_density = gene_count / length_Mb)

# === 7. Prepare Data with explicit labels for Chr9A ===
summary_long <- summary_df %>%
  select(chrom, gene_density, te_density, gc, E1_mean) %>%
  pivot_longer(-chrom, names_to = "feature", values_to = "value") %>%
  mutate(highlight = ifelse(chrom == "Chr9A", "Chr9A", "Other"))

# === 8. Plot explicitly highlighting Chr9A (with annotations) ===
ggplot(summary_long, aes(x = chrom, y = value, fill = highlight)) +
  geom_col() +
  facet_wrap(~feature, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Chr9A" = "#FF6F61", "Other" = "gray80")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "none"
  ) +
  labs(
    title = "Chr9A is Gene Sparse and TE Enriched",
    subtitle = "Clearly highlighting chromosome Chr9A compared to others",
    x = "Chromosome",
    y = "Value"
  ) +
  geom_text(
    data = summary_long %>% filter(chrom == "Chr9A"),
    aes(label = round(value, 2)),
    position = position_stack(vjust = 1.1),
    size = 3.5,
    color = "black",
    fontface = "bold"
  )
