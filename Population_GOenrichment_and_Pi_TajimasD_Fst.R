tbl=read.table("../evalmix/BLR_all_hap1_pruned.3.Q")
tbl=read.table("../evalmix/BLR_all_hap1_pruned.4.table")

pop<-read.table("..evalmix/BLR_all_hap1_pruned.fam")

pop

barplot(t(as.matrix(Q4)), col=rainbow(3),xlab="Individual #", ylab="Ancestry", border=NA)

dev.off()


library(tidyverse)

plot_data <- tbl %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = V1:V3, names_to = "pop", values_to = "prob") %>%
  group_by(id) %>%
  mutate(
    likely_assignment = pop[which.max(prob)],
    assignment_prob = max(prob)
  ) %>%
  ungroup() %>%
  arrange(likely_assignment, desc(assignment_prob)) %>%
  mutate(id = factor(id, levels = unique(id)))  # THIS fixes the order in the plot

library(tidyverse)

# Read Q matrix (with 3 clusters, e.g., from ADMIXTURE K=3)
plot_data <- tbl %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = V1:V3, names_to = "pop", values_to = "prob") %>%
  group_by(id) %>%
  mutate(
    likely_assignment = pop[which.max(prob)],
    assignment_prob = max(prob)
  ) %>%
  ungroup() %>%
  arrange(likely_assignment, desc(assignment_prob)) %>%
  mutate(individual = factor(id, levels = unique(id)))  # <- this locks the plot order


ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B", "#499894"))

#or
#plot_data$id2 <- pop[ord,1]
matching <- data.frame(id1 = 1:41, id2 = pop[ord,1])
plot_data$id2 <- matching$id2[match(plot_data$id, matching$id1)]
ggplot(plot_data, aes(id2, prob, fill = pop)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free', space = 'free') +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B", "#499894",  "#D4A6C8")) +
  theme(axis.text.x = element_text(angle = 90))


# Second plotting method

source("/Users/yiding/Desktop/Scripts/evalAdmix-master/visFuns.R")


#cols_20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
             "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
             "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
             "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")


cols_20 <- c("#4E79A7", "#F28E2B", "#499894",  "#D4A6C8")


x <- tapply(1:nrow(pop),pop[ord,4],mean)
x <- x[as.vector(unique(pop[ord,4]))]

ord<-orderInds(pop = as.vector(pop[,1]), q = tbl)
barplot(t(tbl)[,ord],col=cols_20[1:4],space=0,border=NA,xlab="Individuals",ylab="Admixture proportions for K=4")
text(x,-0.05,as.vector(unique(pop[ord,1])),xpd=T, cex = 0.5, srt = 90)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)


ord<-orderInds(pop = as.vector(pop[,1]), q = tbl)
ord
mat <- t(tbl)[,ord]
colnames(mat) <- names(ord)
mat
barplot(mat,col=cols_20[1:4],space=0,border=NA,xlab="Individuals",ylab="Admixture proportions for K=4", las = 2, cex.names = 0.5)
text(tapply(1:nrow(pop),pop[ord,1],mean),-0.1,unique(pop[ord,1]),xpd=T)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)


#evalAdmix

r<-as.matrix(read.table("..evalmix/output.corres.txt"))

# Plot correlation of residuals
plotCorRes(cor_mat = r, pop = as.vector(pop[,2]), ord=ord, title="K=5", max_z=0.8, min_z=-0.8)

# from evalmix 2025

source("./Scripts/evalAdmix-master/visFuns.R")

# read population labels and estimated admixture proportions
pop<-read.table("..Hap1_based_pop_info.fam")
pop_hap2<-read.table("..evalmix/pruned_hap2_out.fam")

pop

q3<-read.table("..evalmix/BLR_all_hap1_pruned.3.Q",stringsAsFactors=T)
q4<-read.table("..evalmix/BLR_all_hap1_pruned.4.Q",stringsAsFactors=T)
q5<-read.table("..evalmix/BLR_all_hap1_pruned.5.Q",stringsAsFactors=T)
q6<-read.table("..evalmix/BLR_all_hap1_pruned.6.Q",stringsAsFactors=T)

palette(c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))

# order according to population and plot the ADMIXTURE reults
ord<-orderInds(pop = as.vector(pop[,2]), q = q)
ord_hap2<-orderInds(pop = as.vector(pop_hap2[,2]), q = q)
ord_hap2

#make barplot
plotAdmix(q3,ord=ord,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
plotAdmix(q4,ord=ord,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
plotAdmix(q5,ord=ord,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
plotAdmix(q5,ord=ord,pop=pop[,2], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
plotAdmix(q6,ord=ord,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A","#F1CE63"))
plotAdmix(q6,ord=ord,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A","#F1CE63"))

plotAdmix(q2_2,ord=ord_hap2,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
plotAdmix(q3_2,ord=ord_hap2,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
plotAdmix(q4_2,ord=ord_hap2,inds=pop[,1], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
plotAdmix(q4_2,ord=ord_hap2,pop=pop_hap2[,2], colorpal = c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))


r<-as.matrix(read.table("..evalmix/output.corres.4.txt"))

# Plot correlation of residuals
plotCorRes(cor_mat = r, pop = as.vector(pop[,2]), ord=ord, title="Evaluation of K=4", max_z=1, min_z=-1)


# Compute semantic similarity matrix


library(org.Sc.sgd.db)
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(readr)
library(GO.db)
library(AnnotationDbi)
library(viridis)


enrich <- read.csv("..Hap2_del_gProfiler_calbicans_GO_enriched.csv", header=TRUE)
enrich$logp <- -log10(enrich$adjusted_p_value)

gene2go <- read_table("..gene2go_hap2_deletion.txt", col_names = c("gene_id", "GO_term"))

go_gene_counts <- gene2go %>%
  filter(GO_term %in% enrich$term_id) %>%
  count(GO_term, name = "gene_count")

# `enrich$term_id` should be a vector of GO IDs
go_descriptions <- AnnotationDbi::select(
  GO.db,
  keys = enrich$term_id,
  columns = c("TERM"),
  keytype = "GOID"
)

# Join descriptions to enrichment results
enrich_named <- enrich %>%
  left_join(go_descriptions, by = c("term_id" = "GOID"))

# Combine with gene counts
go_data <- enrich_named %>%
  inner_join(go_gene_counts, by = c("term_id" = "GO_term"))

# Pull GO term names from GO.db
enrich$term_desc <- Term(enrich$term_id)

# Merge everything
go_data <- enrich %>%
  inner_join(go_gene_counts, by = c("term_id" = "GO_term"))

# Use semantic similarity (optional): here just build links if they share > x genes
go_links <- gene2go %>%
  filter(GO_term %in% go_data$term_id) %>%
  group_by(GO_term) %>%
  summarise(genes = list(unique(gene_id))) %>%
  tidyr::expand_grid(GO_term2 = .$GO_term) %>%
  filter(GO_term < GO_term2) %>%
  rowwise() %>%
  mutate(overlap = length(intersect(unlist(genes), unlist(go_data$genes[go_data$term_id == GO_term2])))) %>%
  filter(overlap >= 3) # Threshold to reduce noise

edges <- go_links %>%
  select(from = GO_term, to = GO_term2, weight = overlap)

nodes <- go_data %>%
  select(term_id, term_desc, logp, gene_count) %>%
  rename(name = term_id)

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.5, color = "gray50") +
  geom_node_point(aes(size = gene_count, color = logp)) +
  geom_node_text(aes(label = term_desc), repel = TRUE, size = 3) +
  scale_color_viridis_c(option = "plasma", name = "-log10(p adj)") +
  scale_size_continuous(name = "# genes") +
  theme_void() +
  labs(title = "GO Term Enrichment Network")

ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 1, color = "gray50") +
  geom_node_point(aes(size = gene_count, color = logp)) +
  geom_node_text(aes(label = term_desc), repel = TRUE, size = 3) +
  scale_color_viridis_c(option = "plasma", name = "-log10(p adj)") +
  scale_size_continuous(name = "# genes") +
  labs(
    title = "GO Term Enrichment Network",
    x = "X axis semantic",
    y = "Y axis semantic"
  ) +
  theme_bw(base_size = 14) + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(),
    axis.ticks = element_line(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) 


# Load libraries
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(GO.db)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(viridis)
library(purrr)

# 1. Load data
enrich <- read_csv("..Hap2_dup_gProfiler_calbicans_GO_enriched.csv")
gene2go <- read_table("..gene2go_hap2_duplication.txt", col_names = c("gene_id", "GO_term"))


# Get ontology info for each term
go_ontologies <- AnnotationDbi::select(
  GO.db,
  keys = enrich$term_id,
  columns = c("ONTOLOGY"),
  keytype = "GOID"
)

# Filter for MF and BP terms only
valid_ontologies <- go_ontologies %>%
  filter(ONTOLOGY %in% c("MF", "BP"))

# Subset enrich and related data
enrich_filtered <- enrich %>%
  semi_join(valid_ontologies, by = c("term_id" = "GOID"))

# Optional: check how many left
cat("Remaining GO terms after filtering:", nrow(enrich_filtered), "\n")


# 2. Calculate -log10(p adj)
str(enrich$adjusted_p_value)
enrich$adjusted_p_value <- as.numeric(enrich$adjusted_p_value)
enrich$logp <- -log10(enrich$adjusted_p_value)

# 3. Get GO descriptions
enrich$term_desc <- Term(enrich$term_id)

# 4. Count genes per GO term
go_gene_counts <- gene2go %>%
  filter(GO_term %in% enrich$term_id) %>%
  count(GO_term, name = "gene_count")


# 5. Merge enrichment data with gene counts
go_data <- enrich %>%
  inner_join(go_gene_counts, by = c("term_id" = "GO_term"))

# Step 1: Prepare gene sets per GO term
go_gene_list <- gene2go %>%
  filter(GO_term %in% go_data$term_id) %>%
  group_by(GO_term) %>%
  summarise(genes = list(unique(gene_id)), .groups = "drop")

# Step 2: Cross join with unique column names
pairwise_combos <- tidyr::crossing(go_gene_list, go_gene_list, .name_repair = "unique") %>%
  rename(
    GO_term1 = GO_term...1,
    genes1   = genes...2,
    GO_term2 = GO_term...3,
    genes2   = genes...4
  ) %>%
  filter(GO_term1 < GO_term2) %>%
  mutate(overlap = purrr::map2_int(genes1, genes2, ~length(intersect(.x, .y)))) %>%
  filter(overlap >= 3)

library(dplyr)

# Step 3: Edges (character IDs only!)
edges <- pairwise_combos %>%
  dplyr::select(from = GO_term1, to = GO_term2, weight = overlap)

# 1. Create the full node list (including those with no edges)
nodes <- go_data %>%
  dplyr::select(term_id, term_desc, logp, gene_count) %>%
  rename(name = term_id)

# 2. Create the edge list only for overlapping terms
edges <- pairwise_combos %>%
  select(from = GO_term1, to = GO_term2, weight = overlap)

# ✅ 3. Build graph using ALL nodes
graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(component = group_components())  # Optional: label disconnected nodes

ggraph(graph, layout = "fr") +
  geom_node_point(aes(size = gene_count, color = logp)) +
  geom_node_text(
    aes(label = term_desc),
    repel = TRUE,
    size = 3,
    color = "black",
    fontface = "bold"
  ) +
  scale_color_gradient(
    low = "#7B8DBF",  # Blue
    high = "#57B893", # Green
    name = "-log10(p adj)"
  ) +
  scale_size_continuous(name = "# genes", range = c(5, 16)) +
  labs(
    title = "GO Term Enrichment Network",
    x = "Semantic X-axis",
    y = "Semantic Y-axis"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

###############################

library(tidyverse)

# Read data
pi_pop12 <- read_tsv("..pop1_2_hap1.windowed.pi", show_col_types = FALSE) %>%
  mutate(pop = "Population 1-2")

pi_pop345 <- read_tsv("..pop3_4_5_hap1.windowed.pi", show_col_types = FALSE) %>%
  mutate(pop = "Population 3-5")

pi_data <- bind_rows(pi_pop12, pi_pop345)

# -----------------------------
# 2. Clean and format
# -----------------------------
pi_data <- pi_data %>%
  filter(!is.na(PI)) %>%
  mutate(
    mid = (BIN_START + BIN_END) / 2,
    mid_Mb = mid / 1e6   # convert to Mb
  )

# Order chromosomes
chr_order <- pi_data %>%
  distinct(CHROM) %>%
  mutate(chr_num = as.numeric(str_extract(CHROM, "\\d+"))) %>%
  arrange(chr_num) %>%
  pull(CHROM)

pi_data$CHROM <- factor(pi_data$CHROM, levels = chr_order)

# -----------------------------
# 3. Plot
# -----------------------------
p <- ggplot(pi_data, aes(x = mid_Mb, y = PI)) +
  geom_line(linewidth = 0.4) +
  facet_grid(pop ~ CHROM, scales = "free_x") +
  labs(
    x = "Position (Mb)",
    y = expression(paste("Nucleotide diversity (", pi, ")")),
    title = "Chromosome-level genetic diversity"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

p <- ggplot(pi_data, aes(x = mid_Mb, y = PI, colour = pop)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ CHROM, scales = "free_x", ncol = 6) +
  labs(
    x = "Position (Mb)",
    y = expression(paste("Nucleotide diversity (", pi, ")")),
    colour = "Population"
  ) +
  theme_bw()

print(p)

library(tidyverse)

# Read data
pi_pop12 <- read_tsv("..pop1_2_hap1.windowed.pi", show_col_types = FALSE) %>%
  mutate(pop = "Population 1-2")

pi_pop345 <- read_tsv("..pop3_4_5_hap1.windowed.pi", show_col_types = FALSE) %>%
  mutate(pop = "Population 3-5")

pi_data <- bind_rows(pi_pop12, pi_pop345) %>%
  filter(!is.na(PI))

# -----------------------------
# 2. Clean data
# -----------------------------
pi_data <- pi_data %>%
  mutate(
    PI = as.numeric(PI),
    BIN_START = as.numeric(BIN_START),
    BIN_END = as.numeric(BIN_END)
  ) %>%
  filter(
    !is.na(CHROM),
    !is.na(PI),
    !is.na(BIN_START),
    !is.na(BIN_END)
  ) %>%
  mutate(
    mid_bp = (BIN_START + BIN_END) / 2,
    mid_Mb = mid_bp / 1e6
  )

# -----------------------------
# 3. Order chromosomes
# -----------------------------
chr_order <- pi_data %>%
  distinct(CHROM) %>%
  mutate(
    chr_num = as.numeric(str_extract(CHROM, "\\d+")),
    chr_suffix = str_extract(CHROM, "[A-Za-z]+$")
  ) %>%
  arrange(chr_num, chr_suffix) %>%
  pull(CHROM)

pi_data <- pi_data %>%
  mutate(CHROM = factor(CHROM, levels = chr_order))

# -----------------------------
# 4. Chromosome-level faceted plot
#    both populations shown, x-axis in Mb
# -----------------------------
p_faceted <- ggplot(pi_data, aes(x = mid_Mb, y = PI, colour = pop, group = pop)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ CHROM, scales = "free_x", ncol = 6) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(
    x = "Position (Mb)",
    y = expression(paste("Nucleotide diversity (", pi, ")")),
    colour = "Population",
    title = "Chromosome-level nucleotide diversity"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey85", colour = "grey40"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_faceted)

# -----------------------------
# 10. chr_9A vs all other chromosomes in Population 1-2
# -----------------------------
df_test_pop12 <- pi_data %>%
  filter(pop == "Population 1-2") %>%
  mutate(
    target = ifelse(as.character(CHROM) == "chr_9A", "chr_9A", "Other"),
    target = factor(target, levels = c("chr_9A", "Other"))
  )

# Summary
test_summary_pop12 <- df_test_pop12 %>%
  group_by(target) %>%
  summarise(
    mean_pi = mean(PI, na.rm = TRUE),
    median_pi = median(PI, na.rm = TRUE),
    sd_pi = sd(PI, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

print(test_summary_pop12)
write_tsv(test_summary_pop12,
          "..population_1_2_chr9A_vs_other_summary.tsv")

# Wilcoxon test
wilcox_res_pop12 <- wilcox.test(
  x = df_test_pop12$PI[df_test_pop12$target == "chr_9A"],
  y = df_test_pop12$PI[df_test_pop12$target == "Other"],
  alternative = "greater"
)

print(wilcox_res_pop12)
capture.output(
  wilcox_res_pop12,
  file = "..population_1_2_chr9A_vs_other_wilcox.txt"
)

p_value_pop12 <- wilcox_res_pop12$p.value

p_star_pop12 <- case_when(
  p_value_pop12 < 0.001 ~ "***",
  p_value_pop12 < 0.01  ~ "**",
  p_value_pop12 < 0.05  ~ "*",
  TRUE ~ "ns"
)

p_label_pop12 <- paste0("p = ", signif(p_value_pop12, 3), " ", p_star_pop12)

y_min_pop12 <- min(df_test_pop12$PI, na.rm = TRUE)
y_max_pop12 <- max(df_test_pop12$PI, na.rm = TRUE)
y_bracket_pop12 <- y_max_pop12 * 1.05
y_tick_pop12 <- y_max_pop12 * 1.03
y_text_pop12 <- y_max_pop12 * 1.10

p_box_pop12 <- ggplot(df_test_pop12, aes(x = target, y = PI, fill = target)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  scale_fill_manual(values = c("chr_9A" = "red", "Other" = "grey70")) +
  annotate("segment", x = 1, xend = 2, y = y_bracket_pop12, yend = y_bracket_pop12) +
  annotate("segment", x = 1, xend = 1, y = y_tick_pop12, yend = y_bracket_pop12) +
  annotate("segment", x = 2, xend = 2, y = y_tick_pop12, yend = y_bracket_pop12) +
  annotate("text", x = 1.5, y = y_text_pop12, label = p_label_pop12, size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(
    x = NULL,
    y = expression(paste("Windowed nucleotide diversity (", pi, ")")),
    title = "Population 1-2: chr_9A versus all other chromosomes"
  ) +
  coord_cartesian(ylim = c(y_min_pop12, y_max_pop12 * 1.15)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_box_pop12)


# -----------------------------
# 5. Summarize chromosome-wide pi
# -----------------------------
chr_summary <- pi_data %>%
  group_by(pop, CHROM) %>%
  summarise(
    mean_pi = mean(PI, na.rm = TRUE),
    median_pi = median(PI, na.rm = TRUE),
    max_pi = max(PI, na.rm = TRUE),
    sd_pi = sd(PI, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

# -----------------------------
# 6. Rank chromosomes within Population 3-5
# -----------------------------
pop35_rank <- chr_summary %>%
  filter(pop == "Population 3-5") %>%
  arrange(desc(mean_pi)) %>%
  mutate(rank = row_number())

print(pop35_rank)

# Show top chromosome in Population 3-5
top_chr_pop35 <- pop35_rank %>% slice(1)
print(top_chr_pop35)

# -----------------------------
# 7. Barplot of chromosome-wide mean pi
#    for Population 3-5 only
# -----------------------------
pop35_plot_data <- pop35_rank %>%
  mutate(CHROM = fct_reorder(as.character(CHROM), mean_pi, .desc = TRUE))

p_bar_pop35 <- ggplot(pop35_plot_data, aes(x = CHROM, y = mean_pi, fill = CHROM == "chr_9A")) +
  geom_col(width = 0.8) +
  geom_text(aes(label = sprintf("%.5f", mean_pi)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  labs(
    x = "Chromosome",
    y = expression(paste("Mean nucleotide diversity (", pi, ")")),
    title = "Chromosome-wide nucleotide diversity in Population 3-5"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar_pop35)

# -----------------------------
# 8. Statistical comparison:
#    chr_9A vs all other chromosomes in Population 3-5
# -----------------------------
df_test <- pi_data %>%
  filter(pop == "Population 3-5") %>%
  mutate(target = ifelse(CHROM == "chr_9A", "chr_9A", "Other"))

# Check summary first
test_summary <- df_test %>%
  group_by(target) %>%
  summarise(
    mean_pi = mean(PI, na.rm = TRUE),
    median_pi = median(PI, na.rm = TRUE),
    sd_pi = sd(PI, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

print(test_summary)

# Wilcoxon test
wilcox_res <- wilcox.test(
  x = df_test$PI[df_test$target == "chr_9A"],
  y = df_test$PI[df_test$target == "Other"],
  alternative = "greater"
)

print(wilcox_res)

# -----------------------------
# 9. Boxplot: chr_9A vs other chromosomes in Population 3-5
# -----------------------------
# Make sure df_test has the target column
df_test <- pi_data %>%
  filter(pop == "Population 3-5") %>%
  mutate(
    target = ifelse(as.character(CHROM) == "chr_9A", "chr_9A", "Other"),
    target = factor(target, levels = c("chr_9A", "Other"))
  )

# Wilcoxon test
wilcox_res <- wilcox.test(
  x = df_test$PI[df_test$target == "chr_9A"],
  y = df_test$PI[df_test$target == "Other"],
  alternative = "greater"
)

p_value <- wilcox_res$p.value

# Convert to stars
p_star <- case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  TRUE ~ "ns"
)

# Label text
p_label <- paste0("p = ", signif(p_value, 3), " ", p_star)

y_max <- max(df_test$PI, na.rm = TRUE)

p_box <- ggplot(df_test, aes(x = target, y = PI, fill = target)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  scale_fill_manual(values = c("chr_9A" = "red", "Other" = "grey70")) +
  
  # Add significance bar
  geom_segment(aes(x = 1, xend = 2, y = y_max*1.05, yend = y_max*1.05), inherit.aes = FALSE) +
  geom_segment(aes(x = 1, xend = 1, y = y_max*1.03, yend = y_max*1.05), inherit.aes = FALSE) +
  geom_segment(aes(x = 2, xend = 2, y = y_max*1.03, yend = y_max*1.05), inherit.aes = FALSE) +
  
  # Add p-value + stars
  annotate(
    "text",
    x = 1.5,
    y = y_max*1.08,
    label = p_label,
    size = 5
  ) +
  
  labs(
    x = NULL,
    y = expression(paste("Windowed nucleotide diversity (", pi, ")")),
    title = "Population 3-5: chr_9A versus all other chromosomes"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none"
  )

print(p_box)


##############################Tajima's D

library(tidyverse)
library(forcats)
library(patchwork)

# -----------------------------
# 1. Read Tajima's D files
# -----------------------------
taj_pop12 <- read_tsv("..tajimasd_pop12_10kb.Tajima.D",
                      show_col_types = FALSE) %>%
  mutate(pop = "Population 1-2")

taj_pop35 <- read_tsv("..tajimasd_pop345_10kb.Tajima.D",
                      show_col_types = FALSE) %>%
  mutate(pop = "Population 3-5")

taj_data <- bind_rows(taj_pop12, taj_pop35)

# -----------------------------
# 2. Clean data
# -----------------------------
taj_data <- taj_data %>%
  mutate(
    BIN_START = as.numeric(BIN_START),
    N_SNPS = as.numeric(N_SNPS),
    TajimaD = as.numeric(TajimaD)
  ) %>%
  filter(
    !is.na(CHROM),
    !is.na(BIN_START),
    !is.na(TajimaD)
  ) %>%
  mutate(
    mid_bp = BIN_START + 50000,
    mid_Mb = mid_bp / 1e6
  )

# -----------------------------
# 3. Natural chromosome order
# -----------------------------
chr_order <- taj_data %>%
  distinct(CHROM) %>%
  mutate(
    chr_num = as.numeric(stringr::str_extract(CHROM, "\\d+")),
    chr_suffix = stringr::str_extract(CHROM, "[A-Za-z]+$")
  ) %>%
  arrange(chr_num, chr_suffix) %>%
  pull(CHROM)

taj_data <- taj_data %>%
  mutate(CHROM = factor(CHROM, levels = chr_order))

# -----------------------------
# 4. Chromosome summaries
# -----------------------------
taj_summary <- taj_data %>%
  group_by(pop, CHROM) %>%
  summarise(
    mean_tajD   = mean(TajimaD, na.rm = TRUE),
    median_tajD = median(TajimaD, na.rm = TRUE),
    sd_tajD     = sd(TajimaD, na.rm = TRUE),
    n_windows   = n(),
    se_tajD     = sd_tajD / sqrt(n_windows),
    prop_gt2    = mean(TajimaD > 2, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 5. Panel A:
#    Ranked chromosome medians by population
# -----------------------------
rank_plot_data <- taj_summary %>%
  group_by(pop) %>%
  arrange(desc(median_tajD), .by_group = TRUE) %>%
  mutate(rank = row_number(),
         highlight = ifelse(as.character(CHROM) == "chr_9A", "chr_9A", "Other")) %>%
  ungroup()

pA <- ggplot(rank_plot_data, aes(x = rank, y = median_tajD)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(fill = highlight), shape = 21, size = 4, colour = "black", stroke = 0.4) +
  geom_text(aes(label = CHROM), nudge_y = 0.08, size = 3) +
  scale_fill_manual(values = c("Other" = "grey75", "chr_9A" = "red")) +
  facet_wrap(~ pop, ncol = 1, scales = "free_y") +
  labs(
    x = "Chromosome rank within population",
    y = "Median Tajima's D",
    fill = NULL,
    title = "A. Chromosome ranking by Tajima's D"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 6. Panel B:
#    Windowed Tajima's D along chr_9A
# -----------------------------
chr9_data <- taj_data %>%
  filter(as.character(CHROM) == "chr_9A")

pB <- ggplot(chr9_data, aes(x = mid_Mb, y = TajimaD, colour = pop)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) +
  labs(
    x = "Position on chr_9A (Mb)",
    y = "Tajima's D",
    colour = "Population",
    title = "B. Windowed Tajima's D along chr_9A"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 7. Panel C:
#    Population shift per chromosome
# -----------------------------
taj_diff <- taj_summary %>%
  select(pop, CHROM, median_tajD) %>%
  pivot_wider(names_from = pop, values_from = median_tajD) %>%
  mutate(
    delta = `Population 3-5` - `Population 1-2`,
    highlight = ifelse(as.character(CHROM) == "chr_9A", "chr_9A", "Other")
  ) %>%
  arrange(desc(delta)) %>%
  mutate(CHROM = factor(CHROM, levels = CHROM))

pC <- ggplot(taj_diff, aes(x = CHROM, y = delta, fill = highlight)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_col() +
  scale_fill_manual(values = c("Other" = "grey75", "chr_9A" = "red")) +
  labs(
    x = "Chromosome",
    y = expression(Delta*" median Tajima's D (Pop 3-5 - Pop 1-2)"),
    fill = NULL,
    title = "C. Population shift in Tajima's D by chromosome"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 8. Combine
# -----------------------------
final_plot <- pA / pB / pC +
  plot_layout(heights = c(1.2, 1, 1)) +
  plot_annotation(
    title = "Tajima's D places chr_9A among elevated chromosomes, but not as the strongest population-shift signal"
  )

print(final_plot)

taj_wide <- taj_data %>%
  select(CHROM, BIN_START, pop, TajimaD) %>%
  pivot_wider(names_from = pop, values_from = TajimaD) %>%
  mutate(abs_diff = abs(`Population 3-5` - `Population 1-2`))

abs_diff_summary <- taj_wide %>%
  group_by(CHROM) %>%
  summarise(
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    median_abs_diff = median(abs_diff, na.rm = TRUE),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_diff)) %>%
  mutate(highlight = ifelse(CHROM == "chr_9A", "chr_9A", "Other"),
         CHROM = factor(CHROM, levels = CHROM))

ggplot(abs_diff_summary, aes(x = CHROM, y = mean_abs_diff, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Other" = "grey75", "chr_9A" = "red")) +
  labs(
    x = "Chromosome",
    y = "Mean absolute difference in Tajima's D",
    fill = NULL,
    title = "Population difference in Tajima's D by chromosome"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

shift_variation <- taj_wide %>%
  group_by(CHROM) %>%
  summarise(
    sd_diff = sd(`Population 3-5` - `Population 1-2`, na.rm = TRUE),
    iqr_diff = IQR(`Population 3-5` - `Population 1-2`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sd_diff)) %>%
  mutate(highlight = ifelse(CHROM == "chr_9A", "chr_9A", "Other"),
         CHROM = factor(CHROM, levels = CHROM))

ggplot(shift_variation, aes(x = CHROM, y = sd_diff, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Other" = "grey75", "chr_9A" = "red")) +
  labs(
    x = "Chromosome",
    y = "SD of population difference in Tajima's D",
    fill = NULL,
    title = "Heterogeneity of Tajima's D shift by chromosome"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

chr9_diff <- taj_wide %>%
  filter(CHROM == "chr_9A") %>%
  mutate(
    mid_Mb = (BIN_START + 50000) / 1e6,
    diff = `Population 3-5` - `Population 1-2`
  )

ggplot(chr9_diff, aes(x = mid_Mb, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.9, colour = "red") +
  labs(
    x = "Position on chr_9A (Mb)",
    y = "Tajima's D difference (Pop 3-5 - Pop 1-2)",
    title = "Windowed Tajima's D shift along chr_9A"
  ) +
  theme_bw()


# -----------------------------
# Prepare wide format
# -----------------------------
taj_wide <- taj_data %>%
  select(CHROM, BIN_START, pop, TajimaD) %>%
  pivot_wider(names_from = pop, values_from = TajimaD)

# -----------------------------
# Compute difference
# -----------------------------
taj_diff_all <- taj_wide %>%
  mutate(
    mid_Mb = (BIN_START + 50000) / 1e6,
    diff = `Population 3-5` - `Population 1-2`
  )

# -----------------------------
# Plot (all chromosomes)
# -----------------------------
p_diff_genome <- ggplot(taj_diff_all, aes(x = mid_Mb, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.6, colour = "#7B8DBF") +

  
  
#palette(c("#7B8DBF", "#F87850", "#57B893","#D771B6","#B2DF8A"))
  
  # highlight chr_9A
  geom_line(
    data = subset(taj_diff_all, CHROM == "chr_9A"),
    linewidth = 0.8,
    colour = "#57B893"
  ) +
  
  facet_wrap(~ CHROM, scales = "free_x", ncol = 6) +
  
  labs(
    x = "Position (Mb)",
    y = expression(Delta*" Tajima's D (Pop 3-5 - Pop 1-2)"),
    title = "Genome-wide windowed Tajima's D difference"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_diff_genome)

library(tidyverse)
library(patchwork)

# -----------------------------
# 1. Read Tajima's D files
# -----------------------------
taj_pop12 <- read_tsv(
  "..tajimasd_pop12.Tajima.D",
  show_col_types = FALSE
) %>%
  mutate(pop = "Population 1-2")

taj_pop35 <- read_tsv(
  "..tajimasd_pop345.Tajima.D",
  show_col_types = FALSE
) %>%
  mutate(pop = "Population 3-5")

taj_data <- bind_rows(taj_pop12, taj_pop35)

# -----------------------------
# 2. Clean data
# -----------------------------
taj_data <- taj_data %>%
  mutate(
    BIN_START = as.numeric(BIN_START),
    N_SNPS = as.numeric(N_SNPS),
    TajimaD = as.numeric(TajimaD)
  ) %>%
  filter(
    !is.na(CHROM),
    !is.na(BIN_START),
    !is.na(TajimaD)
  ) %>%
  mutate(
    mid_bp = BIN_START + 50000,   # 100 kb windows
    mid_Mb = mid_bp / 1e6
  )

# -----------------------------
# 3. Chromosome order
# -----------------------------
chr_order <- taj_data %>%
  distinct(CHROM) %>%
  mutate(
    chr_num = as.numeric(stringr::str_extract(CHROM, "\\d+")),
    chr_suffix = stringr::str_extract(CHROM, "[A-Za-z]+$")
  ) %>%
  arrange(chr_num, chr_suffix) %>%
  pull(CHROM)

taj_data <- taj_data %>%
  mutate(CHROM = factor(CHROM, levels = chr_order))

# -----------------------------
# 4. Wide table for delta Tajima's D
# -----------------------------
taj_wide <- taj_data %>%
  select(CHROM, BIN_START, mid_Mb, pop, TajimaD) %>%
  pivot_wider(names_from = pop, values_from = TajimaD) %>%
  mutate(
    diff = `Population 3-5` - `Population 1-2`
  )

# -----------------------------
# 5. Optional smoothing for cleaner lines
# -----------------------------
# Change k if you want more/less smoothing
k_smooth <- 5

taj_wide <- taj_wide %>%
  group_by(CHROM) %>%
  arrange(BIN_START, .by_group = TRUE) %>%
  mutate(
    diff_smooth = zoo::rollmean(diff, k = k_smooth, fill = NA, align = "center")
  ) %>%
  ungroup()

taj_data_smooth <- taj_data %>%
  group_by(CHROM, pop) %>%
  arrange(BIN_START, .by_group = TRUE) %>%
  mutate(
    TajimaD_smooth = zoo::rollmean(TajimaD, k = k_smooth, fill = NA, align = "center")
  ) %>%
  ungroup()

# -----------------------------
# 6. Color palette matching manuscript style
# -----------------------------
pop_colors <- c(
  "Population 1-2" = "#E76F51",
  "Population 3-5" = "#2A9D8F"
)

chr9_color <- "#D62828"
bg_line_color <- "grey45"
bg_facet_fill <- "grey88"

# -----------------------------
# 7. Panel A: genome-wide delta Tajima's D
# -----------------------------
pA <- ggplot(taj_wide, aes(x = mid_Mb, y = diff_smooth)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_line(colour = "grey55", linewidth = 0.55, na.rm = TRUE) +
  geom_line(
    data = subset(taj_wide, CHROM == "chr_9A"),
    colour = chr9_color,
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  facet_wrap(~ CHROM, scales = "free_x", ncol = 6) +
  labs(
    x = "Position (Mb)",
    y = expression(Delta * " Tajima's D (Pop 3-5 - Pop 1-2)"),
    title = "A. Genome-wide windowed Tajima's D difference"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = bg_facet_fill, colour = "grey30"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# 8. chr_9A data only
# -----------------------------
chr9_data <- taj_data_smooth %>%
  filter(CHROM == "chr_9A")

chr9_diff <- taj_wide %>%
  filter(CHROM == "chr_9A")

# -----------------------------
# 9. Panel B: chr_9A Tajima's D by population
# -----------------------------
pB <- ggplot(chr9_data, aes(x = mid_Mb, y = TajimaD_smooth, colour = pop)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 1.0, na.rm = TRUE) +
  scale_colour_manual(values = pop_colors) +
  labs(
    x = "Position on chr_9A (Mb)",
    y = "Tajima's D",
    colour = "Population",
    title = "B. Windowed Tajima's D along chr_9A"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# 10. Panel C: chr_9A delta Tajima's D
# -----------------------------
pC <- ggplot(chr9_diff, aes(x = mid_Mb, y = diff_smooth)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_line(colour = chr9_color, linewidth = 1.0, na.rm = TRUE) +
  labs(
    x = "Position on chr_9A (Mb)",
    y = expression(Delta * " Tajima's D"),
    title = "C. Localized Tajima's D shift along chr_9A"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# 11. OPTIONAL Panel D template for SV/CNV/breakpoints
# Replace this with your real track later
# -----------------------------
# Example placeholder:
# sv_track <- your_sv_window_table %>%
#   filter(CHROM == "chr_9A")
#
# pD <- ggplot(sv_track, aes(x = mid_Mb, y = sv_density)) +
#   geom_col(fill = "#4C78A8", width = 0.08) +
#   labs(
#     x = "Position on chr_9A (Mb)",
#     y = "SV density",
#     title = "D. Structural variation along chr_9A"
#   ) +
#   theme_bw(base_size = 13) +
#   theme(
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(face = "bold")
#   )

# -----------------------------
# 12. Combine figure
# -----------------------------
final_plot <- pA / pB / pC +
  plot_layout(heights = c(2.2, 1, 1)) +
  plot_annotation(
    title = "Localized Tajima's D divergence highlights a hotspot on chr_9A"
  )

print(final_plot)


# -----------------------------
# Prepare wide format
# -----------------------------
taj_wide <- taj_data %>%
  select(CHROM, BIN_START, pop, TajimaD) %>%
  pivot_wider(names_from = pop, values_from = TajimaD)

# -----------------------------
# Compute difference
# -----------------------------
taj_diff_all <- taj_wide %>%
  mutate(
    mid_Mb = (BIN_START + 50000) / 1e6,
    diff = `Population 3-5` - `Population 1-2`
  )

# -----------------------------
# Prepare wide format
# -----------------------------
taj_wide <- taj_data %>%
  select(CHROM, BIN_START, pop, TajimaD) %>%
  pivot_wider(names_from = pop, values_from = TajimaD)

# -----------------------------
# Compute difference
# -----------------------------
taj_diff_all <- taj_wide %>%
  mutate(
    mid_Mb = (BIN_START + 50000) / 1e6,
    diff = `Population 3-5` - `Population 1-2`
  )

taj_wide <- taj_data %>%
  select(CHROM, BIN_START, pop, TajimaD) %>%
  pivot_wider(names_from = pop, values_from = TajimaD) %>%
  mutate(abs_diff = abs(`Population 3-5` - `Population 1-2`))

abs_diff_summary <- taj_wide %>%
  group_by(CHROM) %>%
  summarise(
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    median_abs_diff = median(abs_diff, na.rm = TRUE),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_diff)) %>%
  mutate(highlight = ifelse(CHROM == "chr_9A", "chr_9A", "Other"),
         CHROM = factor(CHROM, levels = CHROM))

ggplot(abs_diff_summary, aes(x = CHROM, y = mean_abs_diff, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Other" = "grey75", "chr_9A" = "red")) +
  labs(
    x = "Chromosome",
    y = "Mean absolute difference in Tajima's D",
    fill = NULL,
    title = "Population difference in Tajima's D by chromosome"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


# -----------------------------
# 10. Boxplot: chr_9A vs Other in Population 3-5
# -----------------------------
y_min_taj <- min(df_taj_test$TajimaD, na.rm = TRUE)
y_max_taj <- max(df_taj_test$TajimaD, na.rm = TRUE)
y_bracket_taj <- y_max_taj * 1.05
y_tick_taj <- y_max_taj * 1.02
y_text_taj <- y_max_taj * 1.10

p_taj_box <- ggplot(df_taj_test, aes(x = target, y = TajimaD, fill = target)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  scale_fill_manual(values = c("chr_9A" = "red", "Other" = "grey70")) +
  annotate("segment", x = 1, xend = 2, y = y_bracket_taj, yend = y_bracket_taj) +
  annotate("segment", x = 1, xend = 1, y = y_tick_taj, yend = y_bracket_taj) +
  annotate("segment", x = 2, xend = 2, y = y_tick_taj, yend = y_bracket_taj) +
  annotate("text", x = 1.5, y = y_text_taj, label = p_label_taj, size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(
    x = NULL,
    y = "Windowed Tajima's D",
    title = "Population 3-5: chr_9A versus all other chromosomes"
  ) +
  coord_cartesian(ylim = c(y_min_taj, y_max_taj * 1.15)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none"
  )

print(p_taj_box)


df_taj_test_pop12 <- taj_data %>%
  filter(pop == "Population 1-2") %>%
  mutate(
    target = ifelse(as.character(CHROM) == "chr_9A", "chr_9A", "Other"),
    target = factor(target, levels = c("chr_9A", "Other"))
  )

wilcox_taj_pop12 <- wilcox.test(
  x = df_taj_test_pop12$TajimaD[df_taj_test_pop12$target == "chr_9A"],
  y = df_taj_test_pop12$TajimaD[df_taj_test_pop12$target == "Other"],
  alternative = "two.sided"
)

p_value_taj_pop12 <- wilcox_taj_pop12$p.value

p_star_taj_pop12 <- case_when(
  p_value_taj_pop12 < 0.001 ~ "***",
  p_value_taj_pop12 < 0.01 ~ "**",
  p_value_taj_pop12 < 0.05 ~ "*",
  TRUE ~ "ns"
)

p_label_taj_pop12 <- paste0("p = ", signif(p_value_taj_pop12, 3), " ", p_star_taj_pop12)

y_min_taj12 <- min(df_taj_test_pop12$TajimaD, na.rm = TRUE)
y_max_taj12 <- max(df_taj_test_pop12$TajimaD, na.rm = TRUE)

p_taj_box_pop12 <- ggplot(df_taj_test_pop12, aes(x = target, y = TajimaD, fill = target)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  scale_fill_manual(values = c("chr_9A" = "red", "Other" = "grey70")) +
  
  # significance bracket
  annotate("segment", x = 1, xend = 2, y = y_max_taj12*1.05, yend = y_max_taj12*1.05) +
  annotate("segment", x = 1, xend = 1, y = y_max_taj12*1.03, yend = y_max_taj12*1.05) +
  annotate("segment", x = 2, xend = 2, y = y_max_taj12*1.03, yend = y_max_taj12*1.05) +
  
  # p-value
  annotate("text", x = 1.5, y = y_max_taj12*1.10, label = p_label_taj_pop12, size = 5) +
  
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  
  labs(
    x = NULL,
    y = "Windowed Tajima's D",
    title = "Population 1-2: chr_9A versus all other chromosomes"
  ) +
  
  coord_cartesian(ylim = c(y_min_taj12, y_max_taj12 * 1.15)) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none"
  )

print(p_taj_box_pop12)



#########################Fst analysis
library(tidyverse)

fst <- read_tsv("..pop12_345_hap1_fst.windowed.weir.fst", show_col_types = FALSE)

fst <- fst %>%
  mutate(
    WEIGHTED_FST = as.numeric(WEIGHTED_FST),
    BIN_START = as.numeric(BIN_START),
    BIN_END = as.numeric(BIN_END)
  ) %>%
  filter(!is.na(WEIGHTED_FST)) %>%
  mutate(
    mid_bp = (BIN_START + BIN_END) / 2,
    mid_Mb = mid_bp / 1e6
  )

chr_order <- fst %>%
  distinct(CHROM) %>%
  mutate(
    chr_num = as.numeric(str_extract(CHROM, "\\d+")),
    chr_suffix = str_extract(CHROM, "[A-Za-z]+$")
  ) %>%
  arrange(chr_num, chr_suffix) %>%
  pull(CHROM)

fst <- fst %>%
  mutate(CHROM = factor(CHROM, levels = chr_order))

p_fst <- ggplot(fst, aes(x = mid_Mb, y = WEIGHTED_FST)) +
  geom_line(color = "grey40", linewidth = 0.5) +
  
  # Highlight chr9A
  geom_line(
    data = fst %>% filter(CHROM == "chr_9A"),
    aes(x = mid_Mb, y = WEIGHTED_FST),
    color = "#d73027", linewidth = 0.7
  ) +
  
  facet_wrap(~ CHROM, scales = "free_x", ncol = 6) +
  
  labs(
    x = "Position (Mb)",
    y = "FST (Pop 1–2 vs Pop 3–5)",
    title = "Genome-wide population differentiation highlights chr_9A"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_fst)

fst_summary <- fst %>%
  group_by(CHROM) %>%
  summarise(
    mean_fst = mean(WEIGHTED_FST, na.rm = TRUE),
    median_fst = median(WEIGHTED_FST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_fst))

print(fst_summary)

fst_summary <- fst_summary %>%
  mutate(CHROM = fct_reorder(CHROM, mean_fst, .desc = TRUE))

p_bar_fst <- ggplot(fst_summary, aes(x = CHROM, y = mean_fst,
                                     fill = CHROM == "chr_9A")) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey70")) +
  labs(
    x = "Chromosome",
    y = "Mean FST",
    title = "Chromosome-wide differentiation between populations"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar_fst)


library(tidyverse)
library(forcats)
library(patchwork)

# -----------------------------
# 1. Read FST files
# -----------------------------
fst_hap1 <- read_tsv("..pop12_345_hap1_fst.windowed.weir.fst", show_col_types = FALSE) %>%
  mutate(haplotype = "hap1")

fst_hap2 <- read_tsv("..pop12_345_hap2_fst.windowed.weir.fst", show_col_types = FALSE) %>%
  mutate(haplotype = "hap2")

fst_all <- bind_rows(fst_hap1, fst_hap2)

# -----------------------------
# 2. Clean
# -----------------------------
fst_all <- fst_all %>%
  mutate(
    WEIGHTED_FST = as.numeric(WEIGHTED_FST),
    BIN_START = as.numeric(BIN_START),
    BIN_END = as.numeric(BIN_END)
  ) %>%
  filter(
    !is.na(CHROM),
    !is.na(WEIGHTED_FST),
    !is.na(BIN_START),
    !is.na(BIN_END)
  ) %>%
  mutate(
    mid_bp = (BIN_START + BIN_END) / 2,
    mid_Mb = mid_bp / 1e6
  )

# -----------------------------
# 3. Chromosome order
# -----------------------------
chr_order <- fst_all %>%
  distinct(CHROM) %>%
  mutate(
    chr_num = as.numeric(str_extract(CHROM, "\\d+")),
    chr_suffix = str_extract(CHROM, "[A-Za-z]+$")
  ) %>%
  arrange(chr_num, chr_suffix) %>%
  pull(CHROM)

fst_all <- fst_all %>%
  mutate(CHROM = factor(CHROM, levels = chr_order))

# -----------------------------
# 4. Peak on chr_9A for each haplotype
# -----------------------------
chr9_peaks <- fst_all %>%
  filter(CHROM == "chr_9A") %>%
  group_by(haplotype) %>%
  slice_max(WEIGHTED_FST, n = 1, with_ties = FALSE) %>%
  ungroup()

print(chr9_peaks)

# -----------------------------
# 5. Genome-wide faceted plot for each haplotype
# -----------------------------
plot_fst_haplotype <- function(df, hap_label, chr9_color = "#D62828") {
  peak_df <- df %>%
    filter(CHROM == "chr_9A") %>%
    slice_max(WEIGHTED_FST, n = 1, with_ties = FALSE)
  
  ggplot(df, aes(x = mid_Mb, y = WEIGHTED_FST)) +
    geom_line(color = "grey45", linewidth = 0.5) +
    geom_line(
      data = df %>% filter(CHROM == "chr_9A"),
      color = chr9_color,
      linewidth = 0.8
    ) +
    geom_point(
      data = peak_df,
      aes(x = mid_Mb, y = WEIGHTED_FST),
      color = "black",
      size = 2.2,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ CHROM, scales = "free_x", ncol = 6) +
    labs(
      x = "Position (Mb)",
      y = "FST (Pop 1-2 vs Pop 3-5)",
      title = paste0("Genome-wide FST: ", hap_label)
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey85", colour = "grey40"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

p_hap1 <- plot_fst_haplotype(
  fst_all %>% filter(haplotype == "hap1"),
  "hap1"
)

p_hap2 <- plot_fst_haplotype(
  fst_all %>% filter(haplotype == "hap2"),
  "hap2",
  chr9_color = "#2A9D8F"
)

print(p_hap1)
print(p_hap2)

# save separate
ggsave("genomewide_FST_hap1.pdf", p_hap1, width = 16, height = 8)
ggsave("genomewide_FST_hap1.png", p_hap1, width = 16, height = 8, dpi = 300)

ggsave("genomewide_FST_hap2.pdf", p_hap2, width = 16, height = 8)
ggsave("genomewide_FST_hap2.png", p_hap2, width = 16, height = 8, dpi = 300)

# -----------------------------
# 6. Combined figure: hap1 vs hap2
# -----------------------------
combined_fst <- p_hap1 / p_hap2 +
  plot_annotation(
    title = "Genome-wide FST comparison across haplotypes"
  )

print(combined_fst)

ggsave("genomewide_FST_hap1_vs_hap2.pdf", combined_fst, width = 16, height = 14)
ggsave("genomewide_FST_hap1_vs_hap2.png", combined_fst, width = 16, height = 14, dpi = 300)

# -----------------------------
# 7. Chromosome-wide summary for each haplotype
# -----------------------------
fst_summary <- fst_all %>%
  group_by(haplotype, CHROM) %>%
  summarise(
    mean_fst = mean(WEIGHTED_FST, na.rm = TRUE),
    median_fst = median(WEIGHTED_FST, na.rm = TRUE),
    max_fst = max(WEIGHTED_FST, na.rm = TRUE),
    .groups = "drop"
  )

print(fst_summary)

write_tsv(fst_summary, "chromosome_FST_summary_hap1_hap2.tsv")

# -----------------------------
# 8. Barplot of chromosome-wide mean FST
# -----------------------------
fst_bar_data <- fst_summary %>%
  mutate(
    highlight = case_when(
      CHROM == "chr_9A" & haplotype == "hap1" ~ "chr9A_hap1",
      CHROM == "chr_9A" & haplotype == "hap2" ~ "chr9A_hap2",
      TRUE ~ "Other"
    )
  )

p_bar <- ggplot(fst_bar_data,
                aes(x = CHROM, y = mean_fst, fill = highlight)) +
  geom_col() +
  facet_wrap(~ haplotype, ncol = 1) +
  scale_fill_manual(values = c(
    "Other" = "grey70",
    "chr9A_hap1" = "#D62828",
    "chr9A_hap2" = "#2A9D8F"
  )) +
  labs(
    x = "Chromosome",
    y = "Mean FST",
    title = "Chromosome-wide mean FST across haplotypes"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_bar)

ggsave("chromosome_mean_FST_hap1_hap2.pdf", p_bar, width = 12, height = 8)
ggsave("chromosome_mean_FST_hap1_hap2.png", p_bar, width = 12, height = 8, dpi = 300)

# -----------------------------
# 9. chr_9A zoom for each haplotype
# -----------------------------
p_chr9_zoom <- ggplot(
  fst_all %>% filter(CHROM == "chr_9A"),
  aes(x = mid_Mb, y = WEIGHTED_FST, colour = haplotype)
) +
  geom_line(linewidth = 0.9) +
  geom_point(data = chr9_peaks, size = 2) +
  scale_colour_manual(values = c("hap1" = "#D62828", "hap2" = "#2A9D8F")) +
  labs(
    x = "Position on chr_9A (Mb)",
    y = "FST",
    colour = "Haplotype",
    title = "chr_9A FST peak comparison between haplotypes"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_chr9_zoom)
