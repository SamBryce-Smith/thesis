library(tidyverse)
source("scripts/comparison_utils.R")

# processed PAU dfs produced by calculate_qapa_means_psis.R
pau_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.mean_pau.dosage_labels.tsv")
pau_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.mean_pau.dosage_labels.tsv")

# psi tables containing condition means and deltas
psi_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.psi.mean_median.dosage_labels.tsv")
psi_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.psi.mean_median.dosage_labels.tsv")

# processed DEXSeq results tables filtered to representative PAS produced by mef_representative_pas.R
dexseq_rep_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")
dexseq_rep_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")

# processed DEXSeq results tables produced by QAPA_snakemake pipeline
dexseq_f210i <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")
dexseq_m323k <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")

outdir <- "processed/2024-02-29/analysis/mef_comparisons/comparison"
if (!dir.exists(outdir)) { dir.create(outdir, recursive = T)}


# add in PAU information
dexseq_m323k <- left_join(dexseq_m323k, distinct(pau_m323k, APA_ID, Gene_Name, .keep_all = T) %>% select(-condition), 
                              by = c("APA_ID", "Gene_Name")) %>%
  mutate(gene_chr = paste(Gene, Chr, sep = "_"))
dexseq_f210i <- left_join(dexseq_f210i, distinct(pau_f210i, APA_ID, Gene_Name, .keep_all = T) %>% select(-condition), 
                              by = c("APA_ID", "Gene_Name")) %>%
  mutate(gene_chr = paste(Gene, Chr, sep = "_"))

dexseq_rep_m323k <- left_join(dexseq_rep_m323k, distinct(pau_m323k, APA_ID, Gene_Name, .keep_all = T) %>% select(-condition), 
                              by = c("APA_ID", "Gene_Name")) %>%
  mutate(gene_chr = paste(Gene, Chr, sep = "_"))
dexseq_rep_f210i <- left_join(dexseq_rep_f210i, distinct(pau_f210i, APA_ID, Gene_Name, .keep_all = T) %>% select(-condition), 
                              by = c("APA_ID", "Gene_Name")) %>%
  mutate(gene_chr = paste(Gene, Chr, sep = "_"))

# Add in psi information
dexseq_f210i <- left_join(dexseq_f210i, psi_f210i, by = c("Gene", "Chr"))
dexseq_m323k <- left_join(dexseq_m323k, psi_m323k, by = c("Gene", "Chr"))
dexseq_rep_f210i <- left_join(dexseq_rep_f210i, psi_f210i, by = c("Gene", "Chr"))
dexseq_rep_m323k <- left_join(dexseq_rep_m323k, psi_m323k, by = c("Gene", "Chr"))

dexseq_rep_m323k_targets <- define_apa_targets(dexseq_rep_m323k)
dexseq_rep_f210i_targets <- define_apa_targets(dexseq_rep_f210i)



# have df defined for each 'subset' of targets (proximal, distal, either prox/distal padj < 0.05 + gene (return all evens))
# produce bar plot for all combinations of f210i + m323k for these subsets
psi_bar_plots <- pmap(.l = list(f210i = dexseq_rep_f210i_targets,
               m323k = dexseq_rep_m323k_targets,
               nm = names(dexseq_rep_f210i_targets)),
     .f = function(f210i, m323k, nm) psi_direction_bar_plot(list("F210I" = f210i,
                                                                 "M323K" = m323k),
                                                            plot_title = nm,
                                                            padj_method = "bonferroni")
     ) 

# pull out plots for each subset
psi_bar_plots_p <- map(psi_bar_plots, "plot")

# binomial test results for each subset. Note input counts are always for shortening category
map(psi_bar_plots, "binom_test") %>%
  bind_rows(.id = "target_subset")

psi_bar_plots_p


# Two way comparison of delta usage values

# Compute intersection counts + values for pairwise comparison of target list types (e.g. site-level/gene-level)
m323k_f210i_overlap <- compute_overlap_statistics(dexseq_rep_m323k_targets, dexseq_rep_f210i_targets)
m323k_f210i_overlap

# pull out vectors of intersecting IDs
m323k_f210i_overlap$intersecting_elements %>%
  set_names(paste(m323k_f210i_overlap$comparison_type, m323k_f210i_overlap$comparison_value, sep = "__"))


# test run - APA_IDs for distal
dist_overlap_ids <- m323k_f210i_overlap %>%
  filter(comparison_type == "distal" & comparison_value == "APA_ID") %>%
  pull(intersecting_elements) %>%
  unlist()

gene_overlap_ids <- m323k_f210i_overlap %>%
  filter(comparison_type == "gene" & comparison_value == "Gene_Chr") %>%
  pull(intersecting_elements) %>%
  unlist()

all_overlap_ids <- m323k_f210i_overlap %>%
  filter(comparison_type == "all" & comparison_value == "Gene_Chr") %>%
  pull(intersecting_elements) %>%
  unlist()

# assign labels to sig targets based on overlap

# Compute Euler plots of overlapping target genes

# Need to take initial dfs, label with target (T/F) based on list from *_rep_targets, then pass to euler_plot_wrapper
dexseq_rep_m323k <- left_join(dexseq_rep_m323k,
          distinct(dexseq_rep_m323k_targets$all, Gene, Chr) %>% mutate(target = T),
          by = c("Gene", "Chr")
          ) %>%
  mutate(target = replace_na(target, F),
         # gene_chr = paste(Gene, Chr, sep = "_")
         )

dexseq_rep_f210i <- left_join(dexseq_rep_f210i,
          distinct(dexseq_rep_f210i_targets$all, Gene, Chr) %>% mutate(target = T),
          by = c("Gene", "Chr")
          ) %>%
  mutate(target = replace_na(target, F),
         # gene_chr = paste(Gene, Chr, sep = "_")
         )

f210i_m323k_euler <- euler_plot_wrapper(dexseq_rep_f210i, dexseq_rep_m323k, "F210I", "M323K")

f210i_m323k_euler$all_plot
f210i_m323k_euler$target_plot


# Generate delta psi scatter plot of common & mutant-specific targets  (+ chi-squared)

m323k_q331k_comparison_all <- twoway_comparison_scatter(df1 = dexseq_rep_m323k_targets$all,
                          df2 = dexseq_rep_f210i_targets$all,
                          orig_df1 = dexseq_m323k,
                          orig_df2 = dexseq_f210i,
                          df1_lab = "M323K",
                          df2_lab = "F210I",
                          common_ids = all_overlap_ids,
                          id_cols = c("Gene", "Chr", "Gene_Name")
                          )

# How many dropped (missing psi from other MUT)
m323k_q331k_comparison_all$plot_df %>%
  mutate(missing = is.na(deltaPSImean.HOM_WT.M323K) | is.na(deltaPSImean.HOM_WT.F210I)) %>%
  count(overlap_group, missing)

m323k_q331k_comparison_all


# tidy up axis labels
m323k_q331k_comparison_all$scatter_plot <- m323k_q331k_comparison_all$scatter_plot +
  labs(x = "Change in psi (HOM - WT) - M323K",
       y = "Change in psi (HOM - WT) - F210I") +
  theme(legend.position = "top")

m323k_q331k_comparison_all$scatter_plot

# Save plots to disk

# Psi direction bar plots

walk2(.x = psi_bar_plots_p,
      .y = names(psi_bar_plots_p),
      ~ ggsave(filename = file.path(outdir,
                                             paste("2024-06-29_mef_comparison",
                                                   "psi",
                                                   "direction_bar",
                                                   .y,
                                                   "png",
                                                   sep = ".")
                                             ),
      plot = .x + labs(title = NULL, fill = "3'UTR Change") + theme(legend.position = "top"),
      device = "png",
      units = "mm",
      dpi = "retina",
      height = 125,
      width = 125
      )
      )

walk2(.x = psi_bar_plots_p,
      .y = names(psi_bar_plots_p),
      ~ ggsave(filename = file.path(outdir,
                                    paste("2024-06-29_mef_comparison",
                                          "psi",
                                          "direction_bar",
                                          .y,
                                          "pdf",
                                          sep = ".")
      ),
      plot = .x + labs(title = NULL, fill = "3'UTR Change") + theme(legend.position = "top"),
      device = "pdf",
      units = "mm",
      dpi = "retina",
      height = 125,
      width = 125
      )
)
      
# Euler plots

ggsave(filename = file.path(outdir,
                            paste("2024-06-29_mef_comparison",
                                  "euler",
                                  "m323k_f210i_all_genes",
                                  "png",
                                  sep = ".")
),
plot = f210i_m323k_euler$all_plot,
device = "png",
dpi = "retina",
units = "mm",
height = 125,
width = 125
)


ggsave(filename = file.path(outdir,
                            paste("2024-06-29_mef_comparison",
                                  "euler",
                                  "m323k_f210i_all_genes",
                                  "pdf",
                                  sep = ".")
),
plot = f210i_m323k_euler$all_plot,
device = "pdf",
dpi = "retina",
units = "mm",
height = 125,
width = 125
)

# psi targets scatter
ggsave(filename = file.path(outdir,
                            paste("2024-06-29_mef_comparison",
                                  "psi",
                                  "m323k_f210i_target_scatter",
                                  "png",
                                  sep = ".")
                            ),
       plot = m323k_q331k_comparison_all$scatter_plot,
       device = "png",
       dpi = "retina",
       units = "mm",
       height = 125,
       width = 125
       )

ggsave(filename = file.path(outdir,
                            paste("2024-06-29_mef_comparison",
                                  "psi",
                                  "m323k_f210i_target_scatter",
                                  "pdf",
                                  sep = ".")
),
plot = m323k_q331k_comparison_all$scatter_plot,
device = "pdf",
dpi = "retina",
units = "mm",
height = 125,
width = 125
)

save.image(file.path(outdir, "2024-06-29_mef_comparison.RData"))

