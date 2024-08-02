library(tidyverse)
library(patchwork)
source("scripts/dosage_utils.R")

# processed PAU dfs produced by calculate_qapa_means_psis.R
pau_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.mean_pau.dosage_labels.tsv")
pau_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.mean_pau.dosage_labels.tsv")

# processed DEXSeq results tables filtered to representative PAS produced by mef_representative_pas.R
dexseq_rep_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")
dexseq_rep_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")


outdir <- "processed/2024-02-29/analysis/mef_comparisons"

# Join in DEXSeq with PAU results
dexseq_rep_m323k <- pau_m323k %>%
  # already have PAUs as column names
  select(-condition) %>%
  distinct(APA_ID, .keep_all = T) %>%
  left_join(dexseq_rep_m323k, ., by = c("APA_ID", "Gene_Name"))

dexseq_rep_f210i <- pau_f210i %>%
  # already have PAUs as column names
  select(-condition) %>%
  distinct(APA_ID, .keep_all = T) %>%
  left_join(dexseq_rep_f210i, ., by = c("APA_ID", "Gene_Name"))


# Select genes with any sig PAS, then select proximal/distal as representative
dexseq_rep_m323k_padj <- dexseq_rep_m323k %>%
  group_by(Gene, Chr) %>%
  filter(any(padj < 0.05)) %>%
  ungroup()

# now select proximal site per gene
dexseq_rep_m323k_prox_padj <- dexseq_rep_m323k_padj %>%
  group_by(Gene, Chr) %>%
  filter(Pas_Number == min(Pas_Number)) %>%
  ungroup()

# alternatively, do for distal PAS
dexseq_rep_m323k_dist_padj <- dexseq_rep_m323k_padj %>%
  group_by(Gene, Chr) %>%
  filter(Pas_Number != min(Pas_Number)) %>%
  ungroup()

# repeat for F210I
dexseq_rep_f210i_padj <- dexseq_rep_f210i %>%
  group_by(Gene, Chr) %>%
  filter(any(padj < 0.05)) %>%
  ungroup()

# now select proximal site per gene
dexseq_rep_f210i_prox_padj <- dexseq_rep_f210i_padj %>%
  group_by(Gene, Chr) %>%
  filter(Pas_Number == min(Pas_Number)) %>%
  ungroup()

# alternatively, do for distal PAS
dexseq_rep_f210i_dist_padj <- dexseq_rep_f210i_padj %>%
  group_by(Gene, Chr) %>%
  filter(Pas_Number != min(Pas_Number)) %>%
  ungroup()

# construct 'strict' dfs, where all PAS have padj < 0.05
dexseq_rep_m323k_prox_padj_strict <- filter(dexseq_rep_m323k_prox_padj, padj < 0.05)
dexseq_rep_m323k_dist_padj_strict <- filter(dexseq_rep_m323k_dist_padj, padj < 0.05)

dexseq_rep_f210i_prox_padj_strict <- filter(dexseq_rep_f210i_prox_padj, padj < 0.05)
dexseq_rep_f210i_dist_padj_strict <- filter(dexseq_rep_f210i_dist_padj, padj < 0.05)

# construct pvalue histogram of representative PAS for two datasets (TODO: repeat for all PAS)
dexseq_rep_comb <- bind_rows(F210I = dexseq_rep_f210i, M323K = dexseq_rep_m323k, .id = "mutant")

comb_pval_hist <- dexseq_rep_comb %>%
  ggplot(aes(x = pvalue)) +
  facet_wrap("~ mutant", ncol = 1, scales = "free_y") +
  geom_histogram(binwidth = 0.005) +
  theme_bw(base_size = 14)

# Sig gene count with increasing delta PAU cutoffs

sig_gene_delta_count <- c(0,5,10,15,20) %>%
  set_names() %>%
  map(~ dexseq_rep_comb %>%
        filter(padj < 0.05) %>%
        filter(abs(delta_HOM_WT) >= .x) %>%
        distinct(mutant, Gene, Chr)) %>%
  bind_rows(.id = "abs_delta_cutoff") %>%
  mutate(abs_delta_cutoff = as.integer(abs_delta_cutoff)) %>%
  count(mutant, abs_delta_cutoff)


sig_gene_delta_count_plot <- sig_gene_delta_count %>%
  ggplot(aes(x = abs_delta_cutoff, y = n, label = n, fill = mutant, group = mutant)) +
  geom_col(position = "dodge", width = 3) +
  geom_text(position = position_dodge(width = 3), vjust = -0.5) +
  scale_y_continuous(breaks = seq(0,500,50)) +
  scale_fill_manual(values = c("M323K" = "#ef8a62", "F210I" = "#67a9cf", "Both" = "black", "Q331K" = "#7570b3")) +
  theme_bw(base_size = 14) +
  labs(x = "Minimum absolute change in usage (HOM - WT, %)",
       y = "Gene Count",
       fill = "Mutant") +
  theme(legend.position = "bottom")

sig_gene_delta_count_plot


# Repeat counts, but split by proximal/distal

sig_pas_type_count <- dexseq_rep_comb %>%
  group_by(mutant, Gene, Chr) %>%
  mutate(site_type = if_else(Pas_Number == min(Pas_Number), "Proximal", "Distal"),
         site_type = factor(site_type, levels = c("Proximal", "Distal"))) %>%
  ungroup() %>%
  filter(padj < 0.05) %>%
  distinct(mutant, Gene, Chr, site_type) %>%
  count(mutant, site_type)

sig_pas_type_count %>%
  ggplot(aes(x = site_type, y = n, label = n, fill = mutant, group = mutant)) +
  geom_col(position = "dodge") +
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_manual(values = c("M323K" = "#ef8a62", "F210I" = "#67a9cf", "Both" = "black", "Q331K" = "#7570b3")) +
  theme_bw(base_size = 14) +
  labs(x = "Target PAS Type",
       y = "Count",
       fill = "Mutant") +
  theme(legend.position = "bottom")

# Repeat counting as increase delta threshold
sig_pas_type_delta_count <- c(0,5,10,15,20) %>%
  set_names() %>%
  map(~ dexseq_rep_comb %>%
        group_by(mutant, Gene, Chr) %>%
        mutate(site_type = if_else(Pas_Number == min(Pas_Number), "Proximal", "Distal"),
               site_type = factor(site_type, levels = c("Proximal", "Distal"))) %>%
        ungroup() %>%
        filter(padj < 0.05,
               abs(delta_HOM_WT) >= .x) %>%
        distinct(mutant, Gene, Chr, site_type)) %>%
  bind_rows(.id = "abs_delta_cutoff") %>%
  mutate(abs_delta_cutoff = as.integer(abs_delta_cutoff)) %>%
  count(mutant, abs_delta_cutoff, site_type)

sig_pas_type_delta_count

#
# Dosage plots
#

# First, make standard dosage plots with no delta threshold
# here, list names will be passed as subtitles as MUT the same across datasets
plot_df_list <- list("M323K - Proximal PAS, Gene any padj < 0.05" = dexseq_rep_m323k_prox_padj,
                     "F210I - Proximal PAS, Gene any padj < 0.05" = dexseq_rep_f210i_prox_padj,
                     "M323K - Proximal PAS, padj < 0.05" = dexseq_rep_m323k_prox_padj_strict,
                     "F210I - Proximal PAS, padj < 0.05" = dexseq_rep_f210i_prox_padj_strict,
                     "M323K - Distal PAS, Gene any padj < 0.05" = dexseq_rep_m323k_dist_padj,
                     "F210I - Distal PAS, Gene any padj < 0.05" = dexseq_rep_f210i_dist_padj,
                     "M323K - Distal PAS, padj < 0.05" = dexseq_rep_m323k_dist_padj_strict,
                     "F210I - Distal PAS, padj < 0.05" = dexseq_rep_f210i_dist_padj_strict)

# Make plots for all dosage categories
mef_dosage_all <- map2(.x = plot_df_list,
                       .y = names(plot_df_list),
                       ~ .x %>%
                         dosage_analysis_wrapper(line_labs = labs(title = "MEF",
                                                                  subtitle = .y
                         ),
                         bar_labs = labs(title = "MEF",
                                         subtitle = .y
                         )
                         )
)


# Make plots subsetted to 'clear'/bonafide consistent/inconsistent dosage categories
mef_dosage_subset <- map2(.x = plot_df_list,
                          .y = names(plot_df_list),
                          ~ dosage_analysis_wrapper(df = .x,
                                                    patterns_filter = meanPAU.dosage_pattern %in% c("Up+Up", "Down+Down", "Down+Up", "Up+Down"),
                                                    line_labs = labs(title = "MEF",
                                                                     subtitle = .y),
                                                    bar_labs = labs(title = "MEF",
                                                                    subtitle = .y)
                          )
)


# make plots subsetted to same-containing categories
mef_dosage_subset_same <- map2(.x = plot_df_list,
                               .y = names(plot_df_list),
                               ~ dosage_analysis_wrapper(df = .x,
                                                         patterns_filter = meanPAU.dosage_pattern %in% c("Down+Same", "Up+Same", "Same+Up", "Same+Down"),
                                                         line_labs = labs(title = "MEF",
                                                                          subtitle = .y),
                                                         bar_labs = labs(title = "MEF",
                                                                         subtitle = .y)
                               )
)

# print out all bar plots containing event counts

map(mef_dosage_all, ~ .x$bar_plot)
map(mef_dosage_subset, ~ .x$bar_plot)
map(mef_dosage_subset_same, ~ .x$bar_plot)

map(mef_dosage_all, ~ .x$dosage_plot)
map(mef_dosage_subset, ~ .x$dosage_plot)
map(mef_dosage_subset_same, ~ .x$dosage_plot)


# ATTEMPTS TO MAKE COMBINED PANELS

# combined list of individual plots, removing the title + subtitle from count plots (will be aligned/duplicated vertically)
tmp_dosage_subset <- map(mef_dosage_subset, ~.x$dosage_plot) %>%
  set_names(paste(names(.),
                  "dosage_plot", sep=" - ")
  ) %>%
  map(~ .x + theme_bw(base_size = 12) + labs(title = NULL))

tmp_bar_subset <- map(mef_dosage_subset,
    ~ .x$bar_plot + labs(title = NULL, subtitle = NULL) + theme_bw(base_size = 12)
) %>%
  set_names(paste(names(.), "bar_plot", sep = " - "))

comb_dosage_bar_subset <- c(tmp_dosage_subset, tmp_bar_subset)
names(comb_dosage_bar_subset)

# Generate every combo of 4 plots based on 'target set'
target_groups <- unique(str_split_i(names(comb_dosage_bar_subset), " - ", 2))

# convert to list of n=4 character vectors (F210I | M323K horizontally, with dosage plot at top vertically)
target_groups_plot_order <- pmap(list(f210i_dosage = paste0("F210I - ", target_groups, " - dosage_plot"),
          m323k_dosage = paste0("M323K - ", target_groups, " - dosage_plot"), 
          f210i_bar = paste0("F210I - ", target_groups, " - bar_plot"), 
          m323k_bar = paste0("M323K - ", target_groups, " - bar_plot")
          ),
     function(f210i_dosage, m323k_dosage, f210i_bar, m323k_bar) c(f210i_dosage, m323k_dosage, f210i_bar, m323k_bar)
     ) %>%
  set_names(target_groups) 

# Generate 4 panel combined dosage + count figures for each target set
dosage_count_comb_subset_pau <- target_groups_plot_order %>%
  map(~ wrap_plots(comb_dosage_bar_subset[.x],
                   ncol = 2, 
                   nrow = 2,
                   byrow = T,
                   axes = "collect",
                   axis_titles = "collect") +
        plot_annotation(tag_levels = "A"))

dosage_count_comb_subset_pau

# Repeat for same categories
# combined list of individual plots, removing the title + subtitle from count plots (will be aligned/duplicated vertically)
tmp_dosage_subset_same <- map(mef_dosage_subset_same, ~.x$dosage_plot) %>%
  set_names(paste(names(.),
                  "dosage_plot", sep=" - ")
  ) %>%
  map(~ .x + theme_bw(base_size = 12) + labs(title = NULL))

tmp_bar_subset_same <- map(mef_dosage_subset_same,
                      ~ .x$bar_plot + labs(title = NULL, subtitle = NULL) + theme_bw(base_size = 12)
) %>%
  set_names(paste(names(.), "bar_plot", sep = " - "))

comb_dosage_bar_subset_same <- c(tmp_dosage_subset_same, tmp_bar_subset_same)
names(comb_dosage_bar_subset_same)

# Generate 4 panel combined dosage + count figures for each target set
dosage_count_comb_subset_same_pau <- target_groups_plot_order %>%
  map(~ wrap_plots(comb_dosage_bar_subset_same[.x],
                   ncol = 2, 
                   nrow = 2,
                   byrow = T,
                   axes = "collect",
                   axis_titles = "collect") +
        plot_annotation(tag_levels = "A"))

dosage_count_comb_subset_same_pau




# Plot differences between genotypes for inconsistent categories - to waht extent reflect small changes in usage (i.e. noise in means?)

# Plot for both proximal and distal PAS
mef_subset_inconsistent_deltas <- map(mef_dosage_subset[str_ends(names(mef_dosage_subset), "PAS, padj < 0.05")], 
    ~ .x$plot_df %>%
      filter(meanPAU.dosage_pattern_simplified == "not_consistent") %>%
      distinct(Gene, Chr, meanPAU.dosage_pattern, delta_HET_WT, delta_HOM_HET)
) %>%
  bind_rows(.id = "origin") %>%
  pivot_longer(cols = starts_with("delta"), names_to = "comparison", values_to = "delta", names_prefix = "delta_") %>%
  mutate(abs_delta = abs(delta),
         mutant = str_split_i(origin, " - ", 1),
         pas = str_split_i(str_split_i(origin, " PAS,", 1), " - ", 2),
         pas = factor(pas, levels = c("Proximal", "Distal")),
         comparison = str_replace_all(comparison, "_", " - ")
         ) %>%
  ggplot(aes(x = meanPAU.dosage_pattern, y = abs_delta, colour = mutant)) +
  facet_wrap("pas ~ comparison") +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_continuous(breaks = seq(0,80,20)) +
  scale_colour_manual(values = c("M323K" = "#ef8a62", "F210I" = "#67a9cf")) +
  theme_bw(base_size = 14) +
  labs(x = "PAS Usage Pattern (WT - HET - HOM)",
       y = "Absolute Change in Usage (%)",
       colour = "Mutant") +
  theme(legend.position = "top")

mef_subset_inconsistent_deltas


#
# attempt using both proximal and distal PAS, but this time blot the minimum diff between the two comparisons
map(mef_dosage_subset[str_ends(names(mef_dosage_subset), "PAS, padj < 0.05")], 
    ~ .x$plot_df %>%
      filter(meanPAU.dosage_pattern_simplified == "not_consistent") %>%
      distinct(Gene, Chr, meanPAU.dosage_pattern, delta_HET_WT, delta_HOM_HET)
) %>%
  bind_rows(.id = "origin") %>%
  pivot_longer(cols = starts_with("delta"), names_to = "comparison", values_to = "delta", names_prefix = "delta_") %>%
  mutate(abs_delta = abs(delta),
         mutant = str_split_i(origin, " - ", 1),
         pas = str_split_i(str_split_i(origin, " PAS,", 1), " - ", 2),
         pas = factor(pas, levels = c("Proximal", "Distal")),
         comparison = str_replace_all(comparison, "_", " - ")
  ) %>%
  group_by(mutant, pas, Gene, Chr) %>%
  slice_min(abs_delta, n = 1) %>%
  ungroup() %>%
  ggplot(aes(x = meanPAU.dosage_pattern, y = abs_delta, colour = mutant)) +
  facet_wrap("~ pas") +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_continuous(breaks = seq(0,40,5)) +
  scale_colour_manual(values = c("M323K" = "#ef8a62", "F210I" = "#67a9cf")) +
  theme_bw(base_size = 14) +
  labs(x = "PAS Usage Pattern (WT - HET - HOM)",
       y = "Absolute Change in Usage (%)",
       colour = "Mutant") +
  theme(legend.position = "bottom")


## SAVE PLOTS TO DISK

# output to subdirectories for each target + significance type - extract individual tyoes
target_groups_spl <- str_split_fixed(target_groups, ", ", 2)
pas_type <- str_replace_all(target_groups_spl[, 1], pattern = " ", "_")
sig_type <- str_replace_all(target_groups_spl[, 2], pattern = "  ", "") %>% str_replace_all(" ", "_") %>% str_replace_all("\\<", "lt")
pas_type
sig_type

# create output subdirectories for each combo
outdir_subdirs <- file.path(outdir, pas_type, sig_type)
outdir_subdirs
walk(outdir_subdirs, ~ if (!dir.exists(.x)) {dir.create(.x, recursive = T)})

# assign target_groups as names so can match plot to outdir
names(outdir_subdirs) <- target_groups

# sanity check - order should be identical
stopifnot(all(names(outdir_subdirs) == names(dosage_count_comb_subset_pau)))

# save combined plots to disk
walk2(.x = dosage_count_comb_subset_pau,
      .y = outdir_subdirs,
      ~ ggsave(filename = file.path(.y,
                                    "2024-06-25_mef_dosage_counts.pau.combined.png"),
               plot = .x,
               device = "png",
               units = "mm",
               dpi = "retina",
               height = 200,
               width = 200
               )
      )

walk2(.x = dosage_count_comb_subset_pau,
      .y = outdir_subdirs,
      ~ ggsave(filename = file.path(.y,
                                    "2024-06-25_mef_dosage_counts.pau.combined.pdf"),
               plot = .x,
               device = "pdf",
               units = "mm",
               dpi = "retina",
               height = 200,
               width = 200
      )
)


# repeat for same-containing categories

walk2(.x = dosage_count_comb_subset_same_pau,
      .y = outdir_subdirs,
      ~ ggsave(filename = file.path(.y,
                                    "2024-06-25_mef_dosage_counts.pau.combined.same_containing.png"),
               plot = .x,
               device = "png",
               units = "mm",
               dpi = "retina",
               height = 200,
               width = 200
      )
)

walk2(.x = dosage_count_comb_subset_same_pau,
      .y = outdir_subdirs,
      ~ ggsave(filename = file.path(.y,
                                    "2024-06-25_mef_dosage_counts.pau.combined.same_containing.pdf"),
               plot = .x,
               device = "pdf",
               units = "mm",
               dpi = "retina",
               height = 200,
               width = 200
      )
)


# repeat for inconsistent deltas
ggsave(filename =  file.path(outdir, "2024-06-29_mef_dosage.inconsistent_abs_deltas.proximal_distal_facet.png"),
       plot = mef_subset_inconsistent_deltas,
       device = "png",
       units = "mm",
       dpi = "retina",
       height = 150,
       width = 150
       )

ggsave(filename =  file.path(outdir, "2024-06-29_mef_dosage.inconsistent_abs_deltas.proximal_distal_facet.pdf"),
       plot = mef_subset_inconsistent_deltas,
       device = "pdf",
       units = "mm",
       dpi = "retina",
       height = 150,
       width = 150
)

save.image(file.path(outdir, "2024-06-29_mef_dosage.RData"))
