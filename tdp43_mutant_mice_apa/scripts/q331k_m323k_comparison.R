library(tidyverse)
source("scripts/comparison_utils.R")

# processed PAU dfs produced by calculate_qapa_means_psis.R
#   pau_f210i <- read_tsv("processed/2024-02-29/Q331K_3month_Adult_Brain/.mean_pau.dosage_labels.tsv")
# pau_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.mean_pau.dosage_labels.tsv")

# psi tables containing condition means and deltas
# psi_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.psi.mean_median.dosage_labels.tsv")
# psi_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.psi.mean_median.dosage_labels.tsv")

# processed DEXSeq results tables filtered to representative PAS produced by mef_representative_pas.R
# dexseq_rep_q331k <- read_tsv("processed/2024-02-29/Q331K_3month_Adult_Brain/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")
# dexseq_rep_m323k <- read_tsv("processed/2024-02-29/M323K_3month_Adult_Brain/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")

pau_q331k <- read_tsv("processed/2024-02-29/Q331K_3month_Adult_Brain/Q331K_3month_Adult_Brain.mean_pau.tsv")
pau_m323k <- read_tsv("processed/2024-02-29/M323K_3month_Adult_Brain/M323K_3month_Adult_Brain.mean_pau.tsv")

psi_q331k <- read_tsv("processed/2024-02-29/Q331K_3month_Adult_Brain/Q331K_3month_Adult_Brain.psi.mean_median.tsv")
psi_m323k <- read_tsv("processed/2024-02-29/M323K_3month_Adult_Brain/M323K_3month_Adult_Brain.psi.mean_median.tsv")


# processed DEXSeq results tables produced by QAPA_snakemake pipeline
dexseq_q331k <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/Q331K_3month_Adult_Brain/dexseq_apa.HOMvsWT.results.processed.tsv")
dexseq_m323k <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/M323K_3month_Adult_Brain/dexseq_apa.HOMvsWT.results.processed.tsv")

outdir <- "processed/2024-02-29/analysis/spinal_cord_m323k_q331k_comparisons"
if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

# p-value histograms
pval_hist_comb <- bind_rows("M323K" = dexseq_m323k, "Q331K" = dexseq_q331k, .id = "mutant") %>%
  ggplot(aes(x = pvalue)) + 
  facet_wrap("~ mutant", nrow = 2, scales = "free_y") + 
  geom_histogram(bins = 1000) +
  theme_bw(base_size = 14) +
  labs(x = "P-value",
       y = "Count")

pval_hist_comb

# both have peak at low p-values, although stronger for m323k. Both also have weird peak at 1 (?)

# add in PAU information
dexseq_q331k <- left_join(dexseq_q331k, distinct(pau_q331k, APA_ID, Gene_Name, .keep_all = T) %>% select(-condition), 
                          by = c("APA_ID", "Gene_Name"))
dexseq_m323k <- left_join(dexseq_m323k, distinct(pau_m323k, APA_ID, Gene_Name, .keep_all = T) %>% select(-condition), 
                          by = c("APA_ID", "Gene_Name"))

# add in psi information
dexseq_q331k <- left_join(dexseq_q331k, psi_q331k, 
                          by = c("Gene", "Chr"))
dexseq_m323k <- left_join(dexseq_m323k, psi_m323k, 
                          by = c("Gene", "Chr"))



#
# HOW DOES NUMBER OF GENES VARY WITH A MINIMUM EFFECT SIZE THRESHOLD?
#

dexseq_comb <- bind_rows(Q331K = dexseq_q331k, M323K = dexseq_m323k, .id = "mutant")

sig_gene_delta_count <- c(0,5,10,15,20) %>%
  set_names() %>%
  map(~ dexseq_comb %>%
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
  scale_fill_manual(values = c("#1b9e77", "#7570b3")) +
  theme_bw(base_size = 14) +
  labs(x = "Minimum absolute change in usage (HOM - WT, %)",
       y = "Gene Count",
       fill = "Mutant") +
  theme(legend.position = "bottom")

sig_gene_delta_count_plot


# Produce PSI bar plot for all combinations of q331k + m323k and delta cutoffs

# Generate dfs of target PAS at increasing cutoffs

pau_cutoffs <- c(0,5,10,15,20) %>%
  set_names()

dexseq_q331k_targets_pau <- pau_cutoffs %>%
  map(~ dexseq_q331k %>%
        filter(padj < 0.05,
               abs(delta_HOM_WT) >= .x)
      )

dexseq_m323k_targets_pau <- pau_cutoffs %>%
  map(~ dexseq_m323k %>%
        filter(padj < 0.05,
               abs(delta_HOM_WT) >= .x)
  )

# produce bar plot for all delta cut-offs between each of the two MUTs
psi_bar_plots <- pmap(.l = list(q331k = dexseq_q331k_targets_pau,
                                m323k = dexseq_m323k_targets_pau,
                                nm = names(pau_cutoffs)),
                      .f = function(q331k, m323k, nm) psi_direction_bar_plot(list("Q331K" = q331k,
                                                                                  "M323K" = m323k),
                                                                             plot_title = paste("Absolute (HOM - WT) >=", nm, "%", sep = " "),
                                                                             padj_method = "BH")
                      )  

# pull out plots for each subset
psi_bar_plots_p <- map(psi_bar_plots, "plot")

# binomial test results for each subset. Note input counts are always for shortening category
map(psi_bar_plots, "binom_test") %>%
  bind_rows(.id = "abs_pau_cutoff")

psi_bar_plots_p


# compute bar plot of psi dosage for all PAU cutoffs

# extract binomial test df, readjust pvalues for comparisons at each cutoff
psi_comb_binom <- map(psi_bar_plots, "binom_test") %>%
  bind_rows(.id = "abs_pau_cutoff") %>%
  # convert to factor, ordered from low-high cutoff for easier plotting
  mutate(abs_pau_cutoff = fct_inseq(abs_pau_cutoff)) %>%
  adjust_pvalue(output.col = "p.adj.all", method = "BH") %>%
  add_significance(p.col = "p.adj.all", output.col = "p.adj.all.signif") %>%
  # standardise signif label y positions across MUTs (makes sense if horizontally facet, otherwise group)
  mutate(y.position = max(estimate * 100, (1 - estimate) * 100) + 10) %>%
  ungroup()

psi_comb_binom

# extract combined direction counts pplot, again ordering cutoff low-high
psi_comb_counts <- map(psi_bar_plots, "counts_df") %>%
  bind_rows(.id = "abs_pau_cutoff") %>%
  mutate(abs_pau_cutoff = fct_inseq(abs_pau_cutoff))

psi_comb_delta_thresh_bar <- ggplot(psi_comb_counts, aes(x = abs_pau_cutoff, y = perc, label = n)) +
  facet_wrap("~ mutant", ncol = 2) +
  geom_col(aes(fill = dirn_mean), position = "dodge") +
  geom_text(aes(group = dirn_mean), position = position_dodge(width = 0.9), vjust = -1) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
  ggprism::add_pvalue(psi_comb_binom, label = "p.adj.all.signif", x = "abs_pau_cutoff", label.size = 6) +
  scale_fill_manual(values = c("#b2df8a", "#1f78b4")) +
  labs(title = NULL,
       x = "Minimum Absolute (HOM - WT) Threshold (%)",
       y = "Target APAs (%)",
       fill = "3'UTR length change") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

psi_comb_delta_thresh_bar

#
# Compute overlaps between all possible combinations of significance + effect size cutoffs
#
pau_cutoff_combos <- expand_grid(pau_cutoffs, pau_cutoffs)
colnames(pau_cutoff_combos) <- c("pau_cutoff_q331k", "pau_cutoff_m323k")
pau_cutoff_combos

# generate pretty names of combined thresholds
cutoff_combo_nms <- paste("Q331K_", pau_cutoff_combos$pau_cutoff_q331k,"__M323K_", pau_cutoff_combos$pau_cutoff_m323k, sep = "")


# Generate overlap statistics for all possible combinations of targets defined by just PAU cutoff
q331k_m323k_cutoffs_overlap <- map2(.x = pau_cutoff_combos$pau_cutoff_q331k,
     .y = pau_cutoff_combos$pau_cutoff_m323k,
     ~ compute_overlap_statistics(list(dexseq_q331k_targets_pau[[as.character(.x)]]),
                                  list(dexseq_m323k_targets_pau[[as.character(.y)]])
                                  )
     ) %>%
  set_names(cutoff_combo_nms) %>%
  bind_rows(.id = "cutoff_combo")

q331k_m323k_cutoffs_overlap


# Compute Euler plots of overlapping target genes

# Need to take initial dfs, label with target (T/F) based on list from *_rep_targets, then pass to euler_plot_wrapper
dexseq_m323k <- left_join(dexseq_m323k,
                              distinct(dexseq_m323k_targets_pau$`0`, Gene, Chr) %>% mutate(target = T),
                              by = c("Gene", "Chr")
          ) %>%
  mutate(target = replace_na(target, F),
         gene_chr = paste(Gene, Chr, sep = "_")
  )

dexseq_q331k <- left_join(dexseq_q331k,
                              distinct(dexseq_q331k_targets_pau$`0`, Gene, Chr) %>% mutate(target = T),
                              by = c("Gene", "Chr")
) %>%
  mutate(target = replace_na(target, F),
         gene_chr = paste(Gene, Chr, sep = "_")
  )

q331k_m323k_euler <- euler_plot_wrapper(dexseq_q331k, dexseq_m323k, "Q331K", "M323K")

q331k_m323k_euler$all_plot
q331k_m323k_euler$target_plot



# Create two-way scatter plot comparison of gene-level psi values

# Need: vector of common ids, gene_chr column, psi values, subsetted + original dfs

# add gene_chr to target dfs (so aligns with output of compute_overlap_statistics)
dexseq_m323k_targets_pau <- map(dexseq_m323k_targets_pau, ~ mutate(.x, gene_chr = paste(Gene, Chr, sep = "_")))
dexseq_q331k_targets_pau <- map(dexseq_q331k_targets_pau, ~ mutate(.x, gene_chr = paste(Gene, Chr, sep = "_")))

# get vector of common gene ids at padj < 0.05 cutoff
delta_0_common_genes <- q331k_m323k_cutoffs_overlap %>%
  filter(cutoff_combo == "Q331K_0__M323K_0" & comparison_value == "Gene_Chr") %>%
  pull(intersecting_elements) %>%
  unlist()


# Run wrapper function to generate scatter plot + perform chisq test on quadrant count distributions

m323k_q331k_comparison_delta0 <- twoway_comparison_scatter(dexseq_m323k_targets_pau$`0`,
                                                           dexseq_q331k_targets_pau$`0`,
                                                           dexseq_m323k,
                                                           dexseq_q331k,
                                                           "M323K",
                                                           "Q331K",
                                                           common_ids = delta_0_common_genes
                          )


m323k_q331k_comparison_delta0$scatter_plot
m323k_q331k_comparison_delta0$quadrant_counts
m323k_q331k_comparison_delta0$chisq_result
m323k_q331k_comparison_delta0$chisq_descriptives
m323k_q331k_comparison_delta0$chisq_observed
m323k_q331k_comparison_delta0$chisq_expected

# repeat with greater m323k cutoffs - stronger consistency with Q331K?
m323k_q331k_comparison_delta_5m3 <- twoway_comparison_scatter(dexseq_m323k_targets_pau$`5`,
                                                           dexseq_q331k_targets_pau$`0`,
                                                           dexseq_m323k,
                                                           dexseq_q331k,
                                                           "M323K",
                                                           "Q331K",
                                                           common_ids = unlist(pull(filter(q331k_m323k_cutoffs_overlap,
                                                                                           cutoff_combo == "Q331K_0__M323K_5" & comparison_value == "Gene_Chr"),
                                                                                    intersecting_elements)
                                                                               )
                                                           )

m323k_q331k_comparison_delta_5m3$quadrant_counts

m323k_q331k_comparison_delta0$chisq_result
m323k_q331k_comparison_delta_5m3$chisq_result
m323k_q331k_comparison_delta0$chisq_descriptives
m323k_q331k_comparison_delta_5m3$chisq_descriptives


m323k_q331k_comparison_delta_5b <- twoway_comparison_scatter(dexseq_m323k_targets_pau$`5`,
                                                              dexseq_q331k_targets_pau$`5`,
                                                              dexseq_m323k,
                                                              dexseq_q331k,
                                                              "M323K",
                                                              "Q331K",
                                                              common_ids = unlist(pull(filter(q331k_m323k_cutoffs_overlap,
                                                                                              cutoff_combo == "Q331K_5__M323K_5" & comparison_value == "Gene_Chr"),
                                                                                       intersecting_elements)
                                                              )
                                                             )


m323k_q331k_comparison_delta_5b$quadrant_counts
m323k_q331k_comparison_delta0$chisq_descriptives
m323k_q331k_comparison_delta_5b$chisq_descriptives
m323k_q331k_comparison_delta_5m3$chisq_descriptives


# repeat for all possible combos 
m323k_q331k_comparison_deltas_all <- pmap(list(q331k_cutoff = pau_cutoff_combos$pau_cutoff_q331k,
          m323k_cutoff = pau_cutoff_combos$pau_cutoff_m323k,
          name_cutoff = cutoff_combo_nms
          ),
     function(q331k_cutoff, m323k_cutoff, name_cutoff) {
       
       twoway_comparison_scatter(dexseq_m323k_targets_pau[[as.character(m323k_cutoff)]],
                                 dexseq_q331k_targets_pau[[as.character(q331k_cutoff)]],
                                 dexseq_m323k,
                                 dexseq_q331k,
                                 "M323K",
                                 "Q331K",
                                 common_ids = unlist(pull(filter(q331k_m323k_cutoffs_overlap,
                                                                 cutoff_combo == name_cutoff & comparison_value == "Gene_Chr"),
                                                          intersecting_elements)
                                 )
                                 )
       
       
     },
     .progress = T
     ) %>%
  set_names(cutoff_combo_nms)

# print all plots (annotating with thresholds)
# map2(.x = m323k_q331k_comparison_deltas_all,
#      .y = names(m323k_q331k_comparison_deltas_all),
#      ~ .x$scatter_plot + labs(title = .y))

# M323K-specific targets seem to be over-represented for consistent quadrants - is this restricted to low effect size M323K targets?
map(m323k_q331k_comparison_deltas_all[str_subset(names(m323k_q331k_comparison_deltas_all), "^Q331K_0")],
    ~ .x$chisq_result) %>%
  bind_rows(.id = "delta_cutoffs") %>%
  # adjust with respect to all target set + delta cutoffs (just the Q331K_0 subsets)
  adjust_pvalue(p.col = "p", output.col = "p.adj.all", method = "BH") %>%
  filter(target_set == "M323K")

map(m323k_q331k_comparison_deltas_all[str_subset(names(m323k_q331k_comparison_deltas_all), "^Q331K_0")],
    ~ .x$chisq_descriptives) %>%
  bind_rows(.id = "delta_cutoffs") %>%
  filter(target_set == "M323K")


# As focus on stronger M323K targets, M323K specific targets are no longer over-represented for changing in same direction in Q331K
# Is this just due to noisiness around very low changes in Q331K (i.e. barely affected?)



# Pull out quadrant counts for all Q331K_0 + different M4323K effect sizes
map(m323k_q331k_comparison_deltas_all[str_subset(names(m323k_q331k_comparison_deltas_all), "^Q331K_0")],
    ~ filter(.x$quadrant_counts, target_set == "M323K")
    ) %>%
  bind_rows(.id = "delta_cutoffs") %>%
  mutate(m323k_cutoff = str_remove_all(delta_cutoffs, "Q331K_0__M323K_"),
         m323k_cutoff = fct_inseq(m323k_cutoff),
         consistent = dirn_group %in% c("M323Klonger__Q331Klonger", "M323Kshorter__Q331Kshorter")) %>%
  group_by(m323k_cutoff, consistent) %>%
  summarise(n_tot = sum(n_tot),
            ) %>%
  mutate(frac = n_tot / sum(n_tot))
  
# Whilst proportion of consistent M323K-specific targets does not change between min of 0 + 5, 
# subsequently higher thresholds result in reduced proportion of targets with consistent direction of regulation
# expected vecsue if concordant + large change in M323K, more likely to be found sig by Q331K?

binom_test(x = 152, n = 85 + 152)
binom_test(x = 108, n = 108 + 61)
binom_test(x = 63, n = 103)
binom_test(x = 33, n = 57)
binom_test(x = 18, n = 35)


## Does M323K have a stronger effect in shared targets? 

# plot absolute effect sizes gene-by-gene
m323k_q331k_eff_size_comparison_plot <- m323k_q331k_comparison_deltas_all$Q331K_0__M323K_0$plot_df %>%
  filter(overlap_group == "Both") %>%
  filter(sign(deltaPSImean.HOM_WT.M323K) == sign(deltaPSImean.HOM_WT.Q331K)) %>%
  mutate(dirn = if_else(sign(deltaPSImean.HOM_WT.M323K) == 1, "Lengthening", "Shortening")) %>%
  pivot_longer(cols = starts_with("deltaPSImean"), names_to = "mutant", values_to = "deltaPSImean", names_prefix = "deltaPSImean.HOM_WT.") %>%
  # add max effect size for sorting on plot
  group_by(Gene_Name) %>%
  mutate(max_delta = max(abs(deltaPSImean))) %>%
  ungroup() %>%
  mutate(Gene_Name = fct_reorder(Gene_Name, - max_delta, .fun = max)) %>%
  ggplot(aes(x = Gene_Name, y = abs(deltaPSImean), colour = mutant)) +
  geom_point() +
  scale_colour_manual(values = c("#1b9e77", "#7570b3")) +
  theme_bw(base_size = 14) +
  labs(x = "APA Target Gene",
       y = "Absolute Difference in Mean psi (HOM - WT)",
       colour = "Mutant") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90,
                                   hjust = 0.5,
                                   vjust = 0.5))

m323k_q331k_eff_size_comparison_plot


# Do 
# what about for psi analysis?
# same PAS = exactly the same PAS between two comparisons
# same number of PAS = same number of PAS (delta psis are more comparable)


# Subset to Q331K/M323K targets, extract all analysed/evaluated PAS
int_psi_dexseq_m323k <- m323k_q331k_comparison_deltas_all$Q331K_0__M323K_0$plot_df %>%
  distinct(Gene, Chr, overlap_group) %>%
  left_join(dexseq_m323k, by = c("Gene", "Chr")) %>%
  mutate(pas_sig = padj < 0.05) %>%
  distinct(Gene, Chr, APA_ID, overlap_group, pas_sig)

int_psi_dexseq_q331k <- m323k_q331k_comparison_deltas_all$Q331K_0__M323K_0$plot_df %>%
  distinct(Gene, Chr, overlap_group) %>%
  left_join(dexseq_q331k, by = c("Gene", "Chr")) %>%
  mutate(pas_sig = padj < 0.05) %>%
  distinct(Gene, Chr, APA_ID, overlap_group, pas_sig)

# Join all APA_IDs between datasets 
psi_dexseq_joined <- full_join(int_psi_dexseq_m323k,
                               int_psi_dexseq_q331k,
                               by = c("Gene", "Chr", "APA_ID"),
                               suffix = c(".m323k", ".q331k")
) %>%
  mutate(overlap_group = if_else(is.na(overlap_group.m323k), overlap_group.q331k, overlap_group.m323k)) %>%
  select(-overlap_group.m323k, -overlap_group.q331k)

psi_dexseq_joined

# Do genes have same PAS / number of PAS analysed between the two conditions (i.e. directly comparable psi values?)
psi_summary_joined <- psi_dexseq_joined %>%
  mutate(pas_in_both = !is.na(pas_sig.m323k) & !is.na(pas_sig.q331k)) %>%
  group_by(Gene, Chr, overlap_group) %>%
  summarise(n_pas.m323k = sum(!is.na(pas_sig.m323k)),
            n_pas.q331k = sum(!is.na(pas_sig.q331k)),
            same_n_pas = n_pas.m323k == n_pas.q331k,
            same_pas = all(pas_in_both)) %>%
  ungroup()

psi_summary_joined

psi_counts_summary_joined <- psi_summary_joined %>%
  count(overlap_group, same_pas, same_n_pas) %>%
  group_by(overlap_group) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

psi_counts_summary_joined

# For genes without identical analysed PAS, do they at least share some identical PAS (i.e. could use PAU to evaluate)
psi_summary_joined %>%
  filter(!(same_pas & same_n_pas)) %>%
  # remove genes not expressed in one condition
  filter(n_pas.m323k != 0 & n_pas.q331k != 0) %>%
  left_join(select(psi_dexseq_joined, -overlap_group), by = c("Gene", "Chr"))
  
  
psi_summary_joined %>%
  filter(overlap_group == "Both") %>%
  select(-overlap_group) %>%
  left_join(m323k_q331k_comparison_deltas_all$Q331K_0__M323K_0$plot_df, by = c("Gene", "Chr")) %>%
  mutate(dirn = if_else(sign(deltaPSImean.HOM_WT.M323K) == 1, "Lengthening", "Shortening")) %>%
  ggplot(aes(x = deltaPSImean.HOM_WT.M323K, y = deltaPSImean.HOM_WT.Q331K, alpha = plot_alpha, colour = same_pas)) +
  scale_x_continuous(limits = c(-0.4, 0.4)) +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

psi_summary_joined %>%
  filter(overlap_group == "Both") %>%
  select(-overlap_group) %>%
  left_join(m323k_q331k_comparison_deltas_all$Q331K_0__M323K_0$plot_df, by = c("Gene", "Chr")) %>%
  mutate(dirn = if_else(sign(deltaPSImean.HOM_WT.M323K) == 1, "Lengthening", "Shortening")) %>%
  ggplot(aes(x = abs(deltaPSImean.HOM_WT.M323K), y = abs(deltaPSImean.HOM_WT.Q331K), colour = same_pas, shape = dirn)) +
  scale_x_continuous(limits = c(0, 0.25)) +
  scale_y_continuous(limits = c(0, 0.25)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw(base_size = 14)

psi_dexseq_joined %>%
  left_join(distinct(dexseq_m323k, Gene, Chr, APA_ID, Pas_Number, delta_HOM_WT.m323k = delta_HOM_WT), by = c("Gene", "Chr", "APA_ID")) %>%
  left_join(distinct(dexseq_q331k, Gene, Chr, APA_ID, delta_HOM_WT.q331k = delta_HOM_WT), by = c("Gene", "Chr", "APA_ID")) %>%
  filter(overlap_group == "Both",
         sign(delta_HOM_WT.m323k) == sign(delta_HOM_WT.q331k)) %>%
  group_by(Gene, Chr) %>%
  mutate(prox = Pas_Number == min(Pas_Number)) %>%
  ungroup() %>%
  filter(prox) %>%
  ggplot(aes(x = delta_HOM_WT.m323k, y = delta_HOM_WT.q331k, shape = prox)) +
  scale_x_continuous(limits = c(-50, 50)) +
  scale_y_continuous(limits = c(-50, 50)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14)

# Using PAU of prox to assess this
psi_dexseq_joined %>%
  left_join(distinct(dexseq_m323k, Gene, Chr, APA_ID, Pas_Number, delta_HOM_WT.m323k = delta_HOM_WT), by = c("Gene", "Chr", "APA_ID")) %>%
  left_join(distinct(dexseq_q331k, Gene, Chr, APA_ID, delta_HOM_WT.q331k = delta_HOM_WT), by = c("Gene", "Chr", "APA_ID")) %>%
  filter(overlap_group == "Both",
         sign(delta_HOM_WT.m323k) == sign(delta_HOM_WT.q331k)) %>%
  group_by(Gene, Chr) %>%
  mutate(prox = Pas_Number == min(Pas_Number),
         dirn = as.factor(sign(delta_HOM_WT.m323k))) %>%
  ungroup() %>%
  filter(prox) %>%
  ggplot(aes(x = abs(delta_HOM_WT.m323k), y = abs(delta_HOM_WT.q331k), colour = dirn)) +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(0, 50)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14)

#dist
psi_dexseq_joined %>%
  left_join(distinct(dexseq_m323k, Gene, Chr, APA_ID, Pas_Number, delta_HOM_WT.m323k = delta_HOM_WT), by = c("Gene", "Chr", "APA_ID")) %>%
  left_join(distinct(dexseq_q331k, Gene, Chr, APA_ID, delta_HOM_WT.q331k = delta_HOM_WT), by = c("Gene", "Chr", "APA_ID")) %>%
  filter(overlap_group == "Both",
         sign(delta_HOM_WT.m323k) == sign(delta_HOM_WT.q331k)) %>%
  group_by(Gene, Chr) %>%
  mutate(dist = Pas_Number == max(Pas_Number),
         dirn = as.factor(sign(delta_HOM_WT.m323k))) %>%
  ungroup() %>%
  filter(dist) %>%
  ggplot(aes(x = abs(delta_HOM_WT.m323k), y = abs(delta_HOM_WT.q331k), colour = dirn)) +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(0, 50)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(title = "dist")

  
#  
# SAVE PLOTS TO DISK
#

# pvalue histogram
ggsave(filename = file.path(outdir, "2024-06-30_combined_pvalue_histogram.png"),
       plot = pval_hist_comb,
       units = "mm",
       width = 200,
       height = 100
       )

ggsave(filename = file.path(outdir, "2024-06-30_combined_pvalue_histogram.pdf"),
       plot = pval_hist_comb,
       units = "mm",
       width = 200,
       height = 100
       )

# counts of significant PAS
ggsave(filename = file.path(outdir, "2024-06-30_sig_gene_count_delta_thresholds.bar_plot.png"),
       plot = sig_gene_delta_count_plot + theme(legend.position = "top") + labs(x = "Minimum Absolute Change in Usage (HOM - WT, %)"),
       units = "mm",
       width = 150,
       height = 150       
       )

ggsave(filename = file.path(outdir, "2024-06-30_sig_gene_count_delta_thresholds.bar_plot.pdf"),
       plot = sig_gene_delta_count_plot + theme(legend.position = "top") + labs(x = "Minimum Absolute Change in Usage (HOM - WT, %)"),
       units = "mm",
       width = 150,
       height = 150       
)


# psi bar plots for each cutoff
walk2(.x = psi_bar_plots_p,
      .y = names(psi_bar_plots_p),
      ~ ggsave(filename = file.path(outdir,
                                    paste("2024-06-30_psi_dirn_bar",
                                          "min_abs_delta",
                                          .y,
                                          "png",
                                          sep = ".")
                                    ),
               plot = .x + labs(title = NULL) + theme(legend.position = "top"),
               units = "mm",
               height = 150,
               width = 150
               )
      )

walk2(.x = psi_bar_plots_p,
      .y = names(psi_bar_plots_p),
      ~ ggsave(filename = file.path(outdir,
                                    paste("2024-06-30_psi_dirn_bar",
                                          "min_abs_delta",
                                          .y,
                                          "pdf",
                                          sep = ".")
      ),
      plot = .x + labs(title = NULL) + theme(legend.position = "top"),
      units = "mm",
      height = 150,
      width = 150
      )
)

# combined psi across delta thresholds
ggsave(filename = file.path(outdir, "2024-06-30_psi_dirn_bar.all_deltas.png"),
       plot = psi_comb_delta_thresh_bar + theme(legend.position = "top"),
       units = "mm",
       height = 125,
       width = 175
       )

ggsave(filename = file.path(outdir, "2024-06-30_psi_dirn_bar.all_deltas.pdf"),
       plot = psi_comb_delta_thresh_bar + theme(legend.position = "top"),
       units = "mm",
       height = 125,
       width = 175
)

# psi scatter 
ggsave(filename = file.path(outdir,
                            "2024-06-30_psi_comparison_scatter.targets.no_delta_threshold.png"),
       plot = m323k_q331k_comparison_deltas_all$Q331K_0__M323K_0$scatter_plot + labs(x = "Change in psi (HOM - WT) - M323K",
                                                                                     y = "Change in psi (HOM - WT) - Q331K") +
         theme(legend.position = "top"),
       dpi = "retina",
       units = "mm",
       height = 125,
       width = 125
       )

ggsave(filename = file.path(outdir,
                            "2024-06-30_psi_comparison_scatter.targets.no_delta_threshold.pdf"),
       plot = m323k_q331k_comparison_deltas_all$Q331K_0__M323K_0$scatter_plot + labs(x = "Change in psi (HOM - WT) - M323K",
                                                                                     y = "Change in psi (HOM - WT) - Q331K") +
         theme(legend.position = "top"),
       dpi = "retina",
       units = "mm",
       height = 125,
       width = 125
)


# euler plot
ggsave(filename = file.path(outdir,
                            "2024-07-01_q331k_m323k_euler.all.no_delta_threshold.png"),
       plot = q331k_m323k_euler$all_plot,
       dpi = "retina",
       units = "mm",
       height = 125,
       width = 125
       )

ggsave(filename = file.path(outdir,
                            "2024-07-01_q331k_m323k_euler.all.no_delta_threshold.pdf"),
       plot = q331k_m323k_euler$all_plot,
       dpi = "retina",
       units = "mm",
       height = 125,
       width = 125
)


# common target effect size comparison plot
ggsave(filename = file.path(outdir,
                            "2024-07-01_q331k_m323k_common_target.effect_size_comparison_point.png"),
       plot = m323k_q331k_eff_size_comparison_plot,
       dpi = "retina",
       units = "mm",
       height = 125,
       width = 125
)

ggsave(filename = file.path(outdir,
                            "2024-07-01_q331k_m323k_common_target.effect_size_comparison_point.pdf"),
       plot = m323k_q331k_eff_size_comparison_plot,
       dpi = "retina",
       units = "mm",
       height = 125,
       width = 125
)

