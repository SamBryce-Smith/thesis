library(tidyverse)
source("scripts/dosage_utils.R")
source("scripts/comparison_utils.R")

# processed PAU dfs produced by calculate_qapa_means_psis.R
pau_m323k <- read_tsv("processed/2024-02-29/M323K_Adult_Brain/M323K_Adult_Brain.mean_pau.dosage_labels.tsv")

psi_m323k <- read_tsv("processed/2024-02-29/M323K_Adult_Brain/M323K_Adult_Brain.psi.mean_median.dosage_labels.tsv")

# processed DEXSeq results tables filtered to representative PAS produced by mef_representative_pas.R
dexseq_rep_m323k <- read_tsv("processed/2024-02-29/M323K_Adult_Brain/M323K_Adult_Brain.representative.pas.tsv")

outdir <- "processed/2024-02-29/analysis/m323k_adult_brain"
if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

# Join in DEXSeq with PAU results
dexseq_rep_m323k <- pau_m323k %>%
  # already have PAUs as column names
  select(-condition) %>%
  distinct(APA_ID, .keep_all = T) %>%
  left_join(dexseq_rep_m323k, ., by = c("APA_ID", "Gene_Name")) %>%
  left_join(psi_m323k, by = c("Gene", "Chr"))


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

# construct 'strict' dfs, where all PAS have padj < 0.05
dexseq_rep_m323k_prox_padj_strict <- filter(dexseq_rep_m323k_prox_padj, padj < 0.05)
dexseq_rep_m323k_dist_padj_strict <- filter(dexseq_rep_m323k_dist_padj, padj < 0.05)

prox_specific <- setdiff(unique(dexseq_rep_m323k_prox_padj_strict$Gene), unique(dexseq_rep_m323k_dist_padj_strict$Gene))
dist_specific <- setdiff(unique(dexseq_rep_m323k_dist_padj_strict$Gene), unique(dexseq_rep_m323k_prox_padj_strict$Gene))

#
# HOW DOES NUMBER OF GENES VARY WITH A MINIMUM EFFECT SIZE THRESHOLD?
#


sig_gene_delta_count_m323k <- c(0,5,10,15,20) %>%
  set_names() %>%
  map(~ dexseq_rep_m323k_padj %>%
        filter(padj < 0.05) %>%
        filter(abs(delta_HOM_WT) >= .x) %>%
        distinct(Gene, Chr)) %>%
  bind_rows(.id = "abs_delta_cutoff") %>%
  mutate(abs_delta_cutoff = as.integer(abs_delta_cutoff)) %>%
  count(abs_delta_cutoff)

sig_gene_delta_count_plot_m323k <- sig_gene_delta_count_m323k %>%
  ggplot(aes(x = abs_delta_cutoff, y = n, label = n)) +
  geom_col() +
  geom_text(nudge_y = 10) +
  scale_y_continuous(breaks = seq(0,500,50)) +
  theme_bw(base_size = 14) +
  labs(title = "M323K 6 month (Fratta et al.)",
       subtitle = "PAS padj < 0.05", 
       x = "Minimum absolute change in usage (HOM - WT, %)",
       y = "Gene Count")

sig_gene_delta_count_plot_m323k


# First, make standard dosage plots with no delta threshold
# here, list names will be passed as subtitles as MUT the same across datasets
plot_df_list <- list("Proximal PAS, Gene any padj < 0.05" = dexseq_rep_m323k_prox_padj,
               "Proximal PAS, padj < 0.05" = dexseq_rep_m323k_prox_padj_strict,
               "Distal PAS, Gene any padj < 0.05" = dexseq_rep_m323k_dist_padj,
               "Distal PAS, padj < 0.05" = dexseq_rep_m323k_dist_padj_strict)

# Make plots for all dosage categories
m323k_ad6m_all <- map2(.x = plot_df_list,
                          .y = names(plot_df_list),
                          ~ .x %>%
                            dosage_analysis_wrapper(line_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                                                     subtitle = .y
                                                                     ),
                                                    bar_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                                                    subtitle = .y
                                                                    )
                                                    )
)


# Make plots subsetted to 'clear'/bonafide consistent/inconsistent dosage categories
m323k_ad6m_subset <- map2(.x = plot_df_list,
                     .y = names(plot_df_list),
                     ~ dosage_analysis_wrapper(df = .x,
                                               patterns_filter = meanPAU.dosage_pattern %in% c("Up+Up", "Down+Down", "Down+Up", "Up+Down"),
                                               line_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                                                subtitle = .y),
                                               bar_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                                               subtitle = .y)
                                               )
                     )


# make plots subsetted to same-containing categories
m323k_ad6m_subset_same <- map2(.x = plot_df_list,
                          .y = names(plot_df_list),
                          ~ dosage_analysis_wrapper(df = .x,
                                                    patterns_filter = meanPAU.dosage_pattern %in% c("Down+Same", "Up+Same", "Same+Up", "Same+Down"),
                                                    line_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                                                     subtitle = .y),
                                                    bar_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                                                    subtitle = .y)
                          )
)


map(m323k_ad6m_all, ~ .x$bar_plot)
map(m323k_ad6m_subset, ~ .x$bar_plot)
map(m323k_ad6m_subset, ~ .x$dosage_plot)

# so few 
map(m323k_ad6m_subset_same, ~ .x$dosage_plot)
map(m323k_ad6m_subset_same, ~ .x$bar_plot)

# Make plots at different abs delta thresholds
m323k_ad6m_delta_all <- pmap(list(.x = rep(plot_df_list, 5),
          .y = rep(names(plot_df_list), 5),
          thresh = sort(rep(seq(0,20,5), 4))
          ),
     function(.x, .y, thresh) .x %>%
       filter(abs(delta_HOM_WT) > thresh) %>%
       dosage_analysis_wrapper(line_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                                subtitle = paste(.y, " & abs(HOM - WT) > ", thresh, " %", sep = "")),
                               bar_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                               subtitle = paste(.y, " & abs(HOM - WT) > ", thresh, " %", sep = ""))
       )
          
) %>%
  # add delta threshold to 
  set_names(paste(names(.), sort(rep(seq(0,20,5), 4)), sep = " - "))

# Make plots at different abs delta thresholds, this time for clear dosage categories
m323k_ad6m_delta_subset <- pmap(list(.x = rep(plot_df_list, 5),
                 .y = rep(names(plot_df_list), 5),
                 thresh = sort(rep(seq(0,20,5), 4))
),
function(.x, .y, thresh) .x %>%
  filter(abs(delta_HOM_WT) > thresh) %>%
  dosage_analysis_wrapper(patterns_filter = meanPAU.dosage_pattern %in% c("Up+Up", "Down+Down", "Down+Up", "Up+Down"),
                          line_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                           subtitle = paste(.y, " & abs(HOM - WT) > ", thresh, " %", sep = "")),
                          bar_labs = labs(title = "M323K 6 month (Fratta et al.)",
                                          subtitle = paste(.y, " & abs(HOM - WT) > ", thresh, " %", sep = ""))
  )

) %>%
  # add delta threshold to 
  set_names(paste(names(.), sort(rep(seq(0,20,5), 4)), sep = " - "))

map(m323k_ad6m_delta_subset, ~ .x$dosage_plot)$`Distal PAS, padj < 0.05 - 5`
map(m323k_ad6m_delta_subset, ~ .x$bar_plot)$`Distal PAS, padj < 0.05 - 5`

# Does bias for consistency improve/decline with increasing absolute effect size thresholds?

# extract count dfs for proximal + distal strict sig at each usage threshold
m323k_ad6m_delta_all_count_dfs <- map(m323k_ad6m_delta_all[str_subset(names(m323k_ad6m_delta_all),
                                                                      ", padj < 0.05")
                                                           ],
                                      "counts_df")

m323k_ad6m_delta_all_count_dfs

str_subset(names(m323k_ad6m_delta_all), ", padj < 0.05")

# add the pas type + sig stringency pairing as a column to df
# # then add the abs delta cutoff
m323k_ad6m_delta_all_count_dfs_comb <- map2(.x = m323k_ad6m_delta_all_count_dfs,
                                            .y = names(m323k_ad6m_delta_all_count_dfs),
                                            ~ mutate(.x, group = str_split_i(.y, " - ", 1)
                                                     )
                                            ) %>%
  # update names to just the delta cutoff
  set_names(sort(rep(seq(0,20,5), 2))) %>%
  bind_rows(.id = "delta_cutoff") %>%
  mutate(delta_cutoff = fct_inseq(delta_cutoff)
         )

# Reassign consistent vs not, recalculate counts + %
m323k_ad6m_delta_all_count_dfs_comb_summ <- m323k_ad6m_delta_all_count_dfs_comb %>%
  mutate(dosage_pattern_simple = case_when(grepl("Up", meanPAU.dosage_pattern) & grepl("Down", meanPAU.dosage_pattern) ~ "Inconsistent",
                                           meanPAU.dosage_pattern %in% c("Up+Up", "Down+Down") ~ "Consistent",
                                           str_detect(meanPAU.dosage_pattern, "Same") ~ "Same_containing",
                                           T ~ "other"
                                           )
         ) %>%
  group_by(delta_cutoff, group, dosage_pattern_simple) %>%
  summarise(n = sum(n)) %>%
  # dosage pattern dropped from grouping, now calc perc w/ respect to cutoff + pas pairing
  mutate(perc = (n / sum(n))*100) %>%
  ungroup()

# Generate plot
dosage_consistency_delta_bar <- m323k_ad6m_delta_all_count_dfs_comb_summ %>%
  # Sort so Proximal first, followed by distal
  mutate(group = factor(group, levels = sort(unique(group), decreasing = T))) %>%
  filter(dosage_pattern_simple != "Same_containing") %>%
  ggplot(aes(x = dosage_pattern_simple, y = perc, label = n, fill = delta_cutoff)) +
  facet_wrap("~ group", scales = "fixed") +
  geom_col(position = "dodge") +
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  scale_fill_brewer(type = "seq", palette = "Greens") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x = "PAS Usage Pattern (WT - HET - HOM)", y = "Event Frequency (%)",
       fill = "Minimum abs(HOM - WT)"
       )

dosage_consistency_delta_bar


#
# Directionality of target bias using LABRAT psi
#

# Define groups of high confidence targets based on dosage dependency
# (list of dfs)
dexseq_rep_m323k_targets <- define_apa_targets(dexseq_rep_m323k)

# Want to assess bias at increasing absolute effect size thresholds for each set of target types
# want - combined list of all rep tartgets at each threshold,

# first define vector target labels combined with absolute PAU cutoff
target_cutoff_nms <- expand_grid(cutoff = sort(seq(0,20,5)), target = names(dexseq_rep_m323k_targets)) %>%
  reframe(nm = paste(cutoff, target, sep = " - ")) %>%
  pull()

# constuct a list of all possible targets + PAU cutoff thresholds (i.e. subsetted dfs)
dexseq_rep_m323k_targets_pau <- map2(.x = sort(rep(seq(0,20,5),4)),
                                     .y = rep(names(dexseq_rep_m323k_targets), 5),
                                     ~ filter(dexseq_rep_m323k_targets[[.y]], abs(delta_HOM_WT) >= .x)
) %>%
  set_names(target_cutoff_nms)

# Run separately for each target set (produces combined plot with x axis = absolute delta threshold)
psi_pau_thresh_gene <- psi_direction_bar_plot(dexseq_rep_m323k_targets_pau[str_subset(names(dexseq_rep_m323k_targets_pau), "gene")], id_column = "abs_delta_threshold")
psi_pau_thresh_proximal <- psi_direction_bar_plot(dexseq_rep_m323k_targets_pau[str_subset(names(dexseq_rep_m323k_targets_pau), "proximal")], id_column = "abs_delta_threshold")
psi_pau_thresh_distal <- psi_direction_bar_plot(dexseq_rep_m323k_targets_pau[str_subset(names(dexseq_rep_m323k_targets_pau), "distal")], id_column = "abs_delta_threshold")
psi_pau_thresh_all <- psi_direction_bar_plot(dexseq_rep_m323k_targets_pau[str_subset(names(dexseq_rep_m323k_targets_pau), "all")], id_column = "abs_delta_threshold")

# revise plot labels (as typically intended for comparing two MUTs)
psi_pau_thresh_gene$plot <- psi_pau_thresh_gene$plot +
  labs(subtitle = "gene",
       x = "Minimum Absolute (HOM - WT) Threshold (%)") +
  scale_x_discrete(labels = seq(0,20,5))

psi_pau_thresh_proximal$plot <- psi_pau_thresh_proximal$plot +
  labs(subtitle = "proximal",
       x = "Minimum Absolute (HOM - WT) Threshold (%)") +
  scale_x_discrete(labels = seq(0,20,5))

psi_pau_thresh_distal$plot <- psi_pau_thresh_distal$plot +
  labs(subtitle = "distal",
       x = "Minimum Absolute (HOM - WT) Threshold (%)") +
  scale_x_discrete(labels = seq(0,20,5))

psi_pau_thresh_all$plot <- psi_pau_thresh_all$plot +
  labs(subtitle = "all",
       x = "Minimum Absolute (HOM - WT) Threshold (%)") +
  scale_x_discrete(labels = seq(0,20,5))


psi_pau_thresh_gene$plot
psi_pau_thresh_proximal$plot  
psi_pau_thresh_distal$plot
psi_pau_thresh_all$plot


psi_pau_thresh_gene$binom_test

# Generate plot for single threshold + target type, extract the counts df & binom_df, subset and borrow plotting code
all_binom <- filter(psi_pau_thresh_all$binom_test, abs_delta_threshold == "5 - all") %>%
  mutate(mut = "M323K",
         y.position =  max(estimate * 100, (1 - estimate) * 100) + 10) # adjust y position for this cutoff only

all_binom

psi_plot_all_min5 <- psi_pau_thresh_all$counts_df %>%
  filter(abs_delta_threshold == "5 - all") %>%
  mutate(mut = "M323K") %>%
  ggplot(aes(x = mut, y = perc, label = n)) +
  geom_col(aes(fill = dirn_mean), position = "dodge") +
  geom_text(aes(group = dirn_mean), position = position_dodge(width = 0.9), vjust = -1) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
  ggprism::add_pvalue(all_binom, label = "p.adj.signif", x = "mut", label.size = 6) +
  scale_fill_manual(values = c("#b2df8a", "#1f78b4")) +
  labs(title = NULL,
       x = NULL,
       y = "Target APAs (%)",
       fill = "3'UTR length change") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

psi_plot_all_min5  


#
# Save plots to disk
#

# gene count bar
ggsave(filename = file.path(outdir, "2024-06-30_sig_gene_delta_bar.png"),
       plot = sig_gene_delta_count_plot_m323k + labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150
       )

ggsave(filename = file.path(outdir, "2024-06-30_sig_gene_delta_bar.pdf"),
       plot = sig_gene_delta_count_plot_m323k  + labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150
)

# dosage and bar plot at single abs threshold

# proximal - 5 %
ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.line_plot.proximal.png"),
       plot = m323k_ad6m_delta_subset$`Proximal PAS, padj < 0.05 - 5`$dosage_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
       )

ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.line_plot.proximal.pdf"),
       plot = m323k_ad6m_delta_subset$`Proximal PAS, padj < 0.05 - 5`$dosage_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
)

# Repeat for distal
ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.line_plot.distal.png"),
       plot = m323k_ad6m_delta_subset$`Distal PAS, padj < 0.05 - 5`$dosage_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
)

ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.line_plot.distal.pdf"),
       plot = m323k_ad6m_delta_subset$`Distal PAS, padj < 0.05 - 5`$dosage_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
)

# Save bar plots
ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.bar_plot.proximal.png"),
       plot = m323k_ad6m_delta_subset$`Proximal PAS, padj < 0.05 - 5`$bar_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
)

ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.bar_plot.proximal.pdf"),
       plot = m323k_ad6m_delta_subset$`Proximal PAS, padj < 0.05 - 5`$bar_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
)

# Repeat for distal
ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.bar_plot.distal.png"),
       plot = m323k_ad6m_delta_subset$`Distal PAS, padj < 0.05 - 5`$bar_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
)

ggsave(filename = file.path(outdir,
                            "2024-06-30_dosage.category_subset.min_pau_5.bar_plot.distal.pdf"),
       plot = m323k_ad6m_delta_subset$`Distal PAS, padj < 0.05 - 5`$bar_plot +
         labs(title = NULL),
       units = "mm",
       width = 150,
       height = 150       
)

# dosage consistency across pau thresholds
ggsave(filename = file.path(outdir, "2024-06-30_dosage.category_simple.pau_thresholds.bar_plot.png"),
       plot = dosage_consistency_delta_bar,
       units = "mm",
       width = 200,
       height = 200       
       )

ggsave(filename = file.path(outdir, "2024-06-30_dosage.category_simple.pau_thresholds.bar_plot.pdf"),
       plot = dosage_consistency_delta_bar,
       units = "mm",
       width = 200,
       height = 200       
       )


# psi plot at single threshold
ggsave(filename = file.path(outdir, "2024-06-30_psi.min_pau_5.dirn_bar_plot.png"),
       plot = psi_plot_all_min5,
       units = "mm",
       width = 150,
       height = 150       
       )

ggsave(filename = file.path(outdir, "2024-06-30_psi.min_pau_5.dirn_bar_plot.pdf"),
       plot = psi_plot_all_min5,
       units = "mm",
       width = 150,
       height = 150       
)

# psi plot across delta thresholds
ggsave(filename = file.path(outdir, "2024-06-30_psi.pau_thresholds.dirn_bar_plot.png"),
       plot = psi_pau_thresh_all$plot + labs(title = NULL, subtitle = NULL),
       units = "mm",
       width = 150,
       height = 150 
)

ggsave(filename = file.path(outdir, "2024-06-30_psi.pau_thresholds.dirn_bar_plot.pdf"),
       plot = psi_pau_thresh_all$plot + labs(title = NULL, subtitle = NULL),
       units = "mm",
       width = 150,
       height = 150 
)
