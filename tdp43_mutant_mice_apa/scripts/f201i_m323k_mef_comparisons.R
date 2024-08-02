library(tidyverse)
library(ggVennDiagram)


#' Select representative polyA site according to Hallegger et al. criteria
#'
#' Grouping by gene, this approach first selects rows (polyA sites) with the smallest p-value. If there are ties within a group, the row with the 
#' highest usage value is chosen as the representative. 
#'
#' @param df A data frame containing the QAPA results
#' @param group_cols A character vector specifying the column names based on 
#'                   which the grouping is done. Default is c("Gene_Name", "Chr").
#' @param pvalue_col A character string specifying the column name containing 
#'                   the p-values. Default is "padj".
#' @param usage_col A character string specifying the column name containing 
#'                  the polyA usage values. Default is "mean_PAU_HOM".
#'
#' @return A data frame containing a representative row/polyA site for each group.
#' @export
#'
#' @examples
#' df <- data.frame(Gene_Name = c("Gene1", "Gene1", "Gene2", "Gene2"),
#'                  Chr = c("chr1", "chr1", "chr2", "chr2"),
#'                  padj = c(0.01, 0.05, 0.02, 0.03),
#'                  mean_PAU_HOM = c(0.8, 0.9, 0.7, 0.6))
#' select_rep_pas_halleger(df)
select_rep_pas_halleger <- function(df, group_cols = c("Gene_Name", "Chr"), pvalue_col = "padj", usage_col = "mean_PAU_HOM") {
  
  # for each group, select event/row with smallest pvalue
  tmp <- df %>%
    group_by(across(all_of(group_cols))) %>%
    slice_min(!!sym(pvalue_col), with_ties = T)
  
  # check if pvalue ties within group, if so select most used pas as representatives
  if (nrow(group_keys(tmp)) == nrow(tmp)) {
    # All groups have single site selected
    tmp %>%
      ungroup()
  } else {
    # select site with highest usage as representative
    tmp %>%
      slice_max(!!sym(usage_col)) %>%
      ungroup()
  }
}


#' Select representative polyA site according to most proximal **expressed** site 
#'
#' Grouping by gene, this approach, this function selects the row/event with the smallest Pas_Number i.e. the most gene-start proximal site present in the dataframe.
#'
#' @param df A data frame containing the dataset.
#' @param group_cols A character vector specifying the column names based on 
#'                   which the grouping is done. Default is c("Gene_Name", "Chr").
#' @param pas_num_col A character string specifying the column name containing 5'-3' 1..n polyA site annotation
#'                    Default is "Pas_Number".
#'
#' @return A data frame containing representative rows for each group.
#' @export
#'
#' @examples
#' df <- data.frame(Gene_Name = c("Gene1", "Gene1", "Gene2", "Gene2"),
#'                  Chr = c("chr1", "chr1", "chr2", "chr2"),
#'                  Pas_Number = c(1, 2, 1, 2))
#' select_rep_pas_mostproximal(df)
select_rep_pas_mostproximal <- function(df, group_cols = c("Gene_Name", "Chr"), pas_num_col = "Pas_Number") {
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    slice_min(!!sym(pas_num_col)) %>%
    ungroup()
  
}


#' Select representative polyA site according to most proximal annotated site
#'
#' This approach selects rows where Pas_Number == 1 i.e. it is the most gene-start proximal **annotated** event for that gene 
#'
#' @param df A data frame containing the dataset.
#' @param pas_num_col A character string specifying the column name containing 5'-3' 1..n polyA site annotation
#'                    Default is "Pas_Number".
#'
#' @return A data frame containing representative rows for each group.
#' @export
#'
#' @examples
#' df <- data.frame(Gene_Name = c("Gene1", "Gene1", "Gene2", "Gene2"),
#'                  Chr = c("chr1", "chr1", "chr2", "chr2"),
#'                  Pas_Number = c(1, 2, 1, 2))
#' select_rep_pas_proximal(df)
select_rep_pas_proximal <- function(df, pas_num_col = "Pas_Number") {
  
  df %>%
    filter(!!sym(pas_num_col) == 1)
}


#' Prepare dataframe for base_dosage_plot
#'
#' Prepare a wide-format dataframe for generating dosage dependent PAS usage plots (with base_dosage_plot)
#'
#' @param df input dataframe
#' @param col_prefix The prefix of the columns containing per-condition usages.
#' @param names_col The name of the output column containing condition names.
#' @param values_col The name of the output column containing usage values.
#' @param pattern_col The name of the column containing dosage pattern categories.
#' @param dosage_levels A vector specifying the dosage levels/experimental condition keys and their low-high dosage order.
#' @param dosage_patterns A vector specifying the desired plotting order of dosage pattern labels
#' @return A dataframe ready for use with base_dosage_plot
prep_dosage_plot_df <- function(df, col_prefix, names_col, values_col, pattern_col, dosage_levels = c("WT", "HET", "HOM"), 
                                dosage_patterns = c("Up+Up", "Down+Down", "Up+Same", "Down+Same", 
                                                    "Same+Down", "Same+Up", "Up+Down", "Down+Up", "Other", "NA")) {
  df_processed <- df %>%
    pivot_longer(starts_with(col_prefix), names_to = names_col, values_to = values_col, names_prefix = col_prefix) %>%
    mutate({{ names_col }} := factor(!!sym(names_col), levels = dosage_levels),
           {{ pattern_col }} := factor(!!sym(pattern_col), levels = dosage_patterns)
    )
  
  return(df_processed)
}


#' Generate minimal dosage dependent usage line plot
#'
#' Generate a base ggplot object for visualizing dosage-dependent PAS usage patterns
#'
#' @param data A dataframe containing the data to be plotted (e.g. output of prep_dosage_plot_df. Assumes that have 1 row per condition (x_var), and is ordered factor in lowest-highest dosage
#' @param x_var The name of the variable to be plotted on the x-axis.
#' @param y_var The name of the variable to be plotted on the y-axis.
#' @param group_var The name of the variable to be used for grouping data.
#' @param facet_var The name of the variable to be used for facet wrapping.
#' @param alpha_val The transparency of lines (0 to 1).
#' @param base_font_size The base font size for the plot.
#' @return A ggplot object.
base_dosage_plot <- function(data, x_var, y_var, group_var, facet_var, alpha_val, base_font_size) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], group = .data[[group_var]])) +
    facet_wrap(as.formula(paste0("~", facet_var))) +
    geom_line(alpha = alpha_val) +
    theme_bw(base_size = base_font_size)
}



##### ANALYSIS

# processed PAU dfs produced by calculate_qapa_means_psis.R
pau_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.mean_pau.dosage_labels.tsv")
pau_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.mean_pau.dosage_labels.tsv")

# psi tables containing condition means and deltas
psi_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.psi.mean_median.dosage_labels.tsv")
psi_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.psi.mean_median.dosage_labels.tsv")

# processed DEXSeq results tables produced by QAPA_snakemake pipeline
dexseq_f210i <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")
dexseq_m323k <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")


outdir <- "processed/2024-02-29/analysis/mef_comparisons"
if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

# combine PAUs with DEXSeq results

dexseq_m323k <- pau_m323k %>%
  # already have PAUs as column names
  select(-condition) %>%
  distinct(APA_ID, .keep_all = T) %>%
  left_join(dexseq_m323k, ., by = c("APA_ID", "Gene_Name"))
  
dexseq_f210i <- pau_f210i %>%
  # already have PAUs as column names
  select(-condition) %>%
  distinct(APA_ID, .keep_all = T) %>%
  left_join(dexseq_f210i, ., by = c("APA_ID", "Gene_Name"))

# filter for padj < 0.05
dexseq_m323k_padj <- filter(dexseq_m323k, padj < 0.05)
dexseq_f210i_padj <- filter(dexseq_f210i, padj < 0.05)


# create a minimal df of two MUTs combined

# all expressed/evaluated PAS
comb <- bind_rows("F210I" = dexseq_f210i, "M323K" = dexseq_m323k, .id = "mutant") %>%
  distinct(mutant, Gene_Name, APA_ID, Pas_Number)

# filtered for significance
comb_padj <- bind_rows("F210I" = dexseq_f210i_padj, "M323K" = dexseq_m323k_padj, .id = "mutant") %>%
  distinct(mutant, Gene_Name, APA_ID, Pas_Number)

# number of PAS per MUT
count(comb, mutant)
count(comb_padj, mutant)

# Number of genes per MUT
distinct(comb, mutant, Gene_Name) %>%
  count(mutant)
distinct(comb_padj, mutant, Gene_Name) %>%
  count(mutant)

# Irrespective of significance, how many genes have:
# same proximal PAS between two datasets (pas_number == 1)
# same number of PAS between two datasets (delta psis comparable)
# same PAS betweeen two datasets

pas1 <- comb %>%
  filter(Pas_Number == 1) %>%
  group_by(Gene_Name) %>%
  mutate(both_mut = n_distinct(mutant) == 2,
         common_pas = n_distinct(APA_ID) == 1) %>%
  ungroup() 

pas1_padj <- comb_padj %>%
  filter(Pas_Number == 1) %>%
  group_by(Gene_Name) %>%
  mutate(both_mut = n_distinct(mutant) == 2,
         common_pas = n_distinct(APA_ID) == 1) %>%
  ungroup() 
  
pas1_padj %>%
  distinct(Gene_Name, both_mut, common_pas) %>%
  count(both_mut, common_pas)


comb_padj %>%
  group_by(mutant, Gene_Name) %>%
  slice_min(Pas_Number, n = 1) %>%
  group_by(Gene_Name) %>%
  mutate(both_mut = n_distinct(mutant) == 2,
         common_pas = n_distinct(APA_ID) == 1) %>%
  ungroup() %>%
  arrange(desc(both_mut)) %>%
  distinct(Gene_Name, both_mut, common_pas) %>%
  filter(both_mut & !common_pas)
  
# count(both_mut, common_pas)
# A tibble: 3 × 3
# both_mut common_pas     n
# <lgl>    <lgl>      <int>
#   1 FALSE    TRUE         377
# 2 TRUE     FALSE         15
# 3 TRUE     TRUE           3

# 15 genes where significant shared PAS usage, but the most proximal sig site is different between the two datasets


# WHat about number of PAS
comb_npas <- comb %>%
  group_by(mutant, Gene_Name) %>%
  summarise(n_polya = n_distinct(APA_ID)) %>%
  group_by(Gene_Name) %>%
  summarise(both_mut = n_distinct(mutant) == 2,
         common_n_pas = n_distinct(n_polya) == 1) 

count(comb_npas, both_mut, common_n_pas) %>% mutate(frac = n / sum(n))
# A tibble: 3 × 4
# both_mut common_n_pas     n  frac
# <lgl>    <lgl>        <int> <dbl>
#   1 FALSE    TRUE           829 0.166
# 2 TRUE     FALSE         1093 0.219
# 3 TRUE     TRUE          3065 0.615

# around 80 % genes commonly evaluated
# around 3/4 common genes have shared number of PAS

# For shared num PAS, what fraction have exactly the same PAS?
same_npas_gn <- comb_npas %>% filter(both_mut & common_n_pas) %>%
  pull(Gene_Name)

comb_same_npas_pas <- comb %>%
  filter(Gene_Name %in% same_npas_gn) %>%
  arrange(Gene_Name, mutant, Pas_Number) %>%
  # for each gene, does number of total unique IDs == number of PAS?
  group_by(Gene_Name) %>%
  summarise(same_pas = n_distinct(APA_ID) == max(Pas_Number))

# of the commonly evaluated, 94 % have exactly the same PAS evaluated
comb_same_npas_pas %>%
  count(same_pas) %>%
  mutate(frac = n / sum(n))

# fraction as a proportion of total PAS
length(intersect(comb_same_npas_pas %>% filter(same_pas) %>% pull(Gene_Name),
          unique(comb$Gene_Name))) / length(unique(comb$Gene_Name))

# ~ 58 % of genes evaluated across two datasets have exactly the same number of PAS 

# Intersection between genes with significant APA in both datasets
# followed by scatter plot
# optional subsetting for same PAS, same N pas etc.


# any PAS in gene has significant APA
anysig_gl <- list("F210I" = unique(dexseq_f210i_padj$Gene_Name), "M323K" =  unique(dexseq_m323k_padj$Gene_Name))
ggVennDiagram(anysig_gl, label = "count") +
  labs(title = "F210I vs M323K MEF",
       subtitle = "differential APA targets (any PAS = padj < 0.05)")

ggsave(file.path(outdir, "2024-03-04_mef_venn_padj_005.png"),height = 10, width = 10,dpi = "retina")


anysig_venn <- Venn(anysig_gl)
anysig_shared <- overlap(anysig_venn)
anysig_f210i <- discern(anysig_venn, "F210I")
anysig_m323k <- discern(anysig_venn, "M323K")


# plot scatter of psi values between shared and unique PAS
gene2nm <- bind_rows(distinct(dexseq_f210i, Gene, Gene_Name), distinct(dexseq_m323k, Gene, Gene_Name)) %>% distinct(.keep_all = T)

psi_anysig_df <-  bind_rows("F210I" = psi_f210i, "M323K" = psi_m323k, .id = "mutant") %>%
  left_join(gene2nm, by = "Gene") %>%
  filter(Gene_Name %in% c(anysig_shared, anysig_m323k, anysig_f210i)) %>%
  mutate(overlap_group = case_when(Gene_Name %in% anysig_shared ~ "Both",
                                   Gene_Name %in% anysig_f210i ~ "F210I",
                                   Gene_Name %in% anysig_m323k ~ "M323K",
                                   T ~ ""))

plot_psi_anysig_df <- psi_anysig_df %>%
  select(Gene_Name, overlap_group, mutant, deltaPSImean.HOM_WT, deltaPSImedian.HOM_WT) %>%
  pivot_wider(id_cols = all_of(c("Gene_Name", "overlap_group")),
              names_from = mutant,
              values_from = c(deltaPSImean.HOM_WT, deltaPSImedian.HOM_WT), names_sep = ".") %>%
  mutate(plot_alpha = if_else(overlap_group == "Both", 1, 0.9)) 

mef_comparison_scatter_padj_005 <- ggplot() +
  geom_point(data = filter(plot_psi_anysig_df, overlap_group == "M323K"), aes(x = deltaPSImean.HOM_WT.M323K, y = deltaPSImean.HOM_WT.F210I, color = "M323K"), alpha = 0.6) +
  geom_point(data = filter(plot_psi_anysig_df, overlap_group == "F210I"), aes(x = deltaPSImean.HOM_WT.M323K, y = deltaPSImean.HOM_WT.F210I, color = "F210I"), alpha = 0.6) +
  geom_point(data = filter(plot_psi_anysig_df, overlap_group == "Both"), aes(x = deltaPSImean.HOM_WT.M323K, y = deltaPSImean.HOM_WT.F210I, color = "Both"), alpha = 1, size = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_color_manual(name = "Overlap Group", values = c("M323K" = "#ef8a62", "F210I" = "#67a9cf", "Both" = "black")) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_x_continuous(limits = c(-1,1)) +
  labs(title = "M323K vs F210I MEF differential APA",
       subtitle = "Any polyA site padj < 0.05. Overlap at gene-level") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
  


# assign to quadrants based on dirn of change in M323K and F210I
psi_anysig_dirn_groups <- plot_psi_anysig_df %>%
  mutate(dirn_group = case_when(sign(deltaPSImean.HOM_WT.M323K) == 1 & sign(deltaPSImean.HOM_WT.F210I) == 1 ~ "M323Klonger_F210Ilonger",
                                sign(deltaPSImean.HOM_WT.M323K) == 1 & sign(deltaPSImean.HOM_WT.F210I) == -1 ~ "M323Klonger_F210Ishorter",
                                sign(deltaPSImean.HOM_WT.M323K) == -1 & sign(deltaPSImean.HOM_WT.F210I) == 1 ~ "M323Kshorter_F210Ilonger",
                                sign(deltaPSImean.HOM_WT.M323K) == -1 & sign(deltaPSImean.HOM_WT.F210I) == -1 ~ "M323Kshorter_F210Ishorter",
                                )
         ) %>%
  distinct(Gene_Name, overlap_group, dirn_group)

# Get a df of fraction of genes in each quadrant split by overlap group
plot_psi_dirn_anysig_df_split <- count(psi_anysig_dirn_groups, overlap_group, dirn_group) %>%
  drop_na(dirn_group) %>%
  group_by(overlap_group) %>%
  mutate(frac = n / sum(n),
         dirn_group = fct_reorder(dirn_group, n)) %>%
  ungroup() 

# repeat but global (i.e. ignore overlap group)
plot_psi_dirn_anysig_df_all <- count(psi_anysig_dirn_groups, dirn_group) %>%
  drop_na(dirn_group) %>%
  mutate(frac = n / sum(n),
         overlap_group = "NA") %>%
  mutate(dirn_group = fct_reorder(dirn_group, n)) 

quadrant_frac_bar_padj_005 <- bind_rows(plot_psi_dirn_anysig_df_split, plot_psi_dirn_anysig_df_all) %>%
  mutate(dirn_group = fct_reorder(dirn_group, frac)) %>%
  ggplot(aes(y = dirn_group, x = frac, fill = overlap_group, label = n)) + 
  geom_col(position = "dodge") +
  scale_x_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.2)) +
  scale_fill_manual(values = c("M323K" = "#ef8a62", "F210I" = "#67a9cf", "Both" = "black", "NA" = "purple")) +
  labs(title = "M323K vs F210I MEF Quadrant Fractions",
       subtitle = "NA Overlap Group = altogether (ignore overlap group)",
       x = "Fraction of APA genes",
       y = "",
       fill = "Overlap Group") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
  
quadrant_frac_bar_padj_005

ggsave(filename = file.path(outdir, "2024-03-04_mef_mut_scatter_deltapsi_padj_005.png"),
       plot = mef_comparison_scatter_padj_005,
       device = "png",units = "in", dpi = "retina", height = 8, width = 8
       )  

ggsave(filename = file.path(outdir, "2024-03-04_mef_mut_quadrants_bar_deltapsi_padj_005.png"),
       plot = quadrant_frac_bar_padj_005,
       device = "png",units = "in", dpi = "retina", height = 8, width = 8
) 


### Dosage pattern plots

# Plan - for each mutant:
# - try 3 different selection procedures + PAU
# - try psi values


# Get psi dfs in desired format for plots
dexseq_m323k_padj_psi <- select(dexseq_m323k_padj, Gene, Chr) %>%
  distinct() %>%
  left_join(psi_m323k, by = c("Gene", "Chr"))

dexseq_f210i_padj_psi <- select(dexseq_f210i_padj, Gene, Chr) %>%
  distinct() %>%
  left_join(psi_f210i, by = c("Gene", "Chr"))


# get plot ready dfs for psi values
dosage_plot_dfs_psi <- list("M323K psi" = dexseq_m323k_padj_psi,
     "F210I psi" = dexseq_f210i_padj_psi) %>%
  map(~ prep_dosage_plot_df(.x, col_prefix = "PSImean.", names_col = "condition",
                           values_col = "meanPSI", pattern_col = "meanPSI.dosage_pattern") %>%
        mutate(gene_chr = paste(Gene, Chr, sep = "_"))
      )

#base dosage plots for each MUT using psi values
dosage_plots_psi <- map2(.x = dosage_plot_dfs_psi,
                         .y = names(dosage_plot_dfs_psi),
                         ~ base_dosage_plot(
                           data = .x,
                           x_var = "condition",
                           y_var = "meanPSI",
                           group_var = "gene_chr",
                           facet_var = "meanPSI.dosage_pattern",
                           alpha_val = 0.25,
                           base_font_size = 14
                         ) +
                           labs(title = .y,
                                subtitle = "any PAS padj < 0.05",
                                x = "Condition",
                                y = "mean LABRAT psi")
                         )

dosage_plots_psi


# PAU values, want to try 3 selection procedures for each MUT
# Suggest = select + prepare df separately
# combine + make plots

# try 3 selection procedures and generate plot-ready dfs for each selection
dosage_plot_dfs_pau_m323k <- list("M323K Halleger" = select_rep_pas_halleger,
     "M323K MostProximal" = select_rep_pas_mostproximal,
     "M323K Proximal" = select_rep_pas_proximal) %>%
  map(~ .x(dexseq_m323k_padj) %>%
        prep_dosage_plot_df(
          col_prefix = "mean_PAU_", names_col = "condition", values_col = "meanPAU", pattern_col = "meanPAU.dosage_pattern")
      )

dosage_plot_dfs_pau_f210i <- list("F210I Halleger" = select_rep_pas_halleger,
     "F210I MostProximal" = select_rep_pas_mostproximal,
     "F210I Proximal" = select_rep_pas_proximal) %>%
  map(~ .x(dexseq_f210i_padj) %>%
        prep_dosage_plot_df(
          col_prefix = "mean_PAU_", names_col = "condition", values_col = "meanPAU", pattern_col = "meanPAU.dosage_pattern")
      )

# combine the two together
dosage_plot_dfs_pau <- c(dosage_plot_dfs_pau_m323k, dosage_plot_dfs_pau_f210i)

dosage_plots_pau <- map2(.x = dosage_plot_dfs_pau,
                         .y = names(dosage_plot_dfs_pau),
                         ~ base_dosage_plot(
                           data = .x,
                           x_var = "condition",
                           y_var = "meanPAU",
                           group_var = "APA_ID",
                           facet_var = "meanPAU.dosage_pattern",
                           alpha_val = 0.25,
                           base_font_size = 14
                         ) +
                           labs(title = .y,
                                subtitle = "PAS padj < 0.05",
                                x = "Condition",
                                y = "mean PAS usage (%)")
)


dosage_plots_pau

# Get counts/percentages of the number of events in each dosage category
dosage_group_counts_psi <- map(dosage_plot_dfs_psi,
    ~ .x %>%
      distinct(gene_chr, .keep_all = T) %>%
  count(meanPSI.dosage_pattern, sort = T) %>%
  mutate(perc = (n / sum(n)) * 100 )
  )

dosage_group_counts_pau <- map(dosage_plot_dfs_pau,
    ~ .x %>%
      distinct(APA_ID, .keep_all = T) %>%
      count(meanPAU.dosage_pattern, sort = T) %>%
      mutate(perc = (n / sum(n)) * 100 )
    )

dosage_group_simplified_counts_psi <- map(dosage_plot_dfs_psi,
                               ~ .x %>%
                                 distinct(gene_chr, .keep_all = T) %>%
                                 count(meanPSI.dosage_pattern_simplified, sort = T) %>%
                                 mutate(perc = (n / sum(n)) * 100 )
)

dosage_group_simplified_counts_pau <- map(dosage_plot_dfs_pau,
                               ~ .x %>%
                                 distinct(APA_ID, .keep_all = T) %>%
                                 count(meanPAU.dosage_pattern_simplified, sort = T) %>%
                                 mutate(perc = (n / sum(n)) * 100 )
)

bind_rows(dosage_group_counts_psi, .id = "dataset")
bind_rows(dosage_group_counts_pau, .id = "dataset")

# Make dosage pattern frequency bar plots
dosage_count_plots_pau <- map2(.x = dosage_group_counts_pau,
     .y = names(dosage_group_counts_pau),
     ~ .x %>% 
       ggplot(aes(x = meanPAU.dosage_pattern, y = perc, label = n)) +
       geom_col() +
       geom_text(nudge_y = 5) +
       scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
       theme_bw(base_size = 14) +
       labs(title = .y,
            subtitle = "label = number of events",
            x = "PAS Usage Pattern",
            y = "Event Frequency (%)")
     )


dosage_count_plots_psi <- map2(.x = dosage_group_counts_psi,
     .y = names(dosage_group_counts_psi),
     ~ .x %>% 
       ggplot(aes(x = meanPSI.dosage_pattern, y = perc, label = n)) +
       geom_col() +
       geom_text(nudge_y = 5) +
       scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
       theme_bw(base_size = 14) +
       labs(title = .y,
            subtitle = "label = number of events",
            x = "PAS Usage Pattern",
            y = "Event Frequency (%)")
     )


# Output dosage plots

# line plots
walk2(.x = dosage_plots_pau,
      .y = names(dosage_plots_pau),
      ~ ggsave(filename = file.path(outdir,
                                    paste("2024-05-29_mef_dosage_direction_facet",
                                          "pau",
                                          str_replace_all(.y, " ", "_"),
                                          "png",
                                          sep = ".")
      ),
      plot = .x,
      device = "png",
      units = "in",
      dpi = "retina",
      height = 8,
      width = 8
      )
)

walk2(.x = dosage_plots_psi,
      .y = names(dosage_plots_psi),
      ~ ggsave(filename = file.path(outdir,
                                    paste("2024-05-29_mef_dosage_direction_facet",
                                          "psi",
                                          str_replace_all(.y, " ", "_"),
                                          "png",
                                          sep = ".")
      ),
      plot = .x,
      device = "png",
      units = "in",
      dpi = "retina",
      height = 8,
      width = 8
      )
)

# category counts bar plots
walk2(.x = dosage_count_plots_pau,
      .y = names(dosage_count_plots_pau),
      ~ ggsave(filename = file.path(outdir,
                                    paste("2024-05-29_mef_dosage_category_counts",
                                          "pau",
                                          str_replace_all(.y, " ", "_"),
                                          "png",
                                          sep = ".")
      ),
      plot = .x,
      device = "png",
      units = "in",
      dpi = "retina",
      height = 8,
      width = 8
      )
)

walk2(.x = dosage_count_plots_psi,
      .y = names(dosage_count_plots_psi),
      ~ ggsave(filename = file.path(outdir,
                                    paste("2024-05-29_mef_dosage_category_counts",
                                          "psi",
                                          str_replace_all(.y, " ", "_"),
                                          "png",
                                          sep = ".")
      ),
      plot = .x,
      device = "png",
      units = "in",
      dpi = "retina",
      height = 8,
      width = 8
      )
)