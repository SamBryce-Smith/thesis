library(tidyverse)


#' Summarise expression across all  PAS within a gene/group
summarize_expression <- function(data, grouping_cols = c("Gene", "Chr"), exprn_col = "exonBaseMean") {
  
  data %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(total_exprn = sum(!!sym(exprn_col), na.rm = TRUE)) %>%
    ungroup()
}


#' split input processed DEXSeq df into genes with 2 analysed PAS or >2 analysed PAS
split_by_npas <- function(data, grouping_cols = c("Gene", "Chr"), pas_count_col = "Total_Pas", return_intermediate = FALSE) {
  
  npas_data <- data %>%
    distinct(across(all_of(c(grouping_cols, pas_count_col))))
  
  summary_npas <- npas_data %>%
    summarise(
      n_2pas = sum(!!sym(pas_count_col) == 2),
      frac_2pas = n_2pas / n()
    )
  
  npas2_data <- npas_data %>%
    filter(!!sym(pas_count_col) == 2) %>%
    select(-!!sym(pas_count_col))
  
  #
  npas2_out <- left_join(npas2_data, data, by = grouping_cols)
  #
  npasmore_out <- anti_join(data, npas2_data, by = grouping_cols)
  
  if (return_intermediate) {
    return(list(
    npas2 = npas2_out,
    npasmore = npasmore_out,
    summary_npas = summary_npas))
    
  } else {
    return(list(
      npas2 = npas2_out,
      npasmore = npasmore_out
      ))
    
  }
}


#' Split an input processed DEXSeq df into 2, where groups have no significant events or at least 1
split_by_sig <- function(data, grouping_cols = c("Gene", "Chr"), padj_col = "padj", padj_threshold = 0.05) {
  
  ns <- data %>%
    group_by(across(all_of(grouping_cols))) %>%
    filter(all(!!sym(padj_col) > padj_threshold)) %>%
    ungroup()
  
  sig <- data %>%
    group_by(across(all_of(grouping_cols))) %>%
    filter(any(!!sym(padj_col) < padj_threshold)) %>%
    ungroup()
  
  list(
    ns = ns,
    sig = sig
  )
}

#' split dexseq df into genes/groups with 1, 2 or 3+ significant PAS
split_by_nsig <- function(data, grouping_cols = c("Gene", "Chr"), padj_col = "padj", padj_threshold = 0.05, return_intermediate = FALSE) {
  # Calculate the number of significant padj values per group
  sig_count_data <- data %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(nsig = sum(!!sym(padj_col) < padj_threshold, na.rm = TRUE)) %>%
    ungroup()
  
  # Count occurrences of nsig values
  nsig_counts <- sig_count_data %>%
    count(nsig)
  
  # Split data into nsig == 1, nsig == 2, and nsig > 2
  nsig_1 <- sig_count_data %>%
    filter(nsig == 1) %>%
    select(-nsig) %>%
    left_join(data, by = grouping_cols)
  
  nsig_2 <- sig_count_data %>%
    filter(nsig == 2) %>%
    select(-nsig) %>%
    left_join(data, by = grouping_cols)
  
  nsig_3plus <- sig_count_data %>%
    filter(nsig > 2) %>%
    select(-nsig) %>%
    left_join(data, by = grouping_cols)
  
  # Combine results
  result <- list(
    nsig_1 = nsig_1,
    nsig_2 = nsig_2,
    nsig_3plus = nsig_3plus
  )
  
  if (return_intermediate) {
    result$nsig_counts <- nsig_counts
  }
  
  result
}

#' Select representative PAs for genes with all PAS padj > 0.05/threshold (pick two highest expressed PAS)
rep_pas_0_sig <- function(data, grouping_cols = c("Gene", "Chr"), exprn_col = "exonBaseMean") {
  result <- data %>%
    group_by(across(all_of(grouping_cols))) %>%
    slice_max(!!sym(exprn_col), n = 2) %>%
    ungroup()
  
  result
}


#' Select representative PAS for genes with at least 1 sig padj < 0.05 (picking highest expressed ns PAS as alternative)
rep_pas_1_sig <- function(data, grouping_cols = c("Gene", "Chr"), padj_col = "padj", log2fc_col = "log2fold_HOM_WT", exprn_col = "exonBaseMean", padj_threshold = 0.05) {
  
  # Step 1: Filter, mutate and join data
  sig_ns <- data %>%
    filter(!!sym(padj_col) < padj_threshold) %>%
    mutate(sign_sig = sign(!!sym(log2fc_col))) %>%
    distinct(across(all_of(c(grouping_cols, "sign_sig")))) %>%
    left_join(filter(data, !!sym(padj_col) > padj_threshold), ., by = grouping_cols)
  
  # Step 2: Pick next highest expressed PAS
  rephigh_ns <- sig_ns %>%
    group_by(across(all_of(grouping_cols))) %>%
    slice_max(!!sym(exprn_col), n = 1) %>%
    ungroup()
  
  # Step 3: Combine with significant data
  comb <- bind_rows(rephigh_ns, filter(data, !!sym(padj_col) < padj_threshold)) %>%
    select(-sign_sig)
  
  comb
}


#' Calculate proportion of total expression attributed to representative PAS
rep_pas_expression <- function(combined_data, total_exprn_data, grouping_cols = c("Gene", "Chr"), exprn_col = "exonBaseMean", padj_col = "padj") {
#
 
  # Summarize total expression for each group
  rep_total_exprn_data <- combined_data %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(rep_total_exprn = sum(!!sym(exprn_col), na.rm = TRUE)) %>%
    ungroup()
  
  # Join with total expression data and calculate proportion
  result <- rep_total_exprn_data %>%
    left_join(total_exprn_data, by = grouping_cols) %>%
    mutate(propn_total_exprn = rep_total_exprn / total_exprn)
  
  result
}

#' Select 2 rep PAS for a gene, filtering for genes where the two sig PAS change in opposing directions
rep_pas_2_sig <- function(data, grouping_cols = c("Gene", "Chr"), padj_col = "padj", log2fc_col = "log2fold_HOM_WT", padj_threshold = 0.05) {
  
  # find genes/groups with opposing dirns in usage for 2 sig PAS
  result <- data %>%
    filter(!!sym(padj_col) < padj_threshold) %>%
    group_by(across(all_of(grouping_cols))) %>%
    # +1 for +ve, -1 for -ve 
    summarise(sum_dirn = sum(sign(!!sym(log2fc_col)))) %>%
    ungroup() %>%
    # filter for genes/groups with opposing signs for 2 PAS
    filter(sum_dirn == 0) %>%
    select(-sum_dirn)
  
  # join back in PAS + dexseq cols for those passing opposing signs check
  result %>%
    left_join(filter(data, !!sym(padj_col) < padj_threshold), by = grouping_cols)
}


#' Select representative PAS for genes with > 2 significant PAS (select two PAS with largest expression changing in opposing directions)
rep_pas_3plus_sig <- function(data, grouping_cols = c("Gene", "Chr"), padj_col = "padj", log2fc_col = "log2fold_HOM_WT", exprn_col = "exonBaseMean", padj_threshold = 0.05) {
  result <- data %>%
    filter(!!sym(padj_col) < padj_threshold) %>%
    mutate(delta_sign = sign(!!sym(log2fc_col))) %>%
    group_by(across(all_of(grouping_cols)), delta_sign) %>%
    slice_max(!!sym(exprn_col), n = 1) %>%
    ungroup() %>%
    select(-delta_sign)
  
  result
}

#' Select two representative PAS from a DEXSeq APA results table
rep_pas_wrapper <- function(data, grouping_cols = c("Gene", "Chr"), padj_threshold = 0.05, log2fc_col = "log2fold_HOM_WT", exprn_col = "exonBaseMean", pas_count_col = "Total_Pas",  padj_col = "padj") {
  
  # total expression of all analysed PAS
  total_exprn <- summarize_expression(data,grouping_cols,exprn_col)
  
  # split input into dfs of 2 PAS genes (no need for further processing) and >=2 PAS genes
  data_npas_split <- split_by_npas(data, grouping_cols, pas_count_col, return_intermediate = T)
  
  # split input df of > 2 PAS genes into sig and not sig (sig = any PAS passes threshold)
  sig_pas_split <- split_by_sig(data_npas_split$npasmore, grouping_cols, padj_col, padj_threshold)
  
  # Select representative PAS for ns genes
  rep_0_sig <- rep_pas_0_sig(sig_pas_split$ns, grouping_cols, exprn_col)
  
  # split >2 PAS sig genes by number of significant PAS
  npasmore_sig_split <- split_by_nsig(sig_pas_split$sig, grouping_cols, padj_col, padj_threshold)
  
  # Select rep PAS for genes with 1 sig PAS (sig PAS + highest expressed ns PAS)
  rep_1_sig <- rep_pas_1_sig(npasmore_sig_split$nsig_1, grouping_cols, padj_col, log2fc_col, exprn_col, padj_threshold)
  
  # genes with 2 sig PAS (both sig PAS, but must have different directions of change in usage)
  rep_2_sig <- rep_pas_2_sig(npasmore_sig_split$nsig_2, grouping_cols, padj_col, log2fc_col, padj_threshold)
  
  # TODO: check if rep_2sig drops genes compared to npasmore_sig_split$nsig_2
  
  # genes with > 2 sig PAS (2 sig PAS, largest expression in opposite direction)
  rep_3plus_sig <- rep_pas_3plus_sig(npasmore_sig_split$nsig_3plus, grouping_cols, padj_col, log2fc_col, exprn_col, padj_threshold)
  
  # Have selected rep PAS for all splits - combine
  rep_comb <- bind_rows(data_npas_split$npas2,
                        rep_0_sig,
                        rep_1_sig,
                        rep_2_sig,
                        rep_3plus_sig) %>%
    arrange(across(all_of(grouping_cols)))
  
  # Summarise expression of representative PAS
  rep_total_exprn <- summarize_expression(rep_comb, grouping_cols, exprn_col)
  
  # compare summarised expression of representative to total PAS expression for gene
  rep_total_exprn <- rep_total_exprn %>%
    left_join(total_exprn, by = grouping_cols, suffix = c(".rep", ".all")) %>%
    mutate(rep_fracn_total = total_exprn.rep / total_exprn.all)
  
  
  return(list(rep_df = rep_comb, exprn_df = rep_total_exprn))
  
}


##

# processed DEXSeq results tables produced by QAPA_snakemake pipeline
dexseq_f210i <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")
dexseq_m323k <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")
base_outdir <- "processed/2024-02-29"
comparison_outdir <- "processed/2024-02-29/analysis/mef_comparisons"
if (!dir.exists(comparison_outdir)) {dir.create(comparison_outdir, recursive = T)}

genes_init_f210i <- distinct(dexseq_f210i, Gene, Chr) 
genes_init_m323k <- distinct(dexseq_m323k, Gene, Chr)

rep_pas_m323k <- rep_pas_wrapper(dexseq_m323k)
rep_pas_f210i <- rep_pas_wrapper(dexseq_f210i)

genes_rep_f210i <- distinct(rep_pas_f210i$rep_df, Gene, Chr)
genes_rep_m323k <- distinct(rep_pas_m323k$rep_df, Gene, Chr)


# Check that all genes are retained after selecting rep PAS
bind_rows(init = genes_init_f210i, rep = genes_rep_f210i, .id = "source") %>%
  count(Gene, Chr) %>%
  filter(n != 2)

bind_rows(init = genes_init_m323k, rep = genes_rep_m323k, .id = "source") %>%
  count(Gene, Chr) %>%
  filter(n != 2)

# both 0-length, so all genes retained

rep_pas <- list("M323K_MEFs" = rep_pas_m323k,
                "F210I_MEFs" = rep_pas_f210i)
  
walk2(.x = rep_pas,
      .y = names(rep_pas),
      ~ write_tsv(.x$rep_df,
                  file = file.path(base_outdir, .y, "dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv"),
                  col_names = T
                  )
      )

walk2(.x = rep_pas,
      .y = names(rep_pas),
      ~ write_tsv(.x$exprn_df,
                  file = file.path(base_outdir, .y, "dexseq_apa.HOMvsWT.rep_pas_total_expression.tsv"),
                  col_names = T
      )
)

# Comparison of M323K and F210I for representative PAS
int_rep_m323k <- rep_pas_m323k$rep_df %>%
group_by(Gene, Chr) %>%
  mutate(sig = any(padj < 0.05)) %>%
  ungroup() %>%
  distinct(Gene, Chr, APA_ID, sig, padj, log2fold_HOM_WT)

int_rep_f210i <- rep_pas_f210i$rep_df %>%
  group_by(Gene, Chr) %>%
  mutate(sig = any(padj < 0.05)) %>%
  ungroup() %>%
  distinct(Gene, Chr, APA_ID, sig, padj, log2fold_HOM_WT)

# add all analysed combos of gene & APA_ID together (allowing missing in one group)
rep_joined <- full_join(int_rep_m323k,
                        int_rep_f210i,
                        by = c("Gene", "Chr", "APA_ID"),
                        suffix = c(".m323k", ".f210i")
) %>%
  arrange(Gene, Chr)

# Do genes have same PAS / number of PAS analysed between the two conditions?
rep_summary_joined <- rep_joined %>%
  mutate(pas_in_both = !is.na(sig.m323k) & !is.na(sig.f210i)) %>%
  group_by(Gene, Chr) %>%
  summarise(n_pas.m323k = sum(!is.na(sig.m323k)),
            n_pas.f210i = sum(!is.na(sig.f210i)),
            same_n_pas = n_pas.m323k == n_pas.f210i,
            same_pas = all(pas_in_both),
            sig.m323k = replace_na(any(sig.m323k), F),
            sig.f210i = replace_na(any(sig.f210i), F),
            sig.either = sig.m323k | sig.f210i
  ) %>%
  ungroup()

# what fraction of genes have same expressed PAS (and are directly comparable)?
rep_counts_summary_joined <- rep_summary_joined %>%
  count(sig.either, same_pas, same_n_pas, sort = T) %>%
  group_by(sig.either) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# write to file
write_tsv(rep_summary_joined,
          file.path(comparison_outdir, "f210i_vs_m323k.representative_pas_consistency.summary.tsv"),
          col_names = T
          )

write_tsv(rep_counts_summary_joined,
          file.path(comparison_outdir, "f210i_vs_m323k.representative_pas_consistency.counts.tsv"),
          col_names = T
          )


# what about for psi analysis?
# same PAS = exactly the same PAS between two comparisons
# same number of PAS = same number of PAS (delta psis are more comparable)

# assign significance label + collapse metadata to analysed PAS
int_psi_dexseq_m323k <- dexseq_m323k %>%
  group_by(Gene, Chr) %>%
  mutate(sig = any(padj < 0.05)) %>%
  ungroup() %>%
  distinct(Gene, Chr, APA_ID, sig, log2fold_HOM_WT)

int_psi_dexseq_f210i <- dexseq_f210i %>%
  group_by(Gene, Chr) %>%
  mutate(sig = any(padj < 0.05)) %>%
  ungroup() %>%
  distinct(Gene, Chr, APA_ID, sig, log2fold_HOM_WT)

# Join all APA_IDs between datasets 
psi_dexseq_joined <- full_join(int_psi_dexseq_m323k,
                               int_psi_dexseq_f210i,
                               by = c("Gene", "Chr", "APA_ID"),
                               suffix = c(".m323k", ".f210i")
)

# Do genes have same PAS / number of PAS analysed between the two conditions?
psi_summary_joined <- psi_dexseq_joined %>%
  mutate(pas_in_both = !is.na(sig.m323k) & !is.na(sig.f210i)) %>%
  group_by(Gene, Chr) %>%
  summarise(n_pas.m323k = sum(!is.na(sig.m323k)),
            n_pas.f210i = sum(!is.na(sig.f210i)),
            same_n_pas = n_pas.m323k == n_pas.f210i,
            same_pas = all(pas_in_both),
            sig.m323k = replace_na(any(sig.m323k), F),
            sig.f210i = replace_na(any(sig.f210i), F),
            sig.either = sig.m323k | sig.f210i
  ) %>%
  ungroup()

# what fraction of genes have same expressed PAS (and are directly comparable)?
psi_counts_summary_joined <- psi_summary_joined %>%
  count(sig.either, same_pas, same_n_pas, sort = T) %>%
  group_by(sig.either) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

write_tsv(psi_summary_joined,
          file.path(comparison_outdir, "f210i_vs_m323k.psi_consistency.summary.tsv")
          )

write_tsv(psi_counts_summary_joined,
          file.path(comparison_outdir, "f210i_vs_m323k.psi_consistency.counts.tsv")
          )

# How many genes are consisent between the two datasets that would be able to compare?
# Which one is more strict/sensitive?

# here, need to directly compare the same PAS
rep_consistent_counts <- rep_counts_summary_joined %>%
  filter(same_pas) %>%
  group_by(sig.either) %>%
  summarise(total = sum(n))

# with psi approach, can analyse if same PAS / same number of PAS
psi_consistent_counts_relaxed <- psi_counts_summary_joined %>%
  filter(same_pas | same_n_pas) %>%
  group_by(sig.either) %>%
  summarise(total = sum(n))

# if more cautious (only same PAS), how many genes can you analyse?
psi_consistent_counts_strict <- psi_counts_summary_joined %>%
  filter(same_pas) %>%
  group_by(sig.either) %>%
  summarise(total = sum(n))

consistent_counts_comb <- bind_rows("representative_pas" = rep_consistent_counts,
          "psi_strict" = psi_consistent_counts_strict,
          "psi_relaxed" = psi_consistent_counts_relaxed,
          .id = "comparison_approach")

write_tsv(consistent_counts_comb, file.path(comparison_outdir, "f210i_vs_m323k.representative_vs_psi.summary_counts.tsv"))

# psi approach allows us to look at more PAS (whilst also guaranteeing you represent the entire PAS spectrum)
# TODO: rep pas - are matching PAS more likely to be proximal/distal? Is there an approach where pick 1 match and check consistency?

