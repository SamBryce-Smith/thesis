library(tidyverse)
library(ggprism)
library(rstatix)
library(eulerr)

#' Define APA Targets
#'
#' This function filters for gene-level significance, annotates proximal/distal PAS,
#' and applies filters for reliable PAS targets based on signficance and consistency of change in usage with genotype. 
#' Filtering is performed in the following steps:
#' 1. Gene has any PAS passing padj_threshold (reported in 'gene')
#' Separated for proximal and distal:
#' 2. PAS passes padj threshold AND
#' 2A. Change in usage consistent with genotype OR
#' 2B. Inconsistent with genotype, but change in usage in HETs vs WT is small OR
#' 2C. Inconsistent with genotype, change in usage in HETs vs WT is not small but change in HOMs vs HETs is small
#' 
#' 2B = indicative of noise in means + regulation only in HOMs, 2C = indicative of 'maximal regulation' in HETs and noise in means in HOMs vs HETs
#' 
#'
#' @param dexseq_data A data frame containing DEXSeq results. Minimal Required columns: "Gene", "Chr", "padj", "Pas_Number", "APA_ID", "Gene_Name", "meanPAU.dosage_pattern_simplified", "delta_HET_WT", "delta_HOM_HET".
#' @param padj_threshold A numeric value specifying the threshold for adjusted p-value. Default is 0.05.
#' @param wt_inconsistent_threshold A numeric value specifying the threshold for maximum inconsistency in usage between HETs vs WT. Default is 5.
#' @param het_inconsistent_threshold A numeric value specifying the threshold for maximum inconsistency in usage between HOMs vs HETs. Default is 5.
#'
#' @return A list containing filtered data frames: "all" for combined proximal/distal targets, "proximal" and "distal" for specific site types, and "gene" for the gene-level (any passing padj_threshold) data frame.
#'
define_apa_targets <- function(dexseq_data, padj_threshold = 0.05, wt_inconsistent_threshold = 5, het_inconsistent_threshold = 5) {
  
  # Required columns
  required_columns <- c("Gene", "Chr", "padj", "Pas_Number", "APA_ID", "Gene_Name", 
                        "meanPAU.dosage_pattern_simplified", "delta_HET_WT", "delta_HOM_HET")
  
  # Check if all required columns are present
  missing_columns <- setdiff(required_columns, colnames(dexseq_data))
  
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing in dexseq_data:", paste(missing_columns, collapse = ", ")))
  }
  
  
  # filter for gene-level significance (any PAS below threshold)
  dexseq_sig <- dexseq_data %>%
    group_by(Gene, Chr) %>%
    filter(any(padj < padj_threshold)) %>%
    # annotate proximal/distal
    mutate(site_type = if_else(Pas_Number == min(Pas_Number), "proximal", "distal"),
           site_type = factor(site_type, levels = c("proximal", "distal"))) %>%
    ungroup()
  
  # filter for proximal/distal
  grpd <- group_by(dexseq_sig, site_type)
  
  # apply filters for reliable PAS targets (separately for proximal/ distal)
  dexseq_sig_filt_split <- group_split(grpd) %>%
    set_names(pull(group_keys(grpd))) %>%
    map(~ filter(.x,
                 padj < padj_threshold,
                 (meanPAU.dosage_pattern_simplified == "consistent") | 
                   (meanPAU.dosage_pattern_simplified == "not_consistent" & abs(delta_HET_WT) < wt_inconsistent_threshold) | 
                   (meanPAU.dosage_pattern_simplified == "not_consistent" & abs(delta_HET_WT) >= wt_inconsistent_threshold & abs(delta_HOM_HET) < het_inconsistent_threshold)
    )
    )
  
  # combine across proximal vs distal (i.e. 'gene-level' df)
  dexseq_sig_filt_all <- c(list("all" = bind_rows(dexseq_sig_filt_split)), dexseq_sig_filt_split, list("gene" = dexseq_sig))
  
  return(dexseq_sig_filt_all)
}



#' Compute Overlap Statistics for Two Vectors
#'
#' This function computes overlap statistics for two input vectors and returns a 1-row dataframe.
#'
#' @param vec1 The first vector.
#' @param vec2 The second vector.
#' @param comparison_value A key to identify the comparison (e.g., "APA_ID") (fills the comparison_value column)
#'
#' @return A 1-row dataframe containing the overlap statistics.
#'
#' @export
compute_stats <- function(vec1, vec2) {
  unique1 <- length(unique(vec1))
  unique2 <- length(unique(vec2))
  common <- length(intersect(vec1, vec2))
  distinct1 <- length(setdiff(vec1, vec2))
  distinct2 <- length(setdiff(vec2, vec1))
  total <- length(union(vec1, vec2))
  # intersection / union
  jaccard_sim <- common / (unique1 + unique2 - common)
  
  data.frame(
    n_vec1 = unique1,
    n_vec2 = unique2,
    n_unique_both = total,
    n_common = common,
    n_unique_vec1 = distinct1,
    n_unique_vec2 = distinct2,
    jaccard_similarity = jaccard_sim,
    intersecting_elements = I(list(intersect(vec1, vec2))) # prevents vector being expanded
    )
  
}


compute_stats(c(1,2,3), c(2,3,4))


#' Compute Overlap Statistics for APA Targets
#'
#' This function computes the overlap statistics for APA_IDs and Gene, Chr combinations
#' between the results of running `define_apa_targets` on two separate data frames.
#'
#' @param list1 The first list of data frames, output of `define_apa_targets`.
#' @param list2 The second list of data frames, output of `define_apa_targets`.
#'
#' @return A dataframe containing overlap statistics for all comparisons.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' result1 <- define_apa_targets(dexseq_rep_data1, pau_data1)
#' result2 <- define_apa_targets(dexseq_rep_data2, pau_data2)
#' overlap_stats <- compute_overlap_statistics(result1, result2)
#' }
#'
#' @export
compute_overlap_statistics <- function(list1, list2) {
  
  # Ensure both lists have the same length
  if (length(list1) != length(list2)) {
    stop("The two lists must have the same length.")
  }
  
  # Map over pairs of dataframes and compute the statistics
  stats <- map2(list1, list2, ~ {
    
    # Extract unique APA_IDs and Gene, Chr combinations
    apa_ids_list1 <- unique(.x$APA_ID)
    apa_ids_list2 <- unique(.y$APA_ID)
    
    gene_chr_list1 <- unique(paste(.x$Gene, .x$Chr, sep = "_"))
    gene_chr_list2 <- unique(paste(.y$Gene, .y$Chr, sep = "_"))
    
    # Compute statistics for APA_IDs
    apa_id_stats <- compute_stats(apa_ids_list1, apa_ids_list2)
    
    # Compute statistics for Gene, Chr combinations
    gene_chr_stats <- compute_stats(gene_chr_list1, gene_chr_list2)
    
    # Combine APA_ID and Gene_Chr stats into a single dataframe
    bind_rows(list("APA_ID" = apa_id_stats, "Gene_Chr" = gene_chr_stats),
              .id = "comparison_value")
  }
  )
  
  # Combine all pairwise comparisons into a single dataframe
  result <- bind_rows(stats, .id = "comparison_type")
  
  return(result)
}




#' Generate Direction Bar Plot for PSI Data
#'
#' This function processes a named list of dataframes containing LABRAT psi columns, counts the number of shortening/lengthening events (performing binomial test with null p = 0.5) and generates a bar plot visualising the different frquency of shortening/lengthening events for each input dataframe
#'
#' @param data_list Named list of dataframes, each containing the columns `Gene`, `Chr`, `deltaPSImean.HOM_WT`, and `deltaPSImedian.HOM_WT`.
#' @param padj_method Character string specifying the p-value adjustment method for the binomial test (default is "none").
#' @param plot_title Optional title for the plot (default is NULL).
#' @param fill_values Character vector specifying the fill colors (shortening/lengthening) for the bar plot (default is c("#b2df8a", "#1f78b4")).
#' @param signif_label_size Numeric value specifying the size of the significance label (default is 6).
#' @param signif_label_vjust Numeric value specifying the vertical adjustment for the significance label (default is 10).
#' @param plot_base_size Numeric value specifying the base size for the plot text elements (default is 14).
#' @param id_column Character string specifying the name of the column used to identify the dataframes (populated by list names) (default is "mutant").
#' 
#' @return A named list containing the following elements:
#' \item{plot}{The generated ggplot object.}
#' \item{binom_test}{Dataframe with the binomial test results for each input dataframe (differentiated by id_column).}
#' \item{counts_df}{Dataframe with the counts of direction categories for each input dataframe (differentiated by id_column).}
#' \item{dirn_df}{Dataframe with the direction data (i.e. assignment of shortening/lengthening) for each input dataframe (differentiated by id_column).}
#'
#' @import dplyr
#' @import purrr
#' @import rstatix
#' @import ggplot2
#' @import ggprism
#' @export
psi_direction_bar_plot <- function(data_list,
                                   padj_method = "none", 
                                   plot_title = NULL, 
                                   fill_values = c("#b2df8a", "#1f78b4"),
                                   signif_label_size = 6, signif_label_vjust = 10,
                                   plot_base_size = 14,
                                   id_column = "mutant") {

  # Define the minimal required columns
  required_columns <- c("Gene", "Chr", "deltaPSImean.HOM_WT", "deltaPSImedian.HOM_WT")
  
  # Check that all required columns are present in each dataframe
  check_columns <- function(df) {
    missing_columns <- setdiff(required_columns, colnames(df))
    if (length(missing_columns) > 0) {
      stop(paste("The following required columns are missing:", paste(missing_columns, collapse = ", ")))
    }
  }
  
  # Apply the column check to each dataframe in the list
  walk(data_list, check_columns)
  
  
  # Essential - counts for each mutant in shortening/lengthening category (+ calculate percentage)
  # Then perform binomial test on counts (providing shortening counts), add group1, group2 labels
  # bind the p-value cols together, then calc y position as done before
  ## --> have then defined p value df
    
  # Assign shortening/lengthening based on sign of psi for all genes
  dirn_df_list <- map(data_list, ~ {
    .x %>%
      mutate(
        dirn_mean = if_else(sign(deltaPSImean.HOM_WT) == 1, "Lengthening", "Shortening"),
        dirn_median = if_else(sign(deltaPSImedian.HOM_WT) == 1, "Lengthening", "Shortening"),
        across(starts_with("dirn_me"), ~ factor(.x, levels = c("Shortening", "Lengthening"))),
        dirn_consistent = dirn_mean == dirn_median
      ) %>%
      drop_na(dirn_mean) %>%
      distinct(Gene, Chr, dirn_mean, dirn_median, dirn_consistent)
  })
  
  # Count number of shortening/lengthening events
  counts_df_list <- map(dirn_df_list, ~ {
    .x %>%
      count(dirn_mean) %>%
      mutate(perc = (n / sum(n)) * 100)
  })
  
  # Perform binomial test on each counts_df in the list (+ prep for )
  binom_tests_list <- map(counts_df_list, ~ {
    .x %>%
      pull(n) %>%
      binom_test(detailed = TRUE) %>%
      mutate(group1 = "Shortening", group2 = "Lengthening"
             )
  })
  
  # Combine the lists of dataframes into single dataframe with an id column
  combined_dirn_df <- bind_rows(dirn_df_list, .id = id_column) 
  combined_counts_df <- bind_rows(counts_df_list, .id = id_column) %>%
    mutate("{id_column}" := fct_inorder(!!sym(id_column)))

  # combine binomial test results, adjusting for multiple comparisons + adding plot label location
  combined_binom_tests <- bind_rows(binom_tests_list, .id = id_column) %>%
    mutate("{id_column}" := fct_inorder(!!sym(id_column))) %>%
    adjust_pvalue(output.col = "p.adj", method = padj_method) %>%
    add_significance(p.col = "p.adj", output.col = "p.adj.signif") %>%
    # define with respect to all groups (so standard height for significance layers)
    mutate(y.position = max(estimate * 100, (1 - estimate) * 100) + signif_label_vjust)
  
  # Generate the plot
  plot <- ggplot(combined_counts_df, aes(x = !!sym(id_column), y = perc, label = n)) +
    geom_col(aes(fill = dirn_mean), position = "dodge") +
    geom_text(aes(group = dirn_mean), position = position_dodge(width = 0.9), vjust = -1) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
    ggprism::add_pvalue(combined_binom_tests, label = "p.adj.signif", x = id_column, label.size = signif_label_size) +
    scale_fill_manual(values = fill_values) +
    labs(title = plot_title,
         x = "Mutant",
         y = "Target APAs (%)",
         fill = "3'UTR length change") +
    theme_classic(base_size = plot_base_size) +
    theme(legend.position = "bottom")
  
  # Return a named list with the results
  list(
    plot = plot,
    binom_test = combined_binom_tests,
    counts_df = combined_counts_df,
    dirn_df = combined_dirn_df
  )
}


#' Perform Chi-squared test for two-way comparison quadrant counts
#'
#' This function calculates the quadrant counts from the two-way comparison scatter plot and performs a Chi-squared test 
#' to determine if the distribution of counts across quadrants is different from 
#' what would be expected by chance. Note that the expected counts are scaled by the observed frequency of positive/negative changes in each of the group, rather than assuming equal P of all 4 quadrants.
#' The quadrants are defined based on the direction of change (sign) in two delta columns.
#'
#' @param data A data frame containing the input data (plot_df output from twoway_comparison_scatter).
#' @param group1_delta_col A string specifying the name of the first delta column.
#' @param group2_delta_col A string specifying the name of the second delta column.
#' @param group1_lab A string specifying the label for the first group.
#' @param group2_lab A string specifying the label for the second group.
#' @param pos_lab A string specifying the label for positive direction (default is "longer").
#' @param neg_lab A string specifying the label for negative direction (default is "shorter").
#'
#' @return A named list containing:
#' \item{chisq_test}{The Chi-squared test result (dataframe, output of rstatix::chisq_test)}
#' \item{chisq_test_desc}{Descriptives of the Chi-squared test (output of rstatix::chisq_descriptives)}
#' \item{quadrant_counts}{The counts of events falling into each quadrant (dataframe)}
#' \item{obs_quadrant_dirn_counts}{Observed quadrant direction counts (matrix, input to rstatix::chisq_test)}
#' \item{exp_quadrant_dirn_counts}{Expected quadrant direction counts (matrix, input to rstatix::chisq_test)}
quadrant_count_chisq <- function(data, group1_delta_col, group2_delta_col, group1_lab, group2_lab, pos_lab = "longer", neg_lab = "shorter") {
  
  required_columns <- c("gene_chr", "overlap_group", group1_delta_col, group2_delta_col)
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing from the data:", paste(missing_columns, collapse = ", ")))
  }
  
  # define quadrant labels
  # first split by 'half' for each group
  lab_group1_up <- paste(group1_lab, pos_lab, sep = "")
  lab_group1_down <- paste(group1_lab, neg_lab, sep = "")
  lab_group2_up <- paste(group2_lab, pos_lab, sep = "")
  lab_group2_down <- paste(group2_lab, neg_lab, sep = "")
  
  # now define quadrants
  lab_group1_up_group2_up <- paste(lab_group1_up, lab_group2_up, sep = "__")
  lab_group1_up_group2_down <- paste(lab_group1_up, lab_group2_down , sep = "__")
  lab_group1_down_group2_up <- paste(lab_group1_down, lab_group2_up, sep = "__")
  lab_group1_down_group2_down <- paste(lab_group1_down, lab_group2_down, sep = "__")
  
  # count number of all events falling into each quadrant
  int_quadrant_counts <- data %>%
    mutate(dirn_group = case_when(
      sign(!!sym(group1_delta_col)) == 1 & sign(!!sym(group2_delta_col)) == 1 ~ lab_group1_up_group2_up,
      sign(!!sym(group1_delta_col)) == 1 & sign(!!sym(group2_delta_col)) == -1 ~ lab_group1_up_group2_down,
      sign(!!sym(group1_delta_col)) == -1 & sign(!!sym(group2_delta_col)) == 1 ~ lab_group1_down_group2_up,
      sign(!!sym(group1_delta_col)) == -1 & sign(!!sym(group2_delta_col)) == -1 ~ lab_group1_down_group2_down
      )
    ) %>%
    # some genes not expressed/detected in multiple groups
    drop_na(dirn_group) %>%
    distinct(gene_chr, overlap_group, dirn_group) %>%
    count(overlap_group, dirn_group) 
  
  exp_quadrants <- tibble(dirn_group = c(lab_group1_up_group2_up,lab_group1_up_group2_down, lab_group1_down_group2_up, lab_group1_down_group2_down))
  
  quadrant_counts <- exp_quadrants %>%
    # add 0 counts for missing quadrants (prevents warning with chisq test (ensures symmetric matrix))
    left_join(int_quadrant_counts, by = "dirn_group") %>%
    replace_na(list("n" = 0)) %>%
    group_by(dirn_group) %>%
    summarise(n_tot = sum(n)) %>%
    mutate(frac = n_tot / sum(n_tot))
  
  # get sum of counts in each direction for each group (i.e. short counts, long counts)
  quadrant_dirn_counts <- quadrant_counts %>%
    separate(dirn_group, into = c(group1_lab, group2_lab), sep = "__")
  
  group1_dirn_counts <- quadrant_dirn_counts %>%
    group_by(!!sym(group1_lab)) %>%
    summarise(n = sum(n_tot)) %>%
    mutate(propn = n / sum(n)) %>%
    rename(half_label = !!sym(group1_lab))
  
  group2_dirn_counts <- quadrant_dirn_counts %>%
    group_by(!!sym(group2_lab)) %>%
    summarise(n = sum(n_tot)) %>%
    mutate(propn = n / sum(n)) %>%
    rename(half_label = !!sym(group2_lab))
  
  # get proportions for each group + dirn in a single df
  dirn_propn <- bind_rows(group1_dirn_counts, group2_dirn_counts) %>%
    select(-n) %>%
    pivot_wider(names_from = half_label, values_from = propn)
  
  # get total number of genes analysed
  pop_size <- sum(quadrant_dirn_counts$n_tot)
  
  # calculate expected freqs based on frequency of longer/shorter in input
  exp_quadrant_dirn_counts <- matrix(c(
    pop_size * dirn_propn[[lab_group1_up]] * dirn_propn[[lab_group2_up]],
    pop_size * dirn_propn[[lab_group1_up]] * dirn_propn[[lab_group2_down]],
    pop_size * dirn_propn[[lab_group1_down]] * dirn_propn[[lab_group2_up]],
    pop_size * dirn_propn[[lab_group1_down]] * dirn_propn[[lab_group2_down]]
  ), nrow = 2, byrow = TRUE)

  # put observed quadrant counts in wide format (more clearly assign obs count matrix)
  quadrant_counts_wide <- quadrant_counts %>%
    select(dirn_group, n_tot) %>%
    pivot_wider(names_from = dirn_group, values_from = n_tot)

  obs_quadrant_dirn_counts <- matrix(c(
    quadrant_counts_wide[[lab_group1_up_group2_up]],
    quadrant_counts_wide[[lab_group1_up_group2_down]],
    quadrant_counts_wide[[lab_group1_down_group2_up]],
    quadrant_counts_wide[[lab_group1_down_group2_down]]
  ), nrow = 2, byrow = TRUE)
  
  # Perform chi-sq test
  chisq_test_result <- chisq_test(obs_quadrant_dirn_counts,
                                  p = exp_quadrant_dirn_counts / pop_size, # convert expected frequencies to expected proportions
                                  rescale.p = TRUE,
                                  correct = TRUE)
  
  # Output obs, exp freqs + residuals
  chisq_test_desc <- chisq_descriptives(chisq_test_result)
  
  
  return(list(
    chisq_test = chisq_test_result,
    chisq_test_desc = chisq_test_desc,
    quadrant_counts = quadrant_counts,
    obs_quadrant_dirn_counts = obs_quadrant_dirn_counts,
    exp_quadrant_dirn_counts = exp_quadrant_dirn_counts
  ))
  

}


#' Wrapper function for two-way comparison of target genes (scatter plot, quadrant count analysis)
#'
#' @description This function takes in two data frames, annotates overlapping and unique entries, merges additional information from original data frames, and generates a scatter plot. It also performs chi-squared tests to assess the enrichment of points in specific quadrants.
#'
#' @param df1 Data frame 1 of targets only to be compared.
#' @param df2 Data frame 2 of targets only to be compared.
#' @param orig_df1 Original data frame 1 for merging change in usage for df2-specific targets
#' @param orig_df2 Original data frame 1 for merging change in usage for df1-specific targets
#' @param df1_lab Label to identify data frame 1.
#' @param df2_lab Label to identify data frame 2.
#' @param common_ids Vector of common gene_chr identifiers.
#' @param id_col Column name for gene_chr (default is "gene_chr").
#' @param delta_col Column name containing change in usage (default is "deltaPSImean.HOM_WT").
#' @param id_cols Vector of column names used as identifiers (default is c("Gene", "Chr", "Gene_Name", "APA_ID")).
#'
#' @return A list containing:
#' \describe{
#'   \item{scatter_plot}{ggplot2 object representing the scatter plot.}
#'   \item{plot_df}{Data frame used to generate scatter plot.}
#'   \item{quadrant_counts}{Data frame of quadrant counts for each target set. target_set column differentiates different subsets of targets.}
#'   \item{chisq_result}{Data frame of chi-squared test results with adjusted p-values and significance. target_set column differentiates different subsets of targets.}
#'   \item{chisq_descriptives}{Data frame of chi-squared test descriptive statistics. target_set column differentiates different subsets of targets.}
#'   \item{chisq_observed}{Data frame of chi-squared test observed counts. target_set column differentiates different subsets of targets.}
#'   \item{chisq_expected}{Data frame of chi-squared test expected counts (weighted by observed frequencies of direction of change in usage). target_set column differentiates different subsets of targets.}
#' }
twoway_comparison_scatter <- function(df1, df2, orig_df1, orig_df2, df1_lab, df2_lab, common_ids, id_col = "gene_chr", delta_col = "deltaPSImean.HOM_WT", id_cols = c("Gene", "Chr", "Gene_Name", "APA_ID")) {
  
  
  # columns subset for joining targets with original, unfiltered df
  id_plus_delta_cols <- c(id_cols, delta_col)
  all_id_plus_delta_cols <- c(id_cols, id_col, delta_col)
  
  df1_suffix <- paste(".", df1_lab, sep = "")
  df2_suffix <- paste(".", df2_lab, sep = "")
  

  # for both dfs, annotate common ids then add info from original df of other dataset
  olap_df1 <- df1 %>%
    mutate(overlap_group = if_else(gene_chr %in% common_ids, "Both", df1_lab)) %>%
    # add in info from df2
    select(all_of(all_id_plus_delta_cols), overlap_group) %>%
    left_join(select(orig_df2, all_of(id_plus_delta_cols)), by = id_cols,
              suffix = c(df1_suffix, df2_suffix))
  
  olap_df2 <- df2 %>%
    mutate(overlap_group = if_else(gene_chr %in% common_ids, "Both", df2_lab)) %>%
    # add in info from df1
    select(all_of(all_id_plus_delta_cols), overlap_group) %>%
    left_join(select(orig_df1, all_of(id_plus_delta_cols)), by = id_cols,
              suffix = c(df2_suffix, df1_suffix))
  
  # 
  # combine into single df ready for plotting
  olap_comb <- bind_rows(olap_df1, olap_df2) %>%
    # 'both' group will be duplicated
    distinct(!!sym(id_col), .keep_all = T) %>%
    mutate(plot_alpha = if_else(overlap_group == "Both", 1, 0.8))
  
  # return(olap_comb)
  
  delta_df1_col <- paste(delta_col, df1_lab, sep = ".")
  delta_df2_col <- paste(delta_col, df2_lab, sep = ".")
  
  # # Scatter plot
  scatter_plot <- ggplot() +
    geom_point(data = filter(olap_comb, overlap_group == df1_lab), aes(x = !!sym(delta_df1_col), y = !!sym(delta_df2_col), color = df1_lab), alpha = 0.6) +
    geom_point(data = filter(olap_comb, overlap_group == df2_lab), aes(x = !!sym(delta_df1_col), y = !!sym(delta_df2_col), color = df2_lab), alpha = 0.6) +
    geom_point(data = filter(olap_comb, overlap_group == "Both"), aes(x = !!sym(delta_df1_col), y = !!sym(delta_df2_col), color = "Both"), alpha = 1, size = 2)+
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    scale_color_manual(name = "Overlap Group", values = c("M323K" = "#ef8a62", "F210I" = "#67a9cf", "Both" = "black", "Q331K" = "#7570b3")) +
    scale_y_continuous(limits = c(-1,1)) +
    scale_x_continuous(limits = c(-1,1)) +
    labs(title = NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  
  
  # Perform chi-squared test for different target sets to assess enrichment in specific quadrants
  
  # Generate dfs containing quadrant assignments for all target combinations
  olap_comb_sets <- c(df1_lab, df2_lab, "Both") %>%
    set_names() %>%
    map(~ filter(olap_comb, overlap_group == .x))
  # add 
  olap_comb_sets <- c(olap_comb_sets, list("All" = olap_comb))
  
  chisq_comb_sets <- map(olap_comb_sets,
                         ~ quadrant_count_chisq(.x, delta_df1_col, delta_df2_col, df1_lab, df2_lab)
                         )
  
  # Pull out quadrant counts, chisq result and chisq stats, obs_ for each 
  quadrant_counts_comb <- map(chisq_comb_sets, "quadrant_counts") %>%
    bind_rows(.id = "target_set")
  
  chisq_res_comb <- map(chisq_comb_sets, "chisq_test") %>%
    bind_rows(.id = "target_set") %>%
    adjust_pvalue(output.col = "p.adj", method = "BH") %>%
    add_significance(p.col = "p.adj", output.col = "p.adj.signif") %>%
    mutate(method.padj = "BH")
  
  chisq_desc_comb <- map(chisq_comb_sets, "chisq_test_desc") %>%
    bind_rows(.id = "target_set")
  
  chisq_obs_comb <- map(chisq_comb_sets, "obs_quadrant_dirn_counts") %>%
    map(as_tibble) %>%
    bind_rows(.id = "target_set")
  
  chisq_exp_comb <- map(chisq_comb_sets, "exp_quadrant_dirn_counts") %>%
    map(as_tibble) %>%
    bind_rows(.id = "target_set")
  
  
  
  return(list(scatter_plot = scatter_plot,
              plot_df = olap_comb,
              quadrant_counts = quadrant_counts_comb,
              chisq_result = chisq_res_comb,
              chisq_descriptives = chisq_desc_comb,
              chisq_observed = chisq_obs_comb,
              chisq_expected = chisq_exp_comb
              ))

  
}


#' This function generates Euler plots for both the 'target' subsets and all data subsets from two data frames. It checks for the presence of required columns and ensures the `target_col` is boolean.
#'
#' @param df1 A data frame.
#' @param df2 A 2nd data frame to visualise overlap with df1.
#' @param df1_lab A label for the first data frame to be used in the plot.
#' @param df2_lab A label for the second data frame to be used in the plot.
#' @param id_col A string specifying the name of the column containing unique identifiers that wish to compare (default is "gene_chr").
#' @param target_col A string specifying the name of the boolean column used to identify the subset of targets in each df (default is "target").
#' @param targets_prefix A prefix for the 'targets' identifier group in the list & plot (default is "Targets").
#' @param all_prefix A prefix for the all 'identifiers' group in the list & plot (default is "All").
#' @param label_sep A string used to separate the prefixes and labels in the list names and plot (default is " ").
#' @param labels Logical indicating whether to display labels in the Euler plot (default is TRUE).
#' @param quantities Logical indicating whether to display quantities in the Euler plot (default is TRUE).
#' @param fills Logical indicating whether to fill the Euler plot with colors (default is TRUE).
#'
#' @return A list containing the Euler plots and lists of gene_chr values:
#' \describe{
#'   \item{all_plot}{Euler plot for all identifier groups (all + targets).}
#'   \item{target_plot}{Euler plot for targets subsets (targets only).}
#'   \item{all_list}{Named list of character vectors for all data subsets (i.e. input to all_plot).}
#'   \item{targets_list}{Named list of character vectors for targets subsets (i.e. input to target_plot).}
#' }
#'
#' @importFrom dplyr filter distinct pull mutate
#' @importFrom eulerr euler plot
#' @importFrom rlang sym
#' @export
euler_plot_wrapper <- function(df1, df2, df1_lab, df2_lab, id_col = "gene_chr", target_col = "target", targets_prefix = "Targets", all_prefix = "All", label_sep = " ", labels = TRUE, quantities = TRUE, fills = TRUE) {
  # Check if all expected columns are present in df1
  expected_cols <- c(id_col, target_col)
  
  missing_cols_df1 <- setdiff(expected_cols, colnames(df1))
  if (length(missing_cols_df1) > 0) {
    stop(paste("Missing columns in df1:", paste(missing_cols_df1, collapse = ", ")))
  }
  
  # Check if all expected columns are present in df2
  missing_cols_df2 <- setdiff(expected_cols, colnames(df2))
  if (length(missing_cols_df2) > 0) {
    stop(paste("Missing columns in df2:", paste(missing_cols_df2, collapse = ", ")))
  }
  
  # Check if target_col is boolean
  if (!is.logical(df1[[target_col]]) || !is.logical(df2[[target_col]])) {
    stop("target_col must be a boolean column")
  }
  
  # Create character vectors for 'targets' and 'all' in df1
  targets_df1 <- df1 %>%
    filter(!!sym(target_col) == TRUE) %>%
    distinct(!!sym(id_col)) %>%
    pull(!!sym(id_col))
  
  all_df1 <- df1 %>%
    distinct(!!sym(id_col)) %>%
    pull(!!sym(id_col))
  
  # Create character vectors for 'targets' and 'all' in df2
  targets_df2 <- df2 %>%
    filter(!!sym(target_col) == TRUE) %>%
    distinct(!!sym(id_col)) %>%
    pull(!!sym(id_col))
  
  all_df2 <- df2 %>%
    distinct(!!sym(id_col)) %>%
    pull(!!sym(id_col))
  
  # Generate names for the list based on provided arguments
  list_names_all <- c(
    paste(targets_prefix, df1_lab, sep = label_sep),
    paste(all_prefix, df1_lab, sep = label_sep),
    paste(targets_prefix, df2_lab, sep = label_sep),
    paste(all_prefix, df2_lab, sep = label_sep)
  )
  
  list_names_targets <- c(
    paste(targets_prefix, df1_lab, sep = label_sep),
    paste(targets_prefix, df2_lab, sep = label_sep)
  )
  
  # Create the named lists for eulerr
  gene_chr_list <- list(
    targets_df1,
    all_df1,
    targets_df2,
    all_df2
  )
  names(gene_chr_list) <- list_names_all
  
  gene_chr_list_targets <- list(
    targets_df1,
    targets_df2
  )
  names(gene_chr_list_targets) <- list_names_targets
  
  # Generate the Euler plots (one for targets + non-targets, one for just targets)
  euler_fit_all <- euler(gene_chr_list)
  all_plot <- plot(euler_fit_all, labels = labels, quantities = quantities, fills = fills)
  
  euler_fit_targets <- euler(gene_chr_list_targets)
  target_plot <- plot(euler_fit_targets, labels = labels, quantities = quantities, fills = fills)
  
  # Return the plots and the lists of targets
  return(list(all_plot = all_plot, target_plot = target_plot, all_list = gene_chr_list, targets_list = gene_chr_list_targets))
}
