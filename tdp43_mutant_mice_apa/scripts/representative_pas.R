suppressPackageStartupMessages(library(optparse))

# Define the option list
option_list <- list(
  make_option(
    c("-d", "--dexseq"),
    type = "character",
    help = "Processed DEXSeq results tables produced by QAPA_snakemake pipeline",
    metavar = "character"
  ),
  make_option(
    c("-s", "--significance_threshold"),
    type = "double",
    default = 0.05,
    help = "Adjusted p-value significance threshold (default=%default)",
    metavar = "double"
  ),
  make_option(
    c("-f","--foldchange_column"),
    type = "character",
    default = "log2fold_",
    help = "Prefix of column name containing DEXSeq calculated log2 transformed fold change (default=%default)",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output_prefix"),
    type = "character",
    help = "Prefix for output files ('representative.pas.tsv' for chosen representative PAS, 'representative.total_expression.tsv' for fraction of total expression attributed to representative PAS)",
    metavar = "character"
  )
)

# Create the option parser
opt_parser <- OptionParser(option_list = option_list,
                           description = "Select two representative PAS for each gene from DEXSeq analysis of QAPA quantified PAS.",
                           epilogue = paste("Representative PAS selection strategy:",
                                          "Gene has 2 analysed PAS - keep as is",
                                          "Gene has no PAS passing significance threshold - select two highest expressed PAS (by mean count in all samples, 'exonBaseMean' column)",
                                          "Gene has 1 PAS passing significance threshold - select significant PAS & highest expressed not-significant PAS",
                                          "Gene has 2 PAS passing significance threshold - select both significant PAS. If both significant PAS are changing in same direction, select the highest expressed significant and not-significant PAS (by mean count in all samples, 'exonBaseMean' column)",
                                          "Gene has 3+ PAS passing significance threshold - select two PAS with largest expression changing in opposite directions",
                                          sep = "\n")
                           )

# if no arguments are provided, print the help messages and exit the script
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

# Parse the command line arguments
opt <- parse_args(opt_parser)


# Now know script will be executed, define required functions

suppressPackageStartupMessages(library(tidyverse))

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

#' Select 2 rep PAS for a gene. Where genes have two sig PAS changing in same direction, highest expressed sig and not sig are selecteds
rep_pas_2_sig <- function(data, grouping_cols = c("Gene", "Chr"), padj_col = "padj", log2fc_col = "log2fold_HOM_WT", exprn_col = "exonBaseMean", padj_threshold = 0.05) {
  
  # find genes/groups with opposing dirns in usage for 2 sig PAS
  pre_grp_ids <- distinct(data, across(all_of(grouping_cols)))
  
  result <- data %>%
    filter(!!sym(padj_col) < padj_threshold) %>%
    group_by(across(all_of(grouping_cols))) %>%
    # +1 for +ve, -1 for -ve 
    summarise(sum_dirn = sum(sign(!!sym(log2fc_col)))) %>%
    ungroup() %>%
    # filter for genes/groups with opposing signs for 2 PAS
    filter(sum_dirn == 0) %>%
    select(-sum_dirn)
  
  # 
  post_grp_ids <- distinct(result, across(all_of(grouping_cols)))
  
  # join back PAS entries for genes with 2 opposite dirn sig PAS
  rep_pass <- result %>%
    left_join(filter(data, !!sym(padj_col) < padj_threshold), by = grouping_cols)
  
  if (nrow(pre_grp_ids) == nrow(post_grp_ids)) {
    
    # all have 2 sig PAS in opposite directions
    return(rep_pass)
    
  } else {
    
    message("Genes with 2 sig PAS - genes with 2 sig in same direction observed")
    
    # get missing IDs - report count + actual IDs
    missing_ids <- anti_join(pre_grp_ids, post_grp_ids, by = grouping_cols)
    message(paste("Number of genes affected:", nrow(missing_ids)))
    
    x <- missing_ids %>%
      unite(joined_id, all_of(grouping_cols), sep = " - ") %>%
      pull()
    
    message("Affected IDs:")
    message(paste(x, collapse = "\n"))
    
    # return all PAS for missing genes
    missing_data <- left_join(missing_ids, data, by = grouping_cols)
    
    # select highest expressed PAS passing sig threshold
    missing_sig <- missing_data %>%
      filter(!!sym(padj_col) < padj_threshold) %>%
      group_by(across(all_of(grouping_cols))) %>%
      slice_max(!!sym(exprn_col), n=1) %>%
      ungroup()
    
    # select highest expressed PAS failing sig threshold
    missing_ns <- missing_data %>%
      filter(!!sym(padj_col) >= padj_threshold) %>%
      group_by(across(all_of(grouping_cols))) %>%
      slice_max(!!sym(exprn_col), n=1) %>%
      ungroup()
    
    # combine back with other IDs with opposite dirn sig PAS
    missing_rep <- bind_rows(missing_sig, missing_ns) 
    
    rep_final <- bind_rows(rep_pass, missing_rep) %>%
      arrange(across(all_of(grouping_cols)))
    
    return(rep_final)
    
    
  }
  
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
  rep_2_sig <- rep_pas_2_sig(npasmore_sig_split$nsig_2, grouping_cols, padj_col, log2fc_col, exprn_col, padj_threshold)
  
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


# ANALYSIS

cat("DEXSeq file:", opt$dexseq, "\n")
cat("Significance threshold:", opt$significance_threshold, "\n")
cat("Fold change column prefix:", opt$foldchange_column, "\n")
cat("Output prefix:", opt$output_prefix, "\n")

dexseq <- read_tsv(opt$dexseq, show_col_types=F)
output_prefix <- opt$output_prefix

dexseq_colorder <- colnames(dexseq)

# find fold change column using prefix
fc_column <- dexseq_colorder[str_detect(dexseq_colorder, paste("^", opt$foldchange_column, sep = ""))]
stopifnot("Multiple fold change columns identified for provided prefix"=length(fc_column) == 1)
cat("Inferred fold change column:", fc_column, "\n")

# Number of uniuqe genes before selecting representative APS
genes_init <- distinct(dexseq, Gene, Chr) 

cat("Selecting representative PAS...\n")
rep_pas <- rep_pas_wrapper(dexseq, padj_threshold = opt$significance_threshold, log2fc_col = fc_column)

genes_rep <- distinct(rep_pas$rep_df, Gene, Chr)

# Check that all genes are retained after selecting rep PAS
stopifnot("Different number of genes after selecting representative PAS"=bind_rows(init = genes_init, rep = genes_rep, .id = "source") %>%
  count(Gene, Chr) %>%
  filter(n != 2) %>% nrow() == 0
  )

# Output to file
out_rep <- paste(output_prefix, "representative", "pas" , "tsv", sep = ".")
out_exprn <-paste(output_prefix, "representative", "total_expression" , "tsv", sep = ".")


cat("Output TSV of representative PAS:", out_rep, "\n")
write_tsv(rep_pas$rep_df[dexseq_colorder],
            file = out_rep,
            col_names = T
          )

cat("Output TSV of proportion of total expression attributed to representative PAS:", out_exprn, "\n")
write_tsv(rep_pas$exprn_df,
          file = out_exprn,
          col_names = T
          )


# output to file