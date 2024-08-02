library(tidyverse)
library(RedRibbon)
set.seed(123)


#' Analyze Data and Generate RedRibbon Plot
#'
#' This function performs analysis on two input dataframes and generates a RedRibbon plot.
#' 
#' @param a First input dataframe.
#' @param b Second input dataframe.
#' @param a_id_col Column name in dataframe \code{a} containing the feature identifier.
#' @param b_id_col Column name in dataframe \code{b} representing the feature identifier.
#' @param a_score_col column name in dataframe \code{a} containing the 'score' column (against which to rank features).
#' @param b_score_col Score column name in dataframe \code{b} containing the 'score' column (against which to rank features).
#' @param quadrants_niter Number of iterations for quadrant coordinate analysis. Default is 200.
#' @param quadrants_permutation Logical indicating whether to perform permutation for quadrant analysis. Default is TRUE.
#' @param plot_axes Labels for the x and y axes in the RedRibbon grid heatmap. Default is c("A", "B").
#'
#' @return A list containing:
#' \describe{
#'   \item{red_ribbon_result}{RedRibbon analysis result.}
#'   \item{quadrants_object}{Quadrants object obtained from quadrant analysis.}
#'   \item{combined_data}{Combined dataframe resulting from left join of input dataframes \code{a} and \code{b}.}
#'   \item{quadrant_padj}{Dataframe containing adjusted p-values for each quadrant.}
#'   \item{quadrant_intersections}{list of quadrant-by-quadrant dataframes of \code{combined_data} subsetted for identified overlap thresholds.}
#'   \item{plot}{Generated RedRibbon plot.}
#' }
redribbon_wrapper <- function(a,
                              b,
                              a_id_col,
                              b_id_col,
                              a_score_col,
                              b_score_col,
                              quadrants_niter = 200,
                              quadrants_permutation = TRUE,
                              plot_axes = c("A", "B")
) {
  
  
  # Perform left join
  comb_df <- left_join(
    select(a, id = !!sym(a_id_col), a = !!sym(a_score_col), Chr),
    select(b, id = !!sym(b_id_col), b = !!sym(b_score_col), Chr),
    by = c("id", "Chr")
  )
  
  if (any(is.na(comb_df$b))) {
    stop("NA values found in joined df, please ensure input dfs have completely matching entries & score columns do not have NAs")
  }
  
  
  # RedRibbon analysis
  rr <- RedRibbon(comb_df, enrichment_mode = "hyper-two-tailed")
  
  # Quadrant analysis
  quad <- quadrants(rr, algorithm = "ea", permutation = quadrants_permutation, whole = FALSE, niter = quadrants_niter)
  
  # Extract dfs subsetted for overlapping genes and adjusted p-values
  quadrant_intersections <- map(quad, "positions") %>%
    map(~ comb_df[.x,])
  
  padj_df <- map(quad, "padj") %>%
    enframe(name = "quadrant", value = "padj") %>%
    unnest(padj)
  
  # Plotting
  gg <- ggRedRibbon(
    rr,
    quadrants = quad,
    labels = plot_axes,
    repel.force = 250,
    show.quadrants = TRUE,
    show.pval = TRUE
  ) +
    coord_fixed(ratio = 1, clip = "off")
  
  return(list(
    red_ribbon_result = rr,
    quadrants_object = quad,
    combined_data = comb_df,
    quadrant_padj = padj_df,
    quadrant_intersections = quadrant_intersections,
    plot = gg
  ))
}



# psi tables containing condition means and deltas
psi_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.psi.mean_median.tsv")
psi_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.psi.mean_median.tsv")

# Summary of PAS consistency for genes betweeen two datasets
psi_consistency <- read_tsv("processed/2024-02-29/analysis/mef_comparisons/f210i_vs_m323k.psi_consistency.summary.tsv")

# processed DEXSeq results tables produced by QAPA_snakemake pipeline
dexseq_m323k <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")
dexseq_f210i <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")

# Representative PAS
# processed DEXSeq results tables filtered to 2 representative PAS produced by mef_representative_pas.R
dexseq_rep_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")
dexseq_rep_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.rep_pas.tsv")
rep_consistency <- read_tsv("processed/2024-02-29/analysis/mef_comparisons/f210i_vs_m323k.representative_pas_consistency.summary.tsv")

outdir <- "processed/2024-02-29/analysis/mef_comparisons/redribbon"
if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}


# filter for genes with same analysed + rep PAS for the two MUT datasets
rep_consistency_same <- rep_consistency %>%
  filter(same_pas)

count(rep_consistency_same, sig.either)

# subset rep dfs to those with same PAS
dexseq_rep_f210i_same <- rep_consistency_same %>%
  distinct(Gene, Chr) %>%
  left_join(dexseq_rep_f210i, by = c("Gene", "Chr"))

dexseq_rep_m323k_same <- rep_consistency_same %>%
  distinct(Gene, Chr) %>%
  left_join(dexseq_rep_m323k, by = c("Gene", "Chr"))

# assign direction based on proximal site
dexseq_rep_f210i_same_tformed <- dexseq_rep_f210i_same %>%
  group_by(Gene, Chr) %>%
  mutate(min_pval = min(pvalue)) %>%
  slice_min(Pas_Number, n=1) %>%
  ungroup() %>%
  mutate(log10_pvalue = -log10(min_pval),
         log10_pvalue_signed = log10_pvalue * sign(log2fold_HOM_WT)
         ) %>%
  distinct(Gene, Chr, log10_pvalue, log10_pvalue_signed)

dexseq_rep_m323k_same_tformed <- dexseq_rep_m323k_same %>%
  group_by(Gene, Chr) %>%
  mutate(min_pval = min(pvalue)) %>%
  slice_min(Pas_Number, n=1) %>%
  ungroup() %>%
  mutate(log10_pvalue = -log10(min_pval),
         log10_pvalue_signed = log10_pvalue * sign(log2fold_HOM_WT)
  ) %>%
  distinct(Gene, Chr, log10_pvalue, log10_pvalue_signed)

# Run RedRibbon comparing M323K vs F210I
rep_pau_m323k_f210i_rr <- redribbon_wrapper(a = dexseq_rep_m323k_same_tformed,
                  b = dexseq_rep_f210i_same_tformed,
                  a_id_col = "Gene",
                  b_id_col = "Gene",
                  a_score_col ="log10_pvalue_signed",
                  b_score_col = "log10_pvalue_signed",
                  quadrants_niter = 100,
                  quadrants_permutation = TRUE,
                  plot_axes = c("MEF M323K", "MEF F210I"))

# rep_pau_m323k_f210i_rr$red_ribbon_result
rep_pau_m323k_f210i_rr$plot
rep_pau_m323k_f210i_rr$quadrant_padj
rep_pau_m323k_f210i_rr$quadrant_intersections

rep_pau_m323k_f210i_rr$plot +
  theme(legend.position = "top")

ggsave(filename = file.path(outdir, "2024-06-29_mef_f210i_m323k.redribbon.pau_representative_pas.png"),
       plot = rep_pau_m323k_f210i_rr$plot,
       device = "png",
       units = "mm", 
       width = 150,
       height = 150
       )

ggsave(filename = file.path(outdir, "2024-06-29_mef_f210i_m323k.redribbon.pau_representative_pas.pdf"),
       plot = rep_pau_m323k_f210i_rr$plot,
       device = "pdf",
       units = "mm", 
       width = 150,
       height = 150
)

## Repeat for psi analysis
psi_consistency_same <- psi_consistency %>%
  filter(same_pas)

dexseq_psi_f210i_same <- psi_consistency_same %>%
  distinct(Gene, Chr) %>%
  left_join(distinct(dexseq_f210i, Gene, Chr, pvalue, padj), by = c("Gene", "Chr")) %>%
  left_join(psi_f210i, by = c("Gene", "Chr"))
  

dexseq_psi_m323k_same <- psi_consistency_same %>%
  distinct(Gene, Chr) %>%
  left_join(distinct(dexseq_m323k, Gene, Chr, pvalue, padj), by = c("Gene", "Chr")) %>%
  left_join(psi_m323k, by = c("Gene", "Chr"))


# # Calculating -log10 transformed p values and multiplying by sign of deltaPSImean/madian
dexseq_psi_m323k_tform <- dexseq_psi_m323k_same %>%
  mutate(log10_pvalue = -log10(pvalue),
         log10_pvalue_mean_signed = log10_pvalue * sign(deltaPSImean.HOM_WT),
         log10_pvalue_median_signed = log10_pvalue * sign(deltaPSImedian.HOM_WT),
         consistent_sign = sign(log10_pvalue_mean_signed) == sign(log10_pvalue_median_signed))

dexseq_psi_f210i_tform <- dexseq_psi_f210i_same %>%
  mutate(log10_pvalue = -log10(pvalue),
         log10_pvalue_mean_signed = log10_pvalue * sign(deltaPSImean.HOM_WT),
         log10_pvalue_median_signed = log10_pvalue * sign(deltaPSImedian.HOM_WT),
         consistent_sign = sign(log10_pvalue_mean_signed) == sign(log10_pvalue_median_signed))

# Selecting representative PAS/p-value for each gene
psi_rep_m323k <- dexseq_psi_m323k_tform %>%
  group_by(Gene, Chr) %>%
  slice_max(log10_pvalue, n=1) %>%
  ungroup()

psi_rep_f210i <- dexseq_psi_f210i_tform %>%
  group_by(Gene, Chr) %>%
  slice_max(log10_pvalue, n=1) %>%
  ungroup()


tmp_ids <- psi_rep_f210i %>%
  filter(is.na(log10_pvalue_mean_signed)) %>%
  pull(Gene)

# check that IDs have all NA psi values, not just for max log10_pval
dexseq_psi_f210i_tform %>%
  filter(Gene %in% tmp_ids)


# filter out genes with NA score values, assert that all IDs still match
psi_rep_f210i <- psi_rep_f210i %>%
  filter(!is.na(log10_pvalue_mean_signed))

psi_rep_m323k <- psi_rep_m323k %>%
  filter(!is.na(log10_pvalue_mean_signed))

# Find all unique IDs per dataset, filter from further analysis
psi_inconsistent_na_ids <- unique(c(setdiff(psi_rep_f210i$Gene, psi_rep_m323k$Gene), setdiff(psi_rep_m323k$Gene, psi_rep_f210i$Gene)))
message(paste("Number of unique IDs where NA psi values in one of the datasets - ", length(psi_inconsistent_na_ids)))

psi_rep_f210i <- psi_rep_f210i %>%
  filter(!Gene %in% psi_inconsistent_na_ids)

psi_rep_m323k <- psi_rep_m323k %>%
  filter(!Gene %in% psi_inconsistent_na_ids)

stopifnot( length(unique(c(setdiff(psi_rep_f210i$Gene, psi_rep_m323k$Gene),
                           setdiff(psi_rep_m323k$Gene, psi_rep_f210i$Gene))
                         )
                  ) == 0)

# subset for consistent sign between mean and median psi
psi_rep_f210i_signsame <- psi_rep_f210i %>%
  filter(consistent_sign)

psi_rep_m323k_signsame <- psi_rep_m323k %>%
  filter(consistent_sign)

# Find all unique IDs per dataset, filter from further analysis
psi_signsame_inconsistent_na_ids <- unique(c(setdiff(psi_rep_f210i_signsame$Gene, psi_rep_m323k_signsame$Gene), setdiff(psi_rep_m323k_signsame$Gene, psi_rep_f210i_signsame$Gene)))
message(paste("psi same sign analysis - Number of unique IDs where NA psi values in one of the datasets - ", length(psi_signsame_inconsistent_na_ids)))

psi_rep_f210i_signsame <- psi_rep_f210i_signsame %>%
  filter(!Gene %in% psi_signsame_inconsistent_na_ids)

psi_rep_m323k_signsame <- psi_rep_m323k_signsame %>%
  filter(!Gene %in% psi_signsame_inconsistent_na_ids)

stopifnot( length(unique(c(setdiff(psi_rep_f210i_signsame$Gene, psi_rep_m323k_signsame$Gene),
                           setdiff(psi_rep_m323k_signsame$Gene, psi_rep_f210i_signsame$Gene))
)
) == 0)


rep_psi_m323k_f210i_rr <- redribbon_wrapper(a = psi_rep_m323k,
                                            b = psi_rep_f210i,
                                            a_id_col = "Gene",
                                            b_id_col = "Gene",
                                            a_score_col ="log10_pvalue_mean_signed",
                                            b_score_col = "log10_pvalue_mean_signed",
                                            quadrants_niter = 100,
                                            quadrants_permutation = TRUE,
                                            plot_axes = c("MEF M323K (psi)", "MEF F210I (psi)"))


rep_psi_signsame_m323k_f210i_rr <- redribbon_wrapper(a = psi_rep_m323k_signsame,
                                            b = psi_rep_f210i_signsame,
                                            a_id_col = "Gene",
                                            b_id_col = "Gene",
                                            a_score_col ="log10_pvalue_mean_signed",
                                            b_score_col = "log10_pvalue_mean_signed",
                                            quadrants_niter = 100,
                                            quadrants_permutation = TRUE,
                                            plot_axes = c("MEF M323K (psi)", "MEF F210I (psi)"))


rep_psi_m323k_f210i_rr$quadrant_padj
rep_psi_signsame_m323k_f210i_rr$quadrant_padj
rep_psi_signsame_m323k_f210i_rr$plot
