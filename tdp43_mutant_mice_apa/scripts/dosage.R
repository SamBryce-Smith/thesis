library(tidyverse)

#' Assign Dosage Pattern
#'
#' Assigns dosage pattern of PAS usage to a dataframe based on the comparison of WT, HET and HOM columns.
#'
#' @param df DataFrame The input dataframe.
#' @param column_name_prefix Character The common prefix for the column names.
#' @param wt_key Character The key for the "wild type" category.
#' @param het_key Character The key for the "heterozygous" category.
#' @param hom_key Character The key for the "homozygous" category.
#' @param separator Character, optional The separator used to concatenate column names. Defaults to "_".
#' @param pattern_column_name Character, optional The name for the output column containing specific dosage pattern labels. Defaults to "pattern".
#' @param pattern_simplified_column_name Character, optional The name for the output column containing simplified dosage pattern labels. Defaults to "pattern_simplified".
#'
#' @return DataFrame The input dataframe with an additional columns {pattern_column_name} representing the dosage pattern and {pattern_simplified_column_name} representing simple consistent/not consistent categories
#'
#' @examples
#' # Creating a dataframe
#' df <- data.frame(
#'   meanPAU_WT = c(1, 2, 3),
#'   meanPAU_HET = c(2, 3, 1),
#'   meanPAU_HOM = c(3, 1, 2)
#' )
#'
#' # Assigning dosage pattern
#' result_df <- assign_dosage_pattern(df, "meanPAU", "WT", "HET", "HOM")
#' print(result_df)
#'
#' @export
assign_dosage_pattern <- function(df, column_name_prefix, wt_key, het_key, hom_key, separator = "_", pattern_column_name = "dosage_pattern_specific", pattern_simplified_column_name = "dosage_pattern_simplified") {
  
  # Construct column names
  wt_column <- paste(column_name_prefix, wt_key, sep = separator)
  het_column <- paste(column_name_prefix, het_key, sep = separator)
  hom_column <- paste(column_name_prefix, hom_key, sep = separator)
  
  # Check if constructed column names are found in the input dataframe
  if (!all(c(wt_column, het_column, hom_column) %in% colnames(df))) {
    stop("One or more constructed column names not found in the input dataframe.")
  }
  
  # Apply dosage patterns (specific and simplified)
  result <- df %>%
    mutate({{pattern_column_name}} := case_when(
      !!sym(wt_column) < !!sym(het_column) & !!sym(het_column) < !!sym(hom_column) ~ paste("Up", "+", "Up", sep = ""),
      !!sym(wt_column) > !!sym(het_column) & !!sym(het_column) > !!sym(hom_column) ~ paste("Down", "+", "Down", sep = ""),
      !!sym(wt_column) < !!sym(het_column) & !!sym(het_column) > !!sym(hom_column) ~ paste("Up", "+", "Down", sep = ""),
      !!sym(wt_column) > !!sym(het_column) & !!sym(het_column) < !!sym(hom_column) ~ paste("Down", "+", "Up", sep = ""),
      !!sym(wt_column) > !!sym(het_column) & !!sym(het_column) == !!sym(hom_column) ~ paste("Down", "+", "Same", sep = ""),
      !!sym(wt_column) < !!sym(het_column) & !!sym(het_column) == !!sym(hom_column) ~ paste("Up", "+", "Same", sep = ""),
      !!sym(wt_column) == !!sym(het_column) & !!sym(het_column) > !!sym(hom_column) ~ paste("Same", "+", "Down", sep = ""),
      !!sym(wt_column) == !!sym(het_column) & !!sym(het_column) < !!sym(hom_column) ~ paste("Same", "+", "Up", sep = ""),
      TRUE ~ "Other"
    )) %>%
    mutate({{pattern_simplified_column_name}} := case_when(str_detect(!!sym(pattern_column_name), "Same") ~ "same_containing",
                                                           grepl("Up", !!sym(pattern_column_name)) & grepl("Down", !!sym(pattern_column_name)) ~ "not_consistent",
                                                           T ~ "consistent")
             )
  
  return(result)
}

# Creating a test dataframe
test_df <- data.frame(
  meanPAU_WT = c(1, 2, 3, 4, 5, 6, 7, 8, 3),
  meanPAU_HET = c(2, 3, 1, 4, 5, 6, 8, 7, 2),
  meanPAU_HOM = c(3, 1, 2, 5, 4, 9, 7, 8, 1)
)

assign_dosage_pattern(test_df, column_name_prefix = "meanPAU", wt_key = "WT", het_key = "HET", hom_key = "HOM",separator = "_",pattern_column_name = "mean.dosage_pattern", pattern_simplified_column_name = "mean.dosage_pattern_simplified")


#
# processed PAU dfs produced by calculate_qapa_means_psis.R
pau_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.mean_pau.tsv")
pau_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.mean_pau.tsv")
pau_m323k_ad6m <- read_tsv("processed/2024-02-29/M323K_Adult_Brain/M323K_Adult_Brain.mean_pau.tsv")

# psi tables containing condition means and deltas
psi_m323k <- read_tsv("processed/2024-02-29/M323K_MEFs/M323K_MEFs.psi.mean_median.tsv")
psi_f210i <- read_tsv("processed/2024-02-29/F210I_MEFs/F210I_MEFs.psi.mean_median.tsv")
psi_m323k_ad6m <- read_tsv("processed/2024-02-29/M323K_Adult_Brain/M323K_Adult_Brain.psi.mean_median.tsv")

# Assign dosage patterns based on mean PAU values
pau_m323k <- assign_dosage_pattern(pau_m323k,
                      column_name_prefix = "mean_PAU",
                      wt_key = "WT",
                      het_key = "HET",
                      hom_key = "HOM",
                      separator = "_",
                      pattern_column_name = "meanPAU.dosage_pattern",
                      pattern_simplified_column_name = "meanPAU.dosage_pattern_simplified")

pau_f210i <- assign_dosage_pattern(pau_f210i,
                      column_name_prefix = "mean_PAU",
                      wt_key = "WT",
                      het_key = "HET",
                      hom_key = "HOM",
                      separator = "_",
                      pattern_column_name = "meanPAU.dosage_pattern",
                      pattern_simplified_column_name = "meanPAU.dosage_pattern_simplified")

pau_m323k_ad6m <- assign_dosage_pattern(pau_m323k_ad6m,
                                   column_name_prefix = "mean_PAU",
                                   wt_key = "WT",
                                   het_key = "HET",
                                   hom_key = "HOM",
                                   separator = "_",
                                   pattern_column_name = "meanPAU.dosage_pattern",
                                   pattern_simplified_column_name = "meanPAU.dosage_pattern_simplified")

# Assign dosage patterns based on mean PSI values
psi_m323k <- assign_dosage_pattern(psi_m323k,
                      column_name_prefix = "PSImean",
                      wt_key = "WT",
                      het_key = "HET",
                      hom_key = "HOM",
                      separator = ".",
                      pattern_column_name = "meanPSI.dosage_pattern",
                      pattern_simplified_column_name = "meanPSI.dosage_pattern_simplified")


psi_f210i <- assign_dosage_pattern(psi_f210i,
                                   column_name_prefix = "PSImean",
                                   wt_key = "WT",
                                   het_key = "HET",
                                   hom_key = "HOM",
                                   separator = ".",
                                   pattern_column_name = "meanPSI.dosage_pattern",
                                   pattern_simplified_column_name = "meanPSI.dosage_pattern_simplified")

psi_m323k_ad6m <- assign_dosage_pattern(psi_m323k_ad6m,
                                   column_name_prefix = "PSImean",
                                   wt_key = "WT",
                                   het_key = "HET",
                                   hom_key = "HOM",
                                   separator = ".",
                                   pattern_column_name = "meanPSI.dosage_pattern",
                                   pattern_simplified_column_name = "meanPSI.dosage_pattern_simplified")

# Assign dosage patterns based on median PSI values
psi_m323k <- assign_dosage_pattern(psi_m323k,
                                   column_name_prefix = "PSImedian",
                                   wt_key = "WT",
                                   het_key = "HET",
                                   hom_key = "HOM",
                                   separator = ".",
                                   pattern_column_name = "medianPSI.dosage_pattern",
                                   pattern_simplified_column_name = "medianPSI.dosage_pattern_simplified")

psi_f210i <- assign_dosage_pattern(psi_f210i,
                                   column_name_prefix = "PSImedian",
                                   wt_key = "WT",
                                   het_key = "HET",
                                   hom_key = "HOM",
                                   separator = ".",
                                   pattern_column_name = "medianPSI.dosage_pattern",
                                   pattern_simplified_column_name = "medianPSI.dosage_pattern_simplified")

psi_m323k_ad6m <- assign_dosage_pattern(psi_m323k_ad6m,
                                   column_name_prefix = "PSImedian",
                                   wt_key = "WT",
                                   het_key = "HET",
                                   hom_key = "HOM",
                                   separator = ".",
                                   pattern_column_name = "medianPSI.dosage_pattern",
                                   pattern_simplified_column_name = "medianPSI.dosage_pattern_simplified")

# add label for whether dosage patterns are consistent between using mean and median PSI values
psi_m323k <- psi_m323k %>%
  mutate(meanmedianPSI.dosage_pattern_consistency = if_else(meanPSI.dosage_pattern == medianPSI.dosage_pattern, "consistent", "not_consistent"))

psi_f210i <- psi_f210i %>%
  mutate(PSImeanmedian.dosage_pattern_consistency = if_else(meanPSI.dosage_pattern == medianPSI.dosage_pattern, "consistent", "not_consistent"))

psi_m323k_ad6m <- psi_m323k_ad6m %>%
  mutate(meanmedianPSI.dosage_pattern_consistency = if_else(meanPSI.dosage_pattern == medianPSI.dosage_pattern, "consistent", "not_consistent"))

# wtite to TSV
write_tsv(pau_f210i, "processed/2024-02-29/F210I_MEFs/F210I_MEFs.mean_pau.dosage_labels.tsv")
write_tsv(pau_m323k, "processed/2024-02-29/M323K_MEFs/M323K_MEFs.mean_pau.dosage_labels.tsv")
write_tsv(pau_m323k_ad6m, "processed/2024-02-29/M323K_Adult_Brain/M323K_Adult_Brain.mean_pau.dosage_labels.tsv")
write_tsv(psi_f210i, "processed/2024-02-29/F210I_MEFs/F210I_MEFs.psi.mean_median.dosage_labels.tsv")
write_tsv(psi_m323k, "processed/2024-02-29/M323K_MEFs/M323K_MEFs.psi.mean_median.dosage_labels.tsv")
write_tsv(psi_m323k_ad6m, "processed/2024-02-29/M323K_Adult_Brain/M323K_Adult_Brain.psi.mean_median.dosage_labels.tsv")


