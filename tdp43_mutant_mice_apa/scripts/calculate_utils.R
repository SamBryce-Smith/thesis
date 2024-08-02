suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

#' Calculate QAPA results to calculate mean polyA site expression and usage
#'
#' This function takes a wide-format data frame of QAPA results and a sample table, and calculates the mean TPM (Transcripts Per Million) and PAU (PolyA Usage) for each APA_ID and condition. The results are returned in a long-format data frame.
#'
#' @param qapa_df A data frame in wide format corresponding to QAPA output results table (output of qapa quant)
#' @param sample_table A data frame corresponding to sample table used to configure QAPA_snakemake. Must include sample_name column with values corresponding to sample names input to QAPA, and a condition column to group samples by experimental conditionsample information, including sample names and conditions
#'
#' @return A data frame in long format containing the mean TPM and PAU for each APA_ID and condition. The data frame also includes the Gene_Name for each APA_ID.
#'
convert_qapa_results <- function(qapa_df, sample_table) {
  
  # Converting PAU columns to long format
  qapa_pau_long <- qapa_df %>%
    select(-ends_with(".TPM")) %>% 
    filter(!str_ends(APA_ID, "_S")) %>% 
    pivot_longer(cols = ends_with(".PAU"),
                 names_to = "sample_name",
                 values_to = "PAU",
                 names_pattern = "(.*).PAU")
  
  qapa_pau_long_filt <- qapa_pau_long %>%
    group_by(APA_ID) %>% 
    filter(!all(PAU == 0)) %>% 
    ungroup()
  
  # Converting TPM columns to long format
  qapa_tpm_long <- qapa_df %>%
    select(-ends_with(".PAU")) %>%
    filter(!str_ends(APA_ID, "_S")) %>% 
    pivot_longer(cols = ends_with(".TPM"),
                 names_to = "sample_name",
                 values_to = "TPM",
                 names_pattern = "(.*).TPM") %>%
    group_by(APA_ID) %>% 
    filter(!all(TPM == 0)) %>% 
    ungroup()
  
  # Joining the long tables together
  qapa_comb_long <- left_join(x = qapa_pau_long_filt,
                              y = select(qapa_tpm_long, APA_ID, sample_name, TPM),
                              by = c("APA_ID", "sample_name"))
  
  # Calculate mean TPM for each APA_ID and condition, and PAU
  TPM_mean <- left_join(qapa_comb_long, sample_table, by = "sample_name") %>% 
    group_by(APA_ID, condition) %>% 
    summarise(mean_TPM = mean(TPM, na.rm = TRUE),
              mean_PAU = mean(PAU, na.rm = TRUE)) %>% 
    ungroup()
              
               
  
  # Pivot wider to get data in the desired format
  wider_cond_mean <- pivot_wider(TPM_mean,
                                 id_cols = APA_ID,
                                 names_from = condition,
                                 values_from = c(mean_TPM, mean_PAU)) 
  
  # Joining with Gene_Name
  gene_name <- qapa_comb_long %>%
    select("APA_ID", "Gene_Name") %>%
    distinct()
  
  all_data <- left_join(wider_cond_mean, gene_name, by = "APA_ID")  
  
  
  # Adding mean_TPM to the final result
  all_data_plot <- left_join(all_data, TPM_mean, by = "APA_ID") %>%
    relocate(Gene_Name, .after = APA_ID)
  
  return(all_data_plot)
}




#' Add strand-aware polyA site numbers/index to a data frame
#'
#' This function adds polyA site numbers and total polyA sites for a given group of intervals/genes, where 1 = most proximal site and n = most distal site for given gene (and total number of polyA sites)
#'
#' @param df Input data frame/tibble of polyA sites.
#' @param site_outcol The name of the output column for the polyA site numbers. Default is "Pas_Number".
#' @param total_outcol The name of the output column for the total number of polyA sites in the given group. Default is "Total_Pas".
#' @param region_id_col A character vector specifying the columns that identify a group. Default is c("Gene", "Chr").
#' @param event_id_col The name of the column that contains unique polyA site identifiers. Default is "APA_ID".
#' @param strand_col The name of the column that indicates the strand. Default is "Strand".
#' @param start_col The name of the column that indicates the start coordinate of the interval. Default is "LastExon.Start".
#' @param end_col The name of the column that indicates the end coordinate of the interval. Default is "LastExon.End".
#'
#' @return The input data frame with two additional columns: one for the strand-aware polyA site numbers/indices and one for the total number of polyA sites.
#'
#' @examples
#' df <- data.frame(
#'   Gene = c("Gene1", "Gene1", "Gene2", "Gene2"),
#'   Chr = c("Chr1", "Chr1", "Chr2", "Chr2"),
#'   APA_ID = c("APA1", "APA2", "APA3", "APA4"),
#'   Strand = c("+", "+", "-", "-"),
#'   LastExon.Start = c(100, 200, 300, 400),
#'   LastExon.End = c(150, 250, 350, 450)
#' )
#' add_polya_site_numbers(df)
add_polya_site_numbers = function(df,
                                  site_outcol = "Pas_Number",
                                  total_outcol = "Total_Pas",
                                  region_id_col = c("Gene", "Chr"),
                                  event_id_col = "APA_ID",
                                  strand_col = "Strand",
                                  start_col = "LastExon.Start",
                                  end_col = "LastExon.End") {
  
  df %>%
    dplyr::group_by(across(all_of(region_id_col))) %>%
    # group_vars()
    #if + strand then End coord = 3'end - smallest value = left-most = most proximal
    #if - strand then Start coord = 5'end - largest value = left-most = most proximal
    dplyr::mutate("{site_outcol}" := dplyr::if_else(.data[[strand_col]] == '+',
                                                    true = row_number(.data[[end_col]]), # smallest value assigned 1
                                                    false = row_number(desc(.data[[start_col]])) # largest value assigned 1
    ),
    "{total_outcol}" := n_distinct(.data[[event_id_col]])
    ) %>%
    dplyr::ungroup()
  
}

#' Calculate LABRAT's psi metric from proximal-distal ordered vector of polyA site TPM (Transcripts Per Million) values
#'
#' @param tpms A numeric vector of TPM values for each polyA site in a gene. The vector should be ordered such that the first element corresponds to the first/most proximal polyA site, the second element to the second polyA site ...n polyA sites
#'
#' @return The PSI value for the gene, calculated as the sum of scaled TPMs divided by the sum of unscaled TPMs. The scaling factor for each TPM is its position in the gene (0-indexed) divided by the total number of polyA sites minus 1.
#'
#' @examples
#' tpms <- c(100, 200, 300, 400, 500)
#' labrat_psi_vec(tpms)
#' [1] 0.6666667
labrat_psi_vec <- function(tpms) {
  
  # calculate scaling factors - zero-based idx
  m <- seq(0, length(tpms) - 1) # m = scaling factor, position starting from 0 (pos in gene)
  
  n <- length(tpms)
  
  # calculate scaled TPMs
  # unscaled_tpm * (m / (n -1 ))
  scaled_tpms <- tpms * (m / (n - 1))
  
  # calculate psi - scaled TPM / unscaled TPM 
  sum(scaled_tpms) / sum(tpms)
  
}
