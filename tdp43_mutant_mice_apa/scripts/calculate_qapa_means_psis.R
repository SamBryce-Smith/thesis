
suppressPackageStartupMessages(library(optparse))
 
option_list <- list(make_option(c("-s", "--sample-table"),
                                type="character",
                                dest="sample_table",
                                help="Path to sample table CSV file used as input to QAPA_snakemake pipeline"),
                    make_option(c("-c", "--contrast-table"),
                                type="character",
                                dest="contrast_table",
                                help="Path to contrasts CSV file used as input to QAPA_snakemake pipeline"),
                    make_option(c("-i", "--PAU_table"),
                                dest = "PAU_table",
                                type = "character",
                                help = "Path to QAPA PAU results table (output of qapa quant call)"),
                    make_option(c("-o", "--output-prefix"),
                                dest = "output_prefix",
                                default = "pau_results_long",
                                help = "Output prefix for formatted QAPA results TSVs containing condition-wise mean PAU & TPM values (<output_prefix>.mean_pau.tsv) and condition-wise mean LABRAT psi values (<output_prefix>.mean_psi.tsv) ([default= %default])"),
                     make_option(c("-u", "--utils"),
                                dest = "utils",
                                default = "calculate_utils.R",
                                help = "path to R script containing required utility functions ([default= %default])")
)

opt_parser <- OptionParser(option_list = option_list)

# if no arguments are provided, print the help messages and exit the script
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
 print_help(opt_parser)
 stop()
}

opt <- parse_args(opt_parser)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
source(opt$utils)

#### - ANALYSIS

sample_df <- read_csv(opt$sample_table, show_col_types=FALSE)
contrasts_df <- read_csv(opt$contrast_table, show_col_types=FALSE)
pau_df <- read_tsv(opt$PAU_table, show_col_types=FALSE)

# Calculate condition-wise mean TPM and PAU for each APA_ID in qapa quant output TSV, converting to long format
wider_cond_mean <- convert_qapa_results(pau_df, sample_df)

# Calculate differences in mean PAU between conditions according to pre-specified contrasts

# Get a list of vectors of c(contrast_key, base_key) (i.e. identifiers of experimental condition) from the contrasts table
delta_cols <- apply(contrasts_df,
                    MARGIN = 1, # apply function to rows not columns
                    function(x) paste("mean_PAU_", c(x["contrast_key"], x["base_key"]), sep = ""), # function applied to each row to get the column names
                    simplify = F # ensure output is a list
) 

# Calculate differences in mean PAU between specified contrasts
wider_cond_mean <- delta_cols %>%
  map(~ 
        reframe(wider_cond_mean,
                "delta_{str_remove_all(.x[1],  'mean_PAU_')}_{str_remove_all(.x[2],  'mean_PAU_')}" := !!sym(.x[1]) - !!sym(.x[2]))
  ) %>%
  # combine the delta columns together into a single dataframe
  bind_cols() %>%
  # combine the delta columns with the original table
  bind_cols(wider_cond_mean, .)


# Calculate LABRAT psi values for each gene

# Create pas_number col for calculating psi values
combined_df <- pau_df %>% 
  select(Gene, Chr, APA_ID, Strand, LastExon.Start, LastExon.End) %>% # include the minimal cols needed for the pas_number function
  right_join(wider_cond_mean, by = "APA_ID") %>% # join it to the df we use
  add_polya_site_numbers() # apply the pas_number function to create the pas_number column


# calculate LABRAT psi values for each gene (based on mean PAU values)
modified_example_df <-  combined_df %>% 
  group_by(Gene_Name) %>% 
  arrange(Pas_Number, .by_group = TRUE) %>% 
  # summarise applies function to cols starting w/ "mean_TPM_" (using across for multiple cols), .names = "psi_{.col}" gives it the psi_ to new col name
  summarise(across(starts_with("mean_TPM_"), labrat_psi_vec, .names = "psi_{str_remove_all(.col, 'mean_TPM_')}")) 

# Get a list of vectors of c(psi_contrast, psi_base) - expected column names containing psi values for each contrast defined in the contrasts table
psi_cols <- apply(contrasts_df,
                    MARGIN = 1, # apply function to rows not columns
                    function(x) paste("psi_", c(x["contrast_key"], x["base_key"]), sep = ""), # function applied to each row to get the column names
                    simplify = F # ensure output is a list
) 

# Calculate differences in psi_values between specified contrasts
modified_example_df <- psi_cols %>%
  map(~ # like a for loop, in 1 line of code
        reframe(modified_example_df,
                "delta_psi_{str_remove_all(.x[1],  'psi_')}_{str_remove_all(.x[2],  'psi_')}" := !!sym(.x[1]) - !!sym(.x[2]))
  ) %>%
  # combine the delta columns for different contrasts together into a single dataframe
  bind_cols() %>%
  # combine the delta columns with the original table
  bind_cols(modified_example_df, .)

# write the final tables to a CSV/TSV file (using the output file name/prefix option we define at the start of the script)
write_tsv(wider_cond_mean, file = paste(opt$output_prefix, ".mean_pau.tsv", sep = ""))
write_tsv(modified_example_df, file = paste(opt$output_prefix, ".mean_psi.tsv", sep = ""))



#Rscript OptParseFunction.R -s M323K_adult_brain\m323k_adult_brain_sample_tbl.csv -c M323K_adult_brain\contrasts_m323k_adult_brain.csv -i M323K_adult_brain\all_samples.pau_results.txt -o M323K_adult_brain\m323k_adult_output

# Rscript OptParseFunction.R -s MEF\MEF_F210I\f210i_mef_sample_tbl.csv -c MEF\MEF_F210I\contrasts_f210i_mefs.csv -i MEF\MEF_F210I\f210i_mefs_all_samples.pau_results.txt -o MEF\MEF_F210I\mef_f210i_output

# Rscript OptParseFunction.R -s MEF\MEF_M323K\m323k_mef_sample_tbl.csv -c MEF\MEF_M323K\contrasts_m323k_mefs.csv -i MEF\MEF_M323K\mef_m323k_all_samples.pau_results.txt -o MEF\MEF_M323K\mef_m323k_output

# Rscript OptParseFunction.R -s F210I_embryonic_brain\f210i_embryonic_brain_long_sample_tbl.csv -c F210I_embryonic_brain\contrasts_f210i_embryonic_brain.csv -i F210I_embryonic_brain\f210i_embryonic_all_samples.pau_results.txt -o F210I_embryonic_brain\f210i_embryonic_output

# Rscript OptParseFunction.R -s M323K_3month_Adult_Brain\m323k_adult_brain_3m_sample_tbl.csv -c M323K_3month_Adult_Brain\contrasts_m323k_adult_brain_3m.csv -i M323K_3month_Adult_Brain\m323k_adult_brain_3m_all_samples.pau_results.txt -o M323K_3month_Adult_Brain\m323k_3m_adult_output

# Rscript OptParseFunction.R -s Q331K_3month_Adult_Brain\q331k_adult_brain_3m_sample_tbl.csv -c Q331K_3month_Adult_Brain\contrasts_q331k_adult_brain_3m.csv -i Q331K_3month_Adult_Brain\q331k_adult_brain_3m_all_samples.pau_results.txt -o Q331K_3month_Adult_Brain\q331k_3m_adult_output


