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
                                default = "psi_results",
                                help = "Output prefix for formatted QAPA results matrices: per-sample LABRAT psi values (<output_prefix>.psi.per_sample.tsv), condition-wise mean/median LABRAT psi values and deltas (<output_prefix>.psi.mean_median.tsv), PSI matrix of genes containing NA values (<output_prefix>.psi.per_sample_na.tsv), per-sample gene-level TPM matrix (<output_prefix>.gene_tpm.per_sample.tsv), and condition-(un)wise mean/median gene-level TPM matrix (<output_prefix>.gene_tpm.mean_median.tsv) ([default= %default])"),
                    make_option(c("-u", "--utils"),
                                dest = "utils",
                                default = "calculate_utils.R",
                                help = "path to R script containing required utility functions ([default= %default])")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Calculate LABRAT psi values from QAPA quant output")

# if no arguments are provided, print the help messages and exit the script
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

opt <- parse_args(opt_parser)





suppressPackageStartupMessages(library(tidyverse))
source(opt$utils)

sample_df <- read_csv(opt$sample_table, show_col_types = FALSE)
contrasts_df <- read_csv(opt$contrast_table, show_col_types = FALSE)
qapa_df <- read_tsv(opt$PAU_table, show_col_types = FALSE)

# qapa_df <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/all_samples.pau_results.txt", show_col_types = FALSE)
# sample_df <- read_csv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/f210i_mef_sample_tbl.csv", show_col_types = FALSE)
# contrasts_df <- read_csv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/contrasts_f210i_mefs.csv", show_col_types = FALSE)
output_prefix <- opt$output_prefix
# "processed/qapa_psis"


message("Getting per-sample TPM matrix for APA genes in long format...")
# Add strand-aware polyA site numbers
qapa_df <- add_polya_site_numbers(qapa_df)

# get df of TPM values in long format
qapa_tpm_long <- qapa_df %>%
  select(-ends_with(".PAU")) %>%
  # remove single PAS genes
  filter(!str_ends(APA_ID, "_S")) %>% 
  pivot_longer(cols = ends_with(".TPM"),
               names_to = "sample_name",
               values_to = "TPM",
               names_pattern = "(.*).TPM")

message("Calculating per-sample LABRAT psi values...")
# calculate per-sample LABRAT psi values for each gene (based on TPM values)
psi_per_sample <- qapa_tpm_long %>%
  group_by(sample_name, Gene, Chr) %>%
  # within each sample/gene, sort APA_IDs from 1..n in terms of 5'-3' rank
  arrange(Pas_Number, .by_group = TRUE) %>%
  summarise(PSI = labrat_psi_vec(TPM)) %>%
  ungroup()

# make a 'wide-format' psi matrix for exporting
psi_per_sample_wide <- pivot_wider(psi_per_sample, names_from = "sample_name", values_from = "PSI", names_glue = "{.name}.PSI")

# Add in sample grouping labels (ready to calculate means/medians and deltas)
psi_per_sample <- left_join(psi_per_sample, select(sample_df, sample_name, condition), by = "sample_name")

# track any genes that have NA psi value across any sample (i.e. underlying TPM has NA value (for unknown reason...))
psi_na_any <- psi_per_sample %>%
  group_by(Gene, Chr) %>%
  filter(any(is.na(PSI))) %>%
  ungroup()

# wider format for export
psi_na_any_wide <- pivot_wider(psi_na_any, names_from = "sample_name", values_from = "PSI", names_glue = "{.name}.PSI")

message("Calculating condition-wise mean, median PSIs and deltas across specified contrasts...")
# calculate different per-condition summary PSI values
psi_summ <- psi_per_sample %>%
  # remove any genes with NA psi values
  group_by(Gene, Chr) %>%
  filter(!any(is.na(PSI))) %>%
  group_by(condition, Gene, Chr) %>%
  summarise(PSImean = mean(PSI),
         PSImedian = median(PSI)
         ) %>%
  ungroup()

# Put condition-wise summary values into columns (wide format)
psi_summ_wide <- pivot_wider(psi_summ,
                            id_cols = all_of(c("Gene", "Chr")),
            names_from = condition,
            values_from = c(PSImean, PSImedian),
            names_sep = ".") 


# Get a list of vectors of c(contrast_key, base_key) (i.e. identifiers of experimental condition) from the contrasts table
delta_cols_mean <- apply(contrasts_df,
                    MARGIN = 1, # apply function to rows not columns
                    function(x) paste("PSImean", c(x["contrast_key"], x["base_key"]), sep = "."), # function applied to each row to get the column names
                    simplify = F # ensure output is a list
) 

delta_cols_median <- apply(contrasts_df,
                         MARGIN = 1, # apply function to rows not columns
                         function(x) paste("PSImedian", c(x["contrast_key"], x["base_key"]), sep = "."), # function applied to each row to get the column names
                         simplify = F # ensure output is a list
                         ) 

# delta_cols_mean

# Calculate differences in mean psi between specified contrasts
psi_summ_wide <- delta_cols_mean %>%
  map(~ 
        reframe(psi_summ_wide,
                "deltaPSImean.{str_remove_all(.x[1],  'PSImean.')}_{str_remove_all(.x[2],  'PSImean.')}" := !!sym(.x[1]) - !!sym(.x[2]))
  ) %>%
  # combine the delta columns together into a single dataframe
  bind_cols() %>%
  # combine the delta columns with the original table
  bind_cols(psi_summ_wide, .)

# repeat with median psi
psi_summ_wide <- delta_cols_median %>%
  map(~ 
        reframe(psi_summ_wide,
                "deltaPSImedian.{str_remove_all(.x[1],  'PSImedian.')}_{str_remove_all(.x[2],  'PSImedian.')}" := !!sym(.x[1]) - !!sym(.x[2]))
  ) %>%
  # combine the delta columns together into a single dataframe
  bind_cols() %>%
  # combine the delta columns with the original table
  bind_cols(psi_summ_wide, .)

message("Calculating gene-level TPM matrix across samples and mean, median condition-(un)wise...")
# Calculate gene-level TPM matrix
gene_tpm_long <- qapa_tpm_long %>%
  group_by(sample_name, Gene, Chr) %>%
  summarise(geneTPM = sum(TPM)) %>%
  ungroup()

# wider format to saving to disk
gene_tpm_wide <- pivot_wider(gene_tpm_long, names_from = sample_name, values_from = geneTPM, names_prefix = "geneTPM.")

# Calculate mean & median geneTPM for each condition
gene_tpm_summ <- qapa_tpm_long %>%
  left_join(select(sample_df, sample_name, condition), by = "sample_name") %>%
  group_by(condition, Gene, Chr) %>%
  summarise(geneTPMmean = mean(TPM),
            geneTPMmedian = median(TPM)) %>%
  ungroup()

# also calculate mean & median across all samples
gene_tpm_summ_all <- qapa_tpm_long %>%
  group_by(Gene, Chr) %>%
  summarise(geneTPMmean.all = mean(TPM),
            geneTPMmedian.all = median(TPM)) %>%
  ungroup()

# pivot to wider mtx for export
gene_tpm_summ_wide <- pivot_wider(gene_tpm_summ,
           id_cols = all_of(c("Gene", "Chr")),
           names_from = condition,
           values_from = c(geneTPMmean, geneTPMmedian),
           names_sep = ".")

gene_tpm_summ_wide <- left_join(gene_tpm_summ_wide, gene_tpm_summ_all, by = c("Gene", "Chr"))

# TO OUTPUT:
message(paste("Outputting matrices to disk with prefix: ", output_prefix, sep = ""))

# gene-level per-sample TPM matrix
write_tsv(gene_tpm_wide, paste(output_prefix, "gene_tpm.per_sample.tsv", sep = "."))
          
# gene-level condition-wise &  mean/median TPM matrix
write_tsv(gene_tpm_summ_wide, paste(output_prefix, "gene_tpm.mean_median.tsv", sep = "."))

# matrix of condition wise mean/median labrat psi + deltas for specified contrasts
write_tsv(psi_summ_wide, paste(output_prefix, "psi.mean_median.tsv", sep = "."))

# per-sample psi matrix
write_tsv(psi_per_sample_wide, paste(output_prefix, "psi.per_sample.tsv", sep = "."))

# psi matrix of genes containing NA values in at least 1 sample
write_tsv(psi_na_any_wide, paste(output_prefix, "psi.per_sample_na.tsv", sep = "."))