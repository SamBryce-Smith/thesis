library(tidyverse)

# supplementary table 1 in TSV form (with some funky formatting)
apaeval_datasets <- read_tsv("data/APAeval_Datasets_summary_and_comparisons_to_run.tsv",
                             col_select = all_of(c('SRA_accession_rnaseq', 'sample_name', 'strandedness', 'mate_layout', 'organism', 'avg_read_length', 'SRA_accession_gt', 'sequencing_method'))
                             )
# GT metadata from SRA explorer
sra_meta <- read_tsv("data/apaeval_gt_sra_explorer_metadata.tsv") 

# 'samples' BEDs downlaoded from polyA sit
polyasite_hs <- read_lines("data/atlas.clusters.2.0.GRCh38.96.tsv.gz",n_max = 1) %>%
  str_split("\t", simplify = T) 

polyasite_hs[1, 1:14] %>% paste(collapse = " ; ")
# [1] "chrom ; chromStart ; chromEnd ; name ; score ; strand ; rep ; frac_samples ; nr_prots ; annotation ; gene_name ; gene_id ; repSite_signals ; GSM1614163|PAPERCLIP|HeLa cells|PAPERCLIP for HeLa cells, replicate 1|none
 
# convert to df, extract first row from matrix (only 1 header line), then remove first 13 default columns
polyasite_hs <- polyasite_hs[1, 14:ncol(polyasite_hs)] %>%
  as_tibble_col("atlas_sample_id") %>%
  # extract accession
  mutate(accession = str_split_i(atlas_sample_id, "\\|", 1))
  

# repeat all in one go for mouse
polyasite_mm <- read_lines("data/atlas.clusters.2.0.GRCm38.96.tsv.gz",n_max = 1) %>%
  str_split("\t", simplify = T) %>%
  .[1, 14:ncol(.)] %>%
  as_tibble_col("atlas_sample_id") %>%
  # extract accession
  mutate(accession = str_split_i(atlas_sample_id, "\\|", 1))
  

# Extract GSM accession where possible, tidy up df for joining with datasets
sra_meta <- sra_meta %>%
  mutate(gsm_accession = str_split_i(Title, ":", 1),
         gsm_accession = if_else(str_starts(gsm_accession, "GSM", negate=T),
                                 NA,
                                 gsm_accession)
         ) %>%
  relocate(gsm_accession, .after = Accession) %>%
  rename(SRA_accession_gt = Accession)


# append gsm accession where extracted
apaeval_datasets_sra <- apaeval_datasets %>%
  left_join(select(sra_meta, SRA_accession_gt, gsm_accession, title = Title))



apaeval_datasets_sra_hs <- filter(apaeval_datasets_sra, organism == "Hsap") 
apaeval_datasets_sra_mm <- filter(apaeval_datasets_sra, organism == "Mmus") 

#
apaeval_datasets_sra_mm <- apaeval_datasets_sra_mm %>%
  left_join(polyasite_mm, by = c("SRA_accession_gt" = "accession"), suffix = c("", ".sra")) %>%
  left_join(polyasite_mm, by = c("gsm_accession" = "accession"), suffix = c("", ".gsm"))


apaeval_datasets_sra_hs <- apaeval_datasets_sra_hs %>%
  left_join(polyasite_hs, by = c("SRA_accession_gt" = "accession"), suffix = c("", ".sra")) %>%
  left_join(polyasite_hs, by = c("gsm_accession" = "accession"), suffix = c("", ".gsm")) 



# srx IDs extracted from PolyASite cannot be used here - output to text file, 
