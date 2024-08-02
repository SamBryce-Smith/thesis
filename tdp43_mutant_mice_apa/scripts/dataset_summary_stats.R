library(tidyverse)

# read in dfs output by rna_seq_single_steps/samtools_stats.smk
dfs <- list.files("data/samtools_stats",
                  pattern = "samtools_stats_sn_summary.tsv",
                  recursive = T,
                  full.names = T) %>%
  set_names(basename(dirname(.))) %>%
  map(~ read_tsv(.x, show_col_types = F,
                 col_select = all_of(c("sample_name", "metric", "value")))
      ) %>%
   bind_rows(.id = "dataset")

# assign genotypes from sample names
dfs <- dfs %>%
  mutate(genotype = case_when(str_detect(sample_name, "HOM|hom") ~ "HOM",
                              str_detect(sample_name, "HET|het") ~ "HET",
                              str_detect(sample_name, "WT|wt") ~ "WT",
                              T ~ "NA"
                              )
         ) 

distinct(dfs, genotype)


# assign mutant based on sample name (all names start with MUT)
dfs <- dfs %>%
  mutate(mutant = str_split_i(sample_name, "_", 1))

unique(dfs$dataset)

# assign simple, cleaned dataset names
dfs <- dfs %>%
  mutate(dataset_clean = case_when(dataset == "m323k_12month_spinal_cord" ~ "Fratta 2018 12 months",
                                   dataset == "m323k_f210i_mef" ~ "Fratta 2018 MEF",
                                   dataset == "m323k_q331k_3m" ~ "Unpublished 3 months"))

# clean metric names (remove the ':', replace spaces)
dfs <- dfs %>%
  mutate(metric = str_remove_all(metric, ":$"),
         metric = str_replace_all(metric, " ", "_"))


unique(dfs$metric)


# define column names - some to average
metrics_to_average <- c("reads_mapped", "reads_mapped_and_paired")
metrics_averaged <- c("average_length")

# across datasets + muts , calculate averages for mapped reads (per million scale)
dfs_to_ave <- dfs %>%
  filter(metric %in% metrics_to_average) %>%
  group_by(dataset_clean, mutant, metric) %>%
  summarise(mean = mean(value) / 1e6,
            median = median(value) / 1e6,
            min = min(value) / 1e6,
            max = max(value) / 1e6) %>%
  ungroup()

# read length is already averaged - calculate mean, median, min and max per dataset
dfs_ave <- dfs %>%
  filter(metric %in% metrics_averaged) %>%
  group_by(dataset_clean, mutant, metric) %>%
  summarise(mean = mean(value),
            median = median(value),
            min = min(value),
            max = max(value)
            ) %>%
  ungroup()

# calculate sample numbers (per genotype)
dfs_n_geno <- dfs %>%
  distinct(dataset_clean, mutant, genotype, sample_name) %>%
  count(dataset_clean, mutant, genotype) %>%
  # sample counts per genotype as columns
  pivot_wider(id_cols = all_of(c("dataset_clean", "mutant")),
              names_from = "genotype",
              values_from = "n",
              names_prefix = "n_")

# put metrics as columns (for easier joining/presentation in final table)
dfs_to_ave <- dfs_to_ave %>%
  pivot_wider(id_cols = all_of(c("dataset_clean", "mutant")),
              names_from = "metric",
              values_from = all_of(c("mean", "median", "min", "max")),
              )

dfs_ave <- dfs_ave %>%
  pivot_wider(id_cols = all_of(c("dataset_clean", "mutant")),
              names_from = "metric",
              values_from = all_of(c("mean", "median", "min", "max")),
  )

dfs_n_geno
dfs_ave
dfs_to_ave

# combine into a single df
dfs_summ_all <- dfs_n_geno %>%
  left_join(dfs_ave, by = c("dataset_clean", "mutant")) %>%
  left_join(dfs_to_ave, by = c("dataset_clean", "mutant"))

dfs_summ_all  

if (!dir.exists("processed/summary_stats")) {dir.create("processed/summary_stats", recursive = T)}

dfs_summ_all

write_tsv(dfs_summ_all,
          file.path("processed/summary_stats",
                    "2024-07-12_samtools_summary_stats.all.tsv"),
          col_names = T)

# only the 'mapped' depth columns
write_tsv(select(dfs_summ_all, -ends_with("_paired")),
          file.path("processed/summary_stats",
                    "2024-07-12_samtools_summary_stats.mapped.tsv"),
          col_names = T)

# only the 'mapped and paired' depth columns
write_tsv(select(dfs_summ_all, -ends_with("_mapped")),
          file.path("processed/summary_stats",
                    "2024-07-12_samtools_summary_stats.mapped.tsv"),
          col_names = T)

