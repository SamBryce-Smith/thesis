library(tidyverse)

jaccard_all <- read_tsv("processed/gt_consistency/chr_all.jaccard_index.tsv")
jaccard_te <- read_tsv("processed/gt_consistency/chr_te.jaccard_index.tsv")

unique(jaccard_all$group_id)

# add 'broad' group level identifier (e.g. GTEx, Mayr) - same as used in paper
jaccard_all <- jaccard_all %>%
  mutate(group_id_simple = str_split_i(group_id, "_", 1)) %>%
  relocate(group_id_simple, .after = group_id)

jaccard_te <- jaccard_te %>%
  mutate(group_id_simple = str_split_i(group_id, "_", 1)) %>%
  relocate(group_id_simple, .after = group_id)

# transform to long format (1 row per window size for each comparison)
jaccard_all <- jaccard_all %>%
  pivot_longer(cols = starts_with("jaccardindex"),names_to = c("metric","window_size"), names_sep = "_",
               values_to = "value",
               names_transform = list(window_size = as.integer))

jaccard_te <- jaccard_te %>%
  pivot_longer(cols = starts_with("jaccardindex"),names_to = c("metric","window_size"), names_sep = "_",
               values_to = "value",
               names_transform = list(window_size = as.integer))


jaccard_comb <- bind_rows(all = jaccard_all, te = jaccard_te, .id = "pas_subset")

# scatter plot using all both all and TE, aggregating into larger groups where diff tissues/conds
jaccard_comb %>%
  filter(window_size == 0) %>%
  ggplot(aes(x = group_id_simple, y = value, colour = pas_subset)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_dodge2(width = 0.75),size = 0.75) +
  labs(x = "Dataset",
       y = "Jaccard index",
       colour = "PAS Subset") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
  

# Some datasets have multiple 'conditions' - plot those individually
jaccard_comb %>%
  filter(window_size == 0) %>%
  filter(group_id_simple %in% c("GTEXsim", "HEK293", "Mayr")) %>%
  ggplot(aes(x = group_id, y = value, colour = pas_subset)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_dodge2(width = 0.75),size = 0.75) +
  labs(x = "Dataset",
       y = "Jaccard index",
       colour = "PAS Subset") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")


# Is the Mayr performance affected by window size? only assessed 0 so far
jaccard_comb %>%
  filter(group_id_simple == "Mayr") %>%
  mutate(window_size = factor(window_size, levels = c(0,10,25,50))) %>%
  ggplot(aes( x = window_size, y = value)) +
  facet_wrap("~ group_id", ncol = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()


jaccard_comb %>%
  mutate(window_size = factor(window_size, levels = c(0,10,25,50))) %>%
  ggplot(aes( x = window_size, y = value)) +
  facet_wrap("~ group_id", ncol = 2) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

# Key points so far:
## Jaccard indices slightly higher for last-exon overlapping PAS than complete PASome
## Datasets all have median > 0.6, Mayr exhibits the greatest variability
## Calculation was performed 'within condition' for each dataset
## Mayr dataset variability is driven by a single sub-group ('Mayr_NB'). Within that sub-group, approx half comparisons are ~ 0.3, with remaining similar score to other datasets
## Increasing window size has very minor effect on jaccard index (slight increase in most cases)


#  Should do naive cluster analysis, concerned that jaccard index metric is normalised to length of interval in some way

