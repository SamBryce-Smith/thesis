library(tidyverse)
library(tidytext)
library(ggbump)
library(patchwork)

#' Faceted Heatmap of summarised benchmarking metrics across datasets
#'
#' This function generates a faceted heatmap plot based on the provided dataframe and specifications.
#' 
#' @param df dataframe containing APAeval benchmarking results for individual samples.
#' @param metrics A character vector specifying the metrics to include in the plot.
#' @param metric_col The name of the column in \code{df} containing the metric names. Default is "metric".
#' @param value_col The name of the column in \code{df} containing the median values of the metrics (i.e. the value to control fill intensity and rank tool order on the y axis). Default is "metrics.value.median".
#' @param value_range_col The name of the column in \code{df} containing the interquartile range (IQR) values of the metrics (to be added to text in each cell of heatmap) Default is "metrics.value.iqr".
#' @param dataset_col The name of the column in \code{df} containing dataset names. Default is "dataset".
#' @param participant_col The name of the column in \code{df} containing the participating tool names. Default is "participant_id".
#' @param return_plot_df Logical indicating whether to return the processed dataframe used for plotting purposes instead of the plot. Default is \code{FALSE}.
#' @param plot_base_size The base font size for the plot. Default is 14.
#' @param label_size The size value for text label in heatmap cells (passed to geom_text). Default is 3.
#' @param plot_labs Labels for the plot axes and legend (in ggplot2 labs format). Default is \code{labs(x = "Dataset", y = "Participant", fill = "Median across datasets")}.
#' @param plot_legend_position Position of the legend in the plot. Default is "bottom".
#' @param plot_nrow Number of rows in the facet grid. Default is 1.
#' @param plot_ncol Number of columns in the facet grid. Default is 2.
#' 
#' @return a ggplot2 object containing a metric-faceted heatmap across datasets and tools. If \code{return_plot_df} is \code{TRUE}, returns the processed dataframe used for plotting.
#' 
#' @export
plot_faceted_heatmap <- function(df,
                                 metrics,
                                 metric_col = "metric",
                                 value_col = "metrics.value.median",
                                 value_range_col = "metrics.value.iqr",
                                 dataset_col = "dataset",
                                 participant_col = "participant_id",
                                 return_plot_df = FALSE,
                                 plot_base_size = 14,
                                 label_size = 3,
                                 plot_labs =   labs(x = "Dataset",
                                                    y = "Participant",
                                                    fill = "Median across datasets"),
                                 plot_legend_position = "bottom",
                                 plot_nrow = 1,
                                 plot_ncol = 2) {
  
  plot_df <- df %>%
    filter(!!sym(metric_col) %in% metrics) %>%
    mutate(metrics.value.label = paste(round(!!sym(value_col), 2), "\n+/-", round(!!sym(value_range_col), 2), sep = " "),
           {{metric_col}} := as.factor(!!sym(metric_col)),
           {{participant_col}} := reorder_within(!!sym(participant_col), !!sym(value_col), !!sym(metric_col), fun = median)) 
  
  if (return_plot_df) {
    return(plot_df)
  }
  
  plot_df %>%
    ggplot(aes(x = !!sym(dataset_col), y = !!sym(participant_col), fill = !!sym(value_col), label = metrics.value.label)) +
    facet_wrap(paste("~", metric_col), scales = "free_y", nrow = plot_nrow, ncol = plot_ncol) +
    scale_y_reordered() +
    geom_tile() +
    geom_text(size = 3) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_bw(base_size = plot_base_size) +
    plot_labs +
    theme(legend.position = plot_legend_position,
          axis.text.x = element_text(angle = 90),)
  
}

#' Faceted Heatmap of summarised benchmarking metrics across datasets (sorting by performance in a reference dataset)
#'
#' This function generates a faceted heatmap plot based on the provided dataframe and specifications. Is basically plot_faceted_heatmap but sorts by single dataset, not the median across all datasets
#' 
#' @param df dataframe containing APAeval benchmarking results for individual samples.
#' @param metrics A character vector specifying the metrics to include in the plot.
#' @param metrics A single character vector/string specifying the reference dataset for which to sort participants from 1..n.
#' @param metric_col The name of the column in \code{df} containing the metric names. Default is "metric".
#' @param value_col The name of the column in \code{df} containing the median values of the metrics (i.e. the value to control fill intensity and rank tool order on the y axis). Default is "metrics.value.median".
#' @param value_range_col The name of the column in \code{df} containing the interquartile range (IQR) values of the metrics (to be added to text in each cell of heatmap) Default is "metrics.value.iqr".
#' @param dataset_col The name of the column in \code{df} containing dataset names. Default is "dataset".
#' @param participant_col The name of the column in \code{df} containing the participating tool names. Default is "participant_id".
#' @param return_plot_df Logical indicating whether to return the processed dataframe used for plotting purposes instead of the plot. Default is \code{FALSE}.
#' @param plot_base_size The base font size for the plot. Default is 14.
#' @param label_size The size value for text label in heatmap cells (passed to geom_text). Default is 3.
#' @param plot_label_sep The string to put between the value in value_col and value_range_col in heatmap cells. Default is "+/-" (can do e.g. "\n+/-" to put on separate lines)
#' @param plot_labs Labels for the plot axes and legend (in ggplot2 labs format). Default is \code{labs(x = "Dataset", y = "Participant", fill = "Median across datasets")}.
#' @param plot_legend_position Position of the legend in the plot. Default is "bottom".
#' @param plot_nrow Number of rows in the facet grid. Default is 1.
#' @param plot_ncol Number of columns in the facet grid. Default is 2.
#' 
#' @return a ggplot2 object containing a metric-faceted heatmap across datasets and tools. If \code{return_plot_df} is \code{TRUE}, returns the processed dataframe used for plotting.
#' 
#' @export
plot_faceted_heatmap_ref <- function(df,
                                     metrics,
                                     ref_dataset = "AllExperimental",
                                     metric_col = "metric",
                                     value_col = "metrics.value.median",
                                     value_range_col = "metrics.value.iqr",
                                     dataset_col = "dataset",
                                     participant_col = "participant_id",
                                     plot_base_size = 14,
                                     label_size = 3,
                                     plot_label_sep = "+/-", # '\n+/-' to put on separate line
                                     plot_labs = labs(x = "Dataset",
                                                      y = "Participant",
                                                      fill = "Median across samples"),
                                     plot_legend_position = "bottom",
                                     plot_nrow = 2,
                                     plot_ncol = 2) {
  
  
  value_col_select <- c(value_col, value_range_col)
  group_cols <- c(participant_col, metric_col)
  
  # Subset to reference dataset
  df_ref <- df %>%
    filter(dataset == ref_dataset) %>%
    rename_with(.fn = ~ paste0("reference.", .x, recycle0 = TRUE),
                .cols = all_of(value_col_select)) %>%
    select(all_of(group_cols), starts_with("reference."))
  
  
  # Make a dummy column for ordering - contains the reference values for ranking for just single dataset (experimental summary)
  df <- left_join(df, df_ref, by = group_cols) 

  ref_value_col <- paste0("reference.", value_col)
  
  # Prepare df for plotting
  plot_df <- df %>%
    filter(metric %in% metrics) %>%
    mutate(metrics.value.label = paste(round(!!sym(value_col), 2), plot_label_sep, round(!!sym(value_range_col), 2), sep = " "),
           {{metric_col}} := factor(!!sym(metric_col), levels = metrics),
           {{participant_col}} := reorder_within(!!sym(participant_col), !!sym(ref_value_col), !!sym(metric_col), fun = median)
    ) 
  
  
  plot_df %>%
    ggplot(aes(x = !!sym(dataset_col),
               y = !!sym(participant_col),
               fill = !!sym(value_col),
               label = metrics.value.label)
           ) +
    facet_wrap(paste("~", metric_col), scales = "free_y", nrow = plot_nrow, ncol = plot_ncol) +
    scale_y_reordered() +
    geom_tile() +
    geom_text(size = label_size) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_bw(base_size = plot_base_size) +
    plot_labs +
    theme(legend.position = plot_legend_position,
          axis.text.x = element_text(angle = 90)
          )
}


clean_apaeval_df <- function(df, rm_cols = c("...1"), id_col = "challenge_id") {
  
  #
  df <- select(df, -all_of(rm_cols))
  
  df %>%
  mutate(dataset = str_split_i(!!sym(id_col), "_", 1),
         species = str_extract(!!sym(id_col), "hg38|mm10"),
         datatype = if_else(dataset == "GTEXsim", "simulated", "experimental")) 
  
}

#' Calculate median + inter-quartile ranges across datasets and metrics
#'
#' 
#' @param df The dataframe containing the APAeval benchmarking event outputs data to be summarised.
#' @param group_cols A character vector specifying the columns (& order) to group by. Default is c("participant_id", "dataset", "metric").
#' @param value_col The name of the column in \code{df} containing the metric values. Default is "metrics.value".
#' 
#' @return A tibble summarising the median and IQR of the metric values for each group (i.e. each tool, dataset and metric)
#' 
#' @export
summarise_metrics <- function(df,
                              group_cols = c("participant_id", "dataset", "window_size", "metric"),
                              value_col = "metrics.value") {
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(metrics.value.median = median(!!sym(value_col)),
              metrics.value.iqr = IQR(!!sym(value_col))) %>%
    ungroup()
}


#' Rank Metrics values (1..n in highest-lowest) within Data Frame Groups
#'
#' This function ranks rows within groups of a data frame based on a specified column. Mainly intended to be applied to the output of summarise_metrics
#'
#' @param df A data frame containing metrics to be ranked.
#' @param group_cols Character vector specifying the columns used for grouping the data frame. Defaults to c("window_size", "metric", "dataset")
#' @param rank_by A string specifying the column name to use for ranking within each group. Defaults to "metrics.value.median"
#'   This should be a column containing the metric values used for ranking.
#' @param rank_outcol A string specifying the name of the new column to store the ranks.
#'   Defaults to "rank".
#'
#' @return A modified data frame with an additional column containing the ranks within each group.
#'
#' @import dplyr
#'
#' @examples
#' # Example usage:
#' data <- data.frame(
#'   window_size = c(10, 20, 10, 20, 10, 20),
#'   metric = c("A", "A", "B", "B", "C", "C"),
#'   dataset = c("X", "X", "X", "X", "Y", "Y"),
#'   value = c(5, 10, 3, 8, 15, 20)
#' )
#'
#' # Rank by 'value' within groups defined by 'window_size', 'metric', and 'dataset':
#' ranked_data <- rank_metrics(data, group_cols = c("window_size", "metric", "dataset"), rank_by = "value")
#'
#' @export
rank_metrics <- function(df,
                         group_cols = c("window_size", "metric", "dataset"),
                         rank_by = "metrics.value.median",
                         rank_outcol = "rank") {
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    arrange(desc(!!sym(rank_by))) %>%
    mutate("{rank_outcol}":= min_rank(desc(!!sym(rank_by)))) %>%
    ungroup()
}


#' Make a bump plot comparing ranks in experimental data as a function of window size
plot_rank_vs_window <- function(input_df, metric_str) {
  # Filter the dataframe based on the specified metric
  tmp_df <- input_df %>%
    filter(metric == metric_str)
  
  # Define window sizes for plotting
  plot_window_sizes_vec <- c(10, 25, 50, 100) %>%
    set_names(c(25, 50, 75, 100))
  
  plot_window_sizes_df <- enframe(plot_window_sizes_vec, "plot_window_size", "window_size") %>%
    mutate(plot_window_size = as.integer(plot_window_size))
  
  # Join window sizes with filtered dataframe
  tmp_df <- tmp_df %>%
    left_join(plot_window_sizes_df)
  
  # Plotting
  plot <- tmp_df %>%
    ggplot(aes(x = plot_window_size, y = rank.experimental, group = participant_id)) +
    geom_point() +
    geom_bump() +
    scale_y_reverse(breaks = seq(0, 12)) +
    geom_text(data = tmp_df %>%
                filter(window_size == max(window_size)),
              aes(label = participant_id), nudge_x = 50) +
    scale_x_continuous(breaks = c(25, 50, 75, 100),
                       limits = c(25, 200),
                       labels = plot_window_sizes_vec) +
    labs(title = str_to_title(str_replace(metric_str, "_", " ")),
         x = "Window size (nt)",
         y = "Rank (median of experimental samples)") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          axis.line = element_blank())
  
  return(plot)
}




# identification challenge results
id_df <- read_tsv("data/metrics_01.tsv")
# absolute (TPM) quantification challenge results
absquant0_df <- read_tsv("data/metrics_02_0tpm.tsv") 
absquant1_df <- read_tsv("data/metrics_02_1tpm.tsv")

# remove unnecessary cols, extract dataset name, species and datatype
id_df <- clean_apaeval_df(id_df)
absquant0_df <- clean_apaeval_df(absquant0_df)
absquant1_df <- clean_apaeval_df(absquant1_df)

# for each participant, dataset, window_size + metric, calculate the median + IQR. Rank participants according to median within each dataset
id_df_summ <- summarise_metrics(id_df) %>%
  rank_metrics()
absquant0_df_summ  <- summarise_metrics(absquant0_df) %>%
  rank_metrics()
absquant1_df_summ  <- summarise_metrics(absquant1_df) %>%
  rank_metrics()

# subset to experimental datasets only - calculate global median + IQR
id_df_summ_exp <- id_df %>%
  filter(datatype == "experimental") %>%
  summarise_metrics(group_cols = c("participant_id", "window_size", "metric")) %>%
  rank_metrics(group_cols = c("window_size", "metric"))

absquant0_df_summ_exp  <- absquant0_df %>%
  filter(datatype == "experimental") %>%
  summarise_metrics(group_cols = c("participant_id", "window_size", "metric")) %>%
  rank_metrics(group_cols = c("window_size", "metric"))

absquant1_df_summ_exp  <- absquant1_df %>%
  filter(datatype == "experimental") %>%
  summarise_metrics(group_cols = c("participant_id", "window_size", "metric")) %>%
  rank_metrics(group_cols = c("window_size", "metric"))

# For globally ranked metrics, merge dfs together
# i.e. single df with rank.experimental, rank.simulated as columns
# also calcuate a difference in rank

id_df_summ_comb <- full_join(filter(id_df_summ, dataset == "GTEXsim") %>% select(-dataset),
          id_df_summ_exp,
          by = c("participant_id", "window_size", "metric"),
          suffix = c(".simulated", ".experimental"),
          relationship = "one-to-one"
          ) %>%
  relocate(rank.simulated, .before = rank.experimental) %>%
  # negative = poorer relative performance in experimental rel to simulated
  mutate(rank_diff_sim2exp = rank.simulated - rank.experimental)

absquant0_df_summ_comb <- full_join(filter(absquant0_df_summ, dataset == "GTEXsim") %>% select(-dataset),
          absquant0_df_summ_exp,
          by = c("participant_id", "window_size", "metric"),
          suffix = c(".simulated", ".experimental"),
          relationship = "one-to-one"
) %>%
  relocate(rank.simulated, .before = rank.experimental)

absquant1_df_summ_comb <- full_join(filter(absquant1_df_summ, dataset == "GTEXsim") %>% select(-dataset),
          absquant1_df_summ_exp,
          by = c("participant_id", "window_size", "metric"),
          suffix = c(".simulated", ".experimental"),
          relationship = "one-to-one"
) %>%
  relocate(rank.simulated, .before = rank.experimental)


# Where metrics overlap, visualise de-novo and annotation-based tools together
common_metrics_idquant <- intersect(unique(id_df$metric), unique(absquant0_df$metric))

# Combined df for common metrics of de-novo and annot-based tools, removing APAtrap*1 which is evaluated in both challenges
# Add ranks considering combo of de-novo and reference-based tools
id_absquant0_df_summ <- bind_rows(filter(id_df_summ, metric %in% common_metrics_idquant),
                                  filter(absquant0_df_summ, metric %in% common_metrics_idquant & str_starts(participant_id, "APAtrap", negate = T))
                                  ) %>%
  rank_metrics(rank_outcol = "rank.combined")

id_absquant1_df_summ <- bind_rows(filter(id_df_summ, metric %in% common_metrics_idquant),
                                  filter(absquant1_df_summ, metric %in% common_metrics_idquant & str_starts(participant_id, "APAtrap", negate = T))
                                  ) %>%
  rank_metrics(rank_outcol = "rank.combined")


# Repeat combining with median across all experimental data sets

id_absquant0_df_summ_exp <- bind_rows(filter(id_df_summ_exp, metric %in% common_metrics_idquant),
                                  filter(absquant0_df_summ_exp, metric %in% common_metrics_idquant & str_starts(participant_id, "APAtrap", negate = T))
) %>%
  rank_metrics(group_cols = c("window_size", "metric"),
               rank_outcol = "rank.combined")

id_absquant1_df_summ_exp <- bind_rows(filter(id_df_summ_exp, metric %in% common_metrics_idquant),
                                  filter(absquant1_df_summ_exp, metric %in% common_metrics_idquant & str_starts(participant_id, "APAtrap", negate = T))
) %>%
  rank_metrics(group_cols = c("window_size", "metric"),
               rank_outcol = "rank.combined")



#### PLOTS


## How does window size influence the ranking of tool performance?

# bump plot, summarizing across all experimental (+ simulated) - window size on x axis

# Make rank plots across window sizes for selected metrics
window_size_metrics <- c("Precision", "Sensitivity", "F1_score", "Jaccard_index")

window_size_exp_rank_bump_plot <- c("Precision", "Sensitivity", "F1_score", "Jaccard_index") %>%
  map(~ plot_rank_vs_window(id_df_summ_comb, .x)) %>%
  wrap_plots(ncol = 2, axis_titles = "collect")


window_size_exp_rank_bump_plot



# bump plot comparing simulated vs experimental ranks for a fixed window size across window sizes

# Prepare df for plotting exper vs simulated
exper_vs_sim_rank_plot_df <- id_df_summ_comb %>%
  filter(metric %in% window_size_metrics) %>%
  pivot_longer(cols = starts_with("rank."), names_to = "datatype", values_to = "rank",names_prefix = "rank.") %>%
  mutate(datatype = factor(str_to_sentence(datatype), levels = c("Simulated", "Experimental")),
         plot_label = if_else(datatype == "Experimental", participant_id, ""),
         plot_metric = factor(str_to_title(str_replace_all(metric, "_", " ")),
                              levels = str_to_title(str_replace_all(window_size_metrics, "_", " ")))
  )


# Vector of different window sizes for matching PAS
window_sizes <- unique(id_df_summ$window_size) %>%
  .[!is.na(.)] %>%
  set_names()

# Make faceted plot for 4 metrics across window sizes
exper_vs_simm_rank_plots <- map2(.x = window_sizes,
     .y = names(window_sizes),
     ~ exper_vs_sim_rank_plot_df %>%
       filter(window_size == .x) %>%
       ggplot(aes(x = datatype, y = rank, group = participant_id, label = plot_label)) +
       facet_wrap("~ plot_metric", ncol = 2) +
       geom_point() +
       geom_line() +
       scale_y_reverse(breaks = seq(1,12,1)) +
       geom_text(nudge_x = 0.25) +
       labs(subtitle = paste0("Window size = ", .y),
            x = "Dataset Type",
            y = "Rank (by median across samples)") +
       theme_bw(base_size = 14) +
       theme(legend.position = "none",
             panel.grid = element_blank()
       )
)

# Heatmaps of tool median + IQR for each dataset


# Add medians + ranks from combining all experimental datasets
# first, just for de-novo/annotated independently
id_df_summ <- id_df_summ_exp %>%
  mutate(dataset = "AllExperimental") %>%
  bind_rows(., id_df_summ) %>%
  arrange(window_size, dataset, metric, rank)

absquant0_df_summ <- absquant0_df_summ_exp %>%
  mutate(dataset = "AllExperimental") %>%
  bind_rows(., absquant0_df_summ) %>%
  arrange(window_size, dataset, metric, rank)

absquant1_df_summ <- absquant1_df_summ_exp %>%
  mutate(dataset = "AllExperimental") %>%
  bind_rows(., absquant1_df_summ) %>%
  arrange(window_size, dataset, metric, rank)

# now for combo of de-novo + annotated
id_absquant0_df_summ <- id_absquant0_df_summ_exp %>%
  mutate(dataset = "AllExperimental") %>%
  bind_rows(id_absquant0_df_summ, .) %>%
  arrange(window_size, dataset, metric, rank)

id_absquant1_df_summ <- id_absquant1_df_summ_exp %>%
  mutate(dataset = "AllExperimental") %>%
  bind_rows(id_absquant1_df_summ, .) %>%
  arrange(window_size, dataset, metric, rank)


# Across all window sizes, make heatmap for prec, sens, f1 and jaccard (1 plot per window size) (ranking by median of all experimental samples)
id_heatmaps_ref <- map(window_sizes,
                       ~ plot_faceted_heatmap_ref(filter(id_df_summ, window_size == .x),
                                                  metrics = window_size_metrics,
                                                  ref_dataset = "AllExperimental")
                       )

absquant0_heatmaps_ref <- map(window_sizes,
                       ~ plot_faceted_heatmap_ref(filter(absquant0_df_summ, window_size == .x),
                                                  metrics = window_size_metrics,
                                                  ref_dataset = "AllExperimental")
)

absquant1_heatmaps_ref <- map(window_sizes,
                       ~ plot_faceted_heatmap_ref(filter(absquant1_df_summ, window_size == .x),
                                                  metrics = window_size_metrics,
                                                  ref_dataset = "AllExperimental")
)


# repeat for de-novo + reference-based combined
id_absquant0_heatmaps_ref <- map(window_sizes,
                                 ~ plot_faceted_heatmap_ref(filter(id_absquant0_df_summ, window_size == .x),
                                                            metrics = window_size_metrics,
                                                            ref_dataset = "AllExperimental")
                                 )


id_absquant1_heatmaps_ref <- map(window_sizes,
                             ~ plot_faceted_heatmap_ref(filter(id_absquant1_df_summ, window_size == .x),
                                                        metrics = window_size_metrics,
                                                        ref_dataset = "AllExperimental")
)


# id_heatmaps_ref$`50`
# absquant0_heatmaps_ref$`50`
# absquant1_heatmaps_ref$`50`
# id_absquant0_heatmaps_ref$`50`
# id_absquant1_heatmaps_ref$`50`


# # ID tools for precision, sensitivity and F1 score (1 for each window size)
# id_heatmaps <- map(window_sizes,
#                    ~ id_df_summ %>%
#                                    filter(window_size == .x) %>%
#                                    plot_faceted_heatmap(.,
#                                                         c("Precision", "Sensitivity", "F1_score"),
#                                                         plot_ncol = 3,plot_base_size = 14
#                                    )
# )

# ID tools + PAQR & QAPA on one plot for precision, sensitivity and F1 score (1 for each window size)
# id_absquant0_comb_heatmaps <- map(window_sizes,
#     ~ id_absquant0_df_summ %>%
#       filter(window_size == .x) %>%
#       plot_faceted_heatmap(.,
#                            c("Precision", "Sensitivity", "F1_score"),
#                            plot_ncol = 3
#                            )
#     )
# 
# id_absquant1_comb_heatmaps <- map(window_sizes,
#                                   ~ id_absquant1_df_summ %>%
#                                     filter(window_size == .x) %>%
#                                     plot_faceted_heatmap(.,
#                                                          c("Precision", "Sensitivity", "F1_score"),
#                                                          plot_ncol = 3
#                                     )
# )


if (!dir.exists("processed")) {dir.create("processed")}

walk2(id_absquant0_comb_heatmaps,
      names(id_absquant0_comb_heatmaps),
      ~ ggsave(filename = file.path("processed", paste0("2024-04-11_", "denovo_annot_id_heatmap.prec_sens_f1.window_size_", .y, ".pdf")),
               plot = .x,
               width = 11.7,
               height = 8.3, units = "in")
      )


walk2(id_heatmaps,
      names(id_heatmaps),
      ~ ggsave(filename = file.path("processed", paste0("2024-04-09_", "denovo_id_heatmap.prec_sens_f1.window_size_", .y, ".pdf")),
               plot = .x,
               width = 11.7,
               height = 8.3, units = "in")
)


# Output summary TSVs with median and ranks
# by dataset
# summarised by median across experimental
# Combined with reference-based tools

# id_absquant0_df_summ 
# id_absquant1_df_summ
# id_absquant0_df_summ_exp
# id_absquant1_df_summ_exp

write_tsv(id_absquant0_df_summ, file.path("processed", paste0("2024-04-12_identification_denovo_reference",
                                        ".mintpm_0",
                                        ".metric_medians",
                                        ".per_dataset",
                                        ".tsv"))
          )


write_tsv(id_absquant1_df_summ, file.path("processed", paste0("2024-04-12_identification_denovo_reference",
                                                              ".mintpm_1",
                                                              ".metric_medians",
                                                              ".per_dataset",
                                                              ".tsv"))
          )


write_tsv(id_absquant0_df_summ_exp, file.path("processed", paste0("2024-04-12_identification_denovo_reference",
                                                              ".mintpm_0",
                                                              ".metric_medians",
                                                              ".experimental_global",
                                                              ".tsv"))
)

write_tsv(id_absquant1_df_summ_exp, file.path("processed", paste0("2024-04-12_identification_denovo_reference",
                                                                  ".mintpm_1",
                                                                  ".metric_medians",
                                                                  ".experimental_global",
                                                                  ".tsv"))
)


# Above dfs only contain 'matching' metrics i.e. F1-score, precision and sensitivity 
# Want summaries for abs quant (experimental only medians, per experiment medians)

write_tsv(absquant0_df_summ_exp, file.path("processed", paste0("2024-04-12_absolute_quantification",
                                                               ".mintpm_0",
                                                               ".metric_medians",
                                                               ".experimental_global",
                                                               ".tsv"))
          )

write_tsv(absquant0_df_summ, file.path("processed", paste0("2024-04-12_absolute_quantification",
                                                 ".mintpm_0",
                                                 ".metric_medians",
                                                 ".per_dataset",
                                                 ".tsv"))
          )
       

write_tsv(absquant1_df_summ_exp, file.path("processed", paste0("2024-04-12_absolute_quantification",
                                                           ".mintpm_1",
                                                           ".metric_medians",
                                                           ".experimental_global",
                                                           ".tsv"))
)   


write_tsv(absquant1_df_summ, file.path("processed", paste0("2024-04-12_absolute_quantification",
                                                           ".mintpm_1",
                                                           ".metric_medians",
                                                           ".per_dataset",
                                                           ".tsv"))
)

