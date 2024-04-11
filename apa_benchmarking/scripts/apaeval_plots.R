library(tidyverse)
library(tidytext)


#' Faceted Heatmap of summarised benchmarking metrics across datasets
#'
#' This function generates a faceted heatmap plot based on the provided dataframe and specifications.
#' 
#' @param df dataframe containing APAeval benchmarking results for individual samples.
#' @param metrics A character vector specifying the metrics to include in the plot.
#' @param metric_col The name of the column in \code{df} containing the metric names. Default is "metric".
#' @param value_col The name of the column in \code{df} containing the median values of the metrics (i.e. the value to control fill intensity and rank tools). Default is "metrics.value.median".
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


# identification challenge results
id_df <- read_tsv("data/metrics_01.tsv")
# absolute (TPM) quantification challenge results
absquant0_df <- read_tsv("data/metrics_02_0tpm.tsv") 
absquant1_df <- read_tsv("data/metrics_02_1tpm.tsv")

# remove unnecessary cols, extract dataset name, species and datatype
id_df <- clean_apaeval_df(id_df)
absquant0_df <- clean_apaeval_df(absquant0_df)
absquant1_df <- clean_apaeval_df(absquant1_df)

# for each participant, dataset, window_size + metric, calculate the median + IQR 
id_df_summ <- summarise_metrics(id_df)
absquant0_df_summ  <- summarise_metrics(absquant0_df)
absquant1_df_summ  <- summarise_metrics(absquant1_df)

# Where metrics overlap, visualise de-novo and annotation-based tools together
common_metrics_idquant <- intersect(unique(id_df$metric), unique(absquant0_df$metric))

# Combined df for common metrics of de-novo and annot-based tools, removing APAtrap*1 which is evaluated in both challenges
id_absquant0_df_summ <- bind_rows(filter(id_df_summ, metric %in% common_metrics_idquant),
                                  filter(absquant0_df_summ, metric %in% common_metrics_idquant & str_starts(participant_id, "APAtrap", negate = T))
                                  )
id_absquant1_df_summ <- bind_rows(filter(id_df_summ, metric %in% common_metrics_idquant),
                                  filter(absquant1_df_summ, metric %in% common_metrics_idquant & str_starts(participant_id, "APAtrap", negate = T))
)


# Vector of different window sizes for matching PAS
window_sizes <- unique(id_df_summ$window_size) %>%
  .[!is.na(.)] %>%
  set_names()


# ID tools for precision, sensitivity and F1 score (1 for each window size)
id_heatmaps <- map(window_sizes,
                   ~ id_df_summ %>%
                                   filter(window_size == .x) %>%
                                   plot_faceted_heatmap(.,
                                                        c("Precision", "Sensitivity", "F1_score"),
                                                        plot_ncol = 3,plot_base_size = 14
                                   )
)

# ID tools + PAQR & QAPA on one plot for precision, sensitivity and F1 score (1 for each window size)
id_absquant0_comb_heatmaps <- map(window_sizes,
    ~ id_absquant0_df_summ %>%
      filter(window_size == .x) %>%
      plot_faceted_heatmap(.,
                           c("Precision", "Sensitivity", "F1_score"),
                           plot_ncol = 3
                           )
    )

id_absquant1_comb_heatmaps <- map(window_sizes,
                                  ~ id_absquant1_df_summ %>%
                                    filter(window_size == .x) %>%
                                    plot_faceted_heatmap(.,
                                                         c("Precision", "Sensitivity", "F1_score"),
                                                         plot_ncol = 3
                                    )
)


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
