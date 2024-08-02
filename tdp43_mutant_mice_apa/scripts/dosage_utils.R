library(tidyverse)

#' Prepare dataframe for base_dosage_plot
#'
#' Prepare a wide-format dataframe for generating dosage dependent PAS usage plots (with base_dosage_plot)
#'
#' @param df input dataframe
#' @param col_prefix The prefix of the columns containing per-condition usages.
#' @param names_col The name of the output column containing condition names.
#' @param values_col The name of the output column containing usage values.
#' @param pattern_col The name of the column containing dosage pattern categories.
#' @param dosage_levels A vector specifying the dosage levels/experimental condition keys and their low-high dosage order.
#' @param dosage_patterns A vector specifying the desired plotting order of dosage pattern labels
#' @return A dataframe ready for use with base_dosage_plot
prep_dosage_plot_df <- function(df, col_prefix, names_col, values_col, pattern_col, dosage_levels = c("WT", "HET", "HOM"), 
                                dosage_patterns = c("Up+Up", "Down+Down", "Up+Same", "Down+Same", 
                                                    "Same+Down", "Same+Up", "Up+Down", "Down+Up", "Other", "NA")) {
  df_processed <- df %>%
    pivot_longer(starts_with(col_prefix), names_to = names_col, values_to = values_col, names_prefix = col_prefix) %>%
    mutate("{names_col}" := factor(!!sym(names_col), levels = dosage_levels),
           "{pattern_col}" := factor(!!sym(pattern_col), levels = dosage_patterns)
    )
  
  return(df_processed)
}


#' Generate minimal dosage dependent usage line plot
#'
#' Generate a base ggplot object for visualizing dosage-dependent PAS usage patterns
#'
#' @param data A dataframe containing the data to be plotted (e.g. output of prep_dosage_plot_df. Assumes that have 1 row per condition (x_var), and is ordered factor in lowest-highest dosage
#' @param x_var The name of the variable to be plotted on the x-axis.
#' @param y_var The name of the variable to be plotted on the y-axis.
#' @param group_var The name of the variable to be used for grouping data.
#' @param facet_var The name of the variable to be used for facet wrapping.
#' @param alpha_val The transparency of lines (0 to 1).
#' @param base_font_size The base font size for the plot.
#' @return A ggplot object.
base_dosage_plot <- function(data, x_var, y_var, group_var, facet_var, alpha_val, base_font_size) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], group = .data[[group_var]])) +
    facet_wrap(as.formula(paste0("~", facet_var))) +
    geom_line(alpha = alpha_val) +
    theme_bw(base_size = base_font_size)
}


#' Wrapper function to generate dosage line/spaghetti plots and category count bar plots (and associated dataframes
#'
#' @param df Data frame containing the input data.
#' @param patterns_filter Filter condition for the dosage patterns.
#' @param line_labs ggplot2::labs object containing title and subtitle for the dosage plot (set to labs() to remove blank space).
#' @param bar_labs ggplot2::labs object containing title and subtitle for the bar plot (set to labs() to remove blank space).
#' @param col_prefix Prefix for the column names containing PAU values to be used in the prep_dosage_plot_df function.
#' @param names_col Column name for the genotypes/conditions in the prep_dosage_plot_df function (assumes valid values are WT, HET and HOM).
#' @param values_col Column name for the PAU/usage values in the df.
#' @param pattern_col Column name for the dosage pattern in the df.
#' @param group_var Column name of the group variable for the dosage plot (across genotype i.e. APA_ID).
#' @param alpha_val Alpha value for the base dosage plot.
#' @param base_font_size Base font size for the plots.
#' 
#' @return Named list containing the plot_df, base_dosage_plot, counts_df, and bar_plot.
#' 
dosage_analysis_wrapper <- function(df,
                                    pattern_col = "meanPAU.dosage_pattern",
                                    patterns_filter = !is.na(!!sym(pattern_col)),
                                    line_labs = labs(title = "", subtitle = ""),
                                    bar_labs = labs(title = "", subtitle = ""),
                                    col_prefix = "mean_PAU_",
                                    names_col = "condition",
                                    values_col = "meanPAU",
                                    group_var = "APA_ID",
                                    alpha_val = 0.25,
                                    base_font_size = 14
) {
  
  # Step 1: Prepare the dosage plot data frame
  plot_df <- prep_dosage_plot_df(df,
                                 col_prefix = col_prefix,
                                 names_col = names_col,
                                 values_col = values_col,
                                 pattern_col = pattern_col) %>%
    filter({{ patterns_filter }})
  
  # Step 2: Generate the base dosage plot
  dosage_plot <- base_dosage_plot(
    data = plot_df,
    x_var = names_col,
    y_var = values_col,
    group_var = group_var,
    facet_var = pattern_col,
    alpha_val = alpha_val,
    base_font_size = base_font_size
  ) +
    labs(x = "Genotype", y = "mean PAS usage (%)") +
    line_labs
  
  # Step 3: Calculate summary counts
  counts_df <- plot_df %>%
    distinct(APA_ID, .keep_all = TRUE) %>%
    count(meanPAU.dosage_pattern, sort = TRUE) %>%
    mutate(perc = (n / sum(n)) * 100)
  
  # Step 4: Generate the bar plot
  bar_plot <- counts_df %>%
    ggplot(aes(x = meanPAU.dosage_pattern, y = perc, label = n)) +
    geom_col() +
    geom_text(nudge_y = 5) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    theme_bw(base_size = base_font_size) +
    labs(x = "PAS Usage Pattern (WT - HET - HOM)", y = "Event Frequency (%)") +
    bar_labs
  
  # Return the results as a named list
  return(list(plot_df = plot_df,
              dosage_plot = dosage_plot,
              counts_df = counts_df,
              bar_plot = bar_plot))
}
