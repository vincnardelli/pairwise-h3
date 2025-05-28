library(ggplot2)
library(dplyr)
library(readr) 
library(tidyr) 
library(forcats)


# First proposal ----

# --- Define file paths and parameters ---
# Ensure this directory points to the location of your RDS files
input_data_dir <- "output/ConfigDigest-c439487f_run_20250527_115718/data/"
output_plot_dir_v1 <- "plots_from_analysis_R_v1_english/" # Output for these modified plots

# Create the output directory if it doesn't exist
if (!dir.exists(output_plot_dir_v1)) {
  dir.create(output_plot_dir_v1, recursive = TRUE)
}

summary_stats_file <- file.path(input_data_dir, "summary_stats_approxDGP_20250527_115718.rds")

# --- Load data ---
if (!file.exists(summary_stats_file)) {
  stop(paste("Summary statistics file not found:", summary_stats_file))
}
summary_df <- read_rds(summary_stats_file)

# --- Data Preparation ---
# Convert columns to factors for ggplot2 and define order
summary_df <- summary_df %>%
  mutate(
    N_Scenario = factor(TOTAL_POINTS_TO_GENERATE_scenario,
                        levels = sort(unique(TOTAL_POINTS_TO_GENERATE_scenario)),
                        labels = paste0("N = ", sort(unique(TOTAL_POINTS_TO_GENERATE_scenario)))),
    Lambda_True_Factor = factor(LAMBDA_TRUE_scenario,
                                levels = sort(unique(LAMBDA_TRUE_scenario)),
                                labels = paste0("True Lambda = ", sort(unique(LAMBDA_TRUE_scenario)))),
    DATA_TYPE_scenario = factor(DATA_TYPE_scenario),
    # MODIFIED: Relative Computational Time Efficiency (GM Time / PW Time)
    # Values < 1 mean GM is faster
    # Values > 1 mean GM is slower
    Relative_Time_GM_vs_PW = Mean_Time_GM / Mean_Time_PW
  )

# Filter for scenarios where GM was run
summary_df <- summary_df %>%
  filter(RUN_GM_MODEL_scenario == TRUE)

cat("Glimpse of the prepared summary_df for V1 plots:\n")
glimpse(summary_df)

# --- Plotting Function (Modified for English and new Time Ratio) ---
create_efficiency_plot_v1 <- function(data, y_var, y_label, plot_title, filename_suffix) {
  
  if (!y_var %in% names(data)) {
    warning(paste("Y variable '", y_var, "' not found in data. Skipping plot:", plot_title))
    return(NULL)
  }
  if (nrow(filter(data, !is.na(.data[[y_var]]))) == 0) {
    warning(paste("No non-NA data for Y variable '", y_var, "'. Skipping plot:", plot_title))
    return(NULL)
  }
  
  plot <- ggplot(data, aes(x = LAMBDA_TRUE_scenario, y = .data[[y_var]],
                           color = DATA_TYPE_scenario, shape = DATA_TYPE_scenario,
                           group = DATA_TYPE_scenario)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    facet_wrap(~ N_Scenario, scales = "free_x", labeller = labeller(N_Scenario = label_value)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = unique(data$LAMBDA_TRUE_scenario)) +
    labs(
      title = plot_title,
      x = "True Lambda (DGP Spatial Autocorrelation)",
      y = y_label,
      color = "Data Type",
      shape = "Data Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  if (y_var == "Relative_Time_GM_vs_PW") { # Apply log scale for time ratio
    if(any(data[[y_var]] <= 0, na.rm = TRUE) && any(!is.na(data[[y_var]]))) {
      plot <- plot + scale_y_continuous() 
      cat(paste0("Warning: Non-positive values found in '", y_var, "' for plot '", plot_title, "'. Using linear Y scale.\n"))
    } else if (any(!is.na(data[[y_var]]))) {
      plot <- plot + scale_y_log10()
    }
  } else {
    plot <- plot + scale_y_continuous()
  }
  
  plot_filename <- file.path(output_plot_dir_v1, paste0("RelEff_", filename_suffix, ".pdf"))
  
  print(plot)
  ggsave(plot_filename, plot, width = 12, height = 7, device="pdf", bg="white")
  cat(paste0("Plot saved to: ", plot_filename, "\n"))
  
  return(plot)
}

# --- Generate Plots for Parameter Estimate Efficiencies (V1) ---
# Note: RE = MSE_GM / MSE_PW. If RE < 1, GM is better. If RE > 1, PW is better.

# 1. Relative Efficiency for Beta Slope
plot_beta_v1 <- create_efficiency_plot_v1(
  data = summary_df,
  y_var = "RE_Beta_Slope_PW_vs_GM",
  y_label = "Relative Efficiency (MSE_GM / MSE_PW)",
  plot_title = "Relative Efficiency of Beta1 (Slope) Estimate",
  filename_suffix = "BetaSlope_v1"
)

# 2. Relative Efficiency for Spatial Correlation Parameter (Psi vs Lambda)
plot_spatial_corr_v1 <- create_efficiency_plot_v1(
  data = summary_df,
  y_var = "RE_SpatialCorr_PW_vs_GM",
  y_label = "Relative Efficiency (MSE_GM / MSE_PW)",
  plot_title = "Relative Efficiency of Spatial Parameter (Lambda/Psi) Estimate",
  filename_suffix = "SpatialCorr_v1"
)

# 3. Relative Efficiency for Sigma Squared (Error Variance)
plot_sigma_sq_v1 <- create_efficiency_plot_v1(
  data = summary_df,
  y_var = "RE_SigmaSq_PW_vs_GM",
  y_label = "Relative Efficiency (MSE_GM / MSE_PW)",
  plot_title = "Relative Efficiency of Sigma^2 Epsilon Estimate",
  filename_suffix = "SigmaSq_v1"
)

# --- Generate Plot for Relative Computational Time (V1) ---
# Relative_Time_GM_vs_PW = Mean_Time_GM / Mean_Time_PW
# Values < 1: GM faster
# Values > 1: GM slower (PW faster)
plot_comp_time_v1 <- create_efficiency_plot_v1(
  data = summary_df,
  y_var = "Relative_Time_GM_vs_PW",
  y_label = "Relative Time (Time_GM / Time_PW) - Log Scale", # Modified Label
  plot_title = "Computational Time Ratio (GM vs PW)", # Modified Title
  filename_suffix = "CompTime_v1"
)


# Second proposal ----

# --- Data Preparation for Direct Comparison ---
# Select relevant columns and pivot to long format for MSEs
mse_beta_long <- summary_df %>%
  select(scenario_id, N_Scenario, Lambda_True_Factor, LAMBDA_TRUE_scenario, DATA_TYPE_scenario,
         MSE_Beta_Slope_PW, MSE_Beta_Slope_GM) %>%
  pivot_longer(cols = c(MSE_Beta_Slope_PW, MSE_Beta_Slope_GM),
               names_to = "Method_Metric", values_to = "MSE") %>%
  mutate(Method = ifelse(grepl("_PW$", Method_Metric), "PW", "GM")) %>%
  select(-Method_Metric)

mse_lambda_long <- summary_df %>%
  select(scenario_id, N_Scenario, Lambda_True_Factor, LAMBDA_TRUE_scenario, DATA_TYPE_scenario,
         MSE_Psi_PW, MSE_Lambda_GM) %>% # Note: MSE_Psi_PW for PW
  rename(MSE_Lambda_PW = MSE_Psi_PW) %>% # Rename for consistency in pivot
  pivot_longer(cols = c(MSE_Lambda_PW, MSE_Lambda_GM),
               names_to = "Method_Metric", values_to = "MSE") %>%
  mutate(Method = ifelse(grepl("_PW$", Method_Metric), "PW", "GM")) %>%
  select(-Method_Metric)

mse_sigmasq_long <- summary_df %>%
  select(scenario_id, N_Scenario, Lambda_True_Factor, LAMBDA_TRUE_scenario, DATA_TYPE_scenario,
         MSE_SigmaSq_PW, MSE_SigmaSq_GM) %>%
  pivot_longer(cols = c(MSE_SigmaSq_PW, MSE_SigmaSq_GM),
               names_to = "Method_Metric", values_to = "MSE") %>%
  mutate(Method = ifelse(grepl("_PW$", Method_Metric), "PW", "GM")) %>%
  select(-Method_Metric)

# Select relevant columns and pivot to long format for Time
time_long <- summary_df %>%
  select(scenario_id, N_Scenario, Lambda_True_Factor, LAMBDA_TRUE_scenario, DATA_TYPE_scenario,
         Mean_Time_PW, Mean_Time_GM) %>%
  pivot_longer(cols = c(Mean_Time_PW, Mean_Time_GM),
               names_to = "Method_Metric", values_to = "Time") %>%
  mutate(Method = ifelse(grepl("_PW$", Method_Metric), "PW", "GM")) %>%
  select(-Method_Metric)

cat("\nExample of long format data for Beta MSE (V2 plots):\n")
glimpse(mse_beta_long)

# --- Plotting Function for Direct Comparison ---
create_direct_comparison_plot <- function(data_long, y_var, y_label, plot_title, filename_suffix, use_log_scale_y = TRUE) {
  
  if (!y_var %in% names(data_long)) {
    warning(paste("Y variable '", y_var, "' not found in data_long. Skipping plot:", plot_title))
    return(NULL)
  }
  if (nrow(filter(data_long, !is.na(.data[[y_var]]))) == 0) {
    warning(paste("No non-NA data for Y variable '", y_var, "'. Skipping plot:", plot_title))
    return(NULL)
  }
  
  plot <- ggplot(data_long, aes(x = LAMBDA_TRUE_scenario, y = .data[[y_var]],
                                color = Method, shape = Method, linetype = Method,
                                group = Method)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_line(linewidth = 0.7, alpha = 0.7) +
    facet_grid(N_Scenario ~ DATA_TYPE_scenario, labeller = labeller(N_Scenario = label_value, DATA_TYPE_scenario = label_value)) +
    scale_x_continuous(breaks = unique(data_long$LAMBDA_TRUE_scenario)) +
    labs(
      title = plot_title,
      x = "True Lambda (DGP Spatial Autocorrelation)",
      y = y_label,
      color = "Method",
      shape = "Method",
      linetype = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size=8),
      axis.title.y = element_text(size=10),
      legend.position = "top",
      strip.text = element_text(face = "bold", size=9),
      panel.spacing = unit(0.8, "lines")
    )
  
  if (use_log_scale_y) {
    if(any(data_long[[y_var]] <= 0, na.rm = TRUE) && any(!is.na(data_long[[y_var]]))) {
      plot <- plot + scale_y_continuous() 
      cat(paste0("Warning: Non-positive values found in '", y_var, "' for plot '", plot_title, "'. Using linear Y scale.\n"))
    } else if (any(!is.na(data_long[[y_var]]))) {
      plot <- plot + scale_y_log10()
    }
  } else {
    plot <- plot + scale_y_continuous()
  }
  
  plot_filename <- file.path(output_plot_dir_v1, paste0("DirectComp_", filename_suffix, ".pdf"))
  
  print(plot)
  ggsave(plot_filename, plot, width = 11, height = 9, device="pdf", bg="white") # Adjusted size for grid
  cat(paste0("Plot saved to: ", plot_filename, "\n"))
  
  return(plot)
}

# --- Generate Direct Comparison Plots (V2) ---

# 1. Beta1 MSE Comparison
plot_beta_mse_v2 <- create_direct_comparison_plot(
  data_long = mse_beta_long,
  y_var = "MSE",
  y_label = "Mean Squared Error (MSE) - Log Scale",
  plot_title = "Comparison of Beta1 (Slope) MSE: PW vs GM",
  filename_suffix = "BetaMSE_v2"
)

# 2. Spatial Parameter MSE Comparison
plot_lambda_mse_v2 <- create_direct_comparison_plot(
  data_long = mse_lambda_long,
  y_var = "MSE",
  y_label = "Mean Squared Error (MSE) - Log Scale",
  plot_title = "Comparison of Spatial Parameter (Lambda/Psi) MSE: PW vs GM",
  filename_suffix = "LambdaMSE_v2"
)

# 3. SigmaSq Epsilon MSE Comparison
plot_sigmasq_mse_v2 <- create_direct_comparison_plot(
  data_long = mse_sigmasq_long,
  y_var = "MSE",
  y_label = "Mean Squared Error (MSE) - Log Scale",
  plot_title = "Comparison of Sigma^2 Epsilon MSE: PW vs GM",
  filename_suffix = "SigmaSqMSE_v2"
)

# 4. Computation Time Comparison
plot_time_v2 <- create_direct_comparison_plot(
  data_long = time_long,
  y_var = "Time",
  y_label = "Mean Computation Time (seconds) - Log Scale",
  plot_title = "Comparison of Computation Time: PW vs GM",
  filename_suffix = "CompTime_v2"
)

