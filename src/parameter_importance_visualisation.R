 # ==============================================================================
# Parameter Importance Visualisation Functions
# 
# Part of: MACROM: An Optimal Control Model for Balancing Climate Change Abatement 
# and Damage Trade-offs
# Authors: Nina Rynne, Michael Bode, Melanie Roberts, Ryan Heneghan
# Institution: Griffith University
# 
# Copyright (c) 2025 Nina Rynne
# Licensed under CC-BY-4.0 - see LICENSE file for details
# 
# Citation: If you use this code, please cite:
#   [Full citation of your paper]
# 
# Version: 1.0.0
# Last updated: January 2026
# ==============================================================================

#' @title Parameter Importance Visualisation Functions
#' @description
#' Functions to create spaghetti plots showing optimal control results across 
#' multiple parameter sets from Latin Hypercube Sampling. All runs are plotted 
#' with transparency to show the range of solutions, with median values highlighted.
#'
#' @note All required libraries (ggplot2, dplyr, tidyr, patchwork, purrr, 
#' colorspace, scales, cowplot, here) must be loaded before using these functions.
#'
#' @author Nina Rynne
#' @date January 2026

# ============================================================================
# Global definitions
# ============================================================================

# Define consistent theme for all parameter importance plots
parameter_importance_theme <- theme_bw() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Define color palette for spaghetti plots
# Each variable gets a distinct color with transparency for individual runs
# and a darker version for the median line
spaghetti_colors <- list(
  emissions = list(
    individual = "firebrick",
    median = colorspace::darken("firebrick", 0.3)
  ),
  temperature = list(
    individual = "steelblue",
    median = colorspace::darken("steelblue", 0.3)
  ),
  mitigation = list(
    individual = "forestgreen",
    median = colorspace::darken("forestgreen", 0.3)
  ),
  cdr = list(
    individual = "purple4",
    median = colorspace::darken("purple4", 0.3)
  ),
  abatement_cost = list(
    individual = "darkorange",
    median = colorspace::darken("darkorange", 0.3)
  ),
  damage_cost = list(
    individual = "maroon",
    median = colorspace::darken("maroon", 0.3)
  )
)

# ============================================================================
# Individual plot functions
# ============================================================================

#' @title Create Cumulative Emissions Spaghetti Plot
#' @description
#' Creates a spaghetti plot showing cumulative CO2 emissions over time
#' for multiple parameter sets, with all runs shown in transparent lines
#' and the median highlighted
#' 
#' @param solution_list Named list of solution objects
#' @param title Plot title (optional, default: "Cumulative Emissions")
#' @return ggplot object
create_emissions_spaghetti_plot <- function(solution_list, 
                                            title = "Cumulative Emissions") {
  
  # Extract cumulative emissions data from all runs
  plot_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    data.frame(
      run_id = run_name,
      years = sol$years,
      cumulative_emissions = sol$cumulative_emissions
    )
  })
  
  # Calculate median across all runs
  summary_data <- plot_data %>%
    group_by(years) %>%
    summarise(median_val = median(cumulative_emissions), .groups = "drop")
  
  # Create spaghetti plot with median overlay
  ggplot() +
    geom_line(data = plot_data, 
              aes(x = years, y = cumulative_emissions, group = run_id),
              color = spaghetti_colors$emissions$individual, 
              alpha = 0.05, 
              linewidth = 0.2) +
    geom_line(data = summary_data, 
              aes(x = years, y = median_val),
              color = spaghetti_colors$emissions$median, 
              linewidth = 0.7) +
    labs(y = "Cumulative Emissions (GtCO2)", title = title) +
    parameter_importance_theme
}

#' @title Create Temperature Anomaly Spaghetti Plot
#' @description
#' Creates a spaghetti plot showing temperature anomaly trajectories over time
#' for multiple parameter sets, with all runs shown in transparent lines
#' and the median highlighted.
#' 
#' @param solution_list Named list of solution objects
#' @param title Plot title (optional, default: "Temperature Anomaly")
#' @return ggplot object
create_temperature_spaghetti_plot <- function(solution_list, 
                                              title = "Temperature Anomaly") {
  
  # Extract temperature anomaly data from all runs
  plot_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    data.frame(
      run_id = run_name,
      years = sol$years,
      temperature_anomaly = sol$temperature_anomaly
    )
  })
  
  # Calculate median across all runs
  summary_data <- plot_data %>%
    group_by(years) %>%
    summarise(median_val = median(temperature_anomaly), .groups = "drop")
  
  # Create spaghetti plot with median overlay
  ggplot() +
    geom_line(data = plot_data, 
              aes(x = years, y = temperature_anomaly, group = run_id),
              color = spaghetti_colors$temperature$individual, 
              alpha = 0.05, 
              linewidth = 0.2) +
    geom_line(data = summary_data, 
              aes(x = years, y = median_val),
              color = spaghetti_colors$temperature$median, 
              linewidth = 0.7) +
    labs(y = "Temperature Anomaly (°C)", title = title) +
    parameter_importance_theme
}

#' @title Create Mitigation Spaghetti Plot
#' @description
#' Creates a spaghetti plot showing annual mitigation quantities over time
#' for multiple parameter sets, with all runs shown in transparent lines
#' and the median highlighted.
#' 
#' @param solution_list Named list of solution objects
#' @param title Plot title (optional, default: "Mitigation")
#' @return ggplot object
create_mitigation_spaghetti_plot <- function(solution_list, 
                                             title = "Mitigation") {
  
  # Extract mitigation data from all runs
  plot_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    data.frame(
      run_id = run_name,
      years = sol$years,
      qty_mitig = sol$qty_mitig
    )
  })
  
  # Calculate median across all runs
  summary_data <- plot_data %>%
    group_by(years) %>%
    summarise(median_val = median(qty_mitig), .groups = "drop")
  
  # Create spaghetti plot with median overlay
  ggplot() +
    geom_line(data = plot_data, 
              aes(x = years, y = qty_mitig, group = run_id),
              color = spaghetti_colors$mitigation$individual, 
              alpha = 0.05, 
              linewidth = 0.2) +
    geom_line(data = summary_data, 
              aes(x = years, y = median_val),
              color = spaghetti_colors$mitigation$median, 
              linewidth = 0.7) +
    labs(y = "Mitigation (GtCO2/yr)", title = title) +
    parameter_importance_theme
}

#' @title Create CDR Spaghetti Plot
#' @description
#' Creates a spaghetti plot showing annual carbon dioxide removal quantities over time
#' for multiple parameter sets, with all runs shown in transparent lines
#' and the median highlighted.
#' 
#' @param solution_list Named list of solution objects
#' @param title Plot title (optional, default: "CDR")
#' @return ggplot object
create_cdr_spaghetti_plot <- function(solution_list, 
                                      title = "CDR") {
  
  # Extract CDR data from all runs
  plot_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    data.frame(
      run_id = run_name,
      years = sol$years,
      qty_remov = sol$qty_remov
    )
  })
  
  # Calculate median across all runs
  summary_data <- plot_data %>%
    group_by(years) %>%
    summarise(median_val = median(qty_remov), .groups = "drop")
  
  # Create spaghetti plot with median overlay
  ggplot() +
    geom_line(data = plot_data, 
              aes(x = years, y = qty_remov, group = run_id),
              color = spaghetti_colors$cdr$individual, 
              alpha = 0.05, 
              linewidth = 0.2) +
    geom_line(data = summary_data, 
              aes(x = years, y = median_val),
              color = spaghetti_colors$cdr$median, 
              linewidth = 0.7) +
    labs(y = "CDR (GtCO2/yr)", title = title) +
    parameter_importance_theme
}

#' @title Create Abatement Cost Spaghetti Plot
#' @description
#' Creates a spaghetti plot showing abatement costs (mitigation + CDR costs) 
#' as a proportion of baseline GWP over time for multiple parameter sets, 
#' with all runs shown in transparent lines and the median highlighted.
#' 
#' @param solution_list Named list of solution objects
#' @param title Plot title (optional, default: "Abatement Costs")
#' @return ggplot object
create_abatement_cost_spaghetti_plot <- function(solution_list, 
                                                 title = "Abatement Costs") {
  
  # Extract and calculate abatement costs normalised by baseline GWP
  plot_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    
    # Calculate total abatement cost as sum of mitigation and CDR costs
    abatement_cost_annual <- sol$mitig_costs_annual + sol$remov_costs_annual
    
    # Normalise by baseline GWP (as proportion)
    abatement_cost_proportion <- abatement_cost_annual / sol$baseline_annual_gwp
    
    data.frame(
      run_id = run_name,
      years = sol$years,
      abatement_cost_proportion = abatement_cost_proportion
    )
  })
  
  # Calculate median across all runs
  summary_data <- plot_data %>%
    group_by(years) %>%
    summarise(median_val = median(abatement_cost_proportion, na.rm = TRUE), .groups = "drop")
  
  # Create spaghetti plot with median overlay
  ggplot() +
    geom_line(data = plot_data, 
              aes(x = years, y = abatement_cost_proportion, group = run_id),
              color = spaghetti_colors$abatement_cost$individual, 
              alpha = 0.05, 
              linewidth = 0.2) +
    geom_line(data = summary_data, 
              aes(x = years, y = median_val),
              color = spaghetti_colors$abatement_cost$median, 
              linewidth = 0.7) +
    scale_y_continuous(labels = percent) +
    labs(y = "Abatement Costs\n(% of GWP)", title = title) +
    parameter_importance_theme
}

#' @title Create Damage Cost Spaghetti Plot
#' @description
#' Creates a spaghetti plot showing temperature-related damage costs 
#' as a proportion of baseline GWP over time for multiple parameter sets, 
#' with all runs shown in transparent lines and the median highlighted.
#' 
#' @param solution_list Named list of solution objects
#' @param title Plot title (optional, default: "Damage Costs")
#' @return ggplot object
create_damage_cost_spaghetti_plot <- function(solution_list, 
                                              title = "Damage Costs") {
  
  # Extract and calculate damage costs normalised by baseline GWP
  plot_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    
    # Use temperature costs as damage costs
    damage_cost_annual <- sol$temp_costs_annual
    
    # Normalise by baseline GWP (as proportion)
    damage_cost_proportion <- damage_cost_annual / sol$baseline_annual_gwp
    
    data.frame(
      run_id = run_name,
      years = sol$years,
      damage_cost_proportion = damage_cost_proportion
    )
  })
  
  # Calculate median across all runs
  summary_data <- plot_data %>%
    group_by(years) %>%
    summarise(median_val = median(damage_cost_proportion, na.rm = TRUE), .groups = "drop")
  
  # Create spaghetti plot with median overlay
  ggplot() +
    geom_line(data = plot_data, 
              aes(x = years, y = damage_cost_proportion, group = run_id),
              color = spaghetti_colors$damage_cost$individual, 
              alpha = 0.05, 
              linewidth = 0.2) +
    geom_line(data = summary_data, 
              aes(x = years, y = median_val),
              color = spaghetti_colors$damage_cost$median, 
              linewidth = 0.7) +
    scale_y_continuous(labels = percent) +
    labs(y = "Damage Costs\n(% of GWP)", title = title) +
    parameter_importance_theme
}

# ============================================================================
# Dashboard function
# ============================================================================

#' @title Create Parameter Importance Dashboard
#' @description
#' Creates a comprehensive dashboard with 6 spaghetti plots showing optimal control 
#' results across multiple parameter sets from Latin Hypercube Sampling. Plots are 
#' arranged in a 3×2 grid with shared y-axis limits for comparable plots.
#'
#' @param importance_results Output from run_parameter_importance(), or just the 
#'   successful_runs component
#' @param save_plot Whether to save the plot as PDF (default: FALSE)
#' @param filename Custom filename for saving (default: auto-generated with timestamp)
#' @param print_insights Whether to print analytical insights after creating dashboard (default: FALSE)
#' @param verbose Print progress information (default: TRUE)
#'
#' @return Combined ggplot object containing the dashboard
#'
#' @examples
#' # Create dashboard from parameter importance results
#' dashboard <- create_parameter_importance_dashboard(
#'   importance_results = results,
#'   save_plot = TRUE,
#'   print_insights = TRUE
#' )
create_parameter_importance_dashboard <- function(importance_results, 
                                                  save_plot = FALSE,
                                                  filename = NULL,
                                                  print_insights = FALSE,
                                                  verbose = TRUE) {
  
  # Extract successful_runs component from full results object
  if ("successful_runs" %in% names(importance_results)) {
    solution_list <- importance_results$successful_runs
  } else {
    # Assume it's already the successful_runs component
    solution_list <- importance_results
  }
  
  # Validate that we have solution data
  if (length(solution_list) == 0) {
    stop("No successful runs found. Ensure run_parameter_importance() completed successfully.")
  }
  
  if (verbose) {
    cat("Creating parameter importance dashboard for", length(solution_list), "runs\n")
  }
  
  # Calculate shared y-axis limits for mitigation and CDR plots
  mitigation_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    if (!is.null(sol)) {
      data.frame(mitigation = sol$qty_mitig)
    }
  })
  
  cdr_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    if (!is.null(sol)) {
      data.frame(cdr = sol$qty_remov)
    }
  })
  
  shared_abatement_limits <- c(0, max(c(mitigation_data$mitigation, cdr_data$cdr), na.rm = TRUE))
  
  # Calculate shared y-axis limits for cost plots
  abatement_cost_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    if (!is.null(sol)) {
      abatement_cost_annual <- sol$mitig_costs_annual + sol$remov_costs_annual
      abatement_cost_proportion <- abatement_cost_annual / sol$baseline_annual_gwp
      data.frame(abatement_cost = abatement_cost_proportion)
    }
  })
  
  damage_cost_data <- map_dfr(names(solution_list), function(run_name) {
    sol <- solution_list[[run_name]]
    if (!is.null(sol)) {
      damage_cost_annual <- sol$temp_costs_annual
      damage_cost_proportion <- damage_cost_annual / sol$baseline_annual_gwp
      data.frame(damage_cost = damage_cost_proportion)
    }
  })
  
  shared_cost_limits <- c(0, max(c(abatement_cost_data$abatement_cost, damage_cost_data$damage_cost), na.rm = TRUE))
  
  # Create individual plots with panel labels
  p1 <- create_emissions_spaghetti_plot(solution_list, title = "a) Cumulative emissions")
  
  p2 <- create_temperature_spaghetti_plot(solution_list, title = "b) Temperature anomaly")
  
  p3 <- create_mitigation_spaghetti_plot(solution_list, title = "c) Mitigation") +
    coord_cartesian(ylim = shared_abatement_limits)
  
  p4 <- create_cdr_spaghetti_plot(solution_list, title = "d) CDR") +
    coord_cartesian(ylim = shared_abatement_limits)
  
  p5 <- create_abatement_cost_spaghetti_plot(solution_list, title = "e) Abatement costs") +
    coord_cartesian(ylim = shared_cost_limits)
  
  p6 <- create_damage_cost_spaghetti_plot(solution_list, title = "f) Damage costs") +
    coord_cartesian(ylim = shared_cost_limits)
  
  # Combine plots in 3×2 grid layout
  combined <- (p1 + p2) /
    (p3 + p4) /
    (p5 + p6)
  
  # Count successful runs for subtitle
  n_runs <- length(solution_list)
  
  # Add overall title and subtitle
  final_plot <- combined + 
    plot_annotation(
      title = "Parameter Importance: Optimal Control Results",
      subtitle = paste0("Based on ", n_runs, " parameter combinations"),
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        plot.margin = margin(10, 0, 10, 0)
      )
    )
  
  # Save plot if requested
  if (save_plot) {
    if (is.null(filename)) {
      # Extract scenario name from results for filename
      scenario_name <- "unknown_scenario"
      if ("run_info" %in% names(importance_results) && 
          "scenario" %in% names(importance_results$run_info)) {
        scenario_name <- importance_results$run_info$scenario
        # Clean scenario name for filename (remove spaces, special characters)
        scenario_name <- gsub("[^A-Za-z0-9-]", "_", scenario_name)
      }
      
      filename <- paste0("parameter_importance_", scenario_name, "_", 
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    }
    
    filepath <- here::here("figs", filename)
    ggsave(filepath, final_plot, width = 190, height = 260, units = "mm", device = cairo_pdf)
    
    if (verbose) {
      cat("Parameter importance dashboard saved to:", filepath, "\n")
    }
  }
  
  # Print insights if requested
  if (print_insights) {
    print_parameter_importance_insights(importance_results, verbose = verbose)
  }
  
  return(final_plot)
}

# ============================================================================
# Analysis and insights function
# ============================================================================

#' @title Print Parameter Importance Insights
#' @description
#' Prints detailed analytical insights from parameter importance analysis results,
#' including temperature ranges, cost patterns, strategy patterns, convergence 
#' statistics, and parameter sensitivity observations
#' 
#' @param importance_results Full results object from run_parameter_importance()
#' @param verbose Whether to print section headers (default: TRUE)
#' @return Invisibly returns the summary data frame with calculated metrics
print_parameter_importance_insights <- function(importance_results, verbose = TRUE) {
  
  # Extract summary
  if (is.null(importance_results$summary)) {
    if (verbose) {
      cat("No summary available for insights analysis\n")
    }
    return(invisible(NULL))
  }
  
  summary_df <- importance_results$summary
  
  if (verbose) {
    cat("\n=== PARAMETER IMPORTANCE INSIGHTS ===\n")
  }
  
  cat("Parameter sets analysed:", nrow(summary_df), "\n")
  
  # Temperature insights
  cat("\nTemperature patterns:\n")
  cat("- Peak temperature range:", 
      sprintf("%.2f°C to %.2f°C", 
              min(summary_df$peak_temperature, na.rm = TRUE),
              max(summary_df$peak_temperature, na.rm = TRUE)), "\n")
  cat("- Mean peak temperature:", sprintf("%.2f°C", mean(summary_df$peak_temperature, na.rm = TRUE)), "\n")
  cat("- Standard deviation:", sprintf("%.2f°C", sd(summary_df$peak_temperature, na.rm = TRUE)), "\n")
  
  # Final temperature insights
  cat("\nFinal temperature (2100):\n")
  cat("- Final temperature range:", 
      sprintf("%.2f°C to %.2f°C", 
              min(summary_df$final_temperature, na.rm = TRUE),
              max(summary_df$final_temperature, na.rm = TRUE)), "\n")
  runs_above_1p5 <- sum(summary_df$final_temperature > 1.5, na.rm = TRUE)
  cat("- Runs with final temperature > 1.5°C:", runs_above_1p5, 
      sprintf("(%.1f%%)", 100 * runs_above_1p5 / nrow(summary_df)), "\n")
  
  # Cost insights
  cat("\nCost patterns:\n")
  cat("- Total cost range:", 
      sprintf("%.1f to %.1f trillion $", 
              min(summary_df$total_cost, na.rm = TRUE),
              max(summary_df$total_cost, na.rm = TRUE)), "\n")
  cat("- Mean total cost:", sprintf("%.1f trillion $", mean(summary_df$total_cost, na.rm = TRUE)), "\n")
  
  # Cost component insights
  cat("- Mean mitigation cost:", sprintf("%.1f trillion $", mean(summary_df$mitig_cost, na.rm = TRUE)), "\n")
  cat("- Mean CDR cost:", sprintf("%.1f trillion $", mean(summary_df$remov_cost, na.rm = TRUE)), "\n")
  cat("- Mean damage cost:", sprintf("%.1f trillion $", mean(summary_df$temp_cost, na.rm = TRUE)), "\n")
  
  # Strategy insights
  cat("\nStrategy patterns:\n")
  cat("- Mitigation range:", 
      sprintf("%.1f to %.1f GtCO2", 
              min(summary_df$total_mitigation_units, na.rm = TRUE),
              max(summary_df$total_mitigation_units, na.rm = TRUE)), "\n")
  cat("- CDR range:", 
      sprintf("%.1f to %.1f GtCO2", 
              min(summary_df$total_cdr_units, na.rm = TRUE),
              max(summary_df$total_cdr_units, na.rm = TRUE)), "\n")
  
  # Calculate mitigation vs CDR ratios
  summary_df$mitig_cdr_ratio <- summary_df$total_mitigation_units / 
    pmax(summary_df$total_cdr_units, 0.1)  # Avoid division by zero
  
  cat("- Mean mitigation/CDR ratio:", sprintf("%.2f", mean(summary_df$mitig_cdr_ratio, na.rm = TRUE)), "\n")
  cat("- Ratio range:", 
      sprintf("%.2f to %.2f", 
              min(summary_df$mitig_cdr_ratio, na.rm = TRUE),
              max(summary_df$mitig_cdr_ratio, na.rm = TRUE)), "\n")
  
  # Temperature overshoot insights
  if ("years_above_1p5" %in% names(summary_df)) {
    cat("\nTemperature overshoot:\n")
    cat("- Mean years above 1.5°C:", 
        sprintf("%.1f years", mean(summary_df$years_above_1p5, na.rm = TRUE)), "\n")
    cat("- Range:", 
        sprintf("%d to %d years", 
                min(summary_df$years_above_1p5, na.rm = TRUE),
                max(summary_df$years_above_1p5, na.rm = TRUE)), "\n")
  }
  
  # Convergence insights
  if ("converged" %in% names(summary_df)) {
    n_converged <- sum(summary_df$converged, na.rm = TRUE)
    cat("\nConvergence statistics:\n")
    cat("- Converged runs:", n_converged, "out of", nrow(summary_df),
        sprintf("(%.1f%%)", 100 * n_converged / nrow(summary_df)), "\n")
    
    if ("iterations" %in% names(summary_df)) {
      cat("- Mean iterations:", sprintf("%.1f", mean(summary_df$iterations, na.rm = TRUE)), "\n")
      cat("- Iteration range:", 
          sprintf("%d to %d", 
                  min(summary_df$iterations, na.rm = TRUE),
                  max(summary_df$iterations, na.rm = TRUE)), "\n")
    }
  }
  
  # Cost efficiency insights
  summary_df$cost_per_degree <- summary_df$total_cost / 
    pmax(summary_df$peak_temperature - 1.0, 0.1)  # Cost per degree above 1°C baseline
  
  cat("\nCost efficiency:\n")
  cat("- Mean cost per degree of peak warming:", 
      sprintf("%.1f trillion $/°C", mean(summary_df$cost_per_degree, na.rm = TRUE)), "\n")
  
  # Parameter sensitivity observations
  cat("\nParameter sensitivity observations:\n")
  cat("- Wide variation in outcomes indicates strong parameter sensitivity\n")
  cat("- Temperature range of", 
      sprintf("%.2f°C", max(summary_df$peak_temperature, na.rm = TRUE) - min(summary_df$peak_temperature, na.rm = TRUE)),
      "suggests key parameters significantly affect warming trajectory\n")
  cat("- Cost range of", 
      sprintf("%.1f trillion $", max(summary_df$total_cost, na.rm = TRUE) - min(summary_df$total_cost, na.rm = TRUE)),
      "suggests key parameters significantly affect economic outcomes\n")
  cat("- Strategy mix variation indicates parameter influence on optimal mitigation/CDR balance\n")
  
  # Return augmented summary invisibly
  invisible(summary_df)
}

# ============================================================================
# Usage examples
# ============================================================================

# Assuming you have importance_results from run_parameter_importance():
#
# # Create individual plots
# emissions_plot <- create_emissions_spaghetti_plot(importance_results$successful_runs)
# temp_plot <- create_temperature_spaghetti_plot(importance_results$successful_runs)
# mitigation_plot <- create_mitigation_spaghetti_plot(importance_results$successful_runs)
# cdr_plot <- create_cdr_spaghetti_plot(importance_results$successful_runs)
# abatement_cost_plot <- create_abatement_cost_spaghetti_plot(importance_results$successful_runs)
# damage_cost_plot <- create_damage_cost_spaghetti_plot(importance_results$successful_runs)
#
# # View individual plot
# print(temp_plot)
#
# # Create complete dashboard
# dashboard <- create_parameter_importance_dashboard(importance_results)
# print(dashboard)
#
# # Create dashboard with insights printed
# dashboard <- create_parameter_importance_dashboard(
#   importance_results = importance_results,
#   print_insights = TRUE,
#   verbose = TRUE
# )
#
# # Save dashboard to file
# create_parameter_importance_dashboard(
#   importance_results = importance_results,
#   save_plot = TRUE,
#   filename = "my_parameter_importance.pdf"
# )
#
# # Print insights separately
# print_parameter_importance_insights(importance_results)
#
# # Filter results and create dashboard from filtered data
# filtered_results <- filter_importance_results(
#   importance_results = importance_results,
#   max_final_temperature = 1.6
# )
# dashboard <- create_parameter_importance_dashboard(filtered_results)