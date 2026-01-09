# ==============================================================================
# Scenario Comparison Visualisation Functions
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

#' @title Scenario Comparison Visualisation Functions
#' @description
#' Functions to create plots comparing optimal control results across multiple 
#' SSP scenarios. All scenarios are plotted on the same axes with different colors.
#'
#' @note All required libraries (ggplot2, dplyr, patchwork, purrr, scales, cowplot, here) 
#' must be loaded before using these functions.
#'
#' @author Nina Rynne
#' @date January 2026

# ============================================================================
# Global definitions
# ============================================================================

# Define colorblind-friendly palette for scenarios
# Colors selected for accessibility and distinction between 5 SSP scenarios
scenario_colors <- c(
  "#00ADCF",  # Cyan
  "#173C66",  # Dark blue
  "#F0E442",  # Yellow
  "#E71D25",  # Red
  "#951B1E"   # Dark red
)

# Define consistent theme for all scenario comparison plots
scenario_comparison_theme <- theme_bw() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "none"
  )

# ============================================================================
# Helper functions
# ============================================================================

#' @title Clean Scenario Names for Display
#' @description
#' Removes "-Baseline" suffix from scenario names for cleaner plot legends
#' 
#' @param scenario_names Character vector of scenario names
#' @return Character vector with cleaned names
#' 
#' @examples
#' clean_scenario_names(c("SSP1-Baseline", "SSP2-Baseline"))
#' # Returns: c("SSP1", "SSP2")
clean_scenario_names <- function(scenario_names) {
  gsub("-Baseline$", "", scenario_names)
}

# ============================================================================
# Individual plot functions
# ============================================================================

#' @title Plot Cumulative Emissions Across Scenarios
#' @description
#' Creates a line plot showing cumulative CO2 emissions over time
#' for multiple SSP scenarios
#' 
#' @param scenario_results Named list of solution objects from run_scenario_comparison
#' @param title Plot title (optional)
#' @return ggplot object
plot_scenario_comparison_emissions <- function(scenario_results, 
                                               title = "Cumulative Emissions") {
  
  # Extract cumulative emissions data from all scenarios
  plot_data <- map_dfr(names(scenario_results), function(scenario_name) {
    result <- scenario_results[[scenario_name]]
    if (!is.null(result)) {
      data.frame(
        scenario = scenario_name,
        years = result$years,
        cumulative_emissions = result$cumulative_emissions
      )
    }
  })
  
  # Validate that data was found
  if (nrow(plot_data) == 0) {
    stop("No valid emissions data found across scenarios")
  }
  
  # Clean scenario names for display
  plot_data$scenario_clean <- clean_scenario_names(plot_data$scenario)
  
  # Create cumulative emissions plot
  p <- ggplot(plot_data, aes(x = years, y = cumulative_emissions, color = scenario_clean)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      name = "Scenario",
      values = scenario_colors[1:length(unique(plot_data$scenario_clean))]
    ) +
    labs(
      title = title,
      x = "Year",
      y = expression("Cumulative emissions (GtCO"[2]*")")
    ) +
    scenario_comparison_theme
  
  return(p)
}

#' @title Plot Temperature Trajectories Across Scenarios
#' @description
#' Creates a line plot showing temperature anomaly trajectories over time
#' for multiple SSP scenarios with a 1.5°C reference line
#' 
#' @param scenario_results Named list of solution objects from run_scenario_comparison
#' @param title Plot title (optional)
#' @return ggplot object
plot_scenario_comparison_temperature <- function(scenario_results, 
                                                 title = "Temperature Trajectories") {
  
  # Extract temperature data from all scenarios
  plot_data <- map_dfr(names(scenario_results), function(scenario_name) {
    result <- scenario_results[[scenario_name]]
    if (!is.null(result)) {
      data.frame(
        scenario = scenario_name,
        years = result$years,
        temperature = result$temperature_anomaly
      )
    }
  })
  
  # Validate that data was found
  if (nrow(plot_data) == 0) {
    stop("No valid temperature data found across scenarios")
  }
  
  # Clean scenario names for display
  plot_data$scenario_clean <- clean_scenario_names(plot_data$scenario)
  
  # Create temperature trajectory plot
  p <- ggplot(plot_data, aes(x = years, y = temperature, color = scenario_clean)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1.5, linetype = "dashed", alpha = 0.7, color = "red") +
    scale_color_manual(
      name = "Scenario",
      values = scenario_colors[1:length(unique(plot_data$scenario_clean))]
    ) +
    scale_y_continuous(limits = c(NA, 2)) +
    labs(
      title = title,
      x = "Year",
      y = "Temperature anomaly (°C)"
    ) +
    scenario_comparison_theme
  
  return(p)
}

#' @title Plot Mitigation Strategies Across Scenarios
#' @description
#' Creates a line plot showing annual mitigation quantities over time
#' for multiple SSP scenarios
#' 
#' @param scenario_results Named list of solution objects from run_scenario_comparison
#' @param title Plot title (optional)
#' @return ggplot object
plot_scenario_comparison_mitigation <- function(scenario_results, 
                                                title = "Mitigation Strategies") {
  
  # Extract mitigation data from all scenarios
  plot_data <- map_dfr(names(scenario_results), function(scenario_name) {
    result <- scenario_results[[scenario_name]]
    if (!is.null(result)) {
      data.frame(
        scenario = scenario_name,
        years = result$years,
        mitigation = result$qty_mitig
      )
    }
  })
  
  # Validate that data was found
  if (nrow(plot_data) == 0) {
    stop("No valid mitigation data found across scenarios")
  }
  
  # Clean scenario names for display
  plot_data$scenario_clean <- clean_scenario_names(plot_data$scenario)
  
  # Create mitigation strategies plot
  p <- ggplot(plot_data, aes(x = years, y = mitigation, color = scenario_clean)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      name = "Scenario",
      values = scenario_colors[1:length(unique(plot_data$scenario_clean))]
    ) +
    labs(
      title = title,
      x = "Year",
      y = expression("Mitigation (GtCO"[2]*"/yr)")
    ) +
    scenario_comparison_theme
  
  return(p)
}

#' @title Plot CDR Strategies Across Scenarios
#' @description
#' Creates a line plot showing annual carbon dioxide removal (CDR) quantities 
#' over time for multiple SSP scenarios
#' 
#' @param scenario_results Named list of solution objects from run_scenario_comparison
#' @param title Plot title (optional)
#' @return ggplot object
plot_scenario_comparison_cdr <- function(scenario_results, 
                                         title = "CDR Strategies") {
  
  # Extract CDR data from all scenarios
  plot_data <- map_dfr(names(scenario_results), function(scenario_name) {
    result <- scenario_results[[scenario_name]]
    if (!is.null(result)) {
      data.frame(
        scenario = scenario_name,
        years = result$years,
        cdr = result$qty_remov
      )
    }
  })
  
  # Validate that data was found
  if (nrow(plot_data) == 0) {
    stop("No valid CDR data found across scenarios")
  }
  
  # Clean scenario names for display
  plot_data$scenario_clean <- clean_scenario_names(plot_data$scenario)
  
  # Create CDR strategies plot
  p <- ggplot(plot_data, aes(x = years, y = cdr, color = scenario_clean)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      name = "Scenario",
      values = scenario_colors[1:length(unique(plot_data$scenario_clean))]
    ) +
    labs(
      title = title,
      x = "Year",
      y = expression("CDR (GtCO"[2]*"/yr)")
    ) +
    scenario_comparison_theme
  
  return(p)
}

#' @title Plot Abatement Costs Across Scenarios
#' @description
#' Creates a line plot showing annual abatement costs (mitigation + CDR) as a 
#' percentage of baseline GWP over time for multiple SSP scenarios
#' 
#' @param scenario_results Named list of solution objects from run_scenario_comparison
#' @param title Plot title (optional)
#' @return ggplot object
plot_scenario_comparison_abatement_cost <- function(scenario_results, 
                                                    title = "Abatement Costs") {
  
  # Extract abatement cost data from all scenarios
  plot_data <- map_dfr(names(scenario_results), function(scenario_name) {
    result <- scenario_results[[scenario_name]]
    if (!is.null(result)) {
      # Calculate abatement cost as sum of mitigation and removal costs
      abatement_cost_annual <- result$mitig_costs_annual + result$remov_costs_annual
      
      # Normalize by baseline GWP to get proportion
      abatement_cost_proportion <- abatement_cost_annual / result$baseline_annual_gwp
      
      data.frame(
        scenario = scenario_name,
        years = result$years,
        abatement_cost_proportion = abatement_cost_proportion
      )
    }
  })
  
  # Validate that data was found
  if (nrow(plot_data) == 0) {
    stop("No valid abatement cost data found across scenarios")
  }
  
  # Clean scenario names for display
  plot_data$scenario_clean <- clean_scenario_names(plot_data$scenario)
  
  # Create abatement costs plot
  p <- ggplot(plot_data, aes(x = years, y = abatement_cost_proportion, color = scenario_clean)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      name = "Scenario",
      values = scenario_colors[1:length(unique(plot_data$scenario_clean))]
    ) +
    scale_y_continuous(labels = percent) +
    labs(
      title = title,
      x = "Year",
      y = "Abatement Costs\n(% of GWP)"
    ) +
    scenario_comparison_theme
  
  return(p)
}

#' @title Plot Temperature-Related Damage Costs Across Scenarios
#' @description
#' Creates a line plot showing annual temperature-related damage costs as a 
#' percentage of baseline GWP over time for multiple SSP scenarios
#' 
#' @param scenario_results Named list of solution objects from run_scenario_comparison
#' @param title Plot title (optional)
#' @return ggplot object
plot_scenario_comparison_damage_cost <- function(scenario_results, 
                                                 title = "Damage Costs") {
  
  # Extract damage cost data from all scenarios
  plot_data <- map_dfr(names(scenario_results), function(scenario_name) {
    result <- scenario_results[[scenario_name]]
    if (!is.null(result)) {
      # Use temperature costs as damage costs
      damage_cost_annual <- result$temp_costs_annual
      
      # Normalize by baseline GWP to get proportion
      damage_cost_proportion <- damage_cost_annual / result$baseline_annual_gwp
      
      data.frame(
        scenario = scenario_name,
        years = result$years,
        damage_cost_proportion = damage_cost_proportion
      )
    }
  })
  
  # Validate that data was found
  if (nrow(plot_data) == 0) {
    stop("No valid damage cost data found across scenarios")
  }
  
  # Clean scenario names for display
  plot_data$scenario_clean <- clean_scenario_names(plot_data$scenario)
  
  # Create damage costs plot
  p <- ggplot(plot_data, aes(x = years, y = damage_cost_proportion, color = scenario_clean)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      name = "Scenario",
      values = scenario_colors[1:length(unique(plot_data$scenario_clean))]
    ) +
    scale_y_continuous(labels = percent) +
    labs(
      title = title,
      x = "Year",
      y = "Damage Costs\n(% of GWP)"
    ) +
    scenario_comparison_theme
  
  return(p)
}

# ============================================================================
# Dashboard function
# ============================================================================

#' @title Create Scenario Comparison Dashboard
#' @description
#' Creates a comprehensive dashboard with 6 plots comparing optimal control 
#' results across multiple SSP scenarios. Plots are arranged in a 3×2 grid 
#' with shared y-axis limits for comparable plots and a shared legend.
#'
#' @param scenario_results Output from run_scenario_comparison(), or just the 
#'   scenario_results component
#' @param save_plot Whether to save the plot as PDF (default: FALSE)
#' @param filename Custom filename for saving (default: auto-generated with timestamp)
#' @param print_insights Whether to print analytical insights after creating dashboard (default: FALSE)
#' @param verbose Print progress information (default: TRUE)
#'
#' @return Combined ggplot object containing the dashboard
#'
#' @examples
#' # Create dashboard from scenario comparison results
#' dashboard <- create_scenario_comparison_dashboard(
#'   scenario_results = results,
#'   save_plot = TRUE,
#'   print_insights = TRUE)
#' 
create_scenario_comparison_dashboard <- function(scenario_results, 
                                                 save_plot = FALSE,
                                                 filename = NULL,
                                                 print_insights = FALSE,
                                                 verbose = TRUE,
                                                 include_temperature = TRUE,
                                                 include_emissions = TRUE,
                                                 include_controls = TRUE,
                                                 include_costs = TRUE,
                                                 color_palette = "viridis",
                                                 width = 297,      # A4 landscape width (mm)
                                                 height = 210) {     # A4 landscape height (mm)
  
  # Extract scenario results component from full results object
  if ("scenario_results" %in% names(scenario_results)) {
    scenario_data <- scenario_results$scenario_results
  } else {
    # Assume it's already the scenario_results component
    scenario_data <- scenario_results
  }
  
  # Validate that we have scenario data
  if (length(scenario_data) == 0) {
    stop("No scenario results found. Ensure run_scenario_comparison() completed successfully.")
  }
  
  if (verbose) {
    cat("Creating scenario comparison dashboard for", length(scenario_data), "scenarios\n")
  }
  
  # Calculate shared y-axis limits for mitigation and CDR plots
  mitigation_data <- map_dfr(names(scenario_data), function(scenario_name) {
    result <- scenario_data[[scenario_name]]
    if (!is.null(result)) {
      data.frame(mitigation = result$qty_mitig)
    }
  })
  
  cdr_data <- map_dfr(names(scenario_data), function(scenario_name) {
    result <- scenario_data[[scenario_name]]
    if (!is.null(result)) {
      data.frame(cdr = result$qty_remov)
    }
  })
  
  shared_abatement_limits <- c(0, max(c(mitigation_data$mitigation, cdr_data$cdr), na.rm = TRUE))
  
  # Calculate shared y-axis limits for cost plots
  abatement_cost_data <- map_dfr(names(scenario_data), function(scenario_name) {
    result <- scenario_data[[scenario_name]]
    if (!is.null(result)) {
      abatement_cost_annual <- result$mitig_costs_annual + result$remov_costs_annual
      abatement_cost_proportion <- abatement_cost_annual / result$baseline_annual_gwp
      data.frame(abatement_cost = abatement_cost_proportion)
    }
  })
  
  damage_cost_data <- map_dfr(names(scenario_data), function(scenario_name) {
    result <- scenario_data[[scenario_name]]
    if (!is.null(result)) {
      damage_cost_annual <- result$temp_costs_annual
      damage_cost_proportion <- damage_cost_annual / result$baseline_annual_gwp
      data.frame(damage_cost = damage_cost_proportion)
    }
  })
  
  shared_cost_limits <- c(0, max(c(abatement_cost_data$abatement_cost, damage_cost_data$damage_cost), na.rm = TRUE))
  
  # Create individual plots with panel labels
  p1 <- plot_scenario_comparison_emissions(scenario_data) + 
    labs(title = "a) Cumulative emissions")
  
  p2 <- plot_scenario_comparison_temperature(scenario_data) + 
    labs(title = "b) Temperature trajectories")
  
  p3 <- plot_scenario_comparison_mitigation(scenario_data) + 
    labs(title = "c) Mitigation strategies") +
    coord_cartesian(ylim = shared_abatement_limits)
  
  p4 <- plot_scenario_comparison_cdr(scenario_data) + 
    labs(title = "d) CDR strategies") +
    coord_cartesian(ylim = shared_abatement_limits)
  
  p5 <- plot_scenario_comparison_abatement_cost(scenario_data) + 
    labs(title = "e) Abatement costs") +
    coord_cartesian(ylim = shared_cost_limits)
  
  p6 <- plot_scenario_comparison_damage_cost(scenario_data) + 
    labs(title = "f) Temperature-related damages") +
    coord_cartesian(ylim = shared_cost_limits)
  
  # Extract shared legend from one of the plots
  shared_legend <- get_legend(p1 + theme(legend.position = "right"))
  
  # Combine plots in 3×2 grid layout
  combined <- (p1 + p2) /
    (p3 + p4) /
    (p5 + p6)
  
  # Add shared legend to the right of the grid
  final_plots <- combined | shared_legend
  
  # Set relative widths: wide for main plots, narrow for legend
  final_plots <- final_plots + plot_layout(widths = c(10, 1))
  
  # Count successful scenarios for subtitle
  n_scenarios <- length(scenario_data)
  
  # Add overall title and subtitle
  final_plot <- final_plots + 
    plot_annotation(
      title = "Scenario Comparison: Optimal Control Results",
      subtitle = paste0("Comparing ", n_scenarios, " scenarios"),
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        plot.margin = margin(10, 0, 10, 0)
      )
    )
  
  # Save plot if requested
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0("scenario_comparison_dashboard_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    }
    
    filepath <- here::here("figs", filename)
    ggsave(filepath, final_plot, width = 190, height = 260, units = "mm", device = cairo_pdf)
    
    if (verbose) {
      cat("Scenario comparison dashboard saved to:", filepath, "\n")
    }
  }
  
  # Print insights if requested
  if (print_insights) {
    print_scenario_insights(scenario_results, verbose = verbose)
  }
  
  return(final_plot)
}

# ============================================================================
# Analysis and insights function
# ============================================================================

#' @title Print Scenario Comparison Insights
#' @description
#' Prints detailed analytical insights from scenario comparison results,
#' including temperature patterns, cost patterns, strategy patterns, and 
#' policy recommendations
#' 
#' @param scenario_results Full results object from run_scenario_comparison()
#' @param verbose Whether to print section headers (default: TRUE)
#' @return Invisibly returns the summary data frame with calculated metrics
print_scenario_insights <- function(scenario_results, verbose = TRUE) {
  
  # Extract comparison summary
  if (is.null(scenario_results$comparison_summary)) {
    if (verbose) {
      cat("No comparison summary available for insights analysis\n")
    }
    return(invisible(NULL))
  }
  
  summary_df <- scenario_results$comparison_summary
  
  if (verbose) {
    cat("\n=== SCENARIO COMPARISON INSIGHTS ===\n")
  }
  
  cat("Scenarios analysed:", nrow(summary_df), "\n")
  
  # Temperature insights
  cat("\nTemperature patterns:\n")
  cat("- Peak temperature range:", 
      sprintf("%.2f°C (%s) to %.2f°C (%s)", 
              min(summary_df$peak_temperature, na.rm = TRUE),
              summary_df$scenario[which.min(summary_df$peak_temperature)],
              max(summary_df$peak_temperature, na.rm = TRUE),
              summary_df$scenario[which.max(summary_df$peak_temperature)]), "\n")
  
  # Cost insights
  cat("\nCost patterns:\n")
  cat("- Total cost range:", 
      sprintf("%.1f (%s) to %.1f (%s) trillion $", 
              min(summary_df$total_cost, na.rm = TRUE),
              summary_df$scenario[which.min(summary_df$total_cost)],
              max(summary_df$total_cost, na.rm = TRUE),
              summary_df$scenario[which.max(summary_df$total_cost)]), "\n")
  
  # Strategy insights
  cat("\nStrategy patterns:\n")
  
  # Mitigation patterns
  min_mitig_scenario <- summary_df$scenario[which.min(summary_df$total_mitigation_units)]
  max_mitig_scenario <- summary_df$scenario[which.max(summary_df$total_mitigation_units)]
  cat("- Mitigation reliance: Lowest in", min_mitig_scenario, 
      ", Highest in", max_mitig_scenario, "\n")
  
  # CDR patterns
  min_cdr_scenario <- summary_df$scenario[which.min(summary_df$total_cdr_units)]
  max_cdr_scenario <- summary_df$scenario[which.max(summary_df$total_cdr_units)]
  cat("- CDR reliance: Lowest in", min_cdr_scenario, 
      ", Highest in", max_cdr_scenario, "\n")
  
  # Calculate mitigation vs CDR ratios
  summary_df$mitig_cdr_ratio <- summary_df$total_mitigation_units / 
    pmax(summary_df$total_cdr_units, 0.1)  # Avoid division by zero
  
  balanced_scenario <- summary_df$scenario[which.min(abs(summary_df$mitig_cdr_ratio - 1))]
  mitig_heavy_scenario <- summary_df$scenario[which.max(summary_df$mitig_cdr_ratio)]
  cdr_heavy_scenario <- summary_df$scenario[which.min(summary_df$mitig_cdr_ratio)]
  
  cat("- Most balanced strategy (mitigation ≈ CDR):", balanced_scenario, "\n")
  cat("- Most mitigation-heavy strategy:", mitig_heavy_scenario, "\n")
  cat("- Most CDR-heavy strategy:", cdr_heavy_scenario, "\n")
  
  # Temperature overshoot insights
  if ("years_above_1p5" %in% names(summary_df)) {
    cat("\nTemperature overshoot:\n")
    min_overshoot <- summary_df$scenario[which.min(summary_df$years_above_1p5)]
    max_overshoot <- summary_df$scenario[which.max(summary_df$years_above_1p5)]
    cat("- Shortest overshoot period:", min_overshoot, 
        "(", min(summary_df$years_above_1p5), "years )\n")
    cat("- Longest overshoot period:", max_overshoot,
        "(", max(summary_df$years_above_1p5), "years )\n")
  }
  
  # Cost efficiency insights
  summary_df$cost_per_degree_avoided <- summary_df$total_cost / 
    pmax(summary_df$peak_temperature - 1.0, 0.1)  # Cost per degree above 1°C baseline
  
  most_efficient <- summary_df$scenario[which.min(summary_df$cost_per_degree_avoided)]
  least_efficient <- summary_df$scenario[which.max(summary_df$cost_per_degree_avoided)]
  
  cat("\nCost efficiency:\n")
  cat("- Most cost-efficient scenario:", most_efficient, "\n")
  cat("- Least cost-efficient scenario:", least_efficient, "\n")
  
  # Policy recommendations
  cat("\nKey policy insights:\n")
  cat("- Baseline emissions pathway significantly affects optimal strategy\n")
  cat("- Higher baseline emissions generally require more CDR deployment\n")
  cat("- Trade-offs between cost, temperature, and strategy mix vary by scenario\n")
  cat("- Early action scenarios (lower baselines) tend to be more cost-effective\n")
  
  # Return augmented summary invisibly
  invisible(summary_df)
}

# ============================================================================
# Usage examples
# ============================================================================

# Assuming you have scenario_results from run_scenario_comparison():
#
# # Create individual plots
# emissions_plot <- plot_scenario_comparison_emissions(scenario_results$scenario_results)
# temp_plot <- plot_scenario_comparison_temperature(scenario_results$scenario_results)
# mitigation_plot <- plot_scenario_comparison_mitigation(scenario_results$scenario_results)
# cdr_plot <- plot_scenario_comparison_cdr(scenario_results$scenario_results)
# abatement_cost_plot <- plot_scenario_comparison_abatement_cost(scenario_results$scenario_results)
# damage_cost_plot <- plot_scenario_comparison_damage_cost(scenario_results$scenario_results)
#
# # View individual plot
# print(temp_plot)
#
# # Create complete dashboard
# dashboard <- create_scenario_comparison_dashboard(scenario_results$scenario_results)
# print(dashboard)
#
# # Save dashboard to file
# create_scenario_comparison_dashboard(
#   scenario_results = scenario_results$scenario_results,
#   save_plot = TRUE,
#   filename = "my_scenario_comparison.pdf"
# )