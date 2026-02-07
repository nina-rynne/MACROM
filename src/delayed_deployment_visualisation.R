# ==============================================================================
# Delayed Deployment Visualisation Functions
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

#' @title Delayed Deployment Heatmap Dashboard Functions
#' @description
#' Functions to create heatmap dashboards showing how delays in mitigation and 
#' CDR deployment affect climate outcomes and costs across multiple SSP scenarios.
#' Supports visualization of various metrics including peak temperature, years 
#' above 1.5°C, abatement costs, and temperature-related damage costs.
#' 
#' The main function creates flexible dashboards that can display:
#' - Single or multiple variables
#' - Optional contour lines
#' - Optional gradient vector fields (arrows)
#' - Infeasible parameter combinations
#' - Customizable color palettes
#' - Shared or independent scales across variables
#'
#' @note All required libraries (ggplot2, dplyr, patchwork, purrr, viridis, 
#' cowplot, here, metR) must be loaded before using these functions.
#'
#' @author Nina Rynne
#' @date January 2026

# ==============================================================================
# Global definitions
# ==============================================================================

# Define default color palettes for different variables
# These can be overridden when calling the main function
DELAYED_DEPLOYMENT_PALETTES <- list(
  peak_temperature = list(option = "plasma", direction = -1),
  years_above_1p5 = list(option = "plasma", direction = -1),
  abatement_cost = list(option = "viridis", direction = -1),
  temp_cost = list(option = "viridis", direction = -1),
  total_cost = list(option = "viridis", direction = -1),
  mitig_cost = list(option = "viridis", direction = -1),
  remov_cost = list(option = "viridis", direction = -1)
)

# Define default variable labels for legends
# Note: All abbreviations use full stops for scientific publication standards
DELAYED_DEPLOYMENT_LABELS <- list(
  peak_temperature = "Peak\nTemp.\n(°C)",
  years_above_1p5 = "Years\nAbove\n1.5°C",
  abatement_cost = "Costs\n($USD trillion)",
  temp_cost = "Temp.\nDamage\nCost\n($ trillion)",
  total_cost = "Total\nCost\n($ trillion)",
  mitig_cost = "Mitigation\nCost\n($ trillion)",
  remov_cost = "Removal\nCost\n($ trillion)",
  cost = "Cost\n($ trillion)"
)

# Single-line versions for horizontal/bottom legends
DELAYED_DEPLOYMENT_LABELS_SINGLE_LINE <- list(
  peak_temperature = "Peak Temperature (°C)",
  years_above_1p5 = "Years Above 1.5°C",
  abatement_cost = "Costs ($USD trillion)",
  temp_cost = "Temperature Damage Cost ($ trillion)",
  total_cost = "Total Cost ($ trillion)",
  mitig_cost = "Mitigation Cost ($ trillion)",
  remov_cost = "Removal Cost ($ trillion)",
  cost = "Cost ($ trillion)"
)

# Define SSP scenario order for consistent plotting
SSP_SCENARIO_ORDER <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")

# Define base theme for delayed deployment heatmaps
# This will be adjusted based on layout type (single-row vs multi-row)
# Define base theme for delayed deployment heatmaps
# This will be adjusted based on layout type (single-row vs multi-row)
BASE_DELAYED_DEPLOYMENT_THEME <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# ==============================================================================
# Section 1: Data Preparation Helpers
# ==============================================================================
#
# These functions handle the extraction, filtering, validation, and preprocessing
# of delayed deployment data before visualization. They ensure data consistency
# across all plots and calculate derived quantities like gradients when needed
# for arrow visualizations.
#
# Functions in this section:
# - prepare_delayed_deployment_data(): Extract and validate core data
# - calculate_delayed_deployment_gradients(): Compute directional gradients
# - calculate_variable_limits(): Determine shared or independent scales
# ==============================================================================

#' @title Prepare Delayed Deployment Data for Visualization
#' @description
#' Extracts and validates data from delayed deployment results across scenarios.
#' Filters by SSP scenarios, orders data consistently, validates required columns,
#' and removes any rows with missing values in key variables.
#' 
#' @param deployment_results Results object from delayed deployment analysis
#'   Can be either a data frame directly or a list containing 'combined_results'
#'   Must contain columns:
#'   - scenario: Full scenario name
#'   - scenario_short: Short scenario code (e.g., "SSP1", "SSP2")
#'   - mitigation_delay: Years of mitigation deployment delay
#'   - cdr_delay: Years of CDR deployment delay
#'   - feasible: Logical indicating if combination is feasible
#'   - Plus any variables to be plotted (e.g., peak_temperature, abatement_cost)
#' @param variables Character vector of variable names to validate and prepare
#'   Default is NULL, which performs basic validation only
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' 
#' @return Data frame with filtered, ordered, and validated delayed deployment data
#'   Returns data with scenario_short as an ordered factor and sorted by 
#'   scenario_short, mitigation_delay, and cdr_delay
#' 
#' @examples
#' # Basic preparation
#' plot_data <- prepare_delayed_deployment_data(multi_results)
#' 
#' # Prepare for specific variables
#' plot_data <- prepare_delayed_deployment_data(
#'   multi_results, 
#'   variables = c("peak_temperature", "years_above_1p5")
#' )
prepare_delayed_deployment_data <- function(deployment_results,
                                            variables = NULL,
                                            verbose = TRUE) {
  
  # Validate input structure
  if (is.null(deployment_results)) {
    stop("deployment_results cannot be NULL")
  }
  
  # Extract combined results - handle both data frame and list structures
  if (is.data.frame(deployment_results)) {
    # Data is already a data frame
    combined_data <- deployment_results
  } else if (is.list(deployment_results) && "combined_results" %in% names(deployment_results)) {
    # Data is in a list with combined_results element
    combined_data <- deployment_results$combined_results
  } else {
    stop("deployment_results must be either a data frame or a list containing 'combined_results'")
  }
  
  if (is.null(combined_data) || nrow(combined_data) == 0) {
    stop("No valid data found in deployment_results")
  }
  
  # Define required columns for delayed deployment analysis
  required_cols <- c("scenario", "scenario_short", "mitigation_delay", 
                     "cdr_delay", "feasible")
  
  # Check for required columns
  missing_cols <- setdiff(required_cols, names(combined_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # If variables specified, check they exist
  if (!is.null(variables)) {
    missing_vars <- setdiff(variables, names(combined_data))
    if (length(missing_vars) > 0) {
      stop("Missing requested variables: ", paste(missing_vars, collapse = ", "))
    }
  }
  
  # Filter to valid SSP scenarios only
  plot_data <- combined_data %>%
    filter(scenario_short %in% SSP_SCENARIO_ORDER)
  
  if (nrow(plot_data) == 0) {
    stop("No data found for SSP scenarios: ", paste(SSP_SCENARIO_ORDER, collapse = ", "))
  }
  
  # Keep scenario_short as character (don't convert to factor)
  # Factor conversion was causing issues with filtering in downstream functions
  
  # Remove rows with missing values in specified variables
  if (!is.null(variables)) {
    n_before <- nrow(plot_data)
    
    for (var in variables) {
      plot_data <- plot_data %>%
        filter(!is.na(.data[[var]]))
    }
    
    n_after <- nrow(plot_data)
    
    if (verbose && n_after < n_before) {
      cat(sprintf("Removed %d rows with missing values in specified variables\n", 
                  n_before - n_after))
    }
    
    if (n_after == 0) {
      stop("No valid data remaining after removing missing values")
    }
  }
  
  # Sort data for consistent plotting
  plot_data <- plot_data %>%
    arrange(scenario_short, mitigation_delay, cdr_delay)
  
  # Print summary if verbose
  if (verbose) {
    n_scenarios <- length(unique(plot_data$scenario_short))
    n_combinations <- nrow(plot_data)
    n_feasible <- sum(plot_data$feasible, na.rm = TRUE)
    
    cat(sprintf("Prepared delayed deployment data: %d scenarios, %d combinations, %d feasible\n",
                n_scenarios, n_combinations, n_feasible))
  }
  
  return(plot_data)
}

#' @title Calculate Gradients for Arrow Visualization
#' @description
#' Computes directional gradients (dx, dy) and gradient magnitude for a 
#' specified variable using central differences. These gradients are used
#' to create arrow vector fields showing the direction of steepest change
#' in the response surface.
#' 
#' Central differences provide more accurate gradient estimates than forward
#' or backward differences, especially for interior points. Edge points use
#' forward/backward differences as necessary.
#' 
#' @param data Data frame containing delayed deployment results
#'   Must have columns: scenario_short, mitigation_delay, cdr_delay, and the
#'   variable specified in 'variable' parameter
#' @param variable Character string naming the variable for gradient calculation
#'   (e.g., "abatement_cost", "temp_cost")
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' 
#' @return Data frame with three additional columns:
#'   - dx: Gradient in x-direction (mitigation_delay)
#'   - dy: Gradient in y-direction (cdr_delay)
#'   - gradient_mag: Magnitude of gradient vector (sqrt(dx^2 + dy^2))
#'   
#'   Rows with missing gradients (e.g., isolated points) will have NA values
#' 
#' @details
#' Gradients are calculated separately for each scenario using central differences:
#' - dx = (f(x+h) - f(x-h)) / (2h) where h is the step size
#' - dy = (f(y+k) - f(y-k)) / (2k) where k is the step size
#' 
#' @examples
#' # Calculate gradients for abatement cost
#' data_with_gradients <- calculate_delayed_deployment_gradients(
#'   plot_data, 
#'   variable = "abatement_cost"
#' )
#' 
#' # Calculate gradients for temperature damage cost
#' data_with_gradients <- calculate_delayed_deployment_gradients(
#'   plot_data, 
#'   variable = "temp_cost",
#'   verbose = FALSE
#' )
calculate_delayed_deployment_gradients <- function(data,
                                                   variable,
                                                   verbose = TRUE) {
  
  # Validate inputs
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("data must be a non-empty data frame")
  }
  
  if (!variable %in% names(data)) {
    stop("Variable '", variable, "' not found in data")
  }
  
  required_cols <- c("scenario_short", "mitigation_delay", "cdr_delay")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns for gradient calculation: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  if (verbose) {
    cat(sprintf("Calculating gradients for variable: %s\n", variable))
  }
  
  # Calculate gradients using central differences
  # Group by scenario and calculate gradients within each scenario
  
  # First, calculate dx (gradient in mitigation_delay direction)
  # Group by scenario and cdr_delay, then sort by mitigation_delay
  data_with_gradients <- data %>%
    group_by(scenario_short, cdr_delay) %>%
    arrange(mitigation_delay) %>%
    mutate(
      dx = (lead(.data[[variable]]) - lag(.data[[variable]])) / 
        (lead(mitigation_delay) - lag(mitigation_delay))
    ) %>%
    ungroup()
  
  # Second, calculate dy (gradient in cdr_delay direction)
  # Group by scenario and mitigation_delay, then sort by cdr_delay
  data_with_gradients <- data_with_gradients %>%
    group_by(scenario_short, mitigation_delay) %>%
    arrange(cdr_delay) %>%
    mutate(
      dy = (lead(.data[[variable]]) - lag(.data[[variable]])) / 
        (lead(cdr_delay) - lag(cdr_delay))
    ) %>%
    ungroup()
  
  # Calculate gradient magnitude
  data_with_gradients <- data_with_gradients %>%
    mutate(
      gradient_mag = sqrt(dx^2 + dy^2)
    )
  
  # Report summary statistics
  if (verbose) {
    n_total <- nrow(data_with_gradients)
    n_valid_gradients <- sum(!is.na(data_with_gradients$dx) & 
                               !is.na(data_with_gradients$dy))
    n_missing <- n_total - n_valid_gradients
    
    cat(sprintf("Gradient calculation complete: %d valid, %d missing (edge points)\n",
                n_valid_gradients, n_missing))
    
    if (n_valid_gradients > 0) {
      mag_range <- range(data_with_gradients$gradient_mag, na.rm = TRUE)
      cat(sprintf("Gradient magnitude range: %.4f to %.4f\n", 
                  mag_range[1], mag_range[2]))
    }
  }
  
  return(data_with_gradients)
}

#' @title Calculate Variable Limits for Consistent Scales
#' @description
#' Determines the minimum and maximum values for each variable to ensure
#' consistent color scales across scenarios. Can calculate either shared
#' limits (all variables use same scale) or independent limits (each variable
#' has its own scale).
#' 
#' Shared scales are useful when comparing magnitude across variables (e.g.,
#' different cost components). Independent scales maximize color resolution
#' for each variable separately.
#' 
#' @param data Data frame containing delayed deployment results
#' @param variables Character vector of variable names to calculate limits for
#' @param shared_scale Logical indicating whether to use shared scale across
#'   all variables (TRUE) or independent scales (FALSE). Default is FALSE.
#' @param use_scale_limits Logical indicating whether to cap the colour scale at a
#'   percentile threshold (TRUE) or use the full data range (FALSE). Default is FALSE.
#' @param scale_limit_percentile Numeric percentile (0-100) at which to cap the
#'   colour scale when use_scale_limits = TRUE. Default is 95.
#' @param verbose Logical indicating whether to print limit information (default: TRUE)
#' 
#' @return Named list of limits:
#'   - If shared_scale = TRUE: Single element named "shared" with c(min, max)
#'   - If shared_scale = FALSE: One element per variable, each with c(min, max)
#' 
#' @details
#' When shared_scale = TRUE, calculates the overall minimum and maximum across
#' all specified variables. This is particularly useful for multi-row dashboards
#' where direct visual comparison of magnitudes is important (e.g., comparing
#' total costs, abatement costs, and damage costs on the same scale).
#' 
#' @examples
#' # Independent scales (default)
#' limits <- calculate_variable_limits(
#'   plot_data, 
#'   variables = c("peak_temperature", "years_above_1p5")
#' )
#' # Returns: list(peak_temperature = c(min, max), years_above_1p5 = c(min, max))
#' 
#' # Shared scale for cost comparison
#' limits <- calculate_variable_limits(
#'   plot_data,
#'   variables = c("total_cost", "abatement_cost", "temp_cost"),
#'   shared_scale = TRUE
#' )
#' # Returns: list(shared = c(min, max))
calculate_variable_limits <- function(data,
                                      variables,
                                      shared_scale = FALSE,
                                      use_scale_limits = FALSE,
                                      scale_limit_percentile = 95,
                                      verbose = TRUE) {
  
  # Validate inputs
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("data must be a non-empty data frame")
  }
  
  if (length(variables) == 0) {
    stop("variables must contain at least one variable name")
  }
  
  # Check all variables exist
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  # Calculate limits
  if (shared_scale) {
    # Calculate shared limits across all variables
    all_values <- unlist(data[, variables], use.names = FALSE)
    all_values <- all_values[!is.na(all_values)]
    
    if (length(all_values) == 0) {
      stop("No non-missing values found in specified variables")
    }
    
    # Calculate base limits
    min_val <- min(all_values)
    max_val <- max(all_values)
    
    # Apply percentile capping if requested
    if (use_scale_limits) {
      max_val <- quantile(all_values, probs = scale_limit_percentile / 100)
      
      if (verbose) {
        actual_max <- max(all_values)
        cat(sprintf("Shared scale with %g%% percentile cap:\n", scale_limit_percentile))
        cat(sprintf("  Full range: [%.4f, %.4f]\n", min_val, actual_max))
        cat(sprintf("  Capped range: [%.4f, %.4f]\n", min_val, max_val))
        cat(sprintf("  Values above %.4f will be displayed as maximum colour\n", max_val))
      }
    } else if (verbose) {
      cat(sprintf("Shared scale across %d variables: [%.4f, %.4f]\n",
                  length(variables), min_val, max_val))
    }
    
    shared_limits <- c(min_val, max_val)
    limits_list <- list(shared = shared_limits)
    
  } else {
    # Calculate independent limits for each variable
    limits_list <- list()
    
    for (var in variables) {
      var_values <- data[[var]][!is.na(data[[var]])]
      
      if (length(var_values) == 0) {
        warning("No non-missing values found for variable: ", var)
        next
      }
      
      # Calculate base limits
      min_val <- min(var_values)
      max_val <- max(var_values)
      
      # Apply percentile capping if requested
      if (use_scale_limits) {
        max_val <- quantile(var_values, probs = scale_limit_percentile / 100)
        
        if (verbose) {
          actual_max <- max(var_values)
          cat(sprintf("  %s with %g%% cap: [%.4f, %.4f] (full: [%.4f, %.4f])\n", 
                      var, scale_limit_percentile, min_val, max_val, min_val, actual_max))
        }
      } else if (verbose) {
        cat(sprintf("  %s: [%.4f, %.4f]\n", var, min_val, max_val))
      }
      
      var_limits <- c(min_val, max_val)
      limits_list[[var]] <- var_limits
    }
    
    if (verbose && length(limits_list) > 0) {
      scale_type <- if (use_scale_limits) {
        sprintf("percentile-capped (%g%%)", scale_limit_percentile)
      } else {
        "independent"
      }
      cat(sprintf("Calculated %s scales for %d variables\n", 
                  scale_type, length(limits_list)))
    }
  }
  
  # Validate that we got limits
  if (length(limits_list) == 0) {
    stop("Failed to calculate limits for any variables")
  }
  
  return(limits_list)
}

# ==============================================================================
# Section 2: Plot Component Helpers
# ==============================================================================
#
# These functions create and customize individual plot components used in
# delayed deployment heatmaps. They handle theming, color palettes, labels,
# contour lines, arrows, and infeasible markers. All functions return either
# plot elements or modified ggplot objects that can be further customized.
#
# Functions in this section:
# - get_delayed_deployment_theme(): Return consistent ggplot2 theme
# - get_variable_palette(): Return appropriate color palette for variable
# - get_variable_label(): Return formatted legend label for variable
# - calculate_contour_breaks(): Determine intelligent contour break points
# - create_delayed_deployment_base_heatmap(): Create core geom_tile plot
# - add_delayed_deployment_contours(): Add white contour lines
# - add_delayed_deployment_arrows(): Add gradient vector field arrows
# - add_delayed_deployment_infeasible_markers(): Add red X markers
# ==============================================================================

#' @title Get Delayed Deployment Theme
#' @description
#' Returns a consistent ggplot2 theme for delayed deployment heatmaps. Theme
#' parameters are adjusted based on whether the plot is part of a multi-row
#' dashboard (smaller text for compact layouts) or a single-row dashboard
#' (larger text for readability).
#' 
#' The theme removes grid lines for cleaner heatmaps and uses theme_bw() as
#' a base for publication-quality output.
#' 
#' @param multi_row Logical indicating whether this is a multi-row layout.
#'   If TRUE, uses smaller text sizes. If FALSE, uses larger text sizes.
#'   Default is FALSE (single-row layout).
#' 
#' @return A ggplot2 theme object that can be added to plots with `+`
#' 
#' @details
#' Text sizes are scaled based on layout:
#' - Single-row (multi_row = FALSE): Base size 10, suitable for 2-row layouts
#' - Multi-row (multi_row = TRUE): Base size 8, suitable for 3+ row layouts
#' 
#' All plots have legend.position set to "none" by default, as legends are
#' typically extracted and positioned separately in dashboard layouts.
#' 
#' @examples
#' # For single-row dashboard
#' p <- ggplot(data, aes(x, y)) + 
#'   geom_tile() + 
#'   get_delayed_deployment_theme()
#' 
#' # For multi-row dashboard
#' p <- ggplot(data, aes(x, y)) + 
#'   geom_tile() + 
#'   get_delayed_deployment_theme(multi_row = TRUE)
get_delayed_deployment_theme <- function(multi_row = FALSE) {
  
  # Set text sizes based on layout type
  if (multi_row) {
    # Smaller text for multi-row layouts (3+ rows)
    base_size <- 8
    title_size <- 9
    axis_title_size <- 7
    axis_text_size <- 6
    legend_title_size <- 6
    legend_text_size <- 5.5
    legend_key_size <- 0.22
    margin_size <- 2
  } else {
    # Larger text for single-row layouts (1-2 rows)
    base_size <- 10
    title_size <- 10
    axis_title_size <- 9
    axis_text_size <- 8
    legend_title_size <- 8
    legend_text_size <- 7
    legend_key_size <- 0.35
    margin_size <- 3
  }
  
  # Build theme
  theme_object <- BASE_DELAYED_DEPLOYMENT_THEME +
    theme(
      text = element_text(size = base_size),
      plot.title = element_text(size = title_size, hjust = 0.5),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(legend_key_size, "cm"),
      legend.position = "none",  # Legends extracted separately
      plot.margin = margin(1, margin_size, 1, margin_size)  # Reduce top/bottom from margin_size to 1
    )
  
  return(theme_object)
}

#' @title Get Variable Color Palette
#' @description
#' Returns the appropriate viridis color palette information for a given
#' variable. Uses default palettes defined in DELAYED_DEPLOYMENT_PALETTES
#' but can be overridden with custom palette specifications.
#' 
#' Default palette assignments:
#' - Temperature variables (peak_temperature): plasma (reversed)
#' - Years above threshold (years_above_1p5): plasma (reversed)
#' - Cost variables (all cost types): viridis (reversed)
#' 
#' @param variable Character string naming the variable
#' @param custom_palettes Optional named list of custom palette specifications.
#'   Each element should be a list with 'option' and 'direction' components.
#'   If provided, overrides default palettes. Default is NULL.
#' 
#' @return List with two elements:
#'   - option: Character string of viridis palette name 
#'     (e.g., "plasma", "viridis", "magma", "inferno", "cividis")
#'   - direction: Numeric value of 1 (forward) or -1 (reversed)
#' 
#' @details
#' The viridis family of color palettes are perceptually uniform, colorblind-
#' friendly, and print well in grayscale. Direction = -1 (reversed) is used
#' for most variables to map darker colors to higher values, which is
#' conventional for cost and temperature visualizations.
#' 
#' @examples
#' # Get default palette for temperature
#' palette_info <- get_variable_palette("peak_temperature")
#' # Returns: list(option = "plasma", direction = -1)
#' 
#' # Get default palette for costs
#' palette_info <- get_variable_palette("abatement_cost")
#' # Returns: list(option = "viridis", direction = -1)
#' 
#' # Override with custom palette
#' custom <- list(peak_temperature = list(option = "magma", direction = 1))
#' palette_info <- get_variable_palette("peak_temperature", custom)
#' # Returns: list(option = "magma", direction = 1)
get_variable_palette <- function(variable, custom_palettes = NULL) {
  
  # Check for custom palette first
  if (!is.null(custom_palettes) && variable %in% names(custom_palettes)) {
    palette_info <- custom_palettes[[variable]]
    
    # Validate custom palette structure
    if (!is.list(palette_info) || 
        !"option" %in% names(palette_info) || 
        !"direction" %in% names(palette_info)) {
      stop("Custom palette for '", variable, "' must be a list with 'option' and 'direction' elements")
    }
    
    # Validate direction value
    if (!palette_info$direction %in% c(1, -1)) {
      stop("Palette direction must be 1 or -1, got: ", palette_info$direction)
    }
    
    return(palette_info)
  }
  
  # Use default palette
  if (variable %in% names(DELAYED_DEPLOYMENT_PALETTES)) {
    return(DELAYED_DEPLOYMENT_PALETTES[[variable]])
  }
  
  # If variable not found, issue warning and return default viridis
  warning("No palette defined for variable '", variable, "'. Using default viridis (reversed).")
  return(list(option = "viridis", direction = -1))
}

#' @title Get Variable Legend Label
#' @description
#' Returns the appropriate formatted label for a variable's legend. Uses
#' default labels defined in DELAYED_DEPLOYMENT_LABELS but can be overridden
#' with custom label specifications.
#' 
#' Labels use newline characters (\n) to create multi-line legends that fit
#' better in narrow legend spaces. All abbreviations use full stops for
#' scientific publication standards.
#' 
#' @param variable Character string naming the variable
#' @param custom_labels Optional named list of custom labels.
#'   Each element should be a character string with legend label.
#'   If provided, overrides default labels. Default is NULL.
#' 
#' @return Character string with formatted legend label (may contain \n for line breaks)
#' 
#' @details
#' Default labels follow scientific conventions:
#' - Use proper units in parentheses (e.g., "°C", "$ trillion")
#' - Abbreviations include full stops (e.g., "Temp." not "Temp")
#' - Multi-line format for space efficiency
#' 
#' @examples
#' # Get default label for temperature
#' label <- get_variable_label("peak_temperature")
#' # Returns: "Peak\nTemp.\n(°C)"
#' 
#' # Get default label for costs
#' label <- get_variable_label("abatement_cost")
#' # Returns: "Abatement\nCost\n($ trillion)"
#' 
#' # Override with custom label
#' custom <- list(peak_temperature = "Peak\nTemperature\n(°C)")
#' label <- get_variable_label("peak_temperature", custom)
#' # Returns: "Peak\nTemperature\n(°C)"
get_variable_label <- function(variable, custom_labels = NULL, single_line = FALSE) {
  
  # Check for custom label first
  if (!is.null(custom_labels) && variable %in% names(custom_labels)) {
    label <- custom_labels[[variable]]
    
    # Validate custom label is character
    if (!is.character(label) || length(label) != 1) {
      stop("Custom label for '", variable, "' must be a single character string")
    }
    
    return(label)
  }
  
  # Choose appropriate label set based on single_line parameter
  label_set <- if (single_line) {
    DELAYED_DEPLOYMENT_LABELS_SINGLE_LINE
  } else {
    DELAYED_DEPLOYMENT_LABELS
  }
  
  # Use default label
  if (variable %in% names(label_set)) {
    return(label_set[[variable]])
  }
  
  # If variable not found, issue warning and return variable name as label
  warning("No label defined for variable '", variable, "'. Using variable name as label.")
  
  # Create basic label from variable name (replace underscores with spaces, capitalize)
  basic_label <- gsub("_", " ", variable)
  basic_label <- tools::toTitleCase(basic_label)
  
  return(basic_label)
}

#' @title Calculate Contour Breaks
#' @description
#' Determines intelligent contour line break points for a variable based on
#' its data range. Can use automatic break calculation with adaptive intervals
#' or accept manual break specifications.
#' 
#' Automatic calculation uses adaptive intervals based on the data range to
#' ensure contours are neither too sparse nor too dense. Larger ranges use
#' larger intervals (e.g., 100 for ranges > 500) while smaller ranges use
#' finer intervals (e.g., 1 for ranges < 10).
#' 
#' @param data Data frame containing the variable
#' @param variable Character string naming the variable for contour calculation
#' @param breaks Either "auto" for automatic calculation or a numeric vector
#'   of specific break values. Default is "auto".
#' @param verbose Logical indicating whether to print break information (default: FALSE)
#' 
#' @return Numeric vector of contour break values
#' 
#' @details
#' Automatic break intervals are determined by data range:
#' - Range > 500: interval = 100
#' - Range > 200: interval = 50
#' - Range > 100: interval = 25
#' - Range > 50: interval = 10
#' - Range > 20: interval = 5
#' - Range > 10: interval = 2
#' - Range ≤ 10: interval = 1
#' 
#' For integer-valued variables (e.g., years above 1.5°C), automatic mode
#' uses integer breaks spaced by 1 year.
#' 
#' @examples
#' # Automatic breaks for temperature
#' breaks <- calculate_contour_breaks(data, "peak_temperature")
#' 
#' # Manual breaks for specific contour lines
#' breaks <- calculate_contour_breaks(
#'   data, 
#'   "abatement_cost", 
#'   breaks = c(100, 200, 300, 400, 500)
#' )
#' 
#' # Automatic breaks for years (integer variable)
#' breaks <- calculate_contour_breaks(data, "years_above_1p5")
calculate_contour_breaks <- function(data,
                                     variable,
                                     breaks = "auto",
                                     verbose = FALSE) {
  
  # Validate inputs
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("data must be a non-empty data frame")
  }
  
  if (!variable %in% names(data)) {
    stop("Variable '", variable, "' not found in data")
  }
  
  # If manual breaks provided, validate and return
  if (!identical(breaks, "auto")) {
    if (!is.numeric(breaks) || length(breaks) < 2) {
      stop("Manual breaks must be a numeric vector with at least 2 values")
    }
    
    if (verbose) {
      cat(sprintf("Using manual breaks for %s: %s\n", 
                  variable, paste(breaks, collapse = ", ")))
    }
    
    return(breaks)
  }
  
  # Get variable values
  var_values <- data[[variable]][!is.na(data[[variable]])]
  
  if (length(var_values) == 0) {
    stop("No non-missing values found for variable: ", variable)
  }
  
  var_range <- range(var_values)
  
  # Strategy 1: For peak_temperature, return NULL to let geom_contour() auto-calculate
  if (variable == "peak_temperature") {
    if (verbose) {
      cat(sprintf("Using automatic breaks for %s (range: %.4f to %.4f)\n",
                  variable, var_range[1], var_range[2]))
    }
    return(NULL)
  }
  
  # Strategy 2: For years_above_1p5, use interval of 5 years
  if (variable == "years_above_1p5") {
    interval <- 5
    break_min <- ceiling(var_range[1] / interval) * interval
    break_max <- floor(var_range[2] / interval) * interval
    
    # Generate breaks
    if (break_min > break_max) {
      contour_breaks <- c(var_range[1], var_range[2])
    } else {
      contour_breaks <- seq(break_min, break_max, by = interval)
      
      # Ensure we have at least 2 breaks
      if (length(contour_breaks) < 2) {
        contour_breaks <- c(var_range[1], var_range[2])
      }
    }
    
    if (verbose) {
      cat(sprintf("Using 5-year interval breaks for %s (range: %.4f to %.4f): %s\n",
                  variable, var_range[1], var_range[2],
                  paste(contour_breaks, collapse = ", ")))
    }
    
    return(contour_breaks)
  }
  
  # Strategy 3: For cost variables, use adaptive interval logic
  range_size <- var_range[2] - var_range[1]
  
  # Determine appropriate interval based on range size
  if (range_size > 500) {
    interval <- 100
  } else if (range_size > 200) {
    interval <- 50
  } else if (range_size > 100) {
    interval <- 25
  } else if (range_size > 50) {
    interval <- 10
  } else if (range_size > 20) {
    interval <- 5
  } else if (range_size > 10) {
    interval <- 2
  } else {
    interval <- 1
  }
  
  # Calculate break sequence
  # Round min/max to nearest interval for clean breaks
  break_min <- ceiling(var_range[1] / interval) * interval
  break_max <- floor(var_range[2] / interval) * interval
  
  # Ensure break_min <= break_max
  if (break_min > break_max) {
    # Range is smaller than interval, use variable range directly
    contour_breaks <- c(var_range[1], var_range[2])
    
    if (verbose) {
      cat(sprintf("Range (%.4f to %.4f) smaller than interval (%d). Using min/max as breaks for %s\n",
                  var_range[1], var_range[2], interval, variable))
    }
  } else if (break_min == break_max) {
    # Only one break point, use surrounding values
    contour_breaks <- c(var_range[1], break_min, var_range[2])
    
    if (verbose) {
      cat(sprintf("Single interval break at %.4f for %s\n", break_min, variable))
    }
  } else {
    # Normal case: generate sequence
    contour_breaks <- seq(break_min, break_max, by = interval)
    
    # Ensure we have at least 2 breaks
    if (length(contour_breaks) < 2) {
      contour_breaks <- c(var_range[1], var_range[2])
      
      if (verbose) {
        cat(sprintf("Insufficient breaks generated. Using min/max for %s\n", variable))
      }
    }
  }
  
  if (verbose) {
    cat(sprintf("Final breaks for %s (range: %.4f to %.4f, interval: %d): %s\n",
                variable, var_range[1], var_range[2], interval,
                paste(sprintf("%.4f", contour_breaks), collapse = ", ")))
  }
  
  return(contour_breaks)
}

#' @title Create Base Heatmap
#' @description
#' Creates the foundational heatmap plot using geom_tile() with appropriate
#' color scale, coordinate system, and labels. This is the base layer to which
#' other elements (contours, arrows, markers) can be added.
#' 
#' Uses equal coordinate scaling (coord_equal()) to ensure tiles are square,
#' which is essential for accurate visual interpretation of the delay space.
#' 
#' @param scenario_data Data frame containing data for a single scenario
#'   Must have columns: mitigation_delay, cdr_delay, feasible, and the variable
#' @param variable Character string naming the variable to plot
#' @param variable_limits Numeric vector of length 2: c(min, max) for color scale
#' @param palette_info List with 'option' and 'direction' for viridis palette
#' @param variable_label Character string for legend title (supports \n line breaks)
#' @param scenario_name Character string with scenario name for plot title
#' @param show_title Logical indicating whether to show plot title (default: TRUE)
#' @param title_text Optional custom title text. If NULL and show_title = FALSE,
#'   uses scenario_name. Default is NULL.
#' @param x_label Character string for x-axis label (default: "Mitigation Deployment Delay (years)")
#' @param y_label Character string for y-axis label (default: "CDR Deployment Delay (years)")
#' @param show_x_label Logical indicating whether to show x-axis label (default: TRUE)
#' @param show_y_label Logical indicating whether to show y-axis label (default: TRUE)
#' @param theme_object ggplot2 theme object to apply (default: get_delayed_deployment_theme())
#' 
#' @return ggplot object with base heatmap (geom_tile + scales + coord_equal)
#' 
#' @details
#' The function creates a minimal base heatmap that can be extended with
#' additional layers. Legend position is controlled by the theme_object and
#' is typically set to "none" since legends are extracted separately.
#' 
#' @examples
#' # Basic heatmap for one scenario
#' p <- create_delayed_deployment_base_heatmap(
#'   scenario_data = data %>% filter(scenario_short == "SSP1"),
#'   variable = "peak_temperature",
#'   variable_limits = c(1.5, 2.5),
#'   palette_info = list(option = "plasma", direction = -1),
#'   variable_label = "Peak\nTemp.\n(°C)",
#'   scenario_name = "SSP1"
#' )
#' 
#' # Heatmap without axis labels (for multi-plot dashboard)
#' p <- create_delayed_deployment_base_heatmap(
#'   scenario_data = data %>% filter(scenario_short == "SSP2"),
#'   variable = "abatement_cost",
#'   variable_limits = c(0, 500),
#'   palette_info = list(option = "viridis", direction = -1),
#'   variable_label = "Abatement\nCost\n($ trillion)",
#'   scenario_name = "SSP2",
#'   show_x_label = FALSE,
#'   show_y_label = FALSE
#' )
create_delayed_deployment_base_heatmap <- function(scenario_data,
                                                   variable,
                                                   variable_limits,
                                                   palette_info,
                                                   variable_label,
                                                   scenario_name,
                                                   show_title = FALSE,
                                                   title_text = NULL,
                                                   x_label = "Mitigation Deployment Delay (years)",
                                                   y_label = "CDR Deployment Delay (years)",
                                                   show_x_label = TRUE,
                                                   show_y_label = TRUE,
                                                   theme_object = get_delayed_deployment_theme()) {
  
  # Validate inputs
  if (!is.data.frame(scenario_data) || nrow(scenario_data) == 0) {
    stop("scenario_data must be a non-empty data frame")
  }
  
  required_cols <- c("mitigation_delay", "cdr_delay", variable)
  missing_cols <- setdiff(required_cols, names(scenario_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (length(variable_limits) != 2) {
    stop("variable_limits must be a numeric vector of length 2: c(min, max)")
  }
  
  # Determine title
  plot_title <- if (show_title) {
    if (!is.null(title_text)) title_text else scenario_name
  } else {
    ""
  }
  
  # Determine axis labels
  x_axis_label <- if (show_x_label) x_label else ""
  y_axis_label <- if (show_y_label) y_label else ""
  
  # Create base plot
  p <- ggplot(scenario_data, aes(x = mitigation_delay, y = cdr_delay, fill = .data[[variable]])) +
    geom_tile(color = NA) +  # Add color = NA to remove tile borders
    scale_fill_viridis_c(
      name = variable_label,
      option = palette_info$option,
      direction = palette_info$direction,
      limits = variable_limits,
      oob = scales::squish,  # Squish out-of-bounds values to max colour (not grey)
      labels = scales::comma
    ) +
    labs(
      title = plot_title,
      x = x_axis_label,
      y = y_axis_label
    ) +
    coord_equal() +  # Square tiles for accurate spatial interpretation
    theme_object
  
  return(p)
}

#' @title Add Contour Lines to Heatmap
#' @description
#' Adds white contour lines to an existing heatmap plot. Contours help
#' visualize gradients and patterns in the data, making it easier to identify
#' regions of similar values and steep gradients.
#' 
#' Contours are rendered in white with adjustable transparency to ensure they
#' are visible across the full range of the color scale without obscuring
#' underlying data patterns.
#' 
#' @param plot_object ggplot object to add contours to (typically from
#'   create_delayed_deployment_base_heatmap())
#' @param scenario_data Data frame containing data for the scenario
#'   Must have columns: mitigation_delay, cdr_delay, and the variable
#' @param variable Character string naming the variable for contour calculation
#' @param contour_breaks Numeric vector of values at which to draw contours
#' @param contour_color Color for contour lines (default: "white")
#' @param contour_alpha Transparency of contour lines, 0 (invisible) to 1 (opaque).
#'   Default is 0.6
#' @param contour_linewidth Width of contour lines in mm (default: 0.3)
#' 
#' @return ggplot object with contour lines added
#' 
#' @details
#' Contour lines are drawn using geom_contour() which interpolates values
#' between grid points to create smooth isolines. The function requires
#' regularly spaced data points for accurate contour calculation.
#' 
#' White contours with moderate transparency (alpha = 0.6) provide good
#' visibility across both light and dark regions of the color scale.
#' 
#' @examples
#' # Add contours with default settings
#' p <- create_delayed_deployment_base_heatmap(...)
#' p <- add_delayed_deployment_contours(
#'   p, 
#'   scenario_data, 
#'   variable = "peak_temperature",
#'   contour_breaks = c(1.6, 1.8, 2.0, 2.2)
#' )
#' 
#' # Add contours with custom transparency
#' p <- add_delayed_deployment_contours(
#'   p,
#'   scenario_data,
#'   variable = "abatement_cost",
#'   contour_breaks = seq(100, 500, by = 100),
#'   contour_alpha = 0.8
#' )
add_delayed_deployment_contours <- function(plot_object,
                                            scenario_data,
                                            variable,
                                            contour_breaks,
                                            contour_color = "white",
                                            contour_alpha = 0.6,
                                            contour_linewidth = 0.3) {
  
  # Validate inputs
  if (!inherits(plot_object, "gg")) {
    stop("plot_object must be a ggplot object")
  }
  
  if (!is.data.frame(scenario_data) || nrow(scenario_data) == 0) {
    stop("scenario_data must be a non-empty data frame")
  }
  
  if (!variable %in% names(scenario_data)) {
    stop("Variable '", variable, "' not found in scenario_data")
  }
  
  if (contour_alpha < 0 || contour_alpha > 1) {
    stop("contour_alpha must be between 0 and 1")
  }
  
  # Handle NULL contour_breaks (automatic calculation by geom_contour)
  if (is.null(contour_breaks)) {
    # Add contour layer WITHOUT breaks parameter - let geom_contour() auto-calculate
    plot_with_contours <- plot_object +
      geom_contour(
        data = scenario_data,
        aes(x = mitigation_delay, y = cdr_delay, z = .data[[variable]]),
        color = contour_color,
        alpha = contour_alpha,
        linewidth = contour_linewidth,
        inherit.aes = FALSE
      ) 
  } else {
    # Validate explicit breaks
    if (!is.numeric(contour_breaks) || length(contour_breaks) < 2) {
      stop("contour_breaks must be NULL or a numeric vector with at least 2 values")
    }
    
    # Add contour layer WITH explicit breaks parameter
    plot_with_contours <- plot_object +
      geom_contour(
        data = scenario_data,
        aes(x = mitigation_delay, y = cdr_delay, z = .data[[variable]]),
        breaks = contour_breaks,
        color = contour_color,
        alpha = contour_alpha,
        linewidth = contour_linewidth,
        inherit.aes = FALSE
      )
  }
  
  return(plot_with_contours)
}

#' @title Add Gradient Arrow Vector Field to Heatmap
#' @description
#' Adds gradient arrows to an existing heatmap plot showing the direction and
#' magnitude of steepest change in the response surface. Arrows are colored by
#' gradient magnitude to show both direction (arrow orientation) and rate of
#' change (arrow color).
#' 
#' Uses metR::geom_arrow() which provides proper arrow scaling and filtering.
#' Arrows can be thinned using skip parameters and filtered by minimum magnitude
#' to avoid cluttering the plot with negligible gradients.
#' 
#' @param plot_object ggplot object to add arrows to (typically from
#'   create_delayed_deployment_base_heatmap() with optional contours)
#' @param scenario_data Data frame containing data for the scenario
#'   Must have columns: mitigation_delay, cdr_delay, dx, dy, gradient_mag
#'   (from calculate_delayed_deployment_gradients())
#' @param arrow_scale Numeric scaling factor for arrow length. Larger values
#'   produce longer arrows. Default is 3.0
#' @param arrow_skip Integer specifying how many grid points to skip in each
#'   direction. skip=1 draws every arrow, skip=2 draws every other arrow, etc.
#'   Default is 1 (all arrows)
#' @param min_magnitude Minimum gradient magnitude threshold. Arrows with
#'   magnitude below this value are not drawn, reducing visual clutter.
#'   Default is 0 (show all arrows)
#' @param arrow_size Width of arrow lines in mm (default: 0.5)
#' @param arrow_alpha Transparency of arrows, 0 (invisible) to 1 (opaque).
#'   Default is 0.7
#' @param arrow_angle Angle of arrowhead in degrees (default: 20)
#' @param arrow_length Length of arrowhead in cm (default: 0.4)
#' @param mag_limits Numeric vector of length 2: c(min, max) for gradient
#'   magnitude color scale. Should match across scenarios for consistency.
#' @param show_mag_legend Logical indicating whether to show the magnitude
#'   color legend. Default is FALSE (legend extracted separately)
#' 
#' @return ggplot object with arrow vector field added
#' 
#' @details
#' The arrow vector field visualization requires pre-calculated gradients from
#' calculate_delayed_deployment_gradients(). Arrows show the direction of
#' steepest ascent (increase) in the variable.
#' 
#' Arrow parameters:
#' - preserve.dir = TRUE: Maintains correct direction regardless of aspect ratio
#' - arrow.type = "closed": Uses filled arrowheads for better visibility
#' - skip.x/skip.y: Thins arrows by specified factor in each direction
#' - min.mag: Filters small arrows to reduce clutter
#' 
#' @examples
#' # Add arrows with default settings
#' data_with_grad <- calculate_delayed_deployment_gradients(data, "abatement_cost")
#' p <- create_delayed_deployment_base_heatmap(...)
#' p <- add_delayed_deployment_arrows(
#'   p,
#'   scenario_data = data_with_grad,
#'   mag_limits = c(0, 50)
#' )
#' 
#' # Add filtered arrows (only strong gradients)
#' p <- add_delayed_deployment_arrows(
#'   p,
#'   scenario_data = data_with_grad,
#'   arrow_scale = 5.0,
#'   arrow_skip = 2,
#'   min_magnitude = 10,
#'   mag_limits = c(0, 50)
#' )
add_delayed_deployment_arrows <- function(plot_object,
                                          scenario_data,
                                          arrow_scale = 3.0,
                                          arrow_skip = 1,
                                          min_magnitude = 0,
                                          arrow_size = 0.5,
                                          arrow_alpha = 0.7,
                                          arrow_angle = 20,
                                          arrow_length = 0.4,
                                          mag_limits = NULL,
                                          show_mag_legend = FALSE) {
  
  # Validate inputs
  if (!inherits(plot_object, "gg")) {
    stop("plot_object must be a ggplot object")
  }
  
  if (!is.data.frame(scenario_data) || nrow(scenario_data) == 0) {
    stop("scenario_data must be a non-empty data frame")
  }
  
  required_cols <- c("mitigation_delay", "cdr_delay", "dx", "dy", "gradient_mag")
  missing_cols <- setdiff(required_cols, names(scenario_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns for arrows: ", paste(missing_cols, collapse = ", "),
         "\nDid you run calculate_delayed_deployment_gradients() first?")
  }
  
  # Filter to valid gradients only
  arrow_data <- scenario_data %>%
    filter(!is.na(dx), !is.na(dy))
  
  if (nrow(arrow_data) == 0) {
    warning("No valid gradient data found for arrows. Returning plot unchanged.")
    return(plot_object)
  }
  
  # Calculate magnitude limits if not provided
  if (is.null(mag_limits)) {
    mag_limits <- range(arrow_data$gradient_mag, na.rm = TRUE)
  }
  
  # Add arrow layer using metR
  plot_with_arrows <- plot_object +
    metR::geom_arrow(
      data = arrow_data,
      aes(x = mitigation_delay, 
          y = cdr_delay,
          dx = dx * arrow_scale, 
          dy = dy * arrow_scale,
          colour = gradient_mag),
      skip.x = arrow_skip,
      skip.y = arrow_skip,
      min.mag = min_magnitude,
      preserve.dir = TRUE,      # Keep direction correct regardless of aspect ratio
      arrow.angle = arrow_angle,
      arrow.length = arrow_length,
      arrow.type = "closed",    # Filled arrowhead
      linewidth = arrow_size,
      alpha = arrow_alpha,
      inherit.aes = FALSE
    ) +
    scale_colour_viridis_c(
      name = "Gradient\nMagnitude",
      option = "plasma",
      limits = mag_limits,
      guide = if (show_mag_legend) "colorbar" else "none"
    ) +
    metR::scale_mag(guide = "none")  # Required for dx/dy arrows, suppress legend
  
  return(plot_with_arrows)
}

#' @title Add Infeasible Combination Markers to Heatmap
#' @description
#' Adds red X markers to indicate parameter combinations that are infeasible
#' (e.g., where the optimization failed to converge or constraints could not
#' be satisfied). This provides visual feedback about which regions of the
#' delay space are problematic.
#' 
#' Markers use shape = 4 (X) in red to stand out against the heatmap colors
#' while remaining subtle enough not to dominate the visualization.
#' 
#' @param plot_object ggplot object to add markers to
#' @param scenario_data Data frame containing data for the scenario
#'   Must have columns: mitigation_delay, cdr_delay, feasible (logical)
#' @param marker_color Color for infeasible markers (default: "red")
#' @param marker_size Size of marker in mm (default: 1)
#' @param marker_alpha Transparency of markers, 0 (invisible) to 1 (opaque).
#'   Default is 0.8
#' @param marker_stroke Width of marker lines in mm (default: 0.8)
#' 
#' @return ggplot object with infeasible markers added (or unchanged if all feasible)
#' 
#' @details
#' Markers are only added if there are infeasible combinations in the data.
#' If all combinations are feasible, the function returns the plot unchanged
#' (no warning issued as this is expected behavior for some scenarios).
#' 
#' The 'feasible' column should be a logical vector where FALSE indicates
#' infeasible combinations. Markers use inherit.aes = FALSE to prevent
#' interference with the heatmap's fill aesthetic.
#' 
#' @examples
#' # Add markers with default settings
#' p <- create_delayed_deployment_base_heatmap(...)
#' p <- add_delayed_deployment_infeasible_markers(p, scenario_data)
#' 
#' # Add markers with custom appearance
#' p <- add_delayed_deployment_infeasible_markers(
#'   p,
#'   scenario_data,
#'   marker_size = 2,
#'   marker_alpha = 1.0,
#'   marker_stroke = 1.5
#' )
add_delayed_deployment_infeasible_markers <- function(plot_object,
                                                      scenario_data,
                                                      marker_color = "red",
                                                      marker_size = 1,
                                                      marker_alpha = 0.8,
                                                      marker_stroke = 0.8) {
  
  # Validate inputs
  if (!inherits(plot_object, "gg")) {
    stop("plot_object must be a ggplot object")
  }
  
  if (!is.data.frame(scenario_data) || nrow(scenario_data) == 0) {
    stop("scenario_data must be a non-empty data frame")
  }
  
  required_cols <- c("mitigation_delay", "cdr_delay", "feasible")
  missing_cols <- setdiff(required_cols, names(scenario_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.logical(scenario_data$feasible)) {
    stop("'feasible' column must be logical (TRUE/FALSE)")
  }
  
  if (marker_alpha < 0 || marker_alpha > 1) {
    stop("marker_alpha must be between 0 and 1")
  }
  
  # Filter to infeasible combinations
  infeasible_data <- scenario_data %>%
    filter(!feasible)
  
  # Only add markers if there are infeasible combinations
  if (nrow(infeasible_data) == 0) {
    # All combinations feasible - return plot unchanged
    return(plot_object)
  }
  
  # Add infeasible markers
  plot_with_markers <- plot_object +
    geom_point(
      data = infeasible_data,
      aes(x = mitigation_delay, y = cdr_delay),
      shape = 4,              # X marker
      size = marker_size,
      color = marker_color,
      alpha = marker_alpha,
      stroke = marker_stroke,
      inherit.aes = FALSE     # Don't inherit fill aesthetic from heatmap
    )
  
  return(plot_with_markers)
}

# ==============================================================================
# Section 3: Single Plot Creation
# ==============================================================================
#
# This section contains the main workhorse function that orchestrates the
# creation of a complete heatmap plot for a single scenario. It combines all
# the helper functions from Section 2 to build a publication-ready plot with
# the requested features (contours, arrows, markers, etc.).
#
# The function handles all the conditional logic for what layers to include,
# manages axis labels for dashboard positioning, and ensures consistent
# appearance across all scenarios.
#
# Functions in this section:
# - create_single_delayed_deployment_plot(): Build complete plot for one scenario
# ==============================================================================

#' @title Create Complete Plot for Single Scenario
#' @description
#' Orchestrates the creation of a complete heatmap plot for a single SSP
#' scenario by combining all plot component helpers. This is the main workhorse
#' function that builds publication-ready plots with all requested features.
#' 
#' The function intelligently handles:
#' - Conditional layer addition (contours, arrows, markers)
#' - Axis label positioning for dashboard layouts
#' - Legend visibility control
#' - Theme selection based on layout type
#' - Title formatting based on position
#' 
#' @param scenario_data Data frame containing data for a single scenario
#'   Must have columns: mitigation_delay, cdr_delay, feasible, and the variable
#'   If arrows requested, must also have: dx, dy, gradient_mag
#' @param scenario_name Character string with scenario name (e.g., "SSP1", "SSP2")
#' @param variable Character string naming the variable to plot
#' @param variable_limits Numeric vector c(min, max) for color scale
#' @param variable_label Character string for legend title
#' @param palette_info List with 'option' and 'direction' for viridis palette
#' @param add_contours Logical indicating whether to add contour lines (default: TRUE)
#' @param contour_breaks Either "auto" or numeric vector of contour break values
#' @param contour_alpha Transparency of contour lines (default: 0.6)
#' @param add_arrows Logical indicating whether to add gradient arrows (default: FALSE)
#' @param arrow_scale Numeric scaling factor for arrow length (default: 3.0)
#' @param arrow_skip Integer for arrow thinning (default: 1)
#' @param min_magnitude Minimum gradient magnitude for arrow display (default: 0)
#' @param arrow_size Width of arrow lines in mm (default: 0.5)
#' @param arrow_alpha Transparency of arrows (default: 0.7)
#' @param mag_limits Numeric vector c(min, max) for gradient magnitude scale
#'   Only used if add_arrows = TRUE. If NULL, calculated from data.
#' @param show_infeasible Logical indicating whether to mark infeasible combinations
#'   (default: TRUE)
#' @param show_legend Logical indicating whether to show legend (default: FALSE)
#' @param is_top_row Logical indicating if this is top row in dashboard.
#'   If TRUE, shows column header in title. Default is FALSE.
#' @param is_bottom_row Logical indicating if this is bottom row in dashboard.
#'   If TRUE, shows x-axis label. Default is FALSE.
#' @param is_leftmost_col Logical indicating if this is leftmost column in dashboard.
#'   If TRUE, shows y-axis label. Default is TRUE.
#' @param column_header Optional character string for column header (e.g., "Peak Temperature").
#'   Only used if is_top_row = TRUE. Default is NULL.
#' @param panel_label Optional character string for panel label (e.g., "A)", "B)").
#'   Prepended to scenario name in title. Default is NULL.
#' @param multi_row Logical indicating if part of multi-row dashboard (affects theme).
#'   Default is FALSE.
#' @param verbose Logical indicating whether to print progress messages (default: FALSE)
#' 
#' @return Complete ggplot object for the scenario
#' 
#' @details
#' The function builds plots in layers:
#' 1. Base heatmap (always)
#' 2. Contour lines (if requested)
#' 3. Gradient arrows (if requested)
#' 4. Infeasible markers (if requested)
#' 
#' Axis labels are controlled by position flags:
#' - Top row: Shows column header in title
#' - Bottom row: Shows x-axis label
#' - Leftmost column: Shows y-axis label
#' - Other positions: Minimal labels for compact dashboard
#' 
#' @examples
#' # Basic plot with contours
#' p <- create_single_delayed_deployment_plot(
#'   scenario_data = data %>% filter(scenario_short == "SSP1"),
#'   scenario_name = "SSP1",
#'   variable = "peak_temperature",
#'   variable_limits = c(1.5, 2.5),
#'   variable_label = "Peak\nTemp.\n(°C)",
#'   palette_info = list(option = "plasma", direction = -1),
#'   add_contours = TRUE,
#'   contour_breaks = "auto"
#' )
#' 
#' # Plot with arrows for cost analysis
#' p <- create_single_delayed_deployment_plot(
#'   scenario_data = data_with_gradients %>% filter(scenario_short == "SSP2"),
#'   scenario_name = "SSP2",
#'   variable = "abatement_cost",
#'   variable_limits = c(0, 500),
#'   variable_label = "Abatement\nCost\n($ trillion)",
#'   palette_info = list(option = "viridis", direction = -1),
#'   add_arrows = TRUE,
#'   arrow_scale = 5.0,
#'   min_magnitude = 10,
#'   mag_limits = c(0, 50)
#' )
create_single_delayed_deployment_plot <- function(scenario_data,
                                                  scenario_name,
                                                  variable,
                                                  variable_limits,
                                                  variable_label,
                                                  palette_info,
                                                  add_contours = TRUE,
                                                  contour_breaks = "auto",
                                                  contour_alpha = 0.6,
                                                  add_arrows = FALSE,
                                                  arrow_scale = 3.0,
                                                  arrow_skip = 1,
                                                  min_magnitude = 0,
                                                  arrow_size = 0.5,
                                                  arrow_alpha = 0.7,
                                                  mag_limits = NULL,
                                                  show_infeasible = TRUE,
                                                  show_legend = FALSE,
                                                  is_top_row = FALSE,
                                                  is_bottom_row = FALSE,
                                                  is_leftmost_col = TRUE,
                                                  column_header = NULL,
                                                  panel_label = NULL,
                                                  multi_row = FALSE,
                                                  verbose = FALSE) {
  
  # Validate scenario_data
  if (nrow(scenario_data) == 0) {
    # Return empty plot if no data
    if (verbose) {
      cat(sprintf("No data for scenario %s\n", scenario_name))
    }
    
    empty_plot <- ggplot() + 
      theme_void() + 
      labs(title = paste(scenario_name, "- No Data"))
    
    return(empty_plot)
  }
  
  if (verbose) {
    cat(sprintf("Creating plot for scenario: %s, variable: %s\n", 
                scenario_name, variable))
  }
  
  # Determine title text
  title_text <- ""
  
  if (!is.null(panel_label)) {
    title_text <- paste0(panel_label, " ", scenario_name)
  } else {
    title_text <- scenario_name
  }
  
  # Add column header if top row
  if (is_top_row && !is.null(column_header)) {
    title_text <- column_header
  }
  
  # Determine axis labels based on position
  show_x_label <- is_bottom_row
  show_y_label <- is_leftmost_col
  
  # Get appropriate theme
  theme_obj <- get_delayed_deployment_theme(multi_row = multi_row)
  
  # Adjust legend position if showing legend
  if (show_legend) {
    theme_obj <- theme_obj + theme(legend.position = "right")
  }
  
  # Step 1: Create base heatmap
  p <- create_delayed_deployment_base_heatmap(
    scenario_data = scenario_data,
    variable = variable,
    variable_limits = variable_limits,
    palette_info = palette_info,
    variable_label = variable_label,
    scenario_name = scenario_name,
    show_title = TRUE,
    title_text = title_text,
    show_x_label = show_x_label,
    show_y_label = show_y_label,
    theme_object = theme_obj
  )
  
  if (verbose) {
    cat("  - Base heatmap created\n")
  }
  
  # Step 2: Add contour lines if requested
  if (add_contours) {
    # Calculate contour breaks
    if (identical(contour_breaks, "auto")) {
      contour_breaks_calc <- calculate_contour_breaks(
        data = scenario_data,
        variable = variable,
        breaks = "auto",
        verbose = verbose
      )
    } else {
      contour_breaks_calc <- contour_breaks
    }
    
    # Add contours - handle both NULL (automatic) and explicit breaks
    if (is.null(contour_breaks_calc)) {
      # NULL means automatic calculation by geom_contour()
      p <- add_delayed_deployment_contours(
        plot_object = p,
        scenario_data = scenario_data,
        variable = variable,
        contour_breaks = NULL,
        contour_alpha = contour_alpha
      )
      
      if (verbose) {
        cat("  - Contours added (automatic breaks)\n")
      }
    } else if (length(contour_breaks_calc) >= 2) {
      # Explicit breaks provided
      p <- add_delayed_deployment_contours(
        plot_object = p,
        scenario_data = scenario_data,
        variable = variable,
        contour_breaks = contour_breaks_calc,
        contour_alpha = contour_alpha
      )
      
      if (verbose) {
        cat(sprintf("  - Contours added (%d breaks)\n", length(contour_breaks_calc)))
      }
    } else {
      if (verbose) {
        cat("  - Skipping contours (insufficient breaks)\n")
      }
    }
  }
  
  # Step 3: Add gradient arrows if requested
  if (add_arrows) {
    # Check that gradient columns exist
    required_gradient_cols <- c("dx", "dy", "gradient_mag")
    missing_gradient_cols <- setdiff(required_gradient_cols, names(scenario_data))
    
    if (length(missing_gradient_cols) > 0) {
      warning(sprintf("Cannot add arrows for %s: missing gradient columns: %s", 
                      scenario_name, paste(missing_gradient_cols, collapse = ", ")))
      
      if (verbose) {
        cat("  - Skipping arrows (missing gradient data)\n")
      }
    } else {
      # Filter to rows with valid gradient data
      arrow_data <- scenario_data %>%
        filter(!is.na(dx), !is.na(dy), !is.na(gradient_mag))
      
      if (nrow(arrow_data) > 0) {
        # Calculate magnitude limits if not provided
        if (is.null(mag_limits)) {
          mag_limits <- range(arrow_data$gradient_mag, na.rm = TRUE)
        }
        
        p <- add_delayed_deployment_arrows(
          plot_object = p,
          scenario_data = arrow_data,
          arrow_scale = arrow_scale,
          arrow_skip = arrow_skip,
          min_magnitude = min_magnitude,
          arrow_size = arrow_size,
          arrow_alpha = arrow_alpha,
          mag_limits = mag_limits,
          show_mag_legend = FALSE
        )
        
        if (verbose) {
          cat("  - Arrows added\n")
        }
      } else {
        if (verbose) {
          cat("  - Skipping arrows (no valid gradient data)\n")
        }
      }
    }
  }
  
  # Step 4: Add infeasible markers if requested
  if (show_infeasible) {
    infeasible_data <- scenario_data %>% filter(!feasible)
    
    if (nrow(infeasible_data) > 0) {
      p <- p + 
        geom_point(
          data = infeasible_data, 
          aes(x = mitigation_delay, y = cdr_delay), 
          shape = 4, 
          size = 1, 
          color = "red", 
          alpha = 0.8, 
          stroke = 0.8, 
          inherit.aes = FALSE
        )
      
      if (verbose) {
        cat(sprintf("  - Infeasible markers added (%d points)\n", nrow(infeasible_data)))
      }
    }
  }
  
  return(p)
}

# ==============================================================================
# Section 4: Dashboard Assembly
# ==============================================================================
#
# These functions handle the assembly of individual scenario plots into
# complete dashboards across multiple scenarios. They manage the creation of plot grids,
# legend extraction and positioning, layout arrangement using patchwork, and
# the addition of overall titles and annotations.
#
# The section supports both single-row layouts (1 variable, 5 scenarios in
# 3×2 grid) and multi-row layouts (multiple variables, each variable gets its
# own row of 5 scenarios).
#
# Functions in this section:
# - create_delayed_deployment_plot_grid(): Create all plots for all scenarios
# - extract_delayed_deployment_legend(): Extract legend for separate positioning
# - assemble_delayed_deployment_dashboard(): Arrange plots with patchwork
# - add_delayed_deployment_annotations(): Add title and summary statistics
# - save_delayed_deployment_dashboard(): Save plot to file
# ==============================================================================

#' @title Create Grid of Plots for All Scenarios and Variables
#' @description
#' Creates all individual scenario plots for the requested variables, organizing
#' them into a structured list ready for dashboard assembly. Handles panel
#' labeling (A, B, C, etc.), position-based axis label control, and ensures
#' consistent appearance across all plots.
#' 
#' For multi-variable dashboards, each variable gets its own row of scenarios.
#' For single-variable dashboards, scenarios are arranged in a 3×2 grid.
#' 
#' @param data Prepared data frame with all scenarios (from prepare_delayed_deployment_data())
#'   If arrows requested, must include gradient columns (from calculate_delayed_deployment_gradients())
#' @param variables Character vector of variable names to plot
#' @param variable_limits Named list of limits for each variable (or list with "shared" element)
#'   From calculate_variable_limits()
#' @param palette_info Named list of palette specifications for each variable
#' @param variable_labels Named list of legend labels for each variable
#' @param add_contours Logical indicating whether to add contour lines
#' @param contour_breaks Either "auto" or named list of breaks for each variable
#' @param contour_alpha Transparency of contour lines
#' @param add_arrows Logical indicating whether to add gradient arrows
#' @param arrow_scale Numeric scaling factor for arrow length
#' @param arrow_skip Integer for arrow thinning
#' @param min_magnitude Minimum gradient magnitude for arrow display
#' @param arrow_size Width of arrow lines in mm
#' @param arrow_alpha Transparency of arrows
#' @param mag_limits Numeric vector c(min, max) for gradient magnitude scale
#'   Only used if add_arrows = TRUE
#' @param show_infeasible Logical indicating whether to mark infeasible combinations
#' @param verbose Logical indicating whether to print progress messages
#' 
#' @return Named list with two levels:
#'   - First level: variable names
#'   - Second level: scenario names
#'   - Values: ggplot objects
#'   
#'   Structure: list(
#'     variable1 = list(SSP1 = plot, SSP2 = plot, ...),
#'     variable2 = list(SSP1 = plot, SSP2 = plot, ...)
#'   )
#' 
#' @details
#' Panel labeling follows alphabetical sequence:
#' - Single variable: A) SSP1, B) SSP2, C) SSP3, D) SSP4, E) SSP5
#' - Multiple variables: A-E for first variable row, F-J for second, K-O for third
#' 
#' Axis labels are positioned based on grid location:
#' - Top row: Shows column header (variable name)
#' - Bottom row: Shows x-axis label
#' - Leftmost column: Shows y-axis label
#' - Interior plots: Minimal labels
#' 
#' @examples
#' # Create grid for single variable
#' plot_grid <- create_delayed_deployment_plot_grid(
#'   data = prepared_data,
#'   variables = c("peak_temperature"),
#'   variable_limits = limits_list,
#'   palette_info = palette_list,
#'   variable_labels = label_list
#' )
#' 
#' # Create grid for multiple variables
#' plot_grid <- create_delayed_deployment_plot_grid(
#'   data = prepared_data,
#'   variables = c("total_cost", "abatement_cost", "temp_cost"),
#'   variable_limits = list(shared = c(0, 5000)),
#'   palette_info = palette_list,
#'   variable_labels = label_list,
#'   add_contours = TRUE
#' )
create_delayed_deployment_plot_grid <- function(data,
                                                variables,
                                                variable_limits,
                                                palette_info,
                                                variable_labels,
                                                add_contours = TRUE,
                                                contour_breaks = "auto",
                                                contour_alpha = 0.6,
                                                add_arrows = FALSE,
                                                arrow_scale = 3.0,
                                                arrow_skip = 1,
                                                min_magnitude = 0,
                                                arrow_size = 0.5,
                                                arrow_alpha = 0.7,
                                                mag_limits = NULL,
                                                show_infeasible = TRUE,
                                                verbose = TRUE) {
  
  # Validate inputs
  if (length(variables) == 0) {
    stop("variables must contain at least one variable name")
  }
  
  # Get available scenarios (using unique since scenario_short is character)
  available_scenarios <- unique(data$scenario_short)
  # Ensure they're in the correct order
  available_scenarios <- available_scenarios[available_scenarios %in% SSP_SCENARIO_ORDER]
  available_scenarios <- available_scenarios[order(match(available_scenarios, SSP_SCENARIO_ORDER))]
  
  n_scenarios <- length(available_scenarios)
  n_variables <- length(variables)
  
  if (verbose) {
    cat(sprintf("Creating plot grid: %d variable(s) × %d scenario(s) = %d plots\n",
                n_variables, n_scenarios, n_variables * n_scenarios))
  }
  
  # Determine if multi-row layout
  multi_row <- (n_variables > 1)
  
  # Initialize plot grid structure
  plot_grid <- list()
  
  # Panel label counter (A, B, C, ...)
  panel_counter <- 1
  
  # Loop through each variable
  for (var_idx in seq_along(variables)) {
    variable <- variables[var_idx]
    
    if (verbose) {
      cat(sprintf("\nCreating plots for variable: %s\n", variable))
    }
    
    # Get limits for this variable
    if ("shared" %in% names(variable_limits)) {
      var_limits <- variable_limits$shared
    } else if (variable %in% names(variable_limits)) {
      var_limits <- variable_limits[[variable]]
    } else {
      stop("No limits found for variable: ", variable)
    }
    
    # Get palette for this variable
    var_palette <- palette_info[[variable]]
    
    # Get label for this variable
    var_label <- variable_labels[[variable]]
    
    # Get contour breaks for this variable
    if (is.list(contour_breaks) && variable %in% names(contour_breaks)) {
      var_contour_breaks <- contour_breaks[[variable]]
    } else {
      var_contour_breaks <- contour_breaks  # Use default "auto" or single spec
    }
    
    # Determine position flags
    is_top_row <- (var_idx == 1)
    is_bottom_row <- (var_idx == n_variables)
    
    # For multi-row layouts, use variable name as column header for top row
    column_header <- if (multi_row && is_top_row) NULL else NULL
    
    # Initialize list for this variable
    plot_grid[[variable]] <- list()
    
    # Loop through each scenario
    for (scen_idx in seq_along(available_scenarios)) {
      scenario <- available_scenarios[scen_idx]
      
      # Filter data for this scenario using base R (dplyr filter fails in this context)
      scenario_data <- data[data$scenario_short == scenario, ]
      
      if (verbose) {
        cat(sprintf("  Filtered for %s: %d rows\n", scenario, nrow(scenario_data)))
      }
      
      # Determine if leftmost column
      is_leftmost <- (scen_idx == 1)
      
      # Create panel label (A, B, C, ...)
      panel_label <- paste0(LETTERS[panel_counter], ")")
      panel_counter <- panel_counter + 1
      
      if (verbose) {
        cat(sprintf("  Creating plot %s %s (%d/%d)\n", 
                    panel_label, scenario, scen_idx, n_scenarios))
      }
      
      # Create the plot
      p <- create_single_delayed_deployment_plot(
        scenario_data = scenario_data,
        scenario_name = scenario,
        variable = variable,
        variable_limits = var_limits,
        variable_label = var_label,
        palette_info = var_palette,
        add_contours = add_contours,
        contour_breaks = var_contour_breaks,
        contour_alpha = contour_alpha,
        add_arrows = add_arrows,
        arrow_scale = arrow_scale,
        arrow_skip = arrow_skip,
        min_magnitude = min_magnitude,
        arrow_size = arrow_size,
        arrow_alpha = arrow_alpha,
        mag_limits = mag_limits,
        show_infeasible = show_infeasible,
        show_legend = FALSE,  # Legends extracted separately
        is_top_row = is_top_row,
        is_bottom_row = is_bottom_row,
        is_leftmost_col = is_leftmost,
        column_header = column_header,
        panel_label = panel_label,
        multi_row = multi_row,
        verbose = FALSE  # Suppress individual plot verbosity
      )
      
      # Store plot
      plot_grid[[variable]][[scenario]] <- p
    }
  }
  
  if (verbose) {
    cat(sprintf("\nPlot grid complete: %d plots created\n", panel_counter - 1))
  }
  
  return(plot_grid)
}

#' @title Extract Legend from Plot
#' @description
#' Extracts the legend from a ggplot object for separate positioning in the
#' dashboard layout. Uses cowplot::get_legend() to convert the legend into a
#' grob (graphical object) that can be arranged with patchwork.
#' 
#' The function creates a temporary plot with the legend in the desired position,
#' extracts it, and returns it for inclusion in the final dashboard assembly.
#' 
#' @param plot_object ggplot object to extract legend from
#' @param legend_position Character string specifying legend position for extraction.
#'   Either "right" (for single-row layouts) or "bottom" (for multi-row layouts).
#'   Default is "right".
#' 
#' @return A grob (graphical object) containing the extracted legend that can
#'   be combined with plots using patchwork
#' 
#' @details
#' Legend positioning conventions:
#' - Single-row layouts (1-2 variables): Legend on right side
#' - Multi-row layouts (3+ variables): Legend on bottom, centered
#' 
#' The function requires the cowplot package for get_legend().
#' 
#' @examples
#' # Extract legend for right positioning (single-row)
#' legend <- extract_delayed_deployment_legend(
#'   plot_with_legend,
#'   legend_position = "right"
#' )
#' 
#' # Extract legend for bottom positioning (multi-row)
#' legend <- extract_delayed_deployment_legend(
#'   plot_with_legend,
#'   legend_position = "bottom"
#' )
extract_delayed_deployment_legend <- function(plot_object,
                                              legend_position = "right") {
  
  # Validate inputs
  if (!inherits(plot_object, "gg")) {
    stop("plot_object must be a ggplot object")
  }
  
  if (!legend_position %in% c("right", "bottom", "left", "top")) {
    stop("legend_position must be one of: 'right', 'bottom', 'left', 'top'")
  }
  
  # Check cowplot is available
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("cowplot package is required for legend extraction. Please install it.")
  }
  
  # Create temporary plot with legend in desired position
  # Add legend specifications with smaller sizes
  if (legend_position == "bottom") {
    plot_with_legend <- plot_object +
      theme(legend.position = legend_position) +
      guides(fill = guide_colorbar(
        barwidth = 15,      # Reduced from 25 to make narrower
        barheight = 0.4,    # Reduced from 1.2 to make shorter
        title.position = "top",
        title.hjust = 0.5,
        label.theme = element_text(size = 6),    # Smaller text
        title.theme = element_text(size = 6)     # Smaller title
      ))
  } else {
    plot_with_legend <- plot_object +
      theme(legend.position = legend_position) +
      guides(fill = guide_colorbar(
        barwidth = 1.5,     # Keep same for right legends
        barheight = 15,     # Keep same for right legends
        title.position = "top",
        title.hjust = 0.5,
        label.theme = element_text(size = 10),    # Larger text
        title.theme = element_text(size = 11)     # Larger title
      ))
  }
  
  # Extract legend using cowplot
  legend_grob <- cowplot::get_legend(plot_with_legend)
  
  return(legend_grob)
}

#' @title Assemble Dashboard from Plot Grid
#' @description
#' Arranges individual scenario plots into a complete dashboard using patchwork,
#' with legend positioned appropriately. Handles both single-row layouts
#' (3×2 grid with legend on right) and multi-row layouts (N×5 grid with
#' legend on bottom).
#' 
#' The function intelligently determines the layout based on the number of
#' variables and scenarios, adding blank plots where needed to complete the grid.
#' 
#' @param plot_grid Nested list of plots from create_delayed_deployment_plot_grid()
#'   Structure: list(variable1 = list(SSP1 = plot, ...), variable2 = list(...))
#' @param legend_grob Legend grob from extract_delayed_deployment_legend()
#' @param variables Character vector of variable names (in order)
#' @param scenarios Character vector of scenario names (in order)
#' @param legend_position Character string: "right" for single-row, "bottom" for multi-row
#'   Default is "right"
#' @param verbose Logical indicating whether to print assembly information
#' 
#' @return Combined patchwork object with all plots and legend
#' 
#' @details
#' Layout logic:
#' - Single variable (1 row): 3 columns × 2 rows, plots A-E, blank F, legend right
#' - Multiple variables (N rows): 5 columns × N rows, all plots filled, legend bottom
#' 
#' The function uses patchwork operators:
#' - `+` to combine plots horizontally
#' - `/` to combine rows vertically
#' - `|` to add legend beside main grid
#' - `plot_layout()` to control relative sizes
#' 
#' @examples
#' # Assemble single-row dashboard
#' dashboard <- assemble_delayed_deployment_dashboard(
#'   plot_grid = plot_grid,
#'   legend_grob = legend,
#'   variables = c("peak_temperature"),
#'   scenarios = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
#'   legend_position = "right"
#' )
#' 
#' # Assemble multi-row dashboard
#' dashboard <- assemble_delayed_deployment_dashboard(
#'   plot_grid = plot_grid,
#'   legend_grob = legend,
#'   variables = c("total_cost", "abatement_cost", "temp_cost"),
#'   scenarios = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"),
#'   legend_position = "bottom"
#' )
assemble_delayed_deployment_dashboard <- function(plot_grid,
                                                  legend_grob,
                                                  variables,
                                                  scenarios,
                                                  legend_position = "right",
                                                  verbose = FALSE) {
  
  # Check patchwork is available
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for dashboard assembly. Please install it.")
  }
  
  n_variables <- length(variables)
  n_scenarios <- length(scenarios)
  
  if (verbose) {
    cat(sprintf("Assembling dashboard: %d variable(s) × %d scenario(s)\n",
                n_variables, n_scenarios))
  }
  
  # Determine layout type
  if (n_variables == 1) {
    # Single-row layout: 3×2 grid with blank spot
    if (verbose) {
      cat("Layout: 3×2 grid (single variable) with legend on right\n")
    }
    
    variable <- variables[1]
    
    # Get plots for the single variable
    row1_plots <- list()
    row2_plots <- list()
    
    for (i in 1:min(3, n_scenarios)) {
      if (i <= n_scenarios) {
        row1_plots[[i]] <- plot_grid[[variable]][[scenarios[i]]]
      }
    }
    
    for (i in 4:min(6, n_scenarios)) {
      if (i <= n_scenarios) {
        row2_plots[[i-3]] <- plot_grid[[variable]][[scenarios[i]]]
      }
    }
    
    # Add blank plot if needed (position F for 5 scenarios)
    if (n_scenarios == 5) {
      blank_plot <- ggplot() + 
        theme_void() +
        theme(panel.background = element_rect(fill = "white", color = NA))
      row2_plots[[3]] <- blank_plot
    }
    
    # Combine rows
    row1_combined <- Reduce(`+`, row1_plots)
    row2_combined <- Reduce(`+`, row2_plots)
    
    # Combine rows vertically
    main_grid <- row1_combined / row2_combined
    
    # Add legend to right
    combined_plot <- main_grid | legend_grob
    combined_plot <- combined_plot + patchwork::plot_layout(widths = c(10, 1))
    
  } else {
    # Multi-row layout: N rows × 5 columns with legend on bottom
    if (verbose) {
      cat(sprintf("Layout: %d×5 grid (multi-variable) with legend on bottom\n", n_variables))
    }
    
    # Build rows - each row contains all scenarios for one variable
    row_plots <- list()
    
    for (var_idx in seq_along(variables)) {
      variable <- variables[var_idx]
      
      # Get all scenario plots for this variable
      scenario_plots <- lapply(scenarios, function(scen) {
        plot_grid[[variable]][[scen]]
      })
      
      # Combine horizontally using explicit pipe operator
      # For 5 scenarios: p1 | p2 | p3 | p4 | p5
      if (length(scenario_plots) == 1) {
        row_combined <- scenario_plots[[1]]
      } else if (length(scenario_plots) == 2) {
        row_combined <- scenario_plots[[1]] | scenario_plots[[2]]
      } else if (length(scenario_plots) == 3) {
        row_combined <- scenario_plots[[1]] | scenario_plots[[2]] | scenario_plots[[3]]
      } else if (length(scenario_plots) == 4) {
        row_combined <- scenario_plots[[1]] | scenario_plots[[2]] | scenario_plots[[3]] | scenario_plots[[4]]
      } else if (length(scenario_plots) == 5) {
        row_combined <- scenario_plots[[1]] | scenario_plots[[2]] | scenario_plots[[3]] | scenario_plots[[4]] | scenario_plots[[5]]
      }
      
      row_plots[[var_idx]] <- row_combined
    }
    
    # Combine all rows vertically using /
    if (length(row_plots) == 1) {
      main_grid <- row_plots[[1]]
    } else if (length(row_plots) == 2) {
      main_grid <- row_plots[[1]] / row_plots[[2]]
    } else if (length(row_plots) == 3) {
      main_grid <- row_plots[[1]] / row_plots[[2]] / row_plots[[3]]
    } else if (length(row_plots) == 4) {
      main_grid <- row_plots[[1]] / row_plots[[2]] / row_plots[[3]] / row_plots[[4]]
    } else if (length(row_plots) == 5) {
      main_grid <- row_plots[[1]] / row_plots[[2]] / row_plots[[3]] / row_plots[[4]] / row_plots[[5]]
    }
    
    # Apply tight spacing to main grid with extra bottom margin for legend clearance
    main_grid <- main_grid + plot_layout() & 
      theme(plot.margin = unit(c(1, 2, 3, 2), "mm"))  # Extra bottom margin to prevent legend overlap
    
    # Add legend to bottom
    combined_plot <- main_grid / legend_grob
    
    # Set relative heights (main grid gets most space, legend gets small but readable strip)
    height_ratios <- c(rep(1, n_variables), 0.12)
    combined_plot <- combined_plot + patchwork::plot_layout(
      heights = height_ratios,
      design = NULL
    ) +
      theme(plot.margin = margin(0, 0, 0, 0))  # Remove plot margins to reduce spacing
  }
  
  if (verbose) {
    cat("Dashboard assembly complete\n")
  }
  
  return(combined_plot)
}

#' @title Add Title and Annotations to Dashboard
#' @description
#' Adds overall title, subtitle with summary statistics, and formatting to
#' the assembled dashboard. Uses patchwork::plot_annotation() to add top-level
#' text elements with appropriate styling for publication.
#' 
#' Summary statistics include number of scenarios, total combinations analyzed,
#' and number of feasible combinations. Optional custom subtitle can override
#' the automatic summary.
#' 
#' @param dashboard_object Combined patchwork object from assemble_delayed_deployment_dashboard()
#' @param data Full dataset used for creating plots
#' @param title Character string for main title
#' @param subtitle Optional character string for subtitle. If NULL, generates
#'   automatic summary with statistics. Default is NULL.
#' @param add_statistics Logical indicating whether to include statistics in
#'   subtitle. Only used if subtitle is NULL. Default is TRUE.
#' @param title_size Numeric size for title text (default: 12)
#' @param subtitle_size Numeric size for subtitle text (default: 10)
#' @param verbose Logical indicating whether to print annotation information
#' 
#' @return Annotated patchwork object with title and subtitle
#' 
#' @details
#' Default subtitle format (when subtitle = NULL and add_statistics = TRUE):
#' "N scenarios, X combinations, Y feasible"
#' 
#' Additional information can be added to the subtitle based on analysis type,
#' such as arrow filtering thresholds or other relevant parameters.
#' 
#' Title and subtitle are centered (hjust = 0.5) and use appropriate text
#' sizing for readability in both print and digital formats.
#' 
#' @examples
#' # Add title with automatic statistics
#' final_plot <- add_delayed_deployment_annotations(
#'   dashboard_object = dashboard,
#'   data = plot_data,
#'   title = "Delayed Deployment Analysis: Peak Temperature"
#' )
#' 
#' # Add title with custom subtitle
#' final_plot <- add_delayed_deployment_annotations(
#'   dashboard_object = dashboard,
#'   data = plot_data,
#'   title = "Multi-Variable Cost Analysis",
#'   subtitle = "Shared scale enables direct comparison across cost components"
#' )
#' 
#' # Add title without statistics
#' final_plot <- add_delayed_deployment_annotations(
#'   dashboard_object = dashboard,
#'   data = plot_data,
#'   title = "Temperature and Overshoot Analysis",
#'   add_statistics = FALSE
#' )
add_delayed_deployment_annotations <- function(dashboard_object,
                                               data,
                                               title,
                                               subtitle = NULL,
                                               add_statistics = TRUE,
                                               title_size = 12,
                                               subtitle_size = 10,
                                               verbose = TRUE) {
  
  # Validate inputs
  if (!inherits(dashboard_object, "patchwork")) {
    stop("dashboard_object must be a patchwork object")
  }
  
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("data must be a non-empty data frame")
  }
  
  if (verbose) {
    cat("Skipping annotations (titles and subtitles removed per user request)\n")
  }
  
  # Return dashboard without annotations
  return(dashboard_object)
}

#' @title Save Dashboard to File
#' @description
#' Saves the completed dashboard to a PDF file in the figs/ directory using
#' high-quality cairo_pdf device. Handles filename generation if not provided
#' and creates the output directory if it doesn't exist.
#' 
#' Uses here::here() for robust path handling and ggsave() with cairo_pdf for
#' publication-quality vector output with proper font embedding.
#' 
#' @param plot_object Annotated patchwork object to save
#' @param filename Character string for output filename. If NULL, generates
#'   automatic filename with timestamp. Should include .pdf extension.
#'   Default is NULL.
#' @param width Numeric width in mm (default: 297, A4 landscape width)
#' @param height Numeric height in mm (default: 210, A4 landscape height)
#' @param output_dir Character string for output directory relative to project root.
#'   Default is "figs"
#' @param verbose Logical indicating whether to print save information
#' 
#' @return Invisibly returns the filepath where the plot was saved
#' 
#' @details
#' Default dimensions (297 × 210 mm) correspond to A4 landscape orientation,
#' which is suitable for most dashboard layouts. For multi-row dashboards with
#' 3+ variables, consider increasing height proportionally.
#' 
#' The cairo_pdf device ensures:
#' - High-quality vector output
#' - Proper font embedding
#' - Accurate color reproduction
#' - Professional publication quality
#' 
#' Automatic filename format: "delayed_deployment_dashboard_YYYYMMDD_HHMMSS.pdf"
#' 
#' @examples
#' # Save with automatic filename
#' save_delayed_deployment_dashboard(final_plot)
#' 
#' # Save with custom filename
#' save_delayed_deployment_dashboard(
#'   final_plot,
#'   filename = "temperature_delay_analysis.pdf"
#' )
#' 
#' # Save with custom dimensions for multi-row layout
#' save_delayed_deployment_dashboard(
#'   final_plot,
#'   filename = "cost_comparison_3variables.pdf",
#'   width = 297,
#'   height = 260  # Taller for 3 rows
#' )
save_delayed_deployment_dashboard <- function(plot_object,
                                              filename = NULL,
                                              width = 297,
                                              height = 210,
                                              output_dir = "figs",
                                              verbose = TRUE) {
  
  # Validate inputs
  if (!inherits(plot_object, "patchwork") && !inherits(plot_object, "gg")) {
    stop("plot_object must be a patchwork or ggplot object")
  }
  
  # Check here package is available
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("here package is required for file path handling. Please install it.")
  }
  
  # Generate filename if not provided
  if (is.null(filename)) {
    filename <- paste0("delayed_deployment_dashboard_", 
                       format(Sys.time(), "%Y%m%d_%H%M%S"), 
                       ".pdf")
    
    if (verbose) {
      cat(sprintf("Generated filename: %s\n", filename))
    }
  }
  
  # Ensure filename has .pdf extension
  if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".pdf")
    
    if (verbose) {
      cat(sprintf("Added .pdf extension: %s\n", filename))
    }
  }
  
  # Construct full filepath
  filepath <- here::here(output_dir, filename)
  
  # Create output directory if it doesn't exist
  output_dir_full <- here::here(output_dir)
  if (!dir.exists(output_dir_full)) {
    dir.create(output_dir_full, recursive = TRUE)
    if (verbose) {
      cat(sprintf("Created output directory: %s\n", output_dir_full))
    }
  }
  
  if (verbose) {
    cat(sprintf("Saving dashboard to: %s\n", filepath))
    cat(sprintf("Dimensions: %d × %d mm\n", width, height))
  }
  
  # Save plot using ggsave with cairo_pdf
  ggsave(
    filename = filepath,
    plot = plot_object,
    width = width,
    height = height,
    units = "mm",
    device = cairo_pdf,
    dpi = 300
  )
  
  if (verbose) {
    cat(sprintf("Dashboard saved successfully to: %s\n", filepath))
  }
  
  # Return filepath invisibly
  invisible(filepath)
}

# ==============================================================================
# Section 5: Main Function
# ==============================================================================
#
# This section contains the primary user-facing function that orchestrates the
# entire dashboard creation process. This is the main entry point for users and
# handles all aspects of dashboard generation from data preparation through
# final output.
#
# The main function validates inputs, calls all helper functions in the correct
# sequence, manages parameters, and ensures consistent high-quality output. It
# provides sensible defaults while allowing full customization of all aspects.
#
# Functions in this section:
# - create_delayed_deployment_dashboard(): Main user-facing function
# ==============================================================================

#' @title Create Delayed Deployment Dashboard
#' @description
#' Main user-facing function to create comprehensive heatmap dashboards showing
#' how delays in mitigation and CDR deployment affect climate outcomes and costs
#' across multiple SSP scenarios. This function orchestrates the entire workflow
#' from data preparation through final output.
#' 
#' Supports visualization of any combination of variables with optional contours,
#' gradient arrows, and infeasible markers. Handles both single-variable and
#' multi-variable layouts automatically.
#' 
#' @param deployment_results Results object from delayed deployment analysis
#'   Must contain a 'combined_results' data frame with required columns
#' @param variables Character vector of variable names to plot. Examples:
#'   - "peak_temperature": Peak temperature anomaly (°C)
#'   - "years_above_1p5": Years above 1.5°C threshold
#'   - "abatement_cost": Total abatement costs (mitigation + CDR)
#'   - "temp_cost": Temperature-related damage costs
#'   - "total_cost": Sum of abatement and damage costs
#' @param color_palettes Optional named list of palette specifications for each variable.
#'   Each element should be a list with 'option' and 'direction'. If NULL, uses defaults.
#'   Default is NULL.
#' @param custom_labels Optional named list of custom legend labels for each variable.
#'   If NULL, uses defaults from DELAYED_DEPLOYMENT_LABELS. Default is NULL.
#' @param shared_scale Logical indicating whether to use shared color scale across
#'   all variables (TRUE) or independent scales (FALSE). Shared scales enable
#'   direct magnitude comparison. Default is FALSE.
#' @param use_scale_limits Logical indicating whether to cap the colour scale at a
#'   percentile threshold (TRUE) or use the full data range (FALSE). When TRUE,
#'   values above the threshold are displayed as the maximum colour. Useful for
#'   data with extreme outliers. Default is FALSE.
#' @param scale_limit_percentile Numeric percentile (0-100) at which to cap the
#'   colour scale when use_scale_limits = TRUE. For example, 95 means the colour
#'   scale spans from minimum to 95th percentile, with all higher values shown
#'   as the maximum colour. Default is 95.
#' @param add_contours Logical indicating whether to add contour lines (default: TRUE)
#' @param contour_breaks Either "auto" for automatic calculation or a named list
#'   of numeric vectors specifying breaks for each variable. Default is "auto".
#' @param contour_alpha Transparency of contour lines, 0-1 (default: 0.6)
#' @param add_arrows Logical indicating whether to add gradient vector field arrows.
#'   Requires gradient calculation. Default is FALSE.
#' @param arrow_scale Numeric scaling factor for arrow length (default: 3.0)
#' @param arrow_skip Integer for arrow thinning. 1 = all arrows, 2 = every other, etc.
#'   (default: 1)
#' @param arrow_size Width of arrow lines in mm (default: 0.5)
#' @param arrow_alpha Transparency of arrows, 0-1 (default: 0.7)
#' @param min_magnitude Minimum gradient magnitude threshold for arrow display.
#'   Arrows with magnitude below this are hidden. Default is 0 (show all).
#' @param show_infeasible Logical indicating whether to mark infeasible parameter
#'   combinations with red X markers (default: TRUE)
#' @param title Character string for main dashboard title. If NULL, generates
#'   automatic title based on variables. Default is NULL.
#' @param subtitle Optional character string for subtitle. If NULL, generates
#'   automatic summary with statistics. Default is NULL.
#' @param save_plot Logical indicating whether to save plot to file (default: FALSE)
#' @param filename Character string for output filename. If NULL and save_plot = TRUE,
#'   generates automatic filename. Default is NULL.
#' @param width Numeric width in mm for saved plot (default: 297, A4 landscape)
#' @param height Numeric height in mm for saved plot (default: 210, A4 landscape)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' 
#' @return Complete patchwork dashboard object. If save_plot = TRUE, also saves
#'   to file and invisibly returns the filepath.
#' 
#' @details
#' **Workflow:**
#' 1. Validate inputs and check required packages
#' 2. Prepare and filter data
#' 3. Calculate gradients if arrows requested
#' 4. Determine color scale limits (shared or independent)
#' 5. Get palette and label information
#' 6. Create all individual plots
#' 7. Extract and position legend
#' 8. Assemble dashboard with patchwork
#' 9. Add title and annotations
#' 10. Save to file if requested
#' 
#' **Layout:**
#' - Single variable: 3×2 grid (positions A-E + blank F) with legend on right
#' - Multiple variables: N×5 grid (N rows, 5 scenarios per row) with legend on bottom
#' 
#' **Variable calculation:**
#' If arrows requested, the function automatically calculates gradients for the
#' first variable in the list. For multiple variables, only the first gets arrows.
#' 
#' @examples
#' \dontrun{
#' # Basic usage: temperature analysis with contours
#' dashboard <- create_delayed_deployment_dashboard(
#'   deployment_results = results,
#'   variables = c("peak_temperature", "years_above_1p5"),
#'   add_contours = TRUE
#' )
#' print(dashboard)
#' 
#' # Cost analysis with gradient arrows
#' dashboard <- create_delayed_deployment_dashboard(
#'   deployment_results = results,
#'   variables = c("abatement_cost"),
#'   add_arrows = TRUE,
#'   arrow_scale = 5.0,
#'   min_magnitude = 10,
#'   save_plot = TRUE,
#'   filename = "abatement_cost_arrows.pdf"
#' )
#' 
#' # Multi-variable comparison with shared scale
#' dashboard <- create_delayed_deployment_dashboard(
#'   deployment_results = results,
#'   variables = c("total_cost", "abatement_cost", "temp_cost"),
#'   shared_scale = TRUE,
#'   title = "Comprehensive Cost Comparison",
#'   subtitle = "Shared scale enables direct comparison",
#'   save_plot = TRUE,
#'   height = 260  # Taller for 3 rows
#' )
#' 
#' # Custom color palette
#' custom_palettes <- list(
#'   peak_temperature = list(option = "magma", direction = -1)
#' )
#' dashboard <- create_delayed_deployment_dashboard(
#'   deployment_results = results,
#'   variables = c("peak_temperature"),
#'   color_palettes = custom_palettes
#' )
#' 
#' # Use percentile capping for data with extreme outliers
#' dashboard <- create_delayed_deployment_dashboard(
#'   deployment_results = results,
#'   variables = c("total_cost"),
#'   use_scale_limits = TRUE,
#'   scale_limit_percentile = 95,  # Cap at 95th percentile
#'   add_contours = TRUE
#' )
#' }
create_delayed_deployment_dashboard <- function(deployment_results,
                                                variables,
                                                color_palettes = NULL,
                                                custom_labels = NULL,
                                                shared_scale = FALSE,
                                                use_scale_limits = FALSE,
                                                scale_limit_percentile = 95,
                                                add_contours = TRUE,
                                                contour_breaks = "auto",
                                                contour_alpha = 0.6,
                                                add_arrows = FALSE,
                                                arrow_scale = 3.0,
                                                arrow_skip = 1,
                                                arrow_size = 0.5,
                                                arrow_alpha = 0.7,
                                                min_magnitude = 0,
                                                show_infeasible = TRUE,
                                                title = NULL,
                                                subtitle = NULL,
                                                save_plot = FALSE,
                                                filename = NULL,
                                                width = 297,
                                                height = 210,
                                                verbose = TRUE) {
  
  # ============================================================================
  # Step 1: Validate inputs and check required packages
  # ============================================================================
  
  if (verbose) {
    cat("\n=== DELAYED DEPLOYMENT DASHBOARD CREATION ===\n\n")
    cat("Step 1: Validating inputs and checking packages\n")
  }
  
  # Check required packages
  required_packages <- c("ggplot2", "dplyr", "patchwork", "viridis", "cowplot", "here")
  if (add_arrows) {
    required_packages <- c(required_packages, "metR")
  }
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages not installed: ", paste(missing_packages, collapse = ", "),
         "\nPlease install them with: install.packages(c('", 
         paste(missing_packages, collapse = "', '"), "'))")
  }
  
  # Validate variables
  if (length(variables) == 0) {
    stop("variables must contain at least one variable name")
  }
  
  if (verbose) {
    cat(sprintf("  Requested variables: %s\n", paste(variables, collapse = ", ")))
    cat(sprintf("  Features: contours=%s, arrows=%s, infeasible=%s\n",
                add_contours, add_arrows, show_infeasible))
  }
  
  # ============================================================================
  # Step 2: Prepare and filter data
  # ============================================================================
  
  if (verbose) {
    cat("\nStep 2: Preparing data\n")
  }
  
  plot_data <- prepare_delayed_deployment_data(
    deployment_results = deployment_results,
    variables = variables,
    verbose = verbose
  )
  
  # ============================================================================
  # Step 3: Calculate gradients if arrows requested
  # ============================================================================
  
  if (add_arrows) {
    if (verbose) {
      cat("\nStep 3: Calculating gradients for arrows\n")
    }
    
    # Calculate gradients for the first variable (primary variable for arrows)
    arrow_variable <- variables[1]
    
    plot_data <- calculate_delayed_deployment_gradients(
      data = plot_data,
      variable = arrow_variable,
      verbose = verbose
    )
    
    # Calculate magnitude limits for consistent arrow coloring
    mag_limits <- range(plot_data$gradient_mag, na.rm = TRUE)
    
    if (verbose) {
      cat(sprintf("  Gradient magnitude range: %.4f to %.4f\n", 
                  mag_limits[1], mag_limits[2]))
    }
  } else {
    mag_limits <- NULL
    
    if (verbose) {
      cat("\nStep 3: Skipping gradient calculation (arrows not requested)\n")
    }
  }
  
  # ============================================================================
  # Step 4: Determine color scale limits
  # ============================================================================
  
  if (verbose) {
    cat("\nStep 4: Calculating color scale limits\n")
  }
  
  variable_limits <- calculate_variable_limits(
    data = plot_data,
    variables = variables,
    shared_scale = shared_scale,
    use_scale_limits = use_scale_limits,
    scale_limit_percentile = scale_limit_percentile,
    verbose = verbose
  )
  
  # ============================================================================
  # Step 5: Get palette and label information
  # ============================================================================
  
  if (verbose) {
    cat("\nStep 5: Configuring palettes and labels\n")
  }
  
  # Get palette info for each variable
  palette_info <- list()
  for (var in variables) {
    palette_info[[var]] <- get_variable_palette(var, color_palettes)
  }
  
  # Get label info for each variable
  # Use single-line labels for multi-variable layouts (bottom legend)
  use_single_line <- (length(variables) > 1)
  variable_labels <- list()
  for (var in variables) {
    variable_labels[[var]] <- get_variable_label(var, custom_labels, single_line = use_single_line)
  }
  
  if (verbose) {
    cat(sprintf("  Configured %d variable(s)\n", length(variables)))
  }
  
  # ============================================================================
  # Step 6: Create all individual plots
  # ============================================================================
  
  if (verbose) {
    cat("\nStep 6: Creating plot grid\n")
  }
  
  plot_grid <- create_delayed_deployment_plot_grid(
    data = plot_data,
    variables = variables,
    variable_limits = variable_limits,
    palette_info = palette_info,
    variable_labels = variable_labels,
    add_contours = add_contours,
    contour_breaks = contour_breaks,
    contour_alpha = contour_alpha,
    add_arrows = add_arrows,
    arrow_scale = arrow_scale,
    arrow_skip = arrow_skip,
    min_magnitude = min_magnitude,
    arrow_size = arrow_size,
    arrow_alpha = arrow_alpha,
    mag_limits = mag_limits,
    show_infeasible = show_infeasible,
    verbose = verbose
  )
  
  # ============================================================================
  # Step 7: Extract and position legend
  # ============================================================================
  
  if (verbose) {
    cat("\nStep 7: Extracting legend\n")
  }
  
  # Create a plot with legend visible for extraction
  first_variable <- variables[1]
  first_scenario <- SSP_SCENARIO_ORDER[1]
  
  # Get first plot data for legend extraction
  legend_data <- plot_data %>% filter(scenario_short == first_scenario)
  
  # Create temporary plot with legend
  temp_plot <- create_delayed_deployment_base_heatmap(
    scenario_data = legend_data,
    variable = first_variable,
    variable_limits = if ("shared" %in% names(variable_limits)) {
      variable_limits$shared
    } else {
      variable_limits[[first_variable]]
    },
    palette_info = palette_info[[first_variable]],
    variable_label = variable_labels[[first_variable]],
    scenario_name = first_scenario,
    show_title = FALSE,
    theme_object = get_delayed_deployment_theme(multi_row = length(variables) > 1)
  )
  
  # Add legend to theme with appropriate sizing
  if (length(variables) == 1) {
    # Right legend for single variable
    temp_plot <- temp_plot + theme(legend.position = "right")
  } else {
    # Bottom legend for multiple variables - make it smaller
    temp_plot <- temp_plot + theme(
      legend.position = "bottom",
      legend.key.width = unit(2, "cm"),   # Width of the color bar
      legend.key.height = unit(0.3, "cm"), # Height of the color bar (make this smaller)
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5.5)
    )
  }
  
  # Determine legend position for final layout
  legend_position <- if (length(variables) == 1) "right" else "bottom"
  
  # Extract legend
  legend_grob <- extract_delayed_deployment_legend(
    plot_object = temp_plot,
    legend_position = legend_position
  )
  
  if (verbose) {
    cat(sprintf("  Legend position: %s\n", legend_position))
  }
  
  # ============================================================================
  # Step 8: Assemble dashboard with patchwork
  # ============================================================================
  
  if (verbose) {
    cat("\nStep 8: Assembling dashboard\n")
  }
  
  dashboard <- assemble_delayed_deployment_dashboard(
    plot_grid = plot_grid,
    legend_grob = legend_grob,
    variables = variables,
    scenarios = SSP_SCENARIO_ORDER,
    legend_position = legend_position,
    verbose = verbose
  )
  
  # ============================================================================
  # Step 9: Add title and annotations
  # ============================================================================
  
  if (verbose) {
    cat("\nStep 9: Adding title and annotations\n")
  }
  
  # Generate title if not provided
  if (is.null(title)) {
    if (length(variables) == 1) {
      # Create human-readable title from variable name
      title_var <- gsub("_", " ", variables[1])
      title_var <- tools::toTitleCase(title_var)
      title <- paste("Delayed Deployment Analysis:", title_var)
    } else {
      title <- "Delayed Deployment Analysis: Multi-Variable Comparison"
    }
  }
  
  final_plot <- add_delayed_deployment_annotations(
    dashboard_object = dashboard,
    data = plot_data,
    title = title,
    subtitle = subtitle,
    add_statistics = is.null(subtitle),
    verbose = verbose
  )
  
  # ============================================================================
  # Step 10: Save to file if requested
  # ============================================================================
  
  if (save_plot) {
    if (verbose) {
      cat("\nStep 10: Saving dashboard to file\n")
    }
    
    filepath <- save_delayed_deployment_dashboard(
      plot_object = final_plot,
      filename = filename,
      width = width,
      height = height,
      verbose = verbose
    )
    
    if (verbose) {
      cat(sprintf("\n=== DASHBOARD CREATION COMPLETE ===\n"))
      cat(sprintf("Saved to: %s\n\n", filepath))
    }
    
    return(invisible(final_plot))
  } else {
    if (verbose) {
      cat("\nStep 10: Skipping file save (save_plot = FALSE)\n")
      cat("\n=== DASHBOARD CREATION COMPLETE ===\n\n")
    }
    
    return(final_plot)
  }
}

# ==============================================================================
# Section 6: Convenience Wrappers
# ==============================================================================
#
# These functions provide convenient shortcuts for common dashboard types,
# with pre-configured settings for specific analysis scenarios. They wrap the
# main create_delayed_deployment_dashboard() function with appropriate defaults
# for common use cases.
#
# These wrappers reduce boilerplate code and ensure consistent styling for
# standard analyses while still allowing parameter overrides when needed.
#
# Functions in this section:
# - create_temperature_delay_dashboard(): Temperature and overshoot analysis
# - create_cost_delay_dashboard(): Cost analysis with gradient arrows
# - create_comprehensive_cost_delay_dashboard(): Multi-variable cost comparison
# ==============================================================================

#' @title Create Temperature and Overshoot Dashboard
#' @description
#' Convenience wrapper for creating a dashboard showing peak temperature and
#' years above 1.5°C threshold. Pre-configured with appropriate settings for
#' temperature analysis including plasma color palette and contour lines.
#' 
#' This is a common analysis showing both the peak warming achieved and the
#' duration of temperature overshoot above the Paris Agreement target.
#' 
#' @param deployment_results Results object from delayed deployment analysis
#' @param save_plot Logical indicating whether to save plot to file (default: FALSE)
#' @param filename Character string for output filename (default: NULL for auto-generation)
#' @param ... Additional arguments passed to create_delayed_deployment_dashboard()
#' 
#' @return Dashboard with peak temperature (top row) and years above 1.5°C (bottom row)
#' 
#' @details
#' Pre-configured settings:
#' - Variables: "peak_temperature", "years_above_1p5"
#' - Color palette: plasma (reversed) for both variables
#' - Contours: enabled with automatic breaks
#' - Arrows: disabled (not typically used for temperature analysis)
#' - Layout: 3×2 grid with legend on right
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' dashboard <- create_temperature_delay_dashboard(multi_results)
#' print(dashboard)
#' 
#' # Save to file
#' create_temperature_delay_dashboard(
#'   multi_results,
#'   save_plot = TRUE,
#'   filename = "temperature_overshoot_analysis.pdf"
#' )
#' 
#' # Customize with additional arguments
#' create_temperature_delay_dashboard(
#'   multi_results,
#'   contour_alpha = 0.8,
#'   show_infeasible = FALSE
#' )
#' }
create_temperature_delay_dashboard <- function(deployment_results,
                                               save_plot = FALSE,
                                               filename = NULL,
                                               ...) {
  
  # Pre-configure for temperature analysis
  create_delayed_deployment_dashboard(
    deployment_results = deployment_results,
    variables = c("peak_temperature", "years_above_1p5"),
    shared_scale = FALSE,  # Independent scales (different units)
    add_contours = TRUE,
    add_arrows = FALSE,    # Temperature analysis typically doesn't use arrows
    title = "Delayed Deployment Analysis: Temperature and Overshoot",
    save_plot = save_plot,
    filename = filename,
    ...
  )
}

#' @title Create Cost Analysis Dashboard with Gradient Arrows
#' @description
#' Convenience wrapper for creating a dashboard showing cost metrics with
#' gradient vector field arrows. Pre-configured for cost analysis with
#' viridis color palette, contours, and arrows showing direction of steepest
#' cost increase.
#' 
#' Useful for analyzing abatement costs or damage costs individually with
#' visual indication of cost gradients in the delay space.
#' 
#' @param deployment_results Results object from delayed deployment analysis
#' @param cost_variable Character string naming the cost variable to plot.
#'   Options: "abatement_cost", "temp_cost", "total_cost", "mitig_cost", "remov_cost"
#'   Default is "abatement_cost"
#' @param arrow_scale Numeric scaling factor for arrow length (default: 3.0)
#' @param min_magnitude Minimum gradient magnitude for arrow display (default: 0)
#' @param save_plot Logical indicating whether to save plot to file (default: FALSE)
#' @param filename Character string for output filename (default: NULL for auto-generation)
#' @param ... Additional arguments passed to create_delayed_deployment_dashboard()
#' 
#' @return Dashboard with cost heatmap, contours, and gradient arrows
#' 
#' @details
#' Pre-configured settings:
#' - Variables: User-specified cost variable
#' - Color palette: viridis (reversed)
#' - Contours: enabled with automatic breaks
#' - Arrows: enabled, colored by gradient magnitude (plasma palette)
#' - Layout: 3×2 grid with legend on right
#' 
#' Arrows show the direction of steepest cost increase, with color indicating
#' the rate of change (gradient magnitude).
#' 
#' @examples
#' \dontrun{
#' # Abatement cost analysis (default)
#' dashboard <- create_cost_delay_dashboard(multi_results)
#' 
#' # Temperature damage cost analysis
#' dashboard <- create_cost_delay_dashboard(
#'   multi_results,
#'   cost_variable = "temp_cost"
#' )
#' 
#' # Customize arrow appearance
#' dashboard <- create_cost_delay_dashboard(
#'   multi_results,
#'   arrow_scale = 5.0,
#'   min_magnitude = 10,
#'   save_plot = TRUE,
#'   filename = "abatement_cost_arrows.pdf"
#' )
#' }
create_cost_delay_dashboard <- function(deployment_results,
                                        cost_variable = "abatement_cost",
                                        arrow_scale = 3.0,
                                        min_magnitude = 0,
                                        save_plot = FALSE,
                                        filename = NULL,
                                        ...) {
  
  # Validate cost variable
  valid_cost_vars <- c("abatement_cost", "temp_cost", "total_cost", 
                       "mitig_cost", "remov_cost")
  
  if (!cost_variable %in% valid_cost_vars) {
    stop("cost_variable must be one of: ", paste(valid_cost_vars, collapse = ", "))
  }
  
  # Generate title from variable name
  title_var <- gsub("_", " ", cost_variable)
  title_var <- tools::toTitleCase(title_var)
  default_title <- paste("Delayed Deployment Analysis:", title_var, "with Gradient Arrows")
  
  # Pre-configure for cost analysis with arrows
  create_delayed_deployment_dashboard(
    deployment_results = deployment_results,
    variables = cost_variable,
    shared_scale = FALSE,
    add_contours = TRUE,
    add_arrows = TRUE,
    arrow_scale = arrow_scale,
    min_magnitude = min_magnitude,
    title = default_title,
    save_plot = save_plot,
    filename = filename,
    ...
  )
}

#' @title Create Comprehensive Cost Comparison Dashboard
#' @description
#' Convenience wrapper for creating a multi-row dashboard comparing different
#' cost components (total costs, abatement costs, and temperature damage costs).
#' Pre-configured with shared color scale to enable direct visual comparison
#' of cost magnitudes across all three components.
#' 
#' This dashboard format is useful for understanding the relative contributions
#' of different cost components and how they vary across delay scenarios.
#' 
#' @param deployment_results Results object from delayed deployment analysis
#' @param shared_scale Logical indicating whether to use shared color scale.
#'   TRUE (default) enables direct magnitude comparison. FALSE uses independent scales.
#' @param save_plot Logical indicating whether to save plot to file (default: FALSE)
#' @param filename Character string for output filename (default: NULL for auto-generation)
#' @param height Numeric height in mm for saved plot. Default is 260 (taller than
#'   standard A4 landscape to accommodate 3 rows)
#' @param ... Additional arguments passed to create_delayed_deployment_dashboard()
#' 
#' @return Dashboard with three rows: total costs (A-E), abatement costs (F-J),
#'   and temperature damage costs (K-O), with legend on bottom
#' 
#' @details
#' Pre-configured settings:
#' - Variables: "total_cost", "abatement_cost", "temp_cost"
#' - Color palette: viridis (reversed) for all
#' - Shared scale: TRUE by default (can be overridden)
#' - Contours: enabled with automatic breaks
#' - Arrows: disabled (can be enabled via ...)
#' - Layout: 3×5 grid (3 rows × 5 scenarios) with legend on bottom
#' - Height: 260 mm (taller to accommodate 3 rows)
#' 
#' When shared_scale = TRUE, all three cost types use the same color scale,
#' making it easy to visually compare their relative magnitudes. This is
#' particularly useful for identifying which cost component dominates in
#' different regions of the delay space.
#' 
#' @examples
#' \dontrun{
#' # Basic usage with shared scale
#' dashboard <- create_comprehensive_cost_delay_dashboard(multi_results)
#' print(dashboard)
#' 
#' # Save to file
#' create_comprehensive_cost_delay_dashboard(
#'   multi_results,
#'   save_plot = TRUE,
#'   filename = "comprehensive_cost_comparison.pdf"
#' )
#' 
#' # Use independent scales for each cost type
#' dashboard <- create_comprehensive_cost_delay_dashboard(
#'   multi_results,
#'   shared_scale = FALSE
#' )
#' 
#' # Add arrows to show cost gradients
#' dashboard <- create_comprehensive_cost_delay_dashboard(
#'   multi_results,
#'   add_arrows = TRUE,
#'   arrow_scale = 5.0
#' )
#' }
create_comprehensive_cost_delay_dashboard <- function(deployment_results,
                                                      shared_scale = TRUE,
                                                      save_plot = FALSE,
                                                      filename = NULL,
                                                      height = 260,
                                                      ...) {
  
  # Pre-configure for comprehensive cost comparison
  create_delayed_deployment_dashboard(
    deployment_results = deployment_results,
    variables = c("total_cost", "abatement_cost", "temp_cost"),
    shared_scale = shared_scale,
    add_contours = TRUE,
    add_arrows = FALSE,  # Can be overridden via ...
    title = "Comprehensive Cost Comparison: Total, Abatement, and Temperature Damage",
    subtitle = if (shared_scale) {
      "Shared scale enables direct comparison across cost components"
    } else {
      NULL  # Use default subtitle with statistics
    },
    save_plot = save_plot,
    filename = filename,
    height = height,  # Taller for 3 rows
    ...
  )
}


# ==============================================================================
# Section 7: Usage Examples
# ==============================================================================
#
# Below are comprehensive examples demonstrating how to use the delayed
# deployment visualization functions. Examples progress from basic to advanced
# usage, showing the flexibility and power of the refactored system.
# ==============================================================================

# -----------------------------------------------------------------------------
# Example 1: Basic temperature analysis
# -----------------------------------------------------------------------------
# dashboard <- create_temperature_delay_dashboard(
#   deployment_results = multi_results
# )
# print(dashboard)

# -----------------------------------------------------------------------------
# Example 2: Cost analysis with gradient arrows
# -----------------------------------------------------------------------------
# dashboard <- create_cost_delay_dashboard(
#   deployment_results = multi_results,
#   cost_variable = "abatement_cost",
#   arrow_scale = 5.0,
#   min_magnitude = 10
# )
# print(dashboard)

# -----------------------------------------------------------------------------
# Example 3: Multi-variable cost comparison with shared scale
# -----------------------------------------------------------------------------
# dashboard <- create_comprehensive_cost_delay_dashboard(
#   deployment_results = multi_results,
#   shared_scale = TRUE,
#   save_plot = TRUE,
#   filename = "cost_comparison.pdf"
# )

# -----------------------------------------------------------------------------
# Example 4: Custom variable selection
# -----------------------------------------------------------------------------
# dashboard <- create_delayed_deployment_dashboard(
#   deployment_results = multi_results,
#   variables = c("peak_temperature"),
#   add_contours = TRUE,
#   add_arrows = FALSE,
#   save_plot = TRUE,
#   filename = "peak_temperature_analysis.pdf"
# )

# -----------------------------------------------------------------------------
# Example 5: Custom color palettes
# -----------------------------------------------------------------------------
# custom_palettes <- list(
#   peak_temperature = list(option = "magma", direction = -1),
#   years_above_1p5 = list(option = "inferno", direction = -1)
# )
# 
# dashboard <- create_delayed_deployment_dashboard(
#   deployment_results = multi_results,
#   variables = c("peak_temperature", "years_above_1p5"),
#   color_palettes = custom_palettes
# )

# -----------------------------------------------------------------------------
# Example 6: Temperature damage costs with arrows
# -----------------------------------------------------------------------------
# dashboard <- create_delayed_deployment_dashboard(
#   deployment_results = multi_results,
#   variables = c("temp_cost"),
#   add_contours = TRUE,
#   add_arrows = TRUE,
#   arrow_scale = 3.0,
#   arrow_skip = 1,
#   min_magnitude = 5,
#   save_plot = TRUE,
#   filename = "temp_damage_arrows.pdf",
#   width = 297,
#   height = 210
# )

# -----------------------------------------------------------------------------
# Example 7: Multi-variable with independent scales
# -----------------------------------------------------------------------------
# dashboard <- create_delayed_deployment_dashboard(
#   deployment_results = multi_results,
#   variables = c("peak_temperature", "abatement_cost", "temp_cost"),
#   shared_scale = FALSE,  # Each variable gets its own scale
#   add_contours = TRUE,
#   height = 260,  # Taller for 3 rows
#   save_plot = TRUE
# )

# -----------------------------------------------------------------------------
# Example 8: Custom contour breaks
# -----------------------------------------------------------------------------
# custom_breaks <- list(
#   peak_temperature = c(1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2),
#   years_above_1p5 = seq(0, 100, by = 10)
# )
# 
# dashboard <- create_delayed_deployment_dashboard(
#   deployment_results = multi_results,
#   variables = c("peak_temperature", "years_above_1p5"),
#   contour_breaks = custom_breaks,
#   contour_alpha = 0.8
# )

# -----------------------------------------------------------------------------
# Example 9: Minimal styling for presentations
# -----------------------------------------------------------------------------
# dashboard <- create_delayed_deployment_dashboard(
#   deployment_results = multi_results,
#   variables = c("abatement_cost"),
#   add_contours = FALSE,  # Clean look without contours
#   show_infeasible = FALSE,  # Hide infeasible markers
#   title = "Abatement Costs: Deployment Delay Sensitivity",
#   subtitle = NULL  # No subtitle
# )

# -----------------------------------------------------------------------------
# Example 10: Complete custom workflow (advanced)
# -----------------------------------------------------------------------------
# # Step 1: Prepare data
# plot_data <- prepare_delayed_deployment_data(
#   deployment_results = multi_results,
#   variables = c("abatement_cost"),
#   verbose = TRUE
# )
# 
# # Step 2: Calculate gradients
# plot_data <- calculate_delayed_deployment_gradients(
#   data = plot_data,
#   variable = "abatement_cost",
#   verbose = TRUE
# )
# 
# # Step 3: Use main function with pre-processed data
# # (This is just to demonstrate the workflow - the main function
# # handles these steps automatically)
# dashboard <- create_delayed_deployment_dashboard(
#   deployment_results = multi_results,
#   variables = c("abatement_cost"),
#   add_arrows = TRUE,
#   arrow_scale = 5.0,
#   min_magnitude = 15,
#   verbose = TRUE
# )