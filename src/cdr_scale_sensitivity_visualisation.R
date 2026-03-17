# ==============================================================================
# CDR Scale Sensitivity Visualisation Functions
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
#   Rynne, N., Bode, M., Roberts, M., & Heneghan, R. (2025). 
#   MACROM: An Optimal Control Model for Balancing Climate Change 
#   Abatement and Damage Trade-offs.
#   https://doi.org/10.5281/zenodo.18463951
# 
# Description:
# Heatmap dashboard functions for visualising CDR scale sensitivity results.
# Plots K (CDR carrying capacity, GtCO2/year) on the x-axis and r (CDR growth
# rate) on the y-axis, with one column per SSP scenario and one row per
# requested outcome variable. Structured to mirror delayed_deployment_visualisation.R
# and accepts output from run_cdr_scale_sensitivity() in cdr_scale_sensitivity.R.
# 
# Version: 1.1.0
# Last updated: March 2026
#
# Fixes in v1.1.0:
#   - face = "bold" replaces fontface = "bold" in element_text() (Bug 1)
#   - contour break calculation now returns NULL when range < interval step,
#     preventing "wrong sign in by" seq() error (Bug 2)
#   - wrapper functions use explicit arguments instead of ... to prevent
#     duplicate argument errors (Bug 3)
#   - loop variable renamed show_panel_title to avoid shadowing the
#     show_title parameter of create_cdr_scale_base_heatmap() (Bug 4)
#   - guides = "collect" removed from wrap_plots() calls — incompatible with
#     manually extracted legend approach (Bug 5)
#   - plot_grid renamed panel_grid throughout to avoid shadowing
#     cowplot::plot_grid() (Bug 7)
# ==============================================================================


# ==============================================================================
# Global definitions
# ==============================================================================

# Default colour palettes for each outcome variable
CDR_SCALE_PALETTES <- list(
  peak_temperature = list(option = "plasma",  direction = -1),
  years_above_1p5  = list(option = "viridis",  direction = -1),
  abatement_cost   = list(option = "viridis", direction = -1),
  temp_cost        = list(option = "viridis", direction = -1),
  total_cost       = list(option = "viridis", direction = -1),
  mitig_cost       = list(option = "viridis", direction = -1),
  remov_cost       = list(option = "viridis", direction = -1)
)

# Multi-line legend labels (used when legend is positioned on the right)
CDR_SCALE_LABELS <- list(
  peak_temperature = "Peak\nTemp.\n(°C)",
  years_above_1p5  = "Years\nAbove\n1.5°C",
  abatement_cost   = "Costs\n($USD trillion)",
  temp_cost        = "Temp.\nDamage\nCost\n($ trillion)",
  total_cost       = "Total\nCost\n($ trillion)",
  mitig_cost       = "Mitigation\nCost\n($ trillion)",
  remov_cost       = "Removal\nCost\n($ trillion)",
  cost             = "Cost\n($ trillion)"
)

# Single-line legend labels (used when legend is positioned on the bottom)
CDR_SCALE_LABELS_SINGLE_LINE <- list(
  peak_temperature = "Peak Temperature (°C)",
  years_above_1p5  = "Years Above 1.5°C",
  abatement_cost   = "Costs ($USD trillion)",
  temp_cost        = "Temperature Damage Cost ($ trillion)",
  total_cost       = "Total Cost ($ trillion)",
  mitig_cost       = "Mitigation Cost ($ trillion)",
  remov_cost       = "Removal Cost ($ trillion)",
  cost             = "Cost ($ trillion)"
)

# SSP scenario order for consistent column ordering
SSP_SCENARIO_ORDER_SCALE <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")

# Base ggplot2 theme shared by all CDR scale heatmaps
BASE_CDR_SCALE_THEME <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )


# ==============================================================================
# Section 1: Data Preparation
# ==============================================================================

#' @title Prepare CDR Scale Sensitivity Data for Visualisation
#' @description
#' Extracts, validates, and filters the combined_results data frame from a
#' run_cdr_scale_sensitivity() output object. Ensures required columns are
#' present, filters to valid SSP scenarios, removes rows with missing values
#' in the requested variables, and sorts consistently.
#'
#' @param sensitivity_results Results object from run_cdr_scale_sensitivity().
#'   Can also be passed as a bare data frame (the combined_results element).
#' @param variables Character vector of variable names to validate and prepare.
#'   Default is NULL (basic structural validation only).
#' @param verbose Logical: print progress messages (default: TRUE)
#'
#' @return Data frame ready for plotting, sorted by scenario_short, K, r.
prepare_cdr_scale_data <- function(sensitivity_results,
                                   variables = NULL,
                                   verbose   = TRUE) {
  
  # Accept either a full results list or a bare data frame
  if (is.data.frame(sensitivity_results)) {
    combined_data <- sensitivity_results
  } else if (is.list(sensitivity_results) &&
             "combined_results" %in% names(sensitivity_results)) {
    combined_data <- sensitivity_results$combined_results
  } else {
    stop("sensitivity_results must be either a data frame or a list containing ",
         "'combined_results' (output from run_cdr_scale_sensitivity())")
  }
  
  if (is.null(combined_data) || nrow(combined_data) == 0) {
    stop("No valid data found in sensitivity_results")
  }
  
  # Required structural columns
  required_cols <- c("K", "r", "scenario", "scenario_short", "feasible")
  missing_cols  <- setdiff(required_cols, names(combined_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check requested outcome variables exist
  if (!is.null(variables)) {
    missing_vars <- setdiff(variables, names(combined_data))
    if (length(missing_vars) > 0) {
      stop("Requested variables not found in data: ",
           paste(missing_vars, collapse = ", "))
    }
  }
  
  # Filter to valid SSP scenarios
  plot_data <- combined_data %>%
    filter(scenario_short %in% SSP_SCENARIO_ORDER_SCALE)
  
  if (nrow(plot_data) == 0) {
    stop("No data found matching SSP scenarios: ",
         paste(SSP_SCENARIO_ORDER_SCALE, collapse = ", "))
  }
  
  # Remove rows with missing values in requested variables
  if (!is.null(variables)) {
    n_before <- nrow(plot_data)
    for (var in variables) {
      plot_data <- plot_data %>% filter(!is.na(.data[[var]]))
    }
    n_after <- nrow(plot_data)
    if (verbose && n_after < n_before) {
      cat(sprintf("Removed %d rows with missing values in specified variables\n",
                  n_before - n_after))
    }
    if (n_after == 0) {
      stop("No valid data remaining after removing rows with missing values")
    }
  }
  
  # Sort for consistent plotting
  plot_data <- plot_data %>%
    arrange(scenario_short, K, r)
  
  if (verbose) {
    n_scen  <- length(unique(plot_data$scenario_short))
    n_combs <- nrow(plot_data)
    n_feas  <- sum(plot_data$feasible, na.rm = TRUE)
    cat(sprintf("Prepared CDR scale data: %d scenarios, %d rows, %d feasible\n",
                n_scen, n_combs, n_feas))
  }
  
  return(plot_data)
}


# ==============================================================================
# Section 2: Scale Limits
# ==============================================================================

#' @title Calculate Variable Limits for CDR Scale Heatmaps
#' @description
#' Determines colour scale limits for each outcome variable. Mirrors
#' calculate_variable_limits() from delayed_deployment_visualisation.R exactly,
#' adapted for the CDR scale sensitivity context.
#'
#' @param data Prepared data frame from prepare_cdr_scale_data()
#' @param variables Character vector of variable names
#' @param shared_scale Logical: use shared limits across all variables (default: FALSE)
#' @param use_scale_limits Logical: cap colour scale at a percentile (default: FALSE)
#' @param scale_limit_percentile Numeric percentile for capping (default: 95)
#' @param verbose Logical: print limit information (default: TRUE)
#'
#' @return Named list of limits:
#'   - If shared_scale = TRUE: list(shared = c(min, max))
#'   - If shared_scale = FALSE: list(var1 = c(min, max), var2 = c(min, max), ...)
calculate_cdr_scale_limits <- function(data,
                                       variables,
                                       shared_scale           = FALSE,
                                       use_scale_limits       = FALSE,
                                       scale_limit_percentile = 95,
                                       verbose                = TRUE) {
  
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("data must be a non-empty data frame")
  }
  if (length(variables) == 0) {
    stop("variables must contain at least one variable name")
  }
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  if (shared_scale) {
    
    all_values <- unlist(data[, variables], use.names = FALSE)
    all_values <- all_values[!is.na(all_values)]
    if (length(all_values) == 0) stop("No non-missing values found")
    
    min_val <- min(all_values)
    max_val <- max(all_values)
    
    if (use_scale_limits) {
      max_val <- as.numeric(quantile(all_values, probs = scale_limit_percentile / 100))
      if (verbose) {
        cat(sprintf("Shared scale with %g%% cap: [%.4f, %.4f]\n",
                    scale_limit_percentile, min_val, max_val))
      }
    } else if (verbose) {
      cat(sprintf("Shared scale across %d variables: [%.4f, %.4f]\n",
                  length(variables), min_val, max_val))
    }
    
    return(list(shared = c(min_val, max_val)))
    
  } else {
    
    limits_list <- list()
    for (var in variables) {
      var_values <- data[[var]][!is.na(data[[var]])]
      if (length(var_values) == 0) {
        warning("No non-missing values for variable: ", var)
        next
      }
      min_val <- min(var_values)
      max_val <- max(var_values)
      if (use_scale_limits) {
        max_val <- as.numeric(quantile(var_values, probs = scale_limit_percentile / 100))
        if (verbose) {
          cat(sprintf("  %s with %g%% cap: [%.4f, %.4f]\n",
                      var, scale_limit_percentile, min_val, max_val))
        }
      } else if (verbose) {
        cat(sprintf("  %s: [%.4f, %.4f]\n", var, min_val, max_val))
      }
      limits_list[[var]] <- c(min_val, max_val)
    }
    
    if (length(limits_list) == 0) stop("Failed to calculate limits for any variable")
    return(limits_list)
  }
}


# ==============================================================================
# Section 3: Plot Component Helpers
# ==============================================================================

#' @title Get CDR Scale Heatmap Theme
#' @description
#' Returns the ggplot2 theme for CDR scale heatmap panels. Text sizes are
#' adjusted based on whether the panel sits in a multi-row layout.
#'
#' @param multi_row Logical: TRUE for 3+ row layouts (smaller text). Default FALSE.
#' @return ggplot2 theme object
get_cdr_scale_theme <- function(multi_row = FALSE) {
  
  if (multi_row) {
    base_size         <- 8
    title_size        <- 9
    axis_title_size   <- 7
    axis_text_size    <- 6
    legend_title_size <- 6
    legend_text_size  <- 5.5
    legend_key_size   <- 0.22
    margin_size       <- 2
  } else {
    base_size         <- 10
    title_size        <- 10
    axis_title_size   <- 9
    axis_text_size    <- 8
    legend_title_size <- 8
    legend_text_size  <- 7
    legend_key_size   <- 0.35
    margin_size       <- 3
  }
  
  BASE_CDR_SCALE_THEME +
    theme(
      text            = element_text(size = base_size),
      plot.title      = element_text(size = title_size, hjust = 0.5),
      axis.title      = element_text(size = axis_title_size),
      axis.text       = element_text(size = axis_text_size),
      legend.title    = element_text(size = legend_title_size),
      legend.text     = element_text(size = legend_text_size),
      legend.key.size = unit(legend_key_size, "cm"),
      legend.position = "none",
      plot.margin     = margin(1, margin_size, 1, margin_size)
    )
}


#' @title Get Variable Colour Palette for CDR Scale Heatmaps
#' @description
#' Returns palette information for a variable. Checks custom_palettes first,
#' then CDR_SCALE_PALETTES defaults, then falls back to viridis reversed.
#'
#' @param variable Character string: variable name
#' @param custom_palettes Optional named list overriding defaults. Default NULL.
#' @return List with elements 'option' (character) and 'direction' (1 or -1)
get_cdr_scale_palette <- function(variable, custom_palettes = NULL) {
  
  if (!is.null(custom_palettes) && variable %in% names(custom_palettes)) {
    palette_info <- custom_palettes[[variable]]
    if (!is.list(palette_info) ||
        !all(c("option", "direction") %in% names(palette_info))) {
      stop("Custom palette for '", variable,
           "' must be a list with 'option' and 'direction' elements")
    }
    if (!palette_info$direction %in% c(1, -1)) {
      stop("Palette direction must be 1 or -1")
    }
    return(palette_info)
  }
  
  if (variable %in% names(CDR_SCALE_PALETTES)) {
    return(CDR_SCALE_PALETTES[[variable]])
  }
  
  warning("No palette defined for variable '", variable,
          "'. Using default viridis (reversed).")
  return(list(option = "viridis", direction = -1))
}


#' @title Get Legend Label for CDR Scale Variable
#' @description
#' Returns the formatted legend label for a variable. Checks custom_labels
#' first, then CDR_SCALE_LABELS or CDR_SCALE_LABELS_SINGLE_LINE.
#'
#' @param variable Character string: variable name
#' @param custom_labels Optional named list of character labels. Default NULL.
#' @param single_line Logical: use single-line labels (for bottom legend).
#'   Default FALSE.
#' @return Character string label (may contain \n for line breaks)
get_cdr_scale_label <- function(variable,
                                custom_labels = NULL,
                                single_line   = FALSE) {
  
  if (!is.null(custom_labels) && variable %in% names(custom_labels)) {
    lbl <- custom_labels[[variable]]
    if (!is.character(lbl) || length(lbl) != 1) {
      stop("Custom label for '", variable, "' must be a single character string")
    }
    return(lbl)
  }
  
  label_set <- if (single_line) CDR_SCALE_LABELS_SINGLE_LINE else CDR_SCALE_LABELS
  
  if (variable %in% names(label_set)) {
    return(label_set[[variable]])
  }
  
  warning("No label defined for variable '", variable, "'. Using variable name.")
  basic_label <- tools::toTitleCase(gsub("_", " ", variable))
  return(basic_label)
}


#' @title Calculate Contour Breaks for CDR Scale Variable
#' @description
#' Determines intelligent contour break points using adaptive interval logic.
#' Returns NULL safely when the data range is too small to produce at least two
#' valid break values, preventing the seq() "wrong sign in by" error that occurs
#' when start_val > end_val after rounding to the interval step.
#'
#' @param data Data frame containing the variable
#' @param variable Character string: variable name
#' @param breaks Either "auto" or a numeric vector of explicit break values
#' @param verbose Logical: print break information (default: FALSE)
#' @return Numeric vector of contour break values, or NULL if range insufficient
calculate_cdr_scale_contour_breaks <- function(data,
                                               variable,
                                               breaks  = "auto",
                                               verbose = FALSE) {
  
  if (!is.data.frame(data) || nrow(data) == 0) stop("data must be a non-empty data frame")
  if (!variable %in% names(data)) stop("Variable '", variable, "' not found in data")
  
  # Return explicit breaks as-is
  if (!identical(breaks, "auto")) {
    if (!is.numeric(breaks) || length(breaks) < 2) {
      stop("breaks must be 'auto' or a numeric vector with at least 2 values")
    }
    return(breaks)
  }
  
  # Automatic break calculation
  var_values <- data[[variable]][!is.na(data[[variable]])]
  if (length(var_values) == 0) return(NULL)
  
  data_range <- diff(range(var_values))
  
  # Return NULL when range is effectively zero
  if (data_range <= 0) return(NULL)
  
  interval <- if      (data_range > 500) 100
  else if (data_range > 200)  50
  else if (data_range > 100)  25
  else if (data_range > 50)   10
  else if (data_range > 20)    5
  else if (data_range > 10)    2
  else                         1
  
  min_val   <- min(var_values)
  max_val   <- max(var_values)
  start_val <- ceiling(min_val / interval) * interval
  end_val   <- floor(max_val  / interval) * interval
  
  # Guard: rounding can produce start_val > end_val when the data range is
  # smaller than one interval step. Return NULL so the contour layer is skipped
  # gracefully rather than passing a bad sequence to seq().
  if (start_val > end_val) return(NULL)
  
  auto_breaks <- seq(start_val, end_val, by = interval)
  
  # Require at least 2 breaks for a meaningful contour layer
  if (length(auto_breaks) < 2) return(NULL)
  
  if (verbose) {
    cat(sprintf("Auto contour breaks for %s: %d breaks (interval = %g)\n",
                variable, length(auto_breaks), interval))
  }
  
  return(auto_breaks)
}


#' @title Calculate Gradients for CDR Scale Arrow Visualisation
#' @description
#' Computes directional gradients (dx in the K direction, dy in the r direction)
#' and gradient magnitude using central differences.
#'
#' @param data Data frame containing CDR scale sensitivity results
#' @param variable Character string: variable name for gradient calculation
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return Data frame with additional columns dx, dy, gradient_mag
calculate_cdr_scale_gradients <- function(data,
                                          variable,
                                          verbose = TRUE) {
  
  if (!is.data.frame(data) || nrow(data) == 0) stop("data must be a non-empty data frame")
  if (!variable %in% names(data)) stop("Variable '", variable, "' not found in data")
  
  required_cols <- c("scenario_short", "K", "r")
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (verbose) cat(sprintf("Calculating gradients for variable: %s\n", variable))
  
  # dx: gradient in K direction
  data_with_gradients <- data %>%
    group_by(scenario_short, r) %>%
    arrange(K) %>%
    mutate(
      dx = (lead(.data[[variable]]) - lag(.data[[variable]])) /
        (lead(K) - lag(K))
    ) %>%
    ungroup()
  
  # dy: gradient in r direction
  data_with_gradients <- data_with_gradients %>%
    group_by(scenario_short, K) %>%
    arrange(r) %>%
    mutate(
      dy = (lead(.data[[variable]]) - lag(.data[[variable]])) /
        (lead(r) - lag(r))
    ) %>%
    ungroup()
  
  data_with_gradients <- data_with_gradients %>%
    mutate(gradient_mag = sqrt(dx^2 + dy^2))
  
  if (verbose) {
    n_valid <- sum(!is.na(data_with_gradients$dx) & !is.na(data_with_gradients$dy))
    n_total <- nrow(data_with_gradients)
    cat(sprintf("Gradient calculation complete: %d valid, %d missing (edge points)\n",
                n_valid, n_total - n_valid))
  }
  
  return(data_with_gradients)
}


#' @title Create Base CDR Scale Heatmap Panel
#' @description
#' Creates a single geom_tile heatmap panel for one variable and one SSP
#' scenario, with K on the x-axis and r on the y-axis.
#'
#' @param scenario_data Data frame filtered to a single scenario_short value
#' @param variable Character string: variable to fill by
#' @param variable_limits Numeric vector c(min, max) for the colour scale
#' @param palette_info List with 'option' and 'direction'
#' @param variable_label Character string for the legend title
#' @param scenario_name Character string: scenario label for panel title
#' @param show_title Logical: show scenario name as panel title (default: TRUE)
#' @param show_x_label Logical: show x-axis label (default: TRUE)
#' @param show_y_label Logical: show y-axis label (default: TRUE)
#' @param x_label Character string for x-axis
#' @param y_label Character string for y-axis
#' @param title_text Optional override for panel title text
#' @param theme_object ggplot2 theme from get_cdr_scale_theme()
#' @return ggplot object
create_cdr_scale_base_heatmap <- function(scenario_data,
                                          variable,
                                          variable_limits,
                                          palette_info,
                                          variable_label,
                                          scenario_name,
                                          show_title   = TRUE,
                                          show_x_label = TRUE,
                                          show_y_label = TRUE,
                                          x_label      = "K (GtCO\u2082/year)",
                                          y_label      = "r",
                                          title_text   = NULL,
                                          theme_object = get_cdr_scale_theme()) {
  
  if (!is.data.frame(scenario_data) || nrow(scenario_data) == 0) {
    stop("scenario_data must be a non-empty data frame")
  }
  required_cols <- c("K", "r", variable)
  missing_cols  <- setdiff(required_cols, names(scenario_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (length(variable_limits) != 2) {
    stop("variable_limits must be a numeric vector of length 2: c(min, max)")
  }
  
  plot_title   <- if (show_title)   { if (!is.null(title_text)) title_text else scenario_name } else ""
  x_axis_label <- if (show_x_label) x_label else ""
  y_axis_label <- if (show_y_label) y_label else ""
  
  p <- ggplot(scenario_data,
              aes(x = K, y = r, fill = .data[[variable]])) +
    geom_tile(color = NA) +
    scale_fill_viridis_c(
      name      = variable_label,
      option    = palette_info$option,
      direction = palette_info$direction,
      limits    = variable_limits,
      oob       = scales::squish,
      labels    = scales::comma
    ) +
    labs(
      title = plot_title,
      x     = x_axis_label,
      y     = y_axis_label
    ) +
    theme_object
  
  return(p)
}


#' @title Add Contour Lines to CDR Scale Heatmap
#'
#' @param plot_object ggplot object to add contours to
#' @param scenario_data Data frame for the scenario
#' @param variable Character string: variable for contour z-values
#' @param contour_breaks Numeric vector of break values, or NULL for automatic
#' @param contour_color Colour of contour lines (default: "white")
#' @param contour_alpha Transparency 0-1 (default: 0.6)
#' @param contour_linewidth Line width in mm (default: 0.3)
#' @return ggplot object with contours added
add_cdr_scale_contours <- function(plot_object,
                                   scenario_data,
                                   variable,
                                   contour_breaks    = NULL,
                                   contour_color     = "white",
                                   contour_alpha     = 0.6,
                                   contour_linewidth = 0.3) {
  
  if (!inherits(plot_object, "gg")) stop("plot_object must be a ggplot object")
  if (!is.data.frame(scenario_data) || nrow(scenario_data) == 0) {
    stop("scenario_data must be a non-empty data frame")
  }
  if (!variable %in% names(scenario_data)) {
    stop("Variable '", variable, "' not found in scenario_data")
  }
  if (contour_alpha < 0 || contour_alpha > 1) stop("contour_alpha must be between 0 and 1")
  
  if (is.null(contour_breaks)) {
    plot_object +
      geom_contour(
        data        = scenario_data,
        aes(x = K, y = r, z = .data[[variable]]),
        color       = contour_color,
        alpha       = contour_alpha,
        linewidth   = contour_linewidth,
        inherit.aes = FALSE
      )
  } else {
    if (!is.numeric(contour_breaks) || length(contour_breaks) < 2) {
      stop("contour_breaks must be NULL or a numeric vector with at least 2 values")
    }
    plot_object +
      geom_contour(
        data        = scenario_data,
        aes(x = K, y = r, z = .data[[variable]]),
        breaks      = contour_breaks,
        color       = contour_color,
        alpha       = contour_alpha,
        linewidth   = contour_linewidth,
        inherit.aes = FALSE
      )
  }
}


#' @title Add Gradient Arrow Vector Field to CDR Scale Heatmap
#'
#' @param plot_object ggplot object to add arrows to
#' @param scenario_data Data frame with columns K, r, dx, dy, gradient_mag
#' @param arrow_scale Numeric scaling factor for arrow length (default: 3.0)
#' @param arrow_skip Integer: skip every Nth grid point (default: 1)
#' @param min_magnitude Minimum gradient magnitude to display (default: 0)
#' @param arrow_size Line width in mm (default: 0.5)
#' @param arrow_alpha Transparency 0-1 (default: 0.7)
#' @param mag_limits Numeric vector c(min, max) for gradient magnitude scale
#' @param show_mag_legend Logical: show magnitude legend (default: FALSE)
#' @return ggplot object with arrows added
add_cdr_scale_arrows <- function(plot_object,
                                 scenario_data,
                                 arrow_scale     = 3.0,
                                 arrow_skip      = 1,
                                 min_magnitude   = 0,
                                 arrow_size      = 0.5,
                                 arrow_alpha     = 0.7,
                                 mag_limits      = NULL,
                                 show_mag_legend = FALSE) {
  
  if (!inherits(plot_object, "gg")) stop("plot_object must be a ggplot object")
  if (!requireNamespace("metR", quietly = TRUE)) {
    stop("Package 'metR' is required for arrow visualisation.")
  }
  
  required_cols <- c("K", "r", "dx", "dy", "gradient_mag")
  missing_cols  <- setdiff(required_cols, names(scenario_data))
  if (length(missing_cols) > 0) {
    stop("Missing required gradient columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (is.null(mag_limits)) mag_limits <- range(scenario_data$gradient_mag, na.rm = TRUE)
  
  legend_pos <- if (show_mag_legend) "right" else "none"
  
  plot_object +
    metR::geom_arrow(
      data         = scenario_data,
      aes(x        = K,
          y        = r,
          dx       = dx,
          dy       = dy,
          mag      = gradient_mag,
          color    = gradient_mag),
      scale        = arrow_scale,
      skip.x       = arrow_skip,
      skip.y       = arrow_skip,
      min.mag      = min_magnitude,
      size         = arrow_size,
      alpha        = arrow_alpha,
      preserve.dir = TRUE,
      arrow.type   = "closed",
      inherit.aes  = FALSE
    ) +
    scale_color_viridis_c(
      name   = "Gradient\nMagnitude",
      option = "plasma",
      limits = mag_limits,
      guide  = if (show_mag_legend) guide_colorbar() else "none"
    ) +
    theme(legend.position = legend_pos)
}


#' @title Add Infeasible Markers to CDR Scale Heatmap
#'
#' @param plot_object ggplot object
#' @param scenario_data Data frame containing a feasible column
#' @return ggplot object with red x markers on infeasible cells (unchanged if none)
add_cdr_scale_infeasible_markers <- function(plot_object, scenario_data) {
  
  infeasible_data <- scenario_data %>% filter(!feasible)
  
  if (nrow(infeasible_data) > 0) {
    plot_object +
      geom_point(
        data        = infeasible_data,
        aes(x = K, y = r),
        shape       = 4,
        size        = 1,
        color       = "red",
        alpha       = 0.8,
        stroke      = 0.8,
        inherit.aes = FALSE
      )
  } else {
    plot_object
  }
}


#' @title Extract Legend from CDR Scale Plot
#'
#' @param plot_object ggplot object with a visible legend
#' @param legend_position Character: "right" or "bottom"
#' @return Legend grob object
extract_cdr_scale_legend <- function(plot_object, legend_position = "right") {
  
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for legend extraction.")
  }
  if (!inherits(plot_object, "gg")) stop("plot_object must be a ggplot object")
  
  cowplot::get_legend(plot_object)
}


# ==============================================================================
# Section 4: Plot Grid and Dashboard Assembly
# ==============================================================================

#' @title Create Grid of CDR Scale Plots for All Scenarios and Variables
#' @description
#' Creates all individual heatmap panels for the requested variables and
#' scenarios, returning a nested list indexed by variable then scenario.
#'
#' @param data Prepared data frame from prepare_cdr_scale_data()
#' @param variables Character vector of variable names to plot
#' @param variable_limits Named list of limits from calculate_cdr_scale_limits()
#' @param palette_info Named list of palette specifications per variable
#' @param variable_labels Named list of legend labels per variable
#' @param add_contours Logical: add contour lines (default: TRUE)
#' @param contour_breaks Either "auto" or named list of breaks per variable
#' @param contour_alpha Transparency of contours (default: 0.6)
#' @param add_arrows Logical: add gradient arrow field (default: FALSE)
#' @param arrow_scale Numeric arrow scaling factor (default: 3.0)
#' @param arrow_skip Integer for arrow thinning (default: 1)
#' @param min_magnitude Minimum gradient magnitude for arrows (default: 0)
#' @param arrow_size Arrow line width in mm (default: 0.5)
#' @param arrow_alpha Arrow transparency (default: 0.7)
#' @param mag_limits Numeric c(min, max) for gradient magnitude scale
#' @param show_infeasible Logical: mark infeasible combinations (default: TRUE)
#' @param verbose Logical: print progress messages (default: TRUE)
#'
#' @return Nested named list: list(variable = list(scenario_short = ggplot))
create_cdr_scale_plot_grid <- function(data,
                                       variables,
                                       variable_limits,
                                       palette_info,
                                       variable_labels,
                                       add_contours    = TRUE,
                                       contour_breaks  = "auto",
                                       contour_alpha   = 0.6,
                                       add_arrows      = FALSE,
                                       arrow_scale     = 3.0,
                                       arrow_skip      = 1,
                                       min_magnitude   = 0,
                                       arrow_size      = 0.5,
                                       arrow_alpha     = 0.7,
                                       mag_limits      = NULL,
                                       show_infeasible = TRUE,
                                       verbose         = TRUE) {
  
  n_variables  <- length(variables)
  multi_row    <- n_variables > 1
  theme_object <- get_cdr_scale_theme(multi_row = multi_row)
  
  scenarios_present <- intersect(SSP_SCENARIO_ORDER_SCALE,
                                 unique(data$scenario_short))
  
  panel_counter   <- 0
  panel_labels    <- LETTERS
  panel_grid_list <- list()
  
  for (v_idx in seq_along(variables)) {
    variable <- variables[v_idx]
    
    if (verbose) cat(sprintf("  Building panels for variable: %s\n", variable))
    
    var_limits <- if ("shared" %in% names(variable_limits)) {
      variable_limits$shared
    } else {
      variable_limits[[variable]]
    }
    
    if (identical(contour_breaks, "auto")) {
      var_breaks <- "auto"
    } else if (is.list(contour_breaks) && variable %in% names(contour_breaks)) {
      var_breaks <- contour_breaks[[variable]]
    } else {
      var_breaks <- "auto"
    }
    
    panel_grid_list[[variable]] <- list()
    
    for (s_idx in seq_along(scenarios_present)) {
      scenario_name <- scenarios_present[s_idx]
      
      panel_counter <- panel_counter + 1
      panel_label   <- if (panel_counter <= length(panel_labels)) {
        panel_labels[panel_counter]
      } else {
        as.character(panel_counter)
      }
      
      # Axis label visibility rules
      show_x <- (v_idx == n_variables)   # x label on bottom row only
      show_y <- (s_idx == 1)             # y label on leftmost column only
      # FIX: renamed from show_title to show_panel_title to avoid shadowing the
      # show_title parameter of create_cdr_scale_base_heatmap()
      show_panel_title <- (v_idx == 1)   # scenario title on top row only
      
      scenario_data <- data %>% filter(scenario_short == scenario_name)
      
      if (nrow(scenario_data) == 0) {
        if (verbose) cat(sprintf("    Skipping %s (no data)\n", scenario_name))
        next
      }
      
      # Create base heatmap
      p <- create_cdr_scale_base_heatmap(
        scenario_data   = scenario_data,
        variable        = variable,
        variable_limits = var_limits,
        palette_info    = palette_info[[variable]],
        variable_label  = variable_labels[[variable]],
        scenario_name   = scenario_name,
        show_title      = show_panel_title,
        show_x_label    = show_x,
        show_y_label    = show_y,
        theme_object    = theme_object
      )
      
      # Add panel label
      # FIX: use face = "bold" not fontface = "bold" (correct element_text arg)
      p <- p + labs(tag = paste0(panel_label, ")")) +
        theme(plot.tag = element_text(size = if (multi_row) 8 else 10,
                                      face = "bold",
                                      hjust = 0))
      
      # Add contours
      if (add_contours) {
        contour_breaks_calc <- if (identical(var_breaks, "auto")) {
          calculate_cdr_scale_contour_breaks(scenario_data, variable)
        } else {
          var_breaks
        }
        
        # FIX: calculate_cdr_scale_contour_breaks() now returns NULL when the
        # range is too small; only add the layer when we have >= 2 valid breaks
        if (!is.null(contour_breaks_calc) && length(contour_breaks_calc) >= 2) {
          p <- add_cdr_scale_contours(
            plot_object    = p,
            scenario_data  = scenario_data,
            variable       = variable,
            contour_breaks = contour_breaks_calc,
            contour_alpha  = contour_alpha
          )
        } else {
          if (verbose) {
            cat(sprintf("    Panel %s: contours skipped (insufficient range)\n",
                        panel_label))
          }
        }
      }
      
      # Add gradient arrows
      if (add_arrows) {
        required_grad_cols <- c("dx", "dy", "gradient_mag")
        if (all(required_grad_cols %in% names(scenario_data))) {
          arrow_data <- scenario_data %>%
            filter(!is.na(dx), !is.na(dy), !is.na(gradient_mag))
          if (nrow(arrow_data) > 0) {
            p <- add_cdr_scale_arrows(
              plot_object   = p,
              scenario_data = arrow_data,
              arrow_scale   = arrow_scale,
              arrow_skip    = arrow_skip,
              min_magnitude = min_magnitude,
              arrow_size    = arrow_size,
              arrow_alpha   = arrow_alpha,
              mag_limits    = mag_limits
            )
          }
        } else {
          if (verbose) {
            warning(sprintf("Cannot add arrows for %s: missing gradient columns",
                            scenario_name))
          }
        }
      }
      
      # Add infeasible markers
      if (show_infeasible && "feasible" %in% names(scenario_data)) {
        p <- add_cdr_scale_infeasible_markers(p, scenario_data)
      }
      
      panel_grid_list[[variable]][[scenario_name]] <- p
      
      if (verbose) {
        cat(sprintf("    Panel %s: %s / %s complete\n",
                    panel_label, variable, scenario_name))
      }
    } # end scenario loop
  } # end variable loop
  
  return(panel_grid_list)
}


#' @title Assemble CDR Scale Dashboard with Patchwork
#' @description
#' Assembles individual panels from create_cdr_scale_plot_grid() into a complete
#' patchwork dashboard with an extracted legend.
#'
#' @param panel_grid Nested list from create_cdr_scale_plot_grid()
#' @param legend_grob Legend grob from extract_cdr_scale_legend()
#' @param variables Character vector of variable names (defines row order)
#' @param scenarios Character vector of scenario_short values (defines column order)
#' @param legend_position Character: "right" or "bottom"
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return patchwork object
assemble_cdr_scale_dashboard <- function(panel_grid,
                                         legend_grob,
                                         variables,
                                         scenarios,
                                         legend_position = "right",
                                         verbose         = TRUE) {
  
  n_variables <- length(variables)
  n_scenarios <- length(scenarios)
  
  if (verbose) cat("Assembling dashboard panels...\n")
  
  all_panels <- list()
  for (variable in variables) {
    for (scenario in scenarios) {
      panel <- panel_grid[[variable]][[scenario]]
      if (!is.null(panel)) {
        all_panels <- c(all_panels, list(panel))
      } else {
        all_panels <- c(all_panels, list(patchwork::plot_spacer()))
      }
    }
  }
  
  if (length(all_panels) == 0) stop("No panels were created. Check input data and variable names.")
  
  # ============================================================================
  # Single-variable layout: 3 x 2 grid with legend on right
  # ============================================================================
  if (n_variables == 1) {
    
    while (length(all_panels) < 6) {
      all_panels <- c(all_panels, list(patchwork::plot_spacer()))
    }
    
    # FIX: removed guides = "collect" — incompatible with manually extracted legend
    main_grid <- patchwork::wrap_plots(all_panels, ncol = 3, nrow = 2) &
      theme(plot.margin = unit(c(1, 2, 1, 2), "mm"))
    
    combined_plot <- main_grid | legend_grob
    combined_plot <- combined_plot + patchwork::plot_layout(widths = c(1, 0.08))
    
    # ============================================================================
    # Multi-variable layout: N x 5 grid with legend on bottom
    # ============================================================================
  } else {
    
    # FIX: removed guides = "collect" — incompatible with manually extracted legend
    main_grid <- patchwork::wrap_plots(all_panels, ncol = n_scenarios,
                                       nrow = n_variables) &
      theme(plot.margin = unit(c(1, 2, 3, 2), "mm"))
    
    combined_plot <- main_grid / legend_grob
    height_ratios <- c(rep(1, n_variables), 0.12)
    combined_plot <- combined_plot + patchwork::plot_layout(heights = height_ratios) +
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  
  if (verbose) cat("Dashboard assembly complete\n")
  
  return(combined_plot)
}


#' @title Save CDR Scale Dashboard to File
#'
#' @param plot_object patchwork or ggplot object to save
#' @param filename Character string for output filename (NULL for auto)
#' @param width Numeric width in mm (default: 297)
#' @param height Numeric height in mm (default: 210)
#' @param output_dir Character string for output directory (default: "figs")
#' @param verbose Logical: print save information (default: TRUE)
#' @return Invisibly returns the filepath
save_cdr_scale_dashboard <- function(plot_object,
                                     filename   = NULL,
                                     width      = 297,
                                     height     = 210,
                                     output_dir = "figs",
                                     verbose    = TRUE) {
  
  if (!inherits(plot_object, "patchwork") && !inherits(plot_object, "gg")) {
    stop("plot_object must be a patchwork or ggplot object")
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required for file path handling.")
  }
  
  if (is.null(filename)) {
    filename <- paste0("cdr_scale_sensitivity_dashboard_",
                       format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    if (verbose) cat(sprintf("Generated filename: %s\n", filename))
  }
  
  if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".pdf")
  }
  
  filepath        <- here::here(output_dir, filename)
  output_dir_full <- here::here(output_dir)
  if (!dir.exists(output_dir_full)) dir.create(output_dir_full, recursive = TRUE)
  
  if (verbose) {
    cat(sprintf("Saving dashboard to: %s\n", filepath))
    cat(sprintf("Dimensions: %d x %d mm\n", width, height))
  }
  
  ggsave(
    filename = filepath,
    plot     = plot_object,
    width    = width,
    height   = height,
    units    = "mm",
    device   = cairo_pdf,
    dpi      = 300
  )
  
  if (verbose) cat(sprintf("Dashboard saved successfully to: %s\n", filepath))
  invisible(filepath)
}


# ==============================================================================
# Section 5: Main Dashboard Function
# ==============================================================================

#' @title Create CDR Scale Sensitivity Dashboard
#' @description
#' Main user-facing function. Orchestrates the full workflow from data
#' preparation through final output.
#'
#' @param sensitivity_results Results object from run_cdr_scale_sensitivity()
#' @param variables Character vector of variable names to plot
#' @param color_palettes Optional named list of palette overrides (default: NULL)
#' @param custom_labels Optional named list of legend label overrides (default: NULL)
#' @param shared_scale Logical: shared colour scale across variables (default: FALSE)
#' @param use_scale_limits Logical: cap colour scale at a percentile (default: FALSE)
#' @param scale_limit_percentile Numeric percentile for capping (default: 95)
#' @param add_contours Logical: add contour lines (default: TRUE)
#' @param contour_breaks Either "auto" or named list of numeric vectors per variable
#' @param contour_alpha Transparency of contours (default: 0.6)
#' @param add_arrows Logical: add gradient vector field arrows (default: FALSE)
#' @param arrow_scale Numeric arrow scaling factor (default: 3.0)
#' @param arrow_skip Integer for arrow thinning (default: 1)
#' @param arrow_size Arrow line width in mm (default: 0.5)
#' @param arrow_alpha Arrow transparency (default: 0.7)
#' @param min_magnitude Minimum gradient magnitude for arrow display (default: 0)
#' @param show_infeasible Logical: mark infeasible K/r combinations (default: TRUE)
#' @param save_plot Logical: save dashboard to file (default: FALSE)
#' @param filename Character string for output filename (default: NULL for auto)
#' @param width Numeric width in mm (default: 297)
#' @param height Numeric height in mm (default: 210)
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return Complete patchwork dashboard object
create_cdr_scale_dashboard <- function(sensitivity_results,
                                       variables,
                                       color_palettes         = NULL,
                                       custom_labels          = NULL,
                                       shared_scale           = FALSE,
                                       use_scale_limits       = FALSE,
                                       scale_limit_percentile = 95,
                                       add_contours           = TRUE,
                                       contour_breaks         = "auto",
                                       contour_alpha          = 0.6,
                                       add_arrows             = FALSE,
                                       arrow_scale            = 3.0,
                                       arrow_skip             = 1,
                                       arrow_size             = 0.5,
                                       arrow_alpha            = 0.7,
                                       min_magnitude          = 0,
                                       show_infeasible        = TRUE,
                                       save_plot              = FALSE,
                                       filename               = NULL,
                                       width                  = 297,
                                       height                 = 210,
                                       verbose                = TRUE) {
  
  # Step 1: Validate inputs and check packages
  if (verbose) cat("\n=== CDR SCALE SENSITIVITY DASHBOARD CREATION ===\n\n")
  if (verbose) cat("Step 1: Validating inputs and checking packages\n")
  
  required_packages <- c("ggplot2", "dplyr", "patchwork", "viridis",
                         "cowplot", "here", "scales")
  if (add_arrows) required_packages <- c(required_packages, "metR")
  
  missing_packages <- required_packages[
    !sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Required packages not installed: ", paste(missing_packages, collapse = ", "))
  }
  if (length(variables) == 0) stop("variables must contain at least one variable name")
  
  if (verbose) {
    cat(sprintf("  Variables: %s\n", paste(variables, collapse = ", ")))
    cat(sprintf("  Features: contours=%s, arrows=%s, infeasible=%s\n",
                add_contours, add_arrows, show_infeasible))
  }
  
  # Step 2: Prepare data
  if (verbose) cat("\nStep 2: Preparing data\n")
  plot_data <- prepare_cdr_scale_data(sensitivity_results, variables, verbose)
  
  # Step 3: Gradients (if arrows requested)
  if (add_arrows) {
    if (verbose) cat("\nStep 3: Calculating gradients for arrows\n")
    plot_data  <- calculate_cdr_scale_gradients(plot_data, variables[1], verbose)
    mag_limits <- range(plot_data$gradient_mag, na.rm = TRUE)
  } else {
    mag_limits <- NULL
    if (verbose) cat("\nStep 3: Skipping gradient calculation (arrows not requested)\n")
  }
  
  # Step 4: Colour scale limits
  if (verbose) cat("\nStep 4: Calculating colour scale limits\n")
  variable_limits <- calculate_cdr_scale_limits(
    data                   = plot_data,
    variables              = variables,
    shared_scale           = shared_scale,
    use_scale_limits       = use_scale_limits,
    scale_limit_percentile = scale_limit_percentile,
    verbose                = verbose
  )
  
  # Step 5: Palettes and labels
  if (verbose) cat("\nStep 5: Configuring palettes and labels\n")
  use_single_line <- (length(variables) > 1)
  palette_info    <- list()
  variable_labels <- list()
  for (var in variables) {
    palette_info[[var]]    <- get_cdr_scale_palette(var, color_palettes)
    variable_labels[[var]] <- get_cdr_scale_label(var, custom_labels, use_single_line)
  }
  if (verbose) cat(sprintf("  Configured %d variable(s)\n", length(variables)))
  
  # Step 6: Build panel grid
  if (verbose) cat("\nStep 6: Creating plot grid\n")
  panel_grid <- create_cdr_scale_plot_grid(
    data            = plot_data,
    variables       = variables,
    variable_limits = variable_limits,
    palette_info    = palette_info,
    variable_labels = variable_labels,
    add_contours    = add_contours,
    contour_breaks  = contour_breaks,
    contour_alpha   = contour_alpha,
    add_arrows      = add_arrows,
    arrow_scale     = arrow_scale,
    arrow_skip      = arrow_skip,
    min_magnitude   = min_magnitude,
    arrow_size      = arrow_size,
    arrow_alpha     = arrow_alpha,
    mag_limits      = mag_limits,
    show_infeasible = show_infeasible,
    verbose         = verbose
  )
  
  # Step 7: Extract legend
  if (verbose) cat("\nStep 7: Extracting legend\n")
  first_variable  <- variables[1]
  first_scenario  <- intersect(SSP_SCENARIO_ORDER_SCALE,
                               unique(plot_data$scenario_short))[1]
  legend_data     <- plot_data %>% filter(scenario_short == first_scenario)
  legend_position <- if (length(variables) == 1) "right" else "bottom"
  
  temp_plot <- create_cdr_scale_base_heatmap(
    scenario_data   = legend_data,
    variable        = first_variable,
    variable_limits = if ("shared" %in% names(variable_limits)) {
      variable_limits$shared
    } else {
      variable_limits[[first_variable]]
    },
    palette_info    = palette_info[[first_variable]],
    variable_label  = variable_labels[[first_variable]],
    scenario_name   = first_scenario,
    show_title      = FALSE,
    theme_object    = get_cdr_scale_theme(multi_row = length(variables) > 1)
  )
  
  if (length(variables) == 1) {
    temp_plot <- temp_plot + theme(legend.position = "right")
  } else {
    temp_plot <- temp_plot + theme(
      legend.position   = "bottom",
      legend.key.width  = unit(2, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.title      = element_text(size = 6),
      legend.text       = element_text(size = 5.5)
    )
  }
  
  legend_grob <- extract_cdr_scale_legend(temp_plot, legend_position)
  if (verbose) cat(sprintf("  Legend position: %s\n", legend_position))
  
  # Step 8: Assemble dashboard
  if (verbose) cat("\nStep 8: Assembling dashboard\n")
  scenarios_present <- intersect(SSP_SCENARIO_ORDER_SCALE,
                                 unique(plot_data$scenario_short))
  
  dashboard <- assemble_cdr_scale_dashboard(
    panel_grid      = panel_grid,
    legend_grob     = legend_grob,
    variables       = variables,
    scenarios       = scenarios_present,
    legend_position = legend_position,
    verbose         = verbose
  )
  
  # Step 9: Save if requested
  if (save_plot) {
    if (verbose) cat("\nStep 9: Saving dashboard to file\n")
    filepath <- save_cdr_scale_dashboard(dashboard, filename, width, height, verbose = verbose)
    if (verbose) cat(sprintf("\n=== DASHBOARD CREATION COMPLETE ===\nSaved to: %s\n\n", filepath))
    return(invisible(dashboard))
  } else {
    if (verbose) cat("\nStep 9: Skipping file save (save_plot = FALSE)\n=== DASHBOARD CREATION COMPLETE ===\n\n")
    return(dashboard)
  }
  
} # Close create_cdr_scale_dashboard


# ==============================================================================
# Section 6: Convenience Wrappers
# ==============================================================================

#' @title Create CDR Scale Temperature Dashboard
#' @description
#' Convenience wrapper pre-configured for temperature outcome analysis.
#' Plots peak_temperature and years_above_1p5. All arguments are explicit —
#' no ... passthrough — which prevents duplicate argument errors.
#'
#' @param sensitivity_results Results object from run_cdr_scale_sensitivity()
#' @param contour_alpha Transparency of contour lines (default: 0.6)
#' @param show_infeasible Logical: mark infeasible combinations (default: TRUE)
#' @param use_scale_limits Logical: cap colour scale at percentile (default: FALSE)
#' @param scale_limit_percentile Numeric percentile for capping (default: 95)
#' @param save_plot Logical: save to file (default: FALSE)
#' @param filename Character string for output filename (default: NULL for auto)
#' @param width Numeric width in mm (default: 297)
#' @param height Numeric height in mm (default: 210)
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return Patchwork dashboard object
create_cdr_scale_temperature_dashboard <- function(sensitivity_results,
                                                   contour_alpha          = 0.6,
                                                   show_infeasible        = TRUE,
                                                   use_scale_limits       = FALSE,
                                                   scale_limit_percentile = 95,
                                                   save_plot              = FALSE,
                                                   filename               = NULL,
                                                   width                  = 297,
                                                   height                 = 210,
                                                   verbose                = TRUE) {
  create_cdr_scale_dashboard(
    sensitivity_results    = sensitivity_results,
    variables              = c("peak_temperature", "years_above_1p5"),
    shared_scale           = FALSE,
    use_scale_limits       = use_scale_limits,
    scale_limit_percentile = scale_limit_percentile,
    add_contours           = TRUE,
    contour_alpha          = contour_alpha,
    add_arrows             = FALSE,
    show_infeasible        = show_infeasible,
    save_plot              = save_plot,
    filename               = filename,
    width                  = width,
    height                 = height,
    verbose                = verbose
  )
}

#' @title Create CDR Scale Temperature Dashboards (Separate Files)
#' @description
#' Creates two separate dashboard files — one for peak_temperature and one for
#' years_above_1p5 — each with the 3x2 single-variable layout (5 SSP panels
#' plus legend on right), matching the delayed deployment dashboard style.
#'
#' @param sensitivity_results Results object from run_cdr_scale_sensitivity()
#' @param contour_alpha Transparency of contour lines (default: 0.6)
#' @param show_infeasible Logical: mark infeasible combinations (default: TRUE)
#' @param use_scale_limits Logical: cap colour scale at percentile (default: FALSE)
#' @param scale_limit_percentile Numeric percentile for capping (default: 95)
#' @param filename_temp Filename for peak temperature plot (default: NULL for auto)
#' @param filename_years Filename for years above 1.5 plot (default: NULL for auto)
#' @param width Numeric width in mm (default: 297)
#' @param height Numeric height in mm (default: 210)
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return Invisibly returns a list with both dashboard objects:
#'   list(peak_temperature = ..., years_above_1p5 = ...)
create_cdr_scale_temperature_dashboards <- function(sensitivity_results,
                                                    contour_alpha          = 0.6,
                                                    show_infeasible        = TRUE,
                                                    use_scale_limits       = FALSE,
                                                    scale_limit_percentile = 95,
                                                    filename_temp          = NULL,
                                                    filename_years         = NULL,
                                                    width                  = 297,
                                                    height                 = 210,
                                                    verbose                = TRUE) {
  
  if (verbose) cat("Creating peak temperature dashboard...\n")
  
  dash_temp <- create_cdr_scale_dashboard(
    sensitivity_results    = sensitivity_results,
    variables              = "peak_temperature",
    shared_scale           = FALSE,
    use_scale_limits       = use_scale_limits,
    scale_limit_percentile = scale_limit_percentile,
    add_contours           = TRUE,
    contour_alpha          = contour_alpha,
    add_arrows             = FALSE,
    show_infeasible        = show_infeasible,
    save_plot              = TRUE,
    filename               = filename_temp,
    width                  = width,
    height                 = height,
    verbose                = verbose
  )
  
  if (verbose) cat("\nCreating years above 1.5°C dashboard...\n")
  
  dash_years <- create_cdr_scale_dashboard(
    sensitivity_results    = sensitivity_results,
    variables              = "years_above_1p5",
    shared_scale           = FALSE,
    use_scale_limits       = use_scale_limits,
    scale_limit_percentile = scale_limit_percentile,
    add_contours           = TRUE,
    contour_alpha          = contour_alpha,
    add_arrows             = FALSE,
    show_infeasible        = show_infeasible,
    save_plot              = TRUE,
    filename               = filename_years,
    width                  = width,
    height                 = height,
    verbose                = verbose
  )
  
  invisible(list(
    peak_temperature = dash_temp,
    years_above_1p5  = dash_years
  ))
}


#' @title Create CDR Scale Cost Dashboard
#' @description
#' Convenience wrapper pre-configured for cost outcome analysis.
#' Plots total_cost, abatement_cost, and temp_cost. All arguments are explicit —
#' no ... passthrough — which prevents duplicate argument errors.
#'
#' @param sensitivity_results Results object from run_cdr_scale_sensitivity()
#' @param shared_scale Logical: shared colour scale (default: TRUE)
#' @param contour_alpha Transparency of contour lines (default: 0.6)
#' @param show_infeasible Logical: mark infeasible combinations (default: TRUE)
#' @param use_scale_limits Logical: cap colour scale at percentile (default: FALSE)
#' @param scale_limit_percentile Numeric percentile for capping (default: 95)
#' @param save_plot Logical: save to file (default: FALSE)
#' @param filename Character string for output filename (default: NULL for auto)
#' @param width Numeric width in mm (default: 297)
#' @param height Numeric height in mm (default: 260, taller for 3 rows)
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return Patchwork dashboard object
create_cdr_scale_cost_dashboard <- function(sensitivity_results,
                                            shared_scale           = TRUE,
                                            contour_alpha          = 0.6,
                                            show_infeasible        = TRUE,
                                            use_scale_limits       = FALSE,
                                            scale_limit_percentile = 95,
                                            save_plot              = FALSE,
                                            filename               = NULL,
                                            width                  = 297,
                                            height                 = 260,
                                            verbose                = TRUE) {
  create_cdr_scale_dashboard(
    sensitivity_results    = sensitivity_results,
    variables              = c("total_cost", "abatement_cost", "temp_cost"),
    shared_scale           = shared_scale,
    use_scale_limits       = use_scale_limits,
    scale_limit_percentile = scale_limit_percentile,
    add_contours           = TRUE,
    contour_alpha          = contour_alpha,
    add_arrows             = FALSE,
    show_infeasible        = show_infeasible,
    save_plot              = save_plot,
    filename               = filename,
    width                  = width,
    height                 = height,
    verbose                = verbose
  )
}