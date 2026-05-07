# ==============================================================================
# Capacity Scenario Comparison Visualisation
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
# Last updated: March 2026
# ==============================================================================

#' @title Capacity Scenario Comparison Visualisation
#' @description
#' Creates a 10-panel figure comparing temperature trajectories and CDR strategies
#' across SSP scenarios and growth rate assumptions (slow, moderate, fast). 
#' Panels are faceted by SSP (5 columns), with rows for each variable. Within
#' each panel, the three growth rates are shown as solid, dashed, and dotted lines.
#'
#' @note All required libraries (ggplot2, dplyr, patchwork, purrr, scales, cowplot, here) 
#' must be loaded before using these functions.
#'
#' @author Nina Rynne
#' @date March 2026

# ============================================================================
# Global definitions
# ============================================================================

# Reuse scenario colours from scenario_comparison_visualisation.R
# Colors selected for accessibility and distinction between 5 SSP scenarios
scenario_colors <- c(
  "#00ADCF",  # Cyan    — SSP1
  "#173C66",  # Dark blue — SSP2
  "#F0E442",  # Yellow  — SSP3
  "#E71D25",  # Red     — SSP4
  "#951B1E"   # Dark red — SSP5
)

# SSP names used as panel column headers and for matching results lists
ssp_names <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")

# Growth rate line types and labels, ordered slow → moderate → fast
growth_rate_linetypes <- c(
  "slow"     = "dotted",
  "moderate" = "dashed",
  "fast"     = "solid"
)

growth_rate_labels <- c(
  "slow"     = "Slow",
  "moderate" = "Moderate",
  "fast"     = "Fast"
)

# Reuse base theme from scenario_comparison_visualisation.R
scenario_comparison_theme <- theme_bw() +
  theme(
    text                = element_text(size = 10),
    plot.title          = element_text(size = 10),
    axis.title          = element_text(size = 9),
    axis.text           = element_text(size = 8),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.title        = element_text(size = 9),
    legend.text         = element_text(size = 8),
    legend.position     = "none"
  )

# ============================================================================
# Helper functions
# ============================================================================

#' @title Clean Scenario Names for Display
#' @description
#' Removes "-Baseline" suffix from scenario names for cleaner plot legends.
#' Reused from scenario_comparison_visualisation.R.
#' 
#' @param scenario_names Character vector of scenario names
#' @return Character vector with cleaned names
clean_scenario_names <- function(scenario_names) {
  gsub("-Baseline$", "", scenario_names)
}

#' @title Extract Variable Data Across Growth Rates for One SSP
#' @description
#' Pulls a single numeric variable from the scenario_results component of each
#' growth rate entry, returning a tidy data frame ready for plotting.
#'
#' Expects capacity_results to be the object returned by
#' run_capacity_growth_comparison(), i.e. a named list keyed by growth rate
#' label, where each element is the full list returned by
#' run_scenario_comparison(). The SSP solution objects are therefore at:
#'   capacity_results[[rate]]$scenario_results[["SSP1-Baseline"]]
#'
#' @param capacity_results Named list as returned by
#'   run_capacity_growth_comparison().
#' @param ssp Character string identifying the SSP. Both "SSP1" and
#'   "SSP1-Baseline" are accepted — the function tries both.
#' @param variable Character string; name of the field within each solution
#'   object (e.g. "temperature_anomaly", "qty_remov").
#' @return Data frame with columns: growth_rate, years, value
extract_ssp_variable <- function(capacity_results, ssp, variable) {
  
  map_dfr(names(growth_rate_linetypes), function(rate) {
    
    # Navigate to the scenario_results component produced by
    # run_scenario_comparison() for this growth rate
    rate_entry <- capacity_results[[rate]]
    if (is.null(rate_entry)) {
      warning(sprintf("No results found for growth rate '%s'", rate))
      return(NULL)
    }
    
    scenario_results <- rate_entry$scenario_results
    if (is.null(scenario_results)) {
      warning(sprintf(
        "No $scenario_results component found for growth rate '%s'", rate))
      return(NULL)
    }
    
    # Match SSP with or without "-Baseline" suffix
    ssp_key <- if (!is.null(scenario_results[[ssp]])) {
      ssp
    } else if (!is.null(scenario_results[[paste0(ssp, "-Baseline")]])) {
      paste0(ssp, "-Baseline")
    } else {
      warning(sprintf(
        "No results found for SSP '%s' in growth rate '%s'", ssp, rate))
      return(NULL)
    }
    
    result <- scenario_results[[ssp_key]]
    
    data.frame(
      growth_rate = rate,
      years       = result$years,
      value       = result[[variable]]
    )
  })
}

#' @title Extract CDR Capacity Curve Data Across Growth Rates
#' @description
#' Evaluates the logistic capacity curve for each growth rate using the
#' parameters stored in run_info (g_initial, K, r_value, t_start), returning
#' a tidy data frame in the same shape as extract_ssp_variable() output.
#' The years vector is taken from the first available SSP result for each rate,
#' so the capacity curve spans the same time axis as the plotted data.
#'
#' @param capacity_results Named list as returned by
#'   run_capacity_growth_comparison().
#' @return Data frame with columns: growth_rate, years, value
extract_cdr_capacity <- function(capacity_results) {
  
  map_dfr(names(growth_rate_linetypes), function(rate) {
    
    rate_entry <- capacity_results[[rate]]
    if (is.null(rate_entry)) {
      warning(sprintf("No results found for growth rate '%s'", rate))
      return(NULL)
    }
    
    # Pull capacity parameters saved in run_info
    run_info <- rate_entry$run_info
    if (is.null(run_info$g_initial) || is.null(run_info$K) ||
        is.null(run_info$r_value)   || is.null(run_info$t_start)) {
      warning(sprintf(
        "run_info for growth rate '%s' is missing capacity parameters. ",
        "Re-run with the updated scenario_comparison_capacity.R.", rate))
      return(NULL)
    }
    
    g_initial <- run_info$g_initial
    K         <- run_info$K
    r         <- run_info$r_value
    t_start   <- run_info$t_start
    
    # Pre-compute suppression factor (matches make_logistic_from_zero)
    suppression_factor <- (K / g_initial) - 1
    
    # Get years vector from the first available SSP result for this rate
    scenario_results <- rate_entry$scenario_results
    first_result     <- scenario_results[[1]]
    years            <- first_result$years
    
    # Evaluate logistic formula: g(t) = K / (1 + suppression * exp(-r*(t-t_start)))
    # Return g_initial for years before t_start (matches make_logistic_from_zero)
    capacity_values <- ifelse(
      years < t_start,
      g_initial,
      K / (1 + suppression_factor * exp(-r * (years - t_start)))
    )
    
    data.frame(
      growth_rate = rate,
      years       = years,
      value       = capacity_values
    )
  })
}

# ============================================================================
# Individual panel functions
# ============================================================================

#' @title Plot Temperature Trajectories for One SSP — All Growth Rates
#' @description
#' Creates a single panel showing temperature anomaly over time for one SSP
#' across all three growth rates, with a 1.5°C reference line.
#'
#' @param capacity_results Named list with elements "slow", "moderate", "fast".
#' @param ssp Character string identifying the SSP, e.g. "SSP1".
#' @param ssp_colour Hex colour string for this SSP's lines.
#' @param y_limits Numeric vector of length 2 for shared y-axis limits.
#' @param show_y_axis Logical; show y-axis title and text (default TRUE).
#'   Set FALSE for all but the leftmost column to avoid repetition.
#' @param show_x_axis Logical; show x-axis title (default TRUE).
#'   Set FALSE for the temperature row when CDR row is directly below.
#' @return ggplot object
plot_capacity_temperature <- function(capacity_results,
                                      ssp,
                                      ssp_colour,
                                      y_limits    = c(NA, 2.8),
                                      show_y_axis = TRUE,
                                      show_x_axis = TRUE) {
  
  plot_data <- extract_ssp_variable(capacity_results, ssp, "temperature_anomaly")
  
  if (nrow(plot_data) == 0) {
    stop(sprintf("No valid temperature data found for SSP '%s'", ssp))
  }
  
  # Order growth_rate factor for consistent legend ordering
  plot_data$growth_rate <- factor(plot_data$growth_rate,
                                  levels = names(growth_rate_linetypes))
  
  p <- ggplot(plot_data, aes(x = years, y = value, linetype = growth_rate)) +
    geom_hline(yintercept = 1.5, linetype = "dashed", colour = "red",
               alpha = 0.7, linewidth = 0.5) +
    geom_line(colour = ssp_colour, linewidth = 0.8) +
    scale_linetype_manual(
      name   = "Growth rate",
      values = growth_rate_linetypes,
      labels = growth_rate_labels
    ) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      title = NULL,
      x     = if (show_x_axis) "Year" else NULL,
      y     = if (show_y_axis) "Temperature anomaly (°C)" else NULL
    ) +
    scenario_comparison_theme +
    theme(
      axis.title.y = if (!show_y_axis) element_blank() else element_text(size = 9),
      axis.text.y  = if (!show_y_axis) element_blank() else element_text(size = 8),
      axis.ticks.y = if (!show_y_axis) element_blank() else element_line()
    )
  
  return(p)
}

#' @title Plot CDR Strategies for One SSP — All Growth Rates
#' @description
#' Creates a single panel showing annual CDR (carbon dioxide removal) over time
#' for one SSP across all three growth rates.
#'
#' @param capacity_results Named list with elements "slow", "moderate", "fast".
#' @param ssp Character string identifying the SSP, e.g. "SSP1".
#' @param ssp_colour Hex colour string for this SSP's lines.
#' @param y_limits Numeric vector of length 2 for shared y-axis limits.
#' @param show_y_axis Logical; show y-axis title and text (default TRUE).
#' @param show_x_axis Logical; show x-axis title (default TRUE).
#' @param show_capacity_limit Logical; overlay the maximum CDR capacity curve
#'   for each growth rate as a black line matching the growth rate linetype
#'   (default FALSE).
#' @return ggplot object
plot_capacity_cdr <- function(capacity_results,
                              ssp,
                              ssp_colour,
                              y_limits           = c(0, NA),
                              show_y_axis        = TRUE,
                              show_x_axis        = TRUE,
                              show_capacity_limit = FALSE) {
  
  plot_data <- extract_ssp_variable(capacity_results, ssp, "qty_remov")
  
  if (nrow(plot_data) == 0) {
    stop(sprintf("No valid CDR data found for SSP '%s'", ssp))
  }
  
  plot_data$growth_rate <- factor(plot_data$growth_rate,
                                  levels = names(growth_rate_linetypes))
  
  p <- ggplot(plot_data, aes(x = years, y = value, linetype = growth_rate))
  
  # Add capacity lines first so they render behind the results lines
  if (show_capacity_limit) {
    capacity_data <- extract_cdr_capacity(capacity_results)
    capacity_data$growth_rate <- factor(capacity_data$growth_rate,
                                        levels = names(growth_rate_linetypes))
    p <- p +
      geom_line(data        = capacity_data,
                mapping     = aes(x = years, y = value, linetype = growth_rate),
                colour      = "grey70",
                linewidth   = 0.8,
                inherit.aes = FALSE)
  }
  
  p <- p +
    geom_line(colour = ssp_colour, linewidth = 0.8) +
    scale_linetype_manual(
      name   = "Growth rate",
      values = growth_rate_linetypes,
      labels = growth_rate_labels
    ) +
    scale_y_continuous(limits = y_limits) +
    labs(
      title = NULL,
      x     = if (show_x_axis) "Year" else NULL,
      y     = if (show_y_axis) expression("CDR (GtCO"[2]*"/yr)") else NULL
    ) +
    scenario_comparison_theme +
    theme(
      axis.title.y = if (!show_y_axis) element_blank() else element_text(size = 9),
      axis.text.y  = if (!show_y_axis) element_blank() else element_text(size = 8),
      axis.ticks.y = if (!show_y_axis) element_blank() else element_line()
    )
  
  return(p)
}

# ============================================================================
# Dashboard function
# ============================================================================

#' @title Create Capacity Scenario Comparison Dashboard
#' @description
#' Creates a 10-panel figure (2 rows × 5 columns) comparing temperature
#' trajectories and CDR strategies across SSP scenarios and growth rate
#' assumptions. Each column corresponds to one SSP; within each panel the
#' three growth rates are shown as solid (slow), dashed (moderate), and
#' dotted (fast) lines. A shared legend is placed to the right of the grid.
#'
#' Expects capacity_results to be the object returned by
#' run_capacity_growth_comparison(), i.e. a named list keyed by growth rate
#' label where each element is the full output of run_scenario_comparison().
#'
#' @param capacity_results Named list as returned by
#'   run_capacity_growth_comparison().
#' @param save_plot Logical; whether to save the plot to file (default FALSE).
#' @param filename Character; custom filename (default: auto-generated timestamp).
#'   Extension determines format — use ".pdf" or ".png".
#' @param verbose Logical; print progress messages (default TRUE).
#' @param width Plot width in mm (default 297, A4 landscape).
#' @param height Plot height in mm (default 180).
#'
#' @return Combined ggplot object containing the 10-panel dashboard.
#'
#' @examples
#' # Run the multi-rate comparison (see scenario_comparison_analysis_capacity chunk)
#' capacity_growth_results <- run_capacity_growth_comparison(
#'   parameter_df = parameter_df[1, ],
#'   emissions_df = emissions_df,
#'   economic_df  = economic_df,
#'   scenarios    = c("SSP1-Baseline", "SSP2-Baseline", "SSP3-Baseline",
#'                    "SSP4-Baseline", "SSP5-Baseline"),
#'   r_values     = c(slow = 0.07, moderate = 0.11, fast = 0.18)
#' )
#'
#' # Preview in RStudio viewer
#' fig <- create_capacity_scenario_comparison_dashboard(capacity_growth_results)
#' print(fig)
#'
#' # Save to figs/ directory
#' create_capacity_scenario_comparison_dashboard(
#'   capacity_results = capacity_growth_results,
#'   save_plot        = TRUE,
#'   filename         = "capacity_scenario_comparison.pdf"
#' )
#'
create_capacity_scenario_comparison_dashboard <- function(capacity_results,
                                                          save_plot           = FALSE,
                                                          filename            = NULL,
                                                          verbose             = TRUE,
                                                          show_capacity_limit = FALSE,
                                                          width               = 297,
                                                          height              = 180) {
  
  # Validate top-level structure
  required_rates <- names(growth_rate_linetypes)  # "slow", "moderate", "fast"
  missing_rates  <- setdiff(required_rates, names(capacity_results))
  if (length(missing_rates) > 0) {
    stop(sprintf(
      "capacity_results is missing growth rate(s): %s\nExpected names: %s",
      paste(missing_rates, collapse = ", "),
      paste(required_rates, collapse = ", ")
    ))
  }
  
  if (verbose) {
    cat("Creating capacity scenario comparison dashboard\n")
    cat("Growth rates found:", paste(names(capacity_results), collapse = ", "), "\n")
  }
  
  # ── Compute shared y-axis limits ──────────────────────────────────────────
  # Navigate through $scenario_results to reach the per-SSP solution objects
  
  all_temp <- map_dfr(ssp_names, function(ssp) {
    map_dfr(required_rates, function(rate) {
      result <- capacity_results[[rate]]$scenario_results[[ssp]] %||%
        capacity_results[[rate]]$scenario_results[[paste0(ssp, "-Baseline")]]
      if (!is.null(result)) data.frame(temp = result$temperature_anomaly)
    })
  })
  temp_y_limits <- c(
    floor(min(all_temp$temp, na.rm = TRUE) * 10) / 10,
    ceiling(max(all_temp$temp, na.rm = TRUE) * 10) / 10 + 0.1
  )
  
  all_cdr <- map_dfr(ssp_names, function(ssp) {
    map_dfr(required_rates, function(rate) {
      result <- capacity_results[[rate]]$scenario_results[[ssp]] %||%
        capacity_results[[rate]]$scenario_results[[paste0(ssp, "-Baseline")]]
      if (!is.null(result)) data.frame(cdr = result$qty_remov)
    })
  })
  
  # If capacity lines are shown, include their values in the y-axis upper bound
  # so the capacity curve is never clipped
  if (show_capacity_limit) {
    capacity_data    <- extract_cdr_capacity(capacity_results)
    cdr_upper        <- max(max(all_cdr$cdr, na.rm = TRUE),
                            max(capacity_data$value, na.rm = TRUE))
  } else {
    cdr_upper        <- max(all_cdr$cdr, na.rm = TRUE)
  }
  cdr_y_limits <- c(0, ceiling(cdr_upper * 1.05))
  
  if (verbose) {
    cat(sprintf("Temperature y-axis: %.1f – %.1f °C\n",
                temp_y_limits[1], temp_y_limits[2]))
    cat(sprintf("CDR y-axis: 0 – %.0f GtCO2/yr\n", cdr_y_limits[2]))
  }
  
  # ── Build individual panels ───────────────────────────────────────────────
  # Panel letters a-e for CDR row, f-j for temperature row.
  # Subtitles are plain black text — colour identity is carried by line colour.
  
  panel_letters_cdr  <- letters[1:5]   # a, b, c, d, e
  panel_letters_temp <- letters[6:10]  # f, g, h, i, j
  
  temp_panels <- vector("list", length(ssp_names))
  cdr_panels  <- vector("list", length(ssp_names))
  
  for (i in seq_along(ssp_names)) {
    ssp      <- ssp_names[i]
    ssp_col  <- scenario_colors[i]
    show_y   <- (i == 1)   # y-axis labels only on leftmost column
    
    p_temp <- plot_capacity_temperature(
      capacity_results = capacity_results,
      ssp              = ssp,
      ssp_colour       = ssp_col,
      y_limits         = temp_y_limits,
      show_y_axis      = show_y,
      show_x_axis      = TRUE
    )
    temp_panels[[i]] <- p_temp +
      labs(subtitle = bquote(bold(.(paste0(panel_letters_temp[i], ")"))) ~ .(ssp))) +
      theme(plot.subtitle = element_text(size = 9, hjust = 0, colour = "black"))
    
    p_cdr <- plot_capacity_cdr(
      capacity_results    = capacity_results,
      ssp                 = ssp,
      ssp_colour          = ssp_col,
      y_limits            = cdr_y_limits,
      show_y_axis         = show_y,
      show_x_axis         = FALSE,
      show_capacity_limit = show_capacity_limit
    )
    cdr_panels[[i]] <- p_cdr +
      labs(subtitle = bquote(bold(.(paste0(panel_letters_cdr[i], ")"))) ~ .(ssp))) +
      theme(plot.subtitle = element_text(size = 9, hjust = 0, colour = "black"))
  }
  
  # ── Build combined bottom legend ──────────────────────────────────────────
  # Two separate dummy plots, each yielding one legend grob via get_legend().
  # legend.position = "right" is used deliberately — cowplot::get_legend()
  # reliably finds the legend grob in this position. The grobs are then placed
  # side-by-side below the main grid regardless of their internal position.
  # guides(nrow = 1) forces each legend to lay its keys out horizontally.
  
  ssp_legend_df <- data.frame(
    x           = rep(1:2, times = length(ssp_names)),
    y           = rep(seq_along(ssp_names), each = 2),
    scenario    = factor(
      rep(ssp_names, each = 2),
      levels = ssp_names
    )
  )
  ssp_legend_plot <- ggplot(ssp_legend_df,
                            aes(x = x, y = y, colour = scenario, group = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(
      name   = "Scenario",
      values = setNames(scenario_colors, ssp_names)
    ) +
    guides(colour = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(
      legend.position   = "right",
      legend.direction  = "horizontal",
      legend.title      = element_text(size = 9, face = "bold"),
      legend.text       = element_text(size = 8),
      legend.key.width  = unit(1.5, "lines"),
      legend.key.height = unit(0.8, "lines"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  ssp_leg <- cowplot::get_legend(ssp_legend_plot)
  
  rate_legend_df <- data.frame(
    x           = rep(1:2, times = length(growth_rate_linetypes)),
    y           = rep(seq_along(growth_rate_linetypes), each = 2),
    growth_rate = factor(
      rep(names(growth_rate_linetypes), each = 2),
      levels = names(growth_rate_linetypes)
    )
  )
  rate_legend_plot <- ggplot(rate_legend_df,
                             aes(x = x, y = y,
                                 linetype = growth_rate, group = growth_rate)) +
    geom_line(colour = "grey30", linewidth = 0.8) +
    scale_linetype_manual(
      name   = "Growth rate",
      values = growth_rate_linetypes,
      labels = growth_rate_labels
    ) +
    guides(linetype = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(
      legend.position   = "right",
      legend.direction  = "horizontal",
      legend.title      = element_text(size = 9, face = "bold"),
      legend.text       = element_text(size = 8),
      legend.key.width  = unit(1.5, "lines"),
      legend.key.height = unit(0.8, "lines"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  rate_leg <- cowplot::get_legend(rate_legend_plot)
  
  # Place both legend grobs side by side in a narrow strip
  combined_legend <- cowplot::plot_grid(
    ssp_leg, rate_leg,
    nrow       = 1,
    rel_widths = c(1, 0.7)
  )
  
  # ── Row labels ────────────────────────────────────────────────────────────
  # Add a left-aligned row label as the subtitle of the first panel in each row.
  # These sit above the panel and identify the variable for the whole row.
  temp_panels[[1]] <- temp_panels[[1]] +
    labs(tag = "Temperature trajectories") +
    theme(
      plot.tag          = element_text(size = 9, face = "bold", hjust = 0),
      plot.tag.position = "top"
    )
  
  cdr_panels[[1]] <- cdr_panels[[1]] +
    labs(tag = "CDR deployment") +
    theme(
      plot.tag          = element_text(size = 9, face = "bold", hjust = 0),
      plot.tag.position = "top"
    )
  
  # ── Assemble grid with patchwork ──────────────────────────────────────────
  temp_row <- wrap_plots(temp_panels, nrow = 1)
  cdr_row  <- wrap_plots(cdr_panels,  nrow = 1)
  
  main_grid <- cdr_row / temp_row
  
  # Stack: title / main grid / legend strip
  final_plot <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        "Capacity Scenario Comparison: Optimal Control Results",
        fontface = "bold", size = 12, hjust = 0.5
      ),
    main_grid,
    combined_legend,
    ncol        = 1,
    rel_heights = c(0.04, 1, 0.1)
  )
  
  # ── Save if requested ─────────────────────────────────────────────────────
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0(
        "capacity_scenario_comparison_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    }
    
    filepath <- here::here("figs", filename)
    ggsave(filepath, final_plot,
           width = width, height = height, units = "mm",
           device = cairo_pdf)
    
    # PNG version for publication / sharing
    png_filename <- sub("\\.pdf$", ".png", filename)
    png_filepath <- here::here("figs", png_filename)
    ggsave(png_filepath, final_plot,
           width = width, height = height, units = "mm",
           device = "png", dpi = 300, bg = "white")
    
    if (verbose) {
      cat("Dashboard saved to:", filepath, "\n")
      cat("PNG version saved to:", png_filepath, "\n")
    }
  }
  
  return(final_plot)
}

# ============================================================================
# Usage example
# ============================================================================

# Run the multi-rate comparison from the scenario_comparison_analysis_capacity
# chunk in MACROM_workflow.Rmd, then pass the result directly here:
#
# fig <- create_capacity_scenario_comparison_dashboard(capacity_growth_results)
# print(fig)
#
# # Save to figs/ directory
# create_capacity_scenario_comparison_dashboard(
#   capacity_results = capacity_growth_results,
#   save_plot        = TRUE,
#   filename         = "capacity_scenario_comparison.pdf"
# )
#
# # Or reload a previously saved combined RDS and visualise without re-running:
# capacity_growth_results <- readRDS(here::here("output",
#   "capacity_growth_comparison_combined_20260324_120000.rds"))
# fig <- create_capacity_scenario_comparison_dashboard(capacity_growth_results)
# print(fig)