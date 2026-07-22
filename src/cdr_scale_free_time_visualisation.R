# ==============================================================================
# CDR Scale Free-Terminal-Time Visualisation Functions
#
# Part of: MACROM: An Optimal Control Model for Balancing Climate Change Abatement
# and Damage Trade-offs
# Authors: Nina Rynne, Michael Bode, Melanie Roberts, Ryan Heneghan
# Institution: Griffith University
#
# Copyright (c) 2026 Nina Rynne
# Licensed under CC-BY-4.0 - see LICENSE file for details
#
# Citation: If you use this code, please cite:
#   Rynne, N., Bode, M., Roberts, M., & Heneghan, R. (2025).
#   MACROM: An Optimal Control Model for Balancing Climate Change
#   Abatement and Damage Trade-offs.
#   https://doi.org/10.5281/zenodo.18463951
#
# Description:
# Heatmap visualisation for run_free_time_scale_sensitivity() results
# (cdr_scale_free_time_sensitivity.R). Structurally parallel to
# create_cdr_scale_outcome_plot() in cdr_scale_sensitivity_visualisation.R —
# same panel-per-SSP layout, K on the x-axis, r on the y-axis — but instead of
# classifying cells into no-overshoot/recoverable/unrecoverable categories,
# each tile is filled with the continuous return year T_star on a viridis
# scale FIXED to [2100, 2200] (the free-time search horizon), so panels for
# different scenarios or grids remain visually comparable. Cells with no
# return achieved within the horizon (feasible == FALSE) are filled grey and
# marked with a red "x", reusing add_cdr_scale_infeasible_markers().
#
# DEPENDENCIES: source cdr_scale_sensitivity_visualisation.R first — this file
# reuses its shared globals (BASE_CDR_SCALE_THEME, SSP_SCENARIO_ORDER_SCALE)
# and helper (add_cdr_scale_infeasible_markers()) rather than redefining them.
#
# Version: 1.0.0
# Last updated: July 2026
# ==============================================================================


#' @title Create Free-Terminal-Time Return Year Heatmap
#' @description
#' Heatmap of the free-terminal-time return year (T_star) across the K/r
#' CDR-scale grid, one panel per SSP scenario. Colour scale is fixed to
#' year_limits (default c(2100, 2200), matching the free-time search horizon
#' used by run_free_time_scale_sensitivity()) rather than data-driven, so the
#' same shade always means the same calendar year across panels, scenarios,
#' and separate figures. Cells that never achieve a return within the search
#' horizon (feasible == FALSE, T_star == NA) are filled grey and marked with
#' a red "x".
#'
#' @param sensitivity_results Results object from
#'   run_free_time_scale_sensitivity() (or a bare combined_results data frame)
#' @param year_limits Numeric c(min, max) for the viridis colour scale
#'   (default: c(2100, 2200))
#' @param show_infeasible Logical: mark cells with no return within the
#'   search horizon (default: TRUE)
#' @param save_plot Logical: save to file (default: FALSE)
#' @param filename Character string for output filename (default: NULL for auto)
#' @param width Numeric width in mm (default: 297)
#' @param height Numeric height in mm (default: 100)
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return Patchwork multi-panel plot object (one panel per SSP scenario present)
#'
#' @examples
#' \dontrun{
#' p <- create_free_time_year_heatmap(
#'   sensitivity_results = free_time_scale_results,
#'   save_plot           = TRUE
#' )
#' print(p)
#' }
create_free_time_year_heatmap <- function(sensitivity_results,
                                          year_limits     = c(2100, 2200),
                                          show_infeasible = TRUE,
                                          save_plot       = FALSE,
                                          filename        = NULL,
                                          width           = 297,
                                          height          = 100,
                                          verbose         = TRUE) {

  if (verbose) cat("Preparing data for free-terminal-time return year heatmap...\n")

  # Accept either a full results list or a bare data frame
  combined_data <- if (is.data.frame(sensitivity_results)) {
    sensitivity_results
  } else if (is.list(sensitivity_results) &&
             "combined_results" %in% names(sensitivity_results)) {
    sensitivity_results$combined_results
  } else {
    stop("sensitivity_results must be either a data frame or a list containing ",
         "'combined_results' (output from run_free_time_scale_sensitivity())")
  }

  if (is.null(combined_data) || nrow(combined_data) == 0) {
    stop("No valid data found in sensitivity_results")
  }

  required_cols <- c("K", "r", "scenario", "scenario_short", "T_star", "feasible")
  missing_cols  <- setdiff(required_cols, names(combined_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Filter to valid SSP scenarios. Unlike prepare_cdr_scale_data(), rows with
  # T_star == NA are deliberately KEPT (not dropped) so infeasible cells still
  # render — as grey tiles with a red "x" marker rather than being blank.
  plot_data <- combined_data %>%
    dplyr::filter(scenario_short %in% SSP_SCENARIO_ORDER_SCALE) %>%
    dplyr::arrange(scenario_short, K, r)

  if (nrow(plot_data) == 0) {
    stop("No data found matching SSP scenarios: ",
         paste(SSP_SCENARIO_ORDER_SCALE, collapse = ", "))
  }

  # Warn (don't error) if any feasible T_star falls outside year_limits —
  # squish() will clamp it visually, which could mask a real out-of-range value
  out_of_range <- plot_data %>%
    dplyr::filter(feasible, (T_star < year_limits[1] | T_star > year_limits[2]))
  if (nrow(out_of_range) > 0 && verbose) {
    cat(sprintf("WARNING: %d feasible T_star value(s) fall outside year_limits [%d, %d] and will be clamped visually.\n",
                nrow(out_of_range), year_limits[1], year_limits[2]))
  }

  if (verbose) {
    n_infeasible <- sum(!plot_data$feasible, na.rm = TRUE)
    cat(sprintf("Prepared free-time data: %d rows, %d infeasible (no return within search horizon)\n",
                nrow(plot_data), n_infeasible))
  }

  scenarios_present <- intersect(SSP_SCENARIO_ORDER_SCALE,
                                 unique(plot_data$scenario_short))
  panel_labels <- LETTERS

  shared_theme <- BASE_CDR_SCALE_THEME +
    theme(
      text            = element_text(size = 9),
      plot.title      = element_text(size = 9, hjust = 0.5),
      axis.title      = element_text(size = 8),
      axis.text       = element_text(size = 7),
      legend.position = "none",
      plot.margin     = margin(1, 3, 1, 3)
    )

  panels <- list()

  for (s_idx in seq_along(scenarios_present)) {

    scen        <- scenarios_present[s_idx]
    scen_data   <- plot_data %>% dplyr::filter(scenario_short == scen)
    panel_label <- panel_labels[s_idx]
    show_y      <- (s_idx == 1)

    p <- ggplot(scen_data, aes(x = K, y = r, fill = T_star)) +
      geom_tile(color = NA) +
      scale_fill_viridis_c(
        name      = "Return\nYear",
        option    = "viridis",
        direction = -1,     # later (worse) years render darker/more intense,
                            # matching the direction convention used for the
                            # other "higher is worse" variables in this module
        limits    = year_limits,
        oob       = scales::squish,
        na.value  = "grey85",
        labels    = scales::label_number(big.mark = "")
      ) +
      scale_y_continuous(breaks = c(0.02, 0.04, 0.08, 0.12, 0.16, 0.20),
                         labels = function(x) x * 100) +
      labs(
        title = scen,
        x     = expression("Scale (maximum capacity, GtCO"[2]*"/year)"),
        y     = if (show_y) "Speed (intrinsic growth rate, %)" else NULL,
        tag   = paste0(panel_label, ")")
      ) +
      shared_theme +
      theme(
        plot.tag     = element_text(size = 9, face = "bold", hjust = 0),
        axis.title.y = if (show_y) element_text(size = 8) else element_blank(),
        axis.text.y  = if (show_y) element_text(size = 7) else element_blank(),
        axis.ticks.y = if (show_y) element_line() else element_blank()
      )

    if (show_infeasible) {
      p <- add_cdr_scale_infeasible_markers(p, scen_data)
    }

    panels[[scen]] <- p

    if (verbose) cat(sprintf("  Panel %s (%s) complete\n", panel_label, scen))
  }

  # Assemble panels in a single row with one shared, collected legend at the
  # bottom (same pattern as create_cdr_scale_outcome_plot())
  panel_row <- patchwork::wrap_plots(panels, nrow = 1) &
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.4, "cm"),
      legend.key.width = unit(0.8, "cm")
    )
  combined <- panel_row +
    patchwork::plot_layout(guides = "collect")

  if (verbose) cat("Free-terminal-time return year heatmap assembled\n")

  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0("cdr_scale_free_time_return_year_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    }
    if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".pdf")
    }
    filepath        <- here::here("figs", filename)
    output_dir_full <- here::here("figs")
    if (!dir.exists(output_dir_full)) dir.create(output_dir_full, recursive = TRUE)
    ggsave(filename = filepath, plot = combined,
           width = width, height = height, units = "mm",
           device = cairo_pdf, dpi = 300)
    if (verbose) cat(sprintf("Saved to: %s\n", filepath))
    return(invisible(combined))
  }

  return(combined)
}
