# ==============================================================================
# Lowest Scenario Comparison Dashboard
#
# Part of: MACROM: An Optimal Control Model for Balancing Climate Change Abatement
# and Damage Trade-offs
# Authors: Nina Rynne, Michael Bode, Melanie Roberts, Ryan Heneghan
# Institution: Griffith University
#
# Copyright (c) 2025 Nina Rynne
# Licensed under CC-BY-4.0 - see LICENSE file for details
#
# Version: 1.0.0
# Last updated: May 2026
# ==============================================================================

#' @title Lowest Scenario Comparison Dashboard
#' @description
#' Creates a figure comparing temperature trajectories and CDR strategies across
#' SSP scenarios with SSP-specific (r, K) parameter pairs. Layout mirrors
#' create_capacity_scenario_comparison_dashboard(): 2 rows (CDR top, temperature
#' bottom) × N columns (one per SSP in lowest_results). Within each panel,
#' lines are differentiated by pair (pair_1, pair_2, ...) using linetypes.
#'
#' Pair labels in the legend:
#'   pair_1 → "Lowest growth rate"
#'   pair_2 → "Lowest maximum capacity"
#'   pair_3+ → "Pair 3", "Pair 4", ...
#'
#' Expects lowest_results to be the object returned by run_lowest_comparison(),
#' i.e. a nested list keyed lowest_results[[ssp]][[pair_label]], where each
#' leaf is the full output of run_scenario_comparison() for that single SSP.
#'
#' @author Nina Rynne
#' @date May 2026

# ============================================================================
# Global definitions
# ============================================================================

# Reuse scenario colours from scenario_comparison_visualisation.R
scenario_colors <- c(
  "#00ADCF",  # Cyan      — SSP1
  "#173C66",  # Dark blue — SSP2
  "#F0E442",  # Yellow    — SSP3
  "#E71D25",  # Red       — SSP4
  "#951B1E"   # Dark red  — SSP5
)

# Map SSP full names to their colour index
ssp_color_index <- c(
  "SSP1" = 1, "SSP2" = 2, "SSP3" = 3, "SSP4" = 4, "SSP5" = 5
)

# Available linetypes — assigned to pairs in order
available_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

# Legend labels for pairs — pair_1 and pair_2 have semantic names; rest fall back
pair_legend_labels <- c(
  "pair_1" = "Lowest growth rate",
  "pair_2" = "Lowest maximum capacity"
)

# Base theme (matches existing visualisation files)
scenario_comparison_theme <- theme_bw() +
  theme(
    text             = element_text(size = 10),
    plot.title       = element_text(size = 10),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.position  = "none"
  )

# ============================================================================
# Helper functions
# ============================================================================

#' @title Clean Scenario Names for Display
clean_scenario_names <- function(scenario_names) {
  gsub("-Baseline$", "", scenario_names)
}

#' @title Resolve SSP Colour
#' @description Returns the hex colour for an SSP, matching on the short name
#'   ("SSP1") regardless of whether the full "SSP1-Baseline" form is supplied.
resolve_ssp_colour <- function(ssp) {
  short <- clean_scenario_names(ssp)
  idx   <- ssp_color_index[[short]]
  if (is.null(idx)) {
    warning(sprintf("No colour defined for SSP '%s'; using grey.", ssp))
    return("grey50")
  }
  scenario_colors[idx]
}

#' @title Build Pair Linetype and Label Vectors for a Given Set of Pair Labels
#' @description
#' Returns named vectors of linetypes and display labels for the supplied pair
#' labels. Pairs beyond "pair_2" receive generic labels ("Pair 3", ...).
#'
#' @param pair_labels Character vector of pair labels, e.g. c("pair_1", "pair_2").
#' @return List with elements `linetypes` and `labels`, both named by pair_label.
build_pair_aesthetics <- function(pair_labels) {
  n <- length(pair_labels)
  if (n > length(available_linetypes)) {
    stop("More pairs than available linetypes (max ", length(available_linetypes),
         "). Reduce the number of pairs or extend available_linetypes.")
  }

  linetypes <- setNames(available_linetypes[seq_len(n)], pair_labels)

  labels <- vapply(pair_labels, function(lbl) {
    if (!is.null(pair_legend_labels[[lbl]])) {
      pair_legend_labels[[lbl]]
    } else {
      # Extract trailing number from "pair_N" for fallback label
      num <- sub("^pair_", "", lbl)
      paste("Pair", num)
    }
  }, character(1))
  labels <- setNames(labels, pair_labels)

  list(linetypes = linetypes, labels = labels)
}

#' @title Extract Variable Data Across Pairs for One SSP
#' @description
#' Pulls a single numeric variable from the scenario_results component of each
#' pair for a given SSP, returning a tidy data frame ready for plotting.
#'
#' @param lowest_results Named list as returned by run_lowest_comparison().
#' @param ssp Character; SSP key as it appears in lowest_results (e.g.
#'   "SSP3-Baseline"). Both "SSP3" and "SSP3-Baseline" forms are accepted.
#' @param variable Character; field name within the solution object
#'   (e.g. "temperature_anomaly", "qty_remov").
#' @return Data frame with columns: pair, years, value.
extract_lowest_variable <- function(lowest_results, ssp, variable) {

  # Accept both short ("SSP3") and full ("SSP3-Baseline") key forms
  ssp_entry <- lowest_results[[ssp]]
  if (is.null(ssp_entry)) {
    ssp_full <- paste0(ssp, "-Baseline")
    ssp_entry <- lowest_results[[ssp_full]]
  }
  if (is.null(ssp_entry)) {
    warning(sprintf("No results found for SSP '%s' in lowest_results.", ssp))
    return(data.frame(pair = character(), years = numeric(), value = numeric()))
  }

  map_dfr(names(ssp_entry), function(pair) {

    pair_run         <- ssp_entry[[pair]]
    scenario_results <- pair_run$scenario_results

    if (is.null(scenario_results)) {
      warning(sprintf("No $scenario_results found for SSP '%s', pair '%s'.", ssp, pair))
      return(NULL)
    }

    # Match SSP key inside scenario_results with or without "-Baseline"
    ssp_key <- if (!is.null(scenario_results[[ssp]])) {
      ssp
    } else if (!is.null(scenario_results[[paste0(ssp, "-Baseline")]])) {
      paste0(ssp, "-Baseline")
    } else if (!is.null(scenario_results[[clean_scenario_names(ssp)]])) {
      clean_scenario_names(ssp)
    } else {
      warning(sprintf(
        "Cannot find scenario_results key for SSP '%s', pair '%s'.", ssp, pair))
      return(NULL)
    }

    result <- scenario_results[[ssp_key]]

    data.frame(
      pair  = pair,
      years = result$years,
      value = result[[variable]]
    )
  })
}

#' @title Compute CDR Capacity Curve Data Across Pairs for One SSP
#' @description
#' Computes the logistic CDR capacity curve for each pair of a given SSP
#' directly from ssp_params, g_initial, and t_start. The years vector is read
#' from the first scenario_results entry for each pair in lowest_results.
#' This mirrors the approach used by extract_cdr_capacity() in the old
#' scenario_comparison_capacity_visualisation.R, which reads the logistic
#' parameters from run_info and evaluates the formula directly.
#'
#' @param lowest_results Named list as returned by run_lowest_comparison().
#' @param ssp Character; SSP key in lowest_results (e.g. "SSP3-Baseline").
#' @param ssp_params Named list identical to that passed to run_lowest_comparison().
#' @param g_initial Numeric; shared CDR starting level (GtCO2/yr).
#' @param t_start Numeric; shared CDR deployment start year.
#' @return Data frame with columns: pair, years, value.
extract_lowest_cdr_capacity <- function(lowest_results, ssp,
                                        ssp_params,
                                        g_initial,
                                        t_start) {

  ssp_entry <- lowest_results[[ssp]]
  if (is.null(ssp_entry)) {
    warning(sprintf("No results found for SSP '%s' in lowest_results.", ssp))
    return(data.frame(pair = character(), years = numeric(), value = numeric()))
  }

  # Resolve the ssp_params entry — try the key exactly as given
  ssp_params_entry <- ssp_params[[ssp]]
  if (is.null(ssp_params_entry)) {
    warning(sprintf("ssp_params has no entry for SSP '%s'.", ssp))
    return(data.frame(pair = character(), years = numeric(), value = numeric()))
  }

  pair_labels <- names(ssp_entry)

  map_dfr(seq_along(pair_labels), function(i) {
    pair     <- pair_labels[i]
    pair_vec <- ssp_params_entry[[i]]   # c(r = ..., K = ...)

    if (is.null(pair_vec) || !all(c("r", "K") %in% names(pair_vec))) {
      warning(sprintf("ssp_params entry for SSP '%s', pair %d is missing r or K.", ssp, i))
      return(NULL)
    }

    r <- pair_vec[["r"]]
    K <- pair_vec[["K"]]

    # Years from the first scenario_results entry for this pair — same as how
    # extract_cdr_capacity() gets years in the old visualisation function
    years <- ssp_entry[[pair]]$scenario_results[[1]]$years

    if (is.null(years)) {
      warning(sprintf("Cannot read years for SSP '%s', pair '%s'.", ssp, pair))
      return(NULL)
    }

    suppression_factor <- (K / g_initial) - 1
    capacity_values <- ifelse(
      years < t_start,
      g_initial,
      K / (1 + suppression_factor * exp(-r * (years - t_start)))
    )

    data.frame(pair = pair, years = years, value = capacity_values)
  })
}

# ============================================================================
# Individual panel functions
# ============================================================================

#' @title Plot Temperature Trajectories for One SSP Across All Pairs
#' @param lowest_results Named list as returned by run_lowest_comparison().
#' @param ssp Character; SSP key in lowest_results.
#' @param ssp_colour Hex colour for this SSP's lines.
#' @param pair_aesthetics List from build_pair_aesthetics() with $linetypes and $labels.
#' @param y_limits Numeric vector of length 2 for shared y-axis limits.
#' @param show_y_axis Logical; show y-axis title and labels (default TRUE).
#' @param show_x_axis Logical; show x-axis title (default TRUE).
#' @return ggplot object.
plot_lowest_temperature <- function(lowest_results,
                                    ssp,
                                    ssp_colour,
                                    pair_aesthetics,
                                    y_limits    = c(NA, NA),
                                    show_y_axis = TRUE,
                                    show_x_axis = TRUE) {

  plot_data <- extract_lowest_variable(lowest_results, ssp, "temperature_anomaly")

  if (nrow(plot_data) == 0) {
    stop(sprintf("No valid temperature data found for SSP '%s'", ssp))
  }

  plot_data$pair <- factor(plot_data$pair, levels = names(pair_aesthetics$linetypes))

  p <- ggplot(plot_data, aes(x = years, y = value, linetype = pair)) +
    geom_hline(yintercept = 1.5, linetype = "dashed", colour = "red",
               alpha = 0.7, linewidth = 0.5) +
    geom_line(colour = ssp_colour, linewidth = 0.8) +
    scale_linetype_manual(
      name   = "Parameter combination",
      values = pair_aesthetics$linetypes,
      labels = pair_aesthetics$labels
    ) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      x = if (show_x_axis) "Year" else NULL,
      y = if (show_y_axis) "Temperature anomaly (°C)" else NULL
    ) +
    scenario_comparison_theme +
    theme(
      axis.title.y = if (!show_y_axis) element_blank() else element_text(size = 9),
      axis.text.y  = if (!show_y_axis) element_blank() else element_text(size = 8),
      axis.ticks.y = if (!show_y_axis) element_blank() else element_line()
    )

  return(p)
}

#' @title Plot CDR Deployment for One SSP Across All Pairs
#' @param lowest_results Named list as returned by run_lowest_comparison().
#' @param ssp Character; SSP key in lowest_results.
#' @param ssp_colour Hex colour for this SSP's lines.
#' @param pair_aesthetics List from build_pair_aesthetics() with $linetypes and $labels.
#' @param y_limits Numeric vector of length 2 for shared y-axis limits.
#' @param show_y_axis Logical; show y-axis title and labels (default TRUE).
#' @param show_x_axis Logical; show x-axis title (default TRUE).
#' @param show_capacity_limit Logical; overlay the logistic capacity curve for
#'   each pair as a grey line with matching linetype (default FALSE).
#' @param ssp_params Optional; passed to extract_lowest_cdr_capacity() as
#'   fallback when run_info fields are absent (older RDS files).
#' @param g_initial Optional numeric; shared CDR start level fallback.
#' @param t_start Optional numeric; shared CDR start year fallback.
#' @return ggplot object.
plot_lowest_cdr <- function(lowest_results,
                             ssp,
                             ssp_colour,
                             pair_aesthetics,
                             y_limits            = c(0, NA),
                             show_y_axis         = TRUE,
                             show_x_axis         = TRUE,
                             show_capacity_limit = FALSE,
                             ssp_params          = NULL,
                             g_initial           = NULL,
                             t_start             = NULL) {

  plot_data <- extract_lowest_variable(lowest_results, ssp, "qty_remov")

  if (nrow(plot_data) == 0) {
    stop(sprintf("No valid CDR data found for SSP '%s'", ssp))
  }

  plot_data$pair <- factor(plot_data$pair, levels = names(pair_aesthetics$linetypes))

  p <- ggplot(plot_data, aes(x = years, y = value, linetype = pair))

  # Add capacity curves behind CDR lines when requested
  if (show_capacity_limit) {
    capacity_data <- extract_lowest_cdr_capacity(
      lowest_results = lowest_results,
      ssp            = ssp,
      ssp_params     = ssp_params,
      g_initial      = g_initial,
      t_start        = t_start
    )
    if (nrow(capacity_data) > 0) {
      capacity_data$pair <- factor(capacity_data$pair,
                                   levels = names(pair_aesthetics$linetypes))
      p <- p +
        geom_line(data        = capacity_data,
                  mapping     = aes(x = years, y = value, linetype = pair),
                  colour      = "grey70",
                  linewidth   = 0.8,
                  inherit.aes = FALSE)
    }
  }

  p <- p +
    geom_line(colour = ssp_colour, linewidth = 0.8) +
    scale_linetype_manual(
      name   = "Parameter combination",
      values = pair_aesthetics$linetypes,
      labels = pair_aesthetics$labels
    ) +
    scale_y_continuous(limits = y_limits) +
    labs(
      x = if (show_x_axis) "Year" else NULL,
      y = if (show_y_axis) expression("CDR (GtCO"[2]*"/yr)") else NULL
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

#' @title Create Lowest Scenario Comparison Dashboard
#' @description
#' Creates a 2-row × N-column figure (CDR top, temperature bottom) comparing
#' results across all SSPs and parameter pairs in lowest_results. Each column
#' corresponds to one SSP; within each panel, the pairs are shown as different
#' linetypes. The number of columns adapts automatically to the number of SSPs
#' present in lowest_results.
#'
#' @param lowest_results Named list as returned by run_lowest_comparison().
#' @param show_capacity_limit Logical; overlay the logistic CDR capacity curve
#'   for each pair as a grey line with matching linetype on CDR panels
#'   (default FALSE).
#' @param ssp_params Optional named list; the same ssp_params object passed to
#'   run_lowest_comparison(). Required when loading results from an older RDS
#'   where run_info does not contain capacity parameters, and
#'   show_capacity_limit = TRUE.
#' @param g_initial Optional numeric; shared CDR starting level (GtCO2/yr).
#'   Required alongside ssp_params when using the fallback path.
#' @param t_start Optional numeric; shared CDR start year.
#'   Required alongside ssp_params when using the fallback path.
#' @param save_plot Logical; save the figure to file (default FALSE).
#' @param filename Character; custom filename (default: auto-generated timestamp).
#'   Extension determines format — use ".pdf" or ".png".
#' @param width Plot width in mm (default 297, A4 landscape).
#' @param height Plot height in mm (default 180).
#' @param verbose Logical; print progress messages (default TRUE).
#' @return Combined ggplot object containing the dashboard.
#'
#' @examples
#' fig <- create_lowest_comparison_dashboard(lowest_results)
#' print(fig)
#'
#' # With capacity overlay — pass ssp_params as fallback for older RDS files
#' create_lowest_comparison_dashboard(
#'   lowest_results      = lowest_results,
#'   show_capacity_limit = TRUE,
#'   ssp_params          = ssp_params,
#'   g_initial           = cdr_g_initial,
#'   t_start             = cdr_t_start,
#'   save_plot           = TRUE,
#'   filename            = "lowest_comparison.pdf"
#' )
create_lowest_comparison_dashboard <- function(lowest_results,
                                               show_capacity_limit = FALSE,
                                               ssp_params          = NULL,
                                               g_initial           = NULL,
                                               t_start             = NULL,
                                               save_plot           = FALSE,
                                               filename            = NULL,
                                               width               = 297,
                                               height              = 180,
                                               verbose             = TRUE) {

  # ── Derive SSP list and pair aesthetics ──────────────────────────────────────

  ssp_list <- names(lowest_results)
  if (length(ssp_list) == 0) stop("lowest_results is empty.")

  # Union of all pair labels across all SSPs (preserving pair_1, pair_2, ... order)
  all_pairs <- unique(unlist(lapply(lowest_results, names)))
  # Sort so pair_1 < pair_2 < ... regardless of SSP order
  all_pairs <- all_pairs[order(as.integer(sub("^pair_", "", all_pairs)))]

  pair_aes <- build_pair_aesthetics(all_pairs)

  if (verbose) {
    cat("Creating lowest scenario comparison dashboard\n")
    cat("SSPs:  ", paste(clean_scenario_names(ssp_list), collapse = ", "), "\n")
    cat("Pairs: ", paste(all_pairs, "=", pair_aes$labels[all_pairs]), "\n")
  }

  # ── Compute shared y-axis limits across all runs ──────────────────────────

  all_temp <- map_dfr(ssp_list, function(ssp) {
    map_dfr(names(lowest_results[[ssp]]), function(pair) {
      sr <- lowest_results[[ssp]][[pair]]$scenario_results
      result <- sr[[ssp]] %||% sr[[paste0(ssp, "-Baseline")]] %||%
                sr[[clean_scenario_names(ssp)]]
      if (!is.null(result)) data.frame(temp = result$temperature_anomaly)
    })
  })
  temp_y_limits <- c(
    floor(min(all_temp$temp,   na.rm = TRUE) * 10) / 10,
    ceiling(max(all_temp$temp, na.rm = TRUE) * 10) / 10 + 0.1
  )

  all_cdr <- map_dfr(ssp_list, function(ssp) {
    map_dfr(names(lowest_results[[ssp]]), function(pair) {
      sr <- lowest_results[[ssp]][[pair]]$scenario_results
      result <- sr[[ssp]] %||% sr[[paste0(ssp, "-Baseline")]] %||%
                sr[[clean_scenario_names(ssp)]]
      if (!is.null(result)) data.frame(cdr = result$qty_remov)
    })
  })
  # If capacity lines are shown, include their values in the upper bound so the
  # capacity curve is never clipped — mirrors the existing dashboard behaviour
  if (show_capacity_limit) {
    cap_data <- map_dfr(ssp_list, function(ssp) {
      extract_lowest_cdr_capacity(lowest_results, ssp,
                                  ssp_params = ssp_params,
                                  g_initial  = g_initial,
                                  t_start    = t_start)
    })
    cdr_upper <- if (nrow(cap_data) > 0) {
      max(max(all_cdr$cdr, na.rm = TRUE), max(cap_data$value, na.rm = TRUE))
    } else {
      max(all_cdr$cdr, na.rm = TRUE)
    }
  } else {
    cdr_upper <- max(all_cdr$cdr, na.rm = TRUE)
  }
  cdr_y_limits <- c(0, ceiling(cdr_upper * 1.05))

  if (verbose) {
    cat(sprintf("Temperature y-axis: %.1f – %.1f °C\n",
                temp_y_limits[1], temp_y_limits[2]))
    cat(sprintf("CDR y-axis: 0 – %.0f GtCO2/yr\n", cdr_y_limits[2]))
  }

  # ── Build individual panels ───────────────────────────────────────────────

  n_ssp            <- length(ssp_list)
  panel_letters_cdr  <- letters[seq_len(n_ssp)]
  panel_letters_temp <- letters[seq(n_ssp + 1, 2 * n_ssp)]

  temp_panels <- vector("list", n_ssp)
  cdr_panels  <- vector("list", n_ssp)

  for (i in seq_along(ssp_list)) {
    ssp     <- ssp_list[i]
    ssp_col <- resolve_ssp_colour(ssp)
    show_y  <- (i == 1)
    short   <- clean_scenario_names(ssp)

    # Subset pair aesthetics to only the pairs present for this SSP
    ssp_pairs     <- names(lowest_results[[ssp]])
    ssp_pair_aes  <- list(
      linetypes = pair_aes$linetypes[ssp_pairs],
      labels    = pair_aes$labels[ssp_pairs]
    )

    p_temp <- plot_lowest_temperature(
      lowest_results  = lowest_results,
      ssp             = ssp,
      ssp_colour      = ssp_col,
      pair_aesthetics = ssp_pair_aes,
      y_limits        = temp_y_limits,
      show_y_axis     = show_y,
      show_x_axis     = TRUE
    )
    temp_panels[[i]] <- p_temp +
      labs(subtitle = bquote(bold(.(paste0(panel_letters_temp[i], ")"))) ~ .(short))) +
      theme(plot.subtitle = element_text(size = 9, hjust = 0, colour = "black"))

    p_cdr <- plot_lowest_cdr(
      lowest_results      = lowest_results,
      ssp                 = ssp,
      ssp_colour          = ssp_col,
      pair_aesthetics     = ssp_pair_aes,
      y_limits            = cdr_y_limits,
      show_y_axis         = show_y,
      show_x_axis         = FALSE,
      show_capacity_limit = show_capacity_limit,
      ssp_params          = ssp_params,
      g_initial           = g_initial,
      t_start             = t_start
    )
    cdr_panels[[i]] <- p_cdr +
      labs(subtitle = bquote(bold(.(paste0(panel_letters_cdr[i], ")"))) ~ .(short))) +
      theme(plot.subtitle = element_text(size = 9, hjust = 0, colour = "black"))
  }

  # ── Row labels ────────────────────────────────────────────────────────────
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

  # ── Build combined legend ─────────────────────────────────────────────────

  # SSP colour legend
  short_ssp_names <- clean_scenario_names(ssp_list)
  ssp_legend_df <- data.frame(
    x        = rep(1:2, times = n_ssp),
    y        = rep(seq_len(n_ssp), each = 2),
    scenario = factor(rep(short_ssp_names, each = 2), levels = short_ssp_names)
  )
  ssp_colours_named <- setNames(
    vapply(ssp_list, resolve_ssp_colour, character(1)),
    short_ssp_names
  )
  ssp_legend_plot <- ggplot(ssp_legend_df,
                             aes(x = x, y = y, colour = scenario, group = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(name = "Scenario", values = ssp_colours_named) +
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

  # Pair linetype legend (uses global all_pairs for consistency)
  pair_legend_df <- data.frame(
    x    = rep(1:2, times = length(all_pairs)),
    y    = rep(seq_along(all_pairs), each = 2),
    pair = factor(rep(all_pairs, each = 2), levels = all_pairs)
  )
  pair_legend_plot <- ggplot(pair_legend_df,
                              aes(x = x, y = y, linetype = pair, group = pair)) +
    geom_line(colour = "grey30", linewidth = 0.8) +
    scale_linetype_manual(
      name   = "Parameter combination",
      values = pair_aes$linetypes,
      labels = pair_aes$labels
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
  pair_leg <- cowplot::get_legend(pair_legend_plot)

  combined_legend <- cowplot::plot_grid(
    ssp_leg, pair_leg,
    nrow       = 1,
    rel_widths = c(1, 0.9)
  )

  # ── Assemble grid with patchwork ──────────────────────────────────────────
  temp_row  <- wrap_plots(temp_panels, nrow = 1)
  cdr_row   <- wrap_plots(cdr_panels,  nrow = 1)
  main_grid <- cdr_row / temp_row

  final_plot <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        "Lowest Scenario Comparison: Optimal Control Results",
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
        "lowest_scenario_comparison_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    }

    filepath <- here::here("figs", filename)
    ggsave(filepath, final_plot,
           width = width, height = height, units = "mm",
           device = cairo_pdf)

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
