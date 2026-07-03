# ==============================================================================
# CDR Delay x Growth Rate Visualisation Functions
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
# Last updated: June 2026
# ==============================================================================

#' @note All required libraries (ggplot2, dplyr, patchwork, purrr, cowplot,
#' viridis, here) must be loaded before using these functions.

# ==============================================================================
# Global definitions
# ==============================================================================

# SSP identifiers and display order
ssp_names_dg <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")

# Base plot theme
delay_growth_theme <- theme_bw() +
  theme(
    text             = element_text(size = 10),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.position  = "none"
  )

# ==============================================================================
# Helper functions
# ==============================================================================

#' @title Get Delay Values From a Rate Sub-list
#' @description
#' Recovers the numeric delay values from the names of a rate-level sub-list
#' returned by run_cdr_delay_growth_analysis(). Skips the "run_info" entry.
#'
#' @param rate_results Named list — one rate entry from cdr_delay_growth_results.
#' @return Numeric vector of delay values in ascending order.
get_rate_delays <- function(rate_results) {
  keys   <- setdiff(names(rate_results), "run_info")
  delays <- as.numeric(sub("^delay_", "", keys))
  sort(delays)
}


#' @title Extract Variable Across All Delays for One (Rate, SSP) Pair
#' @description
#' Pulls a single numeric variable from every delay entry for a given growth
#' rate and SSP, returning a tidy data frame ready for plotting.
#'
#' @param delay_growth_results Output from run_cdr_delay_growth_analysis().
#' @param rate Character; growth rate label (e.g. "slow").
#' @param ssp Character; SSP name. Both "SSP1" and "SSP1-Baseline" accepted.
#' @param variable Character; field name within each solution object
#'   (e.g. "temperature_anomaly", "qty_remov").
#' @return Data frame with columns: delay (numeric), years (numeric), value.
extract_delay_variable <- function(delay_growth_results, rate, ssp, variable) {

  rate_results <- delay_growth_results[[rate]]
  if (is.null(rate_results)) {
    warning(sprintf("No results found for growth rate '%s'", rate))
    return(data.frame(delay = numeric(), years = numeric(), value = numeric()))
  }

  delays    <- get_rate_delays(rate_results)
  plot_rows <- vector("list", length(delays))

  for (k in seq_along(delays)) {
    delay     <- delays[k]
    delay_key <- paste0("delay_", delay)

    sc_results <- rate_results[[delay_key]]$scenario_results
    if (is.null(sc_results)) next

    # Accept SSP name with or without "-Baseline" suffix
    ssp_key <- if (!is.null(sc_results[[ssp]])) {
      ssp
    } else if (!is.null(sc_results[[paste0(ssp, "-Baseline")]])) {
      paste0(ssp, "-Baseline")
    } else {
      warning(sprintf("No results for SSP '%s', rate '%s', delay %d", ssp, rate, delay))
      next
    }

    result <- sc_results[[ssp_key]]
    vals   <- result[[variable]]
    if (is.null(vals)) next

    plot_rows[[k]] <- data.frame(
      delay = delay,
      years = result$years,
      value = vals
    )
  }

  do.call(rbind, plot_rows[!sapply(plot_rows, is.null)])
}


#' @title Extract CDR Capacity Curve for One Growth Rate
#' @description
#' Evaluates the logistic capacity curve using parameters stored in the
#' rate-level run_info, returning a data frame in the same shape as
#' extract_delay_variable() output. The years axis is taken from the
#' delay_0 entry so the curve spans the full time horizon.
#'
#' @param delay_growth_results Output from run_cdr_delay_growth_analysis().
#' @param rate Character; growth rate label (e.g. "slow").
#' @return Data frame with columns: years (numeric), value (numeric).
extract_cdr_capacity_curve <- function(delay_growth_results, rate) {

  rate_results <- delay_growth_results[[rate]]
  run_info     <- rate_results[["run_info"]]

  if (is.null(run_info)) {
    warning(sprintf("No run_info found for rate '%s'", rate))
    return(NULL)
  }

  # Get the years axis from the delay_0 entry (first available SSP)
  first_sc <- rate_results[["delay_0"]]$scenario_results
  if (is.null(first_sc)) {
    warning(sprintf("No delay_0 results found for rate '%s'", rate))
    return(NULL)
  }
  years <- first_sc[[1]]$years

  g0                 <- run_info$g_initial
  K                  <- run_info$K
  r                  <- run_info$r_value
  t_start            <- run_info$t_start
  suppression_factor <- (K / g0) - 1

  capacity_values <- ifelse(
    years < t_start,
    g0,
    K / (1 + suppression_factor * exp(-r * (years - t_start)))
  )

  data.frame(years = years, value = capacity_values)
}

# ==============================================================================
# Individual panel functions
# ==============================================================================

#' @title Plot Temperature Trajectories for One (Rate, SSP) Panel
#' @description
#' Creates a single panel showing temperature anomaly over time for one SSP
#' and one growth rate, with lines coloured by CDR delay. A 1.5 degree
#' reference line is included.
#'
#' @param delay_growth_results Output from run_cdr_delay_growth_analysis().
#' @param rate Character; growth rate label.
#' @param ssp Character; SSP name (e.g. "SSP1").
#' @param y_limits Numeric vector of length 2 for shared y-axis limits.
#' @param delay_range Numeric vector of length 2 giving the full delay range
#'   for the colour scale (typically c(0, max_delay_years)).
#' @param line_alpha Line transparency (default 0.7).
#' @param show_y_axis Logical; show y-axis title and labels (default TRUE).
#' @param show_x_axis Logical; show x-axis title (default TRUE).
#' @return ggplot object.
plot_delay_temperature <- function(delay_growth_results,
                                   rate,
                                   ssp,
                                   y_limits    = c(NA, NA),
                                   delay_range = c(0, 70),
                                   line_alpha  = 0.7,
                                   show_y_axis = TRUE,
                                   show_x_axis = TRUE) {

  plot_data <- extract_delay_variable(
    delay_growth_results, rate, ssp, "temperature_anomaly"
  )

  if (nrow(plot_data) == 0) {
    warning(sprintf("No temperature data for rate '%s', SSP '%s'", rate, ssp))
    return(ggplot() + theme_void())
  }

  ggplot(plot_data, aes(x = years, y = value, group = delay, colour = delay)) +
    geom_hline(yintercept = 1.5, linetype = "dashed",
               colour = "red", alpha = 0.7, linewidth = 0.5) +
    geom_line(alpha = line_alpha, linewidth = 0.5) +
    scale_colour_viridis_c(
      option    = "plasma",
      direction = 1,
      limits    = delay_range,
      name      = "CDR delay\n(years)"
    ) +
    scale_y_continuous(limits = y_limits,
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      x = if (show_x_axis) "Year" else NULL,
      y = if (show_y_axis) "Temperature anomaly (°C)" else NULL
    ) +
    delay_growth_theme +
    theme(
      axis.title.y = if (!show_y_axis) element_blank() else element_text(size = 9),
      axis.text.y  = if (!show_y_axis) element_blank() else element_text(size = 8),
      axis.ticks.y = if (!show_y_axis) element_blank() else element_line()
    )
}


#' @title Plot CDR Deployment for One (Rate, SSP) Panel
#' @description
#' Creates a single panel showing annual CDR over time for one SSP and one
#' growth rate, with lines coloured by CDR delay. Optionally overlays the
#' logistic CDR capacity ceiling as a grey line.
#'
#' @param delay_growth_results Output from run_cdr_delay_growth_analysis().
#' @param rate Character; growth rate label.
#' @param ssp Character; SSP name (e.g. "SSP1").
#' @param y_limits Numeric vector of length 2 for shared y-axis limits.
#' @param delay_range Numeric vector of length 2 for the colour scale limits.
#' @param line_alpha Line transparency (default 0.7).
#' @param show_y_axis Logical; show y-axis title and labels (default TRUE).
#' @param show_x_axis Logical; show x-axis title (default TRUE).
#' @param show_capacity_limit Logical; overlay the CDR capacity ceiling
#'   as a grey line (default TRUE).
#' @return ggplot object.
plot_delay_cdr <- function(delay_growth_results,
                           rate,
                           ssp,
                           y_limits            = c(0, NA),
                           delay_range         = c(0, 70),
                           line_alpha          = 0.7,
                           show_y_axis         = TRUE,
                           show_x_axis         = TRUE,
                           show_capacity_limit = TRUE) {

  plot_data <- extract_delay_variable(
    delay_growth_results, rate, ssp, "qty_remov"
  )

  if (nrow(plot_data) == 0) {
    warning(sprintf("No CDR data for rate '%s', SSP '%s'", rate, ssp))
    return(ggplot() + theme_void())
  }

  p <- ggplot(plot_data, aes(x = years, y = value, group = delay, colour = delay))

  if (show_capacity_limit) {
    cap_data <- extract_cdr_capacity_curve(delay_growth_results, rate)
    if (!is.null(cap_data)) {
      p <- p +
        geom_line(data        = cap_data,
                  mapping     = aes(x = years, y = value),
                  colour      = "grey60",
                  linewidth   = 0.8,
                  linetype    = "solid",
                  inherit.aes = FALSE)
    }
  }

  p +
    geom_line(alpha = line_alpha, linewidth = 0.5) +
    scale_colour_viridis_c(
      option    = "plasma",
      direction = 1,
      limits    = delay_range,
      name      = "CDR delay\n(years)"
    ) +
    scale_y_continuous(limits = y_limits) +
    labs(
      x = if (show_x_axis) "Year" else NULL,
      y = if (show_y_axis) expression("CDR (GtCO"[2]*"/yr)") else NULL
    ) +
    delay_growth_theme +
    theme(
      axis.title.y = if (!show_y_axis) element_blank() else element_text(size = 9),
      axis.text.y  = if (!show_y_axis) element_blank() else element_text(size = 8),
      axis.ticks.y = if (!show_y_axis) element_blank() else element_line()
    )
}

# ==============================================================================
# Dashboard function
# ==============================================================================

#' @title Create CDR Delay x Growth Rate Dashboard
#' @description
#' Creates a figure comparing CDR deployment and temperature trajectories
#' across SSP scenarios for one or more CDR growth rates. Lines within each
#' panel are coloured by CDR delay using a plasma colour scale.
#'
#' Layout: one pair of rows (CDR top, temperature bottom) per growth rate,
#' with 5 SSP columns. A single shared colour bar legend is placed at the
#' bottom. If a single rate is requested, the figure is 2 rows x 5 columns;
#' for all three rates it is 6 rows x 5 columns.
#'
#' @param delay_growth_results Output from run_cdr_delay_growth_analysis().
#' @param rates Character vector of growth rate labels to include. If NULL
#'   (default), all non-run_info entries in delay_growth_results are used.
#' @param show_capacity_limit Logical; overlay the CDR capacity ceiling on
#'   CDR panels (default TRUE).
#' @param line_alpha Line transparency (default 0.7). Lower values help when
#'   many delay lines overlap at fine step sizes.
#' @param save_plot Logical; save to figs/ directory (default FALSE).
#' @param filename Character; output filename. Extension sets format (.pdf or
#'   .png). Auto-generated timestamp used if NULL.
#' @param width Plot width in mm (default 297, A4 landscape).
#' @param height Plot height in mm. Defaults to 140 mm per growth rate
#'   included (280 for two rates, 420 for all three).
#' @param verbose Logical; print progress messages (default TRUE).
#' @return Combined ggplot / cowplot object.
#'
#' @examples
#' # All three growth rates combined
#' fig <- create_cdr_delay_growth_dashboard(cdr_delay_growth_results)
#' print(fig)
#'
#' # Single growth rate only
#' fig_slow <- create_cdr_delay_growth_dashboard(
#'   cdr_delay_growth_results,
#'   rates = "slow"
#' )
#' print(fig_slow)
create_cdr_delay_growth_dashboard <- function(delay_growth_results,
                                               rates               = NULL,
                                               show_capacity_limit = TRUE,
                                               line_alpha          = 0.7,
                                               save_plot           = FALSE,
                                               filename            = NULL,
                                               width               = 297,
                                               height              = NULL,
                                               verbose             = TRUE) {

  # --------------------------------------------------------------------------
  # Determine which rates to plot
  # --------------------------------------------------------------------------
  available_rates <- setdiff(names(delay_growth_results), "run_info")

  if (is.null(rates)) {
    rates <- available_rates
  } else {
    missing_rates <- setdiff(rates, available_rates)
    if (length(missing_rates) > 0) {
      stop(sprintf(
        "Rates not found in delay_growth_results: %s\nAvailable: %s",
        paste(missing_rates, collapse = ", "),
        paste(available_rates, collapse = ", ")
      ))
    }
  }

  if (is.null(height)) {
    height <- 140 * length(rates)
  }

  if (verbose) {
    cat("Creating CDR delay x growth rate dashboard\n")
    cat("Growth rates:", paste(rates, collapse = ", "), "\n")
  }

  # --------------------------------------------------------------------------
  # Compute shared y-axis limits across all panels
  # --------------------------------------------------------------------------
  all_temp <- do.call(rbind, lapply(rates, function(rate) {
    do.call(rbind, lapply(ssp_names_dg, function(ssp) {
      extract_delay_variable(delay_growth_results, rate, ssp, "temperature_anomaly")
    }))
  }))

  temp_y_limits <- c(
    floor(min(all_temp$value,   na.rm = TRUE) * 10) / 10,
    ceiling(max(all_temp$value, na.rm = TRUE) * 10) / 10 + 0.1
  )

  all_cdr <- do.call(rbind, lapply(rates, function(rate) {
    do.call(rbind, lapply(ssp_names_dg, function(ssp) {
      extract_delay_variable(delay_growth_results, rate, ssp, "qty_remov")
    }))
  }))

  cdr_upper <- max(all_cdr$value, na.rm = TRUE)

  if (show_capacity_limit) {
    cap_vals <- do.call(rbind, lapply(rates, function(rate) {
      extract_cdr_capacity_curve(delay_growth_results, rate)
    }))
    if (!is.null(cap_vals)) {
      cdr_upper <- max(cdr_upper, max(cap_vals$value, na.rm = TRUE))
    }
  }

  cdr_y_limits <- c(0, ceiling(cdr_upper * 1.05))

  # Delay range for colour scale
  run_info    <- delay_growth_results[["run_info"]]
  delay_range <- if (!is.null(run_info$delays)) {
    c(min(run_info$delays), max(run_info$delays))
  } else {
    delays_all <- get_rate_delays(delay_growth_results[[rates[1]]])
    c(min(delays_all), max(delays_all))
  }

  if (verbose) {
    cat(sprintf("Temperature y-axis: %.1f - %.1f C\n",
                temp_y_limits[1], temp_y_limits[2]))
    cat(sprintf("CDR y-axis: 0 - %.0f GtCO2/yr\n", cdr_y_limits[2]))
    cat(sprintf("Delay colour range: %d - %d years\n",
                delay_range[1], delay_range[2]))
  }

  # --------------------------------------------------------------------------
  # Panel letter sequences across all rate groups
  # --------------------------------------------------------------------------
  # Each rate group has: CDR row (5 panels) + temperature row (5 panels)
  # Letters run a-e (CDR), f-j (temp) for rate 1; k-o, p-t for rate 2; etc.
  n_panels_per_rate <- 10  # 5 SSPs x 2 variables
  all_letters       <- letters[seq_len(length(rates) * n_panels_per_rate)]

  # --------------------------------------------------------------------------
  # Build rate group rows and collect them for assembly
  # --------------------------------------------------------------------------
  rate_rows_list <- vector("list", length(rates))
  rate_labels_for_title <- character(length(rates))

  for (ri in seq_along(rates)) {

    rate      <- rates[ri]
    run_info_rate <- delay_growth_results[[rate]][["run_info"]]
    rate_r    <- if (!is.null(run_info_rate)) run_info_rate$r_value else NA

    rate_labels_for_title[ri] <- sprintf(
      "%s growth rate (r = %.2f)",
      paste0(toupper(substr(rate, 1, 1)), substr(rate, 2, nchar(rate))),
      rate_r
    )

    # Letter offsets for this rate group
    offset_cdr  <- (ri - 1) * n_panels_per_rate
    offset_temp <- offset_cdr + 5

    panel_letters_cdr  <- all_letters[(offset_cdr  + 1):(offset_cdr  + 5)]
    panel_letters_temp <- all_letters[(offset_temp + 1):(offset_temp + 5)]

    cdr_panels  <- vector("list", length(ssp_names_dg))
    temp_panels <- vector("list", length(ssp_names_dg))

    for (si in seq_along(ssp_names_dg)) {
      ssp    <- ssp_names_dg[si]
      show_y <- (si == 1)

      p_cdr <- plot_delay_cdr(
        delay_growth_results = delay_growth_results,
        rate                 = rate,
        ssp                  = ssp,
        y_limits             = cdr_y_limits,
        delay_range          = delay_range,
        line_alpha           = line_alpha,
        show_y_axis          = show_y,
        show_x_axis          = FALSE,
        show_capacity_limit  = show_capacity_limit
      )
      cdr_panels[[si]] <- p_cdr +
        labs(subtitle = bquote(
          bold(.(paste0(panel_letters_cdr[si], ")"))) ~ .(ssp)
        )) +
        theme(plot.subtitle = element_text(size = 9, hjust = 0, colour = "black"))

      p_temp <- plot_delay_temperature(
        delay_growth_results = delay_growth_results,
        rate                 = rate,
        ssp                  = ssp,
        y_limits             = temp_y_limits,
        delay_range          = delay_range,
        line_alpha           = line_alpha,
        show_y_axis          = show_y,
        show_x_axis          = TRUE
      )
      temp_panels[[si]] <- p_temp +
        labs(subtitle = bquote(
          bold(.(paste0(panel_letters_temp[si], ")"))) ~ .(ssp)
        )) +
        theme(plot.subtitle = element_text(size = 9, hjust = 0, colour = "black"))
    }

    # Tag the first panel of each row with the variable name
    cdr_panels[[1]] <- cdr_panels[[1]] +
      labs(tag = "CDR deployment") +
      theme(plot.tag          = element_text(size = 9, face = "bold", hjust = 0),
            plot.tag.position = "top")

    temp_panels[[1]] <- temp_panels[[1]] +
      labs(tag = "Temperature trajectories") +
      theme(plot.tag          = element_text(size = 9, face = "bold", hjust = 0),
            plot.tag.position = "top")

    cdr_row  <- wrap_plots(cdr_panels,  nrow = 1)
    temp_row <- wrap_plots(temp_panels, nrow = 1)

    rate_rows_list[[ri]] <- cdr_row / temp_row
  }

  # --------------------------------------------------------------------------
  # Shared colour bar legend
  # --------------------------------------------------------------------------
  dummy_df <- data.frame(
    x     = 1:2,
    y     = 1:2,
    delay = delay_range
  )
  legend_plot <- ggplot(dummy_df, aes(x = x, y = y, colour = delay)) +
    geom_point() +
    scale_colour_viridis_c(
      option    = "plasma",
      direction = 1,
      limits    = delay_range,
      name      = "CDR delay (years)",
      guide     = guide_colourbar(
        title.position = "left",
        barwidth       = unit(6, "cm"),
        barheight      = unit(0.4, "cm"),
        title.vjust    = 0.8
      )
    ) +
    theme_bw() +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8)
    )
  colour_legend <- cowplot::get_legend(legend_plot)

  # --------------------------------------------------------------------------
  # Rate group title labels (one per rate group, above each pair of rows)
  # --------------------------------------------------------------------------
  title_grobs <- lapply(rate_labels_for_title, function(lbl) {
    cowplot::ggdraw() +
      cowplot::draw_label(lbl, fontface = "bold", size = 10, hjust = 0.5)
  })

  # Interleave title grobs with rate row patchworks
  n_rates         <- length(rates)
  all_grobs       <- vector("list", n_rates * 2 + 1)
  rel_heights_vec <- numeric(n_rates * 2 + 1)

  for (ri in seq_len(n_rates)) {
    all_grobs[[ri * 2 - 1]]       <- title_grobs[[ri]]
    all_grobs[[ri * 2]]           <- rate_rows_list[[ri]]
    rel_heights_vec[ri * 2 - 1]   <- 0.04
    rel_heights_vec[ri * 2]       <- 1
  }
  all_grobs[[n_rates * 2 + 1]]       <- colour_legend
  rel_heights_vec[n_rates * 2 + 1]   <- 0.12

  final_plot <- cowplot::plot_grid(
    plotlist    = all_grobs,
    ncol        = 1,
    rel_heights = rel_heights_vec
  )

  # --------------------------------------------------------------------------
  # Save if requested
  # --------------------------------------------------------------------------
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0(
        "cdr_delay_growth_dashboard_",
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

# ==============================================================================
# Shared aesthetics for line/threshold plots
# ==============================================================================

growth_rate_colours_dg <- c(
  "slow"     = "#7B2D8B",
  "moderate" = "#1E8BC3",
  "fast"     = "#E67E22"
)

growth_rate_linetypes_dg <- c(
  "slow"     = "dotted",
  "moderate" = "dashed",
  "fast"     = "solid"
)

growth_rate_labels_dg <- c(
  "slow"     = "Slow",
  "moderate" = "Moderate",
  "fast"     = "Fast"
)

# ==============================================================================
# Heatmap functions
# ==============================================================================

#' @title Extract Delay Analysis Summary Into Tidy Data Frame
#' @description
#' Flattens the nested cdr_delay_growth_results list into a single tidy data
#' frame with one row per (rate, delay, SSP) combination. Derives an outcome
#' category from peak and final temperature using supplied thresholds.
#'
#' @param delay_growth_results Output from run_cdr_delay_growth_analysis().
#' @param peak_temp_threshold Minimum peak temperature to count as an overshoot
#'   (default: 1.51). Below this, outcome is "No overshoot".
#' @param final_temp_threshold Maximum final temperature to count as recovered
#'   (default: 1.501). Above this, outcome is "Unrecoverable overshoot".
#' @return Data frame with columns: rate, delay, ssp, peak_temperature,
#'   final_temperature, years_above_1p5, outcome.
extract_delay_summary <- function(delay_growth_results,
                                   peak_temp_threshold  = 1.51,
                                   final_temp_threshold = 1.501) {

  rate_labels <- setdiff(names(delay_growth_results), "run_info")
  rows        <- list()

  for (rate in rate_labels) {
    rate_results <- delay_growth_results[[rate]]
    delays       <- get_rate_delays(rate_results)

    for (delay in delays) {
      delay_key  <- paste0("delay_", delay)
      run_result <- rate_results[[delay_key]]
      if (is.null(run_result$comparison_summary)) next

      for (k in seq_len(nrow(run_result$comparison_summary))) {
        row     <- run_result$comparison_summary[k, ]
        ssp_clean <- gsub("-Baseline$", "", row$scenario)

        if (row$peak_temperature < peak_temp_threshold) {
          outcome <- "No overshoot"
        } else if (row$final_temperature <= final_temp_threshold) {
          outcome <- "Recoverable overshoot"
        } else {
          outcome <- "Unrecoverable overshoot"
        }

        rows[[length(rows) + 1]] <- data.frame(
          rate              = rate,
          delay             = delay,
          ssp               = ssp_clean,
          peak_temperature  = row$peak_temperature,
          final_temperature = row$final_temperature,
          years_above_1p5   = row$years_above_1p5,
          outcome           = outcome,
          stringsAsFactors  = FALSE
        )
      }
    }
  }

  result <- do.call(rbind, rows)

  result$rate    <- factor(result$rate,
                           levels = rate_labels)
  result$outcome <- factor(result$outcome,
                           levels = c("No overshoot",
                                      "Recoverable overshoot",
                                      "Unrecoverable overshoot"))
  result$ssp     <- factor(result$ssp,
                           levels = paste0("SSP", 1:5))

  result
}


#' @title Create CDR Delay Heatmap
#' @description
#' Creates a heatmap with CDR delay on the x-axis, growth rate on the y-axis,
#' and SSP scenarios as column facets. The fill variable is controlled by the
#' \code{variable} argument:
#' \itemize{
#'   \item \code{"outcome"}: categorical fill (no overshoot / recoverable /
#'     unrecoverable), using SSP-matched colours.
#'   \item \code{"peak_temperature"}: continuous fill showing peak temperature.
#'   \item \code{"years_above_1p5"}: continuous fill showing years above 1.5C.
#' }
#'
#' @param summary_df Tidy data frame from extract_delay_summary().
#' @param variable Character; one of "outcome", "peak_temperature", or
#'   "years_above_1p5".
#' @param peak_temp_threshold Threshold line value overlaid on continuous plots
#'   (default: 1.5). Pass NULL to suppress.
#' @param save_plot Logical; save to figs/ directory (default FALSE).
#' @param filename Character; output filename (auto-generated if NULL).
#' @param width Plot width in mm (default 297).
#' @param height Plot height in mm (default 100).
#' @param verbose Logical; print progress messages (default TRUE).
#' @return ggplot object.
#'
#' @examples
#' summary_df <- extract_delay_summary(cdr_delay_growth_results)
#'
#' # Outcome heatmap
#' p1 <- create_delay_heatmap(summary_df, variable = "outcome", save_plot = TRUE)
#'
#' # Peak temperature heatmap
#' p2 <- create_delay_heatmap(summary_df, variable = "peak_temperature", save_plot = TRUE)
#'
#' # Years above 1.5C heatmap
#' p3 <- create_delay_heatmap(summary_df, variable = "years_above_1p5", save_plot = TRUE)
create_delay_heatmap <- function(summary_df,
                                  variable  = c("outcome",
                                                "peak_temperature",
                                                "years_above_1p5"),
                                  save_plot = FALSE,
                                  filename  = NULL,
                                  width     = 297,
                                  height    = 100,
                                  verbose   = TRUE) {

  variable <- match.arg(variable)

  # Infer tile width from delay step size in the data
  delays     <- sort(unique(summary_df$delay))
  tile_width <- if (length(delays) > 1) delays[2] - delays[1] else 1

  # Growth rate factor ordered slow -> moderate -> fast
  rate_levels <- levels(summary_df$rate)

  # Reverse factor so slow is at the top of the y-axis
  summary_df$rate <- factor(summary_df$rate,
                             levels = rev(rate_levels))

  # SSP colours matching the rest of the codebase
  ssp_colours <- c(
    SSP1 = "#00ADCF",
    SSP2 = "#173C66",
    SSP3 = "#F0E442",
    SSP4 = "#E71D25",
    SSP5 = "#951B1E"
  )

  # Capitalise rate labels for display
  capitalise <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

  # --------------------------------------------------------------------------
  # Build plot
  # --------------------------------------------------------------------------

  base <- ggplot(summary_df,
                 aes(x = delay, y = rate, width = tile_width, height = 0.9))

  if (variable == "outcome") {

    # Outcome colours: derive light/dark shades per SSP from ssp_colours
    # Use a consistent 3-shade scheme so panels are comparable across SSPs
    outcome_colours <- c(
      "No overshoot"           = "#FFFFFF",
      "Recoverable overshoot"  = "#AECDE0",
      "Unrecoverable overshoot" = "#2C6E8A"
    )

    p <- base +
      geom_tile(aes(fill = outcome), colour = NA) +
      scale_fill_manual(
        name   = "Outcome",
        values = outcome_colours,
        drop   = FALSE
      ) +
      labs(
        x     = "CDR delay (years)",
        y     = "Growth rate",
        title = NULL
      )

  } else if (variable == "peak_temperature") {

    p <- base +
      geom_tile(aes(fill = peak_temperature), colour = NA) +
      scale_fill_viridis_c(
        option = "inferno",
        direction = 1,
        name   = "Peak temperature (°C)"
      ) +
      geom_hline(yintercept = -Inf, colour = NA) +
      labs(
        x     = "CDR delay (years)",
        y     = "Growth rate",
        title = NULL
      )

  } else {

    p <- base +
      geom_tile(aes(fill = years_above_1p5), colour = NA) +
      scale_fill_viridis_c(
        option = "mako",
        direction = 1,
        name   = "Years above 1.5°C"
      ) +
      labs(
        x     = "CDR delay (years)",
        y     = "Growth rate",
        title = NULL
      )
  }

  p <- p +
    facet_wrap(~ ssp, nrow = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(
      labels = function(x) capitalise(x)
    ) +
    delay_growth_theme +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      strip.background = element_blank(),
      strip.text       = element_text(size = 9, face = "bold"),
      panel.border     = element_rect(colour = "grey80", fill = NA),
      axis.text.y      = element_text(size = 8),
      axis.title       = element_text(size = 9)
    )

  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0(
        "cdr_delay_heatmap_", variable, "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    }

    filepath     <- here::here("figs", filename)
    png_filename <- sub("\\.pdf$", ".png", filename)
    png_filepath <- here::here("figs", png_filename)

    ggsave(filepath, p,
           width = width, height = height, units = "mm",
           device = cairo_pdf)
    ggsave(png_filepath, p,
           width = width, height = height, units = "mm",
           device = "png", dpi = 300, bg = "white")

    if (verbose) {
      cat("Saved:", filepath, "\n")
      cat("PNG:", png_filepath, "\n")
    }
  }

  return(p)
}

# ==============================================================================
# Line plot and threshold plot functions
# ==============================================================================

#' @title Create Delay Line Plot
#' @description
#' Plots a continuous outcome metric (peak temperature or years above 1.5C)
#' against CDR delay for all three growth rates. One line per growth rate,
#' distinguished by colour and linetype. SSP scenarios are column facets.
#'
#' @param summary_df Tidy data frame from extract_delay_summary().
#' @param variable Character; "peak_temperature" or "years_above_1p5".
#' @param save_plot Logical; save to figs/ (default FALSE).
#' @param filename Character; output filename (auto-generated if NULL).
#' @param width Plot width in mm (default 297).
#' @param height Plot height in mm (default 100).
#' @param verbose Logical; print progress (default TRUE).
#' @return ggplot object.
create_delay_line_plot <- function(summary_df,
                                    variable  = c("peak_temperature",
                                                  "years_above_1p5"),
                                    save_plot = FALSE,
                                    filename  = NULL,
                                    width     = 297,
                                    height    = 100,
                                    verbose   = TRUE) {

  variable <- match.arg(variable)

  # Ensure growth rate factor is in the correct order for legend
  summary_df$rate <- factor(summary_df$rate,
                             levels = names(growth_rate_colours_dg))

  if (variable == "peak_temperature") {
    y_var   <- "peak_temperature"
    y_label <- "Peak temperature (°C)"
    ref_y   <- 1.5
    file_id <- "peak_temperature"
  } else {
    y_var   <- "years_above_1p5"
    y_label <- "Years above 1.5°C"
    ref_y   <- 0
    file_id <- "years_above_1p5"
  }

  p <- ggplot(summary_df,
              aes(x        = delay,
                  y        = .data[[y_var]],
                  colour   = rate,
                  linetype = rate,
                  group    = rate)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ ssp, nrow = 1) +
    scale_colour_manual(
      name   = "Growth rate",
      values = growth_rate_colours_dg,
      labels = growth_rate_labels_dg
    ) +
    scale_linetype_manual(
      name   = "Growth rate",
      values = growth_rate_linetypes_dg,
      labels = growth_rate_labels_dg
    ) +
    scale_x_continuous(
      name   = "CDR delay (years)",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      name   = y_label,
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    delay_growth_theme +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      strip.background = element_blank(),
      strip.text       = element_text(size = 9, face = "bold"),
      panel.border     = element_rect(colour = "grey80", fill = NA)
    )

  # Add 1.5C reference line for peak temperature
  if (variable == "peak_temperature") {
    p <- p +
      geom_hline(yintercept = 1.5, linetype = "dashed",
                 colour = "red", alpha = 0.7, linewidth = 0.5)
  }

  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0(
        "cdr_delay_lineplot_", file_id, "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    }
    filepath     <- here::here("figs", filename)
    png_filepath <- here::here("figs", sub("\\.pdf$", ".png", filename))

    ggsave(filepath,     p, width = width, height = height,
           units = "mm", device = cairo_pdf)
    ggsave(png_filepath, p, width = width, height = height,
           units = "mm", device = "png", dpi = 300, bg = "white")

    if (verbose) {
      cat("Saved:", filepath, "\n")
      cat("PNG:", png_filepath, "\n")
    }
  }

  return(p)
}


#' @title Compute Last Feasible Delay Thresholds
#' @description
#' For each (rate, SSP) combination, finds the last CDR delay (in years) at
#' which the outcome is still recoverable or no overshoot. This is the latest
#' a decision-maker can act before recovery by 2100 becomes impossible.
#'
#' Combinations that are always feasible (never unrecoverable within the tested
#' range) are flagged with always_feasible = TRUE. Combinations that are never
#' feasible (unrecoverable even at delay = 0) are flagged with never_feasible.
#'
#' @param summary_df Tidy data frame from extract_delay_summary().
#' @return Data frame with columns: rate, ssp, last_feasible_delay,
#'   always_feasible, never_feasible.
compute_delay_thresholds <- function(summary_df) {

  rate_levels <- levels(summary_df$rate)
  ssp_levels  <- levels(summary_df$ssp)
  rows        <- list()

  for (rate in rate_levels) {
    for (ssp in ssp_levels) {
      sub <- summary_df[summary_df$rate == rate & summary_df$ssp == ssp, ]
      if (nrow(sub) == 0) next

      sub <- sub[order(sub$delay), ]

      feasible_delays   <- sub$delay[sub$outcome != "Unrecoverable overshoot"]
      always_feasible   <- length(feasible_delays) == nrow(sub)
      never_feasible    <- length(feasible_delays) == 0

      last_feasible_delay <- if (!never_feasible) max(feasible_delays) else NA_real_

      rows[[length(rows) + 1]] <- data.frame(
        rate                = rate,
        ssp                 = ssp,
        last_feasible_delay = last_feasible_delay,
        always_feasible     = always_feasible,
        never_feasible      = never_feasible,
        stringsAsFactors    = FALSE
      )
    }
  }

  result      <- do.call(rbind, rows)
  result$rate <- factor(result$rate, levels = rate_levels)
  result$ssp  <- factor(result$ssp,  levels = ssp_levels)
  result
}


#' @title Create Delay Threshold Plot
#' @description
#' Dot plot showing the last CDR delay at which recovery to 1.5C by 2100 is
#' still possible, for each SSP x growth rate combination. Points are grouped
#' by growth rate (colour + shape) and positioned on the x-axis by SSP.
#'
#' Combinations that are always feasible within the tested range are shown
#' with an upward-pointing triangle at max_delay and labelled accordingly.
#' Combinations that are never feasible are shown at y = 0 with a downward
#' triangle.
#'
#' @param summary_df Tidy data frame from extract_delay_summary().
#' @param save_plot Logical; save to figs/ (default FALSE).
#' @param filename Character; output filename (auto-generated if NULL).
#' @param width Plot width in mm (default 180).
#' @param height Plot height in mm (default 100).
#' @param verbose Logical; print progress (default TRUE).
#' @return ggplot object.
create_delay_threshold_plot <- function(summary_df,
                                         save_plot = FALSE,
                                         filename  = NULL,
                                         width     = 180,
                                         height    = 100,
                                         verbose   = TRUE) {

  threshold_df <- compute_delay_thresholds(summary_df)

  max_delay <- max(summary_df$delay, na.rm = TRUE)

  # Shape mapping: normal = circle, always feasible = triangle up,
  # never feasible = triangle down
  threshold_df$point_shape <- 16L  # filled circle
  threshold_df$point_shape[threshold_df$always_feasible]  <- 24L  # triangle up
  threshold_df$point_shape[threshold_df$never_feasible]   <- 25L  # triangle down

  # For plotting, cap always_feasible at max_delay and never_feasible at 0
  threshold_df$plot_y <- threshold_df$last_feasible_delay
  threshold_df$plot_y[threshold_df$always_feasible] <- max_delay
  threshold_df$plot_y[threshold_df$never_feasible]  <- 0

  threshold_df$rate <- factor(threshold_df$rate,
                               levels = names(growth_rate_colours_dg))

  p <- ggplot(threshold_df,
              aes(x      = ssp,
                  y      = plot_y,
                  colour = rate,
                  fill   = rate,
                  shape  = factor(point_shape),
                  group  = rate)) +
    geom_hline(yintercept = seq(0, max_delay, by = 10),
               colour = "grey90", linewidth = 0.3) +
    geom_point(size     = 3,
               position = position_dodge(width = 0.4),
               stroke   = 0.5) +
    scale_colour_manual(
      name   = "Growth rate",
      values = growth_rate_colours_dg,
      labels = growth_rate_labels_dg
    ) +
    scale_fill_manual(
      name   = "Growth rate",
      values = growth_rate_colours_dg,
      labels = growth_rate_labels_dg
    ) +
    scale_shape_manual(
      values = c("16" = 16, "24" = 24, "25" = 25),
      guide  = "none"
    ) +
    scale_y_continuous(
      name   = "Last feasible CDR delay (years)",
      limits = c(-2, max_delay + 2),
      breaks = seq(0, max_delay, by = 10),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      x       = NULL,
      caption = paste0(
        "Point at ", max_delay,
        " yrs (▲) = still feasible at maximum tested delay. ",
        "Point at 0 yrs (▼) = never feasible."
      )
    ) +
    delay_growth_theme +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      panel.border     = element_rect(colour = "grey80", fill = NA),
      plot.caption     = element_text(size = 7, colour = "grey40", hjust = 0)
    )

  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0(
        "cdr_delay_threshold_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    }
    filepath     <- here::here("figs", filename)
    png_filepath <- here::here("figs", sub("\\.pdf$", ".png", filename))

    ggsave(filepath,     p, width = width, height = height,
           units = "mm", device = cairo_pdf)
    ggsave(png_filepath, p, width = width, height = height,
           units = "mm", device = "png", dpi = 300, bg = "white")

    if (verbose) {
      cat("Saved:", filepath, "\n")
      cat("PNG:", png_filepath, "\n")
    }
  }

  return(p)
}
