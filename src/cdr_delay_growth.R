# ==============================================================================
# CDR Delay x Growth Rate Analysis Functions
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

#' @title Run CDR Delay x Growth Rate Analysis
#' @description
#' Sweeps over a named set of CDR logistic growth rates (r) and a range of CDR
#' deployment delays, running the full SSP scenario comparison once per
#' (rate, delay) pair. Mitigation is eliminated via make_zero_capacity() so CDR
#' is the sole control throughout.
#'
#' Results are stored as a nested list: results[[rate_label]][[delay_key]], where
#' delay_key is e.g. "delay_0", "delay_10". Each entry is the full output of
#' run_scenario_comparison(). A "run_info" element is added at both the rate
#' level and the top level to record metadata.
#'
#' Intermediate per-rate RDS files are saved after each rate completes, so a
#' crash part-way through does not lose completed work.
#'
#' @param parameter_df Single-row data frame of model parameters.
#' @param emissions_df Data frame with emissions data for all SSP scenarios.
#' @param economic_df Data frame with economic data for all SSP scenarios.
#' @param scenarios Character vector of SSP scenario names to compare.
#' @param r_values Named numeric vector of CDR logistic growth rates to sweep.
#'   Names become list keys in the returned object and in saved filenames
#'   (e.g. c(slow = 0.05, moderate = 0.07, fast = 0.11)).
#' @param K CDR logistic carrying capacity in GtCO2/year (default: 100).
#' @param g_initial Starting CDR deployment level in GtCO2/year (default: 2).
#' @param t_start Year CDR deployment begins (default: 2025).
#' @param max_delay_years Maximum CDR delay to test in years (default: 70).
#' @param delay_step_size Step between tested delay values in years (default: 1).
#'   Use a larger value (e.g. 5 or 10) for fast exploratory runs; set to 1 for
#'   publication-quality results.
#' @param use_mitigation_capacity_limit Logical; set TRUE to apply the
#'   mitigation capacity constraint (default: TRUE). Should remain TRUE when
#'   using make_zero_capacity() to eliminate mitigation.
#' @param mitigation_capacity_function Capacity function for mitigation. Use
#'   make_zero_capacity() to eliminate mitigation entirely (default).
#' @param use_parallel Logical; use parallel processing across SSPs within each
#'   run_scenario_comparison() call (default: TRUE).
#' @param save_results Logical; save intermediate per-rate RDS files after each
#'   rate completes, and a combined RDS on completion (default: TRUE).
#' @param output_dir Directory for saved files (default: "output").
#' @param output_prefix Filename prefix for all saved files
#'   (default: "cdr_delay_growth").
#' @param verbose Logical; print progress messages (default: TRUE).
#'
#' @return Named list. Rate-level entries are keyed by the names of r_values
#'   (e.g. "slow", "moderate", "fast"); each contains:
#'   \describe{
#'     \item{delay_N}{Full run_scenario_comparison() output for delay N years.}
#'     \item{run_info}{List with rate_label, r_value, K, g_initial, t_start, delays.}
#'   }
#'   A top-level "run_info" entry records analysis-wide metadata.
#'
#' @examples
#' cdr_delay_growth_results <- run_cdr_delay_growth_analysis(
#'   parameter_df    = parameter_df[1, ],
#'   emissions_df    = emissions_df,
#'   economic_df     = economic_df,
#'   scenarios       = c("SSP1-Baseline", "SSP2-Baseline", "SSP3-Baseline",
#'                       "SSP4-Baseline", "SSP5-Baseline"),
#'   r_values        = c(slow = 0.05, moderate = 0.07, fast = 0.11),
#'   K               = 100,
#'   max_delay_years = 70,
#'   delay_step_size = 10
#' )
#'
#' # Access temperature trajectory for SSP2, slow rate, 20-year delay
#' cdr_delay_growth_results$slow$delay_20$scenario_results[["SSP2-Baseline"]]$temperature_anomaly

run_cdr_delay_growth_analysis <- function(parameter_df,
                                           emissions_df,
                                           economic_df,
                                           scenarios,
                                           r_values,
                                           K                             = 100,
                                           g_initial                     = 2,
                                           t_start                       = 2025,
                                           max_delay_years               = 70,
                                           delay_step_size               = 1,
                                           use_mitigation_capacity_limit = TRUE,
                                           mitigation_capacity_function  = make_zero_capacity(),
                                           use_parallel                  = TRUE,
                                           save_results                  = TRUE,
                                           output_dir                    = "output",
                                           output_prefix                 = "cdr_delay_growth",
                                           verbose                       = TRUE) {

  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  if (!is.data.frame(parameter_df) || nrow(parameter_df) != 1) {
    stop("parameter_df must be a single-row data frame")
  }
  if (is.null(names(r_values)) || any(names(r_values) == "")) {
    stop("r_values must be a fully named numeric vector ",
         "(e.g., c(slow = 0.05, moderate = 0.07, fast = 0.11))")
  }
  if (any(r_values <= 0)) {
    stop("All r values must be positive")
  }
  if (K <= g_initial) {
    stop("K must be greater than g_initial")
  }
  if (delay_step_size < 1) {
    stop("delay_step_size must be >= 1")
  }
  if (max_delay_years < 0) {
    stop("max_delay_years must be >= 0")
  }

  # --------------------------------------------------------------------------
  # Setup
  # --------------------------------------------------------------------------
  delays      <- seq(0, max_delay_years, by = delay_step_size)
  rate_labels <- names(r_values)
  start_time  <- Sys.time()

  if (verbose) {
    cat("=== CDR DELAY x GROWTH RATE ANALYSIS ===\n")
    cat("Growth rates:", paste(
      paste0(rate_labels, " (r=", r_values, ")"), collapse = ", "), "\n")
    cat("K:", K, "GtCO2/yr | g_initial:", g_initial, "\n")
    cat("CDR delay range: 0 -", max_delay_years,
        "years (step =", delay_step_size, ")\n")
    cat("Delays to test:", length(delays), "\n")
    cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
    cat("Total scenario comparisons:",
        length(r_values) * length(delays), "\n\n")
  }

  # Top-level list: rate entries + run_info
  all_results <- vector("list", length(r_values) + 1)
  names(all_results) <- c(rate_labels, "run_info")

  # --------------------------------------------------------------------------
  # Iterate over growth rates
  # --------------------------------------------------------------------------
  for (i in seq_along(r_values)) {

    rate_label <- rate_labels[i]
    r_val      <- r_values[i]

    if (verbose) {
      cat(strrep("=", 60), "\n")
      cat("Growth rate ", i, " of ", length(r_values),
          ": '", rate_label, "' (r = ", r_val, ")\n", sep = "")
      cat(strrep("=", 60), "\n")
    }

    # Rate-level list: delay entries + run_info
    rate_results <- vector("list", length(delays) + 1)
    names(rate_results) <- c(paste0("delay_", delays), "run_info")

    # -------------------------------------------------------------------------
    # Iterate over delays
    # -------------------------------------------------------------------------
    for (j in seq_along(delays)) {

      delay     <- delays[j]
      delay_key <- paste0("delay_", delay)

      if (verbose) {
        elapsed <- difftime(Sys.time(), start_time, units = "mins")
        cat("  [", i, "/", length(r_values), "] ",
            "Delay ", j, "/", length(delays),
            " (", delay, " yrs) — Elapsed: ",
            sprintf("%.1f", as.numeric(elapsed)), " min\n", sep = "")
      }

      # Capacity curve starts at t_start + delay so the logistic grows from
      # g_initial at the same year CDR deployment actually begins. This
      # represents a delay due to lack of funding/development/political will —
      # no capacity build-up occurs during the delay period.
      cdr_cap_fn <- make_logistic_from_zero(
        g_initial = g_initial,
        K         = K,
        r         = r_val,
        t_start   = t_start + delay
      )

      rate_results[[delay_key]] <- run_scenario_comparison(
        parameter_df                  = parameter_df,
        emissions_df                  = emissions_df,
        economic_df                   = economic_df,
        scenarios                     = scenarios,
        mitigation_delay_years        = 0,
        cdr_delay_years               = delay,
        use_mitigation_capacity_limit = use_mitigation_capacity_limit,
        mitigation_capacity_function  = mitigation_capacity_function,
        use_cdr_capacity_limit        = TRUE,
        cdr_capacity_function         = cdr_cap_fn,
        use_parallel                  = use_parallel,
        save_results                  = FALSE,
        verbose                       = FALSE
      )
    }

    rate_results[["run_info"]] <- list(
      rate_label = rate_label,
      r_value    = r_val,
      K          = K,
      g_initial  = g_initial,
      t_start    = t_start,
      delays     = delays
    )

    all_results[[rate_label]] <- rate_results

    # Save intermediate per-rate RDS after each rate completes
    if (save_results) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      rate_file <- here::here(
        output_dir,
        paste0(output_prefix, "_", rate_label, "_", timestamp, ".rds")
      )
      saveRDS(rate_results, rate_file)
      if (verbose) cat("  Saved:", rate_file, "\n")
    }
  }

  # --------------------------------------------------------------------------
  # Top-level metadata and combined save
  # --------------------------------------------------------------------------
  total_time <- difftime(Sys.time(), start_time, units = "mins")

  all_results[["run_info"]] <- list(
    r_values        = r_values,
    K               = K,
    g_initial       = g_initial,
    t_start         = t_start,
    delays          = delays,
    max_delay_years = max_delay_years,
    delay_step_size = delay_step_size,
    scenarios       = scenarios,
    start_time      = start_time,
    end_time        = Sys.time(),
    total_time_mins = as.numeric(total_time)
  )

  if (verbose) {
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("COMPLETE — Total time: ",
        sprintf("%.1f", as.numeric(total_time)), " minutes\n", sep = "")
    cat(strrep("=", 60), "\n")
  }

  if (save_results) {
    timestamp     <- format(Sys.time(), "%Y%m%d_%H%M%S")
    combined_file <- here::here(
      output_dir,
      paste0(output_prefix, "_combined_", timestamp, ".rds")
    )
    saveRDS(all_results, combined_file)
    if (verbose) cat("Combined results saved to:", combined_file, "\n")
  }

  return(all_results)
}
