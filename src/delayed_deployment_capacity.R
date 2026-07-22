# ==============================================================================
# Delayed Deployment Analysis with CDR Capacity Constraint
#
# Part of: MACROM: An Optimal Control Model for Balancing Climate Change Abatement
# and Damage Trade-offs
# Authors: Nina Rynne, Michael Bode, Melanie Roberts, Ryan Heneghan
# Institution: Griffith University
#
# Copyright (c) 2025 Nina Rynne
# Licensed under CC-BY-4.0 - see LICENSE file for details
#
# Description:
# Combines the delayed deployment grid (mitigation_delay x cdr_delay, as in
# delayed_deployment.R) with the CDR logistic capacity cap (as in
# multi_scenario_comparison_analysis_capacity). Mitigation remains a free,
# unconstrained control that is only delayed. CDR is both delayed AND capped
# by a logistic growth curve, so its ramp-up speed is limited by the chosen
# capacity level (e.g. slow/moderate/fast).
#
# Unlike delayed_deployment.R, no bespoke fast solver exists for the combined
# delay + capacity case, so this uses the general capacity-aware solver
# (optimal_control_shooting via run_scenario_comparison). This is considerably
# slower per combination — see the runtime guidance in the calling Rmd chunk.
#
# Version: 1.0.0
# Last updated: July 2026
# ==============================================================================

#' @title Run Delayed Deployment Analysis with CDR Capacity Constraint
#' @description
#' Sweeps the full grid of mitigation_delay x cdr_delay combinations for each
#' of a named set of CDR logistic growth rates (capacity levels), running the
#' full SSP scenario comparison once per (capacity level, mitigation_delay,
#' cdr_delay) combination. Mitigation is left unconstrained (only delayed);
#' CDR is delayed AND capped by a logistic capacity curve that itself starts
#' ramping only once CDR deployment begins (t_start + cdr_delay), so no
#' capacity builds during the CDR delay period.
#'
#' @param parameter_df Single-row data frame of model parameters.
#' @param emissions_df Data frame with emissions data for all SSP scenarios.
#' @param economic_df Data frame with economic data for all SSP scenarios.
#' @param scenarios Character vector of SSP scenario names to compare.
#' @param r_values Named numeric vector of CDR logistic growth rates defining
#'   the capacity levels to compare. Names become the capacity_level labels
#'   (e.g. c(slow = 0.05, moderate = 0.07, fast = 0.11)).
#' @param K CDR logistic carrying capacity in GtCO2/year (default: 100).
#' @param g_initial Starting CDR deployment level in GtCO2/year (default: 2).
#' @param t_start Year CDR deployment nominally begins (default: 2025).
#' @param max_delay_years Maximum delay to test in years, applied to BOTH
#'   mitigation_delay and cdr_delay (default: 70).
#' @param delay_step_size Step size for the delay grid in years (default: 10).
#'   Use 1 for publication-quality results, larger values (5-10) for fast
#'   exploratory runs — see runtime guidance in the calling Rmd chunk.
#' @param target_emissions Target cumulative emissions constraint in GtCO2
#'   (default: uses co2_target_2100 from parameter_df).
#' @param feasibility_tolerance Tolerance for emission gap to be considered
#'   feasible, in GtCO2 (default: 50, matching delayed_deployment.R).
#' @param use_parallel Logical; use parallel processing across SSPs within
#'   each run_scenario_comparison() call (default: TRUE).
#' @param save_results Logical; save intermediate per-capacity-level RDS files
#'   after each level completes, and a combined RDS/CSV on completion
#'   (default: TRUE).
#' @param output_dir Directory for saved files (default: "output").
#' @param output_prefix Filename prefix for saved files
#'   (default: "delayed_deployment_capacity").
#' @param verbose Logical; print progress messages (default: TRUE).
#'
#' @return List containing:
#'   - results_by_capacity: Named list keyed by capacity level, each containing
#'     combined_results (data frame) and run_info
#'   - combined_results: Data frame with all capacity levels combined, with
#'     columns scenario, scenario_short, mitigation_delay, cdr_delay,
#'     capacity_level, feasible, peak_temperature, total_cost, mitig_cost,
#'     remov_cost, temp_cost, years_above_1p5, and other comparison_summary
#'     columns
#'   - run_info: Metadata about the analysis run
#'
#' @examples
#' delayed_deployment_capacity_results <- run_delayed_deployment_capacity_analysis(
#'   parameter_df    = parameter_df[1, ],
#'   emissions_df    = emissions_df,
#'   economic_df     = economic_df,
#'   scenarios       = c("SSP1-Baseline", "SSP2-Baseline", "SSP3-Baseline",
#'                       "SSP4-Baseline", "SSP5-Baseline"),
#'   r_values        = c(slow = 0.05, moderate = 0.07, fast = 0.11),
#'   max_delay_years = 70,
#'   delay_step_size = 10
#' )
run_delayed_deployment_capacity_analysis <- function(parameter_df,
                                                      emissions_df,
                                                      economic_df,
                                                      scenarios,
                                                      r_values               = c(slow = 0.05, moderate = 0.07, fast = 0.11),
                                                      K                      = 100,
                                                      g_initial              = 2,
                                                      t_start                = 2025,
                                                      max_delay_years        = 70,
                                                      delay_step_size        = 10,
                                                      target_emissions       = NULL,
                                                      feasibility_tolerance  = 50,
                                                      use_parallel           = TRUE,
                                                      save_results           = TRUE,
                                                      output_dir             = "output",
                                                      output_prefix          = "delayed_deployment_capacity",
                                                      verbose                = TRUE) {

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

  if (is.null(target_emissions)) {
    target_emissions <- parameter_df$co2_target_2100
  }

  # --------------------------------------------------------------------------
  # Setup
  # --------------------------------------------------------------------------
  delay_sequence <- seq(0, max_delay_years, by = delay_step_size)
  rate_labels    <- names(r_values)
  n_combinations_per_level <- length(delay_sequence)^2
  start_time     <- Sys.time()

  clean_scenario_names <- function(x) gsub("SSP([0-9])-Baseline", "SSP\\1", x)

  if (verbose) {
    cat("=== DELAYED DEPLOYMENT x CDR CAPACITY ANALYSIS ===\n")
    cat("Capacity levels:", paste(
      paste0(rate_labels, " (r=", r_values, ")"), collapse = ", "), "\n")
    cat("K:", K, "GtCO2/yr | g_initial:", g_initial, "\n")
    cat("Delay range: 0 -", max_delay_years, "years (step:", delay_step_size, ")\n")
    cat("Combinations per capacity level:", n_combinations_per_level, "\n")
    cat("Total scenario comparisons:",
        length(r_values) * n_combinations_per_level, "\n\n")
  }

  results_by_capacity <- vector("list", length(rate_labels))
  names(results_by_capacity) <- rate_labels
  combined_rows <- list()

  # --------------------------------------------------------------------------
  # Iterate over capacity levels
  # --------------------------------------------------------------------------
  for (i in seq_along(r_values)) {

    rate_label <- rate_labels[i]
    r_val      <- r_values[i]

    if (verbose) {
      cat(strrep("=", 60), "\n")
      cat("Capacity level ", i, " of ", length(r_values),
          ": '", rate_label, "' (r = ", r_val, ")\n", sep = "")
      cat(strrep("=", 60), "\n")
    }

    level_rows <- list()

    # ------------------------------------------------------------------------
    # Iterate over the mitigation_delay x cdr_delay grid
    # ------------------------------------------------------------------------
    for (mitig_delay in delay_sequence) {
      for (cdr_delay in delay_sequence) {

        if (verbose) {
          elapsed <- difftime(Sys.time(), start_time, units = "mins")
          cat("  [", rate_label, "] mitigation_delay=", mitig_delay,
              ", cdr_delay=", cdr_delay,
              " — Elapsed: ", sprintf("%.1f", as.numeric(elapsed)), " min\n", sep = "")
        }

        # Capacity curve starts at t_start + cdr_delay so the logistic grows
        # from g_initial at the same year CDR deployment actually begins —
        # no capacity build-up occurs during the CDR delay period.
        cdr_cap_fn <- make_logistic_from_zero(
          g_initial = g_initial,
          K         = K,
          r         = r_val,
          t_start   = t_start + cdr_delay
        )

        comparison <- run_scenario_comparison(
          parameter_df                  = parameter_df,
          emissions_df                  = emissions_df,
          economic_df                   = economic_df,
          scenarios                     = scenarios,
          mitigation_delay_years        = mitig_delay,
          cdr_delay_years               = cdr_delay,
          use_mitigation_capacity_limit = FALSE,
          use_cdr_capacity_limit        = TRUE,
          cdr_capacity_function         = cdr_cap_fn,
          use_parallel                  = use_parallel,
          save_results                  = FALSE,
          verbose                       = FALSE
        )

        if (!is.null(comparison$comparison_summary)) {
          summary_df <- comparison$comparison_summary
          summary_df$scenario_short   <- clean_scenario_names(summary_df$scenario)
          summary_df$mitigation_delay <- mitig_delay
          summary_df$cdr_delay        <- cdr_delay
          summary_df$capacity_level   <- rate_label
          summary_df$feasible         <- summary_df$converged &
            abs(summary_df$final_emissions - target_emissions) <= feasibility_tolerance

          level_rows[[length(level_rows) + 1]] <- summary_df
        }
      }
    }

    level_df <- do.call(rbind, level_rows)

    results_by_capacity[[rate_label]] <- list(
      combined_results = level_df,
      run_info = list(
        rate_label = rate_label,
        r_value    = r_val,
        K          = K,
        g_initial  = g_initial,
        t_start    = t_start,
        delays     = delay_sequence
      )
    )

    combined_rows[[rate_label]] <- level_df

    # Save intermediate per-level RDS after each capacity level completes
    if (save_results) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      level_file <- here::here(
        output_dir,
        paste0(output_prefix, "_", rate_label, "_", timestamp, ".rds")
      )
      saveRDS(results_by_capacity[[rate_label]], level_file)
      if (verbose) cat("  Saved:", level_file, "\n")
    }
  }

  # --------------------------------------------------------------------------
  # Combine and finalise
  # --------------------------------------------------------------------------
  combined_results <- do.call(rbind, combined_rows)
  rownames(combined_results) <- NULL

  total_time <- difftime(Sys.time(), start_time, units = "mins")

  n_combinations <- nrow(combined_results)
  n_feasible     <- sum(combined_results$feasible, na.rm = TRUE)

  if (verbose) {
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("COMPLETE — Total time: ",
        sprintf("%.1f", as.numeric(total_time)), " minutes\n", sep = "")
    cat("Feasible:", n_feasible, sprintf("(%.1f%%)", 100 * n_feasible / n_combinations), "\n")
    cat(strrep("=", 60), "\n")
  }

  results_list <- list(
    results_by_capacity = results_by_capacity,
    combined_results     = combined_results,
    summary_stats = list(
      n_combinations   = n_combinations,
      n_feasible       = n_feasible,
      feasibility_rate = n_feasible / n_combinations
    ),
    run_info = list(
      r_values              = r_values,
      K                     = K,
      g_initial             = g_initial,
      t_start               = t_start,
      max_delay_years       = max_delay_years,
      delay_step_size       = delay_step_size,
      delays                = delay_sequence,
      scenarios             = scenarios,
      target_emissions      = target_emissions,
      feasibility_tolerance = feasibility_tolerance,
      total_time_minutes    = as.numeric(total_time)
    )
  )

  if (save_results) {
    timestamp     <- format(Sys.time(), "%Y%m%d_%H%M%S")
    combined_rds  <- here::here(output_dir, paste0(output_prefix, "_combined_", timestamp, ".rds"))
    combined_csv  <- here::here(output_dir, paste0(output_prefix, "_combined_", timestamp, ".csv"))

    saveRDS(results_list, combined_rds)
    write.csv(combined_results, combined_csv, row.names = FALSE)

    if (verbose) {
      cat("Combined RDS saved to:", combined_rds, "\n")
      cat("Combined CSV saved to:", combined_csv, "\n")
    }
  }

  return(results_list)
}
