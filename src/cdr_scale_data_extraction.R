# ==============================================================================
# CDR Scale Sensitivity Data Extraction Functions
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
# Extracts key summary metrics from CDR scale sensitivity results for reporting
# and downstream analysis. Produces three outputs:
#   - feasible_min: minimum r and minimum K that yield a feasible recovery per SSP
#   - cdr_threshold: boundary CDR units separating recoverable/unrecoverable runs
#   - min_unrecoverable: minimum peak temperature at which recovery becomes impossible
#
# Version: 1.0.0
# Last updated: June 2026
# ==============================================================================


#' @title Logistic CDR Growth Function
#' @description
#' Evaluates the logistic CDR deployment curve at a given time. Used internally
#' to compute area-under-curve (AUC) for feasible_min.
#'
#' @param t Time elapsed since t_start (years)
#' @param r CDR growth rate
#' @param K CDR carrying capacity (GtCO2/year)
#' @param g_initial Starting CDR level (GtCO2/year)
#' @param t0 Reference time offset (default 0)
#'
#' @return CDR deployment level at time t
.cdr_logistic <- function(t, r, K, g_initial = 2, t0 = 0) {
  K / (1 + (K / g_initial - 1) * exp(-r * (t - t0)))
}


#' @title Extract CDR Scale Sensitivity Metrics
#' @description
#' Computes three summary tables from CDR scale sensitivity results:
#'
#' \describe{
#'   \item{feasible_min}{For each SSP, the (r, K) pair with the lowest r and the
#'     (r, K) pair with the lowest K that still achieves a feasible recovery
#'     (final temperature <= 1.5°C) from a genuine overshoot (peak >= threshold).
#'     Includes area-under-curve (AUC) for total CDR deployment over 75 years.}
#'   \item{cdr_threshold}{For each SSP, the maximum CDR units observed in failed
#'     runs and the minimum CDR units observed in successful runs — defining the
#'     boundary between recoverable and unrecoverable scenarios.}
#'   \item{min_unrecoverable}{For each SSP, the minimum rounded peak temperature
#'     at which ALL higher peak temperatures are also unrecoverable.}
#' }
#'
#' Results are printed and saved as CSV files to the output/ directory.
#'
#' @param sensitivity_results Output from run_cdr_scale_sensitivity()
#' @param emissions_df Emissions data frame (from interpolate_ssp_emissions())
#' @param peak_temp_min Minimum peak temperature for a run to count as a genuine
#'   overshoot in feasible_min (default 1.51°C)
#' @param final_temp_max Maximum final temperature for a run to be considered
#'   feasible in feasible_min (default 1.501°C)
#' @param g_initial Starting CDR level used in the sensitivity run (GtCO2/year)
#' @param auc_horizon Years over which to integrate CDR deployment for AUC
#'   (default 75, matching a 2025-2100 horizon)
#' @param verbose Print progress messages (default TRUE)
#'
#' @return Named list with elements: feasible_min, cdr_threshold, min_unrecoverable
#'
#' @examples
#' cdr_metrics <- extract_cdr_scale_metrics(
#'   sensitivity_results = cdr_scale_sensitivity_results,
#'   emissions_df        = emissions_df,
#'   g_initial           = 2
#' )
extract_cdr_scale_metrics <- function(sensitivity_results,
                                      emissions_df,
                                      peak_temp_min   = 1.51,
                                      final_temp_max  = 1.501,
                                      g_initial       = 2,
                                      auc_horizon     = 75,
                                      verbose         = TRUE) {

  results <- sensitivity_results$combined_results

  # ---------------------------------------------------------------------------
  # 1. feasible_min: lowest r and lowest K achieving feasible recovery per SSP
  # ---------------------------------------------------------------------------
  feasible_runs <- results |>
    dplyr::filter(feasible == TRUE,
                  final_temperature <= final_temp_max,
                  peak_temperature  >= peak_temp_min)

  feasible_min <- dplyr::bind_rows(
    feasible_runs |>
      dplyr::group_by(scenario) |>
      dplyr::slice_min(order_by = r, n = 1, with_ties = FALSE) |>
      dplyr::mutate(selected_by = "min_r"),
    feasible_runs |>
      dplyr::group_by(scenario) |>
      dplyr::slice_min(order_by = K, n = 1, with_ties = FALSE) |>
      dplyr::mutate(selected_by = "min_K")
  ) |>
    dplyr::select(scenario, selected_by, r, K,
                  final_temperature, peak_temperature, years_above_1p5,
                  total_cdr_units, mitig_cost, remov_cost, temp_cost, total_cost) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      AUC = integrate(.cdr_logistic, lower = 0, upper = auc_horizon,
                      r = r, K = K, g_initial = g_initial)$value
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(scenario, selected_by)

  # ---------------------------------------------------------------------------
  # 2. cdr_threshold: CDR units at boundary between recoverable/unrecoverable
  # ---------------------------------------------------------------------------
  cdr_threshold <- dplyr::bind_rows(
    results |>
      dplyr::filter(feasible == FALSE) |>
      dplyr::group_by(scenario) |>
      dplyr::summarise(
        max_cdr_unrecoverable = max(total_cdr_units, na.rm = TRUE),
        .groups = "drop"
      ),
    results |>
      dplyr::filter(feasible == TRUE) |>
      dplyr::group_by(scenario) |>
      dplyr::summarise(
        min_cdr_recoverable = min(total_cdr_units, na.rm = TRUE),
        .groups = "drop"
      )
  ) |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ dplyr::first(stats::na.omit(.))),
                     .groups = "drop") |>
    dplyr::arrange(scenario)

  # Cumulative emissions context (2025 onward)
  cumulative_emissions <- emissions_df |>
    dplyr::filter(Year >= 2025) |>
    dplyr::group_by(Scenario) |>
    dplyr::summarise(
      cumulative_emissions = sum(Value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(Scenario)

  cdr_threshold_context <- cdr_threshold |>
    dplyr::left_join(cumulative_emissions, by = c("scenario" = "Scenario")) |>
    dplyr::mutate(
      expected_cdr = cumulative_emissions - 650,
      difference   = max_cdr_unrecoverable - expected_cdr
    )

  # ---------------------------------------------------------------------------
  # 3. min_unrecoverable: lowest peak temp at which all higher temps also fail
  # ---------------------------------------------------------------------------
  min_unrecoverable <- results |>
    dplyr::mutate(peak_temp_rounded = round(peak_temperature, 2)) |>
    dplyr::group_by(scenario, peak_temp_rounded) |>
    dplyr::summarise(any_feasible = any(feasible), .groups = "drop") |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(
      threshold = min(peak_temp_rounded[
        !any_feasible &
          sapply(peak_temp_rounded, function(pt) {
            all(!any_feasible[peak_temp_rounded >= pt])
          })
      ], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(scenario)

  # ---------------------------------------------------------------------------
  # Print and save
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("\n=== FEASIBLE MINIMUM (r and K) PER SSP ===\n")
    print(feasible_min)
    cat("\n=== CDR THRESHOLD (recoverable vs unrecoverable) ===\n")
    print(cdr_threshold_context)
    cat("\n=== MINIMUM UNRECOVERABLE PEAK TEMPERATURE PER SSP ===\n")
    print(min_unrecoverable)
  }

  here::here("output") |> dir.create(showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(feasible_min,     here::here("output", "feasible_min.csv"),    row.names = FALSE)
  utils::write.csv(cdr_threshold,    here::here("output", "cdr_threshold.csv"),   row.names = FALSE)
  utils::write.csv(min_unrecoverable,here::here("output", "min_unrecoverable.csv"),row.names = FALSE)

  if (verbose) cat("\nCSV files saved to output/\n")

  invisible(list(
    feasible_min       = feasible_min,
    cdr_threshold      = cdr_threshold_context,
    min_unrecoverable  = min_unrecoverable
  ))
}
