# ==============================================================================
# CDR Scale Free-Terminal-Time Sensitivity Analysis Functions
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
# Free-terminal-time variant of cdr_scale_sensitivity.R. Sweeps the CDR
# logistic scale parameters K (carrying capacity, GtCO2/year) and r (growth
# rate) across a grid. The analysis is TWO-STAGE per K/r/scenario cell:
#
#   1. Solve the standard fixed-time T = 2100 problem. If the cell is
#      feasible by the existing convention (solver converged and
#      final_temperature <= 1.5 + 0.05°C, exactly as in
#      cdr_scale_sensitivity.R), record T_star = 2100 with
#      optimum_type = "feasible_by_2100" and years_beyond_2100 = 0.
#   2. Only for cells INFEASIBLE at 2100, run
#      optimal_control_free_terminal_time() over [2100, t_max] to find the
#      earliest year cumulative emissions can first be brought down to the
#      1.5°C-consistent target (the year the historical overshoot backlog
#      is cleared — see optimal_control_free_time.R for what this does and
#      doesn't imply about years after T*). The headline output is
#      years_beyond_2100 — how much longer these combinations need beyond
#      the conventional deadline.
#
# WHY THE SEARCH STARTS AT 2100 (not the deployment start year): cells
# infeasible at 2100 are, by construction, still above budget at that point
# — no control path can bring cumulative emissions down to target that
# early. Searching for the EARLIEST reachable year at or after 2100
# therefore always lands on a genuine descent through the target: an
# instant earlier the target wasn't reachable at all, so the earliest
# reachable year cannot be a later, already-past, drifted-back-up touch of
# the target (which is possible at OTHER years further out, once CDR has
# overshot the target and is easing off — see optimal_control_free_time.R
# for the full argument and why no extra "arrived from above" filter is
# needed once the search asks specifically for the earliest year).
#
# Feasibility philosophy (matching cdr_scale_sensitivity.R):
# Every K/r/scenario combination always produces a result row. feasible =
# TRUE when the cell is feasible by 2100 OR the free-time search found a
# genuine, earliest-reachable return year within [2100, t_max]; otherwise
# T_star = NA and metric columns carry the best-effort diagnostics
# available.
#
# Parallel strategy (matching cdr_scale_sensitivity.R): parLapply() across
# K/r combinations, one data frame per combination, bind at the end.
#
# DEPENDENCIES: source these first (see the free_time_scale_sensitivity
# workflow chunk): optimal_control_capacity.R (NOT optimal_control_core.R),
# capacity_helpers.R, optimal_control_free_time.R, cdr_scale_sensitivity.R
# (for build_cdr_scale_grid(), reused unchanged). emissions_df/economic_df
# must extend to t_max.
#
# Version: 2.0.0
# Last updated: July 2026
# ==============================================================================


# ==============================================================================
# Section 1: Results Packagers
# ==============================================================================

#' @title Package Free-Terminal-Time Result for One K/r/Scenario Combination
#' @description
#' Converts one optimal_control_free_terminal_time() result into a
#' standardised single-row data frame for combining across the grid sweep.
#' Metric columns (peak temperature, costs, control totals, ...) are the same
#' set as package_scale_result() in cdr_scale_sensitivity.R but evaluated at
#' T_star rather than 2100, with free-time columns (T_star, optimum_type,
#' hamiltonian_at_T_star, years_beyond_2100) prepended.
#'
#' Column set is fixed regardless of convergence so that bind_rows() can
#' always combine result rows safely.
#'
#' @param ft_result List returned by optimal_control_free_terminal_time()
#' @param K_val Numeric CDR carrying capacity used for this run (GtCO2/year)
#' @param r_val Numeric CDR growth rate used for this run
#' @param scenario Scenario name for this run
#'
#' @return Single-row data frame with the fixed free-time column set
package_free_time_result <- function(ft_result, K_val, r_val, scenario) {

  sol <- ft_result$solution

  # Metrics at T_star from the full inner solution; NA when no solution at
  # all was recovered (solution can be a best-effort attempt even when the
  # search did not converge — keep its diagnostics)
  if (!is.null(sol)) {
    peak_temperature       <- max(sol$temperature_anomaly)
    final_temperature      <- sol$final_temperature
    years_above_1p5        <- sum(sol$temperature_anomaly > 1.5)
    final_emissions        <- sol$final_emissions
    emission_gap           <- sol$emission_gap
    total_cost             <- sol$total_cost
    mitig_cost             <- sol$mitig_cost
    remov_cost             <- sol$remov_cost
    temp_cost              <- sol$temp_cost
    total_mitigation_units <- sum(sol$qty_mitig)
    total_cdr_units        <- sum(sol$qty_remov)
  } else {
    peak_temperature       <- NA_real_
    final_temperature      <- NA_real_
    years_above_1p5        <- NA_real_
    final_emissions        <- NA_real_
    emission_gap           <- NA_real_
    total_cost             <- NA_real_
    mitig_cost             <- NA_real_
    remov_cost             <- NA_real_
    temp_cost              <- NA_real_
    total_mitigation_units <- NA_real_
    total_cdr_units        <- NA_real_
  }

  data.frame(
    K                      = K_val,
    r                      = r_val,
    scenario               = scenario,
    scenario_short         = gsub("SSP([0-9])-Baseline", "SSP\\1", scenario),
    T_star                 = ft_result$T_star,
    optimum_type           = if (is.na(ft_result$optimum_type)) NA_character_
                             else ft_result$optimum_type,
    hamiltonian_at_T_star  = ft_result$hamiltonian_at_T_star,
    years_beyond_2100      = if (!is.na(ft_result$T_star)) {
                               max(0, ft_result$T_star - 2100)
                             } else NA_real_,
    peak_temperature       = peak_temperature,
    final_temperature      = final_temperature,
    years_above_1p5        = years_above_1p5,
    final_emissions        = final_emissions,
    emission_gap           = emission_gap,
    total_cost             = total_cost,
    mitig_cost             = mitig_cost,
    remov_cost             = remov_cost,
    temp_cost              = temp_cost,
    total_mitigation_units = total_mitigation_units,
    total_cdr_units        = total_cdr_units,
    outer_iterations       = ft_result$outer_iterations,
    converged              = ft_result$converged,
    feasible               = ft_result$feasible,
    infeasible_reason      = if (is.na(ft_result$infeasible_reason)) NA_character_
                             else ft_result$infeasible_reason,
    stringsAsFactors       = FALSE
  )
}


#' @title Package Error Result for a Failed Free-Time K/r Combination
#' @description
#' Creates a standardised NA-filled result data frame when a K/r combination
#' fails entirely due to an unhandled solver error. Column set is identical
#' to package_free_time_result() so bind_rows() always combines cleanly.
#'
#' @param K_val Numeric CDR carrying capacity
#' @param r_val Numeric CDR growth rate
#' @param scenarios Character vector of scenario names that were attempted
#' @param error_message Optional error text recorded in infeasible_reason
#'
#' @return Data frame with one row per scenario, all metric columns NA
package_free_time_error_result <- function(K_val, r_val, scenarios,
                                           error_message = NA_character_) {

  data.frame(
    K                      = K_val,
    r                      = r_val,
    scenario               = scenarios,
    scenario_short         = gsub("SSP([0-9])-Baseline", "SSP\\1", scenarios),
    T_star                 = NA_real_,
    optimum_type           = NA_character_,
    hamiltonian_at_T_star  = NA_real_,
    years_beyond_2100      = NA_real_,
    peak_temperature       = NA_real_,
    final_temperature      = NA_real_,
    years_above_1p5        = NA_real_,
    final_emissions        = NA_real_,
    emission_gap           = NA_real_,
    total_cost             = NA_real_,
    mitig_cost             = NA_real_,
    remov_cost             = NA_real_,
    temp_cost              = NA_real_,
    total_mitigation_units = NA_real_,
    total_cdr_units        = NA_real_,
    outer_iterations       = NA_integer_,
    converged              = FALSE,
    feasible               = FALSE,
    infeasible_reason      = error_message,
    stringsAsFactors       = FALSE
  )
}


#' @title Package Feasible-by-2100 Result
#' @description
#' Builds the standard result row for a cell that already meets the 1.5°C
#' target under the conventional fixed T = 2100 solve (stage 1), so the
#' free-terminal-time search is not needed. T_star is recorded as 2100 with
#' optimum_type = "feasible_by_2100" and years_beyond_2100 = 0; metric
#' columns come from the fixed-time solution, matching what the existing
#' fixed-time sensitivity reports for the same cell. hamiltonian_at_T_star
#' is still computed (H at 2100) as a diagnostic.
#'
#' @param solution_2100 Solution list from the fixed-time T = 2100 shooting solve
#' @param K_val Numeric CDR carrying capacity
#' @param r_val Numeric CDR growth rate
#' @param scenario Scenario name
#'
#' @return Single-row data frame with the fixed free-time column set
package_by2100_result <- function(solution_2100, K_val, r_val, scenario) {

  sol <- solution_2100

  data.frame(
    K                      = K_val,
    r                      = r_val,
    scenario               = scenario,
    scenario_short         = gsub("SSP([0-9])-Baseline", "SSP\\1", scenario),
    T_star                 = 2100,
    optimum_type           = "feasible_by_2100",
    hamiltonian_at_T_star  = compute_terminal_hamiltonian(sol),
    years_beyond_2100      = 0,
    peak_temperature       = max(sol$temperature_anomaly),
    final_temperature      = sol$final_temperature,
    years_above_1p5        = sum(sol$temperature_anomaly > 1.5),
    final_emissions        = sol$final_emissions,
    emission_gap           = sol$emission_gap,
    total_cost             = sol$total_cost,
    mitig_cost             = sol$mitig_cost,
    remov_cost             = sol$remov_cost,
    temp_cost              = sol$temp_cost,
    total_mitigation_units = sum(sol$qty_mitig),
    total_cdr_units        = sum(sol$qty_remov),
    outer_iterations       = 0L,
    converged              = TRUE,
    feasible               = TRUE,
    infeasible_reason      = NA_character_,
    stringsAsFactors       = FALSE
  )
}


# ==============================================================================
# Section 2: Per-Cell Two-Stage Solve
# ==============================================================================

#' @title Solve One K/r/Scenario Cell (Fixed 2100 Check, Then Free-Time Search)
#' @description
#' Stage 1 solves the conventional fixed-time T = 2100 problem and applies
#' the same feasibility test as cdr_scale_sensitivity.R (solver converged AND
#' final_temperature <= 1.5 + feasibility_temp_tolerance). Feasible cells are
#' packaged directly. Stage 2, for infeasible cells only, searches
#' [2100, t_max] with optimal_control_free_terminal_time() for the earliest
#' year the target becomes reachable at all. See the file header for why the
#' search is restricted to years beyond 2100.
#'
#' @param parameter_df Single-row parameter data frame
#' @param emissions_df Emissions data frame extending to t_max
#' @param economic_df Economic data frame extending to t_max
#' @param scenario Scenario name
#' @param cdr_capacity_function CDR capacity function for this K/r cell
#' @param K_val,r_val Grid coordinates (recorded in the result row)
#' @param t_max Latest candidate return year (default: 2200)
#' @param mitigation_delay_years,cdr_delay_years Deployment delays
#' @param use_mitigation_capacity_limit,mitigation_capacity_function
#'   Mitigation capacity settings passed to both stages
#' @param feasibility_temp_tolerance Temperature tolerance for the stage-1
#'   feasibility test in °C (default: 0.05, matching package_scale_result())
#'
#' @return Single-row data frame with the fixed free-time column set
solve_free_time_cell <- function(parameter_df,
                                 emissions_df,
                                 economic_df,
                                 scenario,
                                 cdr_capacity_function,
                                 K_val,
                                 r_val,
                                 t_max                         = 2200,
                                 mitigation_delay_years        = 0,
                                 cdr_delay_years               = 0,
                                 use_mitigation_capacity_limit = FALSE,
                                 mitigation_capacity_function  = NULL,
                                 feasibility_temp_tolerance    = 0.05) {

  # ----------------------------------------------------------------------------
  # Stage 1: conventional fixed-time T = 2100 solve
  # ----------------------------------------------------------------------------

  truncated <- truncate_scenario_data(emissions_df, economic_df, t_max = 2100)

  sol_2100 <- tryCatch(
    optimal_control_shooting(
      parameter_df                  = parameter_df,
      emissions_df                  = truncated$emissions_df,
      economic_df                   = truncated$economic_df,
      scenario                      = scenario,
      mitigation_delay_years        = mitigation_delay_years,
      cdr_delay_years               = cdr_delay_years,
      use_mitigation_capacity_limit = use_mitigation_capacity_limit,
      mitigation_capacity_function  = mitigation_capacity_function,
      use_cdr_capacity_limit        = TRUE,
      cdr_capacity_function         = cdr_capacity_function,
      verbose                       = FALSE
    ),
    error = function(e) NULL
  )

  feasible_2100 <- !is.null(sol_2100) &&
    isTRUE(sol_2100$converged) &&
    is.finite(sol_2100$final_temperature) &&
    sol_2100$final_temperature <= (1.5 + feasibility_temp_tolerance)

  if (feasible_2100) {
    return(package_by2100_result(sol_2100, K_val, r_val, scenario))
  }

  # ----------------------------------------------------------------------------
  # Stage 2: free-terminal-time search over [2100, t_max]
  # ----------------------------------------------------------------------------

  ft <- optimal_control_free_terminal_time(
    parameter_df                  = parameter_df,
    emissions_df                  = emissions_df,
    economic_df                   = economic_df,
    scenario                      = scenario,
    t_min                         = 2100,
    t_max                         = t_max,
    mitigation_delay_years        = mitigation_delay_years,
    cdr_delay_years               = cdr_delay_years,
    use_mitigation_capacity_limit = use_mitigation_capacity_limit,
    mitigation_capacity_function  = mitigation_capacity_function,
    use_cdr_capacity_limit        = TRUE,
    cdr_capacity_function         = cdr_capacity_function,
    verbose                       = FALSE
  )

  package_free_time_result(ft, K_val, r_val, scenario)
}


# ==============================================================================
# Section 3: Runtime Estimator
# ==============================================================================

#' @title Estimate Free-Time Sweep Runtime from a Small Timed Subset
#' @description
#' Per-cell cost is bimodal: cells feasible by 2100 need one fixed-time solve
#' (stage 1 only), while infeasible cells additionally run the free-time
#' search (~10-15 inner shooting solves on horizons up to 2200). Before
#' committing to a full grid, this times a few representative K/r cells
#' serially (one scenario each, spread from weakest to strongest CDR so both
#' stages are sampled) and extrapolates to the full sweep, mirroring how
#' runtime is estimated for the other analysis chunks.
#'
#' @param parameter_df Single-row data frame of model parameters
#' @param emissions_df Emissions data frame extending to t_max
#' @param economic_df Economic data frame extending to t_max
#' @param scenarios Character vector of scenario names for the full sweep
#' @param cdr_grid Data frame of K/r combinations for the full sweep
#' @param n_sample Number of grid cells to time (default: 3, spread across
#'   the grid from weakest to strongest CDR)
#' @param n_cores Cores assumed for the parallel sweep (default: detected - 1)
#' @param ... Further arguments passed to optimal_control_free_terminal_time()
#'   (capacity settings are built per-cell as in the sweep itself)
#' @param g_initial Starting CDR level for the logistic curve (default: 2)
#' @param t_start CDR deployment start year (default: 2025)
#' @param use_mitigation_capacity_limit Logical (default: TRUE)
#' @param mitigation_capacity_function Mitigation capacity function
#'   (default: make_zero_capacity())
#' @param verbose Print per-sample timing (default: TRUE)
#'
#' @return Invisibly, a list with mean_minutes_per_cell and
#'   estimated_total_minutes; prints a summary when verbose
estimate_free_time_sweep_runtime <- function(parameter_df,
                                             emissions_df,
                                             economic_df,
                                             scenarios,
                                             cdr_grid,
                                             n_sample                      = 3,
                                             n_cores                       = NULL,
                                             g_initial                     = 2,
                                             t_start                       = 2025,
                                             use_mitigation_capacity_limit = TRUE,
                                             mitigation_capacity_function  = make_zero_capacity(),
                                             verbose                       = TRUE,
                                             ...) {

  n_cores_actual <- if (is.null(n_cores)) {
    max(1, parallel::detectCores() - 1)
  } else {
    as.integer(n_cores)
  }

  # Sample rows spread across the grid ordering (weak -> strong CDR) so both
  # the cheap stage-1-only cells and the expensive stage-2 cells are sampled
  sample_idx <- unique(round(seq(1, nrow(cdr_grid), length.out = n_sample)))
  times <- numeric(length(sample_idx))

  for (i in seq_along(sample_idx)) {
    row <- cdr_grid[sample_idx[i], ]
    cdr_fn <- make_logistic_from_zero(
      g_initial = g_initial, K = row$K, r = row$r, t_start = t_start
    )
    t0 <- Sys.time()
    cell <- solve_free_time_cell(
      parameter_df                  = parameter_df,
      emissions_df                  = emissions_df,
      economic_df                   = economic_df,
      scenario                      = scenarios[1],
      cdr_capacity_function         = cdr_fn,
      K_val                         = row$K,
      r_val                         = row$r,
      use_mitigation_capacity_limit = use_mitigation_capacity_limit,
      mitigation_capacity_function  = mitigation_capacity_function,
      ...
    )
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    if (verbose) {
      cat(sprintf("  Sample %d/%d (K = %.1f, r = %.4f): %.2f min (%s)\n",
                  i, length(sample_idx), row$K, row$r, times[i],
                  ifelse(is.na(cell$optimum_type), "no return by t_max",
                         cell$optimum_type)))
    }
  }

  mean_minutes <- mean(times)
  total_cells  <- nrow(cdr_grid) * length(scenarios)
  est_total    <- mean_minutes * total_cells / n_cores_actual

  if (verbose) {
    cat(sprintf("\nMean per cell:        %.2f min (one scenario)\n", mean_minutes))
    cat(sprintf("Full sweep:           %d cells (%d combinations x %d scenarios)\n",
                total_cells, nrow(cdr_grid), length(scenarios)))
    cat(sprintf("Estimated total:      %.0f min (~%.1f hours) on %d cores\n",
                est_total, est_total / 60, n_cores_actual))
  }

  invisible(list(
    mean_minutes_per_cell   = mean_minutes,
    sample_minutes          = times,
    estimated_total_minutes = est_total,
    n_cores_assumed         = n_cores_actual
  ))
}


# ==============================================================================
# Section 4: Main Sweep Function
# ==============================================================================

#' @title Run CDR Scale Free-Terminal-Time Sensitivity Analysis
#' @description
#' Sweeps the CDR logistic scale parameters K and r across a grid. Each
#' K/r/scenario cell is solved with the two-stage logic of
#' solve_free_time_cell(): feasible-by-2100 cells are recorded with
#' T_star = 2100 / years_beyond_2100 = 0; infeasible cells get the
#' free-terminal-time search over [2100, t_max] for the earliest achievable
#' return year. Returns a combined data frame with T_star, optimum_type,
#' hamiltonian_at_T_star, years_beyond_2100 and the standard outcome metrics
#' evaluated at T_star.
#'
#' optimum_type values in the output:
#'   - "feasible_by_2100":  cell meets the target under the standard fixed
#'                          T = 2100 solve (stage 1)
#'   - "earliest_reachable": earliest year in (2100, t_max] the target
#'                          becomes reachable at all -- see
#'                          optimal_control_free_time.R for why this is also
#'                          always the cost-optimal year and is guaranteed
#'                          to arrive via a genuine descent from above
#'   - NA:                  no return achievable by t_max (T_star = NA,
#'                          feasible = FALSE) or unhandled solver error (see
#'                          infeasible_reason)
#'
#' Structure mirrors run_cdr_scale_sensitivity(): same grid input (use
#' build_cdr_scale_grid() from cdr_scale_sensitivity.R), same always-produce-
#' a-row philosophy.
#'
#' CHECKPOINTING, PROGRESS, AND LOAD BALANCING (added after a 2500-year,
#' 51x51, 5-scenario sweep ran 30+ hours with zero visibility into progress,
#' an accidental laptop sleep raised the possibility of losing all of it,
#' and process monitoring showed only 2-3 of 13 worker processes were still
#' active hours into the run -- the rest had gone idle long ago): the grid
#' is processed in batches of checkpoint_every combinations, dispatched with
#' clusterApplyLB() rather than one static parLapply() call over the whole
#' grid. Per-cell cost here is wildly non-uniform (a narrow
#' feasibility-boundary band in (K, r) costs 100-300x an easy cell -- see
#' the Aug 2026 boundary probe), and the grid is sorted by K then r, so
#' those expensive cells cluster together in the sort order; a static,
#' contiguous per-worker chunk can land almost entirely on that band while
#' other workers' chunks are all easy cells and finish in seconds.
#' clusterApplyLB dispatches one task at a time and gives an idle worker the
#' next queued task immediately, so no worker can get stuck holding a
#' disproportionate share of the hard cells while the rest sit idle.
#' checkpoint_every defaults to 5x the core count when parallel specifically
#' so each worker has multiple tasks to rebalance across per batch (a batch
#' of exactly n_cores tasks gives the scheduler nothing to rebalance -- one
#' task per worker either way).
#'
#' After each batch, progress (combinations done / total, elapsed time, a
#' rate-based ETA) is printed and, if checkpoint_path is supplied, partial
#' results are saved to that file. If checkpoint_path already exists when
#' the function is called again (e.g. after a crash or a deliberate
#' restart), previously-completed K/r combinations are loaded and skipped --
#' only the remaining grid is processed. Worst-case loss on an interruption
#' is one in-flight batch (checkpoint_every combinations), not the whole
#' sweep. CAVEAT: the
#' checkpoint file does not record or validate scenarios/t_max/grid against
#' the current call -- reusing checkpoint_path across genuinely different
#' analysis configurations will silently mix incompatible results together.
#' Use a distinct checkpoint_path per distinct configuration.
#'
#' @param parameter_df Single-row data frame containing model parameters
#' @param emissions_df Emissions data frame extending to t_max
#'   (interpolate_ssp_emissions(..., end_year = 2200))
#' @param economic_df Economic data frame extending to t_max
#' @param scenarios Character vector of scenario names to compare
#' @param cdr_grid Data frame of K/r combinations from build_cdr_scale_grid()
#' @param g_initial Starting CDR level for the logistic curve (GtCO2/year,
#'   default: 2)
#' @param t_start CDR deployment start year (default: 2025)
#' @param t_max Latest candidate return year (default: 2200). emissions_df /
#'   economic_df must extend at least this far.
#' @param mitigation_delay_years Years to delay mitigation start (default: 0)
#' @param cdr_delay_years Years to delay CDR deployment start (default: 0)
#' @param use_mitigation_capacity_limit Logical: activate mitigation capacity
#'   constraint (default: FALSE). Set TRUE with make_zero_capacity() to
#'   eliminate mitigation entirely.
#' @param mitigation_capacity_function Capacity function for mitigation.
#'   Required when use_mitigation_capacity_limit = TRUE.
#' @param feasibility_temp_tolerance Temperature tolerance (°C) for the
#'   stage-1 feasible-by-2100 test (default: 0.05, matching
#'   package_scale_result() in cdr_scale_sensitivity.R)
#' @param use_parallel Logical: enable parallel processing (default: TRUE)
#' @param n_cores Integer number of cores, or NULL for auto-detection
#' @param checkpoint_path Character file path (e.g.
#'   "output/checkpoint_2500_51x51.rds") to save partial results to after
#'   every batch, and to resume from if it already exists. Default NULL
#'   disables checkpointing entirely (original one-shot behaviour).
#' @param checkpoint_every Integer: number of grid combinations per batch /
#'   checkpoint save. Default NULL uses n_cores_actual when parallel (one
#'   full round of work per core between saves) or 5 when serial.
#' @param progress_log Character file path for a plain-text, human-tailable
#'   progress log (one line per batch: timestamp, combinations done,
#'   elapsed time, ETA). Default NULL prints the same line to console only
#'   (governed by verbose) without writing a file.
#' @param save_results Logical: save results to output/ (default: TRUE)
#' @param verbose Logical: print progress information (default: TRUE)
#' @param output_dir Character output directory (default: "output")
#' @param output_prefix Prefix for output filenames
#'   (default: "cdr_scale_free_time_sensitivity")
#'
#' @return List containing:
#'   - combined_results: one row per K/r/scenario with columns K, r,
#'     scenario, scenario_short, T_star, optimum_type,
#'     hamiltonian_at_T_star, years_beyond_2100, peak_temperature,
#'     final_temperature, years_above_1p5, final_emissions, emission_gap,
#'     total_cost, mitig_cost, remov_cost, temp_cost,
#'     total_mitigation_units, total_cdr_units, outer_iterations, converged,
#'     feasible, infeasible_reason
#'   - summary_stats: counts and feasibility rate
#'   - run_info: sweep metadata including runtime and saved file paths
#'
#' @examples
#' \dontrun{
#' grid <- build_cdr_scale_grid(25, 200, 11, 0.02, 0.20, 11)
#' results <- run_free_time_scale_sensitivity(
#'   parameter_df                  = parameter_df[1, ],
#'   emissions_df                  = emissions_df_ext,
#'   economic_df                   = economic_df_ext,
#'   scenarios                     = c("SSP2-Baseline"),
#'   cdr_grid                      = grid,
#'   use_mitigation_capacity_limit = TRUE,
#'   mitigation_capacity_function  = make_zero_capacity()
#' )
#' }
run_free_time_scale_sensitivity <- function(parameter_df,
                                            emissions_df,
                                            economic_df,
                                            scenarios,
                                            cdr_grid,
                                            g_initial                     = 2,
                                            t_start                       = 2025,
                                            t_max                         = 2200,
                                            mitigation_delay_years        = 0,
                                            cdr_delay_years               = 0,
                                            use_mitigation_capacity_limit = FALSE,
                                            mitigation_capacity_function  = NULL,
                                            feasibility_temp_tolerance    = 0.05,
                                            use_parallel                  = TRUE,
                                            n_cores                       = NULL,
                                            checkpoint_path               = NULL,
                                            checkpoint_every              = NULL,
                                            progress_log                  = NULL,
                                            save_results                  = TRUE,
                                            verbose                       = TRUE,
                                            output_dir                    = "output",
                                            output_prefix                 = "cdr_scale_free_time_sensitivity") {

  # ============================================================================
  # Input validation
  # ============================================================================

  if (!is.data.frame(parameter_df)) {
    stop("parameter_df must be a data frame (received ", class(parameter_df)[1], ")")
  }
  if (nrow(parameter_df) != 1) {
    stop("parameter_df must contain exactly one row (received ",
         nrow(parameter_df), " rows)")
  }

  if (!is.data.frame(cdr_grid) || !all(c("K", "r") %in% names(cdr_grid))) {
    stop("cdr_grid must be a data frame with columns 'K' and 'r'. ",
         "Use build_cdr_scale_grid() to create it.")
  }
  if (nrow(cdr_grid) < 1) {
    stop("cdr_grid must contain at least one row")
  }

  if (length(scenarios) == 0) {
    stop("At least one scenario must be specified")
  }

  available_scenarios <- unique(emissions_df$Scenario)
  missing_scenarios   <- setdiff(scenarios, available_scenarios)
  if (length(missing_scenarios) > 0) {
    warning("The following scenarios are not available in emissions_df and will be skipped:\n",
            paste(" -", missing_scenarios, collapse = "\n"))
    scenarios <- intersect(scenarios, available_scenarios)
    if (length(scenarios) == 0) {
      stop("No valid scenarios remain after filtering. Available: ",
           paste(available_scenarios, collapse = ", "))
    }
  }

  if (max(emissions_df$Year) < t_max || max(economic_df$Year) < t_max) {
    stop("emissions_df/economic_df end before t_max = ", t_max,
         ". Rebuild them with interpolate_ssp_emissions()/",
         "interpolate_ssp_economic() using end_year = ", t_max, ".")
  }

  if (use_mitigation_capacity_limit && is.null(mitigation_capacity_function)) {
    stop("mitigation_capacity_function must be supplied when ",
         "use_mitigation_capacity_limit = TRUE")
  }
  if (use_mitigation_capacity_limit && !is.function(mitigation_capacity_function)) {
    stop("mitigation_capacity_function must be a function")
  }

  if (!is.numeric(g_initial) || g_initial < 0) {
    stop("g_initial must be a non-negative number")
  }
  if (!is.numeric(t_start)) {
    stop("t_start must be numeric")
  }

  # ============================================================================
  # Setup and initialisation
  # ============================================================================

  n_combinations     <- nrow(cdr_grid)
  n_scenarios        <- length(scenarios)
  overall_start_time <- Sys.time()

  n_cores_actual <- if (is.null(n_cores)) {
    max(1, parallel::detectCores() - 1)
  } else {
    as.integer(n_cores)
  }

  if (is.null(checkpoint_every)) {
    # A batch of exactly n_cores tasks gives clusterApplyLB nothing to
    # rebalance -- each worker gets exactly one task regardless of
    # scheduler. Multiple tasks per worker per batch (5x here) is what lets
    # a worker that finishes an easy cell immediately pick up more queued
    # work instead of sitting idle until the whole batch closes.
    checkpoint_every <- if (use_parallel) n_cores_actual * 5 else 5
  }
  checkpoint_every <- max(1L, as.integer(checkpoint_every))

  if (verbose) {
    cat("=== CDR SCALE FREE-TERMINAL-TIME SENSITIVITY ANALYSIS ===\n")
    cat("Grid size:            ", n_combinations, "combinations (",
        length(unique(cdr_grid$K)), "K x",
        length(unique(cdr_grid$r)), "r )\n")
    cat("K range:              ",
        sprintf("%.1f - %.1f GtCO2/year", min(cdr_grid$K), max(cdr_grid$K)), "\n")
    cat("r range:              ",
        sprintf("%.4f - %.4f", min(cdr_grid$r), max(cdr_grid$r)), "\n")
    cat("Scenarios:            ", paste(scenarios, collapse = ", "), "\n")
    cat("Free-T search range:  ( 2100,", t_max, "] for cells infeasible at 2100\n")
    cat("g_initial:            ", g_initial, "\n")
    cat("t_start:              ", t_start, "\n")
    cat("Mitigation zeroed:    ", use_mitigation_capacity_limit, "\n")
    cat("Parallel processing:  ", use_parallel, "\n")
    cat("Checkpoint file:      ", if (is.null(checkpoint_path)) "(none)" else checkpoint_path, "\n")
    cat("Checkpoint batch size:", checkpoint_every, "combinations\n")
    cat("Total cells:          ", n_combinations * n_scenarios,
        "(1 fixed solve each + ~10-15 inner solves where infeasible at 2100)\n\n")
  }

  # ============================================================================
  # Checkpoint resume: load prior partial results if present
  # ============================================================================
  # A checkpoint holds the raw per-combination result data frames already
  # computed (checkpoint_results) and a K/r lookup of which grid rows they
  # cover (completed_combinations), used to skip those rows this run.

  checkpoint_results     <- list()
  completed_combinations <- data.frame(K = numeric(0), r = numeric(0))

  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    if (verbose) cat("Found existing checkpoint:", checkpoint_path, "-- resuming\n")
    prior <- readRDS(checkpoint_path)
    checkpoint_results     <- prior$checkpoint_results
    completed_combinations <- prior$completed_combinations
    if (verbose) {
      cat(sprintf("  %d of %d combinations already completed; skipping those\n\n",
                  nrow(completed_combinations), n_combinations))
    }
  }

  remaining_grid <- dplyr::anti_join(cdr_grid, completed_combinations, by = c("K", "r"))
  rownames(remaining_grid) <- NULL
  n_remaining <- nrow(remaining_grid)

  # ============================================================================
  # Define per-combination worker function
  # ============================================================================
  # Processes one row of cdr_grid (one K/r combination) by running the
  # two-stage cell solve serially across all scenarios. Returns a data frame
  # with one row per scenario.

  run_one_combination <- function(K_val, r_val) {

    # Capacity curve starts at t_start + cdr_delay_years so no capacity
    # builds during the delay period (same convention as the fixed-time sweep)
    cdr_fn <- make_logistic_from_zero(
      g_initial = g_initial,
      K         = K_val,
      r         = r_val,
      t_start   = t_start + cdr_delay_years
    )

    scenario_rows <- lapply(scenarios, function(sc) {
      tryCatch({
        solve_free_time_cell(
          parameter_df                  = parameter_df,
          emissions_df                  = emissions_df,
          economic_df                   = economic_df,
          scenario                      = sc,
          cdr_capacity_function         = cdr_fn,
          K_val                         = K_val,
          r_val                         = r_val,
          t_max                         = t_max,
          mitigation_delay_years        = mitigation_delay_years,
          cdr_delay_years               = cdr_delay_years,
          use_mitigation_capacity_limit = use_mitigation_capacity_limit,
          mitigation_capacity_function  = mitigation_capacity_function,
          feasibility_temp_tolerance    = feasibility_temp_tolerance
        )
      }, error = function(e) {
        package_free_time_error_result(K_val, r_val, sc,
                                       error_message = conditionMessage(e))
      })
    })

    dplyr::bind_rows(scenario_rows)
  }

  # ============================================================================
  # Execute the remaining grid in checkpointed batches
  # ============================================================================
  # Unlike a single parLapply() over the whole grid (which returns nothing
  # until every combination is done), the grid is processed in
  # checkpoint_every-sized batches with one shared cluster reused across
  # batches -- periodic progress output and a save point after every batch,
  # with no per-batch cluster startup cost.

  if (n_remaining == 0) {

    if (verbose) cat("All combinations already completed in checkpoint -- nothing to run\n")

  } else {

    cl <- NULL
    if (use_parallel && n_remaining > 1) {
      if (verbose) cat("Starting parallel sweep with", n_cores_actual, "cores...\n")
      cl <- tryCatch({
        cl0 <- parallel::makeCluster(n_cores_actual)
        parallel::clusterEvalQ(cl0, { library(dplyr) })
        parallel::clusterExport(
          cl0,
          varlist = c(
            # Functions defined in this file
            "run_one_combination",
            "solve_free_time_cell",
            "package_free_time_result",
            "package_by2100_result",
            "package_free_time_error_result",
            # Free-time solver chain (must be sourced before calling this)
            "optimal_control_free_terminal_time",
            "compute_terminal_hamiltonian",
            "truncate_scenario_data",
            "optimal_control_shooting",
            "optimal_control_solve",
            # Capacity helper functions
            "make_logistic_from_zero",
            "make_zero_capacity",
            "make_exponential_capacity",
            "make_logistic_capacity",
            "make_linear_capacity",
            "make_piecewise_capacity",
            "make_power_capacity",
            # Data objects
            "parameter_df",
            "emissions_df",
            "economic_df",
            "scenarios",
            # Scalar parameters
            "g_initial",
            "t_start",
            "t_max",
            "mitigation_delay_years",
            "cdr_delay_years",
            "use_mitigation_capacity_limit",
            "mitigation_capacity_function",
            "feasibility_temp_tolerance"
          ),
          envir = environment()
        )
        cl0
      }, error = function(e) {
        if (verbose) {
          cat("Parallel cluster setup failed:", e$message, "\n")
          cat("Falling back to serial processing...\n")
        }
        NULL
      })
      if (is.null(cl)) use_parallel <- FALSE
    }

    if (verbose && !use_parallel) {
      cat("Running serial sweep across", n_remaining, "remaining combinations...\n")
    }

    batch_starts    <- seq(1L, n_remaining, by = checkpoint_every)
    sweep_start     <- Sys.time()
    n_done_this_run <- 0L

    for (batch_start in batch_starts) {

      batch_end  <- min(batch_start + checkpoint_every - 1L, n_remaining)
      batch_rows <- remaining_grid[batch_start:batch_end, , drop = FALSE]

      if (use_parallel && !is.null(cl)) {
        parallel::clusterExport(cl, "batch_rows", envir = environment())
        # clusterApplyLB(), not parLapply(): parLapply splits the batch into
        # one static, contiguous chunk per worker up front and never
        # rebalances. Per-cell cost here is wildly non-uniform (a narrow
        # feasibility-boundary band costs 100-300x an easy cell -- see the
        # Aug 2026 boundary probe), and the grid is sorted by K then r, so
        # boundary cells cluster together in the sort order. A worker whose
        # static chunk happens to cover that band gets stuck for hours while
        # workers with easy chunks finish in seconds and sit permanently
        # idle -- this was directly observed on the live 2500-year sweep
        # (the same 2-3 of 13 processes showed CPU activity across every
        # check spanning 4+ hours; the rest were long since idle).
        # clusterApplyLB dispatches one task at a time and gives an idle
        # worker the next queued task immediately, so no single worker can
        # get stuck with a disproportionate share of the hard cells. Slightly
        # more communication overhead per task than parLapply, which is
        # irrelevant here given individual task costs range from ~1 sec to
        # several minutes.
        batch_results <- parallel::clusterApplyLB(
          cl, seq_len(nrow(batch_rows)),
          function(i) run_one_combination(batch_rows$K[i], batch_rows$r[i])
        )
      } else {
        batch_results <- lapply(
          seq_len(nrow(batch_rows)),
          function(i) run_one_combination(batch_rows$K[i], batch_rows$r[i])
        )
      }

      checkpoint_results     <- c(checkpoint_results, batch_results)
      completed_combinations <- dplyr::bind_rows(completed_combinations,
                                                 batch_rows[, c("K", "r")])
      n_done_this_run <- n_done_this_run + nrow(batch_rows)

      elapsed_min <- as.numeric(difftime(Sys.time(), sweep_start, units = "mins"))
      rate        <- if (elapsed_min > 0) n_done_this_run / elapsed_min else NA_real_
      eta_min     <- if (!is.na(rate) && rate > 0) {
        (n_remaining - n_done_this_run) / rate
      } else NA_real_

      progress_msg <- sprintf(
        "[%s] %d/%d this run (%d/%d overall) | %.1f min elapsed this run | ETA ~%s\n",
        format(Sys.time(), "%H:%M:%S"), n_done_this_run, n_remaining,
        nrow(completed_combinations), n_combinations, elapsed_min,
        if (is.na(eta_min)) "unknown" else sprintf("%.0f min (%.1f hr)", eta_min, eta_min / 60)
      )

      if (verbose) cat(progress_msg)
      if (!is.null(progress_log)) cat(progress_msg, file = progress_log, append = TRUE)

      if (!is.null(checkpoint_path)) {
        checkpoint_dir <- dirname(checkpoint_path)
        if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
        saveRDS(
          list(checkpoint_results = checkpoint_results,
               completed_combinations = completed_combinations),
          checkpoint_path
        )
      }
    }

    if (!is.null(cl)) parallel::stopCluster(cl)
  }

  # ============================================================================
  # Combine results
  # ============================================================================

  if (verbose) cat("\nCombining results...\n")

  results_list <- checkpoint_results[!sapply(checkpoint_results, is.null)]

  if (length(results_list) == 0) {
    stop("No valid results were produced. Check that the solver functions are ",
         "correctly sourced and that the parameter grid is valid.")
  }

  combined_results <- dplyr::bind_rows(results_list)
  rownames(combined_results) <- NULL

  combined_results <- combined_results[
    order(combined_results$K, combined_results$r, combined_results$scenario), ]
  rownames(combined_results) <- NULL

  # ============================================================================
  # Summary statistics
  # ============================================================================

  n_total        <- nrow(combined_results)
  n_feasible     <- sum(combined_results$feasible,  na.rm = TRUE)
  n_converged    <- sum(combined_results$converged, na.rm = TRUE)
  n_failed       <- n_total - n_converged
  expected_total <- n_combinations * n_scenarios
  n_by_2100      <- sum(combined_results$feasible &
                          combined_results$T_star <= 2100, na.rm = TRUE)

  if (n_total < expected_total && verbose) {
    cat(sprintf(
      "WARNING: Expected %d rows (%d combinations x %d scenarios) but got %d.\n",
      expected_total, n_combinations, n_scenarios, n_total))
    cat("Some combinations may have failed entirely. Check error output above.\n")
  }

  total_runtime <- difftime(Sys.time(), overall_start_time, units = "mins")

  summary_stats <- list(
    n_combinations           = n_combinations,
    n_scenarios              = n_scenarios,
    expected_total_rows      = expected_total,
    n_total_rows             = n_total,
    n_converged              = n_converged,
    n_failed                 = n_failed,
    n_feasible               = n_feasible,
    n_feasible_by_2100       = n_by_2100,
    overall_feasibility_rate = if (n_total > 0) n_feasible / n_total else NA_real_
  )

  # ============================================================================
  # Report summary
  # ============================================================================

  if (verbose) {
    cat("\n=== CDR SCALE FREE-TERMINAL-TIME ANALYSIS COMPLETE ===\n")
    cat("Total runtime:        ",
        sprintf("%.1f", as.numeric(total_runtime)), "minutes\n")
    cat("Expected rows:        ", expected_total, "\n")
    cat("Total rows:           ", n_total, "\n")
    cat("Converged:            ", n_converged, "/", n_total, "\n")
    cat("Feasible return:      ", n_feasible, "/", n_total, "\n")
    cat("  ...by 2100:         ", n_by_2100, "\n")
    cat("  ...only after 2100: ", n_feasible - n_by_2100, "\n")
  }

  # ============================================================================
  # Assemble results list
  # ============================================================================

  run_info <- list(
    K_range                       = range(cdr_grid$K),
    r_range                       = range(cdr_grid$r),
    n_K                           = length(unique(cdr_grid$K)),
    n_r                           = length(unique(cdr_grid$r)),
    n_combinations                = n_combinations,
    g_initial                     = g_initial,
    t_start                       = t_start,
    t_min                         = 2100,
    t_max                         = t_max,
    scenarios                     = scenarios,
    mitigation_delay_years        = mitigation_delay_years,
    cdr_delay_years               = cdr_delay_years,
    use_mitigation_capacity_limit = use_mitigation_capacity_limit,
    use_parallel                  = use_parallel,
    checkpoint_path               = checkpoint_path,
    checkpoint_every              = checkpoint_every,
    start_time                    = overall_start_time,
    end_time                      = Sys.time(),
    total_runtime_minutes         = as.numeric(total_runtime)
  )

  results_list_out <- list(
    combined_results = combined_results,
    summary_stats    = summary_stats,
    run_info         = run_info
  )

  # ============================================================================
  # Save results if requested
  # ============================================================================

  if (save_results) {

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    timestamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
    rds_filename <- paste0(output_prefix, "_", timestamp, ".rds")
    csv_filename <- paste0(output_prefix, "_", timestamp, ".csv")
    rds_filepath <- file.path(output_dir, rds_filename)
    csv_filepath <- file.path(output_dir, csv_filename)

    results_list_out$run_info$saved_files <- list(
      rds = rds_filepath,
      csv = csv_filepath
    )

    saveRDS(results_list_out, rds_filepath)
    write.csv(combined_results, csv_filepath, row.names = FALSE)

    if (verbose) {
      cat("\n=== RESULTS SAVED ===\n")
      cat("RDS file:  ", rds_filepath, "\n")
      cat("CSV file:  ", csv_filepath, "\n")
    }
  }

  return(results_list_out)

} # Close run_free_time_scale_sensitivity
