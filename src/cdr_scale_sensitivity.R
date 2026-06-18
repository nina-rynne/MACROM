# ==============================================================================
# CDR Scale Sensitivity Analysis Functions
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
# Sensitivity analysis sweeping the CDR logistic scale function parameters K
# (carrying capacity, GtCO2/year) and r (growth rate) across a user-defined
# grid. For each K/r combination the full scenario comparison is run across all
# requested SSP scenarios, and summary metrics are extracted for heatmap
# visualisation. Designed to be used with cdr_scale_sensitivity_visualisation.R.
#
# Feasibility philosophy (matching delayed_deployment.R):
# Every K/r/scenario combination always produces a result row. The solver saves
# its best solution even when the 1.5°C target is unachievable. The feasible
# flag is set after the fact based on whether final_temperature <= 1.5°C
# (within a small numerical tolerance), not on whether the solver converged.
# Infeasible cells are plotted with their actual outcome values and optionally
# marked in the visualisation, exactly as in the delayed deployment heatmaps.
#
# Parallel strategy (matching delayed_deployment.R):
# Uses parLapply() to distribute K/r combinations across workers. Each worker
# runs run_one_combination() for one grid row and returns a data frame.
# parLapply() always returns a clean flat list with one element per input,
# avoiding the nested list / .combine complexity of foreach.
# 
# Version: 1.3.0
# Last updated: March 2026
# ==============================================================================


# ==============================================================================
# Section 1: Grid Builder
# ==============================================================================

#' @title Build CDR Scale Parameter Grid
#' @description
#' Constructs a data frame of all K x r combinations to sweep across. K values
#' are always linearly spaced. r values are linearly spaced by default, or
#' log-spaced when r_log_scale = TRUE (useful when r spans more than one order
#' of magnitude). n_K and n_r are independent so non-square grids are supported.
#'
#' @param K_min Numeric minimum CDR carrying capacity in GtCO2/year (e.g. 20)
#' @param K_max Numeric maximum CDR carrying capacity in GtCO2/year (e.g. 200)
#' @param n_K Integer number of K values to generate (e.g. 40)
#' @param r_min Numeric minimum CDR growth rate (e.g. 0.01)
#' @param r_max Numeric maximum CDR growth rate (e.g. 0.3)
#' @param n_r Integer number of r values to generate (e.g. 10)
#' @param r_log_scale Logical. If TRUE, r values are log-spaced; if FALSE
#'   (default) they are linearly spaced.
#'
#' @return Data frame with n_K * n_r rows and two columns:
#'   - K: CDR carrying capacity (GtCO2/year)
#'   - r: CDR growth rate
#'
#' @examples
#' # Non-square grid: 40 K values x 10 r values = 400 combinations
#' grid <- build_cdr_scale_grid(20, 200, 40, 0.01, 0.3, 10)
#'
#' # Log-spaced r values for wide r range
#' grid <- build_cdr_scale_grid(20, 200, 10, 0.001, 0.3, 10, r_log_scale = TRUE)
build_cdr_scale_grid <- function(K_min,
                                 K_max,
                                 n_K,
                                 r_min,
                                 r_max,
                                 n_r,
                                 r_log_scale = FALSE) {
  
  # ============================================================================
  # Input validation
  # ============================================================================
  
  if (!is.numeric(K_min) || !is.numeric(K_max) || K_min <= 0 || K_max <= K_min) {
    stop("K_min must be positive and K_max must be greater than K_min")
  }
  if (!is.numeric(r_min) || !is.numeric(r_max) || r_min <= 0 || r_max <= r_min) {
    stop("r_min must be positive and r_max must be greater than r_min")
  }
  if (!is.numeric(n_K) || n_K < 2) {
    stop("n_K must be an integer >= 2")
  }
  if (!is.numeric(n_r) || n_r < 2) {
    stop("n_r must be an integer >= 2")
  }
  
  # ============================================================================
  # Generate sequences
  # ============================================================================
  
  K_values <- seq(K_min, K_max, length.out = as.integer(n_K))
  
  if (r_log_scale) {
    r_values <- exp(seq(log(r_min), log(r_max), length.out = as.integer(n_r)))
  } else {
    r_values <- seq(r_min, r_max, length.out = as.integer(n_r))
  }
  
  # ============================================================================
  # Build full combination grid
  # ============================================================================
  
  grid_df <- expand.grid(K = K_values, r = r_values)
  grid_df <- grid_df[order(grid_df$K, grid_df$r), ]
  rownames(grid_df) <- NULL
  
  return(grid_df)
}


# ==============================================================================
# Section 2: Results Packagers
# ==============================================================================

#' @title Package CDR Scale Sensitivity Result
#' @description
#' Converts the comparison_summary data frame returned by run_scenario_comparison()
#' into a standardised flat data frame for combining across the grid sweep.
#' Appends K, r, and scenario_short columns. One row is returned per scenario.
#'
#' Feasibility: feasible = TRUE only when the solver converged AND
#' final_temperature <= 1.5 + FEASIBILITY_TEMP_TOLERANCE. This correctly
#' identifies combinations where CDR scale was insufficient to return to the
#' climate target, even when the solver converged to its best-effort solution.
#' Matches the delayed deployment approach where feasibility is assessed on
#' outcomes, not solver convergence status.
#'
#' Column set: always outputs the same fixed columns regardless of what is
#' present in comparison_summary. Missing columns are filled with NA to
#' guarantee that bind_rows() can safely combine all result data frames.
#'
#' @param comparison_summary Data frame from run_scenario_comparison()$comparison_summary
#' @param K_val Numeric CDR carrying capacity used for this run (GtCO2/year)
#' @param r_val Numeric CDR growth rate used for this run
#'
#' @return Data frame with one row per scenario and a fixed complete column set
package_scale_result <- function(comparison_summary,
                                 K_val,
                                 r_val) {
  
  # Temperature tolerance for feasibility: solutions within this many degrees
  # of 1.5°C are considered to have met the target (allows for small numerical
  # imprecision in the shooting method convergence)
  FEASIBILITY_TEMP_TOLERANCE <- 0.05  # degrees C
  
  result_df <- comparison_summary
  
  # Add grid parameter and identifier columns
  result_df$K            <- K_val
  result_df$r            <- r_val
  result_df$scenario_short <- gsub("SSP([0-9])-Baseline", "SSP\\1",
                                   result_df$scenario)
  
  # Feasibility based on final_temperature outcome, not solver convergence.
  # converged is a necessary but not sufficient condition: the solver may
  # converge perfectly to a solution that simply cannot reach 1.5°C because
  # the CDR scale parameters are too small.
  result_df$feasible <- result_df$converged &
    !is.na(result_df$final_temperature) &
    result_df$final_temperature <= (1.5 + FEASIBILITY_TEMP_TOLERANCE)
  
  # Fixed canonical column order — must exactly match package_scale_error_result()
  # Any columns missing from comparison_summary are added as NA
  all_metric_cols <- c(
    "peak_temperature",
    "final_temperature",
    "years_above_1p5",
    "final_emissions",
    "emission_gap",
    "total_cost",
    "mitig_cost",
    "remov_cost",
    "temp_cost",
    "total_mitigation_units",
    "total_cdr_units",
    "converged",
    "feasible"
  )
  
  for (col in all_metric_cols) {
    if (!col %in% names(result_df)) result_df[[col]] <- NA
  }
  
  keep_cols <- c("K", "r", "scenario", "scenario_short", all_metric_cols)
  return(result_df[, keep_cols, drop = FALSE])
}


#' @title Package Error Result for Failed K/r Combination
#' @description
#' Creates a standardised NA-filled result data frame when a K/r combination
#' fails entirely due to an unhandled solver error. Column set is identical to
#' package_scale_result() so that bind_rows() can always combine them cleanly.
#'
#' Note: this is only called on unhandled errors. When the solver simply cannot
#' reach 1.5°C, run_scenario_comparison() still returns a best-effort result
#' and package_scale_result() is used with feasible = FALSE.
#'
#' @param K_val Numeric CDR carrying capacity
#' @param r_val Numeric CDR growth rate
#' @param scenarios Character vector of scenario names that were attempted
#'
#' @return Data frame with one row per scenario, all metric columns NA
package_scale_error_result <- function(K_val, r_val, scenarios) {
  
  data.frame(
    K                      = K_val,
    r                      = r_val,
    scenario               = scenarios,
    scenario_short         = gsub("SSP([0-9])-Baseline", "SSP\\1", scenarios),
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
    converged              = FALSE,
    feasible               = FALSE,
    stringsAsFactors       = FALSE
  )
}


# ==============================================================================
# Section 3: Main Sweep Function
# ==============================================================================

#' @title Run CDR Scale Sensitivity Analysis
#' @description
#' Sweeps the CDR logistic scale function parameters K (carrying capacity) and r
#' (growth rate) across a user-defined grid and runs the full scenario comparison
#' for each combination. Returns a combined data frame of summary metrics suitable
#' for heatmap visualisation with cdr_scale_sensitivity_visualisation.R.
#'
#' Every K/r/scenario combination always produces a result row. When the CDR
#' scale is insufficient to return temperature to 1.5°C, the solver returns its
#' best-effort solution and feasible = FALSE. This matches the delayed deployment
#' analysis philosophy.
#'
#' Parallel strategy: uses parLapply() to distribute K/r combinations across
#' workers. Each worker calls run_one_combination() for one grid row and returns
#' a data frame. parLapply() always returns a clean flat list with one element
#' per input — there is no foreach .combine complexity and no nested lists.
#' This directly mirrors how the delayed deployment analysis distributes work.
#'
#' @param parameter_df Single-row data frame containing model parameters
#' @param emissions_df Data frame with emissions data for multiple scenarios
#' @param economic_df Data frame with economic data for multiple scenarios
#' @param scenarios Character vector of scenario names to compare
#'   (e.g. c("SSP1-Baseline", ..., "SSP5-Baseline"))
#' @param cdr_grid Data frame of K/r combinations from build_cdr_scale_grid().
#'   Must contain columns K and r. Supports non-square grids (n_K != n_r).
#' @param g_initial Numeric starting CDR level for the logistic function
#'   (GtCO2/year). Passed to make_logistic_from_zero(). Default is 1.
#' @param t_start Numeric CDR deployment start year. Default is 2025.
#' @param mitigation_delay_years Years to delay mitigation start (default: 0)
#' @param cdr_delay_years Years to delay CDR deployment start (default: 0)
#' @param use_mitigation_capacity_limit Logical: activate mitigation capacity
#'   constraint (default: FALSE). Set TRUE with make_zero_capacity() to
#'   eliminate mitigation entirely.
#' @param mitigation_capacity_function Capacity function for mitigation with
#'   signature function(year). Required when use_mitigation_capacity_limit = TRUE.
#' @param use_parallel Logical: enable parallel processing (default: TRUE)
#' @param n_cores Integer number of cores, or NULL for auto-detection (default: NULL)
#' @param save_results Logical: save results to output/ directory (default: TRUE)
#' @param verbose Logical: print progress information (default: TRUE)
#' @param output_dir Character string for output directory (default: "output")
#' @param output_prefix Prefix for output filenames (default: "cdr_scale_sensitivity")
#'
#' @return List containing:
#'   - combined_results: Data frame with one row per K/r/scenario combination.
#'     Columns: K, r, scenario, scenario_short, peak_temperature,
#'     final_temperature, years_above_1p5, final_emissions, emission_gap,
#'     total_cost, mitig_cost, remov_cost, temp_cost,
#'     total_mitigation_units, total_cdr_units, converged, feasible.
#'   - summary_stats: List with n_combinations, n_scenarios, expected_total_rows,
#'     n_total_rows, n_converged, n_failed, n_feasible, overall_feasibility_rate
#'   - run_info: Metadata including K/r ranges, g_initial, t_start, scenarios,
#'     runtime, and saved file paths
#'
#' @examples
#' \dontrun{
#' grid <- build_cdr_scale_grid(20, 200, 40, 0.01, 0.3, 10)
#'
#' results <- run_cdr_scale_sensitivity(
#'   parameter_df                  = parameter_df[1, ],
#'   emissions_df                  = emissions_df,
#'   economic_df                   = economic_df,
#'   scenarios                     = c("SSP1-Baseline", "SSP2-Baseline",
#'                                     "SSP3-Baseline", "SSP4-Baseline",
#'                                     "SSP5-Baseline"),
#'   cdr_grid                      = grid,
#'   g_initial                     = 1,
#'   t_start                       = 2025,
#'   use_mitigation_capacity_limit = TRUE,
#'   mitigation_capacity_function  = make_zero_capacity(),
#'   use_parallel                  = TRUE,
#'   save_results                  = TRUE,
#'   verbose                       = TRUE
#' )
#' }
run_cdr_scale_sensitivity <- function(parameter_df,
                                      emissions_df,
                                      economic_df,
                                      scenarios,
                                      cdr_grid,
                                      g_initial                     = 1,
                                      t_start                       = 2025,
                                      mitigation_delay_years        = 0,
                                      cdr_delay_years               = 0,
                                      use_mitigation_capacity_limit = FALSE,
                                      mitigation_capacity_function  = NULL,
                                      use_parallel                  = TRUE,
                                      n_cores                       = NULL,
                                      save_results                  = TRUE,
                                      verbose                       = TRUE,
                                      output_dir                    = "output",
                                      output_prefix                 = "cdr_scale_sensitivity") {
  
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
  
  if (verbose) {
    cat("=== CDR SCALE SENSITIVITY ANALYSIS ===\n")
    cat("Grid size:            ", n_combinations, "combinations (",
        length(unique(cdr_grid$K)), "K x",
        length(unique(cdr_grid$r)), "r )\n")
    cat("K range:              ",
        sprintf("%.1f - %.1f GtCO2/year", min(cdr_grid$K), max(cdr_grid$K)), "\n")
    cat("r range:              ",
        sprintf("%.4f - %.4f", min(cdr_grid$r), max(cdr_grid$r)), "\n")
    cat("Scenarios:            ", paste(scenarios, collapse = ", "), "\n")
    cat("g_initial:            ", g_initial, "\n")
    cat("t_start:              ", t_start, "\n")
    cat("Mitigation zeroed:    ", use_mitigation_capacity_limit, "\n")
    cat("Parallel processing:  ", use_parallel, "\n")
    cat("Total solver runs:    ", n_combinations * n_scenarios, "\n")
    cat("Expected result rows: ", n_combinations * n_scenarios,
        "(one per K/r/scenario)\n\n")
  }
  
  # ============================================================================
  # Define per-combination worker function
  # ============================================================================
  # Processes one row of cdr_grid (one K/r combination) by running
  # run_scenario_comparison() serially across all scenarios. Each call returns
  # a data frame with one row per scenario. Defined here so it captures the
  # enclosing environment and can be exported to cluster workers cleanly.
  
  run_one_combination <- function(K_val, r_val) {
    
    # Build CDR capacity function for this K/r combination
    cdr_fn <- make_logistic_from_zero(
      g_initial = g_initial,
      K         = K_val,
      r         = r_val,
      t_start   = t_start
    )
    
    # Run scenario comparison serially across scenarios
    tryCatch({
      
      sc_result <- run_scenario_comparison(
        parameter_df                  = parameter_df,
        emissions_df                  = emissions_df,
        economic_df                   = economic_df,
        scenarios                     = scenarios,
        mitigation_delay_years        = mitigation_delay_years,
        cdr_delay_years               = cdr_delay_years,
        use_mitigation_capacity_limit = use_mitigation_capacity_limit,
        mitigation_capacity_function  = mitigation_capacity_function,
        use_cdr_capacity_limit        = TRUE,
        cdr_capacity_function         = cdr_fn,
        use_parallel                  = FALSE,  # No nested parallelism
        save_results                  = FALSE,  # Central save at the end
        verbose                       = FALSE
      )
      
      # run_scenario_comparison() always returns its best-effort result even
      # when the temperature target cannot be met, so comparison_summary should
      # always be non-empty. Fall back to error rows only if something went
      # genuinely wrong (NULL or empty summary).
      if (!is.null(sc_result$comparison_summary) &&
          nrow(sc_result$comparison_summary) > 0) {
        package_scale_result(
          comparison_summary = sc_result$comparison_summary,
          K_val              = K_val,
          r_val              = r_val
        )
      } else {
        package_scale_error_result(K_val, r_val, scenarios)
      }
      
    }, error = function(e) {
      # Unhandled solver error: return NA rows so this combination still
      # contributes rows to the output (feasible = FALSE, all metrics NA)
      package_scale_error_result(K_val, r_val, scenarios)
    })
  }
  
  # ============================================================================
  # Execute grid sweep: parallel or serial
  # ============================================================================
  # Pattern mirrors delayed_deployment.R: distribute work across a simple
  # worker pool, collect a flat list of data frames, bind at the end.
  # parLapply() is used for parallel execution because it always returns a
  # clean flat list with exactly one element per input index — there is no
  # foreach .combine logic and no nested list structure to untangle.
  
  if (use_parallel && n_combinations > 1) {
    
    n_cores_actual <- if (is.null(n_cores)) {
      max(1, parallel::detectCores() - 1)
    } else {
      as.integer(n_cores)
    }
    
    if (verbose) {
      cat("Starting parallel sweep with", n_cores_actual, "cores...\n")
    }
    
    # Attempt parallel execution; fall back to serial on any failure
    parallel_success <- tryCatch({
      
      cl <- parallel::makeCluster(n_cores_actual)
      
      # Load dplyr on each worker
      parallel::clusterEvalQ(cl, { library(dplyr) })
      
      # Export all objects each worker needs
      parallel::clusterExport(
        cl,
        varlist = c(
          # Functions defined in this file
          "run_one_combination",
          "package_scale_result",
          "package_scale_error_result",
          # Solver chain (must be sourced in the main session before calling this)
          "run_scenario_comparison",
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
          "mitigation_delay_years",
          "cdr_delay_years",
          "use_mitigation_capacity_limit",
          "mitigation_capacity_function"
        ),
        envir = environment()
      )
      
      # parLapply distributes one combination per worker call and returns a
      # flat list — one data frame per element, no nesting, no batching.
      # cdr_grid must be exported separately because it is accessed inside
      # the anonymous function via the index i.
      parallel::clusterExport(cl, "cdr_grid", envir = environment())
      
      results_list <- parallel::parLapply(
        cl  = cl,
        X   = seq_len(n_combinations),
        fun = function(i) {
          run_one_combination(cdr_grid$K[i], cdr_grid$r[i])
        }
      )
      
      parallel::stopCluster(cl)
      TRUE
      
    }, error = function(e) {
      if (verbose) {
        cat("Parallel execution failed:", e$message, "\n")
        cat("Falling back to serial processing...\n")
      }
      # Ensure cluster is stopped even on failure
      if (exists("cl") && !is.null(cl)) {
        tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
      }
      FALSE
    })
    
    if (!parallel_success) {
      use_parallel <- FALSE
    }
  }
  
  # Serial execution: used when use_parallel = FALSE, n_combinations = 1,
  # or parallel failed and fell back. Loops over grid rows one at a time,
  # calls run_one_combination(), and stores each result in a list.
  # This is the simplest possible approach and is guaranteed to work correctly.
  if (!use_parallel || n_combinations == 1) {
    
    if (verbose) {
      cat("Running serial sweep across", n_combinations, "combinations...\n")
    }
    
    results_list <- vector("list", n_combinations)
    
    for (i in seq_len(n_combinations)) {
      K_i <- cdr_grid$K[i]
      r_i <- cdr_grid$r[i]
      
      if (verbose) {
        elapsed <- round(
          as.numeric(difftime(Sys.time(), overall_start_time, units = "mins")), 1)
        cat(sprintf("  Combination %d / %d: K = %.1f, r = %.4f  (%.1f min elapsed)\n",
                    i, n_combinations, K_i, r_i, elapsed))
      }
      
      results_list[[i]] <- run_one_combination(K_i, r_i)
    }
  }
  
  # ============================================================================
  # Combine results
  # ============================================================================
  # At this point results_list is a flat list of data frames — one per
  # combination — regardless of whether parallel or serial was used.
  # parLapply guarantees this structure; the serial loop also guarantees it.
  # Remove any NULLs (should not occur given error handling above, but guard
  # defensively), then bind into a single data frame.
  
  if (verbose) cat("\nCombining results...\n")
  
  results_list <- results_list[!sapply(results_list, is.null)]
  
  if (length(results_list) == 0) {
    stop("No valid results were produced. Check that the solver functions are ",
         "correctly sourced and that the parameter grid is valid.")
  }
  
  # dplyr::bind_rows handles any residual column mismatches by filling with NA
  combined_results <- dplyr::bind_rows(results_list)
  rownames(combined_results) <- NULL
  
  # Sort for consistent output: K ascending, r ascending, scenario alphabetical
  combined_results <- combined_results[
    order(combined_results$K, combined_results$r, combined_results$scenario), ]
  rownames(combined_results) <- NULL
  
  # ============================================================================
  # Summary statistics
  # ============================================================================
  
  n_total          <- nrow(combined_results)
  n_feasible       <- sum(combined_results$feasible,  na.rm = TRUE)
  n_converged      <- sum(combined_results$converged, na.rm = TRUE)
  n_failed         <- n_total - n_converged
  expected_total   <- n_combinations * n_scenarios
  
  # Warn if the row count is lower than expected
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
    overall_feasibility_rate = if (n_total > 0) n_feasible / n_total else NA_real_
  )
  
  # ============================================================================
  # Report summary
  # ============================================================================
  
  if (verbose) {
    cat("\n=== CDR SCALE SENSITIVITY ANALYSIS COMPLETE ===\n")
    cat("Total runtime:        ",
        sprintf("%.1f", as.numeric(total_runtime)), "minutes\n")
    cat("Expected rows:        ", expected_total, "\n")
    cat("Total rows:           ", n_total, "\n")
    cat("Converged:            ", n_converged, "/", n_total, "\n")
    cat("Feasible (<=1.5°C):   ", n_feasible, "/", n_total, "\n")
    cat("Overall feasibility:  ",
        sprintf("%.1f%%", 100 * summary_stats$overall_feasibility_rate), "\n")
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
    scenarios                     = scenarios,
    mitigation_delay_years        = mitigation_delay_years,
    cdr_delay_years               = cdr_delay_years,
    use_mitigation_capacity_limit = use_mitigation_capacity_limit,
    use_parallel                  = use_parallel,
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
    
    # Record file paths in run_info before saving the RDS
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
  
} # Close run_cdr_scale_sensitivity


# ==============================================================================
# Section 4: Convenience Wrapper
# ==============================================================================

#' @title Run CDR Scale Sensitivity from User Parameters
#' @description
#' Convenience wrapper that builds the K × r parameter grid and runs the full
#' sensitivity sweep in a single call. Accepts the same user-facing parameters
#' as the workflow chunk and returns the same output as run_cdr_scale_sensitivity().
#'
#' @param K_min Minimum CDR carrying capacity (GtCO2/year)
#' @param K_max Maximum CDR carrying capacity (GtCO2/year)
#' @param n_K Number of K values to test
#' @param r_min Minimum CDR growth rate
#' @param r_max Maximum CDR growth rate
#' @param n_r Number of r values to test
#' @param r_log_scale Logical: use log-spaced r values (default FALSE)
#' @param g_initial Starting CDR deployment level (GtCO2/year)
#' @param t_start Year CDR deployment begins
#' @param scenarios Character vector of SSP scenario names to compare
#' @param parameter_df Single-row data frame of model parameters
#' @param emissions_df Emissions data frame from interpolate_ssp_emissions()
#' @param economic_df Economic data frame from interpolate_ssp_economic()
#' @param use_mitigation_capacity_limit Logical: activate mitigation capacity
#'   constraint (default TRUE). Pair with make_zero_capacity() to eliminate
#'   mitigation so CDR is the sole control.
#' @param mitigation_capacity_function Capacity function for mitigation with
#'   signature function(year). Required when use_mitigation_capacity_limit = TRUE.
#' @param use_parallel Logical: enable parallel processing (default TRUE)
#' @param save_results Logical: save results to output/ directory (default TRUE)
#' @param verbose Logical: print progress and summary (default TRUE)
#'
#' @return List from run_cdr_scale_sensitivity(): combined_results, summary_stats,
#'   run_info. See run_cdr_scale_sensitivity() for full column descriptions.
#'
#' @examples
#' \dontrun{
#' results <- run_cdr_scale_sensitivity_from_params(
#'   K_min        = 25, K_max = 200, n_K = 51,
#'   r_min        = 0.02, r_max = 0.20, n_r = 51,
#'   g_initial    = 2, t_start = 2025,
#'   scenarios    = c("SSP1-Baseline", "SSP5-Baseline"),
#'   parameter_df = parameter_df[1, ],
#'   emissions_df = emissions_df,
#'   economic_df  = economic_df
#' )
#' }
run_cdr_scale_sensitivity_from_params <- function(K_min,
                                                   K_max,
                                                   n_K,
                                                   r_min,
                                                   r_max,
                                                   n_r,
                                                   r_log_scale                   = FALSE,
                                                   g_initial                     = 2,
                                                   t_start                       = 2025,
                                                   scenarios,
                                                   parameter_df,
                                                   emissions_df,
                                                   economic_df,
                                                   use_mitigation_capacity_limit = TRUE,
                                                   mitigation_capacity_function  = make_zero_capacity(),
                                                   use_parallel                  = TRUE,
                                                   save_results                  = TRUE,
                                                   verbose                       = TRUE) {

  cdr_grid <- build_cdr_scale_grid(
    K_min       = K_min,
    K_max       = K_max,
    n_K         = n_K,
    r_min       = r_min,
    r_max       = r_max,
    n_r         = n_r,
    r_log_scale = r_log_scale
  )

  run_cdr_scale_sensitivity(
    parameter_df                  = parameter_df,
    emissions_df                  = emissions_df,
    economic_df                   = economic_df,
    scenarios                     = scenarios,
    cdr_grid                      = cdr_grid,
    g_initial                     = g_initial,
    t_start                       = t_start,
    mitigation_delay_years        = 0,
    cdr_delay_years               = 0,
    use_mitigation_capacity_limit = use_mitigation_capacity_limit,
    mitigation_capacity_function  = mitigation_capacity_function,
    use_parallel                  = use_parallel,
    save_results                  = save_results,
    verbose                       = verbose
  )
}
