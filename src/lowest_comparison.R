# ==============================================================================
# Lowest Comparison Analysis Functions
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

#' @title Run Lowest Comparison Analysis
#' @description
#' Runs the optimal control algorithm for SSP-specific (r, K) parameter pairs.
#' Unlike run_capacity_growth_comparison(), which applies every r value uniformly
#' across all SSPs, this function allows each SSP to have its own set of CDR
#' logistic parameter combinations. Each combination is a specific (r, K) pair —
#' not all combinations of the supplied values.
#'
#' Pairs are auto-labelled pair_1, pair_2, ... in the order provided.
#' Parallelism is applied across all (SSP, pair) combinations simultaneously.
#'
#' @param parameter_df Single-row data frame containing model parameters.
#' @param emissions_df Data frame with emissions data for multiple scenarios.
#' @param economic_df Data frame with economic data for multiple scenarios.
#' @param ssp_params Named list. Each element name is an SSP scenario string
#'   (e.g., "SSP3-Baseline") and each element is a list of numeric vectors of
#'   the form c(r = <value>, K = <value>), one per parameter combination.
#'   Pairs are taken as-is — no cross-product is computed.
#'   Example:
#'     list(
#'       "SSP3-Baseline" = list(c(r = 0.045, K = 95), c(r = 0.07, K = 25)),
#'       "SSP5-Baseline" = list(c(r = 0.05, K = 100))
#'     )
#' @param g_initial Starting CDR deployment level in GtCO2/year, shared across
#'   all SSPs and pairs. Passed to make_logistic_from_zero() (default: 2).
#' @param t_start Year CDR deployment begins, shared across all SSPs and pairs.
#'   Passed to make_logistic_from_zero() (default: 2025).
#' @param mitigation_delay_years Years to delay mitigation start (default: 0).
#' @param cdr_delay_years Years to delay CDR deployment start (default: 0).
#' @param use_mitigation_capacity_limit Logical; activate mitigation capacity
#'   constraint (default: TRUE).
#' @param mitigation_capacity_function Capacity function for mitigation
#'   (default: make_zero_capacity()).
#' @param use_parallel Logical; parallelise across all (SSP, pair) combinations
#'   (default: TRUE).
#' @param save_results Logical; save combined RDS on completion (default: TRUE).
#' @param output_dir Directory for saved files (default: "output").
#' @param output_prefix Prefix for saved filename (default: "lowest_comparison").
#' @param verbose Logical; print progress messages (default: TRUE).
#'
#' @return Nested list with structure results[[ssp_name]][[pair_label]], where
#'   each leaf is the full output of run_scenario_comparison() for that single
#'   SSP and parameter pair. pair_label is "pair_1", "pair_2", etc.
#'
#' @examples
#' lowest_results <- run_lowest_comparison(
#'   parameter_df = parameter_df[1, ],
#'   emissions_df = emissions_df,
#'   economic_df  = economic_df,
#'   ssp_params   = list(
#'     "SSP3-Baseline" = list(c(r = 0.045, K = 95), c(r = 0.07, K = 25)),
#'     "SSP5-Baseline" = list(c(r = 0.05,  K = 100), c(r = 0.09, K = 60))
#'   ),
#'   g_initial = 2,
#'   t_start   = 2025
#' )
#'
#' # Access a specific SSP/pair result
#' lowest_results[["SSP3-Baseline"]][["pair_1"]]$comparison_summary
#' lowest_results[["SSP3-Baseline"]][["pair_2"]]$scenario_results[["SSP3-Baseline"]]

run_lowest_comparison <- function(parameter_df,
                                  emissions_df,
                                  economic_df,
                                  ssp_params,
                                  g_initial                     = 2,
                                  t_start                       = 2025,
                                  mitigation_delay_years        = 0,
                                  cdr_delay_years               = 0,
                                  use_mitigation_capacity_limit = TRUE,
                                  mitigation_capacity_function  = make_zero_capacity(),
                                  use_parallel                  = TRUE,
                                  save_results                  = TRUE,
                                  output_dir                    = "output",
                                  output_prefix                 = "lowest_comparison",
                                  verbose                       = TRUE) {

  # ============================================================================
  # Input validation
  # ============================================================================

  if (!is.data.frame(parameter_df)) {
    stop("parameter_df must be a data frame (received ", class(parameter_df)[1], ")")
  }
  if (nrow(parameter_df) == 0) {
    stop("parameter_df must contain at least one row (received empty data frame)")
  }
  if (nrow(parameter_df) > 1) {
    stop("parameter_df must contain exactly one row (received ", nrow(parameter_df),
         " rows). Multi-parameter sensitivity analysis is not yet implemented.")
  }

  if (!is.list(ssp_params) || is.null(names(ssp_params)) || any(names(ssp_params) == "")) {
    stop("ssp_params must be a fully named list with one element per SSP scenario")
  }

  # Verify requested SSPs exist in emissions data
  available_scenarios <- unique(emissions_df$Scenario)
  missing_ssps        <- setdiff(names(ssp_params), available_scenarios)
  if (length(missing_ssps) > 0) {
    stop("The following SSPs in ssp_params are not in emissions_df: ",
         paste(missing_ssps, collapse = ", "),
         "\nAvailable: ", paste(available_scenarios, collapse = ", "))
  }

  # Validate each SSP's pair list and auto-label pairs
  labelled_params <- vector("list", length(ssp_params))
  names(labelled_params) <- names(ssp_params)

  for (ssp in names(ssp_params)) {
    pairs <- ssp_params[[ssp]]

    if (!is.list(pairs) || length(pairs) == 0) {
      stop("ssp_params[['", ssp, "']] must be a non-empty list of c(r = ..., K = ...) vectors")
    }

    pair_labels <- paste0("pair_", seq_along(pairs))
    named_pairs <- setNames(pairs, pair_labels)

    for (lbl in pair_labels) {
      pair <- named_pairs[[lbl]]
      if (!all(c("r", "K") %in% names(pair))) {
        stop("Each pair must be a named vector with 'r' and 'K'. ",
             "Problem in ssp_params[['", ssp, "']][['", lbl, "']]")
      }
      if (pair[["r"]] <= 0) {
        stop("r must be positive. Problem in ssp_params[['", ssp, "']][['", lbl, "']]")
      }
      if (pair[["K"]] <= g_initial) {
        stop("K must be greater than g_initial (", g_initial, "). ",
             "Problem in ssp_params[['", ssp, "']][['", lbl, "']]")
      }
    }

    labelled_params[[ssp]] <- named_pairs
  }

  if (use_mitigation_capacity_limit && is.null(mitigation_capacity_function)) {
    stop("mitigation_capacity_function must be supplied when use_mitigation_capacity_limit = TRUE")
  }
  if (use_mitigation_capacity_limit && !is.function(mitigation_capacity_function)) {
    stop("mitigation_capacity_function must be a function")
  }

  # ============================================================================
  # Build flat task list for parallel dispatch
  # ============================================================================
  # Each task is one (SSP, pair_label) combination.

  tasks <- list()
  for (ssp in names(labelled_params)) {
    for (lbl in names(labelled_params[[ssp]])) {
      tasks[[length(tasks) + 1]] <- list(
        ssp   = ssp,
        label = lbl,
        r     = labelled_params[[ssp]][[lbl]][["r"]],
        K     = labelled_params[[ssp]][[lbl]][["K"]]
      )
    }
  }

  n_tasks <- length(tasks)

  if (verbose) {
    cat("=== LOWEST COMPARISON ANALYSIS ===\n")
    cat("SSPs:         ", paste(names(labelled_params), collapse = ", "), "\n")
    cat("Total runs:   ", n_tasks, "\n")
    cat("g_initial:    ", g_initial, "\n")
    cat("t_start:      ", t_start, "\n")
    cat("Parallel:     ", use_parallel, "\n\n")
  }

  start_time <- Sys.time()

  # ============================================================================
  # Helper: run one (SSP, pair) task
  # ============================================================================

  run_one_task <- function(task) {
    tryCatch({
      cdr_cap_fn <- make_logistic_from_zero(
        g_initial = g_initial,
        K         = task$K,
        r         = task$r,
        t_start   = t_start
      )

      result <- run_scenario_comparison(
        parameter_df                  = parameter_df,
        emissions_df                  = emissions_df,
        economic_df                   = economic_df,
        scenarios                     = task$ssp,
        mitigation_delay_years        = mitigation_delay_years,
        cdr_delay_years               = cdr_delay_years,
        use_mitigation_capacity_limit = use_mitigation_capacity_limit,
        mitigation_capacity_function  = mitigation_capacity_function,
        use_cdr_capacity_limit        = TRUE,
        cdr_capacity_function         = cdr_cap_fn,
        use_parallel                  = FALSE,   # single SSP — nothing to parallelise within
        save_results                  = FALSE,   # combined save handled below
        verbose                       = FALSE
      )

      # Return r and K alongside the result so the main process can attach
      # run_info and the capacity curve after assembly (mirrors the approach
      # used in run_capacity_growth_comparison, where all post-processing
      # happens in the main process loop, not inside parallel workers)
      list(success = TRUE, ssp = task$ssp, label = task$label,
           r = task$r, K = task$K, result = result)

    }, error = function(e) {
      list(success = FALSE, ssp = task$ssp, label = task$label, error = e$message)
    })
  }

  # ============================================================================
  # Execution: parallel or serial
  # ============================================================================

  if (use_parallel && n_tasks > 1) {

    parallel_ok <- tryCatch({
      n_cores <- max(1, parallel::detectCores() - 1)
      cl      <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)

      parallel::clusterEvalQ(cl, { library(dplyr) })

      parallel::clusterExport(
        cl,
        c("run_one_task",
          "run_scenario_comparison",
          "optimal_control_shooting",
          "optimal_control_solve",
          "make_logistic_from_zero",
          "make_zero_capacity",
          "parameter_df",
          "emissions_df",
          "economic_df",
          "g_initial",
          "t_start",
          "mitigation_delay_years",
          "cdr_delay_years",
          "use_mitigation_capacity_limit",
          "mitigation_capacity_function"),
        envir = environment()
      )

      task_results_list <- foreach::foreach(
        task = tasks,
        .packages = "dplyr"
      ) %dopar% {
        run_one_task(task)
      }

      parallel::stopCluster(cl)
      TRUE

    }, error = function(e) {
      if (verbose) {
        cat("Parallel execution failed:", e$message, "\n")
        cat("Falling back to serial processing...\n")
      }
      FALSE
    })

    if (!parallel_ok) {
      use_parallel <- FALSE
    }
  }

  if (!use_parallel || n_tasks == 1) {
    if (verbose) cat("Running", n_tasks, "task(s) in serial...\n")

    task_results_list <- vector("list", n_tasks)
    for (i in seq_along(tasks)) {
      task <- tasks[[i]]
      if (verbose) {
        elapsed <- difftime(Sys.time(), start_time, units = "mins")
        cat("Task", i, "of", n_tasks,
            "| SSP:", task$ssp, "| pair:", task$label,
            "| r =", task$r, "| K =", task$K,
            "| Elapsed:", round(elapsed, 1), "min\n")
      }
      task_results_list[[i]] <- run_one_task(task)
    }
  }

  # ============================================================================
  # Assemble nested results list
  # ============================================================================

  all_results    <- vector("list", length(labelled_params))
  names(all_results) <- names(labelled_params)

  for (ssp in names(all_results)) {
    all_results[[ssp]] <- vector("list", length(labelled_params[[ssp]]))
    names(all_results[[ssp]]) <- names(labelled_params[[ssp]])
  }

  failed_tasks <- list()

  for (tr in task_results_list) {
    if (tr$success) {
      result <- tr$result

      # Attach capacity metadata and pre-evaluated capacity curve here, in the
      # main process — mirrors run_capacity_growth_comparison, which does the
      # same thing after run_scenario_comparison() returns in its main loop.
      result$run_info$r_value   <- tr$r
      result$run_info$K         <- tr$K
      result$run_info$g_initial <- g_initial
      result$run_info$t_start   <- t_start

      if (length(result$scenario_results) > 0) {
        years_vec <- result$scenario_results[[1]]$years
        cdr_cap_fn <- make_logistic_from_zero(
          g_initial = g_initial,
          K         = tr$K,
          r         = tr$r,
          t_start   = t_start
        )
        result$cdr_capacity_curve <- vapply(years_vec, cdr_cap_fn, numeric(1))
      }

      all_results[[tr$ssp]][[tr$label]] <- result
    } else {
      key <- paste0(tr$ssp, " / ", tr$label)
      failed_tasks[[key]] <- tr$error
    }
  }

  # ============================================================================
  # Final reporting
  # ============================================================================

  total_time <- difftime(Sys.time(), start_time, units = "mins")

  if (verbose) {
    cat("\n=== LOWEST COMPARISON COMPLETE ===\n")
    cat("Total time:    ", sprintf("%.1f", total_time), "minutes\n")
    cat("Successful:    ", n_tasks - length(failed_tasks), "of", n_tasks, "\n")
    cat("Failed:        ", length(failed_tasks), "\n")

    if (length(failed_tasks) > 0) {
      cat("Failed tasks:\n")
      for (k in names(failed_tasks)) cat("  -", k, ":", failed_tasks[[k]], "\n")
    }
  }

  # ============================================================================
  # Save combined RDS
  # ============================================================================

  if (save_results) {
    timestamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
    rds_filename <- paste0(output_prefix, "_", timestamp, ".rds")
    rds_path     <- here::here(output_dir, rds_filename)
    saveRDS(all_results, rds_path)

    if (verbose) cat("Combined results saved to:", rds_path, "\n")
  }

  return(all_results)
}
