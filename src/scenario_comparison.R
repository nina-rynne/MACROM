# ==============================================================================
# Scenario Comparison Analysis Functions
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
# Last updated: January 2026
# ==============================================================================

#' @title Run Scenario Comparison Analysis
#' @description
#' Runs the optimal control algorithm across multiple SSP scenarios to compare
#' the effects of different baseline emission pathways. This provides the
#' multi-scenario comparison functionality needed for publication figures.
#' Note: Currently only supports single parameter sets (single-row parameter_df).
#'
#' @param parameter_df Single-row data frame containing model parameters
#' @param emissions_df Data frame with emissions data for multiple scenarios
#' @param economic_df Data frame with economic data for multiple scenarios  
#' @param scenarios Vector of scenario names to compare (e.g., c("SSP2-Baseline", "SSP3-Baseline"))
#' @param mitigation_delay_years Number of years to delay mitigation start (default: 0)
#' @param cdr_delay_years Number of years to delay carbon dioxide removal start (default: 0)
#' @param use_parallel Whether to use parallel processing across scenarios (default: TRUE)
#' @param verbose Print progress information (default: TRUE)
#' @param save_results Whether to automatically save results (default: TRUE)
#' @param output_dir Directory to save results (default: "output")
#' @param output_prefix Prefix for output filenames (default: "scenario_comparison")
#'
#' @return List containing:
#'   - scenario_results: Named list of results for each scenario
#'   - comparison_summary: Data frame comparing key metrics across scenarios
#'   - year_first_1p5C: Named vector of years when temperature first reaches 1.5°C for each scenario
#'   - year_peak_temp: Named vector of years when peak temperature occurs for each scenario
#'   - year_mitig_capped: Named vector of years when mitigation is first capped by emissions limit
#'   - failed_scenarios: Named list of error messages for scenarios that failed to run
#'   - run_info: List of metadata about the comparison analysis
#'
#' @examples
#' # Compare across 5 SSP scenarios with single parameter set
#' comparison <- run_scenario_comparison(
#'   parameter_df = my_params,
#'   emissions_df = all_emissions_data,
#'   economic_df = all_economic_data,
#'   scenarios = c("SSP1-Baseline", "SSP2-Baseline", "SSP3-Baseline", 
#'                 "SSP4-Baseline", "SSP5-Baseline"),
#'   mitigation_delay_years = 5,
#'   use_parallel = TRUE
#' )

# Main function definition
run_scenario_comparison <- function(parameter_df,
                                    emissions_df,
                                    economic_df,
                                    scenarios,
                                    mitigation_delay_years = 0,  # Years to delay mitigation
                                    cdr_delay_years = 0,         # Years to delay CDR
                                    use_parallel = TRUE,         # Enable parallel processing
                                    verbose = TRUE,              # display progress logging
                                    save_results = TRUE,         # save results as RDS and csv
                                    output_dir = "output",
                                    output_prefix = "scenario_comparison") {
  
  # ============================================================================
  # Input validation
  # ============================================================================
  
  # Validate parameter_df: must be non-empty data frame with exactly one row
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
  
  # Validate scenarios parameter
  if (length(scenarios) == 0) {
    stop("At least one scenario must be specified")
  }
  
  # Verify requested scenarios exist in emissions data
  available_scenarios <- unique(emissions_df$Scenario)
  missing_scenarios <- setdiff(scenarios, available_scenarios)
  
  if (length(missing_scenarios) > 0) {
    if (verbose) {
      cat("⚠ WARNING: The following scenarios are not available in emissions data:\n")
      cat(paste("  -", missing_scenarios, collapse = "\n"), "\n")
    }
    
    # Filter to only available scenarios
    scenarios <- intersect(scenarios, available_scenarios)
    
    if (length(scenarios) == 0) {
      stop("No valid scenarios available for comparison. Available scenarios: ", 
           paste(available_scenarios, collapse = ", "))
    }
    
    if (verbose) {
      cat("Proceeding with available scenarios:", paste(scenarios, collapse = ", "), "\n\n")
    }
  }
  
  # ============================================================================
  # Setup and initialisation
  # ============================================================================
  
  # Print analysis configuration if verbose mode enabled
  if (verbose) {
    cat("=== SCENARIO COMPARISON ANALYSIS ===\n")
    cat("Number of scenarios:", length(scenarios), "\n")
    cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
    cat("Parameter sets:", nrow(parameter_df), "\n")
    cat("Mitigation delay:", mitigation_delay_years, "years\n")
    cat("CDR delay:", cdr_delay_years, "years\n")
    cat("Using parallel processing:", use_parallel, "\n\n")
  }
  
  # Record start time for performance tracking
  start_time <- Sys.time()
  
  # ============================================================================
  # Helper function definition
  # ============================================================================
  
  # Define helper function to run a single scenario with error handling
  # This allows individual scenario failures without stopping the entire comparison
  run_single_scenario <- function(scenario_name) {
    tryCatch({
      # Run optimal control algorithm for this scenario
      result <- optimal_control_shooting(
        parameter_df = parameter_df,
        emissions_df = emissions_df,
        economic_df = economic_df,
        scenario = scenario_name,
        mitigation_delay_years = mitigation_delay_years,
        cdr_delay_years = cdr_delay_years,
        verbose = FALSE  # Suppress individual scenario output
      )
      
      # Return success with result
      return(list(
        success = TRUE,
        scenario = scenario_name,
        result = result
      ))
      
    }, error = function(e) {
      # Return failure with error message
      return(list(
        success = FALSE,
        scenario = scenario_name,
        error = e$message
      ))
    })
  }
  
  # ============================================================================
  # Scenario execution
  # ============================================================================
  
  # Execute scenario comparisons with parallel or serial processing
  if (use_parallel && length(scenarios) > 1) {
    # Parallel execution across scenarios
    if (verbose) cat("Running", length(scenarios), "scenarios in parallel...\n")
    
    parallel_success <- tryCatch({
      # Create cluster using all cores except one (leave one for system)
      n_cores <- max(1, parallel::detectCores() - 1)
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      
      # Load required packages on each worker node
      # Note: Workers start with fresh R sessions and need explicit package loading
      parallel::clusterEvalQ(cl, {
        library(dplyr)  # Required by optimal_control_shooting and optimal_control_solve
      })
      
      # Export functions and data to cluster workers
      parallel::clusterExport(cl, 
                              c("optimal_control_shooting", 
                                "optimal_control_solve",
                                "run_single_scenario",
                                "parameter_df", 
                                "emissions_df", 
                                "economic_df",
                                "mitigation_delay_years",
                                "cdr_delay_years"),
                              envir = environment())
      
      # Run scenarios in parallel using foreach
      scenario_results_list <- foreach::foreach(
        scenario = scenarios, 
        .packages = "dplyr"
      ) %dopar% {
        run_single_scenario(scenario)
      }
      
      # Clean up: stop cluster to free resources
      parallel::stopCluster(cl)
      
      # Return TRUE to indicate parallel execution succeeded
      TRUE
      
    }, error = function(e) {
      # If parallel execution fails, fall back to serial processing
      if (verbose) {
        cat("Parallel execution failed:", e$message, "\n")
        cat("Falling back to serial processing...\n")
      }
      # Return FALSE to trigger serial execution below
      FALSE
    })
    
    # If parallel execution failed, fall back to serial
    if (!parallel_success) {
      use_parallel <- FALSE
    }
  }
  
  # Serial execution (used when use_parallel=FALSE, only 1 scenario, or parallel failed)
  if (!use_parallel || length(scenarios) == 1) {
    if (verbose) cat("Running", length(scenarios), "scenarios in serial...\n")
    
    scenario_results_list <- list()
    for (i in seq_along(scenarios)) {
      scenario <- scenarios[i]
      
      # Print progress with elapsed time
      if (verbose) {
        elapsed <- difftime(Sys.time(), start_time, units = "mins")
        cat("Scenario", i, "of", length(scenarios), 
            "(", scenario, ") - Elapsed:", round(elapsed, 1), "minutes\n")
      }
      
      # Run scenario and store result
      scenario_results_list[[scenario]] <- run_single_scenario(scenario)
    }
  }
  
  # ============================================================================
  # Results processing
  # ============================================================================
  
  # Separate successful results from failures
  # Note: success/failure is determined by run_single_scenario based on 
  # whether optimal_control_shooting throws an error
  scenario_results <- list()
  failed_scenarios <- list()
  
  for (i in seq_along(scenario_results_list)) {
    result <- scenario_results_list[[i]]
    
    if (result$success) {
      # Store successful result under scenario name
      scenario_results[[result$scenario]] <- result$result
    } else {
      # Store error message for failed scenario
      failed_scenarios[[result$scenario]] <- result$error
    }
  }
  
  # ============================================================================
  # Metrics extraction and summary
  # ============================================================================
  
  # Extract metrics and create comparison summary (only if scenarios succeeded)
  if (length(scenario_results) > 0) {
    n_scenarios <- length(scenario_results)
    summary_data <- vector("list", n_scenarios)  # Pre-allocate list
    
    # Initialise metric vectors with scenario names to be added later
    year_first_1p5C <- numeric(n_scenarios)      # Year each scenario first reaches 1.5°C threshold
    year_peak_temp <- numeric(n_scenarios)       # Year of peak temperature for each scenario
    year_mitig_capped <- numeric(n_scenarios)    # Year mitigation hits emissions constraint for each scenario
    
    for (i in seq_along(scenario_results)) {
      scenario_name <- names(scenario_results)[i]
      run <- scenario_results[[i]]
      
      # Extract temporal metrics
      # 1. Year when temperature first reaches 1.5°C threshold
      temp_1p5_indices <- which(run$temperature_anomaly >= 1.5)
      year_first_1p5 <- if (length(temp_1p5_indices) > 0) {
        run$years[temp_1p5_indices[1]]
      } else {
        NA
      }
      
      # 2. Year of peak temperature
      peak_temp_index <- which.max(run$temperature_anomaly)
      year_peak <- run$years[peak_temp_index]
      
      # 3. Year when mitigation is capped by emissions constraint
      # Use epsilon from solution or default to 0.01
      epsilon <- if ("epsilon_used" %in% names(run)) run$epsilon_used else 0.01
      
      # Detect when mitigation hits constraint: qty_mitig ≈ baseline_emissions - epsilon
      # Tolerance of 0.001 accounts for numerical precision
      tolerance <- 0.001
      max_allowed_mitig <- run$baseline_annual_emissions - epsilon
      capped_indices <- which(abs(run$qty_mitig - max_allowed_mitig) < tolerance)
      year_mitig_cap <- if (length(capped_indices) > 0) {
        run$years[capped_indices[1]]
      } else {
        NA
      }
      
      # Store temporal metrics in pre-allocated vectors
      year_first_1p5C[i] <- year_first_1p5
      year_peak_temp[i] <- year_peak
      year_mitig_capped[i] <- year_mitig_cap
      
      # Create summary row for this scenario
      summary_data[[i]] <- data.frame(
        scenario = scenario_name,
        final_emissions = run$final_emissions,
        final_temperature = run$final_temperature,
        peak_temperature = max(run$temperature_anomaly),
        total_cost = run$total_cost,
        mitig_cost = run$mitig_cost,
        remov_cost = run$remov_cost,
        temp_cost = run$temp_cost,
        total_mitigation_units = sum(run$qty_mitig),
        total_cdr_units = sum(run$qty_remov),
        years_above_1p5 = sum(run$temperature_anomaly > 1.5),
        year_first_1p5C = year_first_1p5,
        year_peak_temp = year_peak,
        year_mitig_capped = year_mitig_cap,
        converged = run$converged,
        iterations = run$iterations,
        stringsAsFactors = FALSE
      )
    }
    
    # Combine all scenario summaries into single data frame
    comparison_summary <- do.call(rbind, summary_data)
    
    # Add scenario names to metric vectors for easy reference
    names(year_first_1p5C) <- names(scenario_results)
    names(year_peak_temp) <- names(scenario_results)
    names(year_mitig_capped) <- names(scenario_results)
  } else {
    # All scenarios failed - return empty/NULL structures
    comparison_summary <- NULL
    year_first_1p5C <- numeric(0)      # Year each scenario first reaches 1.5°C threshold
    year_peak_temp <- numeric(0)       # Year of peak temperature for each scenario
    year_mitig_capped <- numeric(0)    # Year mitigation hits emissions constraint for each scenario
  }
  
  # ============================================================================
  # Final reporting and return
  # ============================================================================
  
  # Calculate total runtime
  total_time <- difftime(Sys.time(), start_time, units = "mins")
  
  # Print summary if verbose mode enabled
  if (verbose) {
    cat("\n=== SCENARIO COMPARISON COMPLETE ===\n")
    cat("Total time:", sprintf("%.1f", total_time), "minutes\n")
    cat("Successful scenarios:", length(scenario_results), "out of", length(scenarios), "\n")
    cat("Failed scenarios:", length(failed_scenarios), "\n")
    
    # List failed scenarios if any
    if (length(failed_scenarios) > 0) {
      cat("Failed scenario names:", paste(names(failed_scenarios), collapse = ", "), "\n")
    }
    
    # Report key metrics if any scenarios succeeded
    if (!is.null(comparison_summary)) {
      cat("\nScenario comparison summary:\n")
      cat("Peak temperature range:", sprintf("%.2f - %.2f°C",
                                             min(comparison_summary$peak_temperature),
                                             max(comparison_summary$peak_temperature)), "\n")
      cat("Total cost range:", sprintf("%.1f - %.1f trillion $",
                                       min(comparison_summary$total_cost),
                                       max(comparison_summary$total_cost)), "\n")
      
      # Report 1.5°C crossing timing
      valid_1p5_years <- year_first_1p5C[!is.na(year_first_1p5C)]
      if (length(valid_1p5_years) > 0) {
        cat("First 1.5°C crossing range:", sprintf("%d - %d",
                                                   min(valid_1p5_years),
                                                   max(valid_1p5_years)), "\n")
      } else {
        cat("No scenarios reach 1.5°C\n")
      }
      
      # Report mitigation capping timing
      valid_cap_years <- year_mitig_capped[!is.na(year_mitig_capped)]
      if (length(valid_cap_years) > 0) {
        cat("Mitigation capping range:", sprintf("%d - %d",
                                                 min(valid_cap_years),
                                                 max(valid_cap_years)), "\n")
      } else {
        cat("No scenarios have mitigation capped by emissions limit\n")
      }
    }
  }
  
  # Save results if requested
  if (save_results) {
    save_scenario_results(
      scenario_results = list(
        scenario_results = scenario_results,
        comparison_summary = comparison_summary,
        year_first_1p5C = year_first_1p5C,
        year_peak_temp = year_peak_temp,
        year_mitig_capped = year_mitig_capped,
        failed_scenarios = failed_scenarios,
        run_info = list(
          scenarios = scenarios,
          n_scenarios = length(scenarios),
          n_successful = length(scenario_results),
          n_failed = length(failed_scenarios),
          parameter_sets = nrow(parameter_df),
          mitigation_delay_years = mitigation_delay_years,
          cdr_delay_years = cdr_delay_years,
          use_parallel = use_parallel,
          start_time = start_time,
          end_time = Sys.time(),
          total_time_minutes = as.numeric(total_time)
        )
      ),
      output_dir = output_dir,
      output_prefix = output_prefix,
      save_csv = TRUE,
      verbose = verbose
    )
  }
  
  # Return comprehensive results list
  return(list(
    scenario_results = scenario_results,
    comparison_summary = comparison_summary,
    year_first_1p5C = year_first_1p5C,
    year_peak_temp = year_peak_temp,
    year_mitig_capped = year_mitig_capped,
    failed_scenarios = failed_scenarios,
    run_info = list(
      scenarios = scenarios,
      n_scenarios = length(scenarios),
      n_successful = length(scenario_results),
      n_failed = length(failed_scenarios),
      parameter_sets = nrow(parameter_df),
      mitigation_delay_years = mitigation_delay_years,
      cdr_delay_years = cdr_delay_years,
      use_parallel = use_parallel,
      start_time = start_time,
      end_time = Sys.time(),
      total_time_minutes = as.numeric(total_time)
    )
  ))
} # Close the run_scenario_comparison function


#' @title Save Scenario Comparison Results
#' @description
#' Saves scenario comparison results to RDS and optionally CSV, with summary reporting
#'
#' @param scenario_results Output from run_scenario_comparison()
#' @param output_dir Directory to save results (default: "output")
#' @param output_prefix Prefix for output filenames (default: "scenario_results")
#' @param save_csv Whether to also save comparison summary as CSV (default: TRUE)
#' @param verbose Print summary information (default: TRUE)
#'
#' @return Invisibly returns a list with filepaths of saved files
save_scenario_results <- function(scenario_results,
                                  output_dir = "output",
                                  output_prefix = "scenario_results",
                                  save_csv = TRUE,
                                  verbose = TRUE) {
  
  # ============================================================================
  # File saving
  # ============================================================================
  
  # Generate timestamp for unique filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Save main results as RDS file
  rds_filename <- paste0(output_prefix, "_", timestamp, ".rds")
  rds_filepath <- here::here(output_dir, rds_filename)
  saveRDS(scenario_results, rds_filepath)
  
  if (verbose) {
    cat("\n=== RESULTS SAVED ===\n")
    cat("RDS file:", rds_filepath, "\n")
  }
  
  # Save comparison summary as CSV if requested and available
  csv_filepath <- NULL
  if (save_csv && !is.null(scenario_results$comparison_summary)) {
    csv_filename <- paste0(output_prefix, "_summary_", timestamp, ".csv")
    csv_filepath <- here::here(output_dir, csv_filename)
    write.csv(scenario_results$comparison_summary, csv_filepath, row.names = FALSE)
    
    if (verbose) {
      cat("CSV summary:", csv_filepath, "\n")
    }
  }
  
  # ============================================================================
  # Extended summary reporting
  # ============================================================================
  
  # Print extended summary statistics if verbose mode enabled
  if (verbose && !is.null(scenario_results$comparison_summary)) {
    cat("\n=== EXTENDED SUMMARY ===\n")
    
    summary_df <- scenario_results$comparison_summary
    
    # Report temperature range across scenarios
    cat("Peak temperature range:", 
        sprintf("%.2f - %.2f°C", 
                min(summary_df$peak_temperature, na.rm = TRUE),
                max(summary_df$peak_temperature, na.rm = TRUE)), "\n")
    
    # Report cost range across scenarios
    cat("Total cost range:", 
        sprintf("%.1f - %.1f trillion $", 
                min(summary_df$total_cost, na.rm = TRUE),
                max(summary_df$total_cost, na.rm = TRUE)), "\n")
    
    # Identify scenarios with extreme values
    min_cost_scenario <- summary_df$scenario[which.min(summary_df$total_cost)]
    max_cost_scenario <- summary_df$scenario[which.max(summary_df$total_cost)]
    min_temp_scenario <- summary_df$scenario[which.min(summary_df$peak_temperature)]
    max_temp_scenario <- summary_df$scenario[which.max(summary_df$peak_temperature)]
    
    cat("Lowest cost scenario:", min_cost_scenario, "\n")
    cat("Highest cost scenario:", max_cost_scenario, "\n")
    cat("Lowest peak temperature scenario:", min_temp_scenario, "\n")
    cat("Highest peak temperature scenario:", max_temp_scenario, "\n")
    
    # Report mitigation and CDR ranges
    cat("Total mitigation range:", 
        sprintf("%.1f - %.1f GtCO2", 
                min(summary_df$total_mitigation_units, na.rm = TRUE),
                max(summary_df$total_mitigation_units, na.rm = TRUE)), "\n")
    
    cat("Total CDR range:", 
        sprintf("%.1f - %.1f GtCO2", 
                min(summary_df$total_cdr_units, na.rm = TRUE),
                max(summary_df$total_cdr_units, na.rm = TRUE)), "\n")
  }
  
  # Return filepaths invisibly for potential further use
  invisible(list(
    rds_filepath = rds_filepath,
    csv_filepath = csv_filepath
  ))
} # Close the save_scenario_results function