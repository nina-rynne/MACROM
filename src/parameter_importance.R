# ==============================================================================
# Parameter Importance Analysis Functions
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

#' @title Run Parameter Importance Analysis with Multiple Parameter Sets
#' @description
#' Runs the optimal control algorithm for multiple parameter sets (typically from
#' Latin Hypercube Sampling) to explore the solution space and assess parameter
#' importance. This provides the parameter sensitivity functionality needed for
#' understanding which model parameters most strongly influence outcomes.
#'
#' @param parameter_df Data frame with multiple parameter sets (one per row)
#' @param emissions_df Data frame with emissions scenario data
#' @param economic_df Data frame with economic scenario data  
#' @param scenario Scenario name to analyse (e.g., "SSP3-Baseline")
#' @param mitigation_delay_years Number of years to delay mitigation start (default: 0)
#' @param cdr_delay_years Number of years to delay carbon dioxide removal start (default: 0)
#' @param use_parallel Whether to use parallel processing (default: TRUE)
#' @param n_cores Number of cores for parallel processing (default: auto-detect)
#' @param save_intermediate Save results periodically during run (default: TRUE)
#' @param save_results Whether to automatically save final results (default: TRUE)
#' @param output_dir Directory to save results (default: "output")
#' @param output_prefix Prefix for output filenames (default: "parameter_importance")
#' @param verbose Print progress information (default: TRUE)
#'
#' @return List containing:
#'   - successful_runs: List of successful optimisation results
#'   - failed_runs: List of failed runs with error information
#'   - summary: Data frame with key metrics from successful runs
#'   - run_info: Metadata about the parameter importance analysis
#'
#' @examples
#' # Standard parameter importance analysis with LHS parameter sets
#' results <- run_parameter_importance(
#'   parameter_df = lhs_parameters,
#'   emissions_df = emissions_data,
#'   economic_df = economic_data,
#'   scenario = "SSP3-Baseline",
#'   mitigation_delay_years = 0,
#'   cdr_delay_years = 0
#' )
#'
#' # Serial processing for debugging without auto-save
#' results <- run_parameter_importance(
#'   parameter_df = lhs_parameters[1:10, ],
#'   emissions_df = emissions_data,
#'   economic_df = economic_data,
#'   scenario = "SSP3-Baseline",
#'   use_parallel = FALSE,
#'   save_results = FALSE
#' )

# Main function definition
run_parameter_importance <- function(parameter_df,
                                     emissions_df,
                                     economic_df,
                                     scenario,
                                     mitigation_delay_years = 0,  # Years to delay mitigation
                                     cdr_delay_years = 0,         # Years to delay CDR
                                     use_parallel = TRUE,         # Enable parallel processing
                                     n_cores = NULL,              # Number of cores (auto-detect if NULL)
                                     save_intermediate = TRUE,    # Save results periodically during run
                                     save_results = TRUE,         # Save final results as RDS and CSV
                                     output_dir = "output",       # Directory to save results
                                     output_prefix = "parameter_importance",  # Prefix for saved files
                                     verbose = TRUE) {            # Display progress logging
  
  # ============================================================================
  # Input validation
  # ============================================================================
  
  # Validate parameter_df: must be non-empty data frame with multiple rows expected
  if (!is.data.frame(parameter_df)) {
    stop("parameter_df must be a data frame (received ", class(parameter_df)[1], ")")
  }
  if (nrow(parameter_df) == 0) {
    stop("parameter_df must contain at least one row (received empty data frame)")
  }
  
  # Validate emissions_df and economic_df
  if (!is.data.frame(emissions_df)) {
    stop("emissions_df must be a data frame (received ", class(emissions_df)[1], ")")
  }
  if (!is.data.frame(economic_df)) {
    stop("economic_df must be a data frame (received ", class(economic_df)[1], ")")
  }
  
  # Count number of parameter sets to run
  n_samples <- nrow(parameter_df)
  
  # Validate scenario exists in emissions data
  available_scenarios <- unique(emissions_df$Scenario)
  if (!scenario %in% available_scenarios) {
    stop("Scenario '", scenario, "' not found in emissions data. ",
         "Available scenarios: ", paste(available_scenarios, collapse = ", "))
  }
  
  # Validate n_cores if provided
  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || n_cores < 1 || n_cores != as.integer(n_cores)) {
      stop("n_cores must be a positive integer (received ", n_cores, ")")
    }
  }
  
  # Warn if only one parameter set (parameter importance needs multiple sets)
  if (n_samples == 1) {
    warning("Only 1 parameter set provided. Parameter importance analysis requires multiple parameter sets. ",
            "Consider increasing the number of LHS samples.")
  }
  
  # Validate delay parameters
  if (mitigation_delay_years < 0 || mitigation_delay_years != as.integer(mitigation_delay_years)) {
    stop("mitigation_delay_years must be a non-negative integer (received ", 
         mitigation_delay_years, ")")
  }
  if (cdr_delay_years < 0 || cdr_delay_years != as.integer(cdr_delay_years)) {
    stop("cdr_delay_years must be a non-negative integer (received ", 
         cdr_delay_years, ")")
  }
  
  # ============================================================================
  # Setup and initialisation
  # ============================================================================
  
  # Print analysis configuration if verbose mode enabled
  if (verbose) {
    cat("=== PARAMETER IMPORTANCE ANALYSIS ===\n")
    cat("Number of parameter samples:", n_samples, "\n")
    cat("Scenario:", scenario, "\n")
    cat("Mitigation delay:", mitigation_delay_years, "years\n")
    cat("CDR delay:", cdr_delay_years, "years\n")
    cat("Using parallel processing:", use_parallel, "\n")
    if (use_parallel) {
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
      }
      cat("Number of cores:", n_cores, "\n")
    }
    cat("\n")
  }
  
  # Auto-detect cores if parallel processing enabled but cores not specified
  if (use_parallel && is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Record start time for performance tracking
  start_time <- Sys.time()
  
  # Initialise storage for results
  successful_runs <- list()
  failed_runs <- list()
  
  # ============================================================================
  # Helper function definition
  # ============================================================================
  
  # Define helper function to run a single parameter set with error handling
  # This allows individual parameter set failures without stopping the entire analysis
  run_single_param_set <- function(i) {
    param_row <- parameter_df[i, ]
    
    tryCatch({
      # Run optimal control algorithm for this parameter set
      result <- optimal_control_shooting(
        parameter_df = param_row,
        emissions_df = emissions_df,
        economic_df = economic_df,
        scenario = scenario,
        mitigation_delay_years = mitigation_delay_years,
        cdr_delay_years = cdr_delay_years,
        verbose = FALSE  # Suppress individual parameter set output
      )
      
      # Return success with result
      return(list(
        success = TRUE,
        param_id = i,
        result = result,
        parameters = param_row
      ))
      
    }, error = function(e) {
      # Return failure with error message
      return(list(
        success = FALSE,
        param_id = i,
        error = e$message,
        parameters = param_row
      ))
    })
  }
  
  # ============================================================================
  # Parameter set execution
  # ============================================================================
  
  # Execute parameter sets with parallel or serial processing
  if (use_parallel && n_samples > 1) {
    # Parallel execution across parameter sets
    if (verbose) cat("Running", n_samples, "parameter sets in parallel...\n")
    
    parallel_success <- tryCatch({
      # Create cluster using all cores except one (leave one for system)
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      
      # Load required packages on each worker node
      parallel::clusterEvalQ(cl, {
        library(dplyr)
      })
      
      # Export necessary functions and variables to workers
      parallel::clusterExport(cl, 
                              c("optimal_control_shooting", 
                                "optimal_control_solve",
                                "emissions_df", 
                                "economic_df", 
                                "scenario",
                                "mitigation_delay_years", 
                                "cdr_delay_years"),
                              envir = environment())
      
      # Run parallel optimisation
      all_results <- foreach::foreach(i = 1:n_samples, .packages = "dplyr") %dopar% {
        run_single_param_set(i)
      }
      
      # Clean up parallel processing
      parallel::stopCluster(cl)
      
      TRUE  # Indicate parallel processing succeeded
      
    }, error = function(e) {
      if (verbose) {
        cat("⚠ Parallel processing failed:", e$message, "\n")
        cat("Falling back to serial processing...\n")
      }
      FALSE  # Indicate parallel processing failed
    })
    
    # Fall back to serial processing if parallel failed
    if (!parallel_success) {
      use_parallel <- FALSE  # Update flag for run_info
    }
  }
  
  # Serial execution (either by choice or fallback from parallel failure)
  if (!use_parallel || n_samples == 1) {
    if (verbose) {
      if (n_samples == 1) {
        cat("Running single parameter set (parallel processing not used)...\n")
      } else {
        cat("Running", n_samples, "parameter sets in serial...\n")
      }
    }
    
    all_results <- list()
    for (i in 1:n_samples) {
      # Progress reporting
      if (verbose && (i %% 10 == 0 || i <= 5)) {
        elapsed <- difftime(Sys.time(), start_time, units = "mins")
        cat("Sample", i, "of", n_samples, "- Elapsed:", round(elapsed, 1), "minutes\n")
      }
      
      # Run single parameter set
      all_results[[i]] <- run_single_param_set(i)
      
      # Intermediate save if requested
      if (save_intermediate && i %% 50 == 0) {
        temp_file <- paste0(output_prefix, "_temp_", i, ".rds")
        saveRDS(all_results[1:i], temp_file)
        if (verbose) cat("Intermediate save:", temp_file, "\n")
      }
    }
  }
  
  # ============================================================================
  # Results processing
  # ============================================================================
  
  # Separate successful and failed runs
  if (verbose) cat("Processing results...\n")
  
  for (i in 1:length(all_results)) {
    result <- all_results[[i]]
    
    if (result$success) {
      # Store successful run with descriptive name
      run_name <- paste0("run_", result$param_id)
      successful_runs[[run_name]] <- result$result
    } else {
      # Store failed run with error information
      failed_name <- paste0("failed_", result$param_id)
      failed_runs[[failed_name]] <- result
    }
  }
  
  # ============================================================================
  # Summary statistics extraction
  # ============================================================================
  
  # Create summary data frame from successful runs
  summary_df <- NULL
  if (length(successful_runs) > 0) {
    summary_data <- list()
    
    for (i in 1:length(successful_runs)) {
      run <- successful_runs[[i]]
      
      # Extract key metrics from this run
      summary_data[[i]] <- data.frame(
        run_id = names(successful_runs)[i],
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
        converged = run$converged,
        iterations = run$iterations,
        stringsAsFactors = FALSE
      )
    }
    
    # Combine all run summaries into single data frame
    summary_df <- do.call(rbind, summary_data)
  }
  
  # ============================================================================
  # Final reporting and return
  # ============================================================================
  
  # Calculate total runtime
  total_time <- difftime(Sys.time(), start_time, units = "mins")
  
  # Print summary report if verbose mode enabled
  if (verbose) {
    cat("\n=== PARAMETER IMPORTANCE ANALYSIS COMPLETE ===\n")
    cat("Total time:", round(total_time, 1), "minutes\n")
    cat("Successful runs:", length(successful_runs), "out of", n_samples, 
        "(", round(100 * length(successful_runs) / n_samples, 1), "%)\n")
    cat("Failed runs:", length(failed_runs), "\n")
    
    # Report key metrics if any runs succeeded
    if (!is.null(summary_df)) {
      cat("\nParameter importance summary:\n")
      cat("Peak temperature range:", sprintf("%.2f - %.2f°C", 
                                             min(summary_df$peak_temperature),
                                             max(summary_df$peak_temperature)), "\n")
      cat("Total cost range:", sprintf("%.1f - %.1f trillion $", 
                                       min(summary_df$total_cost),
                                       max(summary_df$total_cost)), "\n")
      
      # Report mitigation and CDR ranges
      cat("Total mitigation range:", sprintf("%.1f - %.1f GtCO2", 
                                             min(summary_df$total_mitigation_units),
                                             max(summary_df$total_mitigation_units)), "\n")
      cat("Total CDR range:", sprintf("%.1f - %.1f GtCO2", 
                                      min(summary_df$total_cdr_units),
                                      max(summary_df$total_cdr_units)), "\n")
      
      # Report convergence statistics
      n_converged <- sum(summary_df$converged)
      cat("Converged runs:", n_converged, "out of", nrow(summary_df),
          "(", round(100 * n_converged / nrow(summary_df), 1), "%)\n")
    }
  }
  
  # Clean up intermediate files if they were created
  if (save_intermediate && length(successful_runs) > 0) {
    temp_pattern <- paste0(output_prefix, "_temp_.*\\.rds$")
    temp_files <- list.files(pattern = temp_pattern, full.names = TRUE)
    if (length(temp_files) > 0) {
      file.remove(temp_files)
      if (verbose) {
        cat("Cleaned up", length(temp_files), "intermediate temporary files\n")
      }
    }
  }
  
  # Save results if requested
  if (save_results) {
    save_parameter_importance_results(
      parameter_results = list(
        successful_runs = successful_runs,
        failed_runs = failed_runs,
        summary = summary_df,
        run_info = list(
          n_samples = n_samples,
          n_successful = length(successful_runs),
          n_failed = length(failed_runs),
          scenario = scenario,
          mitigation_delay_years = mitigation_delay_years,
          cdr_delay_years = cdr_delay_years,
          use_parallel = use_parallel,
          n_cores = if (use_parallel) n_cores else 1,
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
    successful_runs = successful_runs,
    failed_runs = failed_runs,
    summary = summary_df,
    run_info = list(
      n_samples = n_samples,
      n_successful = length(successful_runs),
      n_failed = length(failed_runs),
      scenario = scenario,
      mitigation_delay_years = mitigation_delay_years,
      cdr_delay_years = cdr_delay_years,
      use_parallel = use_parallel,
      n_cores = if (use_parallel) n_cores else 1,
      start_time = start_time,
      end_time = Sys.time(),
      total_time_minutes = as.numeric(total_time)
    )
  ))
} # Close the run_parameter_importance function


#' @title Save Parameter Importance Results
#' @description
#' Saves parameter importance results to RDS and optionally CSV, with summary reporting
#'
#' @param parameter_results Output from run_parameter_importance()
#' @param output_dir Directory to save results (default: "output")
#' @param output_prefix Prefix for output filenames (default: "parameter_importance")
#' @param save_csv Whether to also save summary as CSV (default: TRUE)
#' @param verbose Print summary information (default: TRUE)
#'
#' @return Invisibly returns a list with filepaths of saved files
save_parameter_importance_results <- function(parameter_results,
                                              output_dir = "output",
                                              output_prefix = "parameter_importance",
                                              save_csv = TRUE,
                                              verbose = TRUE) {
  
  # ============================================================================
  # File saving
  # ============================================================================
  
  # Extract scenario name for filename
  scenario_name <- "unknown_scenario"
  if ("run_info" %in% names(parameter_results) && 
      "scenario" %in% names(parameter_results$run_info)) {
    scenario_name <- parameter_results$run_info$scenario
    # Clean scenario name for filename (remove spaces, special characters)
    scenario_name <- gsub("[^A-Za-z0-9-]", "_", scenario_name)
  }
  
  # Generate timestamp for unique filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Save main results as RDS file with scenario name
  rds_filename <- paste0(output_prefix, "_", scenario_name, "_", timestamp, ".rds")
  rds_filepath <- here::here(output_dir, rds_filename)
  saveRDS(parameter_results, rds_filepath)
  
  if (verbose) {
    cat("\n=== RESULTS SAVED ===\n")
    cat("RDS file:", rds_filepath, "\n")
  }
  
  # Save summary as CSV if requested and available
  csv_filepath <- NULL
  if (save_csv && !is.null(parameter_results$summary)) {
    csv_filename <- paste0(output_prefix, "_summary_", scenario_name, "_", timestamp, ".csv")
    csv_filepath <- here::here(output_dir, csv_filename)
    write.csv(parameter_results$summary, csv_filepath, row.names = FALSE)
    
    if (verbose) {
      cat("CSV summary:", csv_filepath, "\n")
    }
  }
  
  # ============================================================================
  # Extended summary reporting
  # ============================================================================
  
  # Print extended summary statistics if verbose mode enabled
  if (verbose && !is.null(parameter_results$summary)) {
    cat("\n=== EXTENDED SUMMARY ===\n")
    
    summary_df <- parameter_results$summary
    
    # Report temperature range across parameter sets
    cat("Peak temperature range:", 
        sprintf("%.2f - %.2f°C", 
                min(summary_df$peak_temperature, na.rm = TRUE),
                max(summary_df$peak_temperature, na.rm = TRUE)), "\n")
    
    # Report cost range across parameter sets
    cat("Total cost range:", 
        sprintf("%.1f - %.1f trillion $", 
                min(summary_df$total_cost, na.rm = TRUE),
                max(summary_df$total_cost, na.rm = TRUE)), "\n")
    
    # Identify parameter sets with extreme values
    min_cost_run <- summary_df$run_id[which.min(summary_df$total_cost)]
    max_cost_run <- summary_df$run_id[which.max(summary_df$total_cost)]
    min_temp_run <- summary_df$run_id[which.min(summary_df$peak_temperature)]
    max_temp_run <- summary_df$run_id[which.max(summary_df$peak_temperature)]
    
    cat("Lowest cost parameter set:", min_cost_run, "\n")
    cat("Highest cost parameter set:", max_cost_run, "\n")
    cat("Lowest peak temperature parameter set:", min_temp_run, "\n")
    cat("Highest peak temperature parameter set:", max_temp_run, "\n")
    
    # Report mitigation and CDR ranges
    cat("Total mitigation range:", 
        sprintf("%.1f - %.1f GtCO2", 
                min(summary_df$total_mitigation_units, na.rm = TRUE),
                max(summary_df$total_mitigation_units, na.rm = TRUE)), "\n")
    
    cat("Total CDR range:", 
        sprintf("%.1f - %.1f GtCO2", 
                min(summary_df$total_cdr_units, na.rm = TRUE),
                max(summary_df$total_cdr_units, na.rm = TRUE)), "\n")
    
    # Report convergence statistics
    n_converged <- sum(summary_df$converged, na.rm = TRUE)
    cat("Convergence rate:", 
        sprintf("%d / %d (%.1f%%)", 
                n_converged, 
                nrow(summary_df),
                100 * n_converged / nrow(summary_df)), "\n")
  }
  
  # Return filepaths invisibly for potential further use
  invisible(list(
    rds_filepath = rds_filepath,
    csv_filepath = csv_filepath
  ))
} # Close the save_parameter_importance_results function

#' @title Filter Parameter Importance Results
#' @description
#' Filters parameter importance results based on criteria, removing runs that don't
#' meet specified conditions (e.g., fail to return temperature to target)
#'
#' @param parameter_results Output from run_parameter_importance()
#' @param max_final_temperature Maximum acceptable final temperature (default: NULL, no filtering)
#' @param verbose Print information about filtered runs (default: TRUE)
#'
#' @return Filtered parameter_results list with same structure as input
#'
#' @examples
#' # Filter out runs that don't return to 1.5°C
#' cleaned_results <- filter_parameter_results(
#'   parameter_results = results,
#'   max_final_temperature = 1.6
#' )
filter_parameter_results <- function(parameter_results,
                                     max_final_temperature = NULL,
                                     verbose = TRUE) {
  
  # Create a copy of the results
  filtered_results <- parameter_results
  
  # Identify runs to remove based on final temperature
  runs_to_remove <- character(0)
  
  if (!is.null(max_final_temperature) && !is.null(parameter_results$summary)) {
    temp_filter <- parameter_results$summary$final_temperature > max_final_temperature
    runs_to_remove <- parameter_results$summary$run_id[temp_filter]
    
    if (length(runs_to_remove) > 0 && verbose) {
      cat("Identified", length(runs_to_remove), 
          "runs with final temperature >", max_final_temperature, "°C:\n")
      cat(paste(runs_to_remove, collapse = ", "), "\n\n")
    }
  }
  
  # Remove runs if any were identified
  if (length(runs_to_remove) > 0) {
    # Remove from successful_runs list
    filtered_results$successful_runs[runs_to_remove] <- NULL
    
    # Remove from summary dataframe
    filtered_results$summary <- filtered_results$summary %>%
      filter(!run_id %in% runs_to_remove)
    
    # Update run_info
    filtered_results$run_info$n_successful <- length(filtered_results$successful_runs)
    filtered_results$run_info$n_removed <- length(runs_to_remove)
    
    if (verbose) {
      cat("✓ Filtering complete:\n")
      cat("  Removed:", length(runs_to_remove), "runs\n")
      cat("  Remaining:", filtered_results$run_info$n_successful, "runs\n")
    }
  } else {
    if (verbose) {
      cat("No runs removed - all runs meet criteria\n")
    }
    filtered_results$run_info$n_removed <- 0
  }
  
  return(filtered_results)
} # Close the filter_parameter_results function



#' @title Export Parameter Importance Results for Sobol Analysis
#' @description
#' Exports parameter importance results in the format required for external
#' Sobol indices and coefficient of variation analysis. Creates separate CSV
#' files for each output variable (e.g., total cost, peak temperature) with
#' parameter values as columns and output variable as Y column.
#'
#' @param parameter_results Output from run_parameter_importance() or filtered results
#' @param output_dir Directory to save CSV files (default: "output")
#' @param output_variables Named list of functions to extract output variables.
#'   Default includes: totalcost, totalmitigation, totalcdr, totalabatement,
#'   peaktemperature, yearsabove1p5
#' @param parameter_columns Character vector of parameter column names to include.
#'   Default: c("tcre", "costmitiglow", "costremovlow", "costmitigpeak", 
#'              "costremovpeak", "econdampct", "discrate")
#' @param verbose Print progress information (default: TRUE)
#'
#' @return Invisibly returns list of created filepaths
#'
#' @examples
#' # Export with default settings
#' export_for_sobol_analysis(parameter_results)
#' 
#' # Export custom output variables
#' custom_outputs <- list(
#'   finaltemp = function(x) x$final_temperature,
#'   maxcdr = function(x) max(x$qty_remov)
#' )
#' export_for_sobol_analysis(parameter_results, output_variables = custom_outputs)
export_for_sobol_analysis <- function(parameter_results,
                                      output_dir = "output",
                                      output_variables = NULL,
                                      parameter_columns = c("tcre", "costmitiglow", "costremovlow",
                                                            "costmitigpeak", "costremovpeak", 
                                                            "econdampct", "discrate"),
                                      verbose = TRUE) {
  
  # ============================================================================
  # Input validation
  # ============================================================================
  
  # Extract successful_runs component
  if ("successful_runs" %in% names(parameter_results)) {
    successful_runs <- parameter_results$successful_runs
  } else {
    successful_runs <- parameter_results
  }
  
  # Validate that we have successful runs
  if (length(successful_runs) == 0) {
    stop("No successful runs found in parameter_results")
  }
  
  # Extract scenario name for filenames
  scenario_name <- "unknown_scenario"
  if ("run_info" %in% names(parameter_results) && 
      "scenario" %in% names(parameter_results$run_info)) {
    scenario_name <- parameter_results$run_info$scenario
  }
  scenario_clean <- gsub("[^A-Za-z0-9-]", "_", scenario_name)
  
  # ============================================================================
  # Setup default output variables if not provided
  # ============================================================================
  
  if (is.null(output_variables)) {
    output_variables <- list(
      totalcost = function(x) x$total_cost,
      totalmitigation = function(x) sum(x$qty_mitig),
      totalcdr = function(x) sum(x$qty_remov),
      totalabatement = function(x) sum(x$qty_mitig) + sum(x$qty_remov),
      peaktemperature = function(x) max(x$temperature_anomaly),
      yearsabove1p5 = function(x) sum(x$temperature_anomaly > 1.5)
    )
  }
  
  # ============================================================================
  # Build base data frame with parameter values
  # ============================================================================
  
  n_runs <- length(successful_runs)
  
  if (verbose) {
    cat("=== EXPORTING FOR SOBOL ANALYSIS ===\n")
    cat("Scenario:", scenario_name, "\n")
    cat("Number of runs:", n_runs, "\n")
    cat("Output variables:", length(output_variables), "\n")
  }
  
  # Initialise base data frame
  df_base <- data.frame(matrix(nrow = n_runs, ncol = length(parameter_columns)))
  colnames(df_base) <- parameter_columns
  
  # Parameter name mapping (script names -> CSV column names)
  param_mapping <- list(
    tcre = "tcre",
    costmitiglow = "cost_mitig_low",
    costremovlow = "cost_remov_low",
    costmitigpeak = "cost_mitig_peak",
    costremovpeak = "cost_remov_peak",
    econdampct = "econ_dam_pct",
    discrate = "disc_rate"
  )
  
  # Extract parameter values from each run
  for (i in 1:n_runs) {
    params <- successful_runs[[i]]$parameters
    
    for (col_name in parameter_columns) {
      if (col_name %in% names(param_mapping)) {
        param_name <- param_mapping[[col_name]]
        if (param_name %in% names(params)) {
          df_base[i, col_name] <- params[[param_name]]
        } else {
          warning("Parameter '", param_name, "' not found in run ", i)
        }
      }
    }
  }
  
  # ============================================================================
  # Create and save CSV files
  # ============================================================================
  
  # Generate timestamp for filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Store created filepaths
  created_files <- list()
  
  # Create CSV for each output variable
  for (output_name in names(output_variables)) {
    # Create data frame with parameters and output Y
    df_output <- df_base
    df_output$Y <- sapply(successful_runs, output_variables[[output_name]])
    
    # Construct filename with scenario
    filename <- paste0(output_name, "_", scenario_clean, "_", timestamp, ".csv")
    filepath <- here::here(output_dir, filename)
    
    # Save CSV
    write.csv(df_output, filepath, row.names = FALSE)
    
    # Store filepath
    created_files[[output_name]] <- filepath
    
    if (verbose) {
      cat("Created:", filename, "\n")
    }
  }
  
  if (verbose) {
    cat("\n✓ Exported", length(output_variables), "CSV files for Sobol analysis\n")
  }
  
  # Return created filepaths invisibly
  invisible(created_files)
} # Close the export_for_sobol_analysis function