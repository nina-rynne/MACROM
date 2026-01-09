#' @title Optimal Control Core Algorithm for Climate Temperature Overshoot
#' @description
#' This module contains the core optimal control algorithm extracted from 
#' step6_strict_inequality_V2.R. It provides a clean, reusable implementation
#' of the forward-backward sweep with shooting method for climate optimal control.
#'
#' The core algorithm uses a STRICT INEQUALITY constraint for mitigation:
#' 0 <= u_m(t) < E(t) rather than 0 <= u_m(t) <= E(t)
#' This eliminates the need for KKT multipliers at the upper boundary.
#'
#' @author Nina Rynne
#' @project Optimal Control of Temperature Overshoot
#' @date June 2025
#' @license MIT
#' @version 0.1.0
#'
#' @citation If published, include how to cite this work, e.g.:
#' Rynne, N. (2025). Optimal control methods for climate temperature overshoot.
#' Journal of Climate Modeling, XX(X), XXX-XXX. doi:XX.XXXX/XXXXX
#'
#' @contact nina.rynne@griffithuni.edu.au
#' @repository https://github.com/nina-rynne/optimal-control-model
#'
#' @dependencies
#' Required packages: dplyr
#'
#' This module provides the following key functions:
#' - optimal_control_solve: Core forward-backward sweep algorithm
#' - optimal_control_shooting: Shooting method wrapper

#' @title Core Optimal Control Algorithm with Forward-Backward Sweep
#' @description
#' Solves the optimal control problem using forward-backward sweep iteration
#' with strict inequality constraints. This is the core mathematical algorithm
#' extracted from step6_strict_inequality_V2.R.
#'
#' MATHEMATICAL NOTE: Uses strict inequality u_m < E(t) to avoid numerical
#' instabilities at the boundary and eliminate the need for active set detection.
#'
#' @param parameter_df Single-row data frame containing all model parameters
#' @param emissions_df Data frame with emissions scenario data (filtered by scenario)
#' @param economic_df Data frame with economic scenario data (filtered by scenario)
#' @param target_emissions Target cumulative emissions constraint (GtCO2)
#' @param terminal_adjoint Terminal value for the adjoint variable (from shooting method)
#' @param mitigation_delay_years Years to delay mitigation deployment (default: 0)
#' @param cdr_delay_years Years to delay CDR deployment (default: 0)
#' @param max_iterations Maximum forward-backward sweep iterations (default: 500)
#' @param tolerance Convergence tolerance (default: 0.001)
#' @param update_weight Smoothing weight for control updates (default: 0.05)
#' @param epsilon Strict inequality parameter to ensure u_m < E(t) (default: 0.01)
#' @param verbose Print progress information (default: FALSE)
#'
#' @return List containing complete solution with:
#'   - Convergence info: converged, iterations, epsilon_used, constraint_violations
#'   - Time vectors: years, years_rel, n_years
#'   - State variables: cumulative_emissions, temperature_anomaly
#'   - Control variables: qty_mitig, qty_remov
#'   - Adjoint variable: adjoint_var
#'   - Baseline data: baseline_annual_emissions, baseline_annual_gwp
#'   - Cost components: mitig_costs_annual, remov_costs_annual, temp_costs_annual, total_costs_annual
#'   - Summary metrics: final_emissions, emission_gap, total_cost, final_temperature
#'   - Cost totals: mitig_cost, remov_cost, temp_cost
#'   - Diagnostics: max_mitig_fraction
#'
#' @examples
#' # Basic usage with shooting method terminal condition
#' result <- optimal_control_solve(
#'   parameter_df = my_params,
#'   emissions_df = emissions_data,
#'   economic_df = economic_data, 
#'   target_emissions = 650,
#'   terminal_adjoint = 1000
#' )
#'
#' # With delayed deployment
#' result <- optimal_control_solve(
#'   parameter_df = my_params,
#'   emissions_df = emissions_data,
#'   economic_df = economic_data,
#'   target_emissions = 650,
#'   terminal_adjoint = 1000,
#'   mitigation_delay_years = 10,
#'   cdr_delay_years = 5
#' )

optimal_control_solve <- function(parameter_df, 
                                  emissions_df, 
                                  economic_df, 
                                  target_emissions = 650,
                                  terminal_adjoint = -1,
                                  mitigation_delay_years = 0,
                                  cdr_delay_years = 0,
                                  max_iterations = 500,
                                  tolerance = 0.001,
                                  update_weight = 0.05,
                                  epsilon = 0.01,
                                  verbose = FALSE) {

# =============================================================================   
# ============================== SETUP ========================================
  
  # Input validation
  if (!is.data.frame(parameter_df) || nrow(parameter_df) != 1) {
    stop("parameter_df must be a single-row data frame")
  }
  
  if (!is.data.frame(emissions_df) || !is.data.frame(economic_df)) {
    stop("emissions_df and economic_df must be data frames")
  }
  
  # Data should already be filtered by scenario
  emissions_scenario <- emissions_df %>% dplyr::arrange(Year)
  economic_scenario <- economic_df %>% dplyr::arrange(Year)
  
  # Filter data by scenario
  # TESTING. ONLY USE WHEN NOT USING SHOOTING METHOD.
  # emissions_scenario <- emissions_df %>%
  #   dplyr::filter(Scenario == scenario) %>%
  #   dplyr::arrange(Year)
  # 
  # economic_scenario <- economic_df %>%
  #   dplyr::filter(Scenario == scenario) %>%
  #   dplyr::arrange(Year)

  if (nrow(emissions_scenario) == 0 || nrow(economic_scenario) == 0) {
    stop("No data found in emissions or economic data frames")
  }
  
  if (verbose) {
    cat("=== OPTIMAL CONTROL CORE ALGORITHM ===\n")
    cat("Using UPDATED COST FORMULATION with new optimal paths\n")
    cat("Target emissions:", target_emissions, "GtCO2\n")
    cat("Mitigation delay:", mitigation_delay_years, "years\n")
    cat("CDR delay:", cdr_delay_years, "years\n\n")
  }
  
  # Extract parameters from parameter_df
  #cost_mitig_unit <- parameter_df$cost_mitig_unit
  #cost_remov_unit <- parameter_df$cost_remov_unit
  exp_mitig <- parameter_df$exp_mitig
  exp_remov <- parameter_df$exp_remov
  exp_temp_anom <- parameter_df$exp_temp_anom
  clim_temp_init <- parameter_df$clim_temp_init
  tcre <- parameter_df$tcre
  econ_dam_pct <- parameter_df$econ_dam_pct
  disc_rate <- parameter_df$disc_rate
  cost_mitig_low <- parameter_df$cost_mitig_low
  cost_mitig_peak <- parameter_df$cost_mitig_peak
  vol_mitig_low <- parameter_df$vol_mitig_low
  vol_mitig_peak <- parameter_df$vol_mitig_peak
  cost_remov_low <- parameter_df$cost_remov_low
  cost_remov_peak <- parameter_df$cost_remov_peak
  vol_remov_low <- parameter_df$vol_remov_low
  vol_remov_peak <- parameter_df$vol_remov_peak
  
  # calculate integral cost function parameters
  # Mitigation
  b <- (cost_mitig_peak - cost_mitig_low*vol_mitig_peak^2/vol_mitig_low^2)/(vol_mitig_peak - vol_mitig_peak^2/vol_mitig_low)
  a <- 1/vol_mitig_low^2*(cost_mitig_low - vol_mitig_low*b)
  # CDR
  d <- (cost_remov_peak - cost_remov_low*vol_remov_peak^2/vol_remov_low^2)/(vol_remov_peak - vol_remov_peak^2/vol_remov_low)
  c <- 1/vol_remov_low^2*(cost_remov_low - vol_remov_low*d)
  
  # Time and data setup
  years <- emissions_scenario$Year
  years_rel <- years - min(years)
  n_years <- length(years)
  dt <- 1
  baseline_emissions <- emissions_scenario$Value
  baseline_gwp <- economic_scenario$Value
  
  # Calculate deployment start years
  mitigation_start_year <- min(years) + mitigation_delay_years
  cdr_start_year <- min(years) + cdr_delay_years
  
  # Initialize variables
  cumulative_emissions <- rep(0, n_years)
  temperature_anomaly <- rep(clim_temp_init, n_years)
  qty_mitig <- rep(1, n_years)    # Initial guess
  qty_remov <- rep(0.1, n_years)  # Initial guess
  adjoint_var <- rep(0, n_years)
  
  # Set terminal condition from shooting method
  adjoint_var[n_years] <- terminal_adjoint
  
  if (verbose) {
    cat("Initialized with terminal adjoint:", terminal_adjoint, "\n")
    cat("Starting forward-backward sweep iterations...\n")
  }
  
# =============================================================================   
# ============================= FBS LOOP ======================================
  for (iter in 1:max_iterations) {
    
    # Store previous values for convergence checking
    prev_mitig <- qty_mitig
    prev_remov <- qty_remov
    prev_cumulative <- cumulative_emissions
    prev_adjoint <- adjoint_var
    
# =============================================================================   
# ========================== FORWARD SWEEP ====================================
    for (i in 1:n_years) {
      if (i == 1) {
        # First year: start from zero and add first year's net emissions
        annual_net <- baseline_emissions[i] - qty_mitig[i] - qty_remov[i]
        cumulative_emissions[i] <- 0 + annual_net * dt
      } else {
        # Subsequent years: add to previous cumulative total
        annual_net <- baseline_emissions[i] - qty_mitig[i] - qty_remov[i]
        cumulative_emissions[i] <- cumulative_emissions[i-1] + annual_net * dt
      }
    }
    
    # Calculate temperature anomaly based on cumulative emissions
    temperature_anomaly <- clim_temp_init + (cumulative_emissions / 1000) * tcre
    
    # Trying temperature change since start of simulation as damage instead of 
    # including historical increases
    #temperature_anomaly <- (cumulative_emissions / 1000) * tcre
  
    
    # Ensure temperature stays positive for power calculations
    temperature_anomaly <- pmax(temperature_anomaly, 0.1)
    
    # Check for numerical issues
    if (any(!is.finite(temperature_anomaly))) {
      if (verbose) cat("ERROR: Non-finite temperature values at iteration", iter, "\n")
      break
    }
    
# =============================================================================   
# ========================== BACKWARD SWEEP ===================================
    
    # Keep terminal condition fixed (set by shooting method)
    # Integrate backward from terminal condition
    for (i in (n_years-1):1) {
      j <- i  
      
      temp_base <- temperature_anomaly[j]
      discount_factor <- exp(-disc_rate * years_rel[j])
      gwp_value <- baseline_gwp[j]
      
      # Adjoint derivative for quadratic damage function
      adjoint_derivative <- -(exp_temp_anom * gwp_value * econ_dam_pct * 
                                (tcre / 1000) * discount_factor * 
                                (temp_base^(exp_temp_anom - 1)))
      
      if (!is.finite(adjoint_derivative)) {
        adjoint_derivative <- 0
      }
      
      # Backward integration
      adjoint_var[j] <- adjoint_var[j+1] - adjoint_derivative * dt
    }
    
# =============================================================================   
# ========================== UPDATE CONTROLS ==================================
    
    new_mitig <- rep(0, n_years)
    new_remov <- rep(0, n_years)
    
    # Helper function to calculate optimal control from quadratic formula
    calculate_optimal_control <- function(discriminant, adjoint_var_i, exp_term, coeff_a, coeff_b) {
      if (discriminant >= 0 && adjoint_var_i >= 0) {
        return((-coeff_b + sqrt(discriminant)) / (2 * coeff_a))
      } else {
        return(0)
      }
    }
    
    # Helper function to apply bounds to mitigation
    apply_mitigation_bounds <- function(u_optimal, baseline_emiss_i, epsilon, max_growth) {
      if (u_optimal <= 0) {
        return(0)
      } else {
        # Apply both emissions constraint and growth constraint
        upper_bound <- min(baseline_emiss_i - epsilon, max_growth)
        return(min(u_optimal, upper_bound))
      }
    }
    
    # Helper function to apply bounds to removal
    apply_removal_bounds <- function(u_optimal, max_growth) {
      if (u_optimal <= 0) {
        return(0)
      } else {
        return(min(u_optimal, max_growth))
      }
    }
    
    for (i in 1:n_years) {

      current_year <- years[i]
      discount_factor <- exp(-disc_rate * years_rel[i])
      exp_term <- exp(disc_rate * years_rel[i])

      # MITIGATION CONTROL with NEW FORMULATION and DELAYED DEPLOYMENT
      if (current_year >= mitigation_start_year) {
        # Corrected optimal control formula from quadratic solution
        discriminant <- b^2 + 4 * a * adjoint_var[i] * exp_term

        if (discriminant >= 0 && adjoint_var[i] >= 0) {
          u_m_optimal <- (-b + sqrt(discriminant)) / (2 * a)

          # Apply bounds with STRICT INEQUALITY at upper bound
          if (u_m_optimal <= 0) {
            new_mitig[i] <- 0
          } else if (u_m_optimal >= baseline_emissions[i] - epsilon) {
            new_mitig[i] <- baseline_emissions[i] - epsilon
          } else {
            new_mitig[i] <- u_m_optimal
          }
        } else {
          new_mitig[i] <- 0
        }
      } else {
        # Before deployment date - mitigation forced to zero
        new_mitig[i] <- 0
      }

      # CDR CONTROL with NEW FORMULATION and DELAYED DEPLOYMENT
      if (current_year >= cdr_start_year) {
        # New optimal control formula:
        discriminant <- d^2 + 4 * c * adjoint_var[i] * exp_term

        if (discriminant > 0 && adjoint_var[i] > 0) {
          u_r_optimal <- (-d + sqrt(discriminant)) / (2 * c)

          if (u_r_optimal <= 0) {
            new_remov[i] <- 0
          } else {
            new_remov[i] <- u_r_optimal
          }
        } else {
          new_remov[i] <- 0
        }
      } else {
        # Before deployment date - CDR forced to zero
        new_remov[i] <- 0
      }
    }

    # Update with smoothing
    qty_mitig <- update_weight * new_mitig + (1 - update_weight) * prev_mitig
    qty_remov <- update_weight * new_remov + (1 - update_weight) * prev_remov
    
    # Force deployment delays to be exactly zero (smoothing affects this)
    for (i in 1:n_years) {
      current_year <- years[i]
      if (current_year < mitigation_start_year) {
        qty_mitig[i] <- 0
      }
      if (current_year < cdr_start_year) {
        qty_remov[i] <- 0
      }
    }
    
    # Check convergence
    mitig_change <- sum(abs(qty_mitig - prev_mitig))
    remov_change <- sum(abs(qty_remov - prev_remov))
    cumulative_change <- sum(abs(cumulative_emissions - prev_cumulative))
    adjoint_change <- sum(abs(adjoint_var - prev_adjoint))
    
    total_change <- mitig_change + remov_change + cumulative_change + adjoint_change
    
    # Progress reporting
    if (verbose && (iter %% 50 == 0 || iter <= 5)) {
      cat("Iteration", iter, ": Total change =", sprintf("%.6f", total_change), "\n")
    }
    
    if (total_change < tolerance) {
      if (verbose) cat("Converged at iteration", iter, "\n")
      break
    }
  }
  
  # Verify strict inequality constraint is satisfied
  constraint_violations <- sum(qty_mitig >= baseline_emissions)
  if (constraint_violations > 0 && verbose) {
    cat("WARNING: Strict inequality constraint violated at", constraint_violations, "time points\n")
  }

# =============================================================================   
# ========================== CALCULATE RESULTS ================================
  
  # Calculate comprehensive results
  final_emissions <- cumulative_emissions[n_years]
  final_temperature <- temperature_anomaly[n_years]
  discount_factors <- exp(-disc_rate * years_rel)
  
  # Annual costs using NEW FORMULATION
  mitig_costs_annual <- rep(0, n_years)
  remov_costs_annual <- rep(0, n_years)
  temp_costs_annual <- baseline_gwp * econ_dam_pct * (temperature_anomaly^exp_temp_anom) * discount_factors
  
  for (i in 1:n_years) {
    u_m <- qty_mitig[i]
    u_r <- qty_remov[i]
    
    # NEW COST FORMULATION
    mitig_costs_annual[i] <- (a * u_m^3 / 3 + b * u_m^2 / 2) * discount_factors[i]
    
    # NEW COST FORMULATION 
    remov_costs_annual[i] <- (c * u_r^3 / 3 + d * u_r^2 / 2) * discount_factors[i]
  }
  
  total_costs_annual <- mitig_costs_annual + remov_costs_annual + temp_costs_annual
  
  # Total costs (sums for backward compatibility)
  total_mitig_cost <- sum(mitig_costs_annual)
  total_remov_cost <- sum(remov_costs_annual)
  total_temp_cost <- sum(temp_costs_annual)
  total_cost <- sum(total_costs_annual)
  
  if (verbose) {
    cat("\n=== SOLUTION SUMMARY ===\n")
    cat("Final emissions:", round(final_emissions, 1), "GtCO2\n")
    cat("Emission gap:", round(final_emissions - target_emissions, 1), "GtCO2\n")
    cat("Final temperature:", round(final_temperature, 2), "°C\n")
    cat("Total cost:", sprintf("%.1f", total_cost), "trillion $\n")
    cat("Converged:", (total_change < tolerance), "\n")
    cat("Iterations:", iter, "\n")
  }

# =============================================================================   
# ========================== RETURN SOLUTION ==================================
  
  # Return comprehensive solution structure (matching original format)
  return(list(
    # Convergence and diagnostic information
    converged = (total_change < tolerance),
    iterations = iter,
    final_emissions = final_emissions,
    emission_gap = final_emissions - target_emissions,
    epsilon_used = epsilon,
    constraint_violations = constraint_violations,
    max_mitig_fraction = max(qty_mitig / baseline_emissions),
    
    # Time and deployment information
    years = years,
    years_rel = years_rel,
    n_years = n_years,
    mitigation_delay_years = mitigation_delay_years,
    cdr_delay_years = cdr_delay_years,
    mitigation_start_year = mitigation_start_year,
    cdr_start_year = cdr_start_year,
    
    # Baseline data
    baseline_annual_emissions = baseline_emissions,
    baseline_annual_gwp = baseline_gwp,
    
    # State variables
    cumulative_emissions = cumulative_emissions,
    temperature_anomaly = temperature_anomaly,
    
    # Control variables
    qty_mitig = qty_mitig,
    qty_remov = qty_remov,
    
    # Adjoint variable
    adjoint_var = adjoint_var,
    
    # Annual cost vectors (for detailed analysis)
    mitig_costs_annual = mitig_costs_annual,
    remov_costs_annual = remov_costs_annual,
    temp_costs_annual = temp_costs_annual,
    total_costs_annual = total_costs_annual,
    
    # Summary metrics
    final_temperature = final_temperature,
    total_cost = total_cost,
    mitig_cost = total_mitig_cost,
    remov_cost = total_remov_cost,
    temp_cost = total_temp_cost,
    
    # Parameters for reference
    parameters = parameter_df
  ))
}


#' @title Shooting Method for Optimal Control
#' @description
#' Implements the shooting method to find the correct terminal adjoint value
#' that satisfies the cumulative emissions constraint. This wraps the core
#' optimal_control_solve function with iterative boundary value solving.
#'
#' @param parameter_df Single-row data frame containing all model parameters
#' @param emissions_df Data frame with emissions scenario data
#' @param economic_df Data frame with economic scenario data  
#' @param scenario Scenario name to filter data (e.g., "SSP3-Baseline")
#' @param target_emissions Target cumulative emissions constraint (default: uses co2_target_2100 from parameter_df)
#' @param mitigation_delay_years Years to delay mitigation deployment (default: 0)
#' @param cdr_delay_years Years to delay CDR deployment (default: 0)
#' @param max_shooting_iterations Maximum shooting method iterations (default: 100)
#' @param shooting_tolerance Tolerance for emissions constraint (default: 1.0 GtCO2)
#' @param lambda_bounds Initial bounds for terminal adjoint [low, high] (default: c(0, 5000))
#' @param verbose Print progress information (default: FALSE)
#'
#' @return Complete solution list from optimal_control_solve, or best available if not converged
#'
#' @examples
#' # Basic shooting method
#' result <- optimal_control_shooting(
#'   parameter_df = my_params,
#'   emissions_df = emissions_data,
#'   economic_df = economic_data,
#'   scenario = "SSP3-Baseline"
#' )
#'
#' # With delayed deployment
#' result <- optimal_control_shooting(
#'   parameter_df = my_params,
#'   emissions_df = emissions_data,
#'   economic_df = economic_data,
#'   scenario = "SSP3-Baseline",
#'   mitigation_delay_years = 10,
#'   cdr_delay_years = 5
#' )
optimal_control_shooting <- function(parameter_df,
                                     emissions_df,
                                     economic_df, 
                                     scenario,
                                     target_emissions = NULL,
                                     mitigation_delay_years = 0,
                                     cdr_delay_years = 0,
                                     max_shooting_iterations = 1000,
                                     shooting_tolerance = 1.0,
                                     lambda_bounds = c(-50, 50),
                                     verbose = TRUE) {
  
  # Use target from parameter_df if not specified
  if (is.null(target_emissions)) {
    if ("co2_target_2100" %in% names(parameter_df)) {
      target_emissions <- parameter_df$co2_target_2100
    } else {
      stop("target_emissions must be specified or parameter_df must contain co2_target_2100")
    }
  }
  
  # Filter data by scenario
  emissions_scenario <- emissions_df %>%
    dplyr::filter(Scenario == scenario) %>%
    dplyr::arrange(Year)
  
  economic_scenario <- economic_df %>%
    dplyr::filter(Scenario == scenario) %>%
    dplyr::arrange(Year)
  
  if (nrow(emissions_scenario) == 0 || nrow(economic_scenario) == 0) {
    stop("No data found for scenario: ", scenario)
  }
  
  if (verbose) {
    cat("=== SHOOTING METHOD ===\n")
    cat("Scenario:", scenario, "\n")
    cat("Target emissions:", target_emissions, "GtCO2\n")
    cat("Lambda bounds: [", lambda_bounds[1], ",", lambda_bounds[2], "]\n\n")
  }
  
  # Initial bounds for terminal adjoint value
  lambda_low <- lambda_bounds[1]
  lambda_high <- lambda_bounds[2]
  
  # Test lower bound
  result_low <- optimal_control_solve(
    parameter_df = parameter_df,
    emissions_df = emissions_scenario,
    economic_df = economic_scenario,
    target_emissions = target_emissions,
    terminal_adjoint = lambda_low,
    mitigation_delay_years = mitigation_delay_years,
    cdr_delay_years = cdr_delay_years,
    verbose = FALSE
  )
  
  if (result_low$converged) {
    emission_gap_low <- result_low$emission_gap
    if (verbose) {
      cat("Lower bound result: final emissions =", round(result_low$final_emissions, 1), 
          ", gap =", round(emission_gap_low, 1), "\n")
    }
  } else {
    if (verbose) cat("WARNING: Lower bound did not converge\n")
    emission_gap_low <- Inf
  }
  
  # Test upper bound
  result_high <- optimal_control_solve(
    parameter_df = parameter_df,
    emissions_df = emissions_scenario,
    economic_df = economic_scenario,
    target_emissions = target_emissions,
    terminal_adjoint = lambda_high,
    mitigation_delay_years = mitigation_delay_years,
    cdr_delay_years = cdr_delay_years,
    verbose = TRUE
  )
  
  if (result_high$converged) {
    emission_gap_high <- result_high$emission_gap
    if (verbose) {
      cat("Upper bound result: final emissions =", round(result_high$final_emissions, 1), 
          ", gap =", round(emission_gap_high, 1), "\n")
    }
  } else {
    if (verbose) cat("WARNING: Upper bound did not converge\n")
    emission_gap_high <- -Inf
  }
  
  # Check if we have a valid bracket
  if (sign(emission_gap_low) == sign(emission_gap_high)) {
    if (verbose) {
      cat("WARNING: Initial bounds may not bracket the solution\n")
      cat("Adjusting bounds...\n")
    }
    
    # Try to find a better bracket
    if (emission_gap_low > 0 && emission_gap_high > 0) {
      # Both give emissions above target, need lower lambda
      lambda_low <- lambda_low / 2
      lambda_high <- lambda_low + 1000
    } else if (emission_gap_low < 0 && emission_gap_high < 0) {
      # Both give emissions below target, need higher lambda  
      lambda_high <- lambda_high * 2
      lambda_low <- lambda_high - 1000
    }
    
    if (verbose) cat("Adjusted bounds: [", lambda_low, ",", lambda_high, "]\n")
  }
  
  # Main shooting method loop (secant method)
  best_result <- NULL
  best_gap <- Inf
  
  for (iteration in 1:max_shooting_iterations) {
    
    # Progress reporting
    if (verbose && (iteration %% 10 == 1 || iteration <= 5)) {
      cat("Shooting iteration", iteration, ":\n")
    }
    
    # Secant method update
    if (abs(emission_gap_high - emission_gap_low) < 1e-10) {
      if (verbose) cat("WARNING: Numerical issues in secant method\n")
      break
    }
    
    # Calculate new lambda using secant method
    lambda_new <- lambda_low - emission_gap_low * (lambda_high - lambda_low) / 
      (emission_gap_high - emission_gap_low)
    
    if (verbose) cat("Trying lambda =", round(lambda_new, 4), "\n")
    
    # Evaluate at new point
    result_new <- optimal_control_solve(
      parameter_df = parameter_df,
      emissions_df = emissions_scenario,
      economic_df = economic_scenario,
      target_emissions = target_emissions,
      terminal_adjoint = lambda_new,
      mitigation_delay_years = mitigation_delay_years,
      cdr_delay_years = cdr_delay_years,
      verbose = TRUE
    )
    
    if (!result_new$converged) {
      if (verbose) cat("Inner loop did not converge, trying midpoint\n")
      lambda_new <- (lambda_low + lambda_high) / 2
      result_new <- optimal_control_solve(
        parameter_df = parameter_df,
        emissions_df = emissions_scenario,
        economic_df = economic_scenario,
        target_emissions = target_emissions,
        terminal_adjoint = lambda_new,
        mitigation_delay_years = mitigation_delay_years,
        cdr_delay_years = cdr_delay_years,
        verbose = TRUE
      )
    }
    
    emission_gap_new <- result_new$emission_gap
    
    if (verbose) {
      cat("Result: final emissions =", round(result_new$final_emissions, 1), 
          ", gap =", round(emission_gap_new, 1), "\n")
    }
    
    # Track best result
    if (abs(emission_gap_new) < abs(best_gap)) {
      best_result <- result_new
      best_gap <- emission_gap_new
    }
    
    # Check convergence
    if (abs(emission_gap_new) <= shooting_tolerance) {
      if (verbose) cat("Shooting method converged!\n")
      return(result_new)
    }
    
    # Update bounds for next iteration
    if (sign(emission_gap_new) == sign(emission_gap_low)) {
      lambda_low <- lambda_new
      emission_gap_low <- emission_gap_new
    } else {
      lambda_high <- lambda_new  
      emission_gap_high <- emission_gap_new
    }
  }
  
  if (verbose) {
    cat("Shooting method reached max iterations\n")
    cat("Best result found: gap =", round(best_gap, 1), "GtCO2\n")
  }
  
  return(best_result)
}

#' @title Run Optimal Control with Correct Delayed Deployment
#' @description
#' Runs optimal control with deployment delays using the mathematically correct
#' implementation from delayed_deployment.R. Returns full solution objects
#' including complete time series for visualization.
#'
#' @param parameter_df Single-row data frame containing all model parameters
#' @param emissions_df Data frame with emissions scenario data
#' @param economic_df Data frame with economic scenario data
#' @param scenario Scenario name to filter data (e.g., "SSP3-Baseline")
#' @param target_emissions Target cumulative emissions constraint (default: uses co2_target_2100 from parameter_df)
#' @param mitigation_delay_years Years to delay mitigation deployment (default: 0)
#' @param cdr_delay_years Years to delay CDR deployment (default: 0)
#' @param verbose Print progress information (default: FALSE)
#'
#' @return Complete solution list with full time series data
#'
#' @examples
#' # No controls (both delayed for full period)
#' result_no_controls <- run_optimal_control_with_delays(
#'   parameter_df = my_params,
#'   emissions_df = emissions_df,
#'   economic_df = economic_df,
#'   scenario = "SSP3-Baseline",
#'   mitigation_delay_years = 81,
#'   cdr_delay_years = 81
#' )
run_optimal_control_with_delays <- function(parameter_df,
                                            emissions_df,
                                            economic_df,
                                            scenario,
                                            target_emissions = NULL,
                                            mitigation_delay_years = 0,
                                            cdr_delay_years = 0,
                                            verbose = FALSE) {
  
  # Input validation
  if (!is.data.frame(parameter_df) || nrow(parameter_df) != 1) {
    stop("parameter_df must be a single-row data frame")
  }
  
  # Use target from parameter_df if not specified
  if (is.null(target_emissions)) {
    if ("co2_target_2100" %in% names(parameter_df)) {
      target_emissions <- parameter_df$co2_target_2100
    } else {
      stop("target_emissions must be specified or parameter_df must contain co2_target_2100")
    }
  }
  
  # Filter data by scenario (using correct delayed deployment approach)
  emissions_scenario <- emissions_df %>%
    dplyr::filter(Scenario == scenario) %>%
    dplyr::arrange(Year)
  
  economic_scenario <- economic_df %>%
    dplyr::filter(Scenario == scenario) %>%
    dplyr::arrange(Year)
  
  if (nrow(emissions_scenario) == 0 || nrow(economic_scenario) == 0) {
    stop("No data found for scenario: ", scenario)
  }
  
  if (verbose) {
    cat("=== OPTIMAL CONTROL WITH DELAYED DEPLOYMENT ===\n")
    cat("Scenario:", scenario, "\n")
    cat("Mitigation delay:", mitigation_delay_years, "years\n")
    cat("CDR delay:", cdr_delay_years, "years\n")
    cat("Target emissions:", target_emissions, "GtCO2\n\n")
  }
  
  # Extract parameters (following delayed_deployment.R approach)
  cost_mitig_unit <- parameter_df$cost_mitig_unit
  cost_remov_unit <- parameter_df$cost_remov_unit
  exp_mitig <- parameter_df$exp_mitig
  exp_remov <- parameter_df$exp_remov
  exp_temp_anom <- parameter_df$exp_temp_anom
  clim_temp_init <- parameter_df$clim_temp_init
  tcre <- parameter_df$tcre
  econ_dam_pct <- parameter_df$econ_dam_pct
  disc_rate <- parameter_df$disc_rate
  
  # Time and data setup
  years <- emissions_scenario$Year
  years_rel <- years - min(years)
  n_years <- length(years)
  dt <- 1
  baseline_emissions <- emissions_scenario$Value
  baseline_gwp <- economic_scenario$Value
  
  # Calculate deployment start years (following delayed_deployment.R)
  mitigation_start_year <- min(years) + mitigation_delay_years
  cdr_start_year <- min(years) + cdr_delay_years
  
  # Initialize variables
  cumulative_emissions <- rep(0, n_years)
  temperature_anomaly <- rep(clim_temp_init, n_years)
  qty_mitig <- rep(0, n_years)
  qty_remov <- rep(0, n_years)
  adjoint_var <- rep(0, n_years)
  
  # Shooting method parameters
  max_shooting_iterations <- 100
  shooting_tolerance <- 1.0
  lambda_bounds <- c(0, 5000)
  
  # Shooting method implementation (simplified for this specific use)
  lambda_low <- lambda_bounds[1]
  lambda_high <- lambda_bounds[2]
  
  # Test bounds and implement secant method
  for (shooting_iter in 1:max_shooting_iterations) {
    
    # Use midpoint for lambda
    lambda_current <- (lambda_low + lambda_high) / 2
    
    # Set terminal condition
    adjoint_var[n_years] <- lambda_current
    
    # Forward-backward sweep (following delayed_deployment.R structure)
    max_iterations <- 500
    tolerance <- 0.001
    update_weight <- 0.05
    epsilon <- 0.01
    
    for (iter in 1:max_iterations) {
      
      # Store previous values
      prev_mitig <- qty_mitig
      prev_remov <- qty_remov
      prev_cumulative <- cumulative_emissions
      prev_adjoint <- adjoint_var
      
      # === FORWARD SWEEP ===
      for (i in 1:n_years) {
        if (i == 1) {
          annual_net <- baseline_emissions[i] - qty_mitig[i] - qty_remov[i]
          cumulative_emissions[i] <- 0 + annual_net * dt
        } else {
          annual_net <- baseline_emissions[i] - qty_mitig[i] - qty_remov[i]
          cumulative_emissions[i] <- cumulative_emissions[i-1] + annual_net * dt
        }
      }
      
      # Calculate temperature anomaly
      temperature_anomaly <- clim_temp_init + (cumulative_emissions / 1000) * tcre
      temperature_anomaly <- pmax(temperature_anomaly, 0.1)
      
      # === BACKWARD SWEEP ===
      for (i in (n_years-1):1) {
        j <- i
        temp_base <- temperature_anomaly[j]
        discount_factor <- exp(-disc_rate * years_rel[j])
        gwp_value <- baseline_gwp[j]
        
        adjoint_derivative <- -(exp_temp_anom * gwp_value * econ_dam_pct * 
                                  (tcre / 1000) * discount_factor * 
                                  (temp_base^(exp_temp_anom - 1)))
        
        if (!is.finite(adjoint_derivative)) {
          adjoint_derivative <- 0
        }
        
        adjoint_var[j] <- adjoint_var[j+1] - adjoint_derivative * dt
      }
      
      # === UPDATE CONTROLS WITH CORRECT DELAY LOGIC ===
      new_mitig <- rep(0, n_years)
      new_remov <- rep(0, n_years)
      
      for (i in 1:n_years) {
        current_year <- years[i]
        discount_factor <- exp(-disc_rate * years_rel[i])
        
        # MITIGATION CONTROL with DELAYED DEPLOYMENT
        if (current_year >= mitigation_start_year) {
          # Calculate unconstrained optimal mitigation
          mitig_denominator <- exp_mitig * cost_mitig_unit * discount_factor
          if (mitig_denominator > 1e-10 && adjoint_var[i] > 0) {
            um_unconstrained <- adjoint_var[i] / mitig_denominator
          } else {
            um_unconstrained <- 0
          }
          
          # Apply bounds with strict inequality
          if (um_unconstrained <= 0) {
            new_mitig[i] <- 0
          } else if (um_unconstrained >= baseline_emissions[i] - epsilon) {
            new_mitig[i] <- baseline_emissions[i] - epsilon
          } else {
            new_mitig[i] <- um_unconstrained
          }
        } else {
          # Before deployment date - mitigation forced to zero
          new_mitig[i] <- 0
        }
        
        # CDR CONTROL with DELAYED DEPLOYMENT
        if (current_year >= cdr_start_year) {
          remov_denominator <- exp_remov * cost_remov_unit * discount_factor
          if (remov_denominator > 1e-10 && adjoint_var[i] > 0) {
            ur_unconstrained <- adjoint_var[i] / remov_denominator
          } else {
            ur_unconstrained <- 0
          }
          
          if (ur_unconstrained <= 0) {
            new_remov[i] <- 0
          } else {
            new_remov[i] <- ur_unconstrained
          }
        } else {
          # Before deployment date - CDR forced to zero
          new_remov[i] <- 0
        }
      }
      
      # Update with smoothing
      qty_mitig <- update_weight * new_mitig + (1 - update_weight) * prev_mitig
      qty_remov <- update_weight * new_remov + (1 - update_weight) * prev_remov
      
      # Check convergence
      mitig_change <- sum(abs(qty_mitig - prev_mitig))
      remov_change <- sum(abs(qty_remov - prev_remov))
      cumulative_change <- sum(abs(cumulative_emissions - prev_cumulative))
      adjoint_change <- sum(abs(adjoint_var - prev_adjoint))
      
      total_change <- mitig_change + remov_change + cumulative_change + adjoint_change
      
      if (total_change < tolerance) {
        break
      }
    }
    
    # Check emissions constraint
    final_emissions <- cumulative_emissions[n_years]
    emission_gap <- final_emissions - target_emissions
    
    if (abs(emission_gap) <= shooting_tolerance) {
      break
    }
    
    # Update bounds for shooting method
    if (emission_gap > 0) {
      lambda_low <- lambda_current
    } else {
      lambda_high <- lambda_current
    }
  }
  
  # Calculate costs and final metrics
  final_emissions <- cumulative_emissions[n_years]
  final_temperature <- temperature_anomaly[n_years]
  discount_factors <- exp(-disc_rate * years_rel)
  
  # Annual costs
  mitig_costs_annual <- cost_mitig_unit * qty_mitig^exp_mitig * discount_factors
  remov_costs_annual <- cost_remov_unit * qty_remov^exp_remov * discount_factors
  temp_costs_annual <- baseline_gwp * econ_dam_pct * (temperature_anomaly^exp_temp_anom) * discount_factors
  total_costs_annual <- mitig_costs_annual + remov_costs_annual + temp_costs_annual
  
  # Total costs
  total_mitig_cost <- sum(mitig_costs_annual)
  total_remov_cost <- sum(remov_costs_annual)
  total_temp_cost <- sum(temp_costs_annual)
  total_cost <- sum(total_costs_annual)
  
  if (verbose) {
    cat("Solution complete:\n")
    cat("Final emissions:", round(final_emissions, 1), "GtCO2\n")
    cat("Final temperature:", round(final_temperature, 2), "°C\n")
    cat("Total cost:", sprintf("%.1f", total_cost), "trillion $\n")
  }
  
  # Return complete solution structure (matching optimal_control_solve format)
  return(list(
    # Convergence and diagnostic information
    converged = (total_change < tolerance),
    iterations = iter,
    final_emissions = final_emissions,
    emission_gap = emission_gap,
    constraint_violations = 0,
    
    # Time and deployment information
    years = years,
    years_rel = years_rel,
    n_years = n_years,
    mitigation_delay_years = mitigation_delay_years,
    cdr_delay_years = cdr_delay_years,
    mitigation_start_year = mitigation_start_year,
    cdr_start_year = cdr_start_year,
    
    # Baseline data
    baseline_annual_emissions = baseline_emissions,
    baseline_annual_gwp = baseline_gwp,
    
    # State variables
    cumulative_emissions = cumulative_emissions,
    temperature_anomaly = temperature_anomaly,
    
    # Control variables
    qty_mitig = qty_mitig,
    qty_remov = qty_remov,
    
    # Adjoint variable
    adjoint_var = adjoint_var,
    
    # Annual cost vectors
    mitig_costs_annual = mitig_costs_annual,
    remov_costs_annual = remov_costs_annual,
    temp_costs_annual = temp_costs_annual,
    total_costs_annual = total_costs_annual,
    
    # Summary metrics
    final_temperature = final_temperature,
    total_cost = total_cost,
    mitig_cost = total_mitig_cost,
    remov_cost = total_remov_cost,
    temp_cost = total_temp_cost,
    
    # Parameters for reference
    parameters = parameter_df
  ))
}