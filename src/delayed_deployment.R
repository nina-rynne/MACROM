# ==============================================================================
# Delayed Deployment Analysis Functions
# 
# Part of: MACROM: An Optimal Control Model for Balancing Climate Change Abatement and Damage Trade-offs
# Authors: Nina Rynne, Michael Bode, Melanie Roberts, Ryan Heneghan
# Institution: Griffith University
# 
# Copyright (c) 2025 Nina Rynne
# Licensed under CC-BY-4.0 - see LICENSE file for details
# 
# Citation: If you use this code, please cite:
#   [Full citation of your paper]
# 
# Description:
# High-performance implementation of delayed deployment analysis with minimal
# function call overhead and pre-filtered data. Analyses the effects of
# delaying mitigation and/or CDR deployment on optimal control pathways,
# costs, and temperature outcomes across multiple SSP scenarios.
# 
# Version: 2.0.0
# Last updated: January 2026
# ==============================================================================

#' @title Fast Forward-Backward Sweep Algorithm
#' @description
#' Streamlined implementation of the forward-backward sweep method with pre-computed 
#' parameters and coefficients passed directly to eliminate redundant calculations and 
#' extractions. Uses the quadratic cost formulation to solve the optimal control problem
#' via iterative state-adjoint updates until convergence.
#'
#' @param baseline_emissions Pre-filtered baseline emissions vector (GtCO2/year)
#' @param baseline_gwp Pre-filtered baseline economic data vector (trillion $)
#' @param years Pre-filtered years vector
#' @param target_emissions Target cumulative emissions constraint (GtCO2)
#' @param terminal_adjoint Terminal adjoint value from shooting method
#' @param mitigation_start_idx Index where mitigation deployment begins
#' @param cdr_start_idx Index where carbon dioxide removal deployment begins
#' @param a Pre-computed mitigation quadratic coefficient a
#' @param b Pre-computed mitigation quadratic coefficient b
#' @param c Pre-computed CDR quadratic coefficient c
#' @param d Pre-computed CDR quadratic coefficient d
#' @param exp_temp_anom Temperature anomaly exponent parameter
#' @param clim_temp_init Initial climate temperature (°C)
#' @param tcre_scaled Pre-scaled TCRE value (tcre / 1000)
#' @param econ_dam_pct Economic damage percentage
#' @param discount_factors Pre-computed discount factors vector
#' @param exp_factors Pre-computed exponential factors vector
#' @param max_iterations Maximum iterations for convergence (default: 500)
#' @param tolerance Convergence tolerance threshold (default: 0.001)
#' @param epsilon Strict inequality parameter to prevent emissions going negative (default: 0.01)
#'
#' @return List containing:
#'   - converged: Logical indicating if algorithm converged
#'   - iterations: Number of iterations required
#'   - peak_temperature: Maximum temperature anomaly reached (°C)
#'   - final_emissions: Final cumulative emissions (GtCO2)
#'   - emission_gap: Difference between final and target emissions (GtCO2)
#'   - years_above_1p5: Number of years exceeding 1.5°C threshold
#'   - total_cdr_units: Total CDR deployed (GtCO2)
#'   - total_mitig_units: Total mitigation deployed (GtCO2)
#'   - total_cost: Total discounted cost (trillion $)
#'   - mitig_cost: Mitigation component of total cost (trillion $)
#'   - remov_cost: CDR component of total cost (trillion $)
#'   - temp_cost: Temperature damage component of total cost (trillion $)
fast_forward_backward_sweep <- function(baseline_emissions,
                                        baseline_gwp,
                                        years,
                                        target_emissions,
                                        terminal_adjoint,
                                        mitigation_start_idx,
                                        cdr_start_idx,
                                        a, b, c, d,
                                        exp_temp_anom,
                                        clim_temp_init,
                                        tcre_scaled,
                                        econ_dam_pct,
                                        discount_factors,
                                        exp_factors,
                                        max_iterations = 500,
                                        tolerance = 0.001,
                                        epsilon = 0.01) {
  
  # ============================================================================
  # Setup and initialization
  # ============================================================================
  
  # Pre-compute constants for algorithm
  n_years <- length(years)
  dt <- 1                    # Time step (years)
  update_weight <- 0.05      # Weight for smoothed control updates
  
  # Initialize state variables
  cumulative_emissions <- numeric(n_years)                    # Cumulative CO2 emissions (GtCO2)
  temperature_anomaly <- rep(clim_temp_init, n_years)         # Temperature anomaly above pre-industrial (°C)
  
  # Initialize control variables with starting guesses
  qty_mitig <- rep(1, n_years)                                # Mitigation quantity (GtCO2/year)
  qty_remov <- rep(0.1, n_years)                              # CDR quantity (GtCO2/year)
  
  # Initialize adjoint variable with terminal condition
  adjoint_var <- numeric(n_years)                             # Adjoint (co-state) variable
  adjoint_var[n_years] <- terminal_adjoint                    # Set terminal boundary condition
  
  # Pre-compute deployment masks for efficiency
  mitig_active <- seq_len(n_years) >= mitigation_start_idx    # Logical mask: mitigation active years
  cdr_active <- seq_len(n_years) >= cdr_start_idx             # Logical mask: CDR active years
  
  # ============================================================================
  # Helper function definition
  # ============================================================================
  
  # Calculate optimal control using quadratic formula solution
  # Solves the first-order optimality condition derived from Hamiltonian maximization
  # Returns the positive root of the quadratic equation when discriminant is non-negative
  calculate_optimal_control <- function(discriminant, adjoint_var_i, coeff_a, coeff_b) {
    if (discriminant >= 0 && adjoint_var_i >= 0) {
      # Use positive root of quadratic formula
      return((-coeff_b + sqrt(discriminant)) / (2 * coeff_a))
    } else {
      # Return zero control if solution is invalid
      return(0)
    }
  }
  
  # ============================================================================
  # Forward-backward sweep iteration loop
  # ============================================================================
  
  # Main iterative algorithm: alternates between forward state integration
  # and backward adjoint integration until controls converge
  for (iter in seq_len(max_iterations)) {
    
    # Store previous iteration values for convergence check
    prev_mitig <- qty_mitig
    prev_remov <- qty_remov
    prev_cumulative <- cumulative_emissions
    prev_adjoint <- adjoint_var
    
    # ========================================================================
    # Forward sweep: integrate state equations
    # ========================================================================
    
    # Calculate annual net emissions after mitigation and CDR
    annual_net <- baseline_emissions - qty_mitig - qty_remov
    
    # Integrate to get cumulative emissions over time
    cumulative_emissions[1] <- annual_net[1]
    for (i in 2:n_years) {
      cumulative_emissions[i] <- cumulative_emissions[i-1] + annual_net[i]
    }
    
    # Calculate temperature anomaly using TCRE relationship
    # Use pmax to enforce minimum temperature of 0.1°C for numerical stability
    temperature_anomaly <- pmax(clim_temp_init + cumulative_emissions * tcre_scaled, 0.1)
    
    # ========================================================================
    # Backward sweep: integrate adjoint equation
    # ========================================================================
    
    # Integrate adjoint (co-state) variable backward in time
    # The adjoint represents the marginal value of reducing cumulative emissions
    for (i in (n_years-1):1) {
      temp_base <- temperature_anomaly[i]
      
      # Calculate adjoint derivative from temperature damage gradient
      # This is -∂H/∂E where H is the Hamiltonian
      adjoint_derivative <- -(exp_temp_anom * baseline_gwp[i] * econ_dam_pct * 
                                tcre_scaled * discount_factors[i] * 
                                (temp_base^(exp_temp_anom - 1)))
      
      # Update adjoint with backward Euler step, checking for numerical issues
      if (is.finite(adjoint_derivative)) {
        adjoint_var[i] <- adjoint_var[i+1] - adjoint_derivative
      } else {
        # If derivative is non-finite, carry forward previous value
        adjoint_var[i] <- adjoint_var[i+1]
      }
    }
    
    # ========================================================================
    # Update controls using optimality conditions
    # ========================================================================
    
    # Initialize new control values (will be computed from first-order conditions)
    new_mitig <- numeric(n_years)
    new_remov <- numeric(n_years)
    
    # Update mitigation control using quadratic cost formulation
    # Only compute for years where mitigation is active and adjoint is positive
    mitig_indices <- which(mitig_active & adjoint_var > 0)
    
    if (length(mitig_indices) > 0) {
      for (idx in mitig_indices) {
        # Calculate discriminant for quadratic formula
        exp_term <- exp_factors[idx]
        discriminant <- b^2 + 4 * a * adjoint_var[idx] * exp_term
        
        # Solve first-order optimality condition using quadratic formula
        u_m_optimal <- calculate_optimal_control(discriminant, adjoint_var[idx], a, b)
        
        # Apply physical bounds to mitigation quantity
        if (u_m_optimal <= 0) {
          # No mitigation if optimal solution is non-positive
          new_mitig[idx] <- 0
        } else {
          # Cap mitigation to prevent emissions going negative (strict inequality)
          new_mitig[idx] <- min(u_m_optimal, baseline_emissions[idx] - epsilon)
        }
      }
    }
    
    # Update CDR control using quadratic cost formulation
    # Only compute for years where CDR is active and adjoint is positive
    cdr_indices <- which(cdr_active & adjoint_var > 0)
    
    if (length(cdr_indices) > 0) {
      for (idx in cdr_indices) {
        # Calculate discriminant for quadratic formula
        exp_term <- exp_factors[idx]
        discriminant <- d^2 + 4 * c * adjoint_var[idx] * exp_term
        
        # Solve first-order optimality condition using quadratic formula
        u_r_optimal <- calculate_optimal_control(discriminant, adjoint_var[idx], c, d)
        
        # Apply non-negativity constraint (CDR cannot be negative)
        new_remov[idx] <- max(0, u_r_optimal)
      }
    }
    
    # Apply smoothed updates to controls for numerical stability
    # Weighted average between new optimal values and previous iteration
    qty_mitig <- update_weight * new_mitig + (1 - update_weight) * prev_mitig
    qty_remov <- update_weight * new_remov + (1 - update_weight) * prev_remov
    
    # ========================================================================
    # Convergence check
    # ========================================================================
    
    # Calculate total change across all variables (L1 norm)
    total_change <- sum(abs(qty_mitig - prev_mitig)) + 
      sum(abs(qty_remov - prev_remov)) +
      sum(abs(cumulative_emissions - prev_cumulative)) +
      sum(abs(adjoint_var - prev_adjoint))
    
    # Exit loop if converged (total change below tolerance threshold)
    if (total_change < tolerance) {
      break
    }
  } # Close forward-backward sweep loop
  
  # ============================================================================
  # Calculate final metrics and costs
  # ============================================================================
  
  # Extract key outcome metrics from converged solution
  final_emissions <- cumulative_emissions[n_years]
  peak_temperature <- max(temperature_anomaly)
  years_above_1p5 <- sum(temperature_anomaly > 1.5)
  total_cdr_units <- sum(qty_remov)
  total_mitig_units <- sum(qty_mitig)
  
  # Calculate discounted costs using quadratic cost formulation
  # Cost = ∫(a*u³/3 + b*u²/2) * exp(-rt) dt for each control
  mitig_costs <- sum((a * qty_mitig^3 / 3 + b * qty_mitig^2 / 2) * discount_factors)
  remov_costs <- sum((c * qty_remov^3 / 3 + d * qty_remov^2 / 2) * discount_factors)
  temp_costs <- sum(baseline_gwp * econ_dam_pct * (temperature_anomaly^exp_temp_anom) * discount_factors)
  
  # Return comprehensive results list
  return(list(
    converged = (total_change < tolerance),
    iterations = iter,
    final_emissions = final_emissions,
    emission_gap = final_emissions - target_emissions,
    peak_temperature = peak_temperature,
    total_cost = mitig_costs + remov_costs + temp_costs,
    mitig_cost = mitig_costs,
    remov_cost = remov_costs,
    temp_cost = temp_costs,
    years_above_1p5 = years_above_1p5,
    total_cdr_units = total_cdr_units,
    total_mitig_units = total_mitig_units
  ))
} # Close fast_forward_backward_sweep function


#' @title Fast Shooting Method for Optimal Control
#' @description
#' Streamlined shooting method implementation that iteratively searches for the 
#' terminal adjoint value satisfying the cumulative emissions constraint. Uses 
#' bisection to find the adjoint value that makes the optimal control trajectory 
#' meet the target emissions. Pre-computed parameters and coefficients are passed 
#' directly to eliminate redundant calculations.
#'
#' @param baseline_emissions Pre-filtered baseline emissions vector (GtCO2/year)
#' @param baseline_gwp Pre-filtered baseline economic data vector (trillion $)
#' @param years Pre-filtered years vector
#' @param target_emissions Target cumulative emissions constraint (GtCO2)
#' @param mitigation_delay_years Years to delay mitigation deployment
#' @param cdr_delay_years Years to delay carbon dioxide removal deployment
#' @param a Pre-computed mitigation quadratic coefficient a
#' @param b Pre-computed mitigation quadratic coefficient b
#' @param c Pre-computed CDR quadratic coefficient c
#' @param d Pre-computed CDR quadratic coefficient d
#' @param exp_temp_anom Temperature anomaly exponent parameter
#' @param clim_temp_init Initial climate temperature (°C)
#' @param tcre_scaled Pre-scaled TCRE value (tcre / 1000)
#' @param econ_dam_pct Economic damage percentage
#' @param discount_factors Pre-computed discount factors vector
#' @param exp_factors Pre-computed exponential factors vector
#' @param max_iterations Maximum shooting iterations (default: 50)
#' @param tolerance Convergence tolerance for emission gap in GtCO2 (default: 3.0)
#'
#' @return List containing the best result found, with same structure as 
#'   fast_forward_backward_sweep() output
fast_shooting_method <- function(baseline_emissions,
                                 baseline_gwp, 
                                 years,
                                 target_emissions,
                                 mitigation_delay_years,
                                 cdr_delay_years,
                                 a, b, c, d,
                                 exp_temp_anom,
                                 clim_temp_init,
                                 tcre_scaled,
                                 econ_dam_pct,
                                 discount_factors,
                                 exp_factors,
                                 max_iterations = 50,
                                 tolerance = 3.0) {
  
  # ============================================================================
  # Setup and initialization
  # ============================================================================
  
  # Calculate deployment start indices from delay years
  # Index is delay + 1 because R uses 1-based indexing
  mitigation_start_idx <- mitigation_delay_years + 1
  cdr_start_idx <- cdr_delay_years + 1
  
  # Initialize bounds for bisection search on terminal adjoint value
  lambda_low <- -1
  lambda_high <- 50
  
  # ============================================================================
  # Test initial bounds
  # ============================================================================
  
  # Run forward-backward sweep with lower bound of terminal adjoint
  result_low <- fast_forward_backward_sweep(
    baseline_emissions, baseline_gwp, years, target_emissions, 
    lambda_low, mitigation_start_idx, cdr_start_idx,
    a, b, c, d, exp_temp_anom, clim_temp_init, tcre_scaled, 
    econ_dam_pct, discount_factors, exp_factors
  )
  
  # Run forward-backward sweep with upper bound of terminal adjoint
  result_high <- fast_forward_backward_sweep(
    baseline_emissions, baseline_gwp, years, target_emissions,
    lambda_high, mitigation_start_idx, cdr_start_idx,
    a, b, c, d, exp_temp_anom, clim_temp_init, tcre_scaled,
    econ_dam_pct, discount_factors, exp_factors
  )
  
  # Handle convergence failures at boundary points
  # If either bound fails to converge, return the result with smaller emission gap
  if (!result_low$converged || !result_high$converged) {
    best_result <- if (abs(result_low$emission_gap) < abs(result_high$emission_gap)) result_low else result_high
    
    return(list(
      converged = best_result$converged,
      iterations = best_result$iterations,
      final_emissions = best_result$final_emissions,
      emission_gap = best_result$emission_gap,
      peak_temperature = best_result$peak_temperature,
      total_cost = best_result$total_cost,
      mitig_cost = best_result$mitig_cost,
      remov_cost = best_result$remov_cost,
      temp_cost = best_result$temp_cost,
      years_above_1p5 = best_result$years_above_1p5,
      total_cdr_units = best_result$total_cdr_units,
      total_mitig_units = best_result$total_mitig_units
    ))
  }
  
  # Extract emission gaps from boundary results
  emission_gap_low <- result_low$emission_gap
  emission_gap_high <- result_high$emission_gap
  
  # ============================================================================
  # Adjust bounds if necessary
  # ============================================================================
  
  # Check if both bounds produce emission gaps with the same sign
  # This indicates the target is outside the current bracket and bounds need adjustment
  if (sign(emission_gap_low) == sign(emission_gap_high)) {
    if (emission_gap_low > 0 && emission_gap_high > 0) {
      # Both gaps positive: emissions too high, need lower adjoint values
      lambda_low <- lambda_low / 2
      lambda_high <- lambda_low + 1000
    } else {
      # Both gaps negative: emissions too low, need higher adjoint values
      lambda_high <- lambda_high * 2
      lambda_low <- lambda_high - 1000
    }
  }
  
  # ============================================================================
  # Main shooting iteration loop
  # ============================================================================
  
  # Track best result across all iterations in case convergence is not achieved
  best_result <- NULL
  best_gap <- Inf
  
  # Iterate using secant/bisection method to find terminal adjoint satisfying constraint
  for (iteration in seq_len(max_iterations)) {
    
    # Check for numerical stability: exit if gap difference becomes too small
    # This prevents division by near-zero values in secant method update
    if (abs(emission_gap_high - emission_gap_low) < 1e-10) break
    
    # Calculate new terminal adjoint using secant method
    # Linear interpolation between bounds based on emission gap values
    lambda_new <- lambda_low - emission_gap_low * (lambda_high - lambda_low) / 
      (emission_gap_high - emission_gap_low)
    
    # Run forward-backward sweep with new terminal adjoint guess
    result_new <- fast_forward_backward_sweep(
      baseline_emissions, baseline_gwp, years, target_emissions,
      lambda_new, mitigation_start_idx, cdr_start_idx,
      a, b, c, d, exp_temp_anom, clim_temp_init, tcre_scaled,
      econ_dam_pct, discount_factors, exp_factors
    )
    
    # Fallback to bisection if secant method produces non-convergent result
    # Bisection is more robust but slower, so use only when necessary
    if (!result_new$converged) {
      lambda_new <- (lambda_low + lambda_high) / 2
      result_new <- fast_forward_backward_sweep(
        baseline_emissions, baseline_gwp, years, target_emissions,
        lambda_new, mitigation_start_idx, cdr_start_idx,
        a, b, c, d, exp_temp_anom, clim_temp_init, tcre_scaled,
        econ_dam_pct, discount_factors, exp_factors
      )
    }
    
    # Extract emission gap from new result
    emission_gap_new <- result_new$emission_gap
    
    # Update best result tracker if this iteration is closer to target
    if (abs(emission_gap_new) < abs(best_gap)) {
      best_result <- result_new
      best_gap <- emission_gap_new
    }
    
    # Check if emission gap is within tolerance threshold
    # If converged, return immediately with successful result
    if (abs(emission_gap_new) <= tolerance) {
      return(list(
        converged = result_new$converged,
        iterations = result_new$iterations,
        final_emissions = result_new$final_emissions,
        emission_gap = result_new$emission_gap,
        peak_temperature = result_new$peak_temperature,
        total_cost = result_new$total_cost,
        mitig_cost = result_new$mitig_cost,
        remov_cost = result_new$remov_cost,
        temp_cost = result_new$temp_cost,
        years_above_1p5 = result_new$years_above_1p5,
        total_cdr_units = result_new$total_cdr_units,
        total_mitig_units = result_new$total_mitig_units
      ))
    }
    
    # Update bracket bounds based on sign of new emission gap
    # Replace the bound with same sign to maintain opposite signs at boundaries
    if (sign(emission_gap_new) == sign(emission_gap_low)) {
      # New result has same sign as lower bound, so replace lower bound
      lambda_low <- lambda_new
      emission_gap_low <- emission_gap_new
    } else {
      # New result has same sign as upper bound, so replace upper bound
      lambda_high <- lambda_new
      emission_gap_high <- emission_gap_new
    }
  } # Close shooting iteration loop
  
  # ============================================================================
  # Return best result if maximum iterations reached without convergence
  # ============================================================================
  
  return(list(
    converged = best_result$converged,
    iterations = best_result$iterations,
    final_emissions = best_result$final_emissions,
    emission_gap = best_result$emission_gap,
    peak_temperature = best_result$peak_temperature,
    total_cost = best_result$total_cost,
    mitig_cost = best_result$mitig_cost,
    remov_cost = best_result$remov_cost,
    temp_cost = best_result$temp_cost,
    years_above_1p5 = best_result$years_above_1p5,
    total_cdr_units = best_result$total_cdr_units,
    total_mitig_units = best_result$total_mitig_units
  ))
} # Close fast_shooting_method function


# ==============================================================================
# Helper Functions for Delayed Deployment Analysis
# ==============================================================================

#' @title Prepare Scenario Data for Analysis
#' @description
#' Filters and extracts time series data for a single scenario from emissions 
#' and economic data frames. Calculates time-dependent discount and exponential 
#' factors needed for the optimal control algorithm.
#'
#' @param emissions_df Data frame containing emissions data with Scenario and Year columns
#' @param economic_df Data frame containing economic (GWP) data with Scenario and Year columns
#' @param scenario Character string specifying which scenario to extract
#' @param disc_rate Discount rate for calculating present value factors
#'
#' @return List containing:
#'   - years: Vector of years
#'   - baseline_emissions: Emissions values (GtCO2/year)
#'   - baseline_gwp: Economic values (trillion $)
#'   - discount_factors: Present value discount factors exp(-r*t)
#'   - exp_factors: Exponential factors exp(r*t) for control calculations
#'   
#' @return NULL if no data found for the specified scenario
prepare_scenario_data <- function(emissions_df, 
                                  economic_df, 
                                  scenario, 
                                  disc_rate) {
  
  # ============================================================================
  # Filter data for specified scenario
  # ============================================================================
  
  # Extract emissions and economic data for this scenario
  emissions_scenario <- emissions_df %>% 
    filter(Scenario == scenario) %>% 
    arrange(Year)
  
  economic_scenario <- economic_df %>% 
    filter(Scenario == scenario) %>% 
    arrange(Year)
  
  # Check if data exists for this scenario
  if (nrow(emissions_scenario) == 0 || nrow(economic_scenario) == 0) {
    return(NULL)
  }
  
  # ============================================================================
  # Extract time series vectors
  # ============================================================================
  
  # Extract time series vectors for algorithm
  years <- emissions_scenario$Year
  baseline_emissions <- emissions_scenario$Value
  baseline_gwp <- economic_scenario$Value
  
  # ============================================================================
  # Calculate time-dependent factors
  # ============================================================================
  
  # Calculate years relative to start year for discounting
  years_rel <- years - years[1]
  
  # Calculate discount factors for net present value calculations
  discount_factors <- exp(-disc_rate * years_rel)
  
  # Calculate exponential factors used in control optimality conditions
  exp_factors <- exp(disc_rate * years_rel)
  
  # ============================================================================
  # Return prepared data
  # ============================================================================
  
  return(list(
    years = years,
    baseline_emissions = baseline_emissions,
    baseline_gwp = baseline_gwp,
    discount_factors = discount_factors,
    exp_factors = exp_factors
  ))
}


#' @title Package Delay Combination Result
#' @description
#' Converts shooting method result into standardized single-row data frame
#' with all relevant metrics and a feasibility flag.
#'
#' @param mitig_delay Mitigation delay years
#' @param cdr_delay CDR delay years
#' @param result List output from fast_shooting_method()
#' @param feasibility_tolerance Tolerance for emission gap to be considered feasible (GtCO2)
#'
#' @return Single-row data frame with all result columns
package_delay_result <- function(mitig_delay, 
                                 cdr_delay, 
                                 result, 
                                 feasibility_tolerance = 50) {
  
  data.frame(
    mitigation_delay = mitig_delay,
    cdr_delay = cdr_delay,
    peak_temperature = result$peak_temperature,
    final_emissions = result$final_emissions,
    emission_gap = result$emission_gap,
    total_cost = result$total_cost,
    mitig_cost = result$mitig_cost,
    remov_cost = result$remov_cost,
    temp_cost = result$temp_cost,
    converged = result$converged,
    feasible = abs(result$emission_gap) <= feasibility_tolerance,
    years_above_1p5 = result$years_above_1p5,
    total_cdr_units = result$total_cdr_units,
    total_mitig_units = result$total_mitig_units,
    stringsAsFactors = FALSE
  )
}


#' @title Package Error Result for Failed Combination
#' @description
#' Creates standardized error result data frame with NA values for all metrics
#' when a delay combination fails to solve.
#'
#' @param mitig_delay Mitigation delay years
#' @param cdr_delay CDR delay years
#'
#' @return Single-row data frame with NAs for all result columns
package_error_result <- function(mitig_delay, cdr_delay) {
  
  data.frame(
    mitigation_delay = mitig_delay, 
    cdr_delay = cdr_delay,
    peak_temperature = NA, 
    final_emissions = NA, 
    emission_gap = NA,
    total_cost = NA, 
    mitig_cost = NA, 
    remov_cost = NA, 
    temp_cost = NA,
    converged = FALSE, 
    feasible = FALSE, 
    years_above_1p5 = NA,
    total_cdr_units = NA, 
    total_mitig_units = NA, 
    stringsAsFactors = FALSE
  )
}


#' @title Calculate Summary Statistics
#' @description
#' Calculates convergence and feasibility statistics from delay combination results.
#'
#' @param results_df Data frame of results from delay combinations
#'
#' @return List containing:
#'   - n_combinations: Total number of combinations
#'   - n_converged: Number of converged solutions
#'   - n_feasible: Number of feasible solutions
#'   - convergence_rate: Proportion of converged solutions
#'   - feasibility_rate: Proportion of feasible solutions
#'   - feasible_results: Data frame filtered to only feasible results
calculate_summary_stats <- function(results_df) {
  
  # Count successful convergence and feasibility
  n_converged <- sum(results_df$converged, na.rm = TRUE)
  n_feasible <- sum(results_df$feasible, na.rm = TRUE)
  n_combinations <- nrow(results_df)
  
  # Filter to feasible results only (for further analysis)
  feasible_results <- results_df[results_df$feasible & !is.na(results_df$peak_temperature), ]
  
  # Return summary statistics
  return(list(
    n_combinations = n_combinations,
    n_converged = n_converged,
    n_feasible = n_feasible,
    convergence_rate = n_converged / n_combinations,
    feasibility_rate = n_feasible / n_combinations,
    feasible_results = feasible_results
  ))
}


# ==============================================================================
# Main Analysis Function
# ==============================================================================

#' @title Run Unified Delayed Deployment Analysis
#' @description
#' Comprehensive delayed deployment analysis that tests the impact of delaying 
#' mitigation and/or CDR deployment across multiple SSP scenarios. Creates a grid 
#' of delay combinations and uses the shooting method to find optimal control 
#' trajectories for each combination. Supports both single-scenario and multi-scenario 
#' analysis with automatic parallel processing for efficiency.
#'
#' @param parameter_df Single-row data frame containing all model parameters
#' @param emissions_df Data frame with emissions data for multiple scenarios
#' @param economic_df Data frame with economic (GWP) data for multiple scenarios
#' @param scenarios Character vector of scenario names to analyse, or "all" for all SSP 
#'   scenarios, or NULL for all available scenarios (default: NULL)
#' @param target_emissions Target cumulative emissions constraint in GtCO2 
#'   (default: uses co2_target_2100 from parameter_df)
#' @param max_delay_years Maximum delay to test in years (default: 40)
#' @param delay_step_size Step size for delay grid in years (default: 5)
#' @param use_parallel Enable parallel processing across scenarios (default: TRUE)
#' @param parallel_threshold Minimum number of combinations to trigger parallel 
#'   processing (default: 50)
#' @param n_cores Number of cores for parallel processing, or NULL for automatic 
#'   detection (default: NULL)
#' @param verbose Print progress information during execution (default: TRUE)
#'
#' @return List containing:
#'   - results_by_scenario: Named list of results for each scenario
#'   - combined_results: Data frame with all results combined across scenarios
#'   - summary_stats: Summary statistics including feasibility rates
#'   - run_info: Metadata about the analysis run
#'
#' @examples
#' # Run analysis for all SSP scenarios with 5-year delay increments
#' results <- run_unified_delayed_deployment(
#'   parameter_df = my_params,
#'   emissions_df = all_emissions,
#'   economic_df = all_economic,
#'   scenarios = "all",
#'   max_delay_years = 40,
#'   delay_step_size = 5
#' )
run_unified_delayed_deployment <- function(parameter_df,
                                           emissions_df,
                                           economic_df,
                                           scenarios = NULL,
                                           target_emissions = NULL,
                                           max_delay_years = 40,
                                           delay_step_size = 5,
                                           use_parallel = TRUE,
                                           parallel_threshold = 50,
                                           n_cores = NULL,
                                           verbose = TRUE,
                                           save_results = TRUE,      
                                           output_dir = "output",     
                                           output_prefix = "delayed_deployment") {
  
  # ============================================================================
  # Constants
  # ============================================================================
  
  # Feasibility tolerance: emission gap must be within this threshold (GtCO2)
  FEASIBILITY_TOLERANCE <- 50
  
  # ============================================================================
  # Input validation
  # ============================================================================
  
  # Validate parameter_df: must be non-empty data frame with exactly one row
  if (!is.data.frame(parameter_df)) {
    stop("parameter_df must be a data frame (received ", class(parameter_df)[1], ")")
  }
  if (nrow(parameter_df) != 1) {
    stop("parameter_df must contain exactly one row (received ", nrow(parameter_df), " rows)")
  }
  
  # Validate required cost parameters are present
  required_params <- c("cost_mitig_low", "cost_mitig_peak", "vol_mitig_low", "vol_mitig_peak",
                       "cost_remov_low", "cost_remov_peak", "vol_remov_low", "vol_remov_peak")
  missing_params <- setdiff(required_params, names(parameter_df))
  
  if (length(missing_params) > 0) {
    stop("Missing required parameters: ", paste(missing_params, collapse = ", "))
  }
  
  # Set target emissions from parameter_df if not provided
  if (is.null(target_emissions)) {
    target_emissions <- parameter_df$co2_target_2100
  }
  
  # ============================================================================
  # Determine scenarios to process
  # ============================================================================
  
  # Handle scenario selection: "all", NULL, or specific scenario names
  if (is.null(scenarios) || (length(scenarios) == 1 && scenarios == "all")) {
    # Use all standard SSP baseline scenarios
    ssp_scenarios <- c("SSP1-Baseline", "SSP2-Baseline", "SSP3-Baseline", 
                       "SSP4-Baseline", "SSP5-Baseline")
    available_scenarios <- unique(emissions_df$Scenario)
    scenarios_to_process <- intersect(ssp_scenarios, available_scenarios)
    multi_scenario_mode <- TRUE
  } else {
    # Use user-specified scenarios
    scenarios_to_process <- scenarios
    multi_scenario_mode <- length(scenarios) > 1
  }
  
  # Validate that at least one valid scenario is available
  if (length(scenarios_to_process) == 0) {
    stop("No valid scenarios found in emissions data. Available scenarios: ",
         paste(unique(emissions_df$Scenario), collapse = ", "))
  }
  
  # ============================================================================
  # Setup and initialization
  # ============================================================================
  
  # Extract cost and volume parameters from parameter data frame
  cost_mitig_low <- parameter_df$cost_mitig_low
  cost_mitig_peak <- parameter_df$cost_mitig_peak
  vol_mitig_low <- parameter_df$vol_mitig_low
  vol_mitig_peak <- parameter_df$vol_mitig_peak
  cost_remov_low <- parameter_df$cost_remov_low
  cost_remov_peak <- parameter_df$cost_remov_peak
  vol_remov_low <- parameter_df$vol_remov_low
  vol_remov_peak <- parameter_df$vol_remov_peak
  
  # Extract climate and economic parameters
  exp_temp_anom <- parameter_df$exp_temp_anom
  clim_temp_init <- parameter_df$clim_temp_init
  tcre <- parameter_df$tcre
  econ_dam_pct <- parameter_df$econ_dam_pct
  disc_rate <- parameter_df$disc_rate
  
  # Pre-compute quadratic cost function coefficients
  # For mitigation: Cost(u) = a*u³/3 + b*u²/2
  b <- (cost_mitig_peak - cost_mitig_low * vol_mitig_peak^2 / vol_mitig_low^2) / 
    (vol_mitig_peak - vol_mitig_peak^2 / vol_mitig_low)
  a <- (cost_mitig_low - vol_mitig_low * b) / vol_mitig_low^2
  
  # For CDR: Cost(u) = c*u³/3 + d*u²/2
  d <- (cost_remov_peak - cost_remov_low * vol_remov_peak^2 / vol_remov_low^2) / 
    (vol_remov_peak - vol_remov_peak^2 / vol_remov_low)
  c <- (cost_remov_low - vol_remov_low * d) / vol_remov_low^2
  
  # Scale TCRE for use in temperature calculations
  tcre_scaled <- tcre / 1000
  
  # ============================================================================
  # Create delay combination grid
  # ============================================================================
  
  # Generate sequence of delay values to test
  delay_sequence <- seq(0, max_delay_years, by = delay_step_size)
  
  # Calculate total combinations (mitigation delays × CDR delays)
  n_combinations_per_scenario <- length(delay_sequence)^2
  
  # ============================================================================
  # Print analysis configuration
  # ============================================================================
  
  # Report analysis setup if verbose mode enabled
  if (verbose) {
    if (multi_scenario_mode) {
      cat("=== MULTI-SCENARIO DELAYED DEPLOYMENT ANALYSIS ===\n")
      cat("Number of scenarios:", length(scenarios_to_process), "\n")
      cat("Scenarios:", paste(scenarios_to_process, collapse = ", "), "\n")
    } else {
      cat("=== SINGLE SCENARIO DELAYED DEPLOYMENT ANALYSIS ===\n")
      cat("Scenario:", scenarios_to_process[1], "\n")
    }
    cat("Target emissions:", target_emissions, "GtCO2\n")
    cat("Delay range:", 0, "-", max_delay_years, "years (step:", delay_step_size, ")\n")
    cat("Combinations per scenario:", n_combinations_per_scenario, "\n")
    cat("Using parallel processing:", use_parallel, "\n\n")
  }
  
  # Record start time for performance tracking
  overall_start_time <- Sys.time()
  
  # ============================================================================
  # Single scenario execution mode
  # ============================================================================
  
  # For single scenario, parallelize across delay combinations within the scenario
  if (!multi_scenario_mode) {
    scenario <- scenarios_to_process[1]
    
    # Prepare scenario data using helper function
    scenario_data <- prepare_scenario_data(emissions_df, economic_df, scenario, disc_rate)
    
    # Validate that data exists for this scenario
    if (is.null(scenario_data)) {
      stop("No data found for scenario: ", scenario)
    }
    
    # Extract prepared data
    years <- scenario_data$years
    baseline_emissions <- scenario_data$baseline_emissions
    baseline_gwp <- scenario_data$baseline_gwp
    discount_factors <- scenario_data$discount_factors
    exp_factors <- scenario_data$exp_factors
    
    # ========================================================================
    # Parallel processing setup for delay combinations
    # ========================================================================
    
    # Determine whether to use parallel processing based on threshold
    use_parallel_actual <- use_parallel && n_combinations_per_scenario >= parallel_threshold
    
    if (use_parallel_actual) {
      # Determine number of cores: use all except one if not specified
      n_cores_actual <- if (is.null(n_cores)) max(1, parallel::detectCores() - 1) else n_cores
      
      if (verbose) cat("Using parallel processing with", n_cores_actual, "cores\n")
      
      # Create and register parallel cluster
      cl <- makeCluster(n_cores_actual)
      registerDoParallel(cl)
      
      # Export required functions and variables to cluster workers
      clusterExport(cl, c("fast_shooting_method", "fast_forward_backward_sweep", 
                          "package_delay_result", "package_error_result",
                          "baseline_emissions", "baseline_gwp", "years", "target_emissions",
                          "a", "b", "c", "d", "exp_temp_anom", "clim_temp_init", 
                          "tcre_scaled", "econ_dam_pct", "discount_factors", "exp_factors",
                          "FEASIBILITY_TOLERANCE"), 
                    envir = environment())
      
      # Execute parallel loop over all delay combinations
      # Creates full grid: each mitigation delay paired with each CDR delay
      results_df <- foreach(mitig_delay = rep(delay_sequence, each = length(delay_sequence)),
                            cdr_delay = rep(delay_sequence, times = length(delay_sequence)),
                            .combine = rbind,
                            .packages = "dplyr") %dopar% {
                              
                              # Run shooting method with error handling for each combination
                              tryCatch({
                                # Solve optimal control problem for this delay combination
                                result <- fast_shooting_method(
                                  baseline_emissions, baseline_gwp, years, target_emissions,
                                  mitig_delay, cdr_delay, a, b, c, d, exp_temp_anom, 
                                  clim_temp_init, tcre_scaled, econ_dam_pct, discount_factors, exp_factors
                                )
                                
                                # Package results using helper function
                                package_delay_result(mitig_delay, cdr_delay, result, FEASIBILITY_TOLERANCE)
                                
                              }, error = function(e) {
                                # Return error result using helper function
                                package_error_result(mitig_delay, cdr_delay)
                              })
                            }
      
      # Shut down parallel cluster
      stopCluster(cl)
      
    } else {
      # ========================================================================
      # Serial processing alternative (when parallel not used)
      # ========================================================================
      
      # Initialize results list for serial execution
      results_list <- list()
      combination_count <- 0
      
      # Nested loop over all delay combinations
      for (mitig_delay in delay_sequence) {
        for (cdr_delay in delay_sequence) {
          combination_count <- combination_count + 1
          
          # Run shooting method with error handling for each combination
          results_list[[combination_count]] <- tryCatch({
            # Solve optimal control problem for this delay combination
            result <- fast_shooting_method(
              baseline_emissions, baseline_gwp, years, target_emissions,
              mitig_delay, cdr_delay, a, b, c, d, exp_temp_anom, 
              clim_temp_init, tcre_scaled, econ_dam_pct, discount_factors, exp_factors
            )
            
            # Package results using helper function
            package_delay_result(mitig_delay, cdr_delay, result, FEASIBILITY_TOLERANCE)
            
          }, error = function(e) {
            # Return error result using helper function
            package_error_result(mitig_delay, cdr_delay)
          })
        }
      }
      
      # Combine list of data frames into single data frame
      results_df <- do.call(rbind, results_list)
    } # Close parallel vs serial processing conditional
    
    # ========================================================================
    # Calculate summary statistics
    # ========================================================================
    
    # Calculate statistics using helper function
    summary <- calculate_summary_stats(results_df)
    
    # Calculate total runtime
    total_time <- difftime(Sys.time(), overall_start_time, units = "mins")
    
    # ========================================================================
    # Report summary statistics
    # ========================================================================
    
    if (verbose) {
      cat("\n=== ANALYSIS COMPLETE ===\n")
      cat("Total time:", sprintf("%.2f", total_time), "minutes\n")
      cat("Converged:", summary$n_converged, 
          sprintf("(%.1f%%)", 100 * summary$convergence_rate), "\n")
      cat("Feasible:", summary$n_feasible, 
          sprintf("(%.1f%%)", 100 * summary$feasibility_rate), "\n")
    }
    
    # ========================================================================
    # Return results for single scenario
    # ========================================================================
    
    return(list(
      results = results_df,
      summary_stats = list(
        n_combinations = summary$n_combinations,
        n_converged = summary$n_converged,
        n_feasible = summary$n_feasible,
        convergence_rate = summary$convergence_rate,
        feasibility_rate = summary$feasibility_rate
      ),
      feasible_results = summary$feasible_results,
      run_info = list(
        scenario = scenario,
        target_emissions = target_emissions,
        total_time_minutes = as.numeric(total_time)
      )
    ))
  } # Close single scenario mode
  
  # ============================================================================
  # Multi-scenario execution mode
  # ============================================================================
  
  # For multiple scenarios, parallelize across scenarios with optional nested 
  # parallelization within each scenario for delay combinations
  
  if (use_parallel) {
    # ========================================================================
    # Nested parallel processing setup
    # ========================================================================
    
    # Determine total available cores
    total_cores <- if (is.null(n_cores)) parallel::detectCores() - 1 else n_cores
    n_scenarios <- length(scenarios_to_process)
    
    # Distribute cores between scenario-level and combination-level parallelization
    if (n_scenarios >= total_cores) {
      # More scenarios than cores: allocate one core per scenario, no nested parallelism
      cores_per_scenario <- 1
      scenario_cores <- n_scenarios
    } else {
      # Fewer scenarios than cores: distribute remaining cores to combinations within scenarios
      cores_per_scenario <- max(1, floor(total_cores / n_scenarios))
      scenario_cores <- min(total_cores, n_scenarios)
    }
    
    # Report parallelization configuration if verbose mode enabled
    if (verbose) {
      cat("Using nested parallel processing:\n")
      cat("  - Scenarios in parallel:", scenario_cores, "cores\n")
      cat("  - Combinations per scenario:", cores_per_scenario, "cores\n")
      cat("  - Total combinations per scenario:", n_combinations_per_scenario, "\n\n")
    }
    
    # Create and register parallel cluster for scenario-level parallelization
    cl <- makeCluster(scenario_cores)
    registerDoParallel(cl)
    
    # Export required functions, data, and parameters to cluster workers
    clusterExport(cl, c("fast_shooting_method", "fast_forward_backward_sweep", 
                        "prepare_scenario_data", "package_delay_result", 
                        "package_error_result", "calculate_summary_stats",
                        "emissions_df", "economic_df",
                        "target_emissions", "delay_sequence", 
                        "a", "b", "c", "d", 
                        "exp_temp_anom", "clim_temp_init", "tcre_scaled", 
                        "econ_dam_pct", "disc_rate", 
                        "cores_per_scenario", "parallel_threshold",
                        "FEASIBILITY_TOLERANCE"), 
                  envir = environment())
    
    # ========================================================================
    # Execute parallel loop over scenarios
    # ========================================================================
    
    # Each iteration processes one scenario with all its delay combinations
    results_by_scenario <- foreach(scenario = scenarios_to_process,
                                   .packages = c("dplyr", "foreach", "doParallel", "parallel")) %dopar% {
                                     
                                     # Prepare scenario data using helper function
                                     scenario_data <- prepare_scenario_data(emissions_df, economic_df, scenario, disc_rate)
                                     
                                     # Check if data exists for this scenario
                                     if (is.null(scenario_data)) {
                                       return(list(scenario = scenario, error = "No data found"))
                                     }
                                     
                                     # Extract prepared data
                                     years <- scenario_data$years
                                     baseline_emissions <- scenario_data$baseline_emissions
                                     baseline_gwp <- scenario_data$baseline_gwp
                                     discount_factors <- scenario_data$discount_factors
                                     exp_factors <- scenario_data$exp_factors
                                     
                                     # ===========================================================
                                     # Nested parallelization for delay combinations (if enabled)
                                     # ===========================================================
                                     
                                     # Use nested parallel processing if sufficient cores and combinations
                                     if (cores_per_scenario > 1 && length(delay_sequence)^2 >= parallel_threshold) {
                                       
                                       # Create nested parallel cluster for this scenario's combinations
                                       nested_cl <- makeCluster(cores_per_scenario)
                                       registerDoParallel(nested_cl)
                                       
                                       # Export scenario-specific variables to nested cluster workers
                                       clusterExport(nested_cl, c("fast_shooting_method", "fast_forward_backward_sweep",
                                                                  "package_delay_result", "package_error_result",
                                                                  "baseline_emissions", "baseline_gwp", "years", "target_emissions",
                                                                  "a", "b", "c", "d", "exp_temp_anom", "clim_temp_init", 
                                                                  "tcre_scaled", "econ_dam_pct", "discount_factors", "exp_factors",
                                                                  "FEASIBILITY_TOLERANCE"), 
                                                     envir = environment())
                                       
                                       # Execute nested parallel loop over all delay combinations
                                       scenario_results_df <- foreach(mitig_delay = rep(delay_sequence, each = length(delay_sequence)),
                                                                      cdr_delay = rep(delay_sequence, times = length(delay_sequence)),
                                                                      .combine = rbind,
                                                                      .packages = "dplyr") %dopar% {
                                                                        
                                                                        # Run shooting method with error handling
                                                                        tryCatch({
                                                                          # Solve optimal control problem for this delay combination
                                                                          result <- fast_shooting_method(
                                                                            baseline_emissions, baseline_gwp, years, target_emissions,
                                                                            mitig_delay, cdr_delay, a, b, c, d, exp_temp_anom, 
                                                                            clim_temp_init, tcre_scaled, econ_dam_pct, discount_factors, exp_factors
                                                                          )
                                                                          
                                                                          # Package results using helper function
                                                                          package_delay_result(mitig_delay, cdr_delay, result, FEASIBILITY_TOLERANCE)
                                                                          
                                                                        }, error = function(e) {
                                                                          # Return error result using helper function
                                                                          package_error_result(mitig_delay, cdr_delay)
                                                                        })
                                                                      }
                                       
                                       # Shut down nested cluster
                                       stopCluster(nested_cl)
                                       
                                     } else {
                                       # ===========================================================
                                       # Serial processing alternative for combinations
                                       # ===========================================================
                                       
                                       # Process delay combinations serially within this scenario
                                       # (scenario itself still runs in parallel with other scenarios)
                                       scenario_results <- list()
                                       combination_count <- 0
                                       
                                       # Nested loop over all delay combinations
                                       for (mitig_delay in delay_sequence) {
                                         for (cdr_delay in delay_sequence) {
                                           combination_count <- combination_count + 1
                                           
                                           # Run shooting method with error handling for each combination
                                           scenario_results[[combination_count]] <- tryCatch({
                                             # Solve optimal control problem for this delay combination
                                             result <- fast_shooting_method(
                                               baseline_emissions, baseline_gwp, years, target_emissions,
                                               mitig_delay, cdr_delay, a, b, c, d, exp_temp_anom, 
                                               clim_temp_init, tcre_scaled, econ_dam_pct, discount_factors, exp_factors
                                             )
                                             
                                             # Package results using helper function
                                             package_delay_result(mitig_delay, cdr_delay, result, FEASIBILITY_TOLERANCE)
                                             
                                           }, error = function(e) {
                                             # Return error result using helper function
                                             package_error_result(mitig_delay, cdr_delay)
                                           })
                                         }
                                       }
                                       
                                       # Combine list of data frames into single data frame
                                       scenario_results_df <- do.call(rbind, scenario_results)
                                       
                                     } # Close nested parallelization vs serial conditional
                                     
                                     # Calculate summary statistics for this scenario using helper function
                                     summary <- calculate_summary_stats(scenario_results_df)
                                     
                                     # Return results for this scenario
                                     list(
                                       results = scenario_results_df,
                                       summary_stats = list(
                                         n_combinations = summary$n_combinations,
                                         n_converged = summary$n_converged,
                                         n_feasible = summary$n_feasible,
                                         convergence_rate = summary$convergence_rate,
                                         feasibility_rate = summary$feasibility_rate
                                       ),
                                       feasible_results = summary$feasible_results,
                                       run_info = list(
                                         scenario = scenario,
                                         target_emissions = target_emissions
                                       )
                                     )
                                   } # Close parallel foreach loop over scenarios
    
    # Shut down scenario-level parallel cluster
    stopCluster(cl)
    
    # Assign scenario names to results list for easy access
    names(results_by_scenario) <- scenarios_to_process
    
  } else {
    # ========================================================================
    # Serial processing alternative for scenarios
    # ========================================================================
    
    # Process scenarios serially (not in parallel)
    # This mode is used when use_parallel = FALSE
    results_by_scenario <- list()
    
    # Loop through each scenario sequentially
    for (i in seq_along(scenarios_to_process)) {
      scenario <- scenarios_to_process[i]
      
      # Report progress if verbose mode enabled
      if (verbose) {
        cat("Processing scenario", i, "of", length(scenarios_to_process), ":", scenario, "\n")
      }
      
      # Prepare scenario data using helper function
      scenario_data <- prepare_scenario_data(emissions_df, economic_df, scenario, disc_rate)
      
      # Check if data exists for this scenario
      if (is.null(scenario_data)) {
        if (verbose) cat("  ✗ No data found\n")
        next  # Skip to next scenario
      }
      
      # Extract prepared data
      years <- scenario_data$years
      baseline_emissions <- scenario_data$baseline_emissions
      baseline_gwp <- scenario_data$baseline_gwp
      discount_factors <- scenario_data$discount_factors
      exp_factors <- scenario_data$exp_factors
      
      # Process all delay combinations for this scenario
      scenario_results <- list()
      combination_count <- 0
      
      # Nested loop over all delay combinations
      for (mitig_delay in delay_sequence) {
        for (cdr_delay in delay_sequence) {
          combination_count <- combination_count + 1
          
          # Run shooting method with error handling for each combination
          scenario_results[[combination_count]] <- tryCatch({
            # Solve optimal control problem for this delay combination
            result <- fast_shooting_method(
              baseline_emissions, baseline_gwp, years, target_emissions,
              mitig_delay, cdr_delay, a, b, c, d, exp_temp_anom, 
              clim_temp_init, tcre_scaled, econ_dam_pct, discount_factors, exp_factors
            )
            
            # Package results using helper function
            package_delay_result(mitig_delay, cdr_delay, result, FEASIBILITY_TOLERANCE)
            
          }, error = function(e) {
            # Return error result using helper function
            package_error_result(mitig_delay, cdr_delay)
          })
        }
      }
      
      # Combine list of data frames into single data frame
      scenario_results_df <- do.call(rbind, scenario_results)
      
      # Calculate summary statistics for this scenario using helper function
      summary <- calculate_summary_stats(scenario_results_df)
      
      # Store results for this scenario
      results_by_scenario[[scenario]] <- list(
        results = scenario_results_df,
        summary_stats = list(
          n_combinations = summary$n_combinations,
          n_converged = summary$n_converged,
          n_feasible = summary$n_feasible,
          convergence_rate = summary$convergence_rate,
          feasibility_rate = summary$feasibility_rate
        ),
        feasible_results = summary$feasible_results,
        run_info = list(
          scenario = scenario,
          target_emissions = target_emissions
        )
      )
      
      # Report completion if verbose mode enabled
      if (verbose) cat("  ✓ Completed:", summary$n_feasible, "feasible combinations\n")
    } # Close loop over scenarios
    
  } # Close parallel vs serial processing for multi-scenario mode
  
  # ============================================================================
  # Filter and combine results across scenarios
  # ============================================================================
  
  # Remove any scenarios that returned errors
  successful_results <- results_by_scenario[!sapply(results_by_scenario, function(x) !is.null(x$error))]
  results_by_scenario <- successful_results
  
  # Combine all scenario results into single data frame
  successful_scenarios <- names(results_by_scenario)
  combined_results <- map_dfr(successful_scenarios, function(scenario_name) {
    scenario_data <- results_by_scenario[[scenario_name]]$results
    scenario_data$scenario <- scenario_name                                      # Add full scenario name
    scenario_data$scenario_short <- gsub("SSP([0-9])-Baseline", "SSP\\1", scenario_name)  # Add abbreviated name
    return(scenario_data)
  })
  
  # Calculate total runtime
  total_runtime <- difftime(Sys.time(), overall_start_time, units = "mins")
  
  # ============================================================================
  # Report final summary
  # ============================================================================
  
  if (verbose) {
    cat("\n=== MULTI-SCENARIO ANALYSIS COMPLETE ===\n")
    cat("Total runtime:", sprintf("%.1f", total_runtime), "minutes\n")
    cat("Successful scenarios:", length(successful_scenarios), "of", length(scenarios_to_process), "\n")
    cat("Total feasible combinations:", sum(combined_results$feasible, na.rm = TRUE), "\n")
  }
  
  # ============================================================================
  # Prepare results list
  # ============================================================================
  
  results_list <- list(
    results_by_scenario = results_by_scenario,
    combined_results = combined_results,
    summary_stats = list(
      n_scenarios_successful = length(successful_scenarios),
      total_combinations = nrow(combined_results),
      total_feasible = sum(combined_results$feasible, na.rm = TRUE),
      overall_feasibility_rate = sum(combined_results$feasible, na.rm = TRUE) / nrow(combined_results)
    ),
    run_info = list(
      total_runtime_minutes = as.numeric(total_runtime),
      scenarios_processed = successful_scenarios
    )
  )
  
  # ============================================================================
  # Save results if requested
  # ============================================================================
  
  if (save_results) {
    # Generate timestamp for unique filenames
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save main results as RDS file
    rds_filename <- paste0(output_prefix, "_", timestamp, ".rds")
    rds_filepath <- file.path(output_dir, rds_filename)
    
    # Save combined results as CSV
    csv_filename <- paste0(output_prefix, "_", timestamp, ".csv")
    csv_filepath <- file.path(output_dir, csv_filename)
    
    # Add saved file paths to run_info (NOW both paths exist)
    results_list$run_info$saved_files <- list(rds = rds_filepath, csv = csv_filepath)
    
    saveRDS(results_list, rds_filepath)
    
    if (verbose) {
      cat("\n=== RESULTS SAVED ===\n")
      cat("RDS file:", rds_filepath, "\n")
    }
    
    write.csv(combined_results, csv_filepath, row.names = FALSE)
    
    if (verbose) {
      cat("CSV file:", csv_filepath, "\n")
    }
  }
  
  # ============================================================================
  # Return comprehensive results
  # ============================================================================
  
  return(results_list)
} # Close run_unified_delayed_deployment function