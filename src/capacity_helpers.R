# ==============================================================================
# CDR Capacity Constraint Helper Functions
# ==============================================================================
# These functions generate time-varying capacity limits g(t) for CDR deployment
# to enforce realistic ramp-up constraints.
#
# Usage: Pass the generated function to optimal_control_solve via the
# cdr_capacity_function parameter.
#
# Author: Nina Rynne
# Date: February 2026
# ==============================================================================

#' @title Exponential Growth Capacity Function
#' @description
#' Creates a capacity function that grows exponentially from an initial value.
#' Realistic for early-stage technology deployment with learning and scaling.
#'
#' @param u0 Initial capacity in GtCO2/year (e.g., 0.1)
#' @param gamma Annual growth rate as a decimal (e.g., 0.05 = 5% per year)
#' @param t_start Year when growth begins (typically aligns with cdr_delay)
#'
#' @return Function with signature function(year) returning max capacity
#'
#' @examples
#' # 5% annual growth from 0.1 GtCO2/year starting in 2025
#' cap_fn <- make_exponential_capacity(u0 = 0.1, gamma = 0.05, t_start = 2025)
#' cap_fn(2025)  # 0.1
#' cap_fn(2030)  # ~0.128
#' cap_fn(2050)  # ~0.338

make_exponential_capacity <- function(u0, gamma, t_start) {
  function(year) {
    if (year < t_start) {
      return(0)
    }
    u0 * exp(gamma * (year - t_start))
  }
}


#' @title Logistic (S-Curve) Capacity Function
#' @description
#' Creates a capacity function following a logistic curve. Starts slow,
#' accelerates through the midpoint, then saturates at maximum capacity.
#' Realistic for technologies with initial barriers and eventual saturation.
#'
#' @param K Maximum capacity in GtCO2/year (saturation level)
#' @param r Steepness parameter (higher = faster transition)
#' @param t_mid Year at midpoint of transition (inflection point)
#'
#' @return Function with signature function(year) returning max capacity
#'
#' @examples
#' # S-curve reaching 10 GtCO2/year, centered at 2050
#' cap_fn <- make_logistic_capacity(K = 10, r = 0.1, t_mid = 2050)
#' cap_fn(2030)  # ~1.35
#' cap_fn(2050)  # 5.0 (midpoint)
#' cap_fn(2070)  # ~8.65

make_logistic_capacity <- function(K, r, t_mid) {
  function(year) {
    K / (1 + exp(-r * (year - t_mid)))
  }
}


#' @title Linear Ramp with Saturation Capacity Function
#' @description
#' Creates a capacity function that grows linearly until hitting a maximum cap.
#' Simple and interpretable for scenarios with known deployment rate limits.
#'
#' @param rho Linear growth rate in GtCO2/year² (e.g., 0.1 = add 0.1 GtCO2/year
#'   capacity each year)
#' @param K Maximum capacity in GtCO2/year (cap)
#' @param t_start Year when linear growth begins
#'
#' @return Function with signature function(year) returning max capacity
#'
#' @examples
#' # Linear growth of 0.2 GtCO2/year² capped at 15 GtCO2/year
#' cap_fn <- make_linear_capacity(rho = 0.2, K = 15, t_start = 2025)
#' cap_fn(2025)  # 0
#' cap_fn(2030)  # 1.0
#' cap_fn(2100)  # 15 (saturated)

make_linear_capacity <- function(rho, K, t_start) {
  function(year) {
    if (year < t_start) {
      return(0)
    }
    min(rho * (year - t_start), K)
  }
}


#' @title Piecewise Constant Capacity Function
#' @description
#' Creates a capacity function with discrete jumps at specified years.
#' Useful for policy scenarios with specific deployment targets by decade.
#'
#' @param breakpoints Named vector where names are years and values are
#'   capacities. Capacity is constant between breakpoints.
#'
#' @return Function with signature function(year) returning max capacity
#'
#' @examples
#' # Policy targets: 1 GtCO2/year by 2030, 5 by 2050, 10 by 2070
#' cap_fn <- make_piecewise_capacity(c("2025" = 0.5, "2030" = 1, 
#'                                     "2050" = 5, "2070" = 10))

make_piecewise_capacity <- function(breakpoints) {
  years <- as.numeric(names(breakpoints))
  capacities <- as.numeric(breakpoints)
  
  # Sort by year
  ord <- order(years)
  years <- years[ord]
  capacities <- capacities[ord]
  
  function(year) {
    if (year < years[1]) {
      return(0)
    }
    # Find the last breakpoint <= year
    idx <- max(which(years <= year))
    return(capacities[idx])
  }
}