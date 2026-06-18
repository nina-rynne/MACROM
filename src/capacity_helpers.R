# ==============================================================================
# Capacity Constraint Helper Functions
# ==============================================================================
# These functions generate time-varying capacity limits g(t) for mitigation
# and CDR deployment to enforce realistic ramp-up constraints.
#
# Available functional forms:
# - Exponential: Constant percentage growth (e.g., 5% per year)
# - Power Law: Fast early growth with diminishing returns (recommended for simple curves)
# - Logistic (from zero): Three-phase S-curve starting from realistic initial values (recommended)
# - Logistic (standard): S-curve parameterized by midpoint (harder to use)
# - Linear: Constant absolute growth rate with cap
# - Piecewise: Discrete policy targets by time period
#
# Usage: Pass the generated function to optimal_control_solve via the
# mitigation_capacity_function or cdr_capacity_function parameters.
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


#' @title Power Law Capacity Function
#' @description
#' Creates a capacity function following a power law with diminishing returns.
#' Growth is fast initially then slows naturally without requiring an explicit
#' saturation parameter. Realistic for technology deployment with learning
#' curves and resource constraints.
#'
#' The function is: g(t) = a * (t - t_start)^alpha for t >= t_start
#'
#' @param a Scale parameter controlling overall deployment ambition (GtCO2/year).
#'   Represents investment level and policy support. The capacity at year
#'   (t_start + 1) equals a. Higher values indicate more aggressive deployment.
#' @param alpha Exponent controlling the growth pattern (must be 0 < alpha < 1).
#'   Lower values (0.3-0.4) give strong diminishing returns - fast initial growth
#'   followed by rapid slowdown as easy wins are exhausted and constraints bind.
#'   Higher values (0.5-0.6) give more sustained growth - learning effects
#'   outpace resource constraints. Represents the balance between technology
#'   learning rates and physical/economic limits.
#' @param t_start Year when deployment begins
#'
#' @return Function with signature function(year) returning max capacity
#'
#' @examples
#' # Conservative: strong diminishing returns (low-hanging fruit exhausted quickly)
#' cap_fn <- make_power_capacity(a = 1.5, alpha = 0.35, t_start = 2020)
#' cap_fn(2030)  # ~2.4 GtCO2/year
#' cap_fn(2050)  # ~5.0 GtCO2/year
#' cap_fn(2100)  # ~12 GtCO2/year
#'
#' # Moderate: balanced learning and constraints
#' cap_fn <- make_power_capacity(a = 3, alpha = 0.4, t_start = 2020)
#' cap_fn(2030)  # ~6.0 GtCO2/year
#' cap_fn(2050)  # ~12 GtCO2/year
#' cap_fn(2100)  # ~30 GtCO2/year
#'
#' # Aggressive: sustained learning overcomes constraints
#' cap_fn <- make_power_capacity(a = 5, alpha = 0.45, t_start = 2020)
#' cap_fn(2030)  # ~12 GtCO2/year
#' cap_fn(2050)  # ~26 GtCO2/year
#' cap_fn(2100)  # ~70 GtCO2/year

make_power_capacity <- function(a, alpha, t_start) {
  # Input validation
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1 for realistic deployment curves with diminishing returns")
  }
  if (a <= 0) {
    stop("a must be positive (represents deployment scale)")
  }
  
  function(year) {
    if (year < t_start) {
      return(0)
    }
    a * (year - t_start)^alpha
  }
}


#' @title Logistic Growth Starting From Near-Zero
#' @description
#' Creates a logistic S-curve that starts from a specified small initial value
#' and grows to a maximum capacity. This version is easier to parameterize than
#' the standard logistic because you directly specify the starting capacity
#' instead of solving for it via t_mid.
#'
#' Captures three-phase deployment:
#' 1. Slow start (technology development, low investment)
#' 2. Rapid acceleration (technology matures, costs fall, massive scaling)
#' 3. Slowdown to saturation (easier solutions exhausted, physical/economic limits)
#'
#' The function is: g(t) = K / (1 + ((K/g_initial) - 1) * exp(-r*(t - t_start)))
#'
#' @param g_initial Starting capacity at t_start in GtCO2/year (e.g., 0.01 for
#'   current global CDR levels). This is what's actually deployed today.
#' @param K Maximum capacity (saturation level) in GtCO2/year. Represents
#'   physical, economic, or policy limits on deployment (e.g., 30 GtCO2/year).
#' @param r Growth rate controlling how fast the acceleration phase happens
#'   (typically 0.10-0.20). Higher r means faster transition from slow start
#'   to rapid growth. Common values: 0.10 (slow), 0.15 (moderate), 0.20 (fast).
#' @param t_start Year when deployment begins
#'
#' @return Function with signature function(year) returning max capacity
#'
#' @examples
#' # CDR deployment: current levels to realistic 2100 capacity
#' cap_fn <- make_logistic_from_zero(
#'   g_initial = 0.01,  # Current CDR ~0.002-0.01 GtCO2/year
#'   K = 30,            # Saturate at 30 GtCO2/year
#'   r = 0.15,          # Moderate acceleration
#'   t_start = 2020
#' )
#' cap_fn(2030)  # ~0.04 (slow start)
#' cap_fn(2050)  # ~5.4 (acceleration phase)
#' cap_fn(2070)  # ~27.4 (approaching saturation)
#'
#' # Mitigation: current efforts scaling to high capacity
#' cap_fn <- make_logistic_from_zero(
#'   g_initial = 5,     # Current global mitigation ~5 GtCO2/year
#'   K = 50,            # High deployment target
#'   r = 0.12,          # Moderate-slow acceleration
#'   t_start = 2020
#' )
#'
#' # Aggressive scenario: fast acceleration
#' cap_fn <- make_logistic_from_zero(
#'   g_initial = 0.1,
#'   K = 50,
#'   r = 0.20,          # Fast transition to scaling phase
#'   t_start = 2025
#' )

make_logistic_from_zero <- function(g_initial, K, r, t_start) {
  # Input validation
  if (g_initial <= 0) {
    stop("g_initial must be positive (starting capacity)")
  }
  if (K <= g_initial) {
    stop("K must be greater than g_initial (maximum capacity must exceed starting capacity)")
  }
  if (r <= 0) {
    stop("r must be positive (growth rate)")
  }
  
  # Force all parameters into the closure so their values are captured by the
  # enclosing environment rather than looked up by name at call time. This is
  # necessary when the returned function is serialised and sent to parallel
  # workers, which have a fresh environment where the caller's variables (e.g.
  # cdr_t_start) are not in scope.
  force(g_initial)
  force(K)
  force(r)
  force(t_start)
  
  # Pre-calculate the suppression factor for efficiency
  # This represents how much the initial capacity is suppressed relative to K
  suppression_factor <- (K / g_initial) - 1
  
  function(year) {
    if (year < t_start) {
      return(g_initial)
    }
    
    t <- year - t_start
    K / (1 + suppression_factor * exp(-r * t))
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

#' @title Make Zero Capacity Function
#' @description
#' Returns a capacity function that evaluates to zero for all years. Used with
#' use_mitigation_capacity_limit = TRUE to eliminate mitigation as a control
#' entirely, forcing the optimiser to rely solely on CDR. This is the cleanest
#' way to zero mitigation without modifying the solver code.
#'
#' @return Function with signature function(year) returning 0 for all inputs
#'
#' @examples
#' # Create a zero-capacity function
#' zero_fn <- make_zero_capacity()
#' zero_fn(2025)  # Returns 0
#' zero_fn(2060)  # Returns 0
#'
#' # Use in scenario comparison to eliminate mitigation
#' results <- run_scenario_comparison(
#'   ...,
#'   use_mitigation_capacity_limit = TRUE,
#'   mitigation_capacity_function  = make_zero_capacity()
#' )
make_zero_capacity <- function() {
  function(year) 0
}


#' @title Extract Capacity Growth Summary
#' @description
#' Pulls peak temperature and years above 1.5°C for every SSP × growth rate
#' combination from run_capacity_growth_comparison() output into a single tidy
#' data frame, then saves it as a timestamped CSV to the output/ directory.
#'
#' @param capacity_growth_results Output from run_capacity_growth_comparison()
#' @param verbose Print progress messages (default TRUE)
#'
#' @return Tidy data frame with columns: ssp, growth_rate, peak_temperature,
#'   years_above_1p5. Rows ordered by SSP then growth rate.
#'
#' @examples
#' capacity_summary_df <- extract_capacity_summary(capacity_growth_results)
extract_capacity_summary <- function(capacity_growth_results, verbose = TRUE) {

  clean_scenario_names <- function(x) gsub("-Baseline$", "", x)

  capacity_summary_df <- purrr::map_dfr(names(capacity_growth_results), function(rate) {
    summary <- capacity_growth_results[[rate]]$comparison_summary
    if (is.null(summary)) {
      warning(sprintf("No comparison_summary found for growth rate '%s'", rate))
      return(NULL)
    }
    data.frame(
      ssp              = clean_scenario_names(summary$scenario),
      growth_rate      = rate,
      peak_temperature = summary$peak_temperature,
      years_above_1p5  = summary$years_above_1p5,
      stringsAsFactors = FALSE
    )
  })

  growth_rate_order <- names(capacity_growth_results)
  ssp_order         <- sort(unique(capacity_summary_df$ssp))

  capacity_summary_df$growth_rate <- factor(capacity_summary_df$growth_rate,
                                             levels = growth_rate_order)
  capacity_summary_df$ssp         <- factor(capacity_summary_df$ssp,
                                             levels = ssp_order)

  capacity_summary_df <- capacity_summary_df[
    order(capacity_summary_df$ssp, capacity_summary_df$growth_rate), ]

  csv_filename <- paste0("capacity_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  utils::write.csv(capacity_summary_df, here::here("output", csv_filename), row.names = FALSE)

  if (verbose) cat("Capacity summary saved to output/", csv_filename, "\n")

  capacity_summary_df
}