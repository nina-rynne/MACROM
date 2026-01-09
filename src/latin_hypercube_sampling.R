#' @title Experimental Design for Climate Optimal Control
#' @description
#' This script contains functions for generating experimental designs for parameter
#' uncertainty analysis in the optimal control of temperature overshoot model. It includes
#' Latin Hypercube Sampling and other sampling methods for sensitivity analysis.
#'
#' @author Nina Rynne
#' @project Optimal Control of Temperature Overshoot
#' @date March 2025
#' @license MIT
#' @version 0.1.0
#'
#' @citation If published, include how to cite this work, e.g.:
#' Rynne, N. (2025). Optimal control methods for climate temperature overshoot.
#' Journal of Climate Modeling, XX(X), XXX-XXX. doi:XX.XXXX/XXXXX
#' 
#' @funding Include any grant or funding information here
#'
#' @collaborators List any collaborators or contributors
#'
#' @contact nina.rynne@griffithuni.edu.au
#'
#' @repository https://github.com/nina-rynne/optimal-control-model
#'
#' @references
#' - Carnell, R. (2006). lhs: Latin Hypercube Samples (p. 1.2.0) [Computer software]. 
#'   https://CRAN.R-project.org/package=lhs
#'
#' @dependencies
#' Required packages: dplyr, here, lhs
#' 
#' This module provides the following key functions:
#' - generate_lhs_samples: Calculates a latin hypercube sampling of the listed parameters
#' 

#' @title Latin Hypercube Sampling for Model Parameters
#' @description
#' Generates parameter sets using Latin Hypercube Sampling for model sensitivity analysis.
#' The function creates a specified number of parameter combinations where each parameter
#' range is divided into equal probability intervals.
#' 
#' @param n_samples         Number of sample points (default: 50)
#' @param seed              Random seed for reproducibility (default: NULL)
#' @param sampling_method   LHS method: "random", "optimum", or "genetic" (default: "random")
#' @param return_raw        Return raw 0-1 samples alongside scaled samples (default: FALSE)
#' @return                  Data frame of parameter combinations, or if return_raw=TRUE, a list
#'                           containing raw_samples and scaled_samples
#' @examples
#' # Generate 100 samples using default LHS and specific seed
#' lhs_results <- generate_lhs_samples(n_samples = 100, seed = 42)
#' # Generate samples using optimum sampling method and include raw data in return
#' lhs_results <- generate_lhs_samples(sampling_method = "optimum", return_raw = TRUE) 
#'                               

generate_lhs_samples <- function(n_samples = 50, 
                                 seed = NULL, 
                                 sampling_method = "random",
                                 return_raw = FALSE) {
  # Parameter validation
  if (!is.numeric(n_samples) || n_samples <= 0) {
    stop("n_samples must be a positive number")
  }
  
  valid_methods <- c("random", "optimum", "genetic")
  if (!sampling_method %in% valid_methods) {
    stop("sampling_method must be one of: ", paste(valid_methods, collapse=", "))
  }
  
  if (!is.logical(return_raw)) {
    stop("return_raw must be TRUE or FALSE")
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Define parameter ranges
  parameter_ranges <- list(
    tcre = c(min = 0.27, max = 0.63), # transient climate response to cumulative emissions degrees C
    cost_mitig_low = c(min = (0.005), max = (0.02)),
    cost_mitig_peak = c(min = (0.8), max = (1.5)), # cost of deploying mitigation (per Gt (*1e9), converted to trillions (/1e12))
    cost_remov_low = c(min = (0.005), max = (0.05)), # cost of deploying CDR trillion (per Gt (*1e9), converted to trillions (/1e12))
    cost_remov_peak = c(min = (1), max = (2.5)),
    econ_dam_pct = c(min = 0.01, max = 0.2), # proportion of GWP reduced by climate damage
    disc_rate = c(min = 0.01, max = 0.05) # future discounting
  )
  
  # Generate Latin Hypercube samples
  raw_lhs <- NULL
  if (sampling_method == "random") {
    raw_lhs <- randomLHS(n = n_samples, k = length(parameter_ranges))
  } else if (sampling_method == "optimum") {
    raw_lhs <- optimumLHS(n = n_samples, k = length(parameter_ranges))
  } else if (sampling_method == "genetic") {
    raw_lhs <- geneticLHS(n = n_samples, k = length(parameter_ranges))
  }
  
  # Scale the samples to the parameter ranges
  parameter_samples <- data.frame(raw_lhs)
  names(parameter_samples) <- names(parameter_ranges)
  
  # Transform the unit interval values to the actual parameter ranges
  for (param_name in names(parameter_ranges)) {
    param_range <- parameter_ranges[[param_name]]
    parameter_samples[[param_name]] <- param_range["min"] + 
      parameter_samples[[param_name]] * (param_range["max"] - param_range["min"])
  }
  
  if (return_raw) {
    return(list(raw_samples = raw_lhs, 
                scaled_samples = parameter_samples))
  } else {
    return(parameter_samples)
  }
}







