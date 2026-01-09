#' @title Model Parameters for Climate Optimal Control
#' @description
#' This script contains functions for defining, managing, and sampling parameters
#' for the optimal control of temperature overshoot model. These functions handle
#' climate parameters, economic parameters, and numerical solver parameters.
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
#' 
#' @dependencies
#' Required packages: here, yaml
#' 
#' This module provides the following key functions:
#' - get_fixed_parameters: Loads the fixed parameters from a YAML configuration file.
#' - add_fixed_parameters: Takes a dataframe of variable parameters and adds columns 
#'                          for fixed parameters loaded from the YAML configuration
#'  

#' @title Get Fixed Parameters for Climate Model
#' @description
#' Loads fixed parameters from a YAML configuration file. This function reads the 
#' parameter_details.yml file from the project root directory and returns a flat list
#' of all parameters for use in the climate model.
#' @param config_file Path to the YAML configuration file, relative to project root
#' @param flatten Whether to flatten nested parameters into a single-level list (default: TRUE)
#' @param include_sections Which sections to include. Default is all sections.
#' @return List of fixed parameters
#' @examples
#' # Get all fixed parameters
#' fixed_params <- get_fixed_parameters()
#'
#' # Include only physical constants
#' fixed_params <- get_fixed_parameters(include_sections = "physical_constants")
#' 

get_fixed_parameters <- function(config_file = "parameter_details.yml") {

  # Build the full path to the config file using the here package
  config_path <- here::here(config_file)
  
  # Check if the file exists
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path, 
         "\nMake sure parameter_details.yml exists in your project root directory.")
  }
  
  # Read the YAML file
  params <- yaml::read_yaml(config_path)
  
  return(params)
}

#' @title Add Fixed Parameters to latin hypercube sample
#' @description
#' Takes a dataframe of variable parameters (typically from Latin Hypercube
#' Sampling) and adds columns for fixed parameters loaded from the YAML configuration
#' file. This creates a complete set of parameters for each model run.
#' @param param_df Dataframe of parameters (typically from generate_lhs_samples)
#' @param append_position Where to add fixed parameters: "start", "end" (default), or "alphabetical"
#' @param config_file Path to the YAML configuration file, relative to project root
#' @param exclude_params Names of fixed parameters to exclude (if any)
#' @return Dataframe with both variable and fixed parameters
#' @examples
#' # Add fixed parameters to the design
#' full_params <- add_fixed_parameters(variable_params)
#'

add_fixed_parameters_standard_target <- function(param_df,
                                 append_position = "end",
                                 config_file = "parameter_details.yml",
                                 exclude_params = NULL) {
  
  # Parameter validation
  if (!is.data.frame(param_df)) {
    stop("param_df must be a dataframe")
  }
  
  valid_positions <- c("start", "end", "alphabetical")
  if (!append_position %in% valid_positions) {
    stop("append_position must be one of: ", paste(valid_positions, collapse = ", "))
  }
  
  # Get the fixed parameters
  fixed_params <- get_fixed_parameters(config_file)
  
  # Exclude specified parameters if any
  if (!is.null(exclude_params)) {
    fixed_params <- fixed_params[!names(fixed_params) %in% exclude_params]
  }
  
  # Check for duplicate parameter names
  duplicate_params <- intersect(names(param_df), names(fixed_params))
  if (length(duplicate_params) > 0) {
    warning("The following parameters exist in both the input dataframe and fixed parameters: ",
            paste(duplicate_params, collapse = ", "), 
            ". Fixed values will override existing values.")
  }
  
  n_rows <- nrow(param_df)
  
  # Create the output dataframe based on the append_position
  result_df <- param_df
  
  # For each fixed parameter, add a column with the fixed value repeated
  for (param_name in names(fixed_params)) {
    result_df[[param_name]] <- rep(fixed_params[[param_name]], n_rows)
  }
  
  # Reorder columns based on append_position
  if (append_position == "start") {
    # Fixed parameters first, then variable parameters
    col_order <- c(names(fixed_params), setdiff(names(param_df), names(fixed_params)))
    result_df <- result_df[, col_order]
  } else if (append_position == "alphabetical") {
    # All parameters in alphabetical order
    result_df <- result_df[, sort(names(result_df))]
  }
  # For "end" position, no reordering needed as fixed params are already appended at the end
  
  # Add metadata attribute to track parameter source
  attr(result_df, "parameter_source") <- list(
    variable_params = names(param_df),
    fixed_params = names(fixed_params),
    config_file = config_file,
    timestamp = Sys.time()
  )
  
  return(result_df)
}

#' @title Create Parameter Dataframe from YAML with no LHS
#' @description
#' Creates a parameter dataframe with a single row directly from the YAML configuration file.
#' This function reads all parameters from the YAML file and converts them into a single-row
#' dataframe for use in model runs.
#' 
#' @param config_file       Path to the YAML configuration file, relative to project root
#' @param include_sections  Which sections to include from the YAML file (default: all sections)
#' @param exclude_params    Names of parameters to exclude (if any)
#' @param n_rows            Number of identical rows to include in the output dataframe (default: 1)
#' @return                  Dataframe with parameter values
#' @examples
#' # Create a parameter dataframe with all parameters
#' params_df <- create_params_dataframe()
#' 
#' # Create multiple identical rows (useful for compatibility with existing code)
#' params_df <- create_params_dataframe(n_rows = 10)
#' 
#' # Exclude specific parameters
#' params_df <- create_params_dataframe(exclude_params = c("trans_cond", "trans_low"))
#'

create_params_dataframe <- function(config_file = "parameter_details_fixed.yml",
                                    include_sections = NULL,
                                    exclude_params = NULL,
                                    n_rows = 1) {
  
  # Parameter validation
  if (!is.numeric(n_rows) || n_rows <= 0) {
    stop("n_rows must be a positive number")
  }
  
  # Build the full path to the config file using the here package
  config_path <- here::here(config_file)
  
  # Check if the file exists
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path, 
         "\nMake sure parameter_details_fixed.yml exists in your project root directory.")
  }
  
  # Read the YAML file
  yaml_params <- yaml::read_yaml(config_path)
  
  # Filter sections if specified
  if (!is.null(include_sections)) {
    if (!all(include_sections %in% names(yaml_params))) {
      missing_sections <- setdiff(include_sections, names(yaml_params))
      warning("The following requested sections were not found in the YAML file: ", 
              paste(missing_sections, collapse = ", "))
    }
    yaml_params <- yaml_params[names(yaml_params) %in% include_sections]
  }
  
  # Flatten the structure if it's nested
  params_list <- list()
  
  # Check if the YAML structure is nested
  if (any(sapply(yaml_params, is.list))) {
    # Nested structure - need to flatten
    for (section_name in names(yaml_params)) {
      section <- yaml_params[[section_name]]
      if (is.list(section)) {
        # Add each parameter from this section
        for (param_name in names(section)) {
          params_list[[param_name]] <- section[[param_name]]
        }
      } else {
        # Direct parameter
        params_list[[section_name]] <- section
      }
    }
  } else {
    # Already flat structure
    params_list <- yaml_params
  }
  
  # Exclude specified parameters if any
  if (!is.null(exclude_params)) {
    params_list <- params_list[!names(params_list) %in% exclude_params]
  }
  
  # Convert list to dataframe with one row
  params_df <- as.data.frame(params_list, stringsAsFactors = FALSE)
  
  # If n_rows > 1, replicate the row
  if (n_rows > 1) {
    params_df <- params_df[rep(1, n_rows), , drop = FALSE]
    rownames(params_df) <- NULL  # Reset row names
  }
  
  # Add metadata attribute to track parameter source
  attr(params_df, "parameter_source") <- list(
    all_params = names(params_list),
    config_file = config_file,
    timestamp = Sys.time()
  )
  
  return(params_df)
}

#' @title Add Fixed Parameters to latin hypercube sample
#' @description
#' Takes a dataframe of variable parameters (typically from Latin Hypercube
#' Sampling) and adds columns for fixed parameters loaded from the YAML configuration
#' file. This creates a complete set of parameters for each model run.
#' The function also calculates co2_target_2100 as (300/tcre) for each row,
#' overwriting the default value from the YAML file.
#' @param param_df Dataframe of parameters (typically from generate_lhs_samples)
#' @param append_position Where to add fixed parameters: "start", "end" (default), or "alphabetical"
#' @param config_file Path to the YAML configuration file, relative to project root
#' @param exclude_params Names of fixed parameters to exclude (if any)
#' @return Dataframe with both variable and fixed parameters
#' @examples
#' # Add fixed parameters to the design
#' full_params <- add_fixed_parameters(variable_params)
#'

add_fixed_parameters <- function(param_df,
                                 append_position = "end",
                                 config_file = "parameter_details.yml",
                                 exclude_params = NULL) {
  
  # Parameter validation
  if (!is.data.frame(param_df)) {
    stop("param_df must be a dataframe")
  }
  
  valid_positions <- c("start", "end", "alphabetical")
  if (!append_position %in% valid_positions) {
    stop("append_position must be one of: ", paste(valid_positions, collapse = ", "))
  }
  
  # Check if tcre column exists in the input dataframe
  if (!"tcre" %in% names(param_df)) {
    stop("The input dataframe must contain a 'tcre' column to calculate co2_target_2100")
  }
  
  # Get the fixed parameters
  fixed_params <- get_fixed_parameters(config_file)
  
  # Exclude specified parameters if any
  if (!is.null(exclude_params)) {
    fixed_params <- fixed_params[!names(fixed_params) %in% exclude_params]
  }
  
  # Check for duplicate parameter names
  duplicate_params <- intersect(names(param_df), names(fixed_params))
  if (length(duplicate_params) > 0) {
    warning("The following parameters exist in both the input dataframe and fixed parameters: ",
            paste(duplicate_params, collapse = ", "), 
            ". Fixed values will override existing values.")
  }
  
  n_rows <- nrow(param_df)
  
  # Create the output dataframe based on the append_position
  result_df <- param_df
  
  # For each fixed parameter, add a column with the fixed value repeated
  for (param_name in names(fixed_params)) {
    result_df[[param_name]] <- rep(fixed_params[[param_name]], n_rows)
  }
  
  # Calculate co2_target_2100 as (300/tcre) for each row, overwriting the YAML value
  result_df$co2_target_2100 <- 300 / result_df$tcre
  
  # Reorder columns based on append_position
  if (append_position == "start") {
    # Fixed parameters first, then variable parameters
    col_order <- c(names(fixed_params), setdiff(names(param_df), names(fixed_params)))
    result_df <- result_df[, col_order]
  } else if (append_position == "alphabetical") {
    # All parameters in alphabetical order
    result_df <- result_df[, sort(names(result_df))]
  }
  # For "end" position, no reordering needed as fixed params are already appended at the end
  
  # Add metadata attribute to track parameter source
  attr(result_df, "parameter_source") <- list(
    variable_params = names(param_df),
    fixed_params = names(fixed_params),
    calculated_params = "co2_target_2100",
    config_file = config_file,
    timestamp = Sys.time()
  )
  
  return(result_df)
}