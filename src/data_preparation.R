#' @title Data Preparation Functions for Climate Optimal Control
#' @description
#' This script contains functions for loading, processing and preparing data
#' for the optimal control of temperature overshoot model.
#'
#' @author Nina Rynne
#' @project Optimal Control of Temperature Overshoot
#' @date March 2025
#' @license MIT (or whatever license you plan to use)
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
#' Required packages: dplyr, tidyr, readr, here
#' 
#' This module provides the following key functions:
#' - import_ssp_emissions: Imports emissions scenario data from SSP sources
#' - import_ssp_economic: Imports economic data from SSP sources
#' - interpolate_ssp_emissions: Creates a complete emissions time series with interpolated values
#' - interpolate_ssp_economic: Creates a complete economic data time series with interpolated values
#'

# Define the data directory path relative to the project root
DATA_DIR <- "data"

#' @title Import SSP Emissions Data
#' @description
#' Imports emissions scenario data from CSV files for SSP scenarios.
#' The function expects a CSV file with a structure containing columns for
#' Model, Scenario, Region, Variable, Unit, and year columns (2005, 2010, etc.).
#'
#' @param file_name Name of the file within the data directory (e.g., "emissions_data.csv")
#' @param scenarios Vector of scenario names to filter (e.g., c("SSP2-Baseline", "SSP5-Baseline")).
#'        If NULL (default), imports all available scenarios.
#' @param regions Vector of regions to filter (default: "world")
#' @param variables Vector of emission variables to import (e.g., "Emissions|CO2").
#'        If NULL (default), imports all available variables.
#' @return A data frame (tibble) containing filtered emissions data in wide format with
#'         columns for Model, Scenario, Region, Variable, Unit, and year columns
#' @examples
#' # Import all emissions data for all scenarios
#' all_emissions <- import_ssp_emissions("emissions_data.csv")
#'
#' # Import only specific scenarios and variables
#' ssp_emissions <- import_ssp_emissions(
#'   file_name = "emissions_data.csv",
#'   scenarios = c("SSP3-Baseline", "SSP5-Baseline"),
#'   variables = "Emissions|CO2"
#'   )
#' 

import_ssp_emissions <- function(file_name, scenarios = NULL, regions = "world", variables = NULL) {
  # Construct path using here
  file_path <- here::here(DATA_DIR, file_name)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Emissions data file not found: ", file_path)
  }
  
  # Read the CSV file
  data <- readr::read_csv(file_path, show_col_types = FALSE)
  
  # Filter by scenarios if provided
  if (!is.null(scenarios)) {
    data <- data %>% dplyr::filter(Scenario %in% scenarios)
  }
  
  # Filter by regions if provided
  if (regions != "world") {
    data <- data %>% dplyr::filter(Region %in% regions)
  }
  
  # Filter by variables if provided
  if (!is.null(variables)) {
    data <- data %>% dplyr::filter(Variable %in% variables)
  }
  
  # Return the filtered data frame
  return(data)
}

#' @title Import SSP Economic Data
#' @description
#' Imports economic data (GDP, GWP) for SSP scenarios from CSV files.
#' The function expects a CSV file with a structure containing columns for
#' Model, Scenario, Region, Variable, Unit, and year columns (2005, 2010, etc.).
#' Economic variables typically include "GDP" measured in billion US$2005/yr.
#'
#' @param file_name Name of the file within the data directory (e.g., "econ_gwp.csv")
#' @param scenarios Vector of scenario names to filter (e.g., c("SSP2-Baseline", "SSP5-Baseline")).
#'        If NULL (default), imports all available scenarios.
#' @param regions Vector of regions to filter (default: "world")
#' @param variables Vector of economic variables to import (default: NULL)
#' @return A data frame (tibble) containing filtered economic data in wide format with
#'         columns for Model, Scenario, Region, Variable, Unit, and year columns
#' @examples
#' # Import all economic data for all scenarios
#' all_economic <- import_ssp_economic("econ_gwp.csv")
#'
#' # Import only specific scenarios with specific variable GDP
#' ssp_economic <- import_ssp_economic(
#'   file_name = "econ_gwp.csv",
#'   scenarios = c("SSP3-Baseline", "SSP5-Baseline"),
#'   variables = "GDP"
#' )
#'

import_ssp_economic <- function(file_name, scenarios = NULL, regions = "world", variables = NULL) {
  # Construct path using here
  file_path <- here::here(DATA_DIR, file_name)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Economic data file not found: ", file_path)
  }
  
  # Read the CSV file
  data <- readr::read_csv(file_path, show_col_types = FALSE)
  
  # Filter by scenarios if provided
  if (!is.null(scenarios)) {
    data <- data %>% dplyr::filter(Scenario %in% scenarios)
  }
  
  # Filter by regions if provided
  if (regions != "world") {
    data <- data %>% dplyr::filter(Region %in% regions)
  }
  
  # Filter by variables if provided
  if (!is.null(variables)) {
    data <- data %>% dplyr::filter(Variable %in% variables)
  }
  
  # Return the filtered data frame
  return(data)
}

#' @title Interpolate SSP Emissions Time Series
#' @description
#' Creates a complete emissions time series with specified time step by interpolating
#' between available data points in the wide-format SSP data.
#'
#' @param emissions_df Data frame containing emissions data in wide format with year columns
#' @param dt Time step size for interpolation (in years)
#' @param start_year Starting year for the interpolated series
#' @param end_year Ending year for the interpolated series
#' @return A data frame with complete emissions time series at specified time step,
#'         restructured to long format with columns for Scenario, Model, Region, Variable,
#'         Unit, Year, and Value
#' @examples
#' interpolated_emissions <- interpolate_ssp_emissions(
#'   emissions_df = ssp_emissions_data,
#'   dt = 1,
#'   start_year = 2020,
#'   end_year = 2100
#' )

old_interpolate_ssp_emissions <- function(emissions_df, dt = 1, start_year = 2025, end_year = 2100) {
  # Convert wide format to long format first to make interpolation easier
  emissions_long <- emissions_df %>%
    tidyr::pivot_longer(
      cols = where(is.numeric),  # All numeric columns (years)
      names_to = "Year",
      values_to = "Value"
    ) %>%
    dplyr::mutate(Year = as.numeric(Year))  # Ensure Year is numeric
  
  # Get unique combinations of Scenario, Model, Region, and Variable
  series_groups <- emissions_long %>%
    dplyr::select(Scenario, Model, Region, Variable, Unit) %>%
    dplyr::distinct()
  
  # Create an empty list to store results from each iteration
  result_list <- list()
  
  # Process each unique combination
  for (i in 1:nrow(series_groups)) {
    # Extract current group
    current <- series_groups[i, ]
    
    # Filter data for current group
    group_data <- emissions_long %>%
      dplyr::filter(
        Scenario == current$Scenario,
        Model == current$Model,
        Region == current$Region,
        Variable == current$Variable
      )
    
    # Create sequence of years for interpolation
    years_seq <- seq(start_year, end_year, by = dt)
    
    # Create a new dataframe with all years in the sequence
    new_years <- data.frame(
      Scenario = current$Scenario,
      Model = current$Model,
      Region = current$Region,
      Variable = current$Variable,
      Unit = current$Unit,
      Year = years_seq
    )
    
    # Merge with existing data to keep original values
    merged_data <- dplyr::full_join(
      new_years,
      group_data,
      by = c("Scenario", "Model", "Region", "Variable", "Unit", "Year")
    )
    
    # Sort by year
    merged_data <- merged_data %>% dplyr::arrange(Year)
    
    # Get available years and values (non-NA)
    available_years <- merged_data$Year[!is.na(merged_data$Value)]
    available_values <- merged_data$Value[!is.na(merged_data$Value)]
    
    # Only proceed with interpolation if we have at least 2 valid points
    if (length(available_values) >= 2) {
      # Perform linear interpolation
      interp_values <- approx(
        x = available_years,
        y = available_values,
        xout = merged_data$Year,
        method = "linear",
        rule = 2  # rule=2 means extrapolate beyond the data range
      )
      
      merged_data$Value <- interp_values$y
      
      # Store this result in the list
      result_list[[i]] <- merged_data
    } else if (length(available_values) == 1) {
      # If only one valid point, use that value for all years
      merged_data$Value <- available_values[1]
      
      # Store this result in the list
      result_list[[i]] <- merged_data
    }
    # If zero valid points, we don't add anything to the list
  }
  
  # Combine all results from the list using bind_rows
  # This is more reliable than repeatedly using rbind
  if (length(result_list) > 0) {
    interpolated_data <- dplyr::bind_rows(result_list)
    return(interpolated_data)
  } else {
    warning("No valid interpolated data was produced. Check input data structure.")
    return(NULL)
  }
}

#' @title Interpolate SSP Emissions Time Series
#' @description
#' Creates a complete emissions time series with specified time step by interpolating
#' between available data points in the wide-format SSP data. All emission values 
#' are adjusted by a multiplier parameter.
#'
#' @param emissions_df Data frame containing emissions data in wide format with year columns
#' @param dt Time step size for interpolation (in years)
#' @param start_year Starting year for the interpolated series
#' @param end_year Ending year for the interpolated series
#' @param clim_co2_remain Factor to multiply emission values by (default: 0.5)
#' @return A data frame with complete emissions time series at specified time step,
#'         restructured to long format with columns for Scenario, Model, Region, Variable,
#'         Unit, Year, and Value
#' @examples
#' interpolated_emissions <- interpolate_ssp_emissions(
#'   emissions_df = ssp_emissions_data,
#'   dt = 1,
#'   start_year = 2020,
#'   end_year = 2100,
#'   clim_co2_remain = 0.6
#' )

interpolate_ssp_emissions <- function(emissions_df, dt = 1, start_year = 2025, end_year = 2100, clim_co2_remain = 0.5) {
  # Convert wide format to long format first to make interpolation easier
  emissions_long <- emissions_df %>%
    tidyr::pivot_longer(
      cols = where(is.numeric),  # All numeric columns (years)
      names_to = "Year",
      values_to = "Value"
    ) %>%
    dplyr::mutate(Year = as.numeric(Year))  # Ensure Year is numeric
  
  # Get unique combinations of Scenario, Model, Region, and Variable
  series_groups <- emissions_long %>%
    dplyr::select(Scenario, Model, Region, Variable, Unit) %>%
    dplyr::distinct()
  
  # Create an empty list to store results from each iteration
  result_list <- list()
  
  # Process each unique combination
  for (i in 1:nrow(series_groups)) {
    # Extract current group
    current <- series_groups[i, ]
    
    # Filter data for current group
    group_data <- emissions_long %>%
      dplyr::filter(
        Scenario == current$Scenario,
        Model == current$Model,
        Region == current$Region,
        Variable == current$Variable
      )
    
    # Create sequence of years for interpolation
    years_seq <- seq(start_year, end_year, by = dt)
    
    # Create a new dataframe with all years in the sequence
    new_years <- data.frame(
      Scenario = current$Scenario,
      Model = current$Model,
      Region = current$Region,
      Variable = current$Variable,
      Unit = current$Unit,
      Year = years_seq
    )
    
    # Merge with existing data to keep original values
    merged_data <- dplyr::full_join(
      new_years,
      group_data,
      by = c("Scenario", "Model", "Region", "Variable", "Unit", "Year")
    )
    
    # Sort by year
    merged_data <- merged_data %>% dplyr::arrange(Year)
    
    # Get available years and values (non-NA)
    available_years <- merged_data$Year[!is.na(merged_data$Value)]
    available_values <- merged_data$Value[!is.na(merged_data$Value)]
    
    # Only proceed with interpolation if we have at least 2 valid points
    if (length(available_values) >= 2) {
      # Perform linear interpolation
      interp_values <- approx(
        x = available_years,
        y = available_values,
        xout = merged_data$Year,
        method = "linear",
        rule = 2  # rule=2 means extrapolate beyond the data range
      )
      
      merged_data$Value <- interp_values$y
      
      # Store this result in the list
      result_list[[i]] <- merged_data
    } else if (length(available_values) == 1) {
      # If only one valid point, use that value for all years
      merged_data$Value <- available_values[1]
      
      # Store this result in the list
      result_list[[i]] <- merged_data
    }
    # If zero valid points, we don't add anything to the list
  }
  
  # Combine all results from the list using bind_rows
  # This is more reliable than repeatedly using rbind
  if (length(result_list) > 0) {
    interpolated_data <- dplyr::bind_rows(result_list)
    
    # Apply the % emissions remaining after ocean/biosphere uptake
    interpolated_data <- interpolated_data %>%
      dplyr::mutate(Value = Value * clim_co2_remain)
    
    return(interpolated_data)
  } else {
    warning("No valid interpolated data was produced. Check input data structure.")
    return(NULL)
  }
}

#' @title Interpolate SSP Economic Time Series
#' @description
#' Creates a complete economic data time series with specified time step by interpolating
#' between available data points in the wide-format SSP economic data.
#'
#' @param economic_df Data frame containing economic data in wide format with year columns
#' @param dt Time step size for interpolation (in years)
#' @param start_year Starting year for the interpolated series
#' @param end_year Ending year for the interpolated series
#' @return A data frame with complete economic time series at specified time step,
#'         restructured to long format with columns for Scenario, Model, Region, Variable,
#'         Unit, Year, and Value
#' @examples
#' interpolated_economic <- interpolate_ssp_economic(
#'   economic_df = ssp_economic_data,
#'   dt = 1,
#'   start_year = 2020,
#'   end_year = 2100
#' )
#' 
  
interpolate_ssp_economic <- function(economic_df, dt = 1, start_year = 2025, end_year = 2100) {
  # Convert wide format to long format first to make interpolation easier
  economic_long <- economic_df %>%
    tidyr::pivot_longer(
      cols = where(is.numeric),  # All numeric columns (years)
      names_to = "Year",
      values_to = "Value"
    ) %>%
    dplyr::mutate(Year = as.numeric(Year))  # Ensure Year is numeric
  
  # Get unique combinations of Scenario, Model, Region, and Variable
  series_groups <- economic_long %>%
    dplyr::select(Scenario, Model, Region, Variable, Unit) %>%
    dplyr::distinct()
  
  # Create an empty list to store results from each iteration
  result_list <- list()
  
  # Process each unique combination
  for (i in 1:nrow(series_groups)) {
    # Extract current group
    current <- series_groups[i, ]
    
    # Filter data for current group
    group_data <- economic_long %>%
      dplyr::filter(
        Scenario == current$Scenario,
        Model == current$Model,
        Region == current$Region,
        Variable == current$Variable
      )
    
    # Create sequence of years for interpolation
    years_seq <- seq(start_year, end_year, by = dt)
    
    # Create a new dataframe with all years in the sequence
    new_years <- data.frame(
      Scenario = current$Scenario,
      Model = current$Model,
      Region = current$Region,
      Variable = current$Variable,
      Unit = current$Unit,
      Year = years_seq
    )
    
    # Merge with existing data to keep original values
    merged_data <- dplyr::full_join(
      new_years,
      group_data,
      by = c("Scenario", "Model", "Region", "Variable", "Unit", "Year")
    )
    
    # Sort by year
    merged_data <- merged_data %>% dplyr::arrange(Year)
    
    # Get available years and values (non-NA)
    available_years <- merged_data$Year[!is.na(merged_data$Value)]
    available_values <- merged_data$Value[!is.na(merged_data$Value)]
    
    # Only proceed with interpolation if we have at least 2 valid points
    if (length(available_values) >= 2) {
      # Perform linear interpolation
      interp_values <- approx(
        x = available_years,
        y = available_values,
        xout = merged_data$Year,
        method = "linear",
        rule = 2  # rule=2 means extrapolate beyond the data range
      )
      
      merged_data$Value <- interp_values$y
      
      # Store this result in the list
      result_list[[i]] <- merged_data
    } else if (length(available_values) == 1) {
      # If only one valid point, use that value for all years
      merged_data$Value <- available_values[1]
      
      # Store this result in the list
      result_list[[i]] <- merged_data
    }
    # If zero valid points, we don't add anything to the list
  }
  
  # Combine all results from the list using bind_rows
  # This is more reliable than repeatedly using rbind
  if (length(result_list) > 0) {
    interpolated_data <- dplyr::bind_rows(result_list)
    
    # Multiply all values by 1 trillion (1e12)
    #interpolated_data <- interpolated_data %>%
    #  dplyr::mutate(Value = Value * 1e12)
    
    return(interpolated_data)
  } else {
    warning("No valid interpolated data was produced. Check input data structure.")
    return(NULL)
  }
}
  
