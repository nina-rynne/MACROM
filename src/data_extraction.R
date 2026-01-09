# Load required library
library(dplyr)

# Define the vectors you want to extract
target_vectors <- c("baseline_annual_emissions", "baseline_annual_gwp", "cumulative_emissions", 
                    "temperature_anomaly", "qty_mitig", "qty_remov", "mitig_costs_annual", 
                    "remov_costs_annual", "temp_costs_annual", "total_costs_annual")

# Get the scenario results from your nested list
scenario_data <- scenario_results_all_SSPs_20250905_115847$scenario_results

# Get all SSP scenario names (SSP1-Baseline, SSP2-Baseline, etc.)
ssp_names <- names(scenario_data)

# Initialize an empty list to store dataframes for each SSP
all_data_list <- list()

# Loop through each SSP scenario
for (ssp in ssp_names) {
  
  # Get the data for this SSP scenario
  ssp_data <- scenario_data[[ssp]]
  
  # Initialize a dataframe for this SSP
  ssp_df <- data.frame()
  
  # Extract each target vector and add to dataframe
  for (vector_name in target_vectors) {
    if (vector_name %in% names(ssp_data)) {
      vector_data <- ssp_data[[vector_name]]
      
      # Handle different vector lengths by creating a temporary dataframe
      temp_df <- data.frame(
        SSP_Scenario = ssp,
        Variable = vector_name,
        Index = seq_along(vector_data),
        Value = vector_data
      )
      
      # Bind to the SSP dataframe
      ssp_df <- rbind(ssp_df, temp_df)
    } else {
      warning(paste("Vector", vector_name, "not found in", ssp))
    }
  }
  
  # Add this SSP's data to the master list
  all_data_list[[ssp]] <- ssp_df
}

# Combine all SSP dataframes into one master dataframe
final_df <- do.call(rbind, all_data_list)

# Remove row names for cleaner output
rownames(final_df) <- NULL

# Display the structure of the final dataframe
print("Structure of extracted data:")
str(final_df)
print("First few rows:")
head(final_df, 20)

# Export to CSV
output_filename <- paste0("extracted_SSP_data_", Sys.Date(), ".csv")
write.csv(final_df, output_filename, row.names = FALSE)

print(paste("Data exported to:", output_filename))

# Optional: Create a wide format version where each variable becomes a column
# This might be more useful for analysis
wide_df <- final_df %>%
  tidyr::pivot_wider(names_from = Variable, values_from = Value, 
                     id_cols = c(SSP_Scenario, Index))

# Export wide format version too
wide_output_filename <- paste0("extracted_SSP_data_wide_", Sys.Date(), ".csv")
write.csv(wide_df, wide_output_filename, row.names = FALSE)

print(paste("Wide format data exported to:", wide_output_filename))
print("Wide format structure:")
str(wide_df)

################################################################################

# Check what SSP names we actually have
print("Available SSP scenarios:")
print(ssp_names)

# Find max qty_mitig and qty_remov with their years for each SSP
for (ssp in ssp_names) {
  print(paste("Processing:", ssp))
  
  ssp_data <- scenario_data[[ssp]]
  
  # Check if the required vectors exist
  if ("qty_mitig" %in% names(ssp_data) && "qty_remov" %in% names(ssp_data) && "years" %in% names(ssp_data)) {
    
    # Find max qty_mitig and its year
    max_mitig_idx <- which.max(ssp_data$qty_mitig)
    max_mitig_value <- ssp_data$qty_mitig[max_mitig_idx]
    max_mitig_year <- ssp_data$years[max_mitig_idx]
    
    # Find max qty_remov and its year
    max_remov_idx <- which.max(ssp_data$qty_remov)
    max_remov_value <- ssp_data$qty_remov[max_remov_idx]
    max_remov_year <- ssp_data$years[max_remov_idx]
    
    print(sprintf("%s: Max mitig = %.2f in %d, Max remov = %.2f in %d", 
                  ssp, max_mitig_value, max_mitig_year, max_remov_value, max_remov_year))
  } else {
    print(paste("Missing required vectors in", ssp))
    print(paste("Available vectors:", paste(names(ssp_data), collapse = ", ")))
  }
}

################################################################################

# Find max cumulative_emissions with year for each SSP
for (ssp in ssp_names) {
  print(paste("Processing:", ssp))
  
  ssp_data <- scenario_data[[ssp]]
  
  # Check if the required vectors exist
  if ("cumulative_emissions" %in% names(ssp_data) && "years" %in% names(ssp_data)) {
    
    # Find max cumulative_emissions and its year
    min_cumul_idx <- which.min(ssp_data$cumulative_emissions)
    min_cumul_value <- ssp_data$cumulative_emissions[min_cumul_idx]
    min_cumul_year <- ssp_data$years[min_cumul_idx]
    
    print(sprintf("%s: Min cumulative emissions = %.2f in %d", 
                  ssp, min_cumul_value, min_cumul_year))
  } else {
    print(paste("Missing required vectors in", ssp))
    print(paste("Available vectors:", paste(names(ssp_data), collapse = ", ")))
  }
}


################################################################################

# Extract data from combined_results
combined_data <- multi_delay_results_20250905_165048$combined_results

# Define the vectors you want to extract
target_vectors <- c("mitigation_delay", "cdr_delay", "peak_temperature", "years_above_1p5", "scenario_short", "total_mitig_units", "total_cdr_units")

# Create a dataframe from the combined_results data
combined_df <- data.frame(
  Index = seq_along(combined_data$mitigation_delay),
  mitigation_delay = combined_data$mitigation_delay,
  cdr_delay = combined_data$cdr_delay,
  peak_temperature = combined_data$peak_temperature,
  years_above_1p5 = combined_data$years_above_1p5,
  scenario_short = combined_data$scenario_short,
  total_mitig_units = combined_data$total_mitig_units,
  total_cdr_units = combined_data$total_cdr_units
)

# Display the structure and first few rows
print("Structure of extracted combined results data:")
str(combined_df)
print("First few rows:")
head(combined_df, 10)

# Export to CSV
output_filename <- paste0("combined_results_data_", Sys.Date(), ".csv")
write.csv(combined_df, output_filename, row.names = FALSE)

print(paste("Combined results data exported to:", output_filename))
print(paste("Total rows exported:", nrow(combined_df)))

################################################################################


# Extract data from nested list structure and save to CSV

# Initialize empty list to store dataframes from each run
all_runs_data <- list()

# Get the successful_runs from the main object
successful_runs <- `sensitivity_results_SSP5-Baseline_20251029_063633`$successful_runs

# Loop through all available runs in successful_runs
for (run_name in names(successful_runs)) {
  
  # Get the current run data
  current_run <- successful_runs[[run_name]]
  
  # Check if all required variables exist in this run
  required_vars <- c("years", "cumulative_emissions", "temperature_anomaly", 
                     "qty_mitig", "qty_remov", "mitig_costs_annual", 
                     "remov_costs_annual", "temp_costs_annual", "total_costs_annual")
  
  # Skip this run if any required variables are missing
  if (!all(required_vars %in% names(current_run))) {
    cat("Skipping", run_name, "- missing required variables\n")
    next
  }
  
  # Extract the required vectors
  run_df <- data.frame(
    run_id = run_name,
    years = current_run$years,
    cumulative_emissions = current_run$cumulative_emissions,
    temperature_anomaly = current_run$temperature_anomaly,
    qty_mitig = current_run$qty_mitig,
    qty_remov = current_run$qty_remov,
    mitig_costs_annual = current_run$mitig_costs_annual,
    remov_costs_annual = current_run$remov_costs_annual,
    temp_costs_annual = current_run$temp_costs_annual,
    total_costs_annual = current_run$total_costs_annual
  )
  
  # Add to our list
  all_runs_data[[run_name]] <- run_df
  
  cat("Processed", run_name, "with", nrow(run_df), "rows\n")
}

# Combine all runs into one large dataframe
final_df <- do.call(rbind, all_runs_data)

# Create filename with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
filename <- paste0("sensitivity_results_output_values_SSP5_", timestamp, ".csv")

# Save to CSV
write.csv(final_df, filename, row.names = FALSE)

# Print summary
cat("\nExtraction complete!\n")
cat("Total runs processed:", length(all_runs_data), "\n")
cat("Total rows in final dataset:", nrow(final_df), "\n")
cat("File saved as:", filename, "\n")

# Display first few rows
cat("\nFirst few rows of the dataset:\n")
print(head(final_df))

################################################################################

# Extract run #### as its own nested list object

# Check if run exists in the successful_runs
if ("run_1126" %in% names(`sensitivity_results_SSP2-Baseline_20250905_140944`$successful_runs)) {
  
  # Extract run as its own object
  run_1126 <- `sensitivity_results_SSP2-Baseline_20250905_140944`$successful_runs$run_4414
  
  cat("Successfully extracted run_1126\n")
  cat("Structure of run_1126:\n")
  str(run_1126, max.level = 2)
  
  cat("\nAvailable elements in run_1126:\n")
  print(names(run_1126))
  
} else {
  cat("Error: run_1126 not found in successful_runs\n")
  cat("Available runs are:\n")
  print(names(`sensitivity_results_SSP2-Baseline_20250905_140944`$successful_runs))
}

################################################################################

# Create summary dataframe with sums for each run
library(ggplot2)

# Check if all_runs_data exists
if (!exists("all_runs_data")) {
  stop("Error: all_runs_data not found. Please run the extraction code first.")
}

# Create summary by summing qty_mitig and qty_remov for each run
run_summary <- data.frame(
  run_id = character(),
  sum_qty_mitig = numeric(),
  sum_qty_remov = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each run in all_runs_data
for (run_name in names(all_runs_data)) {
  
  # Get the current run dataframe
  current_run_df <- all_runs_data[[run_name]]
  
  # Calculate sums
  sum_mitig <- sum(current_run_df$qty_mitig, na.rm = TRUE)
  sum_remov <- sum(current_run_df$qty_remov, na.rm = TRUE)
  
  # Add to summary dataframe
  run_summary <- rbind(run_summary, data.frame(
    run_id = run_name,
    sum_qty_mitig = sum_mitig,
    sum_qty_remov = sum_remov
  ))
}

# Display summary
cat("Summary of runs:\n")
print(run_summary)

cat("\nSummary statistics:\n")
cat("Number of runs:", nrow(run_summary), "\n")
cat("Mitigation - Min:", round(min(run_summary$sum_qty_mitig), 2), 
    ", Max:", round(max(run_summary$sum_qty_mitig), 2), "\n")
cat("Removal - Min:", round(min(run_summary$sum_qty_remov), 2), 
    ", Max:", round(max(run_summary$sum_qty_remov), 2), "\n")

# Calculate correlation
correlation <- cor(run_summary$sum_qty_mitig, run_summary$sum_qty_remov)
cat("Correlation between mitigation and removal:", round(correlation, 3), "\n")

# Create scatter plot with correlation
scatter_plot <- ggplot(run_summary, aes(x = sum_qty_mitig, y = sum_qty_remov)) +
  geom_point(alpha = 0.7, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
  labs(
    title = "Mitigation vs Removal by Run",
    subtitle = paste("Total runs:", nrow(run_summary), "| Correlation:", round(correlation, 3)),
    x = "Sum of Quantity Mitigation",
    y = "Sum of Quantity Removal"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Display the plot
print(scatter_plot)

# Save the plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
plot_filename <- paste0("mitigation_vs_removal_scatter_", timestamp, ".png")
ggsave(plot_filename, scatter_plot, width = 10, height = 7, dpi = 300)

cat("\nPlot saved as:", plot_filename, "\n")

# Save the summary dataframe as well
summary_filename <- paste0("run_summary_", timestamp, ".csv")
write.csv(run_summary, summary_filename, row.names = FALSE)
cat("Summary dataframe saved as:", summary_filename, "\n")

################################################################################

# Extract vectors from the summary sub-list and save to CSV

# Access the summary data from the main object
summary_data <- `sensitivity_results_SSP5-Baseline_20251029_063633`$summary

# Check if summary exists
if (is.null(summary_data)) {
  stop("Error: summary sub-list not found in the main object")
}

# Check if required vectors exist
required_vectors <- c("run_id", "peak_temperature", "years_above_1p5")
missing_vectors <- setdiff(required_vectors, names(summary_data))

if (length(missing_vectors) > 0) {
  stop(paste("Error: Missing vectors:", paste(missing_vectors, collapse = ", ")))
}

# Extract the required vectors
run_id <- summary_data$run_id
peak_temperature <- summary_data$peak_temperature
years_above_1p5 <- summary_data$years_above_1p5

# Check that all vectors are the same length
vector_lengths <- c(length(run_id), length(peak_temperature), length(years_above_1p5))
if (length(unique(vector_lengths)) > 1) {
  cat("Warning: Vectors have different lengths:", paste(vector_lengths, collapse = ", "), "\n")
}

# Create dataframe
summary_df <- data.frame(
  run_id = run_id,
  peak_temperature = peak_temperature,
  years_above_1p5 = years_above_1p5
)

# Create filename with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
filename <- paste0("summary_extracted_data_SSP5_", timestamp, ".csv")

# Save to CSV
write.csv(summary_df, filename, row.names = FALSE)

# Print summary information
cat("Data extraction complete!\n")
cat("Number of runs:", nrow(summary_df), "\n")
cat("Variables extracted: run_id, peak_temperature, years_above_1p5\n")
cat("File saved as:", filename, "\n")

# Display first few rows
cat("\nFirst few rows of the extracted data:\n")
print(head(summary_df))

# Display some basic statistics
cat("\nSummary statistics:\n")
cat("Peak temperature - Min:", round(min(summary_df$peak_temperature, na.rm = TRUE), 3), 
    ", Max:", round(max(summary_df$peak_temperature, na.rm = TRUE), 3), "\n")
cat("Years above 1.5°C - Min:", min(summary_df$years_above_1p5, na.rm = TRUE), 
    ", Max:", max(summary_df$years_above_1p5, na.rm = TRUE), "\n")


################################################################################

# Analysis of sensitivity results nested list
# Filter runs by final_temperature < 1.6 and extract key metrics

# Note: The top-level list name contains a hyphen, so we use backticks to access it
# Original name: sensitivity_results_SSP1-Baseline_20251029_194435

# Step 1: Access the successful_runs list
# Using backticks to handle the hyphen in the variable name
successful_runs <- `sensitivity_results_SSP5-Baseline_20251029_063633`$successful_runs

# Step 2: Filter runs where final_temperature < 1.6
# Extract run names and their final temperatures
run_names <- names(successful_runs)

# Create a logical vector for filtering
runs_below_threshold <- sapply(run_names, function(run_name) {
  final_temp <- successful_runs[[run_name]]$final_temperature
  # Check if final_temperature exists and is below 1.6
  if (!is.null(final_temp) && length(final_temp) > 0) {
    return(final_temp < 1.6)
  } else {
    return(FALSE)
  }
})

# Filter the runs
filtered_run_names <- run_names[runs_below_threshold]
filtered_runs <- successful_runs[filtered_run_names]

# Step 3: Calculate number of runs remaining
n_runs_remaining <- length(filtered_runs)
cat("Number of runs remaining (final_temperature < 1.6):", n_runs_remaining, "\n\n")

# Step 4: Extract data for analysis
# Create a dataframe with all temperature_anomaly values, qty_mitig, and qty_remov

# Initialize lists to store results
all_temp_anomalies <- list()
all_qty_mitig <- list()
all_qty_remov <- list()
runs_with_temp_above_1_52 <- 0

for (run_name in filtered_run_names) {
  run_data <- filtered_runs[[run_name]]
  
  # Extract temperature anomaly
  temp_anomaly <- run_data$temperature_anomaly
  years <- run_data$years
  
  # Check if any year has temperature_anomaly > 1.52
  if (any(temp_anomaly > 1.52, na.rm = TRUE)) {
    runs_with_temp_above_1_52 <- runs_with_temp_above_1_52 + 1
  }
  
  # Store data with run_name and year information
  for (i in seq_along(temp_anomaly)) {
    all_temp_anomalies[[length(all_temp_anomalies) + 1]] <- data.frame(
      run_name = run_name,
      year = years[i],
      temperature_anomaly = temp_anomaly[i]
    )
  }
  
  # Extract qty_mitig
  qty_mitig <- run_data$qty_mitig
  for (i in seq_along(qty_mitig)) {
    all_qty_mitig[[length(all_qty_mitig) + 1]] <- data.frame(
      run_name = run_name,
      year = years[i],
      qty_mitig = qty_mitig[i]
    )
  }
  
  # Extract qty_remov
  qty_remov <- run_data$qty_remov
  for (i in seq_along(qty_remov)) {
    all_qty_remov[[length(all_qty_remov) + 1]] <- data.frame(
      run_name = run_name,
      year = years[i],
      qty_remov = qty_remov[i]
    )
  }
}

# Combine into dataframes
df_temp <- do.call(rbind, all_temp_anomalies)
df_mitig <- do.call(rbind, all_qty_mitig)
df_remov <- do.call(rbind, all_qty_remov)

# Step 5: Find maximum values and their locations

# Highest temperature_anomaly
max_temp_idx <- which.max(df_temp$temperature_anomaly)
max_temp_info <- df_temp[max_temp_idx, ]

cat("Highest temperature_anomaly:\n")
cat("  Value:", max_temp_info$temperature_anomaly, "\n")
cat("  Run:", max_temp_info$run_name, "\n")
cat("  Year:", max_temp_info$year, "\n\n")

# Highest qty_mitig
max_mitig_idx <- which.max(df_mitig$qty_mitig)
max_mitig_info <- df_mitig[max_mitig_idx, ]

cat("Highest qty_mitig:\n")
cat("  Value:", max_mitig_info$qty_mitig, "\n")
cat("  Run:", max_mitig_info$run_name, "\n")
cat("  Year:", max_mitig_info$year, "\n\n")

# Highest qty_remov
max_remov_idx <- which.max(df_remov$qty_remov)
max_remov_info <- df_remov[max_remov_idx, ]

cat("Highest qty_remov:\n")
cat("  Value:", max_remov_info$qty_remov, "\n")
cat("  Run:", max_remov_info$run_name, "\n")
cat("  Year:", max_remov_info$year, "\n\n")

# Number of runs with at least one year > 1.52
cat("Number of runs with at least one year where temperature_anomaly > 1.52:", 
    runs_with_temp_above_1_52, "\n\n")

# Step 6: Create summary results object
summary_results <- list(
  n_runs_remaining = n_runs_remaining,
  max_temperature_anomaly = list(
    value = max_temp_info$temperature_anomaly,
    run = max_temp_info$run_name,
    year = max_temp_info$year
  ),
  n_runs_temp_above_1_52 = runs_with_temp_above_1_52,
  max_qty_mitig = list(
    value = max_mitig_info$qty_mitig,
    run = max_mitig_info$run_name,
    year = max_mitig_info$year
  ),
  max_qty_remov = list(
    value = max_remov_info$qty_remov,
    run = max_remov_info$run_name,
    year = max_remov_info$year
  )
)

# Print final summary
cat("=" , rep("=", 50), "\n", sep = "")
cat("SUMMARY RESULTS\n")
cat("=" , rep("=", 50), "\n", sep = "")
print(summary_results)

# Optional: Save the filtered dataframes for further analysis
# write.csv(df_temp, "temperature_anomaly_data.csv", row.names = FALSE)
# write.csv(df_mitig, "qty_mitig_data.csv", row.names = FALSE)
# write.csv(df_remov, "qty_remov_data.csv", row.names = FALSE)