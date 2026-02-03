# MACROM User Guide

**Version:** 1.0.0  
**Last Updated:** January 2026  
**For:** MACROM v1.0.0

This guide provides comprehensive instructions for using MACROM (Mitigation and Carbon Removal Optimisation Model).

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation and Setup](#installation-and-setup)
3. [Prerequisites](#prerequisites)
4. [Understanding the Workflow](#understanding-the-workflow)
5. [Quick Start Guide](#quick-start-guide)
6. [Detailed Workflow Instructions](#detailed-workflow-instructions)
7. [Analysis Types](#analysis-types)
8. [Customisation Options](#customisation-options)
9. [Understanding Outputs](#understanding-outputs)
10. [Troubleshooting](#troubleshooting)
12. [Frequently Asked Questions](#frequently-asked-questions)

---

## Introduction

MACROM is an optimal control model that determines the most cost-effective strategies for climate mitigation and carbon dioxide removal (CDR) under different emission scenarios. This guide will help you:

- Set up and run the model
- Choose the appropriate analysis for your research questions
- Interpret and customise outputs
- Troubleshoot common issues

### Who Should Use This Guide

- Researchers analysing climate policy scenarios
- Students learning about optimal control in climate economics
- Policy analysts exploring mitigation and CDR strategies
- Anyone interested in reproducing or extending MACROM analyses

### What You'll Need

- Basic familiarity with R and RStudio
- Understanding of climate change basics (emissions, temperature targets, mitigation)
- Approximately 30 minutes to 2 hours (depending on analysis complexity)

---

## Installation and Setup

### Step 1: Install Required Software

1. **Install R** (version 4.0.0 or higher)
   - Download from: https://cran.r-project.org/
   - Follow installation instructions for your operating system

2. **Install RStudio** (recommended but optional)
   - Download from: https://posit.co/download/rstudio-desktop/
   - RStudio provides a user-friendly interface for R

### Step 2: Download MACROM

**Option A: Using Git (recommended)**
```bash
git clone https://github.com/nina-rynne/MACROM.git
cd MACROM
```

**Option B: Download ZIP**
1. Visit https://github.com/nina-rynne/MACROM
2. Click "Code" → "Download ZIP"
3. Extract the ZIP file to your desired location

### Step 3: Install Required R Packages

Open R or RStudio and run:
```r
# Data manipulation
install.packages(c("dplyr", "tidyr", "purrr", "tidyverse", "readr"))

# File handling
install.packages(c("here", "yaml", "Cairo"))

# Parallel processing
install.packages(c("parallel", "foreach", "doParallel"))

# Analysis
install.packages(c("lhs", "iterators"))

# Visualisation
install.packages(c("ggplot2", "cowplot", "patchwork", "viridisLite", 
                   "viridis", "colorspace", "scales"))
```

**Or install all at once:**
```r
packages <- c("dplyr", "tidyr", "purrr", "tidyverse", "readr", "here", 
              "yaml", "Cairo", "parallel", "foreach", "doParallel", 
              "lhs", "iterators", "ggplot2", "cowplot", "patchwork", 
              "viridisLite", "viridis", "colorspace", "scales")
install.packages(packages)
```

### Step 4: Verify Data Files

Ensure these files are present in the `data/` directory:
- `emissions.csv`: SSP emission scenarios
- `gwp.csv`: Gross World Product projections for SSP scenarios

**Note:** If these files are missing, contact the repository maintainer or check the `data-raw/` directory for instructions on generating them.

### Step 5: Verify Directory Structure

Your MACROM directory should look like this:
```
MACROM/
├── MACROM_workflow.Rmd     # Main user interface
├── README.md               # This file
├── USER_GUIDE.md           # Detailed documentation
├── LICENSE                 # CC_BY_4.0 licence
├── parameter_details.yml
├── parameter_details_fixed.yml
├── data/                   # Processed data files
│   ├── emissions.csv
│   └── gwp.csv
├── data-raw/               # Raw, unprocessed data
├── figs/                   # Generated figures and visualisations
├── output/                 # Analysis outputs (RDS, CSV)
└── src/                    # Source code (R functions)
    ├── data_extraction.R
    ├── data_preparation.R
    ├── delayed_deployment.R
    ├── delayed_deployment_visualisation.R
    ├── latin_hypercube_sampling.R
    ├── model_parameters.R
    ├── optimal_control_core.R
    ├── parameter_importance.R
    ├── parameter_importance_visualisation.R
    ├── scenario_comparison.R
    └── scenario_comparison_visualisation.R
```

---

## Prerequisites

### Technical Requirements

- **Operating System**: Windows, macOS, or Linux
- **R Version**: 4.0.0 or higher
- **RAM**: Minimum 8GB recommended (16GB for large parameter sets)
- **Processor**: Multi-core processor recommended for parallel processing
- **Disk Space**: ~1GB for software and outputs

### Knowledge Requirements

**Essential:**
- Basic R programming (running scripts, understanding data frames)
- Familiarity with RStudio or R command line
- Understanding of file paths and directory structures

**Helpful but not required:**
- Climate change concepts (emissions, temperature anomalies, Paris Agreement)
- Basic economics (cost-benefit analysis, discounting)
- Optimal control theory (helpful for understanding methodology)

### Time Commitment

**Estimated runtimes:**
- Setup and installation: 15-30 minutes (first time only)
- Single scenario with fixed parameters: 5-15 minutes
- All 5 SSP scenarios with fixed parameters: 10-20 minutes
- Parameter sensitivity (20,000 samples, 1 scenario): 1-3 hours
- Delayed deployment analysis (5 scenarios): 30-60 minutes

**Note:** Runtimes vary significantly based on hardware and whether parallel processing is enabled.

---

## Understanding the Workflow

### Overview

MACROM uses a modular workflow structure. The main user interface is `MACROM_workflow.Rmd`, which calls functions from the `src/` directory. You interact only with the workflow file—all complexity is handled by the underlying functions.

### Workflow Structure

The workflow consists of sequential chunks that must be run in order:

1. **Setup** (always required)
   - Libraries
   - Data preparation

2. **Parameter Selection** (choose ONE)
   - Fixed parameters (single parameter set)
   - Latin Hypercube Sampling (multiple parameter sets)

3. **Analysis** (choose one or more)
   - Scenario comparison
   - Parameter importance
   - Delayed deployment

4. **Visualisation** (follows each analysis)
   - Automated figure generation
   - Customisable plotting options

### Decision Tree: Which Analysis Should I Run?
```
START: What is your research question?
│
├─ "How do optimal strategies differ across SSP scenarios?"
│   → Use: Fixed Parameters + Scenario Comparison
│
├─ "Which parameters most influence outcomes?"
│   → Use: Latin Hypercube Sampling + Parameter Importance
│
├─ "What are the costs of delaying action?"
│   → Use: Fixed Parameters + Delayed Deployment
│
└─ "All of the above (comprehensive analysis)"
    → Run all analyses sequentially
```

### Mutually Exclusive Choices

**IMPORTANT:** Some chunks cannot be run together:

❌ **DO NOT RUN BOTH:**
- `latin_hypercube_sampling_call` AND `fixed_model_parameters_call`
- Choose ONE based on your analysis needs

✓ **CAN RUN TOGETHER:**
- Multiple analysis types (scenario comparison, parameter importance, delayed deployment)
- Multiple visualisation chunks

---

## Quick Start Guide

This guide gets you running MACROM with default settings in under 15 minutes.

### Step-by-Step Instructions

#### 1. Open the Workflow

1. Open RStudio
2. Navigate to File → Open File
3. Select `MACROM_workflow.Rmd` from the MACROM directory
4. Click "Open"

#### 2. Run Setup Chunks (REQUIRED)

Run these chunks in order by clicking the green "play" button in the top-right of each chunk:

**Chunk: `libraries`**
- Loads all required R packages
- May take 1-2 minutes on first run
- If errors occur, install missing packages (see Installation section)

**Chunk: `data_preparation_call`**
- Loads and processes emissions and economic data
- Interpolates to annual time steps
- Takes ~30 seconds
- Creates `emissions_df` and `economic_df` objects

#### 3. Choose Parameter Setup (Choose ONE)

**For most users (recommended):**

**Chunk: `fixed_model_parameters_call`**
- Loads single set of default parameter values
- Creates `parameter_df` with 1 row
- Best for: Scenario comparison, Delayed deployment
- Runtime: <1 second

**For sensitivity analysis:**

**Chunk: `latin_hypercube_sampling_call`**
- Generates multiple parameter combinations
- Default: 500 samples
- Best for: Parameter importance analysis
- Runtime: ~5 seconds
- Creates `parameter_df` with 500 rows

#### 4. Run Your First Analysis

**Recommended first analysis: Scenario Comparison**

```r
# In chunk: scenario_comparison_analysis
# Run with default settings (all 5 SSP scenarios)
```

**Expected runtime:** 10-15 minutes with parallel processing
**Output:** `scenario_comparison_results` object

#### 5. View Results

**Chunk: `scenario_comparison_visualisation`**
- Creates publication-quality figures
- Saves to `figs/` directory automatically
- Runtime: ~30 seconds

**You should see:**
- Time series plots of optimal mitigation and CDR, cumulative emissions, temperature, implementation costs and temperature related damages. 

---

## Detailed Workflow Instructions

### 6.1 Libraries

**Chunk:** `libraries`

**Purpose:** Load all required R packages.

**What it does:**
Loads all R packages required for MACROM to function, including packages for data manipulation, file handling, parallel processing, statistical analysis, and visualisation. The libraries are organised into functional groups: data transformation (dplyr, tidyr, purrr), file I/O (readr, here, yaml), parallel processing (parallel, foreach, doParallel), analysis (lhs, iterators), and visualisation (ggplot2, cowplot, patchwork, viridis packages).

**Troubleshooting:**
- If you get "there is no package called X", install it: `install.packages("X")`
- If multiple errors, run the batch install command from Section 2

**Expected output:** Messages showing packages loaded (may include warnings—these are usually harmless)

---

### 6.2 Data Preparation

**Chunk:** `data_preparation_call`

**Purpose:** Load and process emission and economic data.

**What it does:**
1. Loads raw SSP scenario data from `data/emissions.csv` and `data/gwp.csv`
2. Interpolates to annual time steps (default: 2020-2100)
3. Prepares data in format needed by optimal control algorithm

**Customisation:**
```r
# Change time horizon
start_year = 2020   # Model start year
end_year = 2100     # Model end year
dt = 1              # Time step (years)
```

**Output:**
- `emissions_df`: Annual emissions for all SSP scenarios (GtCO₂/year)
- `economic_df`: Annual GWP for all SSP scenarios (trillion $)

**Common modifications:**
- **Shorter horizon:** Set `end_year = 2050`
- **Finer resolution:** Set `dt = 2` (biennial)
- **Different start:** Set `start_year = 2025`

---

### 6.3 Parameter Setup: Fixed Parameters

**Chunk:** `fixed_model_parameters_call`

**Purpose:** Load single set of model parameters from YAML file.

**When to use:**
- Scenario comparison across SSP scenarios
- Delayed deployment analysis
- Any analysis where you want consistent parameters

**What it does:**
1. Reads `parameter_details_fixed.yml`
2. Creates single-row data frame with all model parameters
3. Parameters include:
   - Climate physics (TCRE, initial temperature)
   - Economic parameters (discount rate, damage coefficient, costs)
   - Model configuration (exponents, bounds, delays)

**Output:** `parameter_df` with 1 row

**To modify parameters:**
Edit `parameter_details_fixed.yml`:
```yaml
# Example modifications
disc_rate: 0.05        # Change discount rate from 0.03 to 0.05
tcre: 0.50             # Change climate sensitivity from 0.45 to 0.50
```

Then re-run this chunk.

---

### 6.4 Parameter Setup: Latin Hypercube Sampling

**Chunk:** `latin_hypercube_sampling_call`

**Purpose:** Generate multiple parameter sets for sensitivity analysis.

**When to use:**
- Parameter importance analysis
- Uncertainty quantification
- Understanding model behaviour across parameter space

**What it does:**
1. Reads parameter ranges from `parameter_details.yml`
2. Generates `n_samples` parameter combinations using Latin Hypercube Sampling
3. LHS ensures good coverage of parameter space

**Customisation:**
```r
n_samples = 500      # Number of parameter combinations to test
seed = 12345         # Random seed for reproducibility
```

**Why Latin Hypercube Sampling?**
- More efficient than random sampling
- Ensures good coverage with fewer samples
- Standard method for sensitivity analysis

**Computational impact:**
- 100 samples: ~5-10 minutes for parameter importance
- 5,000 samples: ~1 hours for parameter importance
- 20,000 samples: ~2-3 hours for parameter importance

**Output:** `parameter_df` with multiple rows (e.g., 500 rows)

---

### 6.5 Scenario Comparison Analysis

**Chunk:** `scenario_comparison_analysis`

**Purpose:** Compare optimal strategies across different SSP emission scenarios.

**Requirements:**
- Must use fixed parameters (`parameter_df` with 1 row)
- Run `fixed_model_parameters_call` first

**When to use:**
- Understanding how baseline emissions affect optimal strategies
- Comparing costs and temperature outcomes across futures
- Creating publication figures showing scenario diversity

**What it does:**
1. Runs optimal control algorithm for each specified scenario
2. Finds optimal mitigation and CDR trajectories
3. Calculates resulting temperatures and costs
4. Compares metrics across scenarios

**Customisation options:**
```r
# Select scenarios to compare
scenarios_to_compare <- c("SSP1-Baseline", "SSP2-Baseline", 
                          "SSP3-Baseline", "SSP4-Baseline", 
                          "SSP5-Baseline")

# Processing options
use_parallel = TRUE     # Highly recommended for multiple scenarios
verbose = TRUE          # Show progress messages
```

**Output:**

`scenario_comparison_results` contains:
- `scenario_results`: Full time series for each scenario
- `comparison_summary`: Summary table comparing key metrics
- `year_first_1p5C`: Year when temperature first reaches 1.5°C
- `year_peak_temp`: Year of peak temperature
- Runtime and configuration metadata

**Expected runtime:**
- 5 scenarios (parallel): 10-15 minutes
- 5 scenarios (sequential): 30-45 minutes

---

### 6.6 Parameter Importance Analysis

**Chunk:** `parameter_importance_analysis`

**Purpose:** Run optimal control for multiple parameter combinations to enable sensitivity analysis.

**Requirements:**
- Must use Latin Hypercube Sampling (`parameter_df` with multiple rows)
- Run `latin_hypercube_sampling_call` first
- **Note:** Only parameter importance analysis can use LHS

**When to use:**
- Exploring parameter sensitivity
- Understanding model behaviour across parameter space
- Quantifying uncertainty in outcomes

**What it does:**
1. Runs optimal control for all parameter combinations
2. Collects results for each parameter set
3. Creates summary statistics for all successful runs
4. Separates successful and failed runs

**Customisation options:**
```r
# Select scenario for analysis
scenario_for_importance <- "SSP2-Baseline"  # Usually middle-of-road scenario

# Processing options
use_parallel = TRUE     # Highly recommended for large parameter sets
verbose = FALSE         # Set TRUE to see progress
```

**Why use SSP2-Baseline?**
- Represents "middle-of-road" socio-economic development
- Balances between optimistic (SSP1) and pessimistic (SSP3, SSP5) scenarios
- Most commonly used reference scenario in climate literature

**Output:**

`parameter_importance_results` contains:
- `successful_runs`: List of all successful optimal control solutions
- `failed_runs`: List of failed runs with error information
- `summary`: Data frame with key metrics from all successful runs (peak temperature, costs, emissions, convergence status, etc.)
- `run_info`: Metadata about the analysis run

**Expected runtime:**
- 100 samples: ~5-10 minutes for parameter importance
- 5,000 samples: ~1 hours for parameter importance
- 20,000 samples: ~2-3 hours for parameter importance

**Interpreting results:**
The summary data frame contains outcome metrics for each parameter set, which can then be used to:
- Calculate correlations between parameters and outcomes
- Identify parameter combinations that produce extreme results
- Assess model sensitivity to parameter uncertainty

---

### 6.7 Parameter Importance Filtering

**Chunk:** `parameter_importance_filtering`

**Purpose:** Filter parameter importance results based on temperature criteria.

**Requirements:**
- Must have run `parameter_importance_analysis` first

**What it does:**
Filters the results from parameter importance analysis to retain only parameter sets that meet specified temperature criteria (e.g., maximum acceptable final temperature), which helps focus sensitivity analysis on realistic parameter combinations.

**Customisation options:**
```r
# Set maximum acceptable final temperature
max_acceptable_temperature <- 1.6  # Default: 1.6°C
```

**Output:**
`parameter_importance_results_cleaned` - filtered version of results containing only acceptable parameter combinations.

---

### 6.8 Parameter Importance Visualisation

**Chunk:** `parameter_importance_visualisation`

**Purpose:** Create visualisations showing variation in optimal control trajectories across parameter sets.

**What it does:**
Creates a 3×2 grid of "spaghetti plots" showing how optimal control solutions vary across all parameter combinations. Each plot shows all individual runs as transparent lines with the median trajectory highlighted. The six panels show:
1. **Cumulative emissions** over time
2. **Temperature anomaly** trajectories
3. **Mitigation** deployment over time
4. **CDR** deployment over time
5. **Abatement costs** (as proportion of GWP) over time
6. **Damage costs** (as proportion of GWP) over time

**Customisation options:**
```r
# Styling options
save_plot = TRUE        # Save to file
print_insights = FALSE  # Print detailed statistical insights
verbose = TRUE          # Show progress messages
```

**Output:**
A combined dashboard saved to the `figs/` directory showing the range of possible trajectories and highlighting median behaviour across all parameter sets.

**Understanding the plots:**
- Individual runs shown as transparent lines reveal the spread of possible outcomes
- Median lines (darker, more opaque) show the central tendency
- Wide spread indicates high sensitivity to parameter uncertainty
- Narrow spread indicates robust behaviour across parameter space

---

### 6.9 Parameter Importance Export for Sobol

**Chunk:** `export_parameter_importance_for_sobol`

**Purpose:** Export parameter importance results for external Sobol indices analysis.

**Requirements:**
- Must have run `parameter_importance_filtering` first

**What it does:**
Exports the filtered parameter importance results to CSV format suitable for Sobol sensitivity analysis, which provides variance-based sensitivity indices.

**Output:**
CSV files saved to the `output/` directory containing parameter values and outcome variables formatted for external Sobol analysis tools.

---

### 6.10 Delayed Deployment Analysis

**Chunk:** `delayed_deployment_analysis`

**Purpose:** Analyse effects of delaying mitigation and/or CDR deployment.

**Requirements:**
- Must use single parameter set (`parameter_df` with 1 row)
- Recommended: Run `fixed_model_parameters_call` first

**When to use:**
- Assessing costs of policy delays
- Understanding urgency of climate action
- Exploring trade-offs between mitigation and CDR timing
- Identifying critical deployment windows

**What it does:**
1. Runs optimal control across grid of delay combinations
2. Tests delays from 0 to `max_delay_years` in `delay_step_size` increments
3. Calculates outcomes for all delay combinations
4. Identifies feasible vs infeasible delay scenarios

**Customisation options:**
```r
# Delay range and granularity
max_delay_years <- 70        # Maximum delay to test (years)
delay_step_size <- 10        # Step size (years between tested delays)

# Examples of step size effects:
# step_size = 1:  Tests 0, 1, 2, 3, ... 70 (71 × 71 = 5,041 combinations)
# step_size = 5:  Tests 0, 5, 10, 15, ... 70 (15 × 15 = 225 combinations)
# step_size = 10: Tests 0, 10, 20, 30, ... 70 (8 × 8 = 64 combinations)

# Scenarios to analyse
scenarios_to_analyse <- "all"  # Test all 5 SSP scenarios

# Or select specific scenarios
scenarios_to_analyse <- c("SSP2-Baseline", "SSP3-Baseline")

# Processing options
use_parallel = TRUE     # Highly recommended for this analysis
```

**Computational considerations:**

Total runs = (max_delay / step_size + 1)² × number of scenarios

Examples:
- 70-year max, 10-year steps, 1 scenario: 64 runs (~5 minutes)
- 70-year max, 10-year steps, 5 scenarios: 320 runs (~25-30 minutes)
- 70-year max, 5-year steps, 5 scenarios: 1,125 runs (~1 hour)
- 70-year max, 1-year steps, 5 scenarios: 25,205 runs (~3 hours)

**Recommendation:** Start with `step_size = 10` and `scenarios_to_analyse = "SSP2-Baseline"` for quick testing.

**Output:**

`delayed_deployment_results` contains:
- `results_by_scenario`: Results for each scenario separately
- `combined_results`: All results combined across scenarios
- `summary_stats`: Aggregated metrics and feasibility rates
- Grid of outcomes (temperature, costs) across delay space

**Interpreting results:**

- **Feasible region**: Combinations where temperature target can be met
- **Infeasible region**: Delays too long to achieve target
- **Cost gradients**: How costs increase with delays
- **Critical thresholds**: Maximum delays before infeasibility

---

### 6.11 Delayed Deployment Visualisation

**Chunk:** `delayed_deployment_visualisation`

**Purpose:** Create heatmap visualisations showing how delays affect outcomes.

**What it does:**
- Generates 2D heatmaps: mitigation delay (x-axis) vs CDR delay (y-axis)
- Colour intensity shows outcome magnitude (temperature, cost, etc.)
- Optional contour lines show gradients
- Optional arrows show direction of steepest change
- Marks infeasible combinations

**Customisation options:**
```r
# Select variables to plot
variables = c("peak_temperature")

# Available variables:
#   - peak_temperature: Peak temperature anomaly (°C)
#   - years_above_1p5: Years above 1.5°C threshold
#   - total_cost: Sum of mitigation, CDR, and damage costs
#   - mitig_cost: Mitigation costs only
#   - remov_cost: CDR/removal costs only
#   - temp_cost: Temperature-related damage costs

# Visual features
add_contours = TRUE          # Add white contour lines
contour_alpha = 0.6          # Contour line transparency (0-1)
add_arrows = FALSE           # Add gradient vector field arrows
arrow_scale = 3.0            # Arrow length scaling
arrow_skip = 1               # Arrow thinning (1 = all, 2 = every other)
show_infeasible = FALSE      # Mark infeasible combinations with X

# Colour scales
color_palettes = NULL        # NULL = use defaults (plasma for temp, viridis for cost)
shared_scale = FALSE         # FALSE = independent scales per variable

# Output options
save_plot = TRUE
width = 297                  # A4 landscape width
height = 210                 # Increase for multiple variables
```

**Understanding the heatmaps:**

- **Lower left corner (0,0)**: No delays—immediate action
- **Upper right corner (max,max)**: Maximum delays for both technologies (e.g., 70,70 if max_delay_years = 70)
- **Diagonal**: Equal delays for mitigation and CDR
- **Colour intensity**: Magnitude of outcome variable
- **Contour lines**: Connect points with same values
- **Arrows**: Point in direction of steepest increase

**When to use arrows:**
- Showing direction of greatest cost increase
- Identifying optimal timing strategies
- Visualising trade-offs between mitigation and CDR timing
- **Warning:** Can clutter plots—use sparingly or with `arrow_skip > 1`

---

### 6.12 Pre-built Visualisation Functions

**Chunk:** `delayed_deployment_prebuilt_visualisation`

**Purpose:** Convenience wrappers for common visualisation types.

**What it does:**
Provides three pre-configured dashboard functions for standard analyses.

#### Temperature Dashboard
```r
temperature_dashboard <- create_temperature_delay_dashboard(
  deployment_results = delayed_deployment_results,
  save_plot = FALSE
)
```

**Shows:**
- Peak temperature heatmap
- Years above 1.5°C heatmap
- Pre-configured with appropriate colour scales and labels

**Use when:** Focusing on temperature outcomes and overshoot

#### Cost Dashboard with Arrows
```r
cost_dashboard <- create_cost_delay_dashboard(
  deployment_results = delayed_deployment_results,
  cost_variable = "total_cost",
  arrow_scale = 3.0,
  min_magnitude = 0,
  save_plot = FALSE
)
```

**Shows:**
- Cost heatmap with gradient arrows
- Arrows show direction of steepest cost increase

**Use when:** Understanding cost sensitivity to deployment timing

**Options for `cost_variable`:**
- `"total_cost"`: Sum of mitigation, CDR, and temperature damage costs
- `"temp_cost"`: Temperature damage costs only
- `"mitig_cost"`: Mitigation costs only
- `"remov_cost"`: CDR costs only

#### Comprehensive Cost Comparison
```r
comprehensive_cost_dashboard <- create_comprehensive_cost_delay_dashboard(
  deployment_results = delayed_deployment_results,
  shared_scale = TRUE,
  save_plot = FALSE,
  height = 260
)
```

**Shows:**
- Three-panel comparison: total cost, abatement cost (mitigation + CDR combined), temperature damage cost
- Shared colour scale for direct magnitude comparison

**Use when:** Understanding cost composition and trade-offs

**Note:** 
- Set `shared_scale = FALSE` for independent scales if cost components vary widely in magnitude
- **IMPORTANT**: The code internally uses variables named `"total_cost"`, `"abatement_cost"` (which combines mitigation + CDR costs), and `"temp_cost"` (temperature damage costs). The `abatement_cost` variable represents the sum of mitigation and CDR costs.

---

## Analysis Types

This section summarises the three main analysis types and when to use each.

### 7.1 Scenario Comparison

**Research Question:** "How do optimal mitigation and CDR strategies differ across socio-economic futures?"

**What it does:**
- Compares optimal control solutions across SSP1-5 baseline scenarios
- Shows how baseline emissions affect optimal strategies

**Use for:**
- Understanding range of possible futures
- Identifying common patterns across scenarios
- Comparing costs and temperature outcomes
- Publication figures showing scenario diversity

**Parameters:** Single parameter set (fixed parameters)  
**Outputs:** Time series, summary statistics, multi-panel dashboard  
**Runtime:** 10-20 minutes (5 scenarios, parallel)

---

### 7.2 Parameter Importance

**Research Question:** "Which parameters most influence model outcomes and where should we focus calibration efforts?"

**What it does:**
- Tests model sensitivity across parameter space
- Identifies influential vs negligible parameters
- Quantifies uncertainty propagation

**Use for:**
- Understanding model behaviour
- Prioritising parameter calibration
- Communicating key uncertainties
- Identifying potential model improvements

**Parameters:** Multiple parameter sets (Latin Hypercube Sampling)  
**Outputs:** Sensitivity rankings, correlation matrices, tornado diagrams  
**Runtime:** 2-3 hours (20,000 samples, 1 scenario, parallel)

---

### 7.3 Delayed Deployment

**Research Question:** "What are the consequences of delaying mitigation or CDR implementation?"

**What it does:**
- Tests all combinations of mitigation and CDR delays
- Maps feasible vs infeasible delay regions
- Quantifies cost penalties of delay

**Use for:**
- Assessing urgency of climate action
- Understanding technology timing trade-offs
- Identifying critical deployment windows
- Policy analysis of implementation delays

**Parameters:** Single parameter set (fixed parameters)  
**Outputs:** Heatmaps, feasibility maps, cost gradients  
**Runtime:** 5-10 minutes (5 scenarios, 10-year steps, parallel)

---

## Customisation Options

### 8.1 Time Range

**Where:** `data_preparation_call` chunk
```r
emissions_df <- interpolate_ssp_emissions(
  emissions_df = emissions_imported,
  dt = 1,              # Time step in years
  start_year = 2020,   # Start year
  end_year = 2100      # End year
)
```

**Common modifications:**
- **Shorter horizon:** Set `end_year = 2050` for near-term analysis
- **Mid-century focus:** Set `start_year = 2025`, `end_year = 2075`
- **Coarser resolution:** Set `dt = 2` for biennial time steps

---

### 8.2 CO₂ Target

**Where:** `parameter_details_fixed.yml`
```yaml
co2_target_2100: 650  # CO2 emissions volume at 2100 to reach 1.5°C target
```

**Note:** This represents the target cumulative emissions (GtCO₂) by 2100, which corresponds to a specific temperature target based on the TCRE (Transient Climate Response to Emissions) parameter.

**Common values:**
- Lower targets may be infeasible under some scenarios/delays
- Higher targets correspond to higher temperature outcomes

---

### 8.3 Discount Rate

**Where:** `parameter_details_fixed.yml`
```yaml
disc_rate: 0.03  # Annual discount rate (3%)
```

**Common values:**
- `0.02` (2%): Low discounting—future prioritised
- `0.03` (3%): Middle ground (default)
- `0.05` (5%): Higher discounting—present prioritised

**Effect:** Higher discount rates → less weight on future costs → less aggressive early action

---

### 8.4 Cost Parameters

**Where:** `parameter_details_fixed.yml`

```yaml
# Mitigation cost curve parameters
vol_mitig_low: 1      # Volume at which cost_mitig_low applies
cost_mitig_low: 0.01  # Cost per Gt at low volume ($ trillion)
vol_mitig_peak: 50    # Volume at which cost_mitig_peak applies
cost_mitig_peak: 1    # Cost per Gt at peak volume ($ trillion)

# CDR cost curve parameters  
vol_remov_low: 1      # Volume at which cost_remov_low applies
cost_remov_low: 0.01  # Cost per Gt at low volume ($ trillion)
vol_remov_peak: 50    # Volume at which cost_remov_peak applies
cost_remov_peak: 2    # Cost per Gt at peak volume ($ trillion)

# Economic damage coefficient
econ_dam_pct: 0.05    # Proportion of GWP reduced by climate change
```

**Note:** These parameters define quadratic cost curves. The costs increase from low to peak values as volumes increase. The economic damage coefficient determines the relationship between temperature and economic damages.

---

### 8.5 Parallel Processing

**Where:** Any analysis chunk (`scenario_comparison_analysis`, `parameter_importance_analysis`, `delayed_deployment_analysis`)

```r
use_parallel = TRUE   # Enable parallel processing
n_cores = NULL        # NULL = auto-detect (uses all cores - 1)
```

**When to enable:**
- Multiple scenarios
- Large parameter sets (>100 samples)
- Delayed deployment with fine grid

**When to disable:**
- Debugging
- Single quick test
- Limited RAM

**Note:** Parallel processing requires the `parallel` and `doParallel` packages.

---

### 8.6 Output Saving

**Where:** Analysis chunks

```r
save_results = TRUE           # Automatically save results
output_dir = "output"         # Directory for saved files
output_prefix = "analysis"    # Filename prefix
```

**Saved files:**
- `.rds`: Complete R object (for later analysis)
- `.csv`: Summary table (for Excel/other software)

**File naming:** `{prefix}_{timestamp}.{ext}`

Example: `scenario_comparison_20260202_143022.rds`

---

## Understanding Outputs

### 9.1 Scenario Comparison Outputs

**Main results object:** `scenario_comparison_results`

**Key components:**

#### `scenario_results`
List containing full results for each scenario:
- `time_series`: Year-by-year values
  - `year`: Year
  - `qty_mitig`: Mitigation quantity (GtCO₂/year)
  - `qty_remov`: CDR quantity (GtCO₂/year)
  - `cumulative_emissions`: Cumulative CO₂ (GtCO₂)
  - `temperature_anomaly`: Temperature above pre-industrial (°C)
  - `annual_cost`: Cost that year ($ trillion)

#### `comparison_summary`
Data frame comparing scenarios:
- `scenario`: Scenario name
- `peak_temperature`: Maximum temperature (°C)
- `year_peak_temp`: Year of peak temperature
- `total_cost`: Total discounted cost ($ trillion)
- `mitig_cost`: Mitigation costs ($ trillion)
- `remov_cost`: CDR costs ($ trillion)
- `temp_cost`: Temperature damage costs ($ trillion)
- `total_mitigation`: Total mitigation deployed (GtCO₂)
- `total_cdr`: Total CDR deployed (GtCO₂)

#### `year_first_1p5C`
Named vector: First year each scenario crosses 1.5°C

#### `year_peak_temp`
Named vector: Year of maximum temperature for each scenario

---

### 9.2 Parameter Importance Outputs

**Main results object:** `parameter_importance_results`

**Key components:**

#### `summary`
Data frame with one row per parameter combination:
- Parameter values (one column per parameter)
- Outcome variables (temperature, costs, etc.)
- Convergence flags

#### `importance_summary`
Ranked list of parameter influences:
- Parameter name
- Correlation with each outcome
- Absolute importance score
- Rank

#### `correlation_matrix`
Full correlation matrix between:
- All parameters (inputs)
- All outcome variables (outputs)

#### `successful_runs` and `failed_runs`
Lists of individual parameter set results for detailed inspection

---

### 9.3 Delayed Deployment Outputs

**Main results object:** `delayed_deployment_results`

**Key components:**

#### `results_by_scenario`
List with one element per scenario, each containing:
- `results`: Full grid of delay combinations
- `feasible_results`: Subset of feasible combinations only
- `summary_stats`: Convergence and feasibility rates

#### `combined_results`
Single data frame combining all scenarios:
- `scenario`: Full scenario name (e.g., "SSP2-Baseline")
- `scenario_short`: Abbreviated name (e.g., "SSP2")
- `mitigation_delay`: Years of mitigation delay
- `cdr_delay`: Years of CDR delay
- `feasible`: Logical (TRUE if temperature target met)
- `peak_temperature`: Peak temperature (°C)
- `years_above_1p5`: Years exceeding 1.5°C
- `total_cost`: Total cost ($ trillion)
- `mitig_cost`: Mitigation cost ($ trillion)
- `remov_cost`: CDR cost ($ trillion)
- `temp_cost`: Temperature damage cost ($ trillion)

#### `summary_stats`
- `n_scenarios_successful`: Number of scenarios completed
- `total_combinations`: Total delay combinations tested
- `total_feasible`: Number of feasible combinations
- `overall_feasibility_rate`: Proportion feasible

---

### 9.4 Cost Variable Names in the Code

**IMPORTANT:** Be aware of the following cost variable names used in the actual code:

1. **`mitig_cost`**: Mitigation costs only
2. **`remov_cost`**: CDR/removal costs only  
3. **`temp_cost`**: Temperature-related damage costs
4. **`total_cost`**: Sum of mitigation, CDR, and temperature damage costs
5. **`abatement_cost`**: Combined mitigation + CDR costs (calculated in visualization code)

The visualization functions may calculate `abatement_cost` as the sum of `mitig_cost` and `remov_cost` for convenience in certain plots.

---

## Troubleshooting

### 10.1 Installation Issues

#### Problem: "there is no package called 'X'"
**Solution:**
```r
install.packages("X")
```

#### Problem: Package won't install
**Solutions:**
1. Update R: Check you have R ≥ 4.0.0
2. Try CRAN mirror:
```r
install.packages("X", repos = "https://cloud.r-project.org/")
```
3. Install dependencies first:
```r
install.packages("X", dependencies = TRUE)
```

#### Problem: "non-zero exit status" when installing
**Solution:** On Linux, install system dependencies:
```bash
# Ubuntu/Debian
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev

# Red Hat/Fedora
sudo yum install libcurl-devel openssl-devel libxml2-devel
```

---

### 10.2 Data Issues

#### Problem: "cannot open file 'data/emissions.csv'"
**Solutions:**
1. Verify file exists: `file.exists("data/emissions.csv")`
2. Check working directory: `getwd()`
3. Set correct directory: `setwd("path/to/MACROM")`

#### Problem: Missing emissions scenarios
**Check available scenarios:**
```r
unique(emissions_df$Scenario)
```
**Solution:** Use only scenarios that exist in your data

---

### 10.3 Runtime Issues

#### Problem: Analysis taking too long
**Solutions:**
1. Enable parallel processing: `use_parallel = TRUE`
2. Reduce parameter samples (for parameter importance)
3. Increase delay step size (for delayed deployment)
4. Test with single scenario first

#### Problem: Parallel processing not working
**Solutions:**
1. Check packages installed:
```r
library(parallel)
library(doParallel)
```
2. Check cores available: `parallel::detectCores()`
3. Set cores manually: `n_cores = 4`
4. Try sequential processing: `use_parallel = FALSE`

---

### 10.4 Memory Issues

#### Problem: "cannot allocate vector of size X"
**Solutions:**
1. Close other applications
2. Reduce number of parameter samples
3. Use fewer scenarios simultaneously
4. Increase system RAM (if possible)
5. Process scenarios one at a time:
```r
# Instead of all scenarios at once:
scenarios_to_compare <- "SSP2-Baseline"
# Run, save results, then next scenario
```

---

### 10.5 Convergence Issues

#### Problem: "Optimal control did not converge"
**Possible causes:**
1. Infeasible parameter combination
2. Too stringent temperature target
3. Excessive delays (delayed deployment)

**Solutions:**
1. Check if target is achievable for scenario
2. Relax convergence tolerance
3. Modify parameter bounds
4. Reduce delays or step size

#### Problem: Many failed runs in parameter importance
**Normal if:**
- Using wide parameter ranges
- Some combinations are infeasible
- Typically 10-20% failures acceptable

**Concerning if:**
- >50% failures
- All runs failing

**Solutions:**
1. Check parameter ranges are reasonable
2. Verify scenarios exist
3. Test with fixed parameters first

---

### 10.6 Visualisation Issues

#### Problem: Plots look wrong or empty
**Solutions:**
1. Check results object exists: `str(scenario_comparison_results)`
2. Verify results structure matches expected format
3. Re-run analysis chunk
4. Check for error messages in console

#### Problem: "Cannot find variable 'X' in data"
**Solution:**  
Check variable names match between analysis and visualization:
- Use `names(delayed_deployment_results$combined_results)` to see available variables
- Ensure you're using the correct variable names:
  - `temp_cost` (NOT `temp_damage_cost`)
  - `mitig_cost`, `remov_cost`, `total_cost`

#### Problem: Saved plots have wrong dimensions
**Solutions:**
```r
# Adjust width and height
width = 297   # mm (A4 landscape)
height = 210  # mm

# Or for US Letter
width = 11 * 25.4  # 11 inches to mm
height = 8.5 * 25.4
```

---

### 10.7 Getting More Help

If issues persist:

1. **Check console output** for error messages
2. **Run sessionInfo()** to check R version and packages:
```r
sessionInfo()
```
3. **Create minimal reproducible example:**
```r
# Simplest possible test
scenarios_to_compare <- "SSP2-Baseline"
use_parallel = FALSE
verbose = TRUE
```
---

## Frequently Asked Questions

### General Questions

#### Q: What does MACROM stand for?

**A:** MACROM stands for **Mitigation and Carbon Removal Optimisation Model**. It's an optimal control model for finding cost-effective pathways to climate stabilisation.

#### Q: Is MACROM free to use?

**A:** Yes! MACROM is open source under the CC-BY-4.0 licence. You can use, modify, and distribute it freely, with attribution.

#### Q: What programming experience do I need?

**A:** You need:
- **Basic R skills** (loading data, running scripts, making plots)
- **RStudio familiarity** (or R command line)
- **No optimal control theory required** (helpful but not essential)

Most users succeed with just introductory R knowledge.

#### Q: Can MACROM model specific CDR technologies?

**A:** No. MACROM treats CDR as homogeneous (doesn't distinguish between DACCS, afforestation, etc.). Technology-specific analysis requires code modifications.


#### Q: How accurate are MACROM's cost estimates?

**A:** Cost estimates are **relative, not absolute**. They:
- Show comparative costs across scenarios
- Identify cost-optimal strategies
- Depend heavily on cost function parameters (which are uncertain)
- Useful for comparing strategies, not predicting exact costs
- Should be validated against other models and expert elicitation

Use MACROM for relative comparisons and understanding trade-offs, not absolute cost predictions.

#### Q: Why doesn't MACROM recommend [specific policy]?

**A:** MACROM is a cost-optimisation model, so:
- It finds economically optimal strategies
- Does not consider political feasibility
- Ignores equity, justice, and distributional concerns
- Simplifies technological and social constraints

Results should inform, not dictate, policy decisions. Real policy must consider factors beyond economic optimality.

---

### Analysis and Workflow

#### Q: Which analysis should I run first?

**A:** For most users:
1. **Start with scenario comparison** using fixed parameters and SSP2-Baseline
2. Review outputs to understand model behaviour
3. Then expand to all scenarios or add other analyses

This approach:
- Gives quick results (~10 minutes)
- Helps you understand outputs before tackling longer analyses
- Reveals any data or setup issues early

#### Q: Do I need to run chunks in order?

**A:** Yes, with some flexibility:
- **Always run first**: Libraries, Data Preparation, Parameter Selection
- **Then run in any order**: Analysis chunks (scenario comparison, parameter importance, delayed deployment)
- **After each analysis**: Corresponding visualisation chunk

Dependencies are documented in each chunk.

#### Q: How do I reproduce results from a previous run?

**A:**
1. **Use same random seed** (for LHS):
```r
   generate_lhs_samples(n_samples = 500, seed = 12345)
```

2. **Save parameter values** used:
```r
   write.csv(parameter_df, "output/parameters_used.csv")
```

3. **Use version control** (git) to track code changes

4. **Document R version and packages**:
```r
   writeLines(capture.output(sessionInfo()), "output/session_info.txt")
```

#### Q: Can I stop and restart a long-running analysis?

**A:** Not easily mid-analysis, but you can:
- **Save intermediate results** if modifying source code
- **Run scenarios individually** and combine later
- **Use checkpointing** (requires custom code)

Best practice: Run long analyses overnight or on HPC.

---

### Customisation and Extension

#### Q: How do I change the CO₂ target?

**A:** Edit `parameter_details_fixed.yml`:
```yaml
co2_target_2100: 650  # Change to desired cumulative emissions target (GtCO₂)
```

This target corresponds to a specific temperature outcome based on your TCRE parameter. Then re-run `fixed_model_parameters_call` chunk.

#### Q: Can I model specific CDR technologies (e.g., DACCS vs afforestation)?

**A:** The current version treats CDR as homogeneous. To model specific technologies:
- Modify cost functions in source code
- Add separate state variables for each technology
- Adjust constraints for technology-specific limits

This requires advanced modifications.

#### Q: How do I add a new cost function?

**A:** This requires modifying `optimal_control_core.R`:
1. Define new cost function
2. Add to objective function in solver
3. Add parameters to `parameter_details_fixed.yml`
4. Test thoroughly on simple cases


#### Q: Can MACROM handle regional analysis?

**A:** The current version is global. Regional analysis would require:
- Regional emission and economic data
- Regional climate response (more complex)
- Inter-regional trade and cooperation mechanisms
- Significant model extensions

This is beyond current scope.

---

### Comparison with Other Models

#### Q: How does MACROM compare to IAMs (Integrated Assessment Models)?

**A:** 

**MACROM advantages:**
- Fast computation (minutes vs hours/days)
- Transparent code (fully open source)
- Focused on specific question (optimal mitigation/CDR timing)
- Easy to modify and extend

**IAM advantages:**
- More comprehensive (energy, land use, multiple gases)
- More detailed economic sectors
- Established for policy analysis
- Extensive validation and comparison

**Best use:** MACROM for rapid exploration of optimal control strategies; IAMs for comprehensive policy assessment.

#### Q: Should I use MACROM or DICE/RICE model?

**A:**

**Use MACROM if:**
- Focused on mitigation vs CDR timing
- Need fast iteration and sensitivity analysis
- Want transparent, modifiable code
- Exploring optimal control methodology

**Use DICE/RICE if:**
- Need established model with literature precedent
- Require detailed economic sectors
- Policy analysis requiring model recognition
- Want carbon pricing focus

They're complementary—MACROM can inform DICE/RICE parameter choices and vice versa.

#### Q: Can MACROM results be compared with IPCC scenarios?

**A:** Yes, but carefully:
- MACROM shows optimal strategies (normative)
- IPCC scenarios show possible pathways (descriptive)
- MACROM uses simpler climate model
- Different cost assumptions

Best approach: Compare qualitative patterns (e.g., timing of peak emissions, CDR scale) rather than exact numbers.

---

### Getting Help and Contributing

#### Q: Where can I get help with MACROM?

**A:**
1. This User Guide (you're reading it!)
2. Code comments in `MACROM_workflow.Rmd`
3. GitHub issues: https://github.com/nina-rynne/MACROM/issues
4. Email: nina.rynne@griffithuni.edu.au

Please search existing issues before creating new ones.

#### Q: How can I contribute to MACROM?

**A:** Contributions welcome! See README.md "Contributing" section. Options include:
- Bug reports and fixes
- Documentation improvements
- New features or analysis types
- Testing and validation
- Example analyses and tutorials

Use GitHub pull requests for code contributions.

#### Q: Is MACROM suitable for my undergraduate/graduate research project?

**A:** Yes, if your project involves:
- Climate policy analysis
- Optimal control applications
- Cost-benefit analysis
- Scenario comparison

Requirements:
- Basic R programming skills
- Understanding of climate change basics
- Time for learning curve (~5-10 hours)
- Access to suitable computing resources

Suitable for Master's theses, PhD chapters, coursework projects.

#### Q: Can I use MACROM results in publications?

**A:** Yes! MACROM is open source (CC-BY-4.0 licence). Please:
1. **Cite the MACROM paper** (see README.md)
2. **Acknowledge the model**: "Analysis performed using MACROM"
3. **Document parameters used**: Include in supplementary materials
4. **Share your code** (if possible): Helps reproducibility

#### Q: Will MACROM continue to be developed?

**A:** Yes!

Check GitHub repository for latest updates and roadmap.

---

### Still Have Questions?

If your question isn't answered here:

1. **Search GitHub issues**: https://github.com/nina-rynne/MACROM/issues
2. **Email the developers**: nina.rynne@griffithuni.edu.au
3. **Create a new GitHub issue**: Include:
   - Clear description of question
   - What you've tried
   - Relevant code/outputs
   - `sessionInfo()` output

We aim to respond within 1-2 weeks.

---

**End of User Guide**

Thank you for using MACROM!