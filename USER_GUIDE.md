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
11. [Advanced Usage](#advanced-usage)
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
├── MACROM_workflow.Rmd
├── README.md
├── USER_GUIDE.md
├── LICENSE
├── data/
│   ├── emissions.csv
│   └── gwp.csv
├── figs/              (will be created automatically)
├── output/            (will be created automatically)
└── src/
    ├── data_preparation.R
    ├── model_parameters.R
    ├── latin_hypercube_sampling.R
    ├── optimal_control_core_V3.R
    ├── scenario_comparison.R
    ├── scenario_comparison_visualisation.R
    ├── parameter_importance.R
    ├── parameter_importance_visualisation.R
    ├── delayed_deployment.R
    └── delayed_deployment_visualisation.R
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
- Parameter sensitivity (500 samples, 1 scenario): 1-3 hours
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
├─ "How robust are results to parameter uncertainty?"
│   → Use: Latin Hypercube Sampling + Scenario Comparison
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
- Imports and processes SSP emission and economic data
- Should complete in under 30 seconds
- Creates `emissions_df` and `economic_df` objects

#### 3. Choose Parameter Approach (CHOOSE ONE)

**For most users (RECOMMENDED):**

**Chunk: `fixed_model_parameters_call`**
- Uses single best-estimate parameter set
- Fast and suitable for most analyses
- Creates `parameter_df` object with 1 row

**For sensitivity analysis:**

**Chunk: `latin_hypercube_sampling_call`** (set `eval=TRUE` first)
- Generates 500 parameter combinations
- Much longer runtime (hours instead of minutes)
- Creates `parameter_df` object with 500 rows

#### 4. Run Analysis (CHOOSE ONE OR MORE)

**Option A: Scenario Comparison** (RECOMMENDED FOR FIRST RUN)

**Chunk: `scenario_comparison_analysis`** (set `eval=TRUE` first)
- Compares optimal strategies across all 5 SSP scenarios
- Runtime: ~10-15 minutes with parallel processing
- Creates `scenario_comparison_results` object

**Chunk: `scenario_comparison_visualisation`** (set `eval=TRUE` first)
- Generates publication-quality figures
- Runtime: ~1-2 minutes
- Saves figures to `figs/` directory

#### 5. Check Your Results

After running successfully, you should see:

**In your R environment:**
- `emissions_df`, `economic_df`, `parameter_df`
- `scenario_comparison_results` (or other results objects)

**In the `output/` directory:**
- `scenario_comparison_TIMESTAMP.rds`: Full results
- `scenario_comparison_summary_TIMESTAMP.csv`: Summary table

**In the `figs/` directory:**
- `scenario_dashboard_TIMESTAMP.pdf`: Multi-panel figure

### What Next?

- **Explore results**: Open the PDF figures to view your analysis
- **Load saved data**: Use `readRDS()` to reload results for further analysis
- **Try other analyses**: Run parameter importance or delayed deployment chunks
- **Customise**: Modify parameters in the chunks (see Customisation Options section)

### Quick Troubleshooting

**Problem:** Chunk fails with "could not find function X"  
**Solution:** Ensure Libraries chunk has been run

**Problem:** "object not found" error  
**Solution:** Run previous chunks in order—each chunk depends on the previous ones

**Problem:** Figures don't appear  
**Solution:** Check the `figs/` directory—they're saved automatically even if not displayed

---

## Detailed Workflow Instructions

This section provides in-depth explanations for each workflow component.

### 6.1 Libraries Chunk

**Purpose:** Loads all required R packages for data manipulation, analysis, and visualisation.

**What it does:**
- Loads 20+ R packages
- Organises packages by function (data manipulation, file handling, visualisation, etc.)

**If this chunk fails:**
1. Check which package failed to load (error message will specify)
2. Install the missing package: `install.packages("package_name")`
3. Re-run the chunk

**Common issues:**
- Package not installed → Install it
- Package version too old → Update it: `update.packages("package_name")`

---

### 6.2 Data Preparation Chunk

**Purpose:** Import and process SSP emission and economic data.

**What it does:**
1. Imports raw SSP data from CSV files
2. Interpolates data to annual time steps (2020-2100)
3. Creates `emissions_df` and `economic_df` data frames

**Customisation options:**
```r
# Change time range
emissions_df <- interpolate_ssp_emissions(
  emissions_df = emissions_imported,
  dt = 1,              # Time step: 1 = annual, 0.5 = 6-month, etc.
  start_year = 2020,   # Start year (default: 2020)
  end_year = 2100      # End year (default: 2100)
)
```

**When to customise:**
- Shorter time horizon: Change `end_year` to 2050 or 2075
- Finer temporal resolution: Change `dt` to 0.5 for semi-annual steps
- Different baseline: Change `start_year` (must match data availability)

**Output structure:**

`emissions_df` contains:
- `Year`: Time steps from start_year to end_year
- `Scenario`: SSP scenario name (SSP1-Baseline, SSP2-Baseline, etc.)
- `Emissions`: CO₂ emissions in GtCO₂/year

`economic_df` contains:
- `Year`: Time steps from start_year to end_year
- `Scenario`: SSP scenario name
- `GWP`: Gross World Product in trillion USD

---

### 6.3 Parameter Selection

You must choose ONE of these two approaches.

#### Option A: Fixed Parameters (Recommended)

**Chunk:** `fixed_model_parameters_call`

**Purpose:** Use a single best-estimate parameter set.

**When to use:**
- Comparing across SSP scenarios
- Analysing delayed deployment effects
- Quick test runs
- Generating primary results for publication

**What it does:**
- Reads parameters from `parameter_details_fixed.yml`
- Creates `parameter_df` with 1 row

**How to modify parameters:**
1. Open `parameter_details_fixed.yml` in a text editor
2. Modify parameter values
3. Save the file
4. Re-run this chunk

**Key parameters you might want to modify:**
- `temp_target`: Target temperature (default: 1.5°C)
- `discount_rate`: Discount rate for economic costs (default: 0.02)
- `mitigation_max_rate`: Maximum annual mitigation rate
- `cdr_cumulative_limit`: Total CDR capacity constraint

#### Option B: Latin Hypercube Sampling

**Chunk:** `latin_hypercube_sampling_call`

**Purpose:** Generate multiple parameter combinations for sensitivity analysis.

**When to use:**
- Exploring parameter uncertainty
- Identifying influential parameters
- Generating confidence intervals
- Comprehensive uncertainty quantification

**What it does:**
- Samples parameter space using Latin Hypercube Sampling
- Creates `parameter_df` with `n_samples` rows (default: 500)

**Customisation:**
```r
# Change number of samples
lhs_parameter_df <- generate_lhs_samples(
  n_samples = 500,    # Number of parameter combinations (100-1000 typical)
  seed = 12345        # Random seed for reproducibility
)
```

**Runtime implications:**
- 100 samples: ~30 minutes per scenario
- 500 samples: ~2-3 hours per scenario
- 1000 samples: ~5-6 hours per scenario

**Parameter ranges:**
- Defined in `parameter_details.yml`
- Each parameter has a minimum and maximum value
- LHS samples uniformly across this range

---

### 6.4 Scenario Comparison Analysis

**Chunk:** `scenario_comparison_analysis`

**Purpose:** Compare optimal strategies across different SSP baseline scenarios.

**Requirements:**
- Must use single parameter set (`parameter_df` with 1 row)
- Recommended: Run `fixed_model_parameters_call` first

**What it does:**
1. Runs optimal control algorithm for each SSP scenario
2. Solves for optimal mitigation and CDR trajectories
3. Calculates temperature outcomes and costs
4. Compares results across scenarios

**Customisation options:**
```r
# Select scenarios to compare
scenarios_to_compare <- c("SSP1-Baseline", "SSP2-Baseline", "SSP3-Baseline", 
                          "SSP4-Baseline", "SSP5-Baseline")

# Or compare subset
scenarios_to_compare <- c("SSP2-Baseline", "SSP3-Baseline")

# Delay settings (for testing delayed start)
mitigation_delay <- 0   # Years to delay mitigation (0 = start immediately)
cdr_delay <- 0          # Years to delay CDR (0 = start immediately)

# Processing options
use_parallel = TRUE     # Use parallel processing (faster)
verbose = FALSE         # Print detailed progress (set TRUE for debugging)
```

**Output:**

`scenario_comparison_results` contains:
- `scenario_results`: Full time series for each scenario
- `comparison_summary`: Summary table comparing key metrics
- `year_first_1p5C`: Year when temperature first reaches 1.5°C
- `year_peak_temp`: Year of peak temperature
- Runtime and configuration metadata

**Expected runtime:**
- 1 scenario: 2-3 minutes
- 5 scenarios (parallel): 10-15 minutes
- 5 scenarios (sequential): 30-45 minutes

---

### 6.5 Parameter Importance Analysis

**Chunk:** `parameter_importance_analysis`

**Purpose:** Identify which parameters have the greatest influence on model outcomes.

**Requirements:**
- Must use Latin Hypercube Sampling (`parameter_df` with multiple rows)
- Run `latin_hypercube_sampling_call` first

**When to use:**
- Understanding model sensitivity
- Identifying key uncertainties
- Prioritising parameters for calibration
- Communicating model behavior to stakeholders

**What it does:**
1. Runs optimal control for all parameter combinations
2. Calculates correlations between parameters and outcomes
3. Ranks parameters by influence on key metrics (temperature, costs, etc.)

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
- `importance_summary`: Ranked list of parameter influences
- `correlation_matrix`: Full correlation matrix between parameters and outcomes
- `sensitivity_indices`: Quantitative sensitivity measures
- Individual results for all parameter combinations

**Expected runtime:**
- 100 samples: ~30-45 minutes with parallel processing
- 500 samples: ~2-3 hours with parallel processing
- 1000 samples: ~5-6 hours with parallel processing

**Interpreting results:**

High importance means:
- Changes in this parameter strongly affect outcomes
- Uncertainty in this parameter translates to uncertainty in results
- Careful calibration/estimation of this parameter is critical

Low importance means:
- Model is relatively insensitive to this parameter
- Can use wider uncertainty ranges without major impact
- Less critical for model calibration

---

### 6.6 Parameter Importance Visualisation

**Chunk:** `parameter_importance_visualisation`

**Purpose:** Create visualisations showing parameter sensitivities and relationships.

**What it does:**
- Generates sensitivity plots (tornado diagrams, scatter plots)
- Shows parameter correlations with outcomes
- Identifies nonlinear relationships
- Creates dashboard summarising key findings

**Customisation options:**
```r
# Select outcome variables to analyse
outcome_variables = c("peak_temperature", "total_cost", "years_above_1p5")

# Available outcome variables:
#   - peak_temperature: Maximum temperature anomaly
#   - total_cost: Sum of abatement and damage costs
#   - years_above_1p5: Years exceeding 1.5°C threshold
#   - years_above_2: Years exceeding 2°C threshold
#   - total_mitigation: Cumulative mitigation deployed
#   - total_cdr: Cumulative CDR deployed
#   - abatement_cost: Total mitigation + CDR costs
#   - temp_damage_cost: Temperature-related damages

# Styling options
save_plot = TRUE        # Save to file
width = 297             # Width in mm
height = 210            # Height in mm
```

**Output figures show:**
1. **Tornado diagram**: Parameters ranked by influence
2. **Scatter plots**: Relationships between parameters and outcomes
3. **Correlation heatmap**: Full parameter correlation matrix
4. **Partial dependence plots**: How outcomes change with each parameter

---

### 6.7 Delayed Deployment Analysis

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
max_delay_years <- 50        # Maximum delay to test (years)
delay_step_size <- 10        # Step size (years between tested delays)

# Examples of step size effects:
# step_size = 1:  Tests 0, 1, 2, 3, ... 50 (51 × 51 = 2,601 combinations)
# step_size = 5:  Tests 0, 5, 10, 15, ... 50 (11 × 11 = 121 combinations)
# step_size = 10: Tests 0, 10, 20, 30, 40, 50 (6 × 6 = 36 combinations)

# Scenarios to analyse
scenarios_to_analyse <- "all"  # Test all 5 SSP scenarios

# Or select specific scenarios
scenarios_to_analyse <- c("SSP2-Baseline", "SSP3-Baseline")

# Or single scenario for quick test
scenarios_to_analyse <- "SSP2-Baseline"

# Processing options
use_parallel = TRUE     # Highly recommended for this analysis
```

**Computational considerations:**

Total runs = (max_delay / step_size + 1)² × number of scenarios

Examples:
- 50-year max, 10-year steps, 1 scenario: 36 runs (~5 minutes)
- 50-year max, 10-year steps, 5 scenarios: 180 runs (~30 minutes)
- 50-year max, 5-year steps, 5 scenarios: 605 runs (~2 hours)
- 50-year max, 1-year steps, 5 scenarios: 13,005 runs (~20+ hours)

**Recommendation:** Start with `step_size = 10` and `scenarios_to_analyse = "SSP2-Baseline"` for quick testing.

**Output:**

`delayed_deployment_results` contains:
- `results`: Full results for all delay combinations
- `summary_by_scenario`: Aggregated metrics per scenario
- `feasibility_map`: Which delay combinations are feasible
- Grid of outcomes (temperature, costs) across delay space

**Interpreting results:**

- **Feasible region**: Combinations where temperature target can be met
- **Infeasible region**: Delays too long to achieve target
- **Cost gradients**: How costs increase with delays
- **Critical thresholds**: Maximum delays before infeasibility

---

### 6.8 Delayed Deployment Visualisation

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
#   - abatement_cost: Total abatement costs (mitigation + CDR)
#   - temp_damage_cost: Temperature-related damage costs
#   - total_cost: Sum of abatement and damage costs
#   - mitig_cost: Mitigation costs only
#   - remov_cost: CDR/removal costs only

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
- **Upper right corner (max,max)**: Maximum delays for both technologies
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

### 6.9 Pre-built Visualisation Functions

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
  cost_variable = "abatement_cost",
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
- `"abatement_cost"`: Total mitigation + CDR costs
- `"temp_damage_cost"`: Climate damages
- `"total_cost"`: Sum of abatement and damages

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
- Three-panel comparison: total cost, abatement cost, damage cost
- Shared colour scale for direct magnitude comparison

**Use when:** Understanding cost composition and trade-offs

**Note:** Set `shared_scale = FALSE` for independent scales if cost components vary widely in magnitude

---

## Analysis Types

This section summarises the three main analysis types and when to use each.

### 7.1 Scenario Comparison

**Research Question:** "How do optimal mitigation and CDR strategies differ across socio-economic futures?"

**What it does:**
- Compares optimal control solutions across SSP1-5 baseline scenarios
- Shows how baseline emissions affect optimal strategies
- Identifies robust vs scenario-dependent strategies

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
**Runtime:** 2-3 hours (500 samples, 1 scenario, parallel)

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
**Runtime:** 30-60 minutes (5 scenarios, 10-year steps, parallel)

---

### 7.4 Recommended Analysis Sequences

**For publication (comprehensive analysis):**
1. Scenario Comparison (main results)
2. Parameter Importance (sensitivity/uncertainty)
3. Delayed Deployment (policy implications)

**For quick exploration:**
1. Fixed Parameters + Scenario Comparison (SSP2 only)
2. Review outputs and decide next steps

**For policy analysis:**
1. Fixed Parameters + Delayed Deployment (all SSPs)
2. Focus on cost penalties and feasibility boundaries

**For model development:**
1. Latin Hypercube Sampling + Parameter Importance
2. Identify key parameters for calibration
3. Refine parameter ranges and re-run

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
- **Finer resolution:** Set `dt = 0.5` for semi-annual time steps (increases computation)

---

### 8.2 Temperature Target

**Where:** `parameter_details_fixed.yml`
```yaml
temp_target:
  value: 1.5
  description: "Target temperature anomaly (°C)"
```

**Common values:**
- `1.5`: Paris Agreement aspirational target
- `2.0`: Paris Agreement upper limit
- `1.75`: Intermediate target

**Note:** Lower targets may be infeasible under some scenarios/delays

---

### 8.3 Discount Rate

**Where:** `parameter_details_fixed.yml`
```yaml
discount_rate:
  value: 0.02
  description: "Annual discount rate for economic costs"
```

**Typical range:** 0.01 (1%) to 0.05 (5%)

**Effect:**
- Lower rates → Future costs matter more → Earlier action
- Higher rates → Future costs matter less → Later action

---

### 8.4 Technology Constraints

**Where:** `parameter_details_fixed.yml`

Key constraints:
```yaml
mitigation_max_rate:
  value: 0.10
  description: "Maximum annual mitigation rate (fraction of baseline)"

cdr_cumulative_limit:
  value: 1000
  description: "Total CDR capacity constraint (GtCO2)"

cdr_max_annual_rate:
  value: 20
  description: "Maximum annual CDR deployment (GtCO2/year)"
```

**When to modify:**
- Testing technology scenarios (optimistic vs pessimistic)
- Exploring feasibility boundaries
- Sensitivity to technology assumptions

---

### 8.5 Scenario Selection

**Where:** Analysis chunks (`scenario_comparison_analysis`, `delayed_deployment_analysis`)
```r
# All scenarios
scenarios_to_compare <- c("SSP1-Baseline", "SSP2-Baseline", "SSP3-Baseline", 
                          "SSP4-Baseline", "SSP5-Baseline")

# Subset for testing
scenarios_to_compare <- c("SSP2-Baseline")

# High vs low emissions
scenarios_to_compare <- c("SSP1-Baseline", "SSP5-Baseline")
```

**SSP characteristics:**
- **SSP1:** Sustainable development, low emissions
- **SSP2:** Middle-of-road, moderate emissions
- **SSP3:** Regional rivalry, high emissions
- **SSP4:** Inequality, moderate-high emissions
- **SSP5:** Fossil-fueled development, highest emissions

---

### 8.6 Visualisation Styling

**Colour palettes:**
```r
color_palette = "viridis"    # Default
color_palette = "plasma"     # Alternative (purple-orange)
color_palette = "inferno"    # Alternative (dark-red-yellow)
color_palette = "magma"      # Alternative (dark-purple-pink)
```

**Plot dimensions:**
```r
width = 297     # A4 landscape width (mm)
height = 210    # A4 landscape height (mm)

# For multi-panel layouts:
height = 260    # 3-row layout
height = 297    # A4 portrait height
```

**DPI (resolution):**
```r
dpi = 300       # Publication quality (default)
dpi = 150       # Screen viewing (smaller files)
dpi = 600       # High-resolution print
```

---

### 8.7 Parallel Processing

**Where:** All analysis chunks
```r
use_parallel = TRUE     # Use multiple CPU cores (faster)
use_parallel = FALSE    # Use single core (slower but more stable)
```

**When to disable:**
- Debugging (easier to see error messages)
- Memory constraints (parallel uses more RAM)
- Small problems (overhead not worth it)

**Performance impact:**
- Scenario comparison: 3-5× speedup with 4-8 cores
- Parameter importance: 4-8× speedup with 8+ cores
- Delayed deployment: 3-5× speedup with 4-8 cores

---

### 8.8 Verbosity

**Where:** All analysis chunks
```r
verbose = FALSE     # Minimal output (default)
verbose = TRUE      # Detailed progress messages
```

**Set to TRUE when:**
- Debugging issues
- Long-running jobs (monitor progress)
- First time running analysis (understand what's happening)

**Set to FALSE when:**
- Production runs (cleaner output)
- Generating reports (less clutter)

---

## Understanding Outputs

### 9.1 File Locations

MACROM creates two output directories:

**`output/` directory:**
- Contains data files (RDS, CSV)
- Full results objects with all time series and metadata
- Summary statistics tables

**`figs/` directory:**
- Contains visualisations (PDF, PNG)
- Publication-ready figures
- Dashboards and multi-panel plots

**File naming convention:**
```
{analysis_type}_{timestamp}.{extension}

Examples:
scenario_comparison_20260108_143022.rds
scenario_comparison_summary_20260108_143022.csv
scenario_dashboard_20260108_143045.pdf
```

Timestamp format: `YYYYMMDD_HHMMSS` (year, month, day, hour, minute, second)

---

### 9.2 Scenario Comparison Outputs

#### RDS File (`scenario_comparison_TIMESTAMP.rds`)

Full results object containing:
```r
results <- readRDS("output/scenario_comparison_20260108_143022.rds")

# Structure:
results$scenario_results          # List of results for each scenario
results$comparison_summary        # Summary table (see CSV)
results$year_first_1p5C          # Year reaching 1.5°C (per scenario)
results$year_peak_temp           # Year of peak temperature (per scenario)
results$year_mitig_capped        # Year when mitigation hits emissions limit
results$failed_scenarios         # Any scenarios that failed to solve
results$run_info                 # Metadata (runtime, parameters, etc.)
```

**Per-scenario results structure:**
```r
results$scenario_results$`SSP2-Baseline`$solution

# Contains time series:
$year                  # Time steps
$temperature          # Temperature anomaly (°C)
$emissions            # Controlled emissions (GtCO2/yr)
$mitigation           # Mitigation deployment (GtCO2/yr)
$cdr                  # CDR deployment (GtCO2/yr)
$cumulative_emissions # Cumulative CO2 (GtCO2)
$mitigation_cost      # Annual mitigation cost (trillion USD)
$cdr_cost            # Annual CDR cost (trillion USD)
$damage_cost         # Annual damage cost (trillion USD)
```

#### CSV File (`scenario_comparison_summary_TIMESTAMP.csv`)

Summary statistics table with columns:

- `scenario`: SSP scenario name
- `peak_temperature`: Maximum temperature anomaly (°C)
- `year_peak_temp`: Year of peak temperature
- `years_above_1p5`: Years exceeding 1.5°C
- `years_above_2`: Years exceeding 2°C
- `total_mitigation`: Cumulative mitigation (GtCO2)
- `total_cdr`: Cumulative CDR (GtCO2)
- `total_mitigation_cost`: Sum of mitigation costs (trillion USD)
- `total_cdr_cost`: Sum of CDR costs (trillion USD)
- `total_damage_cost`: Sum of damage costs (trillion USD)
- `total_cost`: Sum of all costs (trillion USD)
- `year_first_1p5C`: Year first reaching 1.5°C (NA if never reached)

**Use for:**
- Quick comparison across scenarios
- Tables in publications
- Further analysis in other software (Excel, Python, etc.)

#### PDF Figure (`scenario_dashboard_TIMESTAMP.pdf`)

Multi-panel visualisation showing:

1. **Temperature trajectories:** Temperature anomaly over time for all scenarios
2. **Emissions:** Baseline vs controlled emissions
3. **Mitigation deployment:** Mitigation efforts over time
4. **CDR deployment:** CDR efforts over time
5. **Cost breakdown:** Mitigation, CDR, and damage costs

**Customisation:** See Section 8.6 (Visualisation Styling)

---

### 9.3 Parameter Importance Outputs

#### RDS File (`parameter_importance_TIMESTAMP.rds`)

Contains:
```r
results <- readRDS("output/parameter_importance_20260108_153045.rds")

results$importance_summary        # Parameter rankings by outcome
results$correlation_matrix        # Full correlation matrix
results$sensitivity_indices       # Quantitative sensitivity measures
results$all_results              # Full results for all parameter combinations
```

**Importance summary structure:**

Ranks parameters by their influence on each outcome variable:

- `parameter`: Parameter name
- `peak_temperature_corr`: Correlation with peak temperature
- `total_cost_corr`: Correlation with total cost
- `importance_score`: Composite importance metric (0-1)

**Interpreting correlations:**
- Values near +1: Strong positive relationship
- Values near -1: Strong negative relationship
- Values near 0: Weak/no relationship

#### Figures

- `parameter_importance_dashboard_TIMESTAMP.pdf`: Multi-panel sensitivity analysis
- `tornado_diagram_TIMESTAMP.pdf`: Parameters ranked by influence
- `scatter_matrix_TIMESTAMP.pdf`: Pairwise relationships

---

### 9.4 Delayed Deployment Outputs

#### RDS File (`delayed_deployment_TIMESTAMP.rds`)

Contains:
```r
results <- readRDS("output/delayed_deployment_20260108_163015.rds")

results$results                  # Full results grid (all delay combinations)
results$summary_by_scenario      # Aggregated metrics per scenario
results$feasibility_map          # Boolean: feasible combinations
results$delay_grid              # Grid of tested delays
```

**Results grid structure:**

For each scenario and delay combination:
```r
results$results[[scenario]][[delay_combo]]

# Contains:
$mitigation_delay_years      # Mitigation delay for this run
$cdr_delay_years            # CDR delay for this run
$peak_temperature           # Peak temperature (°C)
$total_cost                # Total cost (trillion USD)
$feasible                  # Boolean: target achievable?
# ... plus all time series data
```

#### CSV File (`delayed_deployment_summary_TIMESTAMP.csv`)

Tabular format of results grid:

- `scenario`: SSP scenario
- `mitigation_delay`: Mitigation delay (years)
- `cdr_delay`: CDR delay (years)
- `peak_temperature`: Peak temperature (°C)
- `years_above_1p5`: Years above 1.5°C
- `total_cost`: Total cost (trillion USD)
- `abatement_cost`: Mitigation + CDR costs (trillion USD)
- `temp_damage_cost`: Climate damages (trillion USD)
- `feasible`: TRUE/FALSE

**Use for:**
- Creating custom visualisations
- Statistical analysis of delay effects
- Identifying cost-optimal delay strategies

#### Heatmap Figures

- `temperature_delay_dashboard_TIMESTAMP.pdf`: Peak temperature heatmap
- `cost_delay_dashboard_TIMESTAMP.pdf`: Cost heatmaps with gradients
- `comprehensive_cost_delay_dashboard_TIMESTAMP.pdf`: Multi-panel cost comparison

**Reading heatmaps:**
- X-axis: Mitigation delay (years)
- Y-axis: CDR delay (years)
- Colour: Outcome magnitude
- White areas/X marks: Infeasible combinations (target not achievable)
- Contour lines: Connect points with equal values
- Arrows (if enabled): Direction of steepest increase

---

### 9.5 Loading and Analysing Saved Results

**Loading RDS files:**
```r
# Load results
results <- readRDS("output/scenario_comparison_20260108_143022.rds")

# Access specific scenario
ssp2_results <- results$scenario_results$`SSP2-Baseline`

# Extract time series
temperature <- ssp2_results$solution$temperature
years <- ssp2_results$solution$year

# Plot custom visualisation
plot(years, temperature, type = "l", 
     xlab = "Year", ylab = "Temperature Anomaly (°C)",
     main = "SSP2 Temperature Trajectory")
```

**Loading CSV summaries:**
```r
# Load summary table
summary <- read.csv("output/scenario_comparison_summary_20260108_143022.csv")

# Calculate statistics
mean_peak_temp <- mean(summary$peak_temperature)
range_total_cost <- range(summary$total_cost)

# Create custom table
library(dplyr)
summary %>%
  select(scenario, peak_temperature, total_cost, years_above_1p5) %>%
  arrange(total_cost)
```

**Comparing multiple runs:**
```r
# Load different parameter sets or scenarios
run1 <- readRDS("output/scenario_comparison_20260108_143022.rds")
run2 <- readRDS("output/scenario_comparison_20260108_153045.rds")

# Compare key metrics
comparison <- data.frame(
  run = c("Run 1", "Run 2"),
  mean_peak_temp = c(
    mean(run1$comparison_summary$peak_temperature),
    mean(run2$comparison_summary$peak_temperature)
  ),
  mean_total_cost = c(
    mean(run1$comparison_summary$total_cost),
    mean(run2$comparison_summary$total_cost)
  )
)
print(comparison)
```

---

## Troubleshooting

This section addresses common issues and their solutions.

### 10.1 Installation and Setup Issues

#### Problem: Package installation fails

**Error message:** `package 'X' is not available for this version of R`

**Solutions:**
1. Update R to the latest version
2. Check package name spelling
3. Try installing from a different CRAN mirror:
```r
   install.packages("package_name", repos = "https://cloud.r-project.org/")
```
4. For specific packages, check if they require additional system dependencies

#### Problem: "Cannot find function X"

**Error message:** `Error: could not find function "function_name"`

**Solutions:**
1. Ensure Libraries chunk has been run successfully
2. Check that the package containing the function is loaded
3. Restart R session and re-run Libraries chunk
4. Verify package is actually installed: `"package_name" %in% installed.packages()`

#### Problem: Data files not found

**Error message:** `Error: 'emissions.csv' does not exist`

**Solutions:**
1. Verify files are in the `data/` directory
2. Check file names exactly match (case-sensitive on Linux/Mac)
3. Ensure working directory is set to project root:
```r
   getwd()  # Check current directory
   setwd("path/to/MACROM")  # Set if incorrect
```
4. Use RStudio Projects (.Rproj file) to automatically set working directory

---

### 10.2 Runtime and Performance Issues

#### Problem: Analysis runs very slowly

**Symptoms:** Chunks take much longer than expected runtimes

**Solutions:**
1. **Enable parallel processing:**
```r
   use_parallel = TRUE
```

2. **Check CPU usage:**
   - Open Task Manager (Windows) or Activity Monitor (Mac)
   - Verify R is using multiple cores if parallel enabled

3. **Reduce problem size for testing:**
```r
   # Test with single scenario first
   scenarios_to_compare <- c("SSP2-Baseline")
   
   # Reduce LHS samples
   lhs_parameter_df <- generate_lhs_samples(n_samples = 100)  # Instead of 500
   
   # Increase delay step size
   delay_step_size <- 10  # Instead of 5 or 1
```

4. **Close other applications:** Free up RAM and CPU

5. **Check for background processes:** Antivirus, backup software, etc.

#### Problem: Parallel processing fails

**Error message:** `Error in makePSOCKcluster... unable to create cluster`

**Solutions:**
1. **Disable parallel processing:**
```r
   use_parallel = FALSE
```

2. **Reduce number of cores:**
```r
   # In the source code, modify:
   n_cores <- 2  # Instead of detectCores() - 1
```

3. **Check firewall settings:** Some firewalls block inter-process communication

4. **Restart R session:** Clear any hung parallel processes

#### Problem: Out of memory errors

**Error message:** `Error: cannot allocate vector of size X Gb`

**Solutions:**
1. **Reduce problem size:**
   - Fewer LHS samples
   - Fewer scenarios
   - Larger delay step size

2. **Close other applications:** Free up RAM

3. **Disable parallel processing:** Uses less RAM but slower

4. **Increase system RAM:** If working with very large parameter sets

5. **Process scenarios sequentially:**
```r
   # Run scenario comparison one at a time
   for (scenario in scenarios_to_compare) {
     result <- run_scenario_comparison(
       parameter_df = parameter_df[1, ],
       scenarios = scenario,
       ...
     )
     saveRDS(result, paste0("output/", scenario, ".rds"))
   }
```

---

### 10.3 Analysis-Specific Issues

#### Problem: Scenario comparison fails for some scenarios

**Error message:** `Scenario X failed: solver did not converge`

**Possible causes:**
1. Infeasible parameter combination (target unachievable)
2. Numerical instabilities in solver
3. Data issues for that specific scenario

**Solutions:**
1. **Check which scenarios failed:**
```r
   scenario_comparison_results$failed_scenarios
```

2. **Try adjusting solver tolerances** (requires modifying source code in `optimal_control_core_V3.R`)

3. **Relax temperature target or constraints:**
   - Increase `temp_target` from 1.5°C to 1.75°C or 2.0°C
   - Increase `cdr_cumulative_limit`

4. **Run failed scenarios individually with verbose mode:**
```r
   result <- run_scenario_comparison(
     parameter_df = parameter_df[1, ],
     scenarios = "SSP5-Baseline",  # The failed scenario
     verbose = TRUE
   )
```

#### Problem: No convergence in optimal control solver

**Error message:** `Maximum iterations reached without convergence`

**Solutions:**
1. **Check parameter values are reasonable:**
   - Costs not negative
   - Temperature target achievable
   - Technology constraints not too restrictive

2. **Increase maximum iterations** (requires modifying source code)

3. **Try different initial conditions** (requires modifying source code)

4. **Simplify problem:**
   - Remove delays
   - Use less restrictive technology constraints

#### Problem: Parameter importance analysis produces unexpected results

**Symptoms:** 
- All parameters show low importance
- Correlations don't match expectations
- High variance in results

**Solutions:**
1. **Check LHS sample size:**
   - Minimum 100 samples recommended
   - 500+ samples for robust results
```r
   nrow(parameter_df)  # Check sample size
```

2. **Verify parameter ranges:**
   - Check `parameter_details.yml`
   - Ensure ranges are appropriate (not too narrow/wide)

3. **Check for failed runs:**
```r
   # Some parameter combinations may be infeasible
   sum(is.na(parameter_importance_results$all_results$peak_temperature))
```

4. **Examine individual parameter relationships:**
```r
   # Plot parameter vs outcome
   plot(parameter_df$discount_rate, 
        parameter_importance_results$all_results$total_cost)
```

#### Problem: Delayed deployment shows all combinations infeasible

**Symptoms:** Entire heatmap is white/marked infeasible

**Causes:**
- Temperature target too stringent
- Technology constraints too limiting
- Baseline emissions too high

**Solutions:**
1. **Relax temperature target:**
```yaml
   # In parameter_details_fixed.yml
   temp_target:
     value: 2.0  # Instead of 1.5
```

2. **Increase technology capacities:**
```yaml
   cdr_cumulative_limit:
     value: 2000  # Instead of 1000
   
   mitigation_max_rate:
     value: 0.15  # Instead of 0.10
```

3. **Test with SSP1-Baseline:** Lower baseline emissions → more feasible

4. **Reduce maximum delay tested:**
```r
   max_delay_years <- 30  # Instead of 50
```

---

### 10.4 Visualisation Issues

#### Problem: Plots don't appear

**Symptoms:** Code runs without error but no plots displayed

**Solutions:**
1. **Check if plots were saved to file:**
```r
   list.files("figs/", pattern = "*.pdf")
```

2. **Open PDF files directly** from `figs/` directory

3. **Try displaying explicitly:**
```r
   print(dashboard)
```

4. **Check plotting device:**
```r
   dev.cur()  # Check if plot device is open
   dev.off()  # Close current device
```

#### Problem: Figures have wrong dimensions

**Symptoms:** Text too small, plots cut off, poor aspect ratio

**Solutions:**
1. **Adjust width and height:**
```r
   width = 297   # A4 landscape
   height = 210
   
   # Or for taller multi-panel plots:
   height = 260
```

2. **Adjust text size:**
```r
   text_size = 11  # Increase if too small
```

3. **Change DPI for different use cases:**
```r
   dpi = 150  # Screen viewing (larger text relative to plot)
   dpi = 300  # Publication (default)
```

#### Problem: Colours not displaying correctly

**Solutions:**
1. **Try different colour palette:**
```r
   color_palette = "viridis"  # Default
   color_palette = "plasma"   # Alternative
```

2. **Check for colourblind-friendly palettes:**
   - All default palettes are colourblind-safe
   - Avoid custom palettes unless tested

3. **Increase colour contrast:**
   - Use diverging palettes for variables that have meaningful zero-point
   - Use sequential palettes for variables that are always positive

#### Problem: PDF files won't open

**Error message:** File appears corrupted or fails to open

**Solutions:**
1. **Check file was completely written:**
```r
   file.info("figs/scenario_dashboard_20260108_143045.pdf")$size
```
   - Size should be > 0 bytes

2. **Ensure Cairo package is installed and working:**
```r
   library(Cairo)
   capabilities("cairo")  # Should be TRUE
```

3. **Try alternative graphics device:**
```r
   # Modify save_plot code to use:
   pdf(file = "output.pdf")
   print(plot)
   dev.off()
```

---

### 10.5 Data and Results Issues

#### Problem: Results look unrealistic

**Symptoms:**
- Negative costs
- Temperature decreases instantly
- Extremely high mitigation/CDR values

**Solutions:**
1. **Check parameter values:**
```r
   print(parameter_df)
   # Look for unreasonable values
```

2. **Verify data loaded correctly:**
```r
   head(emissions_df)
   summary(emissions_df)
   # Check for NAs, negative values, unrealistic magnitudes
```

3. **Check units:**
   - Emissions: GtCO₂/year
   - Temperature: °C above pre-industrial
   - Costs: Trillion USD
   - GWP: Trillion USD

4. **Examine intermediate results:**
```r
   # Plot baseline emissions
   library(ggplot2)
   ggplot(emissions_df, aes(x = Year, y = Emissions, color = Scenario)) +
     geom_line()
```

#### Problem: Results vary between runs with same parameters

**Symptoms:** Running same analysis twice gives different results

**Possible causes:**
1. Different random seeds (LHS sampling)
2. Numerical instabilities in solver
3. Parallel processing order effects

**Solutions:**
1. **Set random seed for LHS:**
```r
   lhs_parameter_df <- generate_lhs_samples(n_samples = 500, seed = 12345)
```

2. **Use fixed parameters for reproducibility:**
   - Switch to `fixed_model_parameters_call`

3. **Disable parallel processing if order matters:**
```r
   use_parallel = FALSE
```

4. **Save and compare results:**
```r
   result1 <- readRDS("output/run1.rds")
   result2 <- readRDS("output/run2.rds")
   all.equal(result1, result2)
```

---

### 10.6 Getting Help

If you encounter an issue not covered here:

1. **Check error message carefully:** Often contains specific information about the problem

2. **Search GitHub issues:** https://github.com/nina-rynne/MACROM/issues
   - Someone may have encountered the same problem

3. **Create a minimal reproducible example:**
   - Simplify to smallest case that reproduces error
   - Single scenario, fixed parameters, minimal customisation

4. **Report the issue:**
   - Email: nina.rynne@griffithuni.edu.au
   - GitHub: Create new issue with:
     - Error message (full text)
     - Code used (minimal example)
     - R version and operating system
     - Package versions: `sessionInfo()`

5. **Include diagnostic information:**
```r
   # Run this and include output in bug report
   sessionInfo()
```

---

## Advanced Usage

This section covers advanced techniques for experienced users.

### 11.1 Modifying Source Code

MACROM's functionality is contained in R scripts in the `src/` directory. Advanced users can modify these to:

- Change solver algorithms
- Add new cost functions
- Implement different climate models
- Add new analysis types

**Key source files:**

- `optimal_control_core_V3.R`: Core optimal control solver
- `model_parameters.R`: Parameter definitions and defaults
- `data_preparation.R`: Data loading and preprocessing
- `*_visualisation.R`: Plotting functions

**Best practices:**

1. **Make a backup:**
```bash
   cp src/optimal_control_core_V3.R src/optimal_control_core_V3_backup.R
```

2. **Document changes:**
   - Add comments explaining modifications
   - Update function documentation

3. **Test thoroughly:**
   - Run on simple test cases first
   - Verify results are reasonable
   - Compare to original version

4. **Version control:**
   - Use git to track changes
   - Create feature branches for major modifications

### 11.2 Custom Parameter Ranges

To define custom parameter ranges for LHS sampling:

1. **Edit `parameter_details.yml`:**
```yaml
   discount_rate:
     min: 0.01
     max: 0.05
     distribution: uniform
     description: "Annual discount rate"
   
   temp_target:
     min: 1.5
     max: 2.0
     distribution: uniform
     description: "Target temperature (°C)"
```

2. **Add new parameters:**
```yaml
   my_new_parameter:
     min: 0
     max: 100
     distribution: uniform
     description: "Description of new parameter"
```

3. **Modify sampling code if using non-uniform distributions** (requires editing `latin_hypercube_sampling.R`)

### 11.3 Adding New SSP Scenarios

To analyse custom emission scenarios:

1. **Prepare data file:**
   - Format: Same as existing `emissions.csv` and `gwp.csv`
   - Columns: `Year`, `Scenario`, `Emissions` (or `GWP`)
   - Scenario name: Use unique identifier (e.g., "Custom-HighEmissions")

2. **Add to existing data:**
```r
   # Load custom data
   custom_emissions <- read.csv("data/custom_emissions.csv")
   
   # Combine with SSP data
   emissions_df <- rbind(emissions_df, custom_emissions)
```

3. **Run analysis:**
```r
   scenarios_to_compare <- c("SSP2-Baseline", "Custom-HighEmissions")
```

### 11.4 Batch Processing Multiple Configurations

For running many different configurations (e.g., testing multiple temperature targets):
```r
# Define configurations
configs <- list(
  list(temp_target = 1.5, cdr_limit = 1000, name = "config1"),
  list(temp_target = 1.75, cdr_limit = 1000, name = "config2"),
  list(temp_target = 2.0, cdr_limit = 1000, name = "config3"),
  list(temp_target = 1.5, cdr_limit = 1500, name = "config4")
)

# Run all configurations
all_results <- list()
for (config in configs) {
  # Update parameters
  parameter_df$temp_target <- config$temp_target
  parameter_df$cdr_cumulative_limit <- config$cdr_limit
  
  # Run analysis
  results <- run_scenario_comparison(
    parameter_df = parameter_df[1, ],
    emissions_df = emissions_df,
    economic_df = economic_df,
    scenarios = c("SSP2-Baseline"),
    save_results = TRUE,
    output_prefix = paste0("scenario_comparison_", config$name)
  )
  
  all_results[[config$name]] <- results
}

# Compare results
comparison <- data.frame(
  config = names(all_results),
  temp_target = sapply(configs, function(x) x$temp_target),
  cdr_limit = sapply(configs, function(x) x$cdr_limit),
  peak_temp = sapply(all_results, function(x) 
    x$comparison_summary$peak_temperature[1]),
  total_cost = sapply(all_results, function(x) 
    x$comparison_summary$total_cost[1])
)
print(comparison)
```

### 11.5 Custom Visualisations

Creating custom plots from saved results:
```r
library(ggplot2)
library(dplyr)

# Load results
results <- readRDS("output/scenario_comparison_20260108_143022.rds")

# Extract data for custom plot
plot_data <- data.frame()
for (scenario in names(results$scenario_results)) {
  scenario_data <- results$scenario_results[[scenario]]$solution
  scenario_data$Scenario <- scenario
  plot_data <- rbind(plot_data, scenario_data)
}

# Custom temperature trajectory plot
ggplot(plot_data, aes(x = year, y = temperature, color = Scenario)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "orange") +
  labs(
    title = "Temperature Trajectories Under Optimal Control",
    x = "Year",
    y = "Temperature Anomaly (°C above pre-industrial)",
    color = "SSP Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Custom cost comparison
cost_summary <- results$comparison_summary %>%
  select(scenario, total_mitigation_cost, total_cdr_cost, total_damage_cost) %>%
  tidyr::pivot_longer(
    cols = -scenario,
    names_to = "cost_type",
    values_to = "cost"
  )

ggplot(cost_summary, aes(x = scenario, y = cost, fill = cost_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Cost Breakdown by Scenario",
    x = "SSP Scenario",
    y = "Cost (Trillion USD)",
    fill = "Cost Component"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### 11.6 Exporting Results for Other Software

**Export to Excel:**
```r
library(writexl)

# Export scenario comparison summary
write_xlsx(
  list(
    Summary = scenario_comparison_results$comparison_summary,
    Metadata = data.frame(
      Parameter = names(parameter_df),
      Value = as.character(t(parameter_df))
    )
  ),
  path = "output/scenario_comparison_export.xlsx"
)
```

**Export time series data:**
```r
# Extract all time series
for (scenario in names(scenario_comparison_results$scenario_results)) {
  scenario_data <- scenario_comparison_results$scenario_results[[scenario]]$solution
  
  write.csv(
    scenario_data,
    file = paste0("output/timeseries_", scenario, ".csv"),
    row.names = FALSE
  )
}
```

**Export for Python analysis:**
```r
# Save as CSV for pandas
write.csv(
  scenario_comparison_results$comparison_summary,
  "output/for_python_analysis.csv",
  row.names = FALSE
)
```

### 11.7 Integration with Other Models

MACROM can be integrated with other climate or economic models:

**Using MACROM outputs as inputs:**
```r
# Extract optimal mitigation trajectory
optimal_mitigation <- scenario_comparison_results$scenario_results$`SSP2-Baseline`$solution

# Export for use in IAM or Earth System Model
write.csv(
  optimal_mitigation[, c("year", "mitigation", "cdr", "emissions")],
  "output/optimal_trajectory_for_external_model.csv",
  row.names = FALSE
)
```

**Using external model outputs in MACROM:**
```r
# Import custom climate response function
source("path/to/external_climate_model.R")

# Modify optimal_control_core_V3.R to call external function
# (requires editing source code)
```

### 11.8 High-Performance Computing

For very large parameter sets or fine-grained analyses:

**Running on HPC cluster:**
```r
# Save workflow as script
# Create job submission script for your HPC system

# Example PBS script:
# #!/bin/bash
# #PBS -l nodes=1:ppn=16
# #PBS -l walltime=24:00:00
# 
# cd $PBS_O_WORKDIR
# module load R/4.2.0
# Rscript MACROM_batch_run.R
```

**Distributed computing with future package:**
```r
library(future)
library(future.apply)

# Set up distributed backend
plan(multisession, workers = 8)

# Run scenarios in parallel using future
results <- future_lapply(scenarios_to_compare, function(scenario) {
  run_scenario_comparison(
    parameter_df = parameter_df[1, ],
    scenarios = scenario,
    ...
  )
})
```

---

## Frequently Asked Questions

### General Questions

#### Q: What does MACROM stand for?

**A:** MACROM stands for "Mitigation and Carbon Removal Optimisation Model" (or sometimes "Model for Abatement and Climate Response Optimisation Management").

#### Q: What is optimal control and why use it for climate policy?

**A:** Optimal control is a mathematical framework for finding the best way to control a dynamic system over time. For climate policy, it helps identify the most cost-effective path for deploying mitigation and CDR to achieve temperature targets. Unlike simple scenarios or rules-of-thumb, optimal control considers:
- Trade-offs between acting now vs later
- Interactions between mitigation and CDR
- Economic costs vs climate damages
- Physical and technological constraints

#### Q: How long does it take to learn and use MACROM?

**A:** 
- **Basic usage** (running with defaults): 30 minutes - 1 hour
- **Understanding outputs**: 1-2 hours
- **Custom analyses**: 2-5 hours
- **Modifying source code**: 5-10 hours (requires more R experience)

#### Q: Do I need to know optimal control theory to use MACROM?

**A:** No. The workflow handles all the mathematics internally. However, understanding the basic concept helps interpret results:
- The model finds the strategy that minimises total costs (abatement + damages)
- It balances "act now" (higher abatement costs) vs "act later" (higher damages)
- Results represent economically optimal strategies, not necessarily politically feasible ones

---

### Technical Questions

#### Q: What climate model does MACROM use?

**A:** MACROM uses a simplified carbon cycle and climate response model. It includes:
- Cumulative emissions to temperature relationship
- Multi-box carbon cycle (atmosphere, surface ocean, deep ocean)
- Temperature response with thermal inertia
- These simplifications allow fast computation while capturing key dynamics

For detailed equations, see the paper or `optimal_control_core_V3.R`.

#### Q: What cost functions are included?

**A:** MACROM includes three cost components:
1. **Mitigation costs**: Based on marginal abatement cost curves
2. **CDR costs**: Technology-specific removal costs
3. **Climate damages**: Temperature-dependent economic impacts

Cost parameters can be customised in `parameter_details_fixed.yml`.

#### Q: Can MACROM handle other greenhouse gases besides CO₂?

**A:** The current version focuses on CO₂. However, other greenhouse gases can be approximated through:
- CO₂-equivalent emissions in input data
- Adjusted temperature response parameters
- Modified cost functions

Contact the developers if you need multi-gas capability.

#### Q: What's the difference between mitigation and CDR in MACROM?

**A:**
- **Mitigation**: Reduces emissions at source (prevents CO₂ from entering atmosphere)
  - Examples: Renewable energy, efficiency improvements, fuel switching
  - Limited by baseline emission levels
  - Generally cheaper per tonne than CDR
  
- **CDR** (Carbon Dioxide Removal): Removes CO₂ already in the atmosphere
  - Examples: Direct air capture, afforestation, enhanced weathering
  - Not limited by baseline emissions
  - Can achieve net-negative emissions
  - Subject to cumulative capacity constraints

#### Q: How does MACROM handle uncertainty?

**A:** Through Latin Hypercube Sampling (LHS):
- Samples parameter space to explore uncertainty
- Generates distributions of outcomes
- Identifies most influential parameters
- Provides confidence intervals on results

See "Parameter Importance Analysis" section for details.

#### Q: What solver does MACROM use?

**A:** MACROM uses a shooting method to solve the optimal control problem:
- Converts boundary value problem to root-finding
- Uses numerical optimization to find optimal controls
- Fast and robust for this class of problems

---

### Data and Scenarios

#### Q: What are SSP scenarios?

**A:** Shared Socioeconomic Pathways (SSPs) represent different futures of socioeconomic development:

- **SSP1 (Sustainability)**: Low challenges to mitigation and adaptation
- **SSP2 (Middle-of-road)**: Moderate challenges  
- **SSP3 (Regional rivalry)**: High challenges to mitigation and adaptation
- **SSP4 (Inequality)**: Low challenges to mitigation, high to adaptation
- **SSP5 (Fossil-fueled)**: High challenges to mitigation, low to adaptation

MACROM uses "baseline" scenarios (no climate policy) as starting points.

#### Q: Where do the emission and economic data come from?

**A:** The SSP data comes from:
- Emissions: SSP Database (https://tntcat.iiasa.ac.at/SspDb/)
- Economics: SSP GDP projections from integrated assessment models
- Data is preprocessed and included in the repository

#### Q: Can I use my own emission scenarios?

**A:** Yes! See Section 11.3 (Adding New SSP Scenarios) for instructions. Your data must include:
- Annual emissions (GtCO₂/year)
- Gross World Product (trillion USD)
- Years matching your analysis period

#### Q: What time period can I analyse?

**A:** The default is 2020-2100, but you can customise:
- Start year: 2020 or later (limited by data availability)
- End year: Up to 2100 (or later if you extend input data)
- Time step: Annual (1 year) is default, but sub-annual is possible

See Section 8.1 (Time Range customisation).

---

### Results and Interpretation

#### Q: Why do some scenarios show higher costs than others?

**A:** Higher costs typically result from:
- Higher baseline emissions (more abatement needed)
- Tighter temperature targets (requires more aggressive action)
- Delayed deployment (must act faster later, less efficient)
- More restrictive technology constraints (fewer options available)

The scenario comparison analysis helps disentangle these factors.

#### Q: What does "infeasible" mean in delayed deployment results?

**A:** Infeasible means the temperature target cannot be achieved with that delay combination, even with maximum mitigation and CDR deployment. This indicates:
- Delays are too long
- Temperature target is too stringent
- Technology constraints are too limiting
- Need to either act sooner or relax the target

#### Q: Why does temperature overshoot the target temporarily?

**A:** Temperature overshoot occurs because:
- Climate system has inertia (temperatures continue rising after emissions stop)
- CDR takes time to deploy at scale
- Economic optimization may favor some overshoot if damages are manageable
- Need "net-negative" emissions to bring temperature back down

This is a key policy question MACROM helps address.

#### Q: How should I interpret parameter importance results?

**A:** Parameters ranked "high importance" have strong influence on outcomes:
- Prioritise these for careful calibration
- Uncertainty in these parameters translates to outcome uncertainty
- Focus sensitivity discussions on these parameters

Low importance parameters:
- Model is relatively insensitive to these
- Less critical for calibration
- Can use broader uncertainty ranges

#### Q: Are MACROM's cost estimates realistic?

**A:** MACROM's costs are:
- Order-of-magnitude estimates for comparison, not precise forecasts
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
2. Review outputs to understand model behavior
3. Then expand to all scenarios or add other analyses

This approach:
- Gives quick results (~10 minutes)
- Helps you understand outputs before tackling longer analyses
- Reveals any data or setup issues early

#### Q: Can I run multiple analyses on the same parameter set?

**A:** Yes! Once you've run the parameter setup chunks, you can run all analysis types:
```r
# Run fixed parameters ONCE
parameter_df <- create_params_dataframe()

# Then run multiple analyses
run_scenario_comparison(...)
run_delayed_deployment(...)
# etc.
```

Just don't re-run the parameter chunks in between.

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

#### Q: How do I change the temperature target?

**A:** Edit `parameter_details_fixed.yml`:
```yaml
temp_target:
  value: 2.0  # Change from 1.5 to 2.0°C
```

Then re-run `fixed_model_parameters_call` chunk.

#### Q: Can I model specific CDR technologies (e.g., DACCS vs afforestation)?

**A:** The current version treats CDR as homogeneous. To model specific technologies:
- Modify cost functions in source code
- Add separate state variables for each technology
- Adjust constraints for technology-specific limits

This requires advanced modifications (see Section 11.1).

#### Q: How do I add a new cost function?

**A:** This requires modifying `optimal_control_core_V3.R`:
1. Define new cost function
2. Add to objective function in solver
3. Add parameters to `parameter_details_fixed.yml`
4. Test thoroughly on simple cases

Contact developers if you need guidance.

#### Q: Can MACROM handle regional analysis?

**A:** The current version is global. Regional analysis would require:
- Regional emission and economic data
- Regional climate response (more complex)
- Inter-regional trade and cooperation mechanisms
- Significant model extensions

This is beyond current scope but could be future development.

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

**A:** Yes! Planned developments include:
- Additional cost function options
- Multi-gas capability
- Improved solver algorithms
- More visualisation options
- Tutorial materials

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

