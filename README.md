# MACROM: An Optimal Control Model for Balancing Climate Change Abatement and Damage Trade-offs

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)](https://www.r-project.org/)

**Authors:** Nina Rynne, Michael Bode, Melanie Roberts, Ryan Heneghan  
**Institution:** Griffith University  
**Contact:** nina.rynne@griffithuni.edu.au

---

## Overview

MACROM (Mitigation and Carbon Removal Optimisation Model) is a simple, emissions-driven, cost-minimisation climate-economic model that projects cost-optimal deployment of mitigation and carbon dioxide removal (CDR) to achieve specific temperature targets. 

MACROM captures:
1. Trade-offs between the costs of climate action and the economic damages of climate inaction
2. The optimal timing and scale of deploying climate action across the target's time horizon

We use MACROM to assess the mitigation and CDR required to return temperature to the Paris Agreement 1.5°C by 2100 across different socio-economic futures (SSP scenarios).

### Key Features

- **Optimal control framework**: Determines economically optimal mitigation and CDR trajectories
- **Multiple scenarios**: Compatible with all SSP baseline emission pathways
- **Flexible analysis**: Supports scenario comparison, parameter sensitivity, and delayed deployment analysis
- **Publication-ready outputs**: Automated generation of figures and dashboards
- **Open source**: Fully documented and reproducible

## Quick Start

### Prerequisites

- R (≥ 4.0.0)
- RStudio (recommended)
- Required R packages (see `MACROM_workflow.Rmd`)

### Installation

1. Clone the repository:
```bash
   git clone https://github.com/nina-rynne/MACROM.git
   cd MACROM
```

2. Ensure required data files are in the `data/` directory:
   - `emissions.csv`: SSP emission scenarios
   - `gwp.csv`: Gross World Product projections

3. Open `MACROM_workflow.Rmd` in RStudio

### Running the Model

1. Run the Libraries and Data Preparation chunks
2. Choose ONE parameter approach:
   - `fixed_model_parameters_call` for single parameter set (recommended for most users)
   - `latin_hypercube_sampling_call` for sensitivity analysis
3. Run the desired analysis chunks:
   - Scenario comparison
   - Parameter importance (requires LHS)
   - Delayed deployment analysis
4. View results in `output/` and `figs/` directories

**For detailed instructions, see [USER_GUIDE.md](USER_GUIDE.md)**

## Project Structure
```
MACROM/
├── MACROM_workflow.Rmd  # Main user interface
├── README.md                      # This file
├── USER_GUIDE.md                  # Detailed documentation
├── LICENSE                        # CC-BY-4.0 licence
├── data/                          # Processed data files
│   ├── emissions.csv
│   └── gwp.csv
├── data-raw/                      # Raw, unprocessed data
├── figs/                          # Generated figures and visualisations
├── output/                        # Analysis outputs (RDS, CSV)
└── src/                           # Source code (R functions)
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

## Model Details

### What MACROM Does

MACROM simulates the change in cumulative CO₂ emissions over time resulting from:
- **Anthropogenic activity**: Baseline emissions under different SSP scenarios
- **Mitigation**: Reducing emissions at source (preventing release)
- **Carbon Dioxide Removal (CDR)**: Actively removing CO₂ from the atmosphere

The model uses optimal control theory to determine the most cost-effective combination and timing of these interventions to achieve temperature targets while minimising total costs (abatement costs + climate damages).

### Key Model Components

1. **Climate dynamics**: Simple carbon cycle and temperature response model
2. **Economic framework**: Cost functions for mitigation, CDR, and temperature damages
3. **Optimal control**: Shooting method to solve for optimal trajectories
4. **Scenario analysis**: SSP baseline emissions and economic projections

### Analysis Capabilities

- **Scenario Comparison**: Compare optimal strategies across SSP1-5 baseline scenarios
- **Parameter Sensitivity**: Latin Hypercube Sampling for uncertainty quantification
- **Delayed Deployment**: Assess costs of delaying mitigation or CDR implementation
- **Temperature Overshoot**: Analyse peak warming and years above thresholds (1.5°C, 2°C)

## Citation

If you use MACROM in your research, please cite:
```bibtex
@article{rynne2025macrom,
  title={MACROM: An Optimal Control Model for Balancing Climate Change Abatement and Damage Trade-offs},
  author={Rynne, Nina and Bode, Michael and Roberts, Melanie and Heneghan, Ryan},
  journal={[Journal name]},
  year={2025},
  publisher={[Publisher]},
  note={[Add DOI when available]}
}
```

**Plain text citation:**
```
Rynne, N., Bode, M., Roberts, M., & Heneghan, R. (2025). 
MACROM: An Optimal Control Model for Balancing Climate Change 
Abatement and Damage Trade-offs. [Journal details to be added upon publication]
```

## Licence

This project is licensed under the Creative Commons Attribution 4.0 International Licence (CC-BY-4.0).

You are free to:
- **Share**: Copy and redistribute the material in any medium or format
- **Adapt**: Remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- **Attribution**: You must give appropriate credit, provide a link to the licence, and indicate if changes were made

See the [LICENSE](LICENSE) file for the full licence text.

Copyright © 2025 Nina Rynne

