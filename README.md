# Inferring Transition Dynamics from Histogram Time Series

MATLAB code accompanying the manuscript **“Inferring Transition Dynamics from Histogram Time Series: Identifiability at and Near Equilibrium.”**

## Overview

This repository contains the MATLAB code used for the numerical studies in Section 5 of the manuscript. The directory structure follows the organization of the numerical-results section of the paper.

The code includes implementations of the exact histogram-transition likelihood, its Gaussian approximation, simulation routines, parameter estimation, and the numerical experiments used to study identifiability and inference at and near equilibrium.

## Quick start

Clone or download the repository and open MATLAB in the repository root. Run:

```matlab
setup_paths
```

This adds the shared functions in `common/` to the MATLAB path. Then run the relevant experiment from its subsection directory.

Some Monte Carlo experiments can take substantial time. Saved result files are included for selected experiments so that the corresponding plotting scripts can be run without repeating those simulations.

## Numerical experiments

### 5.1 Exact likelihood versus Gaussian approximation

#### 5.1a Population-size sensitivity

Directory: `5_1_exact_vs_gaussian/5_1a_population_size/`

Main files:
- `run_reference_comparison.m` — Monte Carlo comparison of exact and Gaussian likelihood inference across population sizes.
- `plot_reference_comparison.m` — plots the population-size comparison.
- `plot_QP_heatmaps.m` — visualizes the generator and finite-time transition matrices.
- `finalize_reference_checkpoint.m` and `finalize_reference_results.m` — utilities used to finalize the reference Monte Carlo results.

#### 5.1b Observation-length sensitivity

Directory: `5_1_exact_vs_gaussian/5_1b_observation_length/`

Main files:
- `run_T_sensitivity.m` — studies inference as the number of observed transitions increases.
- `plot_T_sensitivity.m` — plots the results.
- `experiment1_T_sensitivity_results.mat` — saved results used by the plotting script.

#### 5.1c Observation-interval sensitivity

Directory: `5_1_exact_vs_gaussian/5_1c_observation_interval/`

Main file:
- `run_delta_sensitivity.m` — studies the effect of the observation interval on recovery of the transition dynamics.

The script can optionally use the Parallel Computing Toolbox when parallel execution is enabled.

#### 5.1d Heterogeneous transition mechanisms

Directory: `5_1_exact_vs_gaussian/5_1d_heterogeneous_mechanisms/`

Main files:
- `run_parameter_robustness.m` — compares performance across heterogeneous stationary distributions and transition-rate parameters.
- `plot_parameter_robustness.m` — plots the resulting error distributions.

### 5.2 Recovery from equilibrium fluctuations

Directory: `5_2_equilibrium_fluctuations/`

Main files:
- `run_experiment2.m` — evaluates recovery of transition dynamics from equilibrium conditional fluctuations.
- `plot_experiment2.m` — plots the recovery results.
- `experiment2_results.mat` — saved simulation results.

### 5.3 Near-equilibrium information transition

Directory: `5_3_near_equilibrium_information/`

Main file:
- `run_experiment3.m` — computes the Gaussian Fisher-information contributions near equilibrium and illustrates the transition between mean-based and fluctuation-based information regimes.

This is a deterministic information calculation rather than a Monte Carlo estimation experiment.

### 5.4 Scalability in large histogram state spaces

Directory: `5_4_scalability/`

Main files:
- `run_experiment4_population.m` — population-size/scalability calculations.
- `run_experiment4_temporal.m` — temporal-data scalability calculations.
- `plot_experiment4.m` — plots the scalability results.
- `experiment4_stage2_temporal_results_R50.mat` — saved temporal experiment results.

These experiments require MATLAB's **Optimization Toolbox** (`fminunc`).

### 5.5 Robustness to transition topology

Directory: `5_5_topology_robustness/`

Main files:
- `run_topology_robustness.m` — compares inference under the path and nonlocal-tree transition topologies used in the manuscript.
- `plot_topology_robustness.m` — plots the comparison.

This experiment requires MATLAB's **Optimization Toolbox** (`fminunc`).

## Shared functions

The `common/` directory contains reusable routines for:

- construction of continuous-time Markov generators and transition matrices;
- simulation of aggregate histogram transitions and time series;
- exact histogram-transition probabilities using dynamic programming;
- exact and Gaussian time-series likelihoods;
- parameter transformations and model fitting; and
- multinomial calculations and sampling.

The exact and Gaussian likelihood experiments in Section 5.1 use these shared routines.

## Software requirements

The code is written in MATLAB.

Most experiments use MATLAB functions available in the base installation. The scalability and topology experiments use `fminunc` and therefore require the **Optimization Toolbox**. The observation-interval experiment can optionally use the **Parallel Computing Toolbox** when parallel execution is enabled.

## Reproducibility

Random-number seeds are specified in the simulation scripts to support reproducible Monte Carlo experiments. Because some experiments involve repeated likelihood optimization and large Monte Carlo studies, runtimes can vary substantially with hardware and MATLAB configuration.

Saved `.mat` files are included where available. For experiments without saved final results, the corresponding `run_*.m` script should be executed before generating the associated plots.

## Repository structure

```text
.
├── setup_paths.m
├── common/
├── 5_1_exact_vs_gaussian/
│   ├── 5_1a_population_size/
│   ├── 5_1b_observation_length/
│   ├── 5_1c_observation_interval/
│   └── 5_1d_heterogeneous_mechanisms/
├── 5_2_equilibrium_fluctuations/
├── 5_3_near_equilibrium_information/
├── 5_4_scalability/
└── 5_5_topology_robustness/
```

## Citation

If you use this code in academic work, please cite the accompanying manuscript:

> I. Nevat, *Inferring Transition Dynamics from Histogram Time Series: Identifiability at and Near Equilibrium*.

The complete bibliographic citation can be added here once the article is published.

## License

No software license is currently specified in this repository. If the code is intended for reuse by others, an explicit license should be added to the repository.
