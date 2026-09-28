# MATLAB simulations — structure aligned with Section 5 of the paper

This folder reorganizes the publication MATLAB code so that the directory
structure follows the numerical-results section of the manuscript.

## Quick start

From the repository root in MATLAB:

```matlab
setup_paths
```

Then enter the relevant subsection folder and run the `run_*.m` script.
Where a separate plotting script is supplied, run it after the result file
has been generated (or use the included saved result file where available).

## Paper-to-code map

### 5.1 Exact likelihood versus Gaussian approximation

**5.1a Population-size sensitivity**
- `5_1_exact_vs_gaussian/5_1a_population_size/run_reference_comparison.m`
- `.../plot_reference_comparison.m`
- `.../plot_QP_heatmaps.m`
- The two `finalize_*.m` scripts are retained because they belong to the
  publication workflow for the reference Monte Carlo run.

**5.1b Observation-length sensitivity**
- `5_1_exact_vs_gaussian/5_1b_observation_length/run_T_sensitivity.m`
- `.../plot_T_sensitivity.m`
- `.../experiment1_T_sensitivity_results.mat`

**5.1c Observation-interval sensitivity**
- `5_1_exact_vs_gaussian/5_1c_observation_interval/run_delta_sensitivity.m`

**5.1d Robustness across heterogeneous transition mechanisms**
- `5_1_exact_vs_gaussian/5_1d_heterogeneous_mechanisms/run_parameter_robustness.m`
- `.../plot_parameter_robustness.m`

### 5.2 Recovery from equilibrium fluctuations
- `5_2_equilibrium_fluctuations/run_experiment2.m`
- `5_2_equilibrium_fluctuations/plot_experiment2.m`
- saved result file included.

### 5.3 Near-equilibrium information transition
- `5_3_near_equilibrium_information/run_experiment3.m`
- This script performs the deterministic Fisher-information calculation and
  produces its own outputs.

### 5.4 Scalability in large histogram state spaces
- `5_4_scalability/run_experiment4_population.m`
- `5_4_scalability/run_experiment4_temporal.m`
- `5_4_scalability/plot_experiment4.m`
- the available final temporal result file is included.

### 5.5 Robustness to transition topology
- `5_5_topology_robustness/run_topology_robustness.m`
- `5_5_topology_robustness/plot_topology_robustness.m`

The publication version of the topology script is set to `R=50`, matching
the manuscript. The uploaded source contained two consecutive `RUN_MODE`
assignments and used `R=100` in full mode; those two inconsistencies were
cleaned up in this consolidated copy.

## Shared functions

`common/` contains the reusable likelihood, simulation, parameterization,
and dynamic-programming functions used by the publication experiments.

## Important reproducibility note

This consolidation reorganizes the code but does not fabricate missing result
files. The uploaded repository itself notes that some final Monte Carlo result
files were not present in the package. Those simulations therefore still need
to be rerun (or their final saved `.mat` files added) for a completely
self-contained archival repository.

## Recommended execution order

1. `setup_paths`
2. Section 5.1 studies
3. Section 5.2
4. Section 5.3
5. Section 5.4
6. Section 5.5

Long Monte Carlo jobs should be run from their own subsection directory so
their checkpoint and result files remain beside the corresponding code.
