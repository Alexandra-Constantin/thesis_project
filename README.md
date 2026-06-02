# Extending the AccelerAge Framework to High-Dimensional Settings Using Bayesian Gompertz AFT Models
Code accompanying my MSc thesis in Statistics and Data Science at Leiden University.

This project extends the AccelerAge framework ([Sluiskes et al., 2024](https://doi.org/10.1007/s10654-024-01114-8)) to high-dimensional predictor settings using Bayesian Gompertz accelerated failure time (AFT) models. The work focuses on estimating biological age from metabolomic biomarkers while accounting for left-truncated and right-censored survival data.

## Project Background
Chronological age is a strong predictor of mortality, but individuals of the same age can differ substantially in their health status and ageing trajectories. Biological age models aim to capture these differences by combining biomarker information into a single measure that better reflects an individual's physiological state.

This thesis builds upon the AccelerAge framework, where biological age is defined through the mean residual life (mrl) function and with respect to some reference population.To enable the use of high-dimensional biomarker data, the original framework was extended using Bayesian Gompertz AFT models implemented in Stan via the brms package.

## Repository Structure
This repository contains the code used for the simulation study presented in the thesis.

### Scripts

**01_dgp_helper_functions.R**

Contains the functions used to generate the simulated datasets and helper functions. This includes the generation of correlated predictors, simulation of survival times under the Gompertz AFT model, construction of the reference life table, and the functions used to obtain mean residual life and biological age estimates.

**02_stan_definitions_priors.R**

Contains the Stan model definitions and prior specifications used in the Bayesian analyses, including both the regularised horseshoe prior and the Gaussian prior.

**03_fit_models_parallel.R**

Fits the Bayesian and frequentist Gompertz AFT models. The simulations are run in parallel to reduce computation time.

**04_run_simulation.R**

Main script used to run the simulation study. It generates the datasets, fits the models, calculates the performance measures, and saves the results.

**05_mcmc_diagnostics.R**

Used to assess the convergence of individual Bayesian model fits. The script extracts and summarises diagnostics such as R-hat, effective sample sizes, divergences, and E-BFMI.

**06_mcmc_diagnostics_grid.R**

Summarises convergence diagnostics across the entire simulation grid; used to identify scenarios where sampling difficulties occurred.

**07_results_plots.R**

Produces the figures used in the Results section of the thesis.

### Workflow
The simulation study was run in the following order:
1. Define the data-generating mechanism and helper functions (`01`).
2. Specify the Bayesian models and priors (`02`).
3. Fit the models (`03`).
4. Run the full simulation study across all scenarios (`04`).
5. Evaluate convergence diagnostics (`05` and `06`).
