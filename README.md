# META-AC: AMOC-Carbon extension of the Model for Economic Tipping point Analysis

This is an extension of the _Model for Economic Tipping point Analysis_ ([META](https://github.com/openmodels/META)) that adds the **AMOC Carbon (AC) feedback**. For the scripts of _Schaumann & Alastrué de Asenjo (2025): Weakening AMOC reduces ocean carbon uptake and increases the social cost of carbon_, see [here for the paper-related code](https://github.com/felixschaumann/AMOC-Carbon). Find the Zenodo archive of this repository here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13830720.svg)](https://doi.org/10.5281/zenodo.13830720)

META is an advanced integrated assessment model (SC-IAM), designed as a model-based meta-analysis of the effects of tipping points on the social cost of carbon (SCC). The model simulates greenhouse gas emissions, temperature and sea-level rise, and market and non-market damages at the country level, and the effects of eight climate tipping points that have been studied in the climate economics literature.

This version of the model is implemented in [Mimi](https://www.mimiframework.org/), an integrated assessment modeling framework developed in [Julia](https://julialang.org/). 

> 📖 For scientific documentation of the core META model, see [here](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/openmodels/META/master/docs.pdf).

## Directories in the repository

The following directories are used for the Mimi model:

 - `data`: Input parameters and validation outputs (under
   `data/benchmark`).
 - `src`: The model code. The scripts directly contained in this directory
   support various types of analyses, with internal model code in
   subdirectories. Specifically, the model components are under
   `src/components` and additional functions are in `src/lib`.
 - `test`: Unit tests, for each component, for the system-wide
   results, and for the Monte Carlo system.
   
## Changes introduced in META-AC 
 
Compared to the core META model, META-AC introduces changes to:

1. [`src/components/AMOC-Carbon.jl`](#1-srccomponentsamoc-carbonjl) (new)
2. [`data/CMIP6_amoc`](#2-datacmip6_amoc) (new)
3. [`src/components/CO2Converter_AMOC.jl`](#3-srccomponentsco2converter_amocjl) (new)
4. [`src/lib/AMOC_Carbon_functions.jl`](#4-srclibamoc_carbon_functionsjl) (new)
5. [`test/test_AMOC_Carbon.jl`](#5-testtest_amoc_carbonjl)
6. [`src/components/AMOC.jl`](#6-srccomponentsamocjl) (modified)
7. [`src/montecarlo.jl`](#7-srcmontecarlojl) (modified)
8. [`src/scc.jl`](#8-srcsccjl) (modified)
 
#### 1. `src/components/AMOC-Carbon.jl`

This new component is the core of this extension. It models the AMOC not as a tipping element whose tipping is stochastically triggered and irreversible. Instead, we introduce explicit variables for `AMOC_strength[tt]` in Sv and `AMOC_decrease[tt]` since 2010 values (which are given by `AMOC_start`).
Further, there is a value for preindustrial AMOC strength `AMOC_pi`, and there are standard deviations `std_AMOC_strength[tt]` and `std_AMOC_pi`.

In every time step, `AMOC_decrease[tt]` is calculated based on the AMOC sensitivity to global temperature, which is parameterized by `beta_AMOC[tt]`. 
This `beta_AMOC[tt]` parameter can be calibrated such that the resulting `AMOC_strength[tt]` exactly corresponds to a given Earth System Model (ESM) projection.

In every time step, a flux of CO$_2$ to the atmosphere is calculated, based on the difference between `AMOC_pi` and `AMOC_strength[tt]`. In order to compare different AMOC projections, the absolute decrease values are scaled by the ratio of preindustrial AMOC strength of MPI-ESM and the ESM that is being calibrated to. The parameter that governs the sensitivity of CO$_2$ release to AMOC weakening is `eta_AMOC`.

The component can be calibrated either to a stylized scenario of AMOC weakening, or to an ESM projection. For a stylized weakening, only the amount of weakening between 2010 and 2100 has to be given to the function `addAMOC_Carbon`.

For calibrating to an ESM, the name of the model has to be given. Subsequently, the `addAMOC_Carbon` will calibrate `AMOC_pi`, `beta_AMOC[tt]` as well as all the uncertainties to the given ESM projection. 

The argument `temp_pattern` of `addAMOC_Carbon` allows to include the AMOC carbon feedback in a model that also features a change in pattern scaling induced by AMOC weakening (as in the standard AMOC component of the META model). If this option is selected, the parameter `pattern_onset` determines the year from which the background temperature pattern is scaled in. Otherwise, if `temp_pattern=false` in the `addAMOC_Carbon` call, the AMOC carbon feedback is included without considering other impacts of AMOC weakening.

#### 2. `data/CMIP6_amoc`

This directory contains AMOC strength projection data for the following ESMs along the SSP2-4.5 scenario:

- MPI-ESM1.2-LR
- ACCESS-ESM1-5
- CanESM5
- GISS-E2-1-G
- MIROC-ES2L
- CESM2
- NorESM2-LM

For each model, the directory also contains data about preindustrial AMOC strengths (`/data/CMIP6_amoc/AMOC_pi_values.csv`). For MPI-ESM1.2-LR and NorESM2-LM, it also contains projection data along SSP5-8.5 and SSP1-2.6, respectively.

#### 3. `src/components/CO2Converter_AMOC.jl`

This component is a copy of `src/components/CO2Converter.jl`, with the only difference being an additional variable for `co2_amoc[tt]`, which is added to annual CO$_2$ emissions `E_co2[tt]`.

#### 4. `src/lib/AMOC_Carbon_functions.jl`

`get_primary_model_name` ensures that the calling a model is not case-sensitive.

`get_projection_ds` loads the correct NCDataSet for a given ESM calibration.

`calibrate_AMOC` loads a given ESM projection NCDataSet and chooses `beta_AMOC[tt]` values such that the `AMOC_Carbon.jl` reproduces the AMOC strength of the given projection for the period 2016-2100. Before that, the AMOC is assumed to be at preindustrial values; after that, it is assumed to continue weakening proportional to the development of global mean temperature.

`add_AMOC` adds the original `AMOC.jl` component.

`get_amoc_carbon_model` builds and runs a `model` with the `AMOC_Carbon.jl` component included. The SSP scenario and the damage mode can be picked; the persistence parameter is fixed at `:damagepersist=0.25`. `temp_pattern` determines whether the AMOC carbon feedback is calculated against the backdrop of AMOC-induced changes in surface temperature patterns.

`get_base_model_w_params` is a function that builds a whole model with the correct set of parameters based on a parameter Dict that is generated by the `DrWatson` function `dict_list`.

`run_amoc_mc_sccs` takes the model from `get_base_model_w_params` and calculates an SCC value for each Monte Carlo sample. A full dictionary of results is returned.

`run_amoc_carbon_mc` is the equivalent function for making projections for each Monte Carlo samples, and saving the results in terms of `cum_CO2_AMOC`, damages, utility, and temperature.

`get_proj_std_ssp245` returns the standard deviation of a given ESM projection of AMOC strength.

#### 5. `test/test_AMOC_Carbon.jl`

This test ensures that the SCC is the same for these following cases:

- calling the `base_model` and calling the `full_model` with all tipping elements switched off
- calling the `base_model` and calling `get_amoc_carbon_model` with both background temperature pattern and the carbon effect switched off
- calling the standard AMOC temperature model through `full_model` and calling `get_amoc_carbon_model` with background temperature pattern switched on and the carbon effect switched off

This test is called by `test/runtests.jl`.

#### 6. `src/components/AMOC.jl`

The only change introduced here is to revert the calculation of `p_AMOC[tt]` back to how it was done in the [META-2021 version](https://github.com/openmodels/META-2021/blob/b9bacaf30403f9086ccea5632d7b96ad00d89de1/src/components/AMOC.jl#L22).

#### 7. `src/montecarlo.jl`

The first change is introduced to allow for repeatable Monte Carlo runs. For this, we comment out the updating of damage parameters within the function `setsim_base`. Additionally, we include the keyword argument `sample_id_subset` in the function `sim_base`, which is then passed on to `create_fair_monte_carlo`. This ensures that one can explicitly set the Monte Carlo samples that are to be used for the FaIR emulator.

The second change consists in passing Boolean keyword arguments to `sim_base` for setting whether to sample over uncertainty in AMOC strength, preindustrial AMOC strength and the sensitivity of carbon fluxes to AMOC weakening (`eta_AMOC`), respectively.
For each of those quantities, there is a keyword argument for an associated standard deviation, for `eta_AMOC`, also the mean value is passed.

#### 8. `src/scc.jl`

Analogous to the changes in `src/montecarlo.jl`, the `sample_id_subset` argument is passed through the relevant function calls. Similarly, the Monte Carlo parameters for uncertainty in AMOC strength, preindustrial AMOC strength and the sensitivity of carbon fluxes to AMOC weakening are passed through.

Additionally, we update the `emuc` parameter of the model's utility component with the `emuc` parameter that is used for SCC calculation in order to avoid a calculation where to different `emuc` values are used alongside, which leads to inconsistencies.
