# META-AC

This is a version of the _Model for Economic Tipping point Analysis_ ([META](https://github.com/openmodels/META)) that includes the **AMOC Carbon (AC) feedback**. For the scripts of _Schaumann & Alastrué de Asenjo (2024): Weakening AMOC reduces carbon drawdown and increases the social cost of carbon_, see [here](https://github.com/felixschaumann/AMOC-Carbon).

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
5. [`src/components/AMOC.jl`](#5-srccomponentsamocjl) (modified)
6. [`src/montecarlo.jl`](#6-srcmontecarlojl) (modified)
7. [`src/scc.jl`](#7-srcsccjl) (modified)
 
#### 1. `src/components/AMOC-Carbon.jl`

This new component 

#### 2. `data/CMIP6_amoc`

This folder contains AMOC strength projection data for the following Earth System Models:

- MIP-ESM
- ACCESS-ESM
- CanESM
- GISS
- MIROC
- CESM2
- NorESM



#### 3. `src/components/CO2Converter_AMOC.jl`

This component is a copy of `src/components/CO2Converter.jl`, with the only difference being an additional variable for `co2_amoc[tt]`, which is added to annual CO$_2$ emissions `E_co2[tt]`.

#### 4. `src/lib/AMOC_Carbon_functions.jl`

#### 5. `src/components/AMOC.jl`

#### 6. `src/montecarlo.jl`

- sample_id_subset

#### 7. `src/scc.jl`

- sample_id_subset