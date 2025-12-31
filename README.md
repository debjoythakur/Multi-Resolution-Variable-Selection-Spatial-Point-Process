# Multi-Resolution-Variable-Selection-Spatial-Point-Process
This project proposes a multi-resolution variable selection method for spatial point process data, motivated by crime and traffic accident occurrences in St. Louis. Using Haar wavelets and penalized likelihood, it identifies which predictors matter, where, and at what spatial scale.

# Abstract
This project proposes a multi-resolution variable selection method for spatial point process data, motivated by crime and traffic accident occurrences in St. Louis. Using Haar wavelets and penalized likelihood, it identifies which predictors matter, where, and at what spatial scale.
Generally, Lasso, Adaptive Lasso, and SCAD are standard approaches in variable selection in the presence of a large number of predictors. In recent years, during intensity function estimation for spatial point processes with a diverging number of predictors, many researchers have considered these penalized methods. But we have discussed a multi-resolution perspective for the variable selection method for spatial point process data. Its advantage is twofold: it not only efficiently selects the predictors but also provides the idea of which points are liable for selecting a predictor at a specific resolution. Actually, our research is motivated by the crime and accident occurrences in  St. Louis and its neighborhoods. It is more relevant to select predictors at the local level, and thus we get the idea of which set of predictors is relevant for the occurrences of crime or accident in which parts of St. Louis. We describe the theoretical and simulation results to justify the accuracy of local-level variable selection during intensity function estimation.


## Overview
This repository implements a **localized spatial point process modeling framework**
for **road-level incident intensity** using **Berman–Turner likelihood approximation**
and **multi-resolution Haar wavelet regularization**.

The methodology integrates:
- Road-mapped accident and crime events
- Neighborhood-scale built-environment covariates
- Spatially varying coefficients via Haar wavelets
- Penalized likelihood estimation (LASSO, SCAD, Adaptive LASSO)
- Post-selection unpenalized inference

The analysis is conducted on **St. Louis, Missouri** road networks.

---

## Methodological Overview

The core components of the framework are:

- Poisson point process modeling of incident locations
- Berman–Turner quadrature approximation of the likelihood
- Spatially varying regression coefficients via 2D Haar wavelets
- Penalized likelihood estimation:
  - Localized LASSO
  - Localized SCAD
  - Global LASSO
  - Global SCAD
  - Global Adaptive LASSO
- Post-selection unpenalized BT refitting for stable inference

This approach allows us to identify:
- **Which covariates matter**
- **Where they matter**
- **At what spatial resolution do they matter**

---

## Repository Structure

├── README.md
├── code/
│ ├── simulation/ # simulation studies
│ └── real_data/ # real-data analysis scripts
├──data/ # Shape Files
└── Output/ # Results
---

## Required R Packages

Core dependencies:
- `sf`
- `spatstat`
- `glmnet`
- `ncvreg`
- `dplyr`
- `ggplot2`

Optional (for multi-panel plots):
- `patchwork`

Installation:
```r
install.packages(c(
  "sf","dplyr","ggplot2",
  "glmnet","ncvreg","patchwork"
))
install.packages("spatstat")
```

flowchart TD
    A[Raw Spatial Data] --> B[CRS Harmonization]
    B --> C[Incident Date Filtering]
    C --> D[Nearest Road Matching]
    D --> E[Road-Level Aggregation]
    E --> F[Neighborhood Buffers]
    F --> G[Covariate Encoding]
    G --> H[Berman–Turner Quadrature]
    H --> I[Haar Basis Expansion]
    I --> J[Penalized Estimation]
    J --> K[Variable Selection]
    K --> L[Post-Selection BT Refit]

