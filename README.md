# geostat-survey
This repository contains the R codes used for calculations and figure plotting related to our research:

> Yalcin, S., Anderson, S., Regular, P., English, P. (2023). Exploring the limits of spatiotemporal and design-based index standardization under reduced survey coverage. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsad155

**Note**: The code and data generated from this research are being archived on Zenodo [DOI: 10.5281/zenodo.8326526].

## Main Scripts:

**1_cod_like_indices_calculations.R**
  - Description: Calculation of cod-like species design and model-based indices.

**2_yellowtail_like_indices_calculations.R**
  - Description: Calculation of yellowtail-like species design and model-based indices.

**3_main_figures.R**
   - Description: Script used for main figure plotting.

## Appendices:

- **4_appendix_1_additional_figures.R**
  - Description: Script containing additional figures for Appendix 1.

- **5_appendix_2_testing_sampling_variation.R**
  - Description: Codes for testing sampling variation in Appendix 2.

- **6_appendix_3_recovery_spillover.R**
  - Description: Script for mapping the base scenario, recovery of the closed area, and recovery with spillover effect in Appendix 3.
 
## Functions:

- **bootstrapping_fn.R**
  - Description: Bootstrap function for design-based indices.

- **data_prep_fn.R**
  - Description: Data preparation functions.

- **model_run_fn.R**
  - Description: Model functions.

- **pop_cod_fn.R**
  - Description: Functions related to cod-like population simulations.
 
- **pop_yellowtail_fn.R**
  - Description: Functions related to yellowtail-like population simulations.
 
## Data folder:
  - This folder includes the data generated through the analysis and can be used for plotting.
