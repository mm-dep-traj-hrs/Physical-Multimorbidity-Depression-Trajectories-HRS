# Physical Multimorbidity Patterns and Depressive Symptom Trajectories Among Older U.S. Adults

Supplementary code for *Prevalent Physical Multimorbidity Patterns and Depressive Symptom Trajectories among Community Dwelling U.S. Older Adults*.

**Author:** Nicholas Bishop, University of Arizona  
**Contact:** njbishop@arizona.edu  
**Last updated:** 2026-05-12

## Overview

This repository contains the R and Mplus output files used to construct the analytic sample, generate inverse probability of dropout weights (IPW), and estimate latent growth models of depressive symptom trajectories across somatic multimorbidity patterns using the Health and Retirement Study (HRS), 2012–2020.

## Data sources

Analyses use two HRS files:

- **RAND HRS Longitudinal File 2022 (v1):** `randhrs1992_2022v1.sav`
- **HRS Childhood Health and Family file:** `AGGCHLDFH2016A_R.sav`

Both files require registration and approval at https://hrsdata.isr.umich.edu/. Update file paths in the R script to match your local directory.

## Repository structure

- **R code for analysis.R**          Sample construction, data management, variable coding, IPW creation
- **Mplus\model fit testing**        Mplus output files for model fit testing of the unconditional growth model using complete analytic sample
- **Mplus\model fit testing ipw**    Mplus output files for model fit testing of the unconditional growth model using IPW-trimmed sample
- **Mplus\growth models**            Mplus output files for primary and IPW-adjusted latent growth models

## Software

- R version 4.5.2, with packages listed at the top of `R code for analysis.R`
- Mplus version 8.8

## Notes

- `MplusAutomation` citation: Hallquist, M. N., & Wiley, J. F. (2018). MplusAutomation: An R package for facilitating large-scale latent variable analyses in Mplus. *Structural Equation Modeling*, *25*(4), 621–638. https://doi.org/10.1080/10705511.2017.1402334
- The IPW workflow uses the `missForest` package for covariate imputation prior to weight estimation. Cores for parallelization are set in the script and may need adjustment for your environment.

## Citation

If you use this code, please cite:

Bishop, N. J., Walker, K. J., Nagel, C. L., Newsom, J. T., Botoseneanu, A., Allore, H. G., Triolo, F., & Quiñones, A. R. (UPDATE). Prevalent physical multimorbidity patterns and depressive symptom trajectories among community dwelling U.S. older adults. *Aging & Mental Health*.
