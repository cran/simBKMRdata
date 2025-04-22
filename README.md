---
title: "simBKMRdata"
author: "Kazi Tanvir Hasan and Gabriel Odom"
---

Simulate multivariate normal or multivariate skewed exposure data for downstream *in silico* experiments with [Bayesian Kernel Machine Regression](https://CRAN.R-project.org/package=bkmr).

## Outline of Functions

1. (Graphically) Assess Skewness of Exposure variables $z_i$: this part we will show as examples in the vignettes, but we won't include functions.
2. Transformation Functions:
    + To Normality: $(z - \bar{z})/\text{sd}(z)$; $\log_{b}\left[(z - \bar{z})/\text{sd}(z)\right]$
    + To Gamma: $z/\text{sd}(z)$; $\log_{b}\left[z/\text{sd}(z)\right]$; $\log_{b}\left[(z+1)/\text{sd}(z+1)\right]$
3. Calculate List of Parameters for Groups $i = 1, \ldots, G$
    + MV Normal: $n_i$, $\hat{\boldsymbol\mu}_i$, $\hat{\boldsymbol\Sigma}_i$
    + MV Skew Gamma: $n_i$, $\hat{\alpha}_i$, $\hat{\beta}_i$, $\hat{\mathbb{P}}_i$ (group Spearman correlation matrix Rho)
4. Simulate List of Exposure Data Sets for Groups $i = 1, \ldots, G$
    + MV Normal: use [`MASS::mvrnorm()`](https://rdrr.io/cran/MASS/man/mvrnorm.html)
    + MV Skew Gamma: copy the internals of [`lcmix::rmvgamma()`](https://rdrr.io/rforge/lcmix/man/mvgamma.html) from RForge, cite it, and include it here
5. Use BKMR to analyze the simulated data (we will show some quick examples in vignettes, but not include any functions)
6. Calculate a PIP threshold that preserves a 5% test size for real or simulated data (as close as we can for now). This function should only depend on the response vector, or summary statistics of it (specifically $|\text{cv}(y)|$ and $n$). 




