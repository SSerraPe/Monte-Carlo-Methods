# Monte Carlo and Importance Sampling Estimators

This repository contains a simulation study of Monte Carlo and Importance Sampling (IS) methods for estimating expectations and tail probabilities under various target distributions using **R**.

## Overview

We analyze the empirical performance of different estimators by computing:

- **Bias**
- **Variance**
- **Mean Squared Error (MSE)**

The simulations cover:

1. **Raw Monte Carlo (MC)** for a standard normal distribution.
2. **Importance Sampling** using a bimodal target:  
   \( p(x) = 0.5 \mathcal{N}(-3, 1) + 0.5 \mathcal{N}(5, 4) \)
3. **Comparison of UIS vs. SNIS** estimators.
4. **Sensitivity analysis** to the mean and variance of the proposal distribution.
5. **Extension to higher dimensions**, examining MSE scaling and weight degeneracy.

## Key Estimands

- \( \mathbb{E}[X] \)
- \( \mathbb{E}[X^2] \)
- \( \Pr(X > \gamma) \)

## Features

- Monte Carlo estimator characterization for multiple sample sizes
- Importance sampling visualization with good and bad proposals
- Empirical comparison between UIS and SNIS
- Grid-based evaluation of proposal parameter impact
- Performance degradation analysis in higher dimensions
- Annotated plots using **ggplot2** and **patchwork**

## Requirements

This project uses the following R packages:

```r
tidyverse
ggplot2
patchwork
scales
gridExtra
flextable
officer
````

Install them via:

```r
install.packages(c("tidyverse", "ggplot2", "patchwork", "scales", "gridExtra", "flextable", "officer"))
```

## Usage

Simply run the main R script to:

* Simulate samples from the target and proposal distributions.
* Compute performance metrics.
* Generate all visualizations and tables.

## Output

The script generates:

* Histograms and trajectories of estimates
* Bias, variance, and MSE vs. sample size
* Comparisons of UIS vs. SNIS
* Sensitivity tables of proposals
* MSE plots by dimension
* Sorted IS weights (to assess degeneracy)

## License

MIT License

## Author

\[Your Name] – 2025

```
