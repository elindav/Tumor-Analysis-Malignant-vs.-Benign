# Breast Nuclear Area: A Gamma Distribution Analysis

### Statistical Modeling, Parameter Estimation, and Confidence Interval Analysis

<p align="center">
  <strong>AUA · CS 108 — Statistics · Fall 2026</strong>
</p>

---

## Overview

This project investigates whether mean cell-nucleus area measurements can help quantify differences between malignant and benign breast-mass observations.

Using the **Wisconsin Diagnostic Breast Cancer Dataset**, we apply Gamma distribution models to analyze the two diagnostic groups separately. The project focuses on distributional modeling, parameter estimation, estimator properties, and confidence interval construction.

In particular, we investigate how accurately the population mean nuclear area among malignant observations can be estimated and evaluate the performance of two confidence interval methods through simulation.

## Research Questions

The analysis addresses the following questions:

- How do mean nuclear area measurements differ between malignant and benign observations?
- How well does the Gamma distribution describe the observed measurements?
- How do the Method of Moments (MoM) and Maximum Likelihood Estimation (MLE) perform as parameter estimation methods?
- How accurately can the population mean nuclear area for malignant observations be estimated?
- How does the empirical coverage of asymptotic confidence intervals change with sample size?

## Dataset

The project uses the **Wisconsin Diagnostic Breast Cancer (WDBC) Dataset**, obtained from the UCI Machine Learning Repository and accessed through the `scikit-learn` distribution.

The dataset contains 569 diagnostic observations, each with 30 numerical features describing cell-nucleus characteristics.

This analysis focuses on `area_mean`, the mean nuclear area per observation.

| Diagnostic group | Observations |
|:---|---:|
| Malignant | 212 |
| Benign | 357 |
| **Total** | **569** |

The measurements describe cell nuclei extracted from fine-needle aspirate (FNA) images, not the physical dimensions of entire tumors.

## Statistical Methodology

### 1. Exploratory Data Analysis

We examine the distribution of mean nuclear area across both diagnostic groups using:

- Descriptive statistics
- Comparative boxplots
- Density histograms
- Theoretical Gamma Q–Q plots
- Empirical and fitted cumulative distribution functions

### 2. Gamma Distribution Modeling

Separate two-parameter Gamma distributions are fitted to malignant and benign observations.

The probability density function is:

\[
f(x;k,\beta)=\frac{\beta^k}{\Gamma(k)}
x^{k-1}e^{-\beta x}, \quad x>0
\]

where:

- \(k\) represents the shape parameter.
- \(\beta\) represents the rate parameter.

The population mean and variance are:

\[
E(X)=\frac{k}{\beta},
\qquad
Var(X)=\frac{k}{\beta^2}
\]

### 3. Parameter Estimation

Two estimation approaches are investigated:

- **Method of Moments (MoM):** Estimates the Gamma parameters by matching theoretical and empirical moments.
- **Maximum Likelihood Estimation (MLE):** Estimates parameters by maximizing the likelihood function.

The MLE shape parameter is obtained numerically using R's `uniroot()` function.

Their theoretical properties and finite-sample performance are also examined through simulation.

### 4. Confidence Interval Construction

Two asymptotic 95% confidence interval methods are constructed for the malignant population mean:

- **Wald (Normal) Interval**
- **Log-Transformed Interval**

The methods are compared in terms of interval endpoints, length, and empirical coverage.

### 5. Simulation Study

Monte Carlo simulation is used to evaluate estimator performance and confidence interval coverage.

The confidence interval simulation uses:

- 5,000 repetitions per sample size.
- Sample sizes of \(n=30,\ 212,\ 500\).
- A fixed random seed: `set.seed(108)`.
- The malignant MLE estimates as fixed simulation parameters.

The analysis examines empirical coverage probability and average interval length.

## Key Results

### Comparison of Diagnostic Groups

| Statistic | Malignant | Benign |
|:---|---:|---:|
| Sample size | 212 | 357 |
| Mean nuclear area | 978.38 | 462.79 |
| Median | 932.00 | 458.40 |
| Standard deviation | 367.94 | 134.29 |
| Estimated Gamma shape (MLE) | 7.65 | 11.31 |
| Estimated Gamma rate (MLE) | 0.00782 | 0.02444 |

The estimated mean nuclear area in the malignant group is approximately **2.11 times** that of the benign group.

The malignant observations also exhibit greater variability and greater relative dispersion under the fitted Gamma model.

### 95% Confidence Intervals

For the malignant population mean:

\[
\hat{\mu}_M=978.38
\]

| Method | Lower bound | Upper bound |
|:---|---:|---:|
| Wald (Normal) | 930.77 | 1025.98 |
| Log-Transformed | 931.91 | 1027.16 |

Both methods produce very similar intervals for the observed sample.

### Simulation: Empirical Coverage

| Sample size | Wald | Log-Transformed |
|---:|---:|---:|
| 30 | 93.24% | 93.62% |
| 212 | 94.66% | 94.74% |
| 500 | 95.34% | 95.48% |

The simulation suggests that both methods approach their nominal 95% coverage as sample size increases. Their performance is similar, while the log-transformed interval additionally guarantees positive endpoints.

## Main Findings

- Separate Gamma distributions provide reasonable approximations to the observed nuclear-area measurements.
- The malignant group has a substantially higher average nuclear area than the benign group.
- MLE generally demonstrates lower empirical variance and mean squared error than MoM in the estimator simulation.
- Both confidence interval methods perform similarly, with coverage approaching 95% as sample size increases.
- The Gamma model fits the central portions of the distributions better than some tail regions, particularly for malignant observations.

## Limitations

The results should be interpreted with several limitations in mind:

- The Gamma distributions are working approximations rather than exact descriptions of the data.
- Some discrepancies remain in the distributional tails.
- The observations are treated as independent and identically distributed.
- The dataset is not documented as a random sample of all breast-mass patients.
- The analysis concerns mean cell-nucleus area, not the dimensions of entire tumors.
- Simulation results assume the fitted Gamma parameters represent the true generating parameters for the simulated samples.

These limitations restrict the extent to which the findings can be generalized to broader clinical populations.

## Technologies and Tools

- **R** — Statistical computing and analysis
- **R Markdown** — Reproducible statistical reporting
- **Numerical Methods** — Numerical solution of the MLE shape equation
- **Monte Carlo Simulation** — Evaluation of estimator properties and confidence interval coverage
- **Statistical Visualization** — Distribution plots, Q–Q plots, CDFs, and simulation graphs

## References

1. Wolberg, W. H., Street, W. N., & Mangasarian, O. L. *Breast Cancer Wisconsin (Diagnostic).* UCI Machine Learning Repository.  
   https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic

2. Scikit-learn Developers. *Breast Cancer Wisconsin (Diagnostic) Dataset.*  
   https://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_breast_cancer.html

## Authors

**Elina Davtyan · Nana Budaghyan · Lusine Norekyan**

---

<p align="center">
  <sub>Academic statistical analysis · Fall 2026</sub>
</p>
