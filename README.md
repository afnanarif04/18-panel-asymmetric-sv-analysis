# Asymmetric Panel Stochastic Volatility with Cross-Sectional Dependence
## A Kalman Filter Approach with Applications

---

## Overview

This repository contains the replication code and data for the paper:

The paper proposes the **Asymmetric Panel Stochastic Volatility (APSV)** model, a two-stage Kalman filter estimator that simultaneously accommodates:
- Leverage and size effects in log-volatility dynamics
- Multi-factor cross-sectional dependence via Common Correlated Effects (CCE) defactoring
- Heterogeneous factor loadings across panel units

---

## Repository Structure

```
├── README.md
├── code/
│   ├── 00_Functions.R                          # Core: C++ Kalman filter, CCE, QML, HPJ
│   ├── 00_Empirical_Utilities.R                # Data download helpers
│   ├── APSV_Demo_Simple.R                      # 5-minute standalone demo
│   ├── MonteCarlo/
│   │   ├── MC_01_Baseline.R                    # φ=0.95, γ₁=−0.15, γ₂=0.10
│   │   ├── MC_02_ConfigA.R                     # High persistence (φ=0.98)
│   │   ├── MC_03_ConfigB.R                     # Low persistence (φ=0.90)
│   │   ├── MC_04_ConfigC.R                     # Strong leverage, no size effect
│   │   └── MC_05_ConfigD.R                     # Strong leverage and size effect
│   ├── Figures/
│   │   ├── Figure_1_Baseline_Performance.R
│   │   ├── Figure_2_CCE_vs_Naive.R
│   │   ├── Figure_3_HPJ_Bias_Correction.R
│   │   └── Figure_4_Convergence_Rate.R
│   └── Empirical/
│       ├── Section_5.1_Crypto.R                # Cryptocurrency panel (N=7, T=2191)
│       ├── Section_5.2_MSCI.R                  # MSCI developed markets (N=12, T=6539)
│       └── Section_5.3_AsaiMcAleer2011.R       # Benchmark comparison
└── data/
    ├── APSV_Crypto_Data.xlsx                   # Cryptocurrency daily returns
    ├── APSV_MSCI_Data.xlsx                     # MSCI ETF daily returns
    └── APSV_AsaiMcAleer2011_Data.xlsx          # A&M 2011 replication dataset
```

---

## Quick Start

### Option 1: Run the 5-minute demo (recommended first step)
```r
# Open in RStudio, then run:
source("code/APSV_Demo_Simple.R")
# Runtime: ~5 minutes on a standard laptop
# Output: APSV_Demo_Result.pdf showing the key finding
```

### Option 2: Reproduce a single Monte Carlo configuration
```r
setwd("code/MonteCarlo")
source("../00_Functions.R")
source("MC_01_Baseline.R")
# Runtime: ~3-6 hours for 500 replications across 20 grid cells
# Output: Results/MC_Baseline.rds
```

### Option 3: Reproduce the empirical applications
```r
setwd("code/Empirical")
source("../00_Functions.R")
source("../00_Empirical_Utilities.R")
source("Section_5.1_Crypto.R")      # Section 5.1 — Cryptocurrency panel
source("Section_5.2_MSCI.R")        # Section 5.2 — MSCI equity panel
source("Section_5.3_AsaiMcAleer2011.R")  # Section 5.3 — Benchmark
# Runtime: 10-30 minutes per script (downloads data from Yahoo Finance)
```

---

## Software Requirements

- **R version:** >= 4.0.0
- **Required packages:**

```r
install.packages(c(
  "Rcpp", "ggplot2", "tidyr", "dplyr", "patchwork",
  "quantmod", "xts", "zoo", "moments", "viridis"
))
```

The C++ Kalman filter in `00_Functions.R` compiles automatically on first run via `Rcpp::sourceCpp()`. This takes approximately 30 seconds the first time and is cached thereafter.

---

## Data Sources

All data are freely available. The Excel files in `data/` contain the exact cleaned datasets used in the paper.

| Application | Source | Period | N | T |
|---|---|---|---|---|
| Cryptocurrency panel | Yahoo Finance | Jan 2019 – Dec 2024 | 7 | 2,191 |
| MSCI equity panel | Yahoo Finance (iShares ETFs) | Jan 1999 – Dec 2024 | 12 | 6,539 |
| A&M 2011 replication | Yahoo Finance + FRED | Jan 1990 – Dec 2007 | 3 | ~4,500 |

To re-download the data from scratch:
```r
source("code/00_Empirical_Utilities.R")
download_crypto_data(start = "2019-01-01", end = "2024-12-31")
download_msci_data(start = "1999-01-04", end = "2024-12-31")
```

---

## Model Summary

The APSV measurement and transition equations are:

```
y_it  = exp(h_it / 2) * ε_it

h_it  = μ_i + φ * h_{i,t-1} + γ₁ * ε_{i,t-1}
              + γ₂ * (|ε_{i,t-1}| − E|ε|)
              + λᵢ' * F_t + η_it
```

**Identified parameters:** `(φ, γ₁, σ²_u)` where `σ²_u = γ₁² + γ₂²(1−2/π) + σ²_e`

**Estimator:** Two-stage CCE-QML with half-panel jackknife bias correction (HPJ)

---

## Key Results

| Application | Pooled φ̂ | Pooled γ̂₁ | Key finding |
|---|---|---|---|
| Cryptocurrencies (N=7) | 0.967 | **+0.023** | Inverted leverage — positive shocks amplify volatility |
| MSCI equities (N=12) | 0.999 | **−0.047** | Standard equity leverage, strong country heterogeneity |
| A&M 2011 replication | 0.991 | **−0.010** | CCE removes >60% of apparent TOPIX leverage |

---
