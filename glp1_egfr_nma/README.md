# GLP-1 RA eGFR NMA — R Analysis Code

Network meta-analysis of GLP-1 receptor agonists and incretin-based therapies
on chronic annual eGFR slope (mL/min/1.73 m²/year) in type 2 diabetes.

**Primary outcome:** Mean annual eGFR slope from random-intercept +
random-coefficient mixed models (mL/min/1.73 m²/year; negative = decline;
less negative = better renal preservation).

**Framework:** Cochrane systematic review methodology; PRISMA-NMA reporting.

---

## Repository structure

```
glp1_egfr_nma/
├── README.md
├── data/
│   └── first_eGFR_NMA_netmeta.xlsx   # master NMA workbook (arm_level +
│                                      # study_covariates + treatment_dictionary)
├── 01_frequentist_nma.R               # netmeta frequentist NMA (full pipeline)
├── 02_bayesian_nma_multinma.R         # Bayesian NMA via multinma + Stan (primary)
├── 03_sensitivity_analyses.R          # double-blind sub-network + other sensitivities
└── 04_figures.R                       # publication-ready figures
```

---

## Installation

### Frequentist (01_frequentist_nma.R)
```r
install.packages(c("netmeta", "readxl", "dplyr", "metafor", "ggplot2"))
```

### Bayesian (02_bayesian_nma_multinma.R)
```r
# Step 1 — cmdstanr (Stan backend; no JAGS required)
install.packages("cmdstanr",
                 repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()   # one-time, ~3–5 min

# Step 2 — multinma
install.packages("multinma")

# Verify
library(multinma)
```

> **Note:** JAGS is NOT required. These scripts use Stan via `cmdstanr`.
> Tested on macOS (Apple Silicon + Intel) and Linux with R ≥ 4.2.

---

## Key methodological notes

- **Node definition:** drug + dose + route + frequency. Doses are never pooled
  across nodes (dulaglutide 0.75 mg ≠ 1.5 mg; SC semaglutide ≠ oral semaglutide).
- **Estimand:** random-coefficient chronic annual slope. Annualised-total
  derivations and 52-week changes are flagged in `include_in_first_NMA`.
- **Heterogeneity:** τ = 1.78 (full network, vague prior); τ = 0.17
  (double-blind sub-network, informative HalfNormal(0.25) prior).
- **Inconsistency:** Bayesian node-split — all 95% CrIs include zero; largest
  ω = −2.09 (Insulin_glargine vs Dulaglutide_1.5, AWARD-7 bridge).
- **Primary interpretable result:** double-blind sub-network (script 03).

---

## Session info (at time of analysis)

- R 4.6
- multinma (version as installed)
- cmdstanr / CmdStan
- netmeta ≥ 2.0
- readxl, dplyr, metafor, ggplot2

---

## Citation

Campos-Garcia FJ et al. [Title]. Cochrane Mexico. [Year in preparation].

Data extraction form: `GLP1_NMA_DataExtractionForm.xlsx`
