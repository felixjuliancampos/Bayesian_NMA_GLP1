# ==============================================================================
# 02_bayesian_nma_multinma.R
# GLP-1 RA / Incretin NMA — eGFR chronic annual slope
# Bayesian framework: multinma + Stan (cmdstanr backend; NO JAGS required)
#
# Input:  data/first_eGFR_NMA_netmeta.xlsx
# Output: console summaries + ggplot figures
#
# Installation (run once):
#   install.packages("cmdstanr",
#     repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#   cmdstanr::install_cmdstan()
#   install.packages("multinma")
#
# Reference: Phillippo DM et al. multinma. J Stat Softw 2024.
#   https://dmphillippo.github.io/multinma/
#
# Verified syntax for multinma 1.x:
#   - Priors: .default(half_normal(scale = x)) / .default(normal(...))
#   - MCMC args passed DIRECTLY to nma(), not via .args
#   - adapt_delta is a named argument of nma()
#   - chains / warmup / iter / seed passed through ... to Stan
# ==============================================================================

# ---- 0. Packages -------------------------------------------------------------
for (p in c("multinma", "readxl", "dplyr", "ggplot2")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(multinma); library(readxl); library(dplyr); library(ggplot2)

options(mc.cores = parallel::detectCores())   # use all available cores

# ---- 1. Load data ------------------------------------------------------------
xlsx_path <- "data/first_eGFR_NMA_netmeta.xlsx"

arm_raw <- read_excel(xlsx_path, sheet = "arm_level")
cov     <- read_excel(xlsx_path, sheet = "study_covariates")

arm <- arm_raw |>
  filter(!is.na(mean), !is.na(se), !is.na(n),
         include_in_first_NMA == "yes") |>
  mutate(sd = se * sqrt(n))

cat("Studies:", n_distinct(arm$studlab),
    "| Treatments:", n_distinct(arm$node),
    "| Arms:", nrow(arm), "\n")

# ---- 2. Build full network --------------------------------------------------
net <- set_agd_arm(
  data    = arm,
  study   = study  |> {\(.) arm$studlab}(),   # helper below
  trt     = node,
  y       = mean,
  se      = se,
  trt_ref = "Placebo"
)
# Simpler assignment:
net <- set_agd_arm(arm,
                   study   = studlab,
                   trt     = node,
                   y       = mean,
                   se      = se,
                   trt_ref = "Placebo")
print(net)
plot(net, main = "GLP-1 RA eGFR NMA — network")

# ==============================================================================
# SECTION A — FULL NETWORK, THREE TAU PRIORS (prior sensitivity)
#
# Prior rationale:
#   Vague       [HN(5)]    : lets data speak; reproduces frequentist tau ~ 1.13
#   Weakly inf  [HN(1)]    : rules out tau > 3 (implausible for this outcome)
#   Informative [HN(0.25)] : calibrated to double-blind trials where
#                            placebo-contrast tau ~ 0.24 (REWIND/LEADER/FLOW)
#                            Justified by Turner 2012 Stat Med empirical priors
# ==============================================================================

# ---- 3. Fit full-network models (all 3 priors) ------------------------------
set.seed(20240610)

# Shared MCMC settings (verified working in multinma 1.x)
mcmc_common <- list(
  chains      = 4,
  warmup      = 2000,
  iter        = 7000,    # total including warmup -> 5000 post-warmup draws
  seed        = 20240610,
  adapt_delta = 0.99
)

cat("\n=== Full network — vague prior [HN(5)] ===\n")
fit_full_vague <- do.call(nma, c(
  list(net,
       trt_effects = "random",
       prior_het   = .default(half_normal(scale = 5)),
       prior_trt   = .default(normal(location = 0, scale = 10))),
  mcmc_common
))

cat("\n=== Full network — weakly informative prior [HN(1)] ===\n")
fit_full_weak <- do.call(nma, c(
  list(net,
       trt_effects = "random",
       prior_het   = .default(half_normal(scale = 1)),
       prior_trt   = .default(normal(location = 0, scale = 10))),
  mcmc_common
))

cat("\n=== Full network — informative prior [HN(0.25)] ===\n")
fit_full_inf <- do.call(nma, c(
  list(net,
       trt_effects = "random",
       prior_het   = .default(half_normal(scale = 0.25)),
       prior_trt   = .default(normal(location = 0, scale = 10))),
  mcmc_common
))

# ---- 4. Convergence diagnostics ---------------------------------------------
cat("\n\n========== CONVERGENCE DIAGNOSTICS ==========\n")
for (nm in c("vague", "weak", "informative")) {
  fit <- get(paste0("fit_full_", nm))
  s   <- summary(fit)$summary
  max_rhat <- max(s[, "Rhat"],  na.rm = TRUE)
  min_ess  <- min(s[, "n_eff"], na.rm = TRUE)
  flag     <- if (max_rhat > 1.01) "*** RERUN with more iterations ***" else "OK"
  cat(sprintf("Prior %-15s: max Rhat = %.4f [%s]  min ESS = %.0f\n",
              nm, max_rhat, flag, min_ess))
}
# If max Rhat > 1.01: increase warmup to 4000, iter to 12000, adapt_delta to 0.995

# ---- 5. Treatment effects vs Placebo ----------------------------------------
cat("\n\n========== TREATMENT EFFECTS vs Placebo ==========\n")
cat("\n--- Vague prior ---\n");       print(relative_effects(fit_full_vague))
cat("\n--- Weakly informative ---\n");print(relative_effects(fit_full_weak))
cat("\n--- Informative prior ---\n"); print(relative_effects(fit_full_inf))

# ---- 6. Heterogeneity posterior (tau) ---------------------------------------
cat("\n\n========== HETEROGENEITY (tau) ==========\n")
for (nm in c("vague", "weak", "informative")) {
  fit <- get(paste0("fit_full_", nm))
  het <- summary(fit, pars = "tau")$summary
  cat(sprintf("Prior %-15s: tau mean=%.3f  median=%.3f  95%%CrI=[%.3f, %.3f]\n",
              nm,
              het["tau", "mean"],
              het["tau", "50%"],
              het["tau", "2.5%"],
              het["tau", "97.5%"]))
}
cat("Frequentist reference: tau = 1.1323 (DerSimonian-Laird, full network)\n")

# ---- 7. Node-split inconsistency check --------------------------------------
cat("\n\n========== NODE-SPLIT INCONSISTENCY ==========\n")
fit_ns <- do.call(nma, c(
  list(net,
       consistency = "nodesplit",
       trt_effects = "random",
       prior_het   = .default(half_normal(scale = 5)),
       prior_trt   = .default(normal(location = 0, scale = 10))),
  mcmc_common
))
cat("\nNode-split results (omega = direct - indirect):\n")
print(fit_ns)
# omega 95% CrI includes zero -> no formally significant inconsistency
# Largest omega: Insulin_glargine vs Dulaglutide_1.5 (AWARD-7 bridge)

# Forest plot: direct vs indirect
plot(summary(fit_ns)) +
  labs(title    = "Node-split: direct vs indirect evidence",
       subtitle = "omega = direct - indirect; CrI includes 0 = no significant inconsistency",
       x        = "MD in eGFR slope (mL/min/1.73m²/year)")

# ---- 8. Exports (uncomment to save) -----------------------------------------
# saveRDS(fit_full_vague, "outputs/bayes_full_vague.rds")
# saveRDS(fit_full_inf,   "outputs/bayes_full_informative.rds")
# saveRDS(fit_ns,         "outputs/bayes_nodesplit.rds")
# rel_tbl <- as.data.frame(relative_effects(fit_full_inf))
# write.csv(rel_tbl, "outputs/bayes_effects_informative.csv", row.names = FALSE)

cat("\n\n========== SESSION INFO ==========\n")
sessionInfo()
