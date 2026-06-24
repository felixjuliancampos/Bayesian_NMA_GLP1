# ==============================================================================
# 01_frequentist_nma.R
# GLP-1 RA / Incretin NMA — eGFR chronic annual slope
# Frequentist framework: netmeta (graph-theoretical NMA)
#
# Input:  data/first_eGFR_NMA_netmeta.xlsx
# Output: console + optional CSV exports (commented at end)
#
# Reference: Rücker G, Schwarzer G. Netmeta. R package. 2023.
# ==============================================================================

# ---- 0. Packages -------------------------------------------------------------
required <- c("netmeta", "readxl", "dplyr", "metafor", "ggplot2")
for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(netmeta); library(readxl); library(dplyr)
library(metafor); library(ggplot2)

# ---- 1. Load data ------------------------------------------------------------
xlsx_path <- "data/first_eGFR_NMA_netmeta.xlsx"   # adjust if needed

arm <- read_excel(xlsx_path, sheet = "arm_level")
cov <- read_excel(xlsx_path, sheet = "study_covariates")

cat("Loaded", nrow(arm), "arm rows from", n_distinct(arm$studlab), "studies.\n")

# ---- 2. Filter to primary NMA rows ------------------------------------------
# include_in_first_NMA == "yes" = primary analysis
# "sensitivity_pooled_node"     = toggle on for sensitivity (pooled-sema node)
arm_pop <- arm |>
  filter(!is.na(mean), !is.na(se), !is.na(n),
         include_in_first_NMA == "yes") |>
  mutate(sd = se * sqrt(n))   # pseudo-SD required by pairwise()

cat("\nArms in primary analysis:\n")
print(arm_pop |> select(studlab, node, n, mean, se, sd))

# ---- 3. Arm-level → pairwise contrast format --------------------------------
pw <- pairwise(
  treat   = node,
  n       = n,
  mean    = mean,
  sd      = sd,
  studlab = studlab,
  data    = arm_pop,
  sm      = "MD"
)

cat("\nPairwise contrasts (TE = MD in eGFR slope, mL/min/1.73m²/year):\n")
print(pw |> select(studlab, treat1, treat2, TE, seTE))

# ---- 3b. Verify connectivity before fitting ---------------------------------
nc <- netconnection(treat1  = pw$treat1,
                    treat2  = pw$treat2,
                    studlab = pw$studlab)
print(nc)
stopifnot("Network not connected — check include_in_first_NMA flags." =
            nc$n.subnets == 1)

# ---- 4. Primary NMA — random effects ----------------------------------------
nma <- netmeta(
  TE              = TE,
  seTE            = seTE,
  treat1          = treat1,
  treat2          = treat2,
  studlab         = studlab,
  data            = pw,
  sm              = "MD",
  common          = FALSE,
  random          = TRUE,
  reference.group = "Placebo",
  details.chkmultiarm = TRUE
)

cat("\n\n========== PRIMARY NMA RESULTS ==========\n")
summary(nma)

# ---- 5. Network graph -------------------------------------------------------
netgraph(
  nma,
  plastic          = FALSE,
  thickness        = "number.of.studies",
  points           = TRUE,
  cex.points       = 4,
  number.of.studies = TRUE,
  main             = "GLP-1 RA eGFR NMA — full network"
)

# ---- 6. Forest plot (all treatments vs Placebo) ----------------------------
forest(
  nma,
  reference.group = "Placebo",
  sortvar         = -Pscore,   # sort by P-score (higher = better eGFR preservation)
  smlab           = "MD in eGFR slope vs Placebo (mL/min/1.73m²/year)",
  xlim            = c(-5, 5),
  main            = "Frequentist NMA — random effects"
)

# ---- 7. League table --------------------------------------------------------
cat("\n\n========== LEAGUE TABLE (random effects) ==========\n")
league <- netleague(nma, common = FALSE, digits = 3)
print(league$random)

# ---- 8. P-scores (ranking) --------------------------------------------------
# Higher slope = less decline = BETTER → small.values = "bad"
cat("\n\n========== TREATMENT RANKINGS (P-score) ==========\n")
ranks <- netrank(nma, common = FALSE, small.values = "bad")
print(ranks)

# ---- 9. Heterogeneity & inconsistency ---------------------------------------
cat("\n\n========== HETEROGENEITY ==========\n")
cat(sprintf("tau  = %.4f\n", nma$tau))
cat(sprintf("tau² = %.4f\n", nma$tau^2))
cat(sprintf("I²   = %.1f%%\n", nma$I2 * 100))

cat("\n\n========== DESIGN-BY-TREATMENT INCONSISTENCY ==========\n")
decomp.design(nma)

cat("\n\n========== NODE-SPLIT (direct vs indirect) ==========\n")
ns <- netsplit(nma)
print(ns)
forest(ns, main = "Node-split: direct vs indirect evidence")

# ---- 10. Meta-regression on effect modifiers --------------------------------
cat("\n\n========== META-REGRESSION ==========\n")
# Merge study-level covariates into pairwise data
pw_cov <- pw |> left_join(cov, by = "studlab")

modifiers <- c(
  "mean_baseline_eGFR_mL.min.1.73m2",
  "log10_mean_baseline_UACR",
  "mean_baseline_HbA1c_pct",
  "mean_baseline_BMI_kg.m2",
  "mean_DM_duration_years",
  "pct_baseline_SGLT2i"
)

for (m in modifiers) {
  cat("\n---------- Modifier:", m, "----------\n")
  n_avail <- sum(!is.na(cov[[m]]))
  if (n_avail < 3) {
    cat("Insufficient data (n =", n_avail, "). Skipping.\n"); next
  }
  fit <- tryCatch(
    metareg(nma, as.formula(paste("~", m))),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    cat("metareg failed:", fit$message, "\n")
  } else {
    print(summary(fit))
  }
}

# ---- 11. Exports (uncomment to save) ----------------------------------------
# write.csv(league$random, "outputs/league_table_random.csv")
# write.csv(as.data.frame(ranks), "outputs/pscores.csv")
# write.csv(pw |> select(studlab, treat1, treat2, TE, seTE),
#           "outputs/pairwise_contrasts.csv", row.names = FALSE)
# saveRDS(nma, "outputs/nma_object.rds")

cat("\n\n========== SESSION INFO ==========\n")
sessionInfo()
