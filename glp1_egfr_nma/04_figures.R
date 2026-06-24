# ==============================================================================
# 04_figures.R
# GLP-1 RA / Incretin NMA — publication-ready figures
#
# Requires objects from 02 and 03. Source those scripts first, OR load RDS:
#   fit_full_vague <- readRDS("outputs/bayes_full_vague.rds")
#   fit_full_inf   <- readRDS("outputs/bayes_full_informative.rds")
#   fit_db         <- readRDS("outputs/bayes_double_blind.rds")
#   fit_ns         <- readRDS("outputs/bayes_nodesplit.rds")
# ==============================================================================

library(ggplot2); library(multinma); library(dplyr); library(patchwork)

# ---- FIGURE 1 — Double-blind NMA forest plot --------------------------------
# Primary result: Bayesian double-blind sub-network, informative tau prior
# Requires: fit_db, het_db from 03_sensitivity_analyses.R

fig1 <- plot(relative_effects(fit_db), ref_line = 0) +
  labs(
    title    = "Figure 1. Bayesian NMA — eGFR slope vs Placebo",
    subtitle = paste0(
      "Double-blind placebo-controlled trials (k=4) | ",
      "Informative \u03c4 prior [HN(0.25)] | ",
      "\u03c4 median = 0.17 (95% CrI 0.01\u20130.56)"
    ),
    x        = "MD in eGFR slope vs Placebo (mL/min/1.73m\u00b2/year)\nPositive = less decline (better)",
    caption  = paste0(
      "Studies: REWIND (Botros 2023), LEADER (Mann 2017), FLOW (Perkovic 2024), ",
      "SURPASS-CVOT (Zoungas 2026).\n",
      "Tirzepatide estimate is indirect via SURPASS-CVOT (Dulaglutide_1.5 bridge)."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(plot.caption = element_text(size = 8, hjust = 0))
print(fig1)

# ---- FIGURE 2 — Tau prior sensitivity ---------------------------------------
# Shows posterior tau across all three priors + frequentist reference

# Build tau summary table from all models
# (run 02 first to have fit_full_vague, fit_full_weak, fit_full_inf)
make_tau_row <- function(fit, label) {
  het <- summary(fit, pars = "tau")$summary
  data.frame(
    prior  = label,
    median = het["tau","50%"],
    lo95   = het["tau","2.5%"],
    hi95   = het["tau","97.5%"]
  )
}

tau_priors <- bind_rows(
  make_tau_row(fit_full_vague, "Full network\nVague [HN(5)]"),
  make_tau_row(fit_full_weak,  "Full network\nWeakly inf. [HN(1)]"),
  make_tau_row(fit_full_inf,   "Full network\nInformative [HN(0.25)]"),
  data.frame(prior  = "Double-blind\nsub-network [HN(0.25)]",
             median = 0.17, lo95 = 0.01, hi95 = 0.56)
)
tau_priors$prior <- factor(tau_priors$prior, levels = rev(tau_priors$prior))

fig2 <- ggplot(tau_priors, aes(x = prior, y = median)) +
  geom_pointrange(aes(ymin = lo95, ymax = hi95),
                  size = 0.9, linewidth = 1,
                  colour = c("#C00000","#C00000","#C00000","#2E75B6")) +
  geom_hline(yintercept = 1.1323, linetype = "dashed",
             colour = "grey50", linewidth = 0.8) +
  annotate("text", x = 1.6, y = 1.25,
           label = "Frequentist tau = 1.13\n(DerSimonian-Laird)",
           size = 3, colour = "grey40", hjust = 0) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 7),
                     breaks = c(0, 0.5, 1, 2, 3, 4, 5, 6)) +
  labs(
    title    = "Figure 2. Posterior heterogeneity (\u03c4) by model",
    subtitle = "Median and 95% CrI | Dashed = frequentist DerSimonian-Laird estimate",
    x        = NULL,
    y        = "\u03c4 (mL/min/1.73m\u00b2/year)",
    caption  = paste0(
      "Blue = double-blind sub-network (primary). Red = full network with three prior choices.\n",
      "Extreme tau in the full network reflects AWARD-7 / GRADE structural heterogeneity",
      " (baseline eGFR 35-95; mixed estimands)."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(plot.caption = element_text(size = 8, hjust = 0))
print(fig2)

# ---- FIGURE 3 — Node-split forest (direct vs indirect) ---------------------
# Requires: fit_ns from 02_bayesian_nma_multinma.R

fig3 <- plot(summary(fit_ns)) +
  labs(
    title    = "Figure 3. Node-split: direct vs indirect evidence",
    subtitle = paste0(
      "\u03c9 = direct \u2212 indirect | 95% CrI including zero = no significant inconsistency\n",
      "Largest \u03c9 = \u22122.09 (Insulin_glargine vs Dulaglutide_1.5): AWARD-7 bridge"
    ),
    x        = "MD in eGFR slope (mL/min/1.73m\u00b2/year)",
    caption  = paste0(
      "Full network with vague tau prior [HN(5)]. ",
      "Frequentist inconsistency test: Q=4.18, df=1, p=0.041.\n",
      "All Bayesian node-split 95% CrIs include zero (underpowered with k=6 studies)."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(plot.caption = element_text(size = 8, hjust = 0))
print(fig3)

# ---- FIGURE 4 — Baseline eGFR gradient (manuscript key signal) -------------
# Bar/dot chart showing the direct placebo-anchored MDs vs baseline eGFR
# This is the key clinical finding independent of NMA heterogeneity

gradient_data <- data.frame(
  trial        = c("LEADER\n(Liraglutide 1.8)",
                   "REWIND\n(Dulaglutide 1.5)",
                   "SUSTAIN-6+PIONEER-6\n(Semaglutide pooled)",
                   "FLOW\n(Semaglutide 1.0)"),
  baseline_egfr = c(80, 77, 75, 47),
  MD            = c(0.127, 0.190, 0.590, 1.165),
  lo95          = c(-0.43, 0.13, 0.27, 0.86),
  hi95          = c(0.68, 0.25, 0.91, 1.47),
  blinding      = c("Double-blind", "Double-blind", "Pooled*", "Double-blind")
)

fig4 <- ggplot(gradient_data,
               aes(x = baseline_egfr, y = MD,
                   colour = blinding, shape = blinding)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
                width = 1.5, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = trial), hjust = -0.12, size = 3, colour = "grey30") +
  scale_x_reverse(limits = c(95, 35),
                  breaks  = c(80, 70, 60, 50, 40)) +
  scale_colour_manual(values = c("Double-blind" = "#2E75B6",
                                 "Pooled*"       = "#ED7D31")) +
  scale_shape_manual( values = c("Double-blind" = 16, "Pooled*" = 17)) +
  labs(
    title    = "Figure 4. Baseline eGFR gradient — direct placebo-anchored evidence",
    subtitle = "Larger GLP-1 RA benefit at lower baseline eGFR (pre-specified effect modifier)",
    x        = "Mean baseline eGFR (mL/min/1.73m\u00b2)",
    y        = "MD in annual eGFR slope vs Placebo\n(mL/min/1.73m\u00b2/year)",
    colour   = NULL, shape = NULL,
    caption  = paste0(
      "FLOW: Perkovic 2024. REWIND: Botros 2023. LEADER: Mann 2017.\n",
      "*Pooled: Tuttle 2023 (SUSTAIN-6 + PIONEER-6; pools SC 0.5, SC 1.0, oral 14 mg — ",
      "distinct from primary Semaglutide_1.0_QW node)."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 8, hjust = 0))
print(fig4)

# ---- FIGURE 5 — Combined manuscript panel (patchwork) ----------------------
# Combines Fig1 (forest) + Fig2 (tau) side by side for a 2-panel manuscript figure
fig_panel <- (fig1 | fig2) +
  plot_annotation(
    title  = "GLP-1 RA / Incretin NMA — Bayesian analysis",
    tag_levels = "A"
  )
print(fig_panel)

# ---- Save figures (uncomment) -----------------------------------------------
# ggsave("outputs/fig1_forest_doubleblind.pdf",     fig1,      w=8, h=5)
# ggsave("outputs/fig2_tau_prior_comparison.pdf",   fig2,      w=8, h=5)
# ggsave("outputs/fig3_nodesplit.pdf",              fig3,      w=8, h=6)
# ggsave("outputs/fig4_egfr_gradient.pdf",          fig4,      w=8, h=5)
# ggsave("outputs/fig_panel_AB.pdf",                fig_panel, w=14,h=6)

cat("Figures 1-5 generated.\n")
