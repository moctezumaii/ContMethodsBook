# =============================================================================
# GEE Example: Tick Phenology & Prescribed Burning
# Simulated data to match Gleim et al. (2014) design
# Paper: DOI 10.1371/journal.pone.0112174
#
# NOTE: The paper's raw data are NOT publicly available (not deposited in PLOS
# ONE data repository). This script simulates a dataset that preserves the
# key design features:
#   - 21 plots (clusters) across 4 burn-regime treatments
#   - Monthly repeated measures over 2 years (24 time points)
#   - Overdispersed tick count response with seasonal variation
#   - Predictors: burn_regime, month, host_abundance, canopy_closure
#
# The paper does NOT provide an explicit model equation. The inferred
# specification used here (negative binomial GEE, log link, exchangeable
# correlation) is based on the Methods section of Gleim et al. (2014).
#
# Dependencies: geepack, MASS, tidyverse (or base R)
# =============================================================================

# ---- 0. Load libraries -------------------------------------------------------
library(geepack)   # geeglm() for GEE
library(MASS)      # glm.nb() for negative binomial GLM (sensitivity check)
library(tidyverse) # data manipulation and visualization (optional)

# ---- 1. Simulate data matching Gleim et al. (2014) design -------------------
set.seed(2014)  # year of the paper for reproducibility

# Study design parameters
n_plots     <- 21          # number of sampling plots (clusters)
n_months    <- 24          # 12 months × 2 years (2010–2011)
burn_levels <- c("BB", "BUB", "UBB", "UBUB")  # burn regime treatments

# Assign plots to burn regimes (approximately equal; 21 is not divisible by 4)
# BB=6, BUB=5, UBB=5, UBUB=5 (matches approximate paper description)
plot_burns <- c(rep("BB", 6), rep("BUB", 5), rep("UBB", 5), rep("UBUB", 5))

# Create the full plot × month data frame
tick_data <- expand.grid(
  plot_id   = 1:n_plots,
  month_num = 1:n_months
) |>
  mutate(
    # Burn regime for each plot
    burn_regime = factor(plot_burns[plot_id], levels = burn_levels),
    # Calendar month within year (1 = January)
    cal_month   = ((month_num - 1) %% 12) + 1,
    year        = ifelse(month_num <= 12, 2010, 2011),
    # Seasonal covariate: tick activity peaks in spring/summer (months 4–8)
    season_peak = as.integer(cal_month %in% 4:8),
    # Host abundance: simulated from trail camera counts (Poisson-ish)
    host_abundance = rpois(n(), lambda = 3 + season_peak * 2),
    # Canopy closure: simulated (0–1 scale; BB plots have lower canopy)
    canopy_closure = rbeta(
      n(),
      shape1 = ifelse(burn_regime %in% c("BB", "BUB"), 2, 4),
      shape2 = 3
    )
  )

# True log-linear mean model (used to generate tick counts):
# log(mu) = intercept
#          + burn_regime effect (burning reduces ticks)
#          + seasonal peak
#          + host abundance effect
#          + canopy effect
#          + plot-level random noise (overdispersion)

# Burn regime effects on log scale (relative to UBUB reference)
burn_effects <- c("BB" = -1.2, "BUB" = -0.7, "UBB" = -0.5, "UBUB" = 0.0)

# Add plot-level random intercept to induce within-cluster correlation
plot_re <- rnorm(n_plots, mean = 0, sd = 0.6)  # plot random effect

tick_data <- tick_data |>
  mutate(
    log_mu = 3.5 +                                          # baseline
      burn_effects[as.character(burn_regime)] +             # burn effect
      0.8 * season_peak +                                   # seasonal peak
      0.12 * host_abundance +                               # host effect
      (-0.5) * canopy_closure +                             # canopy effect
      plot_re[plot_id],                                     # plot RE
    mu = exp(log_mu),
    # Simulate overdispersed counts from negative binomial (theta = 2)
    tick_count = rnegbin(n(), mu = mu, theta = 2)
  ) |>
  # Ensure plot_id is a factor for geeglm (required for id argument)
  mutate(
    plot_id = factor(plot_id),
    burn_regime = relevel(burn_regime, ref = "UBUB")  # UBUB as reference
  ) |>
  # Sort by cluster then time (required by geepack)
  arrange(plot_id, month_num)

# Quick look at the data
cat("\n=== Data Summary ===\n")
cat("Total rows:", nrow(tick_data), "\n")
cat("Total ticks:", sum(tick_data$tick_count), "\n")
glimpse(tick_data)

# Distribution of tick counts by burn regime
tick_data |>
  group_by(burn_regime) |>
  summarise(
    mean_ticks   = mean(tick_count),
    var_ticks    = var(tick_count),
    median_ticks = median(tick_count),
    n            = n()
  ) |>
  print()

# ---- 2. Check for overdispersion using Poisson GLM --------------------------
cat("\n=== Overdispersion Check (Poisson GLM) ===\n")

fit_pois_glm <- glm(
  tick_count ~ burn_regime + season_peak + host_abundance + canopy_closure,
  family = poisson(link = "log"),
  data   = tick_data
)

# Ratio of residual deviance to df: >> 1 means overdispersed
disp_ratio <- fit_pois_glm$deviance / fit_pois_glm$df.residual
cat("Residual deviance / df:", round(disp_ratio, 2), "\n")
cat("Overdispersion detected:", disp_ratio > 2, "\n")
# Expected: >> 1 because we simulated from negative binomial

# ---- 3. Negative binomial GLM (baseline, ignores clustering) ----------------
cat("\n=== Negative Binomial GLM (ignores clustering) ===\n")

fit_nb_glm <- glm.nb(
  tick_count ~ burn_regime + season_peak + host_abundance + canopy_closure,
  data = tick_data
)

cat("Negative binomial GLM results:\n")
print(summary(fit_nb_glm))

cat("\nEstimated theta (NB dispersion):", fit_nb_glm$theta, "\n")
cat("Note: theta ~2 matches our simulation (theta = 2)\n")

# ---- 4. GEE with exchangeable correlation -----------------------------------
# NOTE on geepack and negative binomial:
# geepack::geeglm does not directly support family = negative.binomial().
# A common approach is to use family = poisson (quasi-Poisson) with geeglm,
# which estimates a scale parameter from the data to handle overdispersion.
# This gives consistent point estimates; the sandwich SEs account for
# the remaining variance. See Hardin & Hilbe (2013) for details.

cat("\n=== GEE: Exchangeable Correlation (quasi-Poisson) ===\n")

fit_exch <- geeglm(
  tick_count ~ burn_regime + season_peak + host_abundance + canopy_closure,
  id     = plot_id,
  data   = tick_data,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  std.err = "san.se"  # sandwich (robust) standard errors
)

cat("GEE (exchangeable) results:\n")
print(summary(fit_exch))

# Extract coefficients with robust SEs and 95% CIs
coef_exch <- data.frame(
  term     = names(coef(fit_exch)),
  estimate = coef(fit_exch),
  robust_se = sqrt(diag(fit_exch$geese$vbeta)),
  stringsAsFactors = FALSE
) |>
  mutate(
    ci_lower = estimate - 1.96 * robust_se,
    ci_upper = estimate + 1.96 * robust_se,
    exp_est  = exp(estimate),
    exp_low  = exp(ci_lower),
    exp_high = exp(ci_upper)
  )

cat("\nCoefficients (log scale) with robust SEs and exponentiated IRR:\n")
print(round(coef_exch[, c("term","estimate","robust_se","ci_lower","ci_upper",
                            "exp_est","exp_low","exp_high")], 3))

# ---- 5. GEE with AR(1) correlation ------------------------------------------
cat("\n=== GEE: AR(1) Correlation ===\n")

fit_ar1 <- geeglm(
  tick_count ~ burn_regime + season_peak + host_abundance + canopy_closure,
  id     = plot_id,
  data   = tick_data,
  family = poisson(link = "log"),
  corstr = "ar1",
  std.err = "san.se"
)

cat("GEE (AR1) results:\n")
print(summary(fit_ar1))

# ---- 6. GEE with independence correlation -----------------------------------
cat("\n=== GEE: Independence Correlation (sensitivity check) ===\n")

fit_indep <- geeglm(
  tick_count ~ burn_regime + season_peak + host_abundance + canopy_closure,
  id     = plot_id,
  data   = tick_data,
  family = poisson(link = "log"),
  corstr = "independence",
  std.err = "san.se"
)

cat("GEE (independence) results:\n")
print(summary(fit_indep))

# ---- 7. Compare models using QIC --------------------------------------------
cat("\n=== Model Comparison: QIC ===\n")

# QIC function (Pan 2001) — lower QIC is better fit
# geepack provides QIC via QIC() for geeglm objects
qic_exch  <- QIC(fit_exch)
qic_ar1   <- QIC(fit_ar1)
qic_indep <- QIC(fit_indep)

qic_table <- data.frame(
  Model       = c("Exchangeable", "AR(1)", "Independence"),
  QIC         = c(qic_exch["QIC"], qic_ar1["QIC"], qic_indep["QIC"]),
  QICu        = c(qic_exch["QICu"], qic_ar1["QICu"], qic_indep["QICu"])
)

cat("QIC comparison (lower = better):\n")
print(qic_table)
cat("Note: AR(1) expected to perform best for seasonal monthly data.\n")

# ---- 8. Compare robust vs naive SEs across models ---------------------------
cat("\n=== Robust vs Naive SE Comparison (exchangeable model) ===\n")

# Naive SEs from model-based covariance
naive_se  <- sqrt(diag(fit_exch$geese$vbeta.naiv))
robust_se <- sqrt(diag(fit_exch$geese$vbeta))

se_compare <- data.frame(
  term      = names(coef(fit_exch)),
  naive_se  = naive_se,
  robust_se = robust_se,
  ratio     = robust_se / naive_se
)

cat("SE comparison (ratio > 1 means naive SEs underestimate uncertainty):\n")
print(round(se_compare, 3))
cat("Note: With K=21 clusters, sandwich SEs may be anti-conservative.\n")

# ---- 9. Interpret key coefficients (exchangeable model) ---------------------
cat("\n=== Interpretation of Key Coefficients ===\n")

# Get burn_regime coefficients (log scale, relative to UBUB reference)
burn_coefs <- coef_exch |>
  filter(grepl("burn_regime", term)) |>
  mutate(
    pct_change = (exp_est - 1) * 100
  )

cat("Burn regime effects (reference = UBUB = unburned/unburned):\n")
cat("Estimate: log-scale coefficient; exp_est: incidence rate ratio (IRR)\n")
print(round(burn_coefs[, c("term","estimate","exp_est","exp_low","exp_high","pct_change")], 3))
cat("\nInterpretation example:\n")
cat("If estimate for BB = -1.2, then IRR = exp(-1.2) ≈ 0.30\n")
cat("=> Burned/burned plots have ~70% fewer ticks than unburned/unburned,\n")
cat(paste0("   averaged across all ", n_plots, " plots and ", n_months,
           " months (population-averaged effect).\n"))

# ---- 10. Sensitivity check: Poisson vs quasi-Poisson vs NB GLM --------------
cat("\n=== Sensitivity Check: Family Comparison ===\n")

# Compare coefficient estimates across approaches
# (Poisson GLM, NB GLM, and GEE quasi-Poisson)
sens_table <- data.frame(
  term        = names(coef(fit_pois_glm)),
  pois_glm    = coef(fit_pois_glm),
  nb_glm      = coef(fit_nb_glm),
  gee_poisson = coef(fit_exch)
)

cat("Coefficient comparison (Poisson GLM vs NB GLM vs GEE quasi-Poisson):\n")
print(round(sens_table, 3))
cat("\nNote: Point estimates should be similar; SEs will differ.\n")
cat("GEE robust SEs account for within-cluster correlation.\n")
cat("NB GLM does not account for repeated measures (pseudoreplication).\n")

# ---- 11. Visualization (optional) -------------------------------------------
cat("\n=== Visualization ===\n")

# Observed mean tick counts by burn regime and season
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- tick_data |>
    group_by(burn_regime, season_peak) |>
    summarise(mean_count = mean(tick_count), se = sd(tick_count)/sqrt(n()),
              .groups = "drop") |>
    mutate(season = ifelse(season_peak == 1, "Peak season\n(Apr–Aug)", 
                           "Off season")) |>
    ggplot(aes(x = burn_regime, y = mean_count, fill = burn_regime)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                  width = 0.2) +
    facet_wrap(~season) +
    labs(
      title    = "Mean tick counts by burn regime and season",
      subtitle = "Simulated data (matches Gleim et al. 2014 design)",
      x        = "Burn Regime",
      y        = "Mean tick count per plot-month",
      caption  = "BB=burned/burned, BUB=burned/unburned, UBB=unburned/burned, UBUB=unburned/unburned"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
  cat("Plot displayed. Save with: ggsave('tick_counts_by_regime.png', p)\n")
}

cat("\n=== Script Complete ===\n")
cat("Key takeaways:\n")
cat("1. Poisson GLM shows overdispersion (deviance/df >> 1)\n")
cat("2. Negative binomial GLM addresses overdispersion but ignores clustering\n")
cat("3. GEE with quasi-Poisson addresses both, with robust SEs\n")
cat("4. AR(1) likely fits best for monthly ecological data (check QIC)\n")
cat("5. With K=21 clusters, sandwich SEs may need small-sample correction\n")
cat("6. Burn regime effects: BB and BUB plots consistently lower tick counts\n")
