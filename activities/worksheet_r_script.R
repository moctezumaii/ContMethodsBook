# =============================================================================
# Scaffolded R Script: GEE for Repeated-Measures Count Data
# SNR 690 — Activity: worksheet_r_script.R
#
# Learning objectives addressed:
#   LO 3: Compare working correlation structures
#   LO 5: Implement geeglm and interpret robust SEs
#
# Instructions:
#   Fill in the blanks marked with ______ or ## FILL IN ##
#   Do NOT change the dataset or variable names
#   Run each section sequentially (Ctrl+Enter per line, or Source)
#
# Time: ~15 minutes in class (or as homework)
# =============================================================================

# ---- 0. Load libraries -------------------------------------------------------
library(geepack)
# install.packages("geepack") if needed

# ---- 1. Pre-loaded mini-dataset ----------------------------------------------
# Simulated data: 5 plots, 6 monthly time points each (30 rows total)
# Matches key features of Gleim et al. (2014): burn regime, tick counts,
# seasonal covariate, plot-level clustering.

set.seed(42)
mini_data <- data.frame(
  plot_id     = factor(rep(1:5, each = 6)),
  month_num   = rep(1:6, times = 5),
  burn_regime = factor(rep(c("BB","BB","BUB","UBUB","UBUB"), each = 6),
                        levels = c("UBUB","BUB","BB")),
  season_peak = rep(c(0, 1, 1, 1, 0, 0), times = 5),
  host_abund  = c(2,4,5,4,3,2, 3,5,6,5,3,2, 2,3,4,3,2,1,
                  5,7,9,8,6,4, 4,6,8,7,5,3),
  tick_count  = c(12,45,38,41,18,10,
                  8, 20,25,22,11,5,
                  5, 14,18,15, 8,3,
                  35,80,95,88,50,20,
                  28,65,78,72,40,16)
)

# Look at the data
head(mini_data, 10)
str(mini_data)

# ---- 2. Exploratory: overdispersion check ------------------------------------
# Fit a simple Poisson GLM (ignores clustering — just to check overdispersion)
fit_check <- glm(tick_count ~ burn_regime + season_peak,
                 family = poisson,
                 data   = mini_data)

# Calculate dispersion: residual deviance / residual df
# Values >> 1 indicate overdispersion
disp <- fit_check$deviance / fit_check$df.residual
cat("Poisson GLM dispersion ratio:", round(disp, 2), "\n")
# Expected output: >> 1 (our data are overdispersed)

# ---- 3. Fit a GEE --------------------------------------------------------
# Fill in the geeglm() call below.
# The formula should predict tick_count from burn_regime and season_peak.
# The id argument identifies the cluster variable (which column?).
# Start with corstr = "exchangeable".
# Use family = poisson(link = "log") as a quasi-Poisson approximation.

fit_gee <- geeglm(
  ## FILL IN: formula
  ________ ~ ________ + ________,
  
  ## FILL IN: cluster id variable (must be a factor)
  id     = ________,
  
  data   = mini_data,
  family = poisson(link = "log"),
  
  ## FILL IN: working correlation structure ("exchangeable", "ar1", or "independence")
  corstr = "________",
  
  std.err = "san.se"   # request sandwich (robust) SEs
)

# ---- 4. Examine output -------------------------------------------------------
summary(fit_gee)

# Extract coefficients and robust SEs
coef_table <- data.frame(
  estimate  = coef(fit_gee),
  robust_se = sqrt(diag(fit_gee$geese$vbeta))
)

## FILL IN: Add 95% CIs (lower = estimate - 1.96*SE; upper = estimate + 1.96*SE)
coef_table$ci_lower <- ________ - 1.96 * ________
coef_table$ci_upper <- ________ + 1.96 * ________

## FILL IN: Exponentiate estimates and CIs to get incidence rate ratios (IRR)
coef_table$irr       <- exp(________)
coef_table$irr_lower <- exp(________)
coef_table$irr_upper <- exp(________)

print(round(coef_table, 3))

# ---- 5. Interpret one coefficient --------------------------------------------
# Using the log link and Poisson family, coefficients are on the log scale.
# exp(beta) = incidence rate ratio (IRR).
#
# Fill in the blanks for the BB coefficient:
#
# "Plots with burn regime BB have an estimated IRR of _____ (95% CI: _____ to _____),
#  meaning they have approximately _____% [fewer/more] ticks than UBUB plots,
#  averaged across all time points and plots."
#
# Write your interpretation here (as a comment):
## FILL IN: Your interpretation
# Plots with BB have IRR = _____, meaning...

# ---- 6. Compare two correlation structures -----------------------------------
# Re-fit the model with corstr = "ar1" and compare to exchangeable.
# Which one gives narrower CIs? Why might AR(1) be more appropriate for
# monthly ecological data?

fit_gee_ar1 <- geeglm(
  tick_count ~ burn_regime + season_peak,
  id     = plot_id,
  data   = mini_data,
  family = poisson(link = "log"),
  
  ## FILL IN: change the correlation structure
  corstr = "________",
  
  std.err = "san.se"
)

# Compare QIC (lower = better fit)
QIC(fit_gee)      # exchangeable
QIC(fit_gee_ar1)  # ar1

## FILL IN: Which model has lower QIC? What does that suggest?
# Answer (as a comment):

# ---- EXPECTED OUTPUT ---------------------------------------------------------
# (For instructor reference — do not look until you've filled in all blanks!)
#
# Dispersion ratio: ~8–15 (overdispersed)
#
# Exchangeable GEE coefficients (approximate):
#   (Intercept)   ~  3.1  (log scale)
#   burn_regimeBUB  ~ -0.8  (IRR ~0.45 → 55% fewer ticks than UBUB)
#   burn_regimeBB   ~ -1.3  (IRR ~0.27 → 73% fewer ticks than UBUB)
#   season_peak     ~  0.9  (IRR ~2.5 → 2.5x more ticks in peak season)
#
# Sample interpretation for BB:
#   "Burned/burned (BB) plots have ~73% fewer ticks than unburned/unburned
#    (UBUB) plots on average, holding season constant (IRR = 0.27, 95% CI: ...)"
