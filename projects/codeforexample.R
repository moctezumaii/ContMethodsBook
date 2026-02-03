# ============================================================================
# Week 03: Method-Question-Data Triangle
# Simulation Examples: Mismatch, Good Match, Overcomplicated
# ============================================================================

# Load packages
library(lme4)       # For mixed models
library(ggplot2)    # For plotting
library(dplyr)      # For data manipulation
library(broom)      # For tidy model output
library(broom.mixed) # For tidy mixed model output

set.seed(2024)  # For reproducibility

# ============================================================================
# EXAMPLE 1: MISMATCH
# Ignoring nested/grouped data structure
# ============================================================================

# Scenario: Testing fertilizer effect on crop yield
# - 5 fields (random effect - different baseline fertility)
# - 4 plots per field (2 control, 2 fertilized)
# - Problem: Plots within fields are NOT independent

# --- Simulate the data ---

n_fields <- 5
n_plots_per_field <- 4

# Create data structure
mismatch_data <- expand.grid(
  field = factor(1:n_fields),
  plot = 1:n_plots_per_field
) %>%
  mutate(
    # Assign treatment: 2 control, 2 fertilized per field
    treatment = rep(c("control", "control", "fertilizer", "fertilizer"), n_fields),
    treatment = factor(treatment, levels = c("control", "fertilizer"))
  )

# Simulate field-level random effects (some fields are just better)
field_effects <- data.frame(
  field = factor(1:n_fields),
  field_effect = rnorm(n_fields, mean = 0, sd = 8)  # Large field variation!
)

# True treatment effect
true_effect <- 3  # Fertilizer adds 3 units (small effect)

# Generate yields
mismatch_data <- mismatch_data %>%
  left_join(field_effects, by = "field") %>%
  mutate(
    # Yield = baseline + field effect + treatment effect + noise
    yield = 50 + 
            field_effect + 
            ifelse(treatment == "fertilizer", true_effect, 0) +
            rnorm(n(), mean = 0, sd = 2)  # Small residual error
  )

# --- The WRONG analysis (ignoring field structure) ---
wrong_model <- lm(yield ~ treatment, data = mismatch_data)

# --- The CORRECT analysis (accounting for field) ---
correct_model <- lmer(yield ~ treatment + (1|field), data = mismatch_data)

# --- Compare results ---
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("EXAMPLE 1: MISMATCH - Ignoring Nested Structure\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("TRUE treatment effect:", true_effect, "\n\n")

cat("WRONG MODEL (lm - ignores field structure):\n")
summary(wrong_model)$coefficients

cat("\n\nCORRECT MODEL (lmer - accounts for field):\n")
summary(correct_model)$coefficients

# Extract p-values for comparison
wrong_p <- summary(wrong_model)$coefficients["treatmentfertilizer", "Pr(>|t|)"]
correct_t <- summary(correct_model)$coefficients["treatmentfertilizer", "t value"]
# Approximate p-value for mixed model (Satterthwaite or similar would be better)
correct_p <- 2 * (1 - pnorm(abs(correct_t)))

cat("\n\nP-VALUE COMPARISON:\n")
cat("Wrong model p-value:   ", round(wrong_p, 4), "\n")
cat("Correct model p-value: ", round(correct_p, 4), "\n")
cat("\nThe wrong model may show 'significant' results that aren't real!\n")

# --- Visualization ---
p1 <- ggplot(mismatch_data, aes(x = treatment, y = yield, color = field)) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  geom_line(aes(group = field), alpha = 0.3) +
  stat_summary(aes(group = 1), fun = mean, geom = "crossbar", 
               width = 0.5, color = "black", linewidth = 1) +
  labs(
    title = "Example 1: MISMATCH",
    subtitle = "Ignoring field structure inflates significance",
    x = "Treatment",
    y = "Yield",
    color = "Field"
  ) +
  theme_minimal(base_size = 14)

print(p1)

# ============================================================================
# EXAMPLE 2: GOOD MATCH
# Count data with appropriate Poisson GLM
# ============================================================================

# Scenario: Pollinator visits to flowers
# - 2 treatments: native vs non-native plants
# - 30 plants per treatment
# - Response: count of pollinator visits (0 to ~45)

# --- Simulate the data ---

n_per_group <- 30

goodmatch_data <- data.frame(
  plant_id = 1:(2 * n_per_group),
  treatment = factor(rep(c("native", "non_native"), each = n_per_group),
                     levels = c("non_native", "native"))
)

# True effect: native plants get more visits
baseline_visits <- 8      # non-native mean
native_effect <- 0.5       # log-scale effect (multiplicative)

# Simulate counts from Poisson
goodmatch_data <- goodmatch_data %>%
  mutate(
    log_mu = log(baseline_visits) + ifelse(treatment == "native", native_effect, 0),
    visits = rpois(n(), lambda = exp(log_mu))
  )

# --- The CORRECT analysis (Poisson GLM) ---
poisson_model <- glm(visits ~ treatment, family = poisson, data = goodmatch_data)

# --- For comparison: what if we used lm? ---
wrong_lm <- lm(visits ~ treatment, data = goodmatch_data)

# --- Results ---
cat("\n\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("EXAMPLE 2: GOOD MATCH - Count Data with Poisson GLM\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("TRUE effect (log scale):", native_effect, "\n")
cat("TRUE effect (multiplicative):", round(exp(native_effect), 2), "x more visits\n\n")

cat("POISSON GLM (correct):\n")
summary(poisson_model)$coefficients

cat("\n\nInterpretation:\n")
est <- coef(poisson_model)["treatmentnative"]
cat("Native plants receive", round(exp(est), 2), "times as many visits\n")
cat("That's", round((exp(est) - 1) * 100, 1), "% more visits\n")

# --- Visualization ---
p2 <- ggplot(goodmatch_data, aes(x = treatment, y = visits, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("non_native" = "#E69F00", "native" = "#56B4E9")) +
  labs(
    title = "Example 2: GOOD MATCH",
    subtitle = "Count data → Poisson GLM",
    x = "Plant Type",
    y = "Pollinator Visits"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(p2)

# Check model assumptions with a quick mean-variance check
cat("\n\nModel check - Mean vs Variance by group:\n")
goodmatch_data %>%
  group_by(treatment) %>%
  summarise(
    mean_visits = mean(visits),
    var_visits = var(visits),
    ratio = var_visits / mean_visits
  ) %>%
  print()
cat("(Ratio near 1 suggests Poisson is appropriate)\n")

# ============================================================================
# EXAMPLE 3: OVERCOMPLICATED
# Simple comparison done with unnecessarily complex model
# ============================================================================

# Scenario: Tree height in 2 forest types
# - 30 trees per forest type
# - Response: height in meters (continuous, normal)
# - Simple question: Is there a difference?

# --- Simulate the data ---

n_trees <- 30

overcomp_data <- data.frame(
  tree_id = 1:(2 * n_trees),
  forest_type = factor(rep(c("deciduous", "coniferous"), each = n_trees))
)

# True difference
deciduous_mean <- 18
coniferous_mean <- 22
tree_sd <- 4

overcomp_data <- overcomp_data %>%
  mutate(
    height = ifelse(forest_type == "deciduous",
                    rnorm(n(), deciduous_mean, tree_sd),
                    rnorm(n(), coniferous_mean, tree_sd))
  )

# --- Simple approach (appropriate) ---
simple_ttest <- t.test(height ~ forest_type, data = overcomp_data)
simple_lm <- lm(height ~ forest_type, data = overcomp_data)

# --- Overcomplicated approach (unnecessary) ---
# In practice, you might use brms or Stan for Bayesian model
# Here we simulate what that would give us (same answer, more effort)

# For demonstration, let's use a mixed model with no random effects needed
# (This is silly but illustrates the point)
# We'll also bootstrap for "uncertainty quantification"

bootstrap_diff <- replicate(1000, {
  boot_data <- overcomp_data %>% sample_n(n(), replace = TRUE)
  diff(tapply(boot_data$height, boot_data$forest_type, mean))
})

# --- Results ---
cat("\n\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("EXAMPLE 3: OVERCOMPLICATED - Simple Question, Complex Method\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("TRUE difference:", coniferous_mean - deciduous_mean, "meters\n\n")

cat("SIMPLE T-TEST (appropriate):\n")
cat("  Estimated difference:", round(diff(simple_ttest$estimate), 2), "\n")
cat("  95% CI: [", round(simple_ttest$conf.int[1], 2), ",", 
    round(simple_ttest$conf.int[2], 2), "]\n")
cat("  p-value:", format(simple_ttest$p.value, digits = 3), "\n")
cat("  Time to run: ~0.001 seconds\n\n")

cat("SIMPLE LINEAR MODEL (also appropriate):\n")
print(round(summary(simple_lm)$coefficients, 3))

cat("\n\nOVERCOMPLICATED BOOTSTRAP (unnecessary):\n")
cat("  Estimated difference:", round(mean(bootstrap_diff), 2), "\n")
cat("  95% CI: [", round(quantile(bootstrap_diff, 0.025), 2), ",",
    round(quantile(bootstrap_diff, 0.975), 2), "]\n")
cat("  Time to run: ~1-2 seconds (1000x slower)\n")
cat("  And a full Bayesian model would take minutes!\n\n")

cat("CONCLUSION: Same answer, unnecessary complexity!\n")
cat("Save fancy methods for when you NEED them.\n")

# --- Visualization ---
p3 <- ggplot(overcomp_data, aes(x = forest_type, y = height, fill = forest_type)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = c("deciduous" = "#228B22", "coniferous" = "#006400")) +
  labs(
    title = "Example 3: OVERCOMPLICATED",
    subtitle = "Normal data, simple comparison → t-test is fine!",
    x = "Forest Type",
    y = "Tree Height (m)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(p3)

# ============================================================================
# SUMMARY COMPARISON
# ============================================================================

cat("\n\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SUMMARY: THE THREE SCENARIOS\n
")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

summary_table <- data.frame(
  Example = c("1. Mismatch", "2. Good Match", "3. Overcomplicated"),
  Problem = c("Ignored grouping", "None", "Unnecessary complexity"),
  Consequence = c("False positive risk", "Valid inference", "Wasted effort"),
  Lesson = c("Check independence", "Match data to distribution", "Start simple")
)

print(summary_table, row.names = FALSE)

# ============================================================================
# COMBINED FIGURE FOR SLIDES
# ============================================================================

library(patchwork)

combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Three Scenarios: Mismatch, Good Match, Overcomplicated",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_plot)

# Save for slides
# ggsave("three_scenarios.png", combined_plot, width = 14, height = 5, dpi = 300)