# Formative Quiz: Generalized Estimating Equations (GEE)
**SNR 690 — Contemporary Approaches to Quantitative Methods**  
**Topic:** GEE for Repeated-Measures Count Data (Gleim et al. 2014)

*This quiz is formative — for learning and self-assessment, not for grades. Attempt each question before reading the answer.*

---

## Q1 (Multiple Choice): When to Use GEE vs. GLMM

A researcher has 30 forest plots. Each plot is sampled monthly for 18 months for bird abundance (count data). She wants to know whether habitat type (logged vs. unlogged) affects average bird abundance across all plots.

**Which approach is most appropriate and why?**

A) Repeated-measures ANOVA — handles the repeated monthly observations within plots.  
B) Generalized Linear Mixed Model (GLMM) — provides cluster-specific (conditional) effects of habitat type on each plot.  
C) **Generalized Estimating Equations (GEE)** — provides population-averaged (marginal) effect of habitat type across all plots; directly answers the research question.  
D) Standard Poisson GLM with robust SEs — sufficient because habitat type is a between-plot variable.

**Correct answer: C**

**Explanation:** The research question asks about the *average* effect of habitat type across all plots — a population-averaged (marginal) question. GEE is designed for exactly this: it models the marginal mean E[Y] and produces population-averaged coefficients. GLMM (B) would give a conditional effect — the effect for a specific plot with random effect b_i = 0 — which is a subtly different quantity. For a non-linear link (like log for counts), marginal and conditional effects are numerically different. Standard Poisson GLM (D) ignores the within-plot repeated measures correlation, inflating Type I error.

---

## Q2 (Multiple Choice): What the Sandwich Variance Estimator Does

In a GEE with exchangeable working correlation, the analyst suspects the true within-cluster correlation is actually AR(1), not exchangeable. She reports robust (sandwich) standard errors.

**What is the effect of this misspecification on her inference?**

A) Her β̂ estimates are biased because the working covariance is wrong.  
B) Her β̂ estimates are consistent (unbiased), but her sandwich SEs are invalid because they assume exchangeable correlation is correct.  
C) **Her β̂ estimates are consistent, and her sandwich SEs are valid regardless of the working correlation misspecification** — provided the number of clusters is large enough.  
D) She must re-fit the model with the correct AR(1) correlation to get valid estimates and SEs.

**Correct answer: C**

**Explanation:** This is GEE's central advantage. The estimating equation Σ D_i' V_i⁻¹ (Y_i - μ_i) = 0 gives consistent estimates as long as the mean model (link + family + predictors) is correctly specified — regardless of whether R(α) is right. The sandwich SE Var̂(β̂) = B⁻¹MB⁻¹ is also consistent under misspecification because the "meat" M uses empirical residuals, not model-assumed residuals. **The caveat:** this requires K → ∞. With small K (e.g., K = 21 in Gleim et al.), the sandwich estimator can be anti-conservative.

---

## Q3 (Multiple Choice): Choosing a Working Correlation Structure

A researcher samples 50 soil plots monthly for 12 months, measuring nitrogen concentration (continuous). She expects that measurements in adjacent months are more correlated than measurements 6 months apart.

**Which working correlation structure is most ecologically reasonable for her data?**

A) Independence — GEE is robust to misspecification anyway, so it doesn't matter.  
B) **AR(1)** — assumes correlation decays with lag; appropriate when adjacent observations are more correlated than distant ones.  
C) Exchangeable — assumes equal correlation between all pairs of time points; appropriate for clustered data without temporal ordering.  
D) Unstructured — estimates all pairwise correlations separately; always most accurate.

**Correct answer: B**

**Explanation:** AR(1) specifies Corr(Y_t, Y_s) = α^|t-s|, which decays as the lag |t-s| increases. This directly matches the researcher's expectation that adjacent months are more correlated. Exchangeable (C) assumes equal correlation at all lags — unrealistic for time-series data with seasonal patterns. Independence (A) ignores correlation entirely; while GEE estimates remain consistent, efficiency is reduced. Unstructured (D) is flexible but requires estimating n(n-1)/2 correlation parameters — with 12 time points that is 66 parameters, which may be unstable and is not available in all GEE implementations.

---

## Q4 (Short Answer): Write the GEE Estimating Equation and Define V_i

**Prompt:** Write the GEE estimating equation and define the working covariance matrix V_i. Identify each component.

**Answer:**

The GEE estimating equation is:

$$\sum_{i=1}^{K} \mathbf{D}_i^\top \mathbf{V}_i^{-1}(\mathbf{Y}_i - \boldsymbol{\mu}_i) = \mathbf{0}$$

where:

| Symbol | Definition |
|--------|-----------|
| K | Number of clusters (independent sampling units) |
| Y_i | Observed response vector for cluster i (length n_i) |
| μ_i | Vector of marginal means: μ_it = E[Y_it \| X_it] = h(X_it β) |
| D_i | Jacobian matrix: D_i = ∂μ_i/∂β (n_i × p matrix) |
| V_i | Working covariance matrix (see below) |

**The working covariance matrix:**

$$\mathbf{V}_i = \phi \, \mathbf{A}_i^{1/2} \, R(\alpha) \, \mathbf{A}_i^{1/2}$$

- φ: dispersion parameter (estimated from data)
- A_i: diagonal matrix with Var(μ_it) on the diagonal (depends on family; for Poisson: μ_it; for negative binomial: μ_it + μ_it²/θ)
- R(α): working correlation matrix (exchangeable, AR(1), or independence)

**Key property:** GEE estimates β̂ are consistent if the mean structure is correctly specified, regardless of whether R(α) is the true correlation matrix.

---

## Q5 (Short Answer): Interpret a GEE Coefficient from Negative Binomial Model with Log Link

**Prompt:** A GEE model (negative binomial family, log link, exchangeable correlation) is fit to monthly tick count data from 21 plots in four burn regimes. The output shows:

```
                   Estimate  Std.err  Wald   Pr(>|W|)
(Intercept)          3.451    0.312  122.1   <0.001
burn_regimeBUB      -0.693    0.178   15.2   <0.001
burn_regimeBB       -1.204    0.201   35.9   <0.001
burn_regimeUBB      -0.511    0.163    9.8    0.002
season_peak          0.875    0.134   42.6   <0.001
```

*Reference category: UBUB (unburned/unburned). Standard errors are robust (sandwich).*

**Interpret the coefficient for `burn_regimeBB` in plain language.**

**Answer:**

The coefficient for `burn_regimeBB` is −1.204 on the log scale.

**Step 1:** Exponentiate to get the incidence rate ratio (IRR):  
exp(−1.204) ≈ 0.30

**Step 2:** Interpret as a multiplicative change relative to the reference:  
Burned/burned (BB) plots have, on average, approximately **70% fewer ticks** than unburned/unburned (UBUB) plots (IRR = 0.30; robust 95% CI: approximately 0.20–0.45), holding season and other covariates constant.

**Step 3:** Emphasize the marginal interpretation (because this is GEE):  
This is a **population-averaged effect**: averaged across all 21 plots and all time points in the study, BB plots have 70% fewer ticks than UBUB plots. This is not the effect for any specific plot.

**Step 4:** The robust (sandwich) SEs account for the within-plot temporal correlation, giving valid inference even if the exchangeable working correlation is misspecified.

---

*Aligned to Learning Objectives: Q1→LO1, Q2→LO2+LO5, Q3→LO3, Q4→LO2, Q5→LO5*
