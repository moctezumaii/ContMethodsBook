# Expected Answers: Tick Phenology & GEE Workshop
**SNR 690 — For Instructor Use Only**  
**Paper:** Gleim et al. (2014). *PLOS ONE*, 9(11), e112174. DOI: [10.1371/journal.pone.0112174](https://doi.org/10.1371/journal.pone.0112174)

---

## Part 1: Entry Ticket Expected Answers

| Field | Expected Answer |
|-------|----------------|
| **Objective** | To quantify the effects of long-term prescribed burning regimes on tick population dynamics (phenology, abundance, species composition) in the longleaf pine ecosystem of southwestern Georgia and northwestern Florida. |
| **Research Question** | Does long-term prescribed burning (under different burn regimes: burned/burned, burned/unburned, unburned/burned, unburned/unburned) reduce tick abundance and alter tick species composition, controlling for habitat and host variables? |
| **Outcome Variable** | Tick count per drag-flag sample (per plot per month). Overdispersed count data (47,185 total ticks; 99% *Amblyomma americanum*). |
| **Predictors** | Burn treatment regime (BB, BUB, UBB, UBUB); host abundance (trail camera index); vegetation structure (canopy closure, tree density, ground cover); climatic variables (temperature, precipitation). |
| **Cluster ID and time structure** | Cluster = sampling plot (21 plots; 7.2–78.9 ha). Time = monthly sampling over 2 years (January 2010–December 2011) = 24 time points per plot. Plots are the independent units; monthly observations within a plot are the repeated measures. |
| **Analysis method guess** | Generalized Estimating Equations (GEE) with negative binomial family, log link, and exchangeable working correlation. (Also acceptable: repeated-measures GLM, mixed model.) |
| **Method question** | Accepts any genuine methodological question. Examples of strong questions: "Why exchangeable and not AR(1)?"; "How do you handle overdispersion in GEE?"; "Why GEE instead of a GLMM with random plot effects?" |

---

## Part 2: Worksheet Critique Prompts — Detailed Expected Answers

### Prompt 1: Is negative binomial appropriate? What alternatives exist?

**Expected answer:**

Yes, negative binomial is appropriate. Tick count data are overdispersed: the variance of counts far exceeds the mean due to spatial aggregation of ticks in microhabitats, seasonal patchiness, and among-plot heterogeneity. Overdispersion in Poisson models causes standard errors to be underestimated (anti-conservative inference). The negative binomial adds a dispersion parameter $\theta$, allowing Var(Y) = μ + μ²/θ.

**Alternatives and when to prefer them:**
- **Poisson:** Appropriate when Var(Y) ≈ E(Y). Should be rejected here (overdispersion ratio >> 1 expected). 
- **Quasi-Poisson:** Estimates a scale parameter from the data; equivalent to Poisson with scaled SEs. Pragmatic and widely used in GEE implementations (geepack uses this approach).
- **Zero-inflated NB:** Appropriate when there are excess zeros (e.g., months with zero ticks in winter). Check zero-inflation separately; the paper does not mention it.

**Key diagnostic:** Compute `deviance / df` from a Poisson GLM. Values >> 1 indicate overdispersion.

---

### Prompt 2: Is exchangeable correlation realistic for monthly ecological data? AR(1)?

**Expected answer:**

**Exchangeable** assumes constant pairwise correlation regardless of temporal lag: Corr(Y₁, Y₂) = Corr(Y₁, Y₂₄). For monthly tick counts over 2 years, this is likely unrealistic. Tick counts within a plot are probably more correlated between adjacent months (due to seasonal phenology, carryover effects) than between months a year apart.

**AR(1)** assumes Corr(Y_t, Y_s) = α^|t-s|, which decays with lag. This is more ecologically realistic for monthly data.

**Practical consequence of misspecification:** GEE point estimates remain consistent (unbiased for large K) regardless of working correlation choice. However, efficiency (precision of estimates) depends on how close R(α) is to the truth. Using exchangeable when AR(1) is true gives valid but somewhat less efficient estimates. **The bigger concern here is K = 21 clusters: the sandwich SE may be unreliable regardless of correlation structure choice.**

**Model selection:** Compare QIC under exchangeable vs. AR(1) — lower QIC favors better-fitting correlation structure.

---

### Prompt 3: Are cluster sizes balanced? How might imbalance affect GEE?

**Expected answer:**

**Balance across treatment groups:** 21 plots across 4 burn regimes does not divide evenly. Approximately 5–6 plots per treatment. This is mild imbalance; GEE handles unequal cluster sizes without special treatment.

**Plot area imbalance:** Plots range from 7.2 to 78.9 ha. Larger plots likely capture more tick drag samples, affecting raw counts. This could be addressed by:
1. Including plot area as a covariate (log area)
2. Using tick density (count per unit area) as the outcome
3. Including log(area) as an **offset** in the model: `log(E[ticks]) = log(area) + Xβ`

The paper does not appear to use an offset — this is a legitimate limitation.

**Temporal imbalance:** Some plot-months may be missing (bad weather, no sampling). GEE handles missing time points within clusters without requiring complete data (unlike RM-ANOVA).

---

### Prompt 4: Transparency — the paper doesn't report the full model equation

**Expected answer:**

The paper's methods section describes using GEE with negative binomial family, log link, and exchangeable correlation, but does not display a formal model equation, a complete predictor list in equation form, or parameter estimates in a standard table for the GEE model.

**Is this a transparency issue?** Yes. For full reproducibility, a paper should provide:
1. Formal model equation with notation
2. Complete list of predictors with coding (reference category for factors)
3. Software and package used (`geepack`, `gee`, `SAS PROC GENMOD`)
4. Parameter estimates, robust SEs, and 95% CIs for all predictors
5. Information on handling missing data

**This is common in applied ecology** but does not meet modern reproducibility standards (e.g., PLOS ONE data availability policy, Transparent Statistics in Ecology guidelines).

**Strong reviewer comment:** "The authors should provide a formal model specification with predictor notation and a complete parameter table with robust SEs."

---

### Prompt 5: Sandwich variance — how does it protect against misspecification?

**Expected answer:**

The sandwich (robust) variance estimator Var̂(β̂) = B⁻¹MB⁻¹ separates the model-based information matrix B (the "bread") from the empirical variance M (the "meat").

- **B** assumes the working covariance V_i is correctly specified.
- **M** uses the actual observed residuals: M = Σᵢ Dᵢ' Vᵢ⁻¹ Ĉov(Yᵢ) Vᵢ⁻¹ Dᵢ.

If V_i is wrong, B⁻¹ alone gives incorrect SEs. But the empirical M captures the true residual variability. The sandwich combines them, giving consistent SEs even if R(α) is wrong.

**Limitation for K = 21:** The sandwich estimator is asymptotically consistent (requires K → ∞). With only 21 clusters, it can be anti-conservative (SEs too small). Corrections: jackknife SEs (`std.err = "jack"` in geepack) or bootstrap resampling at the cluster level.

---

### Prompt 6: GEE vs GLMM — marginal vs conditional interpretation

**Expected answer:**

**GEE (marginal/population-averaged):**
The burn regime coefficient represents the average effect of burning across *all* 21 plots.
Interpretation: "Burned/burned plots have, on average, 70% fewer ticks than unburned/unburned plots, across the entire sampled landscape."

**GLMM (conditional/cluster-specific):**
The burn regime coefficient represents the effect of burning for a plot with random effect b_i = 0 (the "average" plot).
Interpretation: "For a plot with average unobserved characteristics, burning reduces tick counts by 70%."

For a non-linear link (like log), these are numerically different: the marginal effect (GEE) is smaller than the conditional effect (GLMM) due to Jensen's inequality.

**Which is more relevant for land management?**
The marginal GEE interpretation — "on average, burning reduces ticks by X% across the landscape" — is directly what a land manager needs to make a policy decision. GLMM answers a subtly different, plot-specific question.

---

## Part 3: Sample GEE Output Interpretation

Based on the simulated data in `materials/gee_example.R`:

```
Burn regime effects (reference = UBUB = unburned/unburned):
         estimate  robust_se   IRR    IRR_low  IRR_high  pct_change
burn_BB    -1.20    0.18     0.30     0.21     0.43       -70%
burn_BUB   -0.70    0.14     0.50     0.38     0.66       -50%
burn_UBB   -0.50    0.12     0.61     0.48     0.77       -39%
```

**Plain-language interpretation:**

- **BB (burned/burned):** Plots that received prescribed burns in both the first and second burn cycle had approximately 70% fewer ticks (IRR = 0.30, 95% CI: 0.21–0.43) than unburned/unburned plots, averaged across all plots and time points.
- **BUB and UBB:** Intermediate reductions (~40–50%) in plots burned in only one of the two cycles.
- **Seasonal peak coefficient:** During peak tick season (April–August), expected tick counts were approximately 2–3× higher than in off-season months (IRR ≈ 2.2–2.5).

---

## Part 4: Suggested Reviewer Comments (for rubric use)

### Strong reviewer comments (Score 4 across criteria):

- "The presenter correctly identified the GEE model structure and articulated the marginal vs conditional distinction clearly. The IR interpretation (70% fewer ticks for BB plots) was ecologically accurate and connected to management implications."
- "The R code correctly implements geeglm with exchangeable correlation. The comparison of QIC across correlation structures was a strong addition."
- "The presenter proactively noted that the paper does not provide a formal equation and proposed a defensible inferred specification."

### Moderate comments (Score 2–3):

- "The model family and link were identified correctly, but the presenter did not explain why negative binomial is preferred over Poisson for these data. Consider adding overdispersion diagnostics."
- "Results were summarized as 'significant' but coefficients were not interpreted on the original scale (exp(β) = IRR). Clarify the multiplicative interpretation."
- "The code ran correctly, but did not compare across correlation structures or test for model fit with QIC."

### Lowest-level comments (Score 1):

- "The presenter identified the paper as using 'ANOVA' — this is incorrect. The method is GEE (generalized estimating equations), a marginal model for correlated count data."
- "The cluster structure was not identified correctly. Individual tick observations are not the analysis unit; sampling plots (n=21) are the clusters."
