# Expanded Lecture Notes: GEE for Repeated-Measures Count Data
## SNR 690 — Tick Phenology & Prescribed Burning (Gleim et al. 2014)

---

## Overview

These notes expand on the slide deck (`slides/slides-ticks-gee.qmd`) and provide:
- Short derivation of the GEE estimating equation from quasi-likelihood
- Sandwich variance derivation (bread and meat)
- Working correlation comparison (pros, cons, guidance)
- Why negative binomial over Poisson for overdispersed counts
- Connection to the specific design and analysis in Gleim et al. (2014)
- Key references

---

## Slide 1 (Title) — Speaker Notes

Introduce the paper briefly: it is an applied ecology paper asking whether long-term prescribed burning reduces tick abundance. The paper is notable for using GEE to account for repeated sampling within plots — a correct and defensible choice, though the methods reporting could be more transparent.

---

## Slide 2 (Agenda) — Speaker Notes

Orient students to the structure of the 75-minute class. Emphasize that the student presenter leads the first half; the second half is instructor-led methods and group activity. The minute paper at the end feeds into next class.

---

## Slide 3 (Entry Ticket Debrief) — Speaker Notes

Fill in the bullet points based on what you actually saw in entry tickets. Key corrections to make:
- Cluster = plot (not individual ticks)
- Outcome = count data → appropriate distribution choices are Poisson or negative binomial
- GEE handles the within-plot repeated measures structure

---

## Slide 4 (Why GEE?) — Extended Notes

### Marginal vs. Conditional Models

The distinction between **marginal** and **conditional** models is foundational.

**Conditional model (GLMM):**
$$
\log(E[Y_{it} \mid b_i]) = \mathbf{X}_{it}^\top \boldsymbol{\beta} + b_i
$$
where $b_i \sim N(0, \sigma^2)$ is a random effect for plot $i$.

The interpretation of $\boldsymbol{\beta}$ here is: *the effect of a one-unit change in a covariate on the log-expected tick count for a plot with random effect $b_i$*. The coefficient is **plot-specific**.

**Marginal model (GEE):**
$$
\log(E[Y_{it}]) = \mathbf{X}_{it}^\top \boldsymbol{\beta}
$$

Here $\boldsymbol{\beta}$ is the **population-averaged** effect — the average effect across all plots.

### Why This Matters Practically

For the tick burning question: "Does prescribed burning reduce tick abundance on average across the landscape?" — the marginal GEE answer is directly what we want. If we asked "what is the tick reduction in plot 7 specifically?", GLMM would be more appropriate.

In large-cohort observational studies (where GEE originated in biostatistics), the population-averaged effect is almost always the policy-relevant quantity.

---

## Slide 5 (GEE Estimating Equation) — Derivation from Quasi-Likelihood

### Setup

Let:
- $K$ = number of clusters (21 plots)
- $n_i$ = number of observations in cluster $i$ (up to 24 months)
- $Y_{it}$ = tick count for plot $i$ at month $t$
- $\mu_{it} = E[Y_{it} \mid \mathbf{X}_{it}]$ = marginal mean
- $\mathbf{Y}_i = (Y_{i1}, \ldots, Y_{in_i})^\top$, $\boldsymbol{\mu}_i = (\mu_{i1}, \ldots, \mu_{in_i})^\top$

### From Score Equations

For a generalized linear model (GLM) with response vector $\mathbf{Y}_i$ and mean $\boldsymbol{\mu}_i$, the quasi-likelihood score equation for a single cluster would be:

$$
\mathbf{S}_i(\boldsymbol{\beta}) = \mathbf{D}_i^\top \mathbf{V}_i^{-1} (\mathbf{Y}_i - \boldsymbol{\mu}_i)
$$

where:
- $\mathbf{D}_i = \frac{\partial \boldsymbol{\mu}_i}{\partial \boldsymbol{\beta}}$ is the $n_i \times p$ Jacobian of the mean w.r.t. parameters
- $\mathbf{V}_i$ is the working covariance matrix of $\mathbf{Y}_i$

Summing over all clusters and setting equal to zero gives the GEE:

$$
\sum_{i=1}^{K} \mathbf{D}_i^\top \mathbf{V}_i^{-1} (\mathbf{Y}_i - \boldsymbol{\mu}_i) = \mathbf{0}
$$

### The Working Covariance

$$
\mathbf{V}_i = \phi \, \mathbf{A}_i^{1/2} \, R(\alpha) \, \mathbf{A}_i^{1/2}
$$

where:
- $\phi$ = dispersion parameter (estimated from data; for negative binomial, this is the overdispersion)
- $\mathbf{A}_i$ = diagonal matrix with $\text{Var}(\mu_{it})$ on the diagonal (depends on the family; for Poisson: $\mu_{it}$; for negative binomial: $\mu_{it} + \mu_{it}^2/\theta$)
- $R(\alpha)$ = working correlation matrix

### Key Result: Consistency Under Misspecification

Liang & Zeger (1986) showed that if the mean structure $\boldsymbol{\mu}_i = h(\mathbf{X}_i \boldsymbol{\beta})$ is correctly specified, then $\hat{\boldsymbol{\beta}}$ is **consistent** regardless of whether $R(\alpha)$ is the true correlation matrix. This is GEE's defining advantage.

---

## Slide 6 (Sandwich Variance) — Extended Notes

### Derivation

The GEE estimator $\hat{\boldsymbol{\beta}}$ solves $\sum_i \mathbf{S}_i(\boldsymbol{\beta}) = \mathbf{0}$.

By a Taylor expansion:

$$
\hat{\boldsymbol{\beta}} - \boldsymbol{\beta}_0 \approx \left[\sum_i \frac{\partial \mathbf{S}_i}{\partial \boldsymbol{\beta}}\right]^{-1} \sum_i \mathbf{S}_i(\boldsymbol{\beta}_0)
$$

The asymptotic variance of $\hat{\boldsymbol{\beta}}$ is:

$$
\text{Var}(\hat{\boldsymbol{\beta}}) = \mathbf{B}^{-1} \, \text{Var}\!\left(\sum_i \mathbf{S}_i\right) \mathbf{B}^{-1}
$$

where $\mathbf{B} = \sum_i \mathbf{D}_i^\top \mathbf{V}_i^{-1} \mathbf{D}_i$ (the **bread**).

Since the clusters are independent, $\text{Var}(\sum_i \mathbf{S}_i) = \sum_i \text{Var}(\mathbf{S}_i)$.

Estimating $\text{Var}(\mathbf{S}_i)$ empirically gives the **meat**:

$$
\mathbf{M} = \sum_i \mathbf{D}_i^\top \mathbf{V}_i^{-1} \widehat{\text{Cov}}(\mathbf{Y}_i) \mathbf{V}_i^{-1} \mathbf{D}_i
$$

where $\widehat{\text{Cov}}(\mathbf{Y}_i) = (\mathbf{Y}_i - \hat{\boldsymbol{\mu}}_i)(\mathbf{Y}_i - \hat{\boldsymbol{\mu}}_i)^\top$.

The **sandwich estimator** is:

$$
\widehat{\text{Var}}(\hat{\boldsymbol{\beta}}) = \mathbf{B}^{-1} \mathbf{M} \mathbf{B}^{-1}
$$

### Small-Sample Warning

The sandwich estimator assumes $K \to \infty$. With $K = 21$ clusters (as in this paper), it can be anti-conservative. Corrections include:
- Jackknife SE: `std.err = "jack"` in `geepack::geeglm`
- Pan (2001) small-sample correction
- Bootstrap SEs (resample entire clusters)

---

## Slide 7 (Paper's Model) — Extended Notes

### What the Paper Reports

Gleim et al. (2014) state (paraphrasing from Methods, p. 3):

> "We used generalized estimating equations (GEE) with a negative binomial distribution and log link function to account for the repeated monthly sampling within plots. The exchangeable working correlation structure was used."

The paper does **not** display a formal model equation, parameter table, or explicit predictor list in a single model specification table. Predictors mentioned in the text include: burn treatment (BB, BUB, UBB, UBUB), host abundance indices (trail camera data), vegetation structure (canopy closure, tree density, ground cover), and climatic variables.

### Inferred Model

A defensible inferred specification is:

$$
\log(E[\text{ticks}_{it}]) = \beta_0 + \beta_1 \cdot \text{BurnRegime}_i + \beta_2 \cdot \text{HostAbundance}_{it} + \beta_3 \cdot \text{Canopy}_{it} + \beta_4 \cdot \text{Season}_{it} + \cdots
$$

Subscript $i$ = plot (cluster), $t$ = month (time point within cluster).

### Interpretation of Coefficients (Log Link, Negative Binomial)

For a log-link model:
- $\exp(\hat{\beta}_1)$ = multiplicative change in expected tick count per unit change in $X_1$
- If $\hat{\beta}_1 = -0.8$ for UBB vs UBUB: $\exp(-0.8) \approx 0.45$ → 55% reduction in tick counts in UBB plots compared to UBUB plots

---

## Slide 8 (Working Correlation) — Extended Notes

### Exchangeable Correlation

$$
R(\alpha) = \begin{pmatrix} 1 & \alpha & \alpha & \cdots \\ \alpha & 1 & \alpha & \cdots \\ \vdots & & \ddots & \\ \alpha & \cdots & \alpha & 1 \end{pmatrix}
$$

**Pros:** Simple; one parameter to estimate; analogous to random intercept model; well-behaved with small $K$.

**Cons:** Assumes January–February correlation = January–December correlation. Unrealistic for time series.

### AR(1) Correlation

$$
R_{ts} = \alpha^{|t-s|}
$$

**Pros:** Ecologically realistic for monthly data (closer months more correlated); one parameter.

**Cons:** Requires equally spaced time points; can behave oddly with missing data.

### Independence

$$
R = I
$$

**Pros:** Simplest; robust with large $K$; equivalent to ignoring clustering (but using sandwich SEs).

**Cons:** Least efficient (widest CIs if true correlation is non-zero).

### Model Selection via QIC

The QIC (Pan 2001) is used to compare GEE models. It is analogous to AIC but accounts for the quasi-likelihood:

$$
\text{QIC}(R) = -2 \hat{Q}(\hat{\boldsymbol{\beta}}; I) + 2 \text{tr}(\hat{\boldsymbol{\Omega}}_I \hat{\mathbf{V}}_R)
$$

where the independence model quasi-likelihood is used as the reference, and $\hat{\mathbf{V}}_R$ is the sandwich variance under working correlation $R$.

---

## Why Negative Binomial Over Poisson for Overdispersed Counts

### The Poisson Constraint

The Poisson distribution assumes $\text{Var}(Y) = E[Y] = \mu$.

For tick count data, this is almost never true:
- Spatial heterogeneity in plots → some plots consistently have very high or low tick counts
- Temporal aggregation → tick counts spike during peak season
- Biological clustering → ticks aggregate in microhabitats

### The Negative Binomial Alternative

The negative binomial distribution with parameters $\mu$ and $\theta$ has:
$$
\text{Var}(Y) = \mu + \frac{\mu^2}{\theta}
$$

As $\theta \to \infty$, NB → Poisson. For small $\theta$ (high overdispersion), variance much exceeds the mean.

### Diagnosing Overdispersion

```r
# Fit Poisson GLM
fit_pois <- glm(tick_count ~ burn_regime + month, 
                family = poisson, data = tick_data)
# Check: residual deviance / df
fit_pois$deviance / fit_pois$df.residual
# Values >> 1 indicate overdispersion
```

In the Gleim et al. data, 47,185 ticks across 21 plots × 24 months × 4 treatments suggests variance >> mean. Negative binomial is appropriate.

### `geepack` Limitation

`geepack::geeglm` does **not** directly support `family = negative.binomial(...)` in the same way as `glm`. Common workarounds:
1. Use `family = poisson` with the scale parameter estimated from data (quasi-Poisson approach)
2. Fit negative binomial GLM (not GEE) with `MASS::glm.nb` and use manual sandwich SEs
3. Use the `gee` package (which does accept `family = "negative.binomial"`)

See `materials/gee_example.R` for implementation of all three approaches.

---

## Connection to Gleim et al. (2014) Design

### Study Design Summary

| Feature | Detail |
|---------|--------|
| Plots (clusters) | 21 |
| Burn treatments | BB (burned/burned), BUB, UBB, UBUB = 4 levels |
| Sampling duration | 2 years × 12 months = 24 time points per plot |
| Total observations | ~21 × 24 = 504 plot-months (some missing due to weather) |
| Outcome | Tick count per drag-flag sample |
| Primary species | *Amblyomma americanum* (99% of ticks) |

### GEE Application in this Context

- **Cluster = plot**: 21 independent sampling units
- **Within-cluster correlation**: Monthly observations within a plot are correlated (seasonal patterns, habitat effects)
- **Exchangeable assumption**: Debatable — monthly data probably have autocorrelation that decays with lag (AR(1) more realistic)
- **Negative binomial**: Appropriate for overdispersed counts
- **Log link**: Natural for multiplicative effects on counts

---

## Key References

1. **Liang, K.-Y., & Zeger, S. L. (1986).** Longitudinal data analysis using generalized linear models. *Biometrika*, 73(1), 13–22.  
   *The foundational GEE paper. Read this.*

2. **Zeger, S. L., & Liang, K.-Y. (1986).** Longitudinal data analysis for discrete and continuous outcomes. *Biometrics*, 42(1), 121–130.  
   *Companion paper with more focus on the discrete (count) case.*

3. **Hardin, J. W., & Hilbe, J. M. (2013).** *Generalized Estimating Equations* (2nd ed.). Chapman & Hall/CRC.  
   *Best textbook treatment; accessible for applied researchers.*

4. **Pan, W. (2001).** Akaike's information criterion in generalized estimating equations. *Biometrics*, 57(1), 120–125.  
   *Introduces QIC for model selection in GEE.*

5. **Højsgaard, S., Halekoh, U., & Yan, J. (2006).** The R package geepack for generalized estimating equations. *Journal of Statistical Software*, 15(2), 1–11.  
   *R implementation reference.*

6. **Gleim, E. R., et al. (2014).** The phenology of ticks and the effects of long-term prescribed burning on tick population dynamics. *PLOS ONE*, 9(11), e112174.  
   DOI: https://doi.org/10.1371/journal.pone.0112174
