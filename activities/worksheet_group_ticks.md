# Group Worksheet: Critiquing a GEE Analysis
**SNR 690 — Contemporary Approaches to Quantitative Methods**  
**Paper:** Gleim et al. (2014). Tick phenology & prescribed burning. *PLOS ONE*, 9(11), e112174. DOI: [10.1371/journal.pone.0112174](https://doi.org/10.1371/journal.pone.0112174)  
**Time limit:** 10 minutes | **Group size:** 3–4 students

---

## The Paper's Model Specification (Inferred)

> **Note:** Gleim et al. (2014) do **not** provide an explicit model equation. The specification below is inferred from the Methods section.

| Element | Specification |
|---------|--------------|
| **Family** | Negative binomial (overdispersed counts) |
| **Link** | Log |
| **Correlation structure** | Exchangeable |
| **Clusters** | 21 sampling plots |
| **Repeated measures** | Monthly sampling over 2 years (up to 24 time points per plot) |
| **Predictor(s)** | Burn regime (BB, BUB, UBB, UBUB), host abundance, canopy closure, tree density, ground cover, temperature, precipitation |

**Inferred equation:**  
log(E[ticks_it]) = β₀ + β₁·BurnRegime_i + β₂·HostAbundance_it + β₃·Canopy_it + β₄·Season_it + …

---

## Critique Prompts

**Answer as many as you can in 10 minutes. Your group's main deliverable is prompt 7 (the challenge question).**

**1. Family choice: Is negative binomial appropriate?**  
What is overdispersion, and how would you check for it? What alternatives exist (Poisson, quasi-Poisson, zero-inflated)? Under what conditions would you prefer each?

*Your group's answer:*  
&nbsp;  
&nbsp;  

---

**2. Correlation structure: Is exchangeable realistic for monthly ecological data?**  
For a 24-month time series, exchangeable correlation assumes Corr(Jan, Feb) = Corr(Jan, Dec). Is that biologically realistic? Would AR(1) be more defensible? What is the practical consequence of misspecifying the correlation structure?

*Your group's answer:*  
&nbsp;  
&nbsp;  

---

**3. Cluster balance: Are the 21 plots balanced across treatments?**  
21 plots across 4 burn regimes does not divide evenly. How might unequal cluster sizes affect GEE? Does unequal plot *area* (7.2–78.9 ha) matter? Should plot area be included as a covariate or offset?

*Your group's answer:*  
&nbsp;  
&nbsp;  

---

**4. Transparency: The paper does not report the full model equation — is this a problem?**  
What information would you need to fully replicate the analysis? Does the absence of a formal model equation meet the standards for reproducible research? What would you require in a peer review?

*Your group's answer:*  
&nbsp;  
&nbsp;  

---

**5. Robust SEs: How does the sandwich variance estimator protect against correlation misspecification?**  
Explain in plain language why GEE point estimates are consistent even if the working correlation is wrong. Under what conditions does the sandwich estimator fail? Is K = 21 clusters sufficient?

*Your group's answer:*  
&nbsp;  
&nbsp;  

---

**6. GEE vs. GLMM: What would change if they used a GLMM instead?**  
Describe one practical difference in how you would interpret the burn regime coefficient from a GEE vs. a GLMM. Which interpretation is more relevant to a land manager deciding whether to prescribe burn?

*Your group's answer:*  
&nbsp;  
&nbsp;  

---

## ⚡ Deliverable: One Challenge Question for the Presenter

Write one methodological question that you want the presenter to answer. This should be a question that requires genuine engagement with the paper's methods — not a factual recall question.

> **Our group's challenge question:**  
&nbsp;  
&nbsp;  
&nbsp;  

---

*Learning objectives addressed: LO 1, 2, 3, 4 (critique), 5 (implementation)*  
*Instructor: Collect challenge questions at 10 min; call on 2–3 groups during group reports (1:05–1:12)*
