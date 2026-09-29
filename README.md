# Bayesian Generalized Hyperbolic Secant GLMs

Project for the Bayesian Statistics course — Politecnico di Milano, A.Y. 2025-2026

**Authors:** Nalini Lorenzo, Conte Francesco Emanuele, Hinterwaldner Nicholas, Pellizzari Sofia, Tolledi Simone, Lombardi Sofia

## What this is about

The Generalized Hyperbolic Secant (GHS) is one of the six distributions that make up the Natural
Exponential Family with Quadratic Variance Function (NEF-QVF), alongside the Gaussian, Poisson,
Binomial, Gamma and Negative Binomial — but it is by far the least used in practice. Unlike the
Gaussian, it has full support on the real line *and* a variance that grows quadratically with the
mean, which makes it a natural candidate for continuous, heteroskedastic, heavy-tailed responses
that can take both positive and negative values (ratios, differences, log-errors, and similar
variation measures).

This project implements a Bayesian Generalized Linear Model with GHS response, develops the MCMC
machinery needed to sample from its posterior (a Metropolis-within-Gibbs sampler with a MALA step
for the regression coefficients), and tests it both on simulated data and on a real dataset of
hourly temperature changes at New York City airports.

The full theoretical derivation and the complete set of results are in the project report:
[`Bayesian_report.pdf`](Bayesian_report.pdf). This README covers implementation details that
complement the report rather than repeating it — see below.

## Repository structure

All files sit in the repository root (no subfolders):

```
.
├── ghs_mala.R                # MalaWG(): Metropolis-within-Gibbs sampler, MALA step for
│                              # beta + Metropolis-Hastings step for r
├── ghs_sampler.R              # f_GHS(): GHS density; rghs(): exact sampler (Ratio-of-Uniforms)
├── simulated_data_analysis.R  # fits the model on simulated data with known beta and r,
│                              # to check that the sampler recovers the true parameters
├── temperature_analysis.R     # fits the model on the NYC airports temperature dataset,
│                              # with full diagnostics and a comparison against a Bayesian LM
└── Bayesian_report.pdf
```

`temperature_analysis.R` and `simulated_data_analysis.R` both `source()` the two library files, so
they must stay in the same folder as the analysis scripts.

## Requirements

```r
install.packages(c("pracma", "progress", "coda", "ggplot2", "patchwork", "nycflights13"))
```

- `pracma` — provides `gammaz()`, the complex Gamma function used in the GHS density
- `progress` — progress bar during MCMC sampling
- `coda` — MCMC diagnostics (effective sample size, autocorrelation plots)
- `ggplot2` / `patchwork` — credible interval plots
- `nycflights13` — source of the `weather` dataset used in `temperature_analysis.R`

## How to run

```r
source("simulated_data_analysis.R")  # a few minutes, n = 2000
source("temperature_analysis.R")     # around 11 mins on Colab, n ≈ 26,000
```

## Implementation notes

A few choices in the code that go beyond what's in the report:

- **Adaptive step size for MALA.** The report explains *that* the MALA step size ε is tuned
  during burn-in to hit a target acceptance rate (0.574, per Roberts & Rosenthal, 1998), but not
  the exact update rule. The code uses a recursion of the form:
  `log(ε) ← log(ε) + γ_k · (α̂ − target)`, where `α̂` is the running acceptance rate and
  `γ_k = k^(-0.7)` shrinks with the iteration count `k`. During the simulations, it was impossible to achieve good mixing     using a fixed ε: it either rejected almost every proposal (if too large) or barely moved the chain (if too small).

- **The MALA step for beta never calls the complex Gamma function.** The full GHS log-likelihood
  includes a term with `gammaz()`, expensive to evaluate on a dataset with tens of thousands of
  points. When updating beta with r held fixed, that term (and the rest of the normalizing
  constant) doesn't depend on beta and therefore it can be dropped from `log_lik()`.
  `gammaz()` is only evaluated in the Metropolis-Hastings step for r, where it can't be avoided.
  
- **Beta is sampled in a whitened parametrization.** The report's formula for the MALA proposal
  (Sec. 3.2/4.2) is written directly on beta. The code instead samples a standard multivariate Gaussian z and maps it into beta via `beta = mu_prior + L %*% z`, where L is the Cholesky factor of the prior covariance matrix. This way, the prior gradient in z-space collapses to -z, so the inverse prior covariance matrix — needed to evaluate the prior gradient directly on beta — never has to be formed at all. Only L is computed, once, outside the MCMC loop, and reused at every iteration to map z into beta.


## Key results

- On simulated data, the sampler recovers the true beta and r within their 95% posterior credible
  intervals, with acceptance rates and effective sample sizes consistent with healthy MALA mixing.
- On the NYC temperature dataset, the GHS GLM outperforms a Bayesian Gaussian linear model in
  terms of WAIC, consistent with the heteroskedasticity and heavy tails detected in exploratory
  analysis.

## Limitations

The GHS GLM only pays off when the data's variance actually grows quadratically with the mean, the
specific relationship the GHS imposes. When that assumption doesn't hold, a standard Bayesian
Gaussian linear model tends to be more robust.
