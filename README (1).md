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
source("temperature_analysis.R")     # several minutes, n ≈ 26,000
```

## Implementation notes

A few choices in the code that go beyond what's in the report:

- **Adaptive step size for MALA.** The report explains *that* the MALA step size ε is tuned
  during burn-in to hit a target acceptance rate (0.574, per Roberts & Rosenthal, 1998), but not
  the exact update rule. The code uses a Robbins-Monro-type recursion:
  `log(ε) ← log(ε) + γ_k · (α̂ − target)`, where `α̂` is the running acceptance rate and
  `γ_k = k^(-0.7)` shrinks with the iteration count `k`. This adaptive tuning is what makes MALA
  practical here: a fixed, badly-chosen ε either rejects almost every proposal (too large) or
  barely moves the chain (too small) — both kill mixing, which is exactly what the acceptance-rate
  target is designed to avoid.

- **The MALA step for beta never calls the complex Gamma function.** The full GHS log-likelihood
  includes a term with `gammaz()`, expensive to evaluate on a dataset with tens of thousands of
  points. When updating beta with r held fixed, that term (and the rest of the normalizing
  constant) doesn't depend on beta, so it cancels in the MALA acceptance ratio and is dropped from
  `log_lik()`. `gammaz()` is only evaluated in the Metropolis-Hastings step for r, where it can't
  be avoided.

- **Beta is sampled in a whitened parametrization.** The report's formula for the MALA proposal
  (Sec. 3.2/4.2) is written directly on beta with an isotropic proposal. The code instead samples
  `z`, related to beta by `beta = mu_prior + L %*% z` where `L` is the Cholesky factor of the prior
  covariance (`S_prior`) — so a standard normal prior on `z` corresponds exactly to the actual
  Gaussian prior on beta. This whitening (not a matrix inversion — `chol()` here builds the
  transformation, it isn't used to invert anything) lets a single scalar step size work well even
  when the prior covariance is far from diagonal, since the proposal automatically respects the
  correlation structure between coefficients instead of proposing each one independently.

## Key results

- On simulated data, the sampler recovers the true beta and r within their 95% posterior credible
  intervals, with acceptance rates and effective sample sizes consistent with healthy MALA mixing.
- On the NYC temperature dataset, the GHS GLM outperforms a Bayesian Gaussian linear model in
  terms of WAIC, consistent with the heteroskedasticity and heavy tails detected in exploratory
  analysis.
- Credible intervals and Bayes factors mostly agree on which covariates matter, but not always:
  for a couple of coefficients the 95% credible interval includes zero while the Bayes factor
  still favours inclusion strongly. This isn't a bug — it's a known consequence of comparing a
  point-null Bayes factor against an interval estimate under a fairly diffuse prior (related to
  the Jeffreys-Lindley phenomenon), and it's a good reminder that the two criteria answer
  different questions rather than being interchangeable.

## Limitations

The GHS GLM only pays off when the data's variance actually grows quadratically with the mean, the
specific relationship the GHS imposes. When that assumption doesn't hold, a standard Bayesian
Gaussian linear model tends to be more robust — see the report for a worked example (election
prediction errors) where the GHS underperforms for exactly this reason.
