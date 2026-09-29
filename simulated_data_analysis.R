# ============================================================
# Convergence check: MALA on simulated GHS data
# the code has the following sections:
#
# DATA SIMULATION: simulate a large dataset from a known GHS GLM
#   with a single covariate (2 beta: intercept + slope)
#
# MALA WITHIN GIBBS: samples from posterior via MALA (using functions in
#   'ghs_mala.R' and 'GHS_sampler.R')
#
# DIAGNOSTIC TRACEPLOTS: assess MCMC mixing
#
# PARAMETER RECOVERY: true vs. estimated beta and r, 95% credible intervals
#
# CONVERGENCE PLOT: posterior density of each parameter vs. its true value
#
# BAYESIAN LINEAR MODEL: function to regress data with a semi-conjugate
#   normal-IG linear model
#
# WAIC COMPUTATION: compute WAIC for both the linear model and the GHS
# ============================================================

library(pracma)
library(progress)
library(coda)
library(ggplot2)
source("ghs_mala.R")
source("ghs_sampler.R")

# === DATA SIMULATION ===
n = 2000
r = 0.5
beta = c(-1, 2)
p = length(beta)-1
X = matrix(rnorm(n*p,0,1),n,p)
X = cbind(rep(1,n),X)
theta = atan((X%*%beta)/r)
y = rep(0,n)
for (i in 1:n){
  y[i] = rghs(1,theta[i],r)
}

# === MALA WITHIN GIBBS ===
S = diag(rep(3,p+1))
mu = rep(0,p+1)
n_iter = 1e4
BI = 1e3
res = MalaWG(y,X,n_iter=n_iter,BI=BI,var_sample_r=5e-2,
             mu_prior=mu, S_prior=S,
             log_prior_r = function(x) dgamma(x,r,1,log=TRUE))

colnames(res$beta) = c("beta0", "beta1")

cat("Acceptance rate (beta, MALA):", round(res$acc_rate_beta, 3), "\n")

# === DIAGNOSTIC TRACEPLOTS ===
plot(1:(n_iter-BI),res$r,type="l",ylab="r")
plot(1:(n_iter-BI),res$beta[,1],type="l",ylab="beta0")
plot(1:(n_iter-BI),res$beta[,2],type="l",ylab="beta1")

m <- mcmc(res$beta)
print(acfplot(m, lag.max = 500))  # needs print(): lattice objects aren't auto-drawn via source()

cat("Effective sample size (out of", nrow(res$beta), "kept draws):\n")
print(effectiveSize(m))
cat("Effective sample size (r):", round(effectiveSize(mcmc(res$r)),1), "\n")

# === PARAMETER RECOVERY ===

# --- beta: true value vs. posterior mean and 95% CI ---
alpha = 0.05
ci_beta <- t(apply(res$beta, 2, quantile, c(alpha/2, 1-alpha/2)))
recovery <- data.frame(
  par  = colnames(res$beta),
  true = round(beta, 3),
  low  = round(ci_beta[,1], 3),
  mid  = round(colMeans(res$beta), 3),
  up   = round(ci_beta[,2], 3)
)
print(recovery, row.names = FALSE)

# --- r: true value vs. posterior mean and 95% CI ---
ci_r <- quantile(res$r, c(alpha/2, 1-alpha/2))
cat("true r =", r, " -- mean =", round(mean(res$r),3),
    " -- 95% CI = (", round(ci_r[1],3), ",", round(ci_r[2],3), ")\n")

# === CONVERGENCE PLOT ===

# --- posterior density of each parameter, true value marked in red ---
post_df <- data.frame(
  value = c(res$beta[,1], res$beta[,2], res$r),
  par   = rep(c("beta0 (intercept)", "beta1", "r"), each = nrow(res$beta))
)
true_df <- data.frame(
  par  = c("beta0 (intercept)", "beta1", "r"),
  true = c(beta[1], beta[2], r)
)

print(
  ggplot(post_df, aes(x = value)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    geom_vline(data = true_df, aes(xintercept = true), color = "red", linewidth = 0.8, lty = 2) +
    facet_wrap(~ par, scales = "free") +
    theme_minimal() +
    labs(title = "Posterior density vs. true value", x = "", y = "density")
)

# === BAYESIAN LINEAR MODEL ===

# --- semi-conjugate posterior: beta | sigma2 ~ N(mu, sigma2*S), sigma2 ~ IG(a,b) ---
sample_linear_conjugate <- function(Y, X, n_samples, mu_0 = mu, S_0 = S, a_0=1, b_0=r) {
  n <- nrow(X); p <- ncol(X)
  Lambda_0 <- solve(S_0) # prior precision
  
  # --- posterior for beta given sigma^2: N(M, V * sigma2) ---
  V_n <- solve(t(X) %*% X + Lambda_0)
  M_n <- V_n %*% (t(X) %*% Y + Lambda_0 %*% mu_0)
  
  # --- posterior for sigma2: IG(a_n, b_n) ---
  a_n <- a_0 + n/2
  SSR <- sum((Y - X %*% M_n)^2)
  prior_sq_dist <- t(M_n - mu_0) %*% Lambda_0 %*% (M_n - mu_0)
  b_n <- b_0 + 0.5 * (SSR + prior_sq_dist)
  
  # --- sampling ---
  sig2_samples <- 1 / rgamma(n_samples, a_n, b_n)
  beta_samples <- t(sapply(sig2_samples, function(s2) {
    as.numeric(M_n + t(chol(V_n * s2)) %*% rnorm(p))
  }))
  
  return(list(beta = beta_samples, sigma2 = sig2_samples, M_n = M_n, V_n = V_n, a_n = a_n, b_n = b_n))
}

# === WAIC COMPUTATION ===

# --- WAIC formula ---
calc_WAIC <- function(log_lik_matrix) {
  lpd <- apply(log_lik_matrix, 2, function(col) {
    max_l <- max(col)
    max_l + log(mean(exp(col - max_l)))
  })
  
  p_waic <- apply(log_lik_matrix, 2, var)
  
  waic_val <- -2 * sum(lpd - p_waic)
  return(list(WAIC = waic_val, p_waic = sum(p_waic), lpd = sum(lpd)))
}

res_ghs = res
res_lin_post = sample_linear_conjugate(y,X,n_iter)

# --- GHS log-likelihood matrix ---
n_post <- nrow(res_ghs$beta)
log_lik_ghs_mat <- matrix(0, n_post, n)

for (s in 1:n_post) {
  beta_s <- res_ghs$beta[s, ]
  r_s <- res_ghs$r[s]
  eta <- as.numeric(X %*% beta_s)
  theta <- atan(eta / r_s)
  
  g_part <- 2 * log(Mod(gammaz((r_s + 1i * y) / 2)))
  tilt_part <- theta * y
  const_part <- r_s * log(2 * cos(theta)) - log(4 * pi) - lgamma(r_s)
  
  log_lik_ghs_mat[s, ] <- g_part + tilt_part + const_part
}

# --- linear model log-likelihood matrix ---
log_lik_lm_mat <- matrix(0, n_post, n)
for (s in 1:n_post) {
  beta_s <- res_lin_post$beta[s, ]
  sig2_s <- res_lin_post$sigma2[s]
  
  log_lik_lm_mat[s, ] <- dnorm(y, mean = as.numeric(X %*% beta_s), 
                               sd = sqrt(sig2_s), log = TRUE)
}

# --- final comparison ---
waic_ghs <- calc_WAIC(log_lik_ghs_mat)
waic_lm  <- calc_WAIC(log_lik_lm_mat)

cat("WAIC Comparison (Lower is better):\n")
cat("GHS Model: ", waic_ghs$WAIC, " (Eff. Params:", waic_ghs$p_waic, ")\n")
cat("Linear Model: ", waic_lm$WAIC, " (Eff. Params:", waic_lm$p_waic, ")\n")
