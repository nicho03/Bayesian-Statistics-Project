# ============================================================
# Bayesian GHS GLM applied to real data: temperature variation
# at New York airports (Sec. 5 of the report)
# the code has the following sections:
#
# DATA PREPARATION: download the dataset and set the temperature variation as
#   target variable
#
# MALA WITHIN GIBBS: samples from posterior via MALA (using function in 
#   'ghs_mala.R')
#
# DIAGNOSTIC TRACEPLOT: assess for MCMC variable mixing
#
# BAYESIAN LINEAR MODEL: function to regress data with a semi-conjugate     
#   normall-IG linear model
#
# WAIC COMUPATION: compute WAIC for both the linear model and the GHS
# 
# POSTERIOR CREDIBLE INTERVALS: credible intervals for GHS model
#
# BAYES FACTORS: BF_10 for each coefficient
#
# BAYESIAN STANDARDIZED RESIDUALS AND QQ PLOT: QQ plot for GHS model
# ============================================================

library(nycflights13)
library(pracma)    # gammaz(), used directly below in the WAIC section
library(ggplot2)
library(patchwork)
source("ghs_mala.R")
source("ghs_sampler.R")

set.seed(1)   # for reproducibility of the MCMC run

# === DATA PREPARATION ===
data("weather")
df <- weather[, c("temp", "dewp", "humid", "wind_speed", "precip", "pressure")]

# --- response variable ---
df$delta_temp <- c(NA, diff(df$temp))  # y_t = temp_t - temp_{t-1}

# --- shift the dataframe by 1 ---
df_lagged <- as.data.frame(lapply(df, function(x) c(NA, x[-length(x)])))

# --- filter out anomalies and get y and X ---
final_data <- cbind(delta_temp = df$delta_temp, df_lagged)
final_data <- final_data[complete.cases(final_data), ]
final_data = final_data[-c(10523,10522),]

y <- final_data$delta_temp
X <- model.matrix(delta_temp ~ ., data = final_data)

# === MALA WITHIN GIBBS ===
ols_model = lm(y ~ X - 1)
col_names = colnames(X)

# --- model parameters ---
n = dim(X)[1]
p = dim(X)[2]

S = 5*vcov(ols_model)
mu = rep(0,p)
n_iter = 1e4
BI = 1e3
r = 4

# --- run MALA within GIBBS ---
res = MalaWG(y,X,
             n_iter=n_iter,
             BI=BI,
             var_sample_r=5e-3,
             mu_prior=mu,
             S_prior=S,
             log_prior_r = function(x) dgamma(x,r,1,log=TRUE)
             )

print('MCMC executed')

# === DIAGNOSTIC TRACEPLOTS ===
indices <- 1:(n_iter - BI)

plot(indices, res$r,type = "l", ylab = "r", xlab = "iteration", main = "r traceplot")

#par(mfrow = c(3, 3), mar = c(3.5, 3.5, 2, 1)) 
for (i in seq_along(col_names)) {
  plot(indices, res$beta[, i],type = "l", ylab = "beta", xlab = "iteration", main = paste(col_names[i], "coefficient traceplot"))
}

# === BAYESIAN LINEAR MODEL ===

# Bayesian semi-conjugate linear model
# Model:
#   beta | sigma2 ~ N(mu, sigma2 * S)
#   sigma2 ~ IG(a, b)
sample_linear_conjugate <- function(Y, X, n_samples, mu_0 = mu, S_0 = S, a_0=1, b_0=r) {
  n <- nrow(X); p <- ncol(X)
  Lambda_0 <- solve(S_0) # Prior precision
  
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


#=== WAIC COMPUTATION ===

# --- WAIC formula ---
calc_WAIC <- function(log_lik_matrix) {
  # log pointwise predictive density (lpd) using log-sum-exp for stability
  lpd <- apply(log_lik_matrix, 2, function(col) {
    max_l <- max(col)
    max_l + log(mean(exp(col - max_l)))
  })
  
  # effective number of parameters (p_waic)
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
  
  # pointwise log-likelihood for each observation i
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
  
  # pointwise gaussian log-likelihood
  log_lik_lm_mat[s, ] <- dnorm(y, mean = as.numeric(X %*% beta_s), 
                               sd = sqrt(sig2_s), log = TRUE)
}

# --- final comparison ---
waic_ghs <- calc_WAIC(log_lik_ghs_mat)
waic_lm  <- calc_WAIC(log_lik_lm_mat)

cat("WAIC Comparison (Lower is better):\n")
cat("GHS Model: ", waic_ghs$WAIC, " (Eff. Params:", waic_ghs$p_waic, ")\n")
cat("Linear Model: ", waic_lm$WAIC, " (Eff. Params:", waic_lm$p_waic, ")\n")


# === POSTERIOR CREDIBLE INTERVALS ===

# --- confidence level ---
alpha = 0.05

# --- confidence intervals ---
ci_beta <- t(apply(res$beta, 2, quantile, c(alpha/2,1-alpha/2)))
mean_beta = colMeans(res$beta)

# named ci_df (not df) since df already holds the raw weather data above --
# reusing the same name for two different objects was confusing to read back
ci_df <- data.frame(
  par = colnames(X),
  low = ci_beta[,1],
  mid  = mean_beta,
  up = ci_beta[,2]
  
)

# --- plots ---
df_row1 <- ci_df[1, , drop = FALSE]
df_row6 <- ci_df[6, , drop = FALSE]
df_rest <- ci_df[-c(1,6), ]


make_plot <- function(data){
  ggplot(data, aes(x = par, y = mid)) +
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = up), width = .2) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.6, lty = 2) + 
    coord_flip() +
    theme_minimal() +
    labs(x = "", y = "")
}

p_rest <- make_plot(df_rest)
p_1    <- make_plot(df_row1)
p_6    <- make_plot(df_row6)

h <- c(nrow(df_rest), nrow(df_row6), nrow(df_row1))

print(
  (p_rest / p_6 / p_1) +
    plot_layout(heights = h) +
    plot_annotation(title = "Bayesian 95% Posterior Credible Intervals") &
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(5.5, 12, 5.5, 5.5)
    )
)


# === BAYES FACTORS ===

# H0: beta_i = 0  vs.  H1: beta_i != 0  (Sec. 5.4 of the report)
# BF_01 = pi(beta_i = 0 | y) / pi(beta_i = 0)   -- posterior density over prior
#         density, both evaluated at zero
# BF_10 = 1 / BF_01                              -- values > 1 favour H1,
#         i.e. evidence in favour of keeping the covariate

BF_10 <- sapply(seq_len(p), function(i) {
  post_dens <- density(res$beta[, i])
  posterior_at_zero <- approx(post_dens$x, post_dens$y, xout = 0)$y
  if (is.na(posterior_at_zero)) {
    posterior_at_zero <- dnorm(0, mean(res$beta[, i]), sd(res$beta[, i]))
  }
  prior_at_zero <- dnorm(0, mu[i], sqrt(S[i, i]))
  prior_at_zero / posterior_at_zero
})

bf_df <- data.frame(
  Coefficient = colnames(X),
  BF_10 = ifelse(BF_10 > 1000, ">1000", formatC(BF_10, digits = 2, format = "f"))
)

cat("\nBayes factors (BF_10, evidence in favour of including each covariate):\n")
print(bf_df, row.names = FALSE)


# === BAYESIAN STANDARDIZED RESIDUALS AND QQ PLOT ===

# --- standization ---
beta_hat <- colMeans(res_ghs$beta)
r_hat <- mean(res_ghs$r)
mu_hat <- X %*% beta_hat
var_hat <- r_hat + mu_hat^2 / r_hat
stand_residuals <- (y - mu_hat) / sqrt(var_hat)

# --- standardized residual scatterplot ---
plot(mu_hat, stand_residuals + rnorm(n, 0, 0.2), xlim=c(-1,1), xlab="y_fitted", ylab="posterior predictive residual", main="Posterior predictive residual scatterplot")

# --- QQ-plot ---
theoretical_ghs <- rghs(100000, theta = 0, r = 1) 
observed_sorted <- sort(stand_residuals)
theoretical_quantiles <- quantile(theoretical_ghs, probs = seq(1/n, 1 - 1/n, length.out = n))
plot(theoretical_quantiles, observed_sorted,
     main = "GHS Q-Q Plot",
     xlab = "Theoretical GHS Quantiles",
     ylab = "Residuals",
     pch = 16, col = rgb(0, 0, 0, 0.3))
# --- 45° line ---
abline(0, 1, col = "red", lwd = 2)
