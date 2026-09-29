# ============================================================
# Bayesian GHS GLM -- Metropolis-within-Gibbs with a MALA step
# ============================================================
#
# Implements the algorithm described in Sec. 4.2 of the report
# ("Adding a MALA step to improve mixing"):
#   - beta is updated with a Metropolis-Adjusted Langevin proposal,
#     computed in a "whitened" parametrization z (so that the
#     Gaussian prior N(mu_prior, S_prior) becomes a standard normal
#     on z, formula (12) of the report);
#   - r is updated with a Metropolis-Hastings step using a Gamma
#     proposal centered on the current value (Sec. 4.1, formula (18)).
#
# Model (Sec. 4, formula (16)):
#   y_i | theta_i, r ~ GHS(theta_i, r)
#   theta_i = atan( (x_i^T beta) / r )
#   beta ~ N(mu_prior, S_prior)
#   r    ~ log_prior_r   (user-supplied, e.g. a Gamma density)
library(pracma)   # gammaz(): complex Gamma function used in the GHS density (formula (7))
library(progress)

MalaWG <- function(Y, X, n_iter, BI = 0,
                    var_sample_r,
                    mu_prior, S_prior,
                    log_prior_r,
                    epsilon = 0.1,
                    target_alpha = 0.574,
                    gamma_exponent = 0.7,
                    add_intercept = FALSE) {

  check <- function(u, log_alpha) {
    if (is.na(log_alpha) || is.nan(log_alpha) || is.infinite(log_alpha)) return(FALSE)
    return(u < log_alpha)
  }

  X <- as.matrix(X)
  Y <- as.numeric(Y)
  if (add_intercept) {
    has_intercept_like <- (ncol(X) >= 1 && all(abs(X[, 1] - 1) < 1e-12))
    if (!has_intercept_like) X <- cbind(1, X)
  }

  n <- nrow(X); p <- ncol(X)
  U <- chol(S_prior)
  L <- t(U)

  beta_from_z <- function(z) as.numeric(mu_prior + L %*% z)

  # --- log-likelihood of beta given r ---
  log_lik <- function(beta, r) {
    eta <- as.numeric(X %*% beta)
    b <- eta / r
    sum(Y * atan(b) - (r / 2) * log(1 + b^2))
  }

  # --- gradient of the log-likelihood wrt beta, formula (13) ---
  grad_log_lik_beta <- function(beta, r) {
    eta <- as.numeric(X %*% beta)
    b <- eta / r
    # d/d eta_i: (Y_i - eta_i) / ( r * (1 + (eta_i/r)^2 ) )
    w <- (Y - eta) / (r * (1 + b^2))   # length n
    as.numeric(crossprod(X, w))        # length p
  }

  log_post_z <- function(z, r) {
    beta <- beta_from_z(z)
    log_lik(beta, r) - 0.5 * sum(z * z)
  }

  grad_log_post_z <- function(z, r) {
    beta <- beta_from_z(z)
    g_beta <- grad_log_lik_beta(beta, r)
    g_z <- as.numeric(t(L) %*% g_beta)  # chain rule d/dz = L^T * grad_beta
    g_z - z
  }

  # --- posterior for r (likelihood of r given beta + prior) ---
  log_posterior_r <- function(beta, r) {
    eta <- as.numeric(X %*% beta)
    theta <- atan(eta / r)  # <-- reparametrization
    g <- gammaz((r + 1i * Y) / 2)
    lg_g <- 2 * sum(log(Mod(g)))
    const_part <- r * sum(log(2 * cos(theta))) - n * lgamma(r)
    const_part + lg_g + log_prior_r(r) + sum(theta * Y)
  }


  pb <- progress_bar$new(
  format = " MWG [:bar] :percent eta: :eta",
  total = n_iter,
  stream = stderr(),
  clear = FALSE,
  width = 60
  )
  
  # --- init
  z_cur <- rep(0, p); r_cur <- 1
  beta_sim <- matrix(NA_real_, n_iter, p); r_sim <- numeric(n_iter)

  # --- adaptive tuning of epsilon (MALA step size)
  log_epsilon <- log(epsilon)
  accepted_total_beta <- 0; accepted_adapt_beta <- 0
  sigma2 <- (exp(log_epsilon))^2 / p^(1 / 3)

  logp_cur <- log_post_z(z_cur, r_cur)
  grad_cur <- grad_log_post_z(z_cur, r_cur)

  for (i in 1:n_iter) {
    # ===== beta step: MALA on z =====
    z_mean <- z_cur + (sigma2 / 2) * grad_cur
    z_new <- as.numeric(z_mean + sqrt(sigma2) * rnorm(p))

    logp_new <- log_post_z(z_new, r_cur)
    grad_new <- grad_log_post_z(z_new, r_cur)

    qold <- -sum((z_cur - z_new - (sigma2 / 2) * grad_new)^2) / (2 * sigma2)
    qnew <- -sum((z_new - z_cur - (sigma2 / 2) * grad_cur)^2) / (2 * sigma2)

    log_alpha_beta <- (logp_new - logp_cur) + (qold - qnew)

    if (check(log(runif(1)), log_alpha_beta)) {
      z_cur <- z_new; logp_cur <- logp_new; grad_cur <- grad_new
      accepted_total_beta <- accepted_total_beta + 1
      if (i <= BI) accepted_adapt_beta <- accepted_adapt_beta + 1
    }

    if (i <= BI && BI > 0) {
      alpha_hat <- accepted_adapt_beta / i
      gamma_k <- 1 / (i^gamma_exponent)
      log_epsilon <- log_epsilon + gamma_k * (alpha_hat - target_alpha)
      sigma2 <- (exp(log_epsilon))^2 / p^(1 / 3)
    }

    # === map z into beta ===
    beta_cur <- beta_from_z(z_cur)

    # ===== r step (gamma proposal)  =====
    a <- r_cur^2 / var_sample_r; b <- r_cur / var_sample_r
    r_new <- rgamma(1, a, b)

    log_alpha_r <- log_posterior_r(beta_cur, r_new) - log_posterior_r(beta_cur, r_cur) +
      dgamma(r_cur, r_new^2 / var_sample_r, r_new / var_sample_r, log = TRUE) -
      dgamma(r_new, a, b, log = TRUE)

    if (check(log(runif(1)), log_alpha_r)) {
      r_cur <- r_new
      # refresh current posterior and gradient
      logp_cur <- log_post_z(z_cur, r_cur)   
      grad_cur <- grad_log_post_z(z_cur, r_cur)
    }

    beta_sim[i, ] <- beta_cur; r_sim[i] <- r_cur
    if (i %% 100 == 0 || i == n_iter) {
      pb$tick(100)
}
  }

  keep <- (BI + 1):n_iter
  list(beta = beta_sim[keep, , drop = FALSE],
       r = r_sim[keep],
       acc_rate_beta = accepted_total_beta / n_iter,
       epsilon_final = exp(log_epsilon))
}
