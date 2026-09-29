# ============================================================
# GHS density and exact sampler (Ratio-of-Uniforms, Appendix C)
# ============================================================

# === GHS DENSITY ===
f_GHS <- function(x, theta, r) {
  # log-space evaluation to avoid overflow/underflow for large r or x
  g_mod_sq <- (Mod(gammaz((r + 1i * x) / 2)))^2
  log_const <- r * log(2 * cos(theta)) - log(4 * pi) - lgamma(r)
  return(exp(log_const + log(g_mod_sq) + theta * x))
}

# === GHS SAMPLER ===
rghs <- function(n, theta, r) {
  # Mode of the density: sets the height of the bounding box
  opt_mode <- optimize(function(x) f_GHS(x, theta, r),
                        interval = c(-20, 20), maximum = TRUE)
  u_max <- sqrt(opt_mode$objective)

  # Horizontal boundaries of the Ratio-of-Uniforms region
  v_max <- optimize(function(x) x * sqrt(f_GHS(x, theta, r)),
                     interval = c(0, 100), maximum = TRUE)$objective
  v_min <- optimize(function(x) x * sqrt(f_GHS(x, theta, r)),
                     interval = c(-100, 0), maximum = FALSE)$objective

  samples <- numeric(n)
  accepted <- 0
  while (accepted < n) {
    u <- runif(1, 0, u_max)
    v <- runif(1, v_min, v_max)
    x <- v / u
    if (u <= sqrt(f_GHS(x, theta, r))) {
      accepted <- accepted + 1
      samples[accepted] <- x
    }
  }
  return(samples)
}
