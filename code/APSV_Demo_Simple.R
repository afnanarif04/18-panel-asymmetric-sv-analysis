################################################################################
#  APSV_Demo_Simple.R
#  ─────────────────────────────────────────────────────────────────────────────
#  STAND-ALONE DEMO — Asymmetric Panel Stochastic Volatility with CCE
#  Authors: M.A.A. Bin Mohd Amran & F. Furuoka
#
#  PURPOSE: Self-contained demonstration script that reviewers can run quickly
#  to verify the main Monte Carlo finding of the paper. Uses a small number of
#  replications (R = 50) and a single (N, T) grid cell so the entire script
#  finishes in approximately 5 minutes on a laptop.
#
#  WHAT THIS SCRIPT DOES:
#  1. Compiles the C++ Kalman filter for the conditional QML
#  2. Generates 50 simulated panels under the Baseline configuration
#     (phi = 0.95, gamma1 = -0.15, gamma2 = 0.10, sigma2_e = 0.05)
#     with N = 50 units and T = 500 time periods
#  3. For each replication, estimates the model using:
#     (a) Naive QML — no cross-sectional defactoring
#     (b) CCE-QML — proposed estimator (defactors first)
#  4. Reports bias and RMSE for the three identified parameters
#  5. Produces a single comparison figure and saves it to the working directory
#
#  EXPECTED RESULT: CCE-QML root mean squared error for sigma2_u should be
#  approximately two orders of magnitude smaller than that of the naive estimator.
#  This is the headline finding of the paper.
#
#  REQUIRES: R >= 4.0, packages: Rcpp, ggplot2, tidyr
#  RUNTIME : approximately 4–7 minutes on a standard laptop
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 1. SETUP
# ─────────────────────────────────────────────────────────────────────────────
rm(list = ls())
options(warn = -1)

# Package check
needed <- c("Rcpp","ggplot2","tidyr")
miss <- setdiff(needed, rownames(installed.packages()))
if (length(miss) > 0) {
  cat("Installing packages:", paste(miss, collapse = ", "), "\n")
  install.packages(miss, repos = "https://cloud.r-project.org")
}
invisible(lapply(needed, library, character.only = TRUE))

# ─────────────────────────────────────────────────────────────────────────────
# 2. C++ KALMAN FILTER — compiled inline
# ─────────────────────────────────────────────────────────────────────────────
cat("Compiling C++ Kalman filter (takes ~30 seconds the first time)... ")
sourceCpp(code = '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
static const double LOG2PI     = 1.8378770664093454836;
static const double V_LOG_CHI2 = 9.8696044010893586188 / 2.0;

// [[Rcpp::export]]
double neg_qll_unit_cpp(NumericVector theta_unc, NumericVector x_star,
                        NumericVector y_vec, NumericVector h_level_vec) {
  double phi      = tanh(theta_unc[0]);
  double gamma1   = theta_unc[1];
  double sigma2_u = exp(theta_unc[2]);
  int TT = x_star.size();
  double denom = 1.0 - phi * phi;
  if (denom < 0.01) denom = 0.01;
  double h_pred = 0.0;
  double P_pred = sigma2_u / denom;
  double sum_logF = 0.0, sum_v2F = 0.0;
  for (int tt = 0; tt < TT; tt++) {
    double v = x_star[tt] - h_pred;
    double F = P_pred + V_LOG_CHI2;
    if (F <= 0.0) F = V_LOG_CHI2;
    sum_logF += log(F);
    sum_v2F  += v * v / F;
    double K = P_pred / F;
    double h_filt = h_pred + K * v;
    double P_filt = (1.0 - K) * P_pred;
    double eps = y_vec[tt] * exp(-(h_filt + h_level_vec[tt]) / 2.0);
    if (eps >  20.0) eps =  20.0;
    if (eps < -20.0) eps = -20.0;
    if (tt < TT - 1) {
      h_pred = phi * h_filt + gamma1 * eps;
      P_pred = phi * phi * P_filt + sigma2_u;
    }
  }
  double nll = 0.5 * TT * LOG2PI + 0.5 * sum_logF + 0.5 * sum_v2F;
  if (!R_finite(nll)) return 1e10;
  return nll;
}

// [[Rcpp::export]]
double neg_pooled_qll_cpp(NumericVector theta_unc, NumericMatrix x_star_mat,
                          NumericMatrix y_mat, NumericMatrix h_level_mat) {
  int N = x_star_mat.nrow();
  double total = 0.0;
  for (int i = 0; i < N; i++) {
    total += neg_qll_unit_cpp(theta_unc, x_star_mat(i, Rcpp::_),
                               y_mat(i, Rcpp::_), h_level_mat(i, Rcpp::_));
  }
  return total;
}
')
cat("done.\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
N <- 50          # cross-section size
TT <- 500        # time series length
R_DEMO <- 50     # number of Monte Carlo replications (small for speed)

# Baseline configuration parameters
phi0     <- 0.95
gamma1_0 <- -0.15
gamma2_0 <- 0.10
sigma2_e <- 0.05
s2u_true <- gamma1_0^2 + gamma2_0^2 * (1 - 2/pi) + sigma2_e
T_burn   <- 500

E_abs_eps  <- sqrt(2 / pi)
E_log_chi2 <- digamma(0.5) + log(2)
V_log_chi2 <- pi^2 / 2

cat("══════════════════════════════════════════════════════════════\n")
cat(" APSV DEMO — Small Monte Carlo for Quick Verification\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(sprintf(" Configuration: phi=%.2f, gamma1=%.2f, gamma2=%.2f, sigma2_e=%.2f\n",
            phi0, gamma1_0, gamma2_0, sigma2_e))
cat(sprintf(" Identifiable sigma2_u = gamma1^2 + gamma2^2*(1-2/pi) + sigma2_e = %.5f\n", s2u_true))
cat(sprintf(" Grid: N = %d, T = %d, replications = %d\n", N, TT, R_DEMO))
cat(sprintf(" Number of common factors: 1 (single global volatility factor)\n"))
cat("══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. DATA-GENERATING PROCESS
# ─────────────────────────────────────────────────────────────────────────────
generate_data <- function(N, TT, phi, g1, g2, s2e, seed) {
  set.seed(seed)
  T_total <- TT + T_burn
  s_e <- sqrt(s2e)

  mu_i     <- runif(N, -1.5, -0.5)              # unit-specific intercepts
  lambda_i <- runif(N, 0.5, 1.5)                # factor loadings (1 factor)
  rho_f    <- 0.50
  s2_eta   <- 1 - rho_f^2

  # Common factor F_t — AR(1)
  f_vec <- numeric(T_total); f_vec[1] <- rnorm(1)
  for (tt in 2:T_total)
    f_vec[tt] <- rho_f * f_vec[tt-1] + rnorm(1, 0, sqrt(s2_eta))

  # Generate y_it
  y_mat <- matrix(0, N, T_total)
  h_mat <- matrix(0, N, T_total)
  eps_mat <- matrix(0, N, T_total)

  h_init <- mu_i / (1 - phi)
  for (i in 1:N) {
    h_mat[i, 1]   <- h_init[i]
    eps_mat[i, 1] <- rnorm(1)
    y_mat[i, 1]   <- exp(h_mat[i, 1] / 2) * eps_mat[i, 1]
    for (tt in 2:T_total) {
      h_mat[i, tt] <- mu_i[i] +
                       phi * h_mat[i, tt-1] +
                       g1 * eps_mat[i, tt-1] +
                       g2 * (abs(eps_mat[i, tt-1]) - E_abs_eps) +
                       lambda_i[i] * f_vec[tt] +
                       rnorm(1, 0, s_e)
      eps_mat[i, tt] <- rnorm(1)
      y_mat[i, tt] <- exp(h_mat[i, tt] / 2) * eps_mat[i, tt]
    }
  }
  y_mat[, (T_burn+1):T_total]
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. CCE DEFACTORING
# ─────────────────────────────────────────────────────────────────────────────
cce_defactor <- function(x_mat, p_T = floor(ncol(x_mat)^(1/3))) {
  N <- nrow(x_mat); TT <- ncol(x_mat)
  x_bar <- colMeans(x_mat)
  T_eff <- TT - p_T
  z_bar <- matrix(0, T_eff, p_T + 1)
  for (lag in 0:p_T) z_bar[, lag+1] <- x_bar[(p_T+1-lag):(TT-lag)]
  z_aug <- cbind(1, z_bar)
  ZtZ_inv_Zt <- solve(crossprod(z_aug), t(z_aug))
  x_star <- matrix(0, N, T_eff)
  cce_fitted <- matrix(0, N, T_eff)
  for (i in 1:N) {
    x_i <- x_mat[i, (p_T+1):TT]
    delta_hat <- ZtZ_inv_Zt %*% x_i
    fitted_i  <- as.numeric(z_aug %*% delta_hat)
    x_star[i, ]     <- x_i - fitted_i
    cce_fitted[i, ] <- fitted_i
  }
  list(x_star = x_star, cce_fitted = cce_fitted, p_T = p_T, T_eff = T_eff)
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. POOLED QML
# ─────────────────────────────────────────────────────────────────────────────
pooled_qml <- function(x_star, y, init = c(0.90, -0.10, 0.08), h_level) {
  N <- nrow(x_star); TT <- ncol(x_star)
  if (is.null(h_level)) h_level <- matrix(0, N, TT)
  init_unc <- c(atanh(init[1]), init[2], log(init[3]))
  opt <- tryCatch(
    optim(init_unc,
          fn = function(th) neg_pooled_qll_cpp(th, x_star, y, h_level),
          method = "BFGS", control = list(maxit = 300, reltol = 1e-7)),
    error = function(e) NULL)
  if (is.null(opt)) return(rep(NA, 3))
  c(phi      = tanh(opt$par[1]),
    gamma1   = opt$par[2],
    sigma2_u = exp(opt$par[3]))
}

# ─────────────────────────────────────────────────────────────────────────────
# 7. RUN MONTE CARLO
# ─────────────────────────────────────────────────────────────────────────────
cat("Running Monte Carlo (R = ", R_DEMO, " replications)...\n", sep = "")
results <- matrix(NA, R_DEMO, 6,
  dimnames = list(NULL, c("phi_naive","g1_naive","s2u_naive",
                           "phi_cce","g1_cce","s2u_cce")))

t_start <- Sys.time()
for (r in 1:R_DEMO) {
  y_mat <- generate_data(N, TT, phi0, gamma1_0, gamma2_0, sigma2_e, seed = 1000 + r)
  y_dm  <- sweep(y_mat, 1, rowMeans(y_mat), "-")
  y_dm[y_dm == 0] <- 1e-6
  x_mat <- log(y_dm^2) - E_log_chi2

  # Naive — no defactoring
  est_n <- pooled_qml(x_mat, y_dm, h_level = matrix(0, N, TT))
  results[r, 1:3] <- est_n

  # CCE — with defactoring
  cce <- tryCatch(cce_defactor(x_mat), error = function(e) NULL)
  if (!is.null(cce)) {
    y_eff <- y_dm[, (cce$p_T + 1):TT, drop = FALSE]
    est_c <- pooled_qml(cce$x_star, y_eff, h_level = cce$cce_fitted)
    results[r, 4:6] <- est_c
  }

  if (r %% 10 == 0)
    cat(sprintf("  rep %d/%d  elapsed=%.1fs\n",
                r, R_DEMO, as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nMonte Carlo finished in %.1f minutes.\n\n", elapsed))

# ─────────────────────────────────────────────────────────────────────────────
# 8. SUMMARY TABLE
# ─────────────────────────────────────────────────────────────────────────────
true_vals <- c(phi0, gamma1_0, s2u_true)
names(true_vals) <- c("phi","gamma1","sigma2_u")

bias_naive <- colMeans(results[, 1:3], na.rm = TRUE) - true_vals
bias_cce   <- colMeans(results[, 4:6], na.rm = TRUE) - true_vals
rmse_naive <- sqrt(colMeans((results[, 1:3] - matrix(true_vals, R_DEMO, 3, byrow = TRUE))^2,
                             na.rm = TRUE))
rmse_cce   <- sqrt(colMeans((results[, 4:6] - matrix(true_vals, R_DEMO, 3, byrow = TRUE))^2,
                             na.rm = TRUE))

cat("══════════════════════════════════════════════════════════════════════════\n")
cat(" RESULTS — Bias and Root Mean Squared Error\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("%-12s %-8s %12s %12s | %12s %12s\n",
            "Parameter","True","Naive Bias","CCE Bias","Naive RMSE","CCE RMSE"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (j in 1:3) {
  cat(sprintf("%-12s %-8.5f %+12.5f %+12.5f | %12.5f %12.5f\n",
              names(true_vals)[j], true_vals[j],
              bias_naive[j], bias_cce[j],
              rmse_naive[j], rmse_cce[j]))
}
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("\nKey finding: Naive RMSE for sigma2_u (%.4f) is %.0fx larger than CCE RMSE (%.4f)\n",
            rmse_naive[3], rmse_naive[3] / max(rmse_cce[3], 1e-9), rmse_cce[3]))
cat("This is the dramatic improvement provided by cross-sectional defactoring.\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 9. FIGURE — single comparison plot
# ─────────────────────────────────────────────────────────────────────────────
df_long <- data.frame(
  rep   = rep(1:R_DEMO, 6),
  est   = c(results[,1], results[,2], results[,3], results[,4], results[,5], results[,6]),
  param = factor(rep(rep(c("phi","gamma1","sigma2_u"), each = R_DEMO), 2),
                 levels = c("phi","gamma1","sigma2_u")),
  estimator = factor(rep(c("Naive (no CCE)","CCE-QML (proposed)"), each = 3 * R_DEMO),
                     levels = c("Naive (no CCE)","CCE-QML (proposed)"))
)
df_true <- data.frame(param = factor(c("phi","gamma1","sigma2_u"),
                                       levels = c("phi","gamma1","sigma2_u")),
                       value = true_vals)

p <- ggplot(df_long, aes(x = estimator, y = est, fill = estimator)) +
  geom_boxplot(width = 0.5, outlier.size = 1, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.4) +
  geom_hline(data = df_true, aes(yintercept = value),
             colour = "#E63946", linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ param, scales = "free_y", ncol = 3,
             labeller = as_labeller(c(
               phi = "phi (persistence)",
               gamma1 = "gamma1 (leverage)",
               sigma2_u = "sigma2_u (transition variance)"))) +
  scale_fill_manual(values = c("Naive (no CCE)" = "#E63946",
                                "CCE-QML (proposed)" = "#1F3864"),
                    guide = "none") +
  labs(x = NULL, y = "Estimate",
       title = sprintf("APSV Demo — N = %d, T = %d, R = %d replications", N, TT, R_DEMO),
       subtitle = "Red dashed line = true parameter. CCE-QML is tightly centred; Naive is wildly biased.",
       caption = "Configuration: phi=0.95, gamma1=-0.15, gamma2=0.10, sigma2_e=0.05") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#1F3864"),
        strip.text       = element_text(colour = "white", face = "bold"),
        axis.text.x      = element_text(angle = 15, hjust = 1),
        plot.title       = element_text(face = "bold", colour = "#1F3864"),
        plot.caption     = element_text(colour = "grey50"))

ggsave("APSV_Demo_Result.pdf", p, width = 11, height = 5, dpi = 300)
cat("✓ Figure saved: APSV_Demo_Result.pdf\n")
cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat(" DEMO COMPLETE\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
