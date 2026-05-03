################################################################################
#  00_Functions.R  (FINAL — Three-Parameter Conditional KF)
#  ═════════════════════════════════════════════════════════════════════════════
#  ESTIMATOR: Pooled CCE-QML with conditional Kalman filter
#
#  IDENTIFIED PARAMETERS:
#    phi      — persistence (AR coefficient in log-volatility)
#    gamma1   — leverage effect (sign-dependent asymmetry)
#    sigma2_u — total transition variance = gamma1^2 + gamma2^2*(1-2/pi) + sigma2_e
#
#  WHY THREE (NOT FOUR):
#    gamma2 (size effect) and sigma2_e (idiosyncratic innovation variance) enter
#    the transition equation symmetrically — both contribute additive variance.
#    The log-squared QML cannot separate them: the QLL surface has a ridge
#    along (gamma2, sigma2_e) pairs with equal sigma2_u. This is structural
#    non-identification, not finite-sample bias.
#
#    gamma1 IS identified because the sign of eps_hat — recovered from
#    y_it * exp(-(h_filt + cce_fitted_t)/2) — provides asymmetric information
#    that breaks the symmetry.
#
#  ESTIMATION STEPS:
#    Step 1: CCE defactoring with time-varying cce_fitted (delta_hat_i' * z_bar_t)
#    Step 2: Pooled conditional KF-QML for (phi, gamma1, sigma2_u)
#            - prediction: h_pred = phi*h_filt + gamma1*eps_hat
#            - P_pred = phi^2*P_filt + sigma2_u
#            - eps_hat = y * exp(-(h_filt + cce_fitted_t)/2)
#    Step 3: ABC bias correction for phi (optional)
#    Step 4: Naive estimator (no CCE) as benchmark
#
#  Paper: "Asymmetric Panel Stochastic Volatility with Cross-Sectional Dependence"
#  Authors: M.A.A. Bin Mohd Amran & F. Furuoka
################################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 0. WINDOWS DLL FIX + RCPP
# ═══════════════════════════════════════════════════════════════════════════════

user_tmp <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R_tmp")
dir.create(user_tmp, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMP = user_tmp, TEMP = user_tmp, TMPDIR = user_tmp)
rcpp_cache <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R_rcpp_cache")
dir.create(rcpp_cache, showWarnings = FALSE, recursive = TRUE)

library(Rcpp)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. C++ KALMAN FILTER — Three-parameter conditional KF
# ═══════════════════════════════════════════════════════════════════════════════

cpp_code_string <- '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static const double PI_VAL     = 3.14159265358979323846;
static const double LOG2PI     = 1.8378770664093454836;
static const double V_LOG_CHI2 = PI_VAL * PI_VAL / 2.0;
static const double E_ABS_EPS  = 0.7978845608028654;

// Conditional KF for one unit: estimates (phi, gamma1, sigma2_u)
// Prediction uses gamma1*eps_hat (sign-dependent leverage)
// Transition variance uses sigma2_u (total, absorbs gamma2 and sigma2_e)
// Level correction from cce_fitted (time-varying vector)
//
// theta_unc = (atanh(phi), gamma1, log(sigma2_u))
//
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

    // Recover standardised return using time-varying level
    double eps = y_vec[tt] * exp(-(h_filt + h_level_vec[tt]) / 2.0);
    if (eps >  20.0) eps =  20.0;
    if (eps < -20.0) eps = -20.0;

    if (tt < TT - 1) {
      // Prediction: only gamma1 (sign-dependent leverage)
      h_pred = phi * h_filt + gamma1 * eps;
      // Transition variance: full sigma2_u
      P_pred = phi * phi * P_filt + sigma2_u;
    }
  }

  double nll = 0.5 * TT * LOG2PI + 0.5 * sum_logF + 0.5 * sum_v2F;
  if (!R_finite(nll)) return 1e10;
  return nll;
}

// Pooled QLL across N units
// [[Rcpp::export]]
double neg_pooled_qll_cpp(NumericVector theta_unc, NumericMatrix x_star_mat,
                          NumericMatrix y_mat, NumericMatrix h_level_mat) {
  int N = x_star_mat.nrow();
  double total = 0.0;
  for (int i = 0; i < N; i++) {
    total += neg_qll_unit_cpp(theta_unc,
                              x_star_mat(i, Rcpp::_),
                              y_mat(i, Rcpp::_),
                              h_level_mat(i, Rcpp::_));
  }
  return total;
}
'

cat("Compiling C++ Kalman filter (3-param conditional) ... ")
sourceCpp(code = cpp_code_string, cacheDir = rcpp_cache, rebuild = TRUE)
cat("done.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. GLOBALS
# ═══════════════════════════════════════════════════════════════════════════════

R       <- 500
N_grid  <- c(30, 50, 100, 200, 400)
T_grid  <- c(250, 500, 1000, 2000)
T_burn  <- 500
run_hpj <- FALSE

r_factors <- 1;    rho_f      <- 0.50
lambda_lo <- 0.5;  lambda_hi  <- 1.5
mu_lo     <- -1.5; mu_hi     <- -0.5

E_abs_eps  <- sqrt(2 / pi)
E_log_chi2 <- digamma(0.5) + log(2)
V_log_chi2 <- pi^2 / 2
sigma2_a   <- 1 - 2 / pi     # Var(|eps| - E|eps|)

configs <- list(
  Baseline = list(phi = 0.95, gamma1 = -0.15, gamma2 = 0.10, sigma2_e = 0.05),
  A        = list(phi = 0.98, gamma1 = -0.15, gamma2 = 0.10, sigma2_e = 0.05),
  B        = list(phi = 0.90, gamma1 = -0.15, gamma2 = 0.10, sigma2_e = 0.05),
  C        = list(phi = 0.95, gamma1 = -0.20, gamma2 = 0.00, sigma2_e = 0.05),
  D        = list(phi = 0.95, gamma1 = -0.20, gamma2 = 0.15, sigma2_e = 0.05))

# ═══════════════════════════════════════════════════════════════════════════════
# 3. DGP (unchanged — still generates all 4 parameters)
# ═══════════════════════════════════════════════════════════════════════════════

generate_apsv_data <- function(N, TT, phi, gamma1, gamma2, sigma2_e,
                                r = 1, rho_f = 0.50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  T_total <- TT + T_burn; sigma_e <- sqrt(sigma2_e)
  mu_i     <- runif(N, mu_lo, mu_hi)
  lambda_i <- matrix(runif(N * r, lambda_lo, lambda_hi), N, r)
  sigma2_eta <- 1 - rho_f^2
  f_mat <- matrix(0, r, T_total); f_mat[, 1] <- rnorm(r)
  for (tt in 2:T_total)
    f_mat[, tt] <- rho_f * f_mat[, tt-1] + rnorm(r, 0, sqrt(sigma2_eta))
  y_mat <- h_mat <- eps_mat <- matrix(0, N, T_total)
  h_init <- mu_i / (1 - phi)
  for (i in 1:N) {
    h_mat[i, 1] <- h_init[i]; eps_mat[i, 1] <- rnorm(1)
    y_mat[i, 1] <- exp(h_mat[i, 1] / 2) * eps_mat[i, 1]
    for (tt in 2:T_total) {
      h_mat[i, tt] <- mu_i[i] + phi * h_mat[i, tt-1] +
        gamma1 * eps_mat[i, tt-1] +
        gamma2 * (abs(eps_mat[i, tt-1]) - E_abs_eps) +
        sum(lambda_i[i, ] * f_mat[, tt]) + rnorm(1, 0, sigma_e)
      eps_mat[i, tt] <- rnorm(1)
      y_mat[i, tt] <- exp(h_mat[i, tt] / 2) * eps_mat[i, tt]
    }
  }
  idx <- (T_burn + 1):T_total
  list(y = y_mat[, idx],
       true_params = c(phi = phi, gamma1 = gamma1, gamma2 = gamma2, sigma2_e = sigma2_e))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. CCE DEFACTORING — returns cce_fitted (time-varying level correction)
# ═══════════════════════════════════════════════════════════════════════════════

cce_defactor <- function(x_mat, p_T = NULL) {
  N <- nrow(x_mat); TT <- ncol(x_mat)
  if (is.null(p_T)) p_T <- floor(TT^(1/3))
  x_bar <- colMeans(x_mat); T_eff <- TT - p_T
  z_bar <- matrix(0, T_eff, p_T + 1)
  for (lag in 0:p_T) z_bar[, lag+1] <- x_bar[(p_T+1-lag):(TT-lag)]
  z_aug <- cbind(1, z_bar)
  ZtZ_inv_Zt <- solve(crossprod(z_aug), t(z_aug))
  x_star     <- matrix(0, N, T_eff)
  cce_fitted <- matrix(0, N, T_eff)
  for (i in 1:N) {
    x_i <- x_mat[i, (p_T+1):TT]
    delta_hat <- ZtZ_inv_Zt %*% x_i
    fitted_i  <- as.numeric(z_aug %*% delta_hat)
    x_star[i, ]     <- x_i - fitted_i
    cce_fitted[i, ] <- fitted_i
  }
  list(x_star = x_star, cce_fitted = cce_fitted, T_eff = T_eff, p_T = p_T)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. POOLED QML — Three-parameter: (phi, gamma1, sigma2_u)
# ═══════════════════════════════════════════════════════════════════════════════

pooled_qml <- function(x_star_mat, y_mat, init = c(0.90, -0.10, 0.08),
                       h_level_mat = NULL, compute_se = TRUE) {
  N <- nrow(x_star_mat); TT <- ncol(x_star_mat)
  if (is.null(h_level_mat)) h_level_mat <- matrix(0, N, TT)

  neg_f <- function(th) neg_pooled_qll_cpp(th, x_star_mat, y_mat, h_level_mat)

  # init = c(phi, gamma1, sigma2_u) in natural space
  init_unc <- c(atanh(init[1]), init[2], log(init[3]))

  opt <- tryCatch(
    optim(init_unc, neg_f, method = "BFGS",
          control = list(maxit = 500, reltol = 1e-8), hessian = compute_se),
    error = function(e) NULL)

  if (is.null(opt))
    return(list(theta = rep(NA, 3), se = rep(NA, 3), converged = FALSE))

  theta_hat <- c(phi      = tanh(opt$par[1]),
                 gamma1   = opt$par[2],
                 sigma2_u = exp(opt$par[3]))

  se_hat <- rep(NA, 3)
  if (compute_se) {
    J <- diag(c(1 - tanh(opt$par[1])^2, 1, exp(opt$par[3])))
    V <- tryCatch(solve(opt$hessian), error = function(e) diag(NA, 3))
    se_hat <- sqrt(abs(diag(J %*% V %*% t(J))))
  }

  list(theta = theta_hat, se = se_hat, loglik = -opt$value,
       converged = (opt$convergence == 0))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. BIAS CORRECTION
# ═══════════════════════════════════════════════════════════════════════════════

abc_correction_phi <- function(phi_hat, sigma2_u_hat, TT) {
  snr <- sigma2_u_hat / (sigma2_u_hat + V_log_chi2)
  phi_abc <- phi_hat + (1 + phi_hat) / TT * snr
  max(min(phi_abc, 0.999), -0.999)
}

hpj_correction <- function(x_star_mat, y_mat, init, h_level_mat) {
  Te <- ncol(x_star_mat); Th <- floor(Te / 2)
  ef <- pooled_qml(x_star_mat, y_mat, init, h_level_mat, FALSE)
  e1 <- pooled_qml(x_star_mat[, 1:Th, drop=F], y_mat[, 1:Th, drop=F],
                   init, h_level_mat[, 1:Th, drop=F], FALSE)
  e2 <- pooled_qml(x_star_mat[, (Th+1):Te, drop=F], y_mat[, (Th+1):Te, drop=F],
                   init, h_level_mat[, (Th+1):Te, drop=F], FALSE)
  if (any(is.na(ef$theta)) || any(is.na(e1$theta)) || any(is.na(e2$theta)))
    return(list(theta = rep(NA, 3)))
  th <- 2 * ef$theta - 0.5 * (e1$theta + e2$theta)
  th[1] <- max(min(th[1], 0.999), -0.999)   # phi
  th[3] <- max(th[3], 1e-6)                  # sigma2_u
  list(theta = th)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 7. SINGLE REPLICATION
# ═══════════════════════════════════════════════════════════════════════════════

run_one_replication <- function(N, TT, cfg, rep_id) {
  # Generate data (full 4-param DGP)
  dat <- generate_apsv_data(N, TT, cfg$phi, cfg$gamma1, cfg$gamma2, cfg$sigma2_e,
                            r = r_factors, rho_f = rho_f,
                            seed = rep_id * 1000 + N + TT)
  y_mat <- dat$y; y_mat[y_mat == 0] <- 1e-6
  x_mat <- log(y_mat^2) - E_log_chi2
  p_T <- floor(TT^(1/3))

  # True sigma2_u (what we're estimating instead of separate gamma2, sigma2_e)
  true_s2u <- cfg$gamma1^2 + cfg$gamma2^2 * sigma2_a + cfg$sigma2_e

  # Step 1: CCE defactoring
  cce <- cce_defactor(x_mat, p_T)
  x_star     <- cce$x_star
  cce_fitted <- cce$cce_fitted
  y_eff <- y_mat[, (p_T+1):TT]

  # Step 2: Pooled conditional KF-QML for (phi, gamma1, sigma2_u)
  init <- c(0.90, cfg$gamma1, true_s2u)
  est_qml <- pooled_qml(x_star, y_eff, init, cce_fitted, TRUE)

  # HPJ bias correction (if enabled)
  if (run_hpj) {
    est_hpj <- hpj_correction(x_star, y_eff, init, cce_fitted)
  } else {
    est_hpj <- list(theta = rep(NA, 3))
  }

  # ABC bias correction (phi only)
  if (any(is.na(est_qml$theta))) {
    est_abc <- list(theta = rep(NA, 3))
  } else {
    phi_abc <- abc_correction_phi(est_qml$theta[1], est_qml$theta[3], ncol(x_star))
    est_abc <- list(theta = c(phi_abc, est_qml$theta[2], est_qml$theta[3]))
  }

  # Naive: no CCE defactoring
  x_naive <- x_mat[, (p_T+1):TT]
  h_level_naive <- matrix(rowMeans(x_mat), nrow = N, ncol = ncol(x_naive))
  est_naive <- pooled_qml(x_naive, y_eff, init, h_level_naive, FALSE)

  # Collect results — output 3 identified params per estimator
  # Map: slot 1 = phi, slot 2 = gamma1, slot 3 = sigma2_u
  #       slot 4 = NA (placeholder for backward compat with 21_Compile)
  q  <- unname(est_qml$theta)
  h  <- unname(est_hpj$theta)
  a  <- unname(est_abc$theta)
  nv <- unname(est_naive$theta)
  se <- unname(est_qml$se)
  if (is.null(se) || length(se) < 3) se <- rep(NA_real_, 3)

  c(QML_phi  = q[1], QML_g1  = q[2], QML_g2 = NA_real_, QML_s2u  = q[3],
    HPJ_phi  = h[1], HPJ_g1  = h[2], HPJ_g2 = NA_real_, HPJ_s2u  = h[3],
    ABC_phi  = a[1], ABC_g1  = a[2], ABC_g2 = NA_real_, ABC_s2u  = a[3],
    Naive_phi = nv[1], Naive_g1 = nv[2], Naive_g2 = NA_real_, Naive_s2u = nv[3],
    QML_conv = ifelse(is.null(est_qml$converged), 0, as.numeric(est_qml$converged)),
    Naive_conv = ifelse(is.null(est_naive$converged), 0, as.numeric(est_naive$converged)),
    QML_se_phi = se[1], QML_se_g1 = se[2], QML_se_s2u = se[3])
}

# ═══════════════════════════════════════════════════════════════════════════════
# 8. SUMMARY STATISTICS
# ═══════════════════════════════════════════════════════════════════════════════

cn21 <- c("QML_phi","QML_g1","QML_g2","QML_s2u",
          "HPJ_phi","HPJ_g1","HPJ_g2","HPJ_s2u",
          "ABC_phi","ABC_g1","ABC_g2","ABC_s2u",
          "Naive_phi","Naive_g1","Naive_g2","Naive_s2u",
          "QML_conv","Naive_conv",
          "QML_se_phi","QML_se_g1","QML_se_s2u")

summarise_mc <- function(results, true_params) {
  ests <- c("QML", "HPJ", "ABC", "Naive")
  pars <- c("phi", "g1", "s2u")  # 3 identified params (skip g2)

  if (is.list(true_params)) {
    true_s2u <- as.numeric(true_params$gamma1)^2 +
                as.numeric(true_params$gamma2)^2 * sigma2_a +
                as.numeric(true_params$sigma2_e)
    pt <- c(as.numeric(true_params$phi),
            as.numeric(true_params$gamma1),
            true_s2u)
  } else {
    true_s2u <- as.numeric(true_params["gamma1"])^2 +
                as.numeric(true_params["gamma2"])^2 * sigma2_a +
                as.numeric(true_params["sigma2_e"])
    pt <- c(as.numeric(true_params["phi"]),
            as.numeric(true_params["gamma1"]),
            true_s2u)
  }

  out <- data.frame()
  for (e in ests) {
    for (j in seq_along(pars)) {
      cn <- paste0(e, "_", pars[j])
      if (cn %in% colnames(results)) {
        v <- results[, cn]; v <- v[is.finite(v)]
        if (length(v) == 0) next
        b <- mean(v) - pt[j]; s <- sd(v); r <- sqrt(b^2 + s^2)
        out <- rbind(out, data.frame(
          Estimator = e, Parameter = pars[j], True = round(pt[j], 5),
          Bias = round(b, 5), SD = round(s, 5), RMSE = round(r, 5),
          stringsAsFactors = FALSE))
      }
    }
  }
  out
}

# ═══════════════════════════════════════════════════════════════════════════════
# 9. CELL RUNNER (sequential, auto-resume, incremental save)
# ═══════════════════════════════════════════════════════════════════════════════

run_cell <- function(N, TT, cfg, cfg_name, R, save_prefix = "cell") {
  key <- paste0("N", N, "_T", TT)
  save_file <- sprintf("%s_%s_%s.RData", save_prefix, cfg_name, key)
  partial_file <- sprintf("partial_%s_%s.RData", cfg_name, key)

  if (file.exists(save_file)) {
    cat(sprintf("  N=%3d T=%4d SKIP (done)\n", N, TT))
    return(invisible(NULL))
  }

  cat(sprintf("  N=%3d T=%4d ", N, TT))
  t0 <- Sys.time()

  results_list <- vector("list", R)
  start_rep <- 1

  if (file.exists(partial_file)) {
    load(partial_file)
    if (!is.null(partial_data) && length(partial_data) > 0) {
      n_done <- length(partial_data)
      results_list[1:n_done] <- partial_data
      start_rep <- n_done + 1
      cat(sprintf("(resume %d) ", start_rep))
    }
  }

  if (start_rep <= R) {
    last_print <- Sys.time()
    for (rep_id in start_rep:R) {
      results_list[[rep_id]] <- tryCatch(
        run_one_replication(N, TT, cfg, rep_id),
        error = function(e) rep(NA, 21))

      if (as.numeric(difftime(Sys.time(), last_print, units = "secs")) > 5) {
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        done <- rep_id - start_rep + 1
        eta <- elapsed * (R - rep_id) / done
        cat(sprintf("[%d/%d ETA:%.0fm] ", rep_id, R, eta))
        last_print <- Sys.time()
      }

      if (rep_id %% 25 == 0) {
        partial_data <- results_list[1:rep_id]
        save(partial_data, file = partial_file)
      }
    }
  }

  results_mat <- do.call(rbind, results_list[1:R])
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  if (file.exists(partial_file)) file.remove(partial_file)

  if (!is.null(results_mat) && nrow(results_mat) > 0) {
    colnames(results_mat) <- cn21
    nv <- sum(is.finite(results_mat[, "QML_phi"]))
    cat(sprintf("done %.1fmin (%d/%d ok)\n", el, nv, R))
    sd <- summarise_mc(results_mat, cfg)
    sd$N <- N; sd$T <- TT; sd$Config <- cfg_name
    cell_result <- list(raw = results_mat, summary = sd)
    save(cell_result, file = save_file)
    print(sd)
  } else {
    cat(sprintf("FAIL %.1fmin\n", el))
  }
}

cat("✓ 00_Functions.R loaded (3-PARAM: phi, gamma1, sigma2_u)\n")
cat(sprintf("  sigma2_u = gamma1^2 + gamma2^2*(1-2/pi) + sigma2_e\n"))
cat(sprintf("  Baseline true sigma2_u = %.5f\n",
    (-0.15)^2 + 0.10^2 * sigma2_a + 0.05))
