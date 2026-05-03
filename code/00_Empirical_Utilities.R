################################################################################
#  00_Empirical_Utilities.R
#  ─────────────────────────────────────────────────────────────────────────────
#  Shared helper functions for all 5 empirical applications.
#  source("00_Empirical_Utilities.R") at the top of each app script.
#
#  REQUIRES: 00_Functions.R (Kalman filter + CCE) in same directory
#            install.packages(c("quantmod","xts","zoo","tidyverse",
#                               "moments","tseries","patchwork","viridis"))
################################################################################

library(quantmod)
library(xts)
library(zoo)
library(tidyverse)
library(moments)
library(patchwork)
library(viridis)

dir.create("data",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

E_log_chi2 <- digamma(0.5) + log(2)
V_log_chi2 <- pi^2 / 2

# ─────────────────────────────────────────────────────────────────────────────
# DATA DOWNLOAD
# Fixed version: aligns via merge.xts on Date index (avoids xts() error)
# ─────────────────────────────────────────────────────────────────────────────
download_returns <- function(tickers, names,
                              date_start = "1999-01-04",
                              date_end   = "2024-12-31",
                              min_obs    = 0.95) {
  cat(sprintf("  Downloading %d series from Yahoo Finance...\n", length(tickers)))
  px_list <- list()

  for (i in seq_along(tickers)) {
    px <- tryCatch(
      suppressWarnings(
        getSymbols(tickers[i], src = "yahoo",
                   from = date_start, to = date_end,
                   auto.assign = FALSE)
      ),
      error = function(e) {
        cat(sprintf("    FAIL: %s (%s)\n", tickers[i], conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(px) && nrow(px) > 50) {
      adj <- tryCatch(Ad(px), error = function(e) Cl(px))
      colnames(adj) <- names[i]
      px_list[[names[i]]] <- adj
    }
    Sys.sleep(0.4)   # be polite to Yahoo rate limiter
  }

  if (length(px_list) == 0) stop("No series downloaded. Check internet / tickers.")

  # Merge all xts objects on their shared Date index
  # merge.xts handles misaligned dates correctly
  merged <- Reduce(function(a, b) merge(a, b, join = "outer"), px_list)
  merged <- na.locf(merged, maxgap = 5)    # forward-fill short gaps
  merged <- na.locf(merged, fromLast = TRUE, maxgap = 5)  # back-fill at start

  # Log-returns in percent: r_t = 100*(log P_t - log P_{t-1})
  ret_xts <- diff(log(merged)) * 100
  ret_xts <- ret_xts[-1, ]

  # Keep columns with sufficient coverage
  ok <- colSums(!is.na(coredata(ret_xts))) >= min_obs * nrow(ret_xts)
  ret_xts <- ret_xts[, ok, drop = FALSE]

  # Convert to plain matrix; remove any residual NA rows/cols
  ret_mat <- coredata(ret_xts)
  rownames(ret_mat) <- as.character(index(ret_xts))
  good_cols <- colSums(is.na(ret_mat)) == 0
  ret_mat   <- ret_mat[, good_cols, drop = FALSE]

  # Replace exact zeros (can appear from fill) with tiny value
  ret_mat[ret_mat == 0] <- 1e-8

  cat(sprintf("  Panel: N = %d, T = %d  [%s to %s]\n",
              ncol(ret_mat), nrow(ret_mat),
              rownames(ret_mat)[1], rownames(ret_mat)[nrow(ret_mat)]))
  list(ret   = ret_mat,
       dates = as.Date(rownames(ret_mat)),
       N     = ncol(ret_mat),
       T     = nrow(ret_mat))
}

# ─────────────────────────────────────────────────────────────────────────────
# BAI-NG (2002) IC1 — number of common factors
# ─────────────────────────────────────────────────────────────────────────────
bai_ng_ic1 <- function(x_mat, r_max = 8) {
  N <- nrow(x_mat); TT <- ncol(x_mat); NT <- N * TT
  # Standardise rows
  x_c <- x_mat - rowMeans(x_mat)
  x_s <- x_c / pmax(apply(x_c, 1, sd), 1e-10)
  V   <- t(x_s) %*% x_s / NT
  eig <- eigen(V, symmetric = TRUE)
  lam <- pmax(eig$values, 0)
  Fh  <- eig$vectors[, seq_len(r_max), drop = FALSE]

  ic <- numeric(r_max + 1)
  for (k in 0:r_max) {
    Vk <- if (k == 0) {
      sum(x_s^2) / NT
    } else {
      Fk <- Fh[, 1:k, drop = FALSE]
      Lk <- x_s %*% Fk / TT
      sum((x_s - Lk %*% t(Fk))^2) / NT
    }
    pen     <- k * (N + TT) / NT * log(NT / (N + TT))
    ic[k+1] <- log(max(Vk, 1e-15)) + pen
  }
  r_hat <- which.min(ic) - 1
  list(r_hat        = r_hat,
       ic           = ic,
       var_explained = cumsum(lam[1:r_max]) / sum(lam))
}

# ─────────────────────────────────────────────────────────────────────────────
# HAUSMAN TEST for slope homogeneity (PMG vs MG)
# ─────────────────────────────────────────────────────────────────────────────
hausman_test <- function(qml_est, mg_mat, N) {
  mg_mean  <- colMeans(mg_mat, na.rm = TRUE)
  mg_se    <- apply(mg_mat, 2, sd, na.rm = TRUE) / sqrt(max(N - 1, 1))
  diff_v   <- mg_mean - unname(qml_est$theta)
  V_mg     <- diag(mg_se^2)
  V_qml    <- diag((unname(qml_est$se))^2 / N)
  V_diff   <- V_mg - V_qml + diag(1e-14, 3)
  # Force PSD
  ev       <- eigen(V_diff, symmetric = TRUE)
  ev$values<- pmax(ev$values, 1e-12)
  V_inv    <- ev$vectors %*% diag(1/ev$values) %*% t(ev$vectors)
  stat     <- as.numeric(t(diff_v) %*% V_inv %*% diff_v)
  stat     <- max(stat, 0)
  pval     <- pchisq(stat, df = 3, lower.tail = FALSE)
  list(stat = stat, pval = pval, df = 3)
}

# ─────────────────────────────────────────────────────────────────────────────
# FULL APSV PIPELINE for one dataset
# ─────────────────────────────────────────────────────────────────────────────
run_apsv <- function(ret_mat, app_name,
                     init_phi   = 0.95,
                     init_g1    = -0.10,
                     init_s2u   = 0.07) {
  init <- c(init_phi, init_g1, init_s2u)
  N    <- ncol(ret_mat)
  TT   <- nrow(ret_mat)

  cat(sprintf("\n┌─ %s ─┐\n", app_name))
  cat(sprintf("│  N = %d, T = %d\n", N, TT))

  # Demean, build log-squared proxy
  ret_dm  <- sweep(ret_mat, 2, colMeans(ret_mat, na.rm = TRUE), "-")
  y_mat   <- t(ret_dm)   # N × T
  y_mat[y_mat == 0] <- 1e-6
  x_mat   <- log(y_mat^2) - E_log_chi2

  # CCE defactoring
  p_T   <- floor(TT^(1/3))
  cce   <- cce_defactor(x_mat, p_T)
  x_star<- cce$x_star
  cce_f <- cce$cce_fitted
  y_eff <- y_mat[, (p_T + 1):TT, drop = FALSE]
  T_eff <- ncol(x_star)

  # Bai-Ng factor selection
  bn <- bai_ng_ic1(x_mat)
  cat(sprintf("│  Bai-Ng IC1: r_hat = %d\n", bn$r_hat))

  # Pooled CCE-QML
  cat("│  QML ... ")
  est_qml <- tryCatch(
    pooled_qml(x_star, y_eff, init, cce_f, compute_se = TRUE),
    error = function(e) {
      cat(sprintf("FAILED (%s)\n", conditionMessage(e)))
      list(theta = rep(NA, 3), se = rep(NA, 3), converged = FALSE)
    }
  )
  cat(sprintf("phi=%.4f  g1=%+.4f  s2u=%.5f\n",
              est_qml$theta[1], est_qml$theta[2], est_qml$theta[3]))

  # HPJ bias correction
  cat("│  HPJ ... ")
  est_hpj <- tryCatch(
    hpj_correction(x_star, y_eff, init, cce_f),
    error = function(e) list(theta = rep(NA, 3))
  )
  if (any(is.na(est_hpj$theta))) est_hpj$theta <- est_qml$theta  # fallback
  cat(sprintf("phi=%.4f  g1=%+.4f  s2u=%.5f\n",
              est_hpj$theta[1], est_hpj$theta[2], est_hpj$theta[3]))

  # MG estimation
  cat("│  MG (unit-by-unit) ... ")
  nm  <- colnames(ret_mat)
  mg_mat <- matrix(NA, N, 3,
                   dimnames = list(nm, c("phi","gamma1","sigma2_u")))
  for (i in seq_len(N)) {
    e_i <- tryCatch(
      pooled_qml(x_star[i,,drop=F], y_eff[i,,drop=F], init,
                 h_level_mat = cce_f[i,,drop=F], compute_se = FALSE),
      error = function(e) list(theta = rep(NA_real_, 3))
    )
    mg_mat[i,] <- e_i$theta
  }
  ok_mg   <- sum(complete.cases(mg_mat))
  mg_mean <- colMeans(mg_mat, na.rm = TRUE)
  mg_sd   <- apply(mg_mat, 2, sd, na.rm = TRUE)
  cat(sprintf("%d/%d units OK\n", ok_mg, N))

  # Hausman test
  hausman <- tryCatch(
    hausman_test(est_qml, mg_mat, N),
    error = function(e) list(stat = NA, pval = NA, df = 3)
  )

  # Derived quantities
  phi_hpj   <- est_hpj$theta[1]
  snr       <- est_hpj$theta[3] / (est_hpj$theta[3] + V_log_chi2)
  half_life <- if (is.finite(phi_hpj) && abs(phi_hpj) < 1 && phi_hpj > 0)
                 -log(2) / log(phi_hpj) else NA_real_

  cat(sprintf("│  Half-life = %.1f days  SNR = %.5f\n", half_life, snr))
  cat(sprintf("│  Hausman chi2(%d) = %.2f  p = %.4f  [%s]\n",
              hausman$df, hausman$stat, hausman$pval,
              ifelse(!is.na(hausman$pval) && hausman$pval < 0.01,
                     "REJECT H0 ***", "not reject")))
  cat(sprintf("└─────────────────────────────────────────────────┘\n"))

  list(app       = app_name,
       N         = N, T = TT, T_eff = T_eff, p_T = p_T,
       QML       = est_qml,
       HPJ       = est_hpj,
       MG_mat    = mg_mat,
       MG_mean   = mg_mean,
       MG_sd     = mg_sd,
       Hausman   = hausman,
       n_factors = bn$r_hat,
       var_exp   = bn$var_explained,
       SNR       = snr,
       half_life = half_life,
       x_star    = x_star,
       cce_f     = cce_f,
       y_eff     = y_eff)
}

# ─────────────────────────────────────────────────────────────────────────────
# SHARED PLOT THEME
# ─────────────────────────────────────────────────────────────────────────────
theme_apsv <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
    strip.background = element_rect(fill = "#1F3864"),
    strip.text       = element_text(colour = "white", face = "bold", size = 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(size = 12, face = "bold", colour = "#1F3864"),
    plot.subtitle    = element_text(size = 9, colour = "grey35"),
    plot.caption     = element_text(size = 8, colour = "grey55", hjust = 0),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9)
  )

cat("✓ 00_Empirical_Utilities.R loaded\n")
