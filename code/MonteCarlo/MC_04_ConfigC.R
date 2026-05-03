################################################################################
#  MC_04_ConfigC.R
#  Monte Carlo simulation — Configuration C: STRONG LEVERAGE, NO SIZE
#  Configuration: φ = 0.95, γ₁ = -0.20, γ₂ = 0.00, σ²_e = 0.05
#  Implied σ²_u = γ₁² + γ₂²(1−2/π) + σ²_e = 0.09000
#  Grid: N ∈ {30,50,100,200,400} × T ∈ {250,500,1000,2000} = 20 cells
#  Replications per cell: R = 500
#  Estimators reported: CCE-QML (proposed), Naive (no defactoring), HPJ (bias-corrected)
################################################################################

rm(list = ls())
if (requireNamespace("rstudioapi", quietly = TRUE))
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../00_Functions.R")

CONFIG_NAME <- "C"
cfg <- configs[[CONFIG_NAME]]
phi0 <- cfg$phi; g1_0 <- cfg$gamma1; g2_0 <- cfg$gamma2; s2e_0 <- cfg$sigma2_e
s2u_true <- g1_0^2 + g2_0^2 * (1 - 2/pi) + s2e_0

cat("══════════════════════════════════════════════════════════════\n")
cat(sprintf(" Monte Carlo: Configuration %s\n", CONFIG_NAME))
cat(sprintf(" True parameters: phi=%.2f gamma1=%.2f gamma2=%.2f sigma2_e=%.2f\n",
            phi0, g1_0, g2_0, s2e_0))
cat(sprintf(" Implied sigma2_u (identifiable) = %.5f\n", s2u_true))
cat(sprintf(" R = %d replications per cell, %d cells total\n",
            R, length(N_grid)*length(T_grid)))
cat("══════════════════════════════════════════════════════════════\n\n")

dir.create("Results", showWarnings = FALSE)
results_all <- list()
t_start <- Sys.time()

for (N in N_grid) {
  for (TT in T_grid) {
    cat(sprintf("\n[%s] N = %d, T = %d  (%s)\n",
                CONFIG_NAME, N, TT,
                format(Sys.time(), "%H:%M:%S")))

    cell_estimates <- matrix(NA, R, 9,
      dimnames = list(NULL, c("phi_qml","g1_qml","s2u_qml",
                               "phi_naive","g1_naive","s2u_naive",
                               "phi_hpj","g1_hpj","s2u_hpj")))

    for (r in 1:R) {
      seed_r <- 1000000 + match(N, N_grid) * 10000 +
                match(TT, T_grid) * 1000 + r
      sim <- generate_apsv_data(N, TT, phi0, g1_0, g2_0, s2e_0,
                                 r = r_factors, rho_f = rho_f, seed = seed_r)

      # 1. Build log-squared proxy
      y_dm   <- sweep(sim$y, 1, rowMeans(sim$y), "-")
      y_dm[y_dm == 0] <- 1e-6
      x_mat  <- log(y_dm^2) - E_log_chi2

      # 2. CCE-QML (proposed)
      cce <- tryCatch(cce_defactor(x_mat), error = function(e) NULL)
      if (!is.null(cce)) {
        y_eff <- y_dm[, (cce$p_T + 1):TT, drop = FALSE]
        est <- tryCatch(
          pooled_qml(cce$x_star, y_eff,
                     init = c(0.90, -0.10, 0.08),
                     h_level_mat = cce$cce_fitted, compute_se = FALSE),
          error = function(e) list(theta = rep(NA, 3)))
        cell_estimates[r, 1:3] <- est$theta
      }

      # 3. Naive (no CCE) — h_level = 0 for all units
      h_zero <- matrix(0, N, TT)
      est_n <- tryCatch(
        pooled_qml(x_mat, y_dm, init = c(0.90, -0.10, 0.08),
                   h_level_mat = h_zero, compute_se = FALSE),
        error = function(e) list(theta = rep(NA, 3)))
      cell_estimates[r, 4:6] <- est_n$theta

      # 4. HPJ bias correction
      if (!is.null(cce)) {
        est_h <- tryCatch(
          hpj_correction(cce$x_star, y_eff,
                         init = c(0.90, -0.10, 0.08),
                         h_level_mat = cce$cce_fitted),
          error = function(e) list(theta = rep(NA, 3)))
        cell_estimates[r, 7:9] <- est_h$theta
      }

      if (r %% 50 == 0)
        cat(sprintf("  rep %d/%d  done\n", r, R))
    }

    # Cell summary
    cell_df <- as.data.frame(cell_estimates)
    summary_cell <- list(
      config = CONFIG_NAME, N = N, T = TT,
      true   = c(phi = phi0, g1 = g1_0, s2u = s2u_true),
      bias_qml   = c(mean(cell_df$phi_qml,  na.rm=T) - phi0,
                     mean(cell_df$g1_qml,   na.rm=T) - g1_0,
                     mean(cell_df$s2u_qml,  na.rm=T) - s2u_true),
      rmse_qml   = c(sqrt(mean((cell_df$phi_qml  - phi0)^2,    na.rm=T)),
                     sqrt(mean((cell_df$g1_qml   - g1_0)^2,    na.rm=T)),
                     sqrt(mean((cell_df$s2u_qml  - s2u_true)^2,na.rm=T))),
      rmse_naive = c(sqrt(mean((cell_df$phi_naive  - phi0)^2,    na.rm=T)),
                     sqrt(mean((cell_df$g1_naive   - g1_0)^2,    na.rm=T)),
                     sqrt(mean((cell_df$s2u_naive  - s2u_true)^2,na.rm=T))),
      rmse_hpj   = c(sqrt(mean((cell_df$phi_hpj  - phi0)^2,    na.rm=T)),
                     sqrt(mean((cell_df$g1_hpj   - g1_0)^2,    na.rm=T)),
                     sqrt(mean((cell_df$s2u_hpj  - s2u_true)^2,na.rm=T))),
      estimates = cell_df
    )
    results_all[[paste0("N",N,"_T",TT)]] <- summary_cell

    cat(sprintf("  CCE-QML bias:  phi=%+.5f  g1=%+.5f  s2u=%+.5f\n",
                summary_cell$bias_qml[1], summary_cell$bias_qml[2],
                summary_cell$bias_qml[3]))
    cat(sprintf("  CCE-QML RMSE:  phi=%.5f  g1=%.5f  s2u=%.5f\n",
                summary_cell$rmse_qml[1], summary_cell$rmse_qml[2],
                summary_cell$rmse_qml[3]))
    cat(sprintf("  Naive   RMSE:  phi=%.5f  g1=%.5f  s2u=%.5f  [ratio %5.1fx]\n",
                summary_cell$rmse_naive[1], summary_cell$rmse_naive[2],
                summary_cell$rmse_naive[3],
                summary_cell$rmse_naive[3] / max(summary_cell$rmse_qml[3], 1e-9)))
    cat(sprintf("  HPJ     RMSE:  phi=%.5f  g1=%.5f  s2u=%.5f\n",
                summary_cell$rmse_hpj[1], summary_cell$rmse_hpj[2],
                summary_cell$rmse_hpj[3]))
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nTotal runtime: %.1f minutes\n", elapsed))

saveRDS(results_all, sprintf("Results/MC_%s.rds", CONFIG_NAME))
cat(sprintf("\n✓ Results saved: Results/MC_%s.rds\n", CONFIG_NAME))
