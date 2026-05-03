################################################################################
#  Figure_3_HPJ_Bias_Correction.R
#  ───────────────────────────────────────────────────────────────────────────
#  Generates Figure 3 of the paper:
#  Bias of CCE-QML vs HPJ across N×T grid. Demonstrates that HPJ exactly
#  cancels the leading O(1/T) Nickell bias for phi.
#  Output: figures/Figure3_HPJ_Bias_Correction.pdf
#
#  REQUIRES: Results/MC_Baseline.rds (run MonteCarlo/MC_01_Baseline.R first)
################################################################################

rm(list = ls())
if (requireNamespace("rstudioapi", quietly = TRUE))
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(dplyr)
library(patchwork)

dir.create("figures", showWarnings = FALSE)

res_path <- "../MonteCarlo/Results/MC_Baseline.rds"
if (!file.exists(res_path))
  stop(paste0("Results not found at: ", res_path))
results <- readRDS(res_path)

# Build comparison data frame
rows <- list()
for (key in names(results)) {
  cell <- results[[key]]
  # Compute bias for HPJ
  est_hpj_phi <- mean(cell$estimates$phi_hpj, na.rm = TRUE) - cell$true["phi"]
  est_hpj_g1  <- mean(cell$estimates$g1_hpj, na.rm = TRUE) - cell$true["g1"]
  est_hpj_s2u <- mean(cell$estimates$s2u_hpj, na.rm = TRUE) - cell$true["s2u"]
  rows[[length(rows)+1]] <- data.frame(
    N = cell$N, T = cell$T,
    param = rep(c("phi","gamma1","sigma2_u"), 2),
    estimator = rep(c("CCE-QML (uncorrected)","HPJ (bias-corrected)"), each = 3),
    bias = c(cell$bias_qml,
             c(est_hpj_phi, est_hpj_g1, est_hpj_s2u)),
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)

df$param <- factor(df$param, levels = c("phi","gamma1","sigma2_u"))
df$N_lbl <- factor(df$N, levels = c(30,50,100,200,400),
                   labels = paste0("N=", c(30,50,100,200,400)))

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#1F3864"),
        strip.text       = element_text(colour = "white", face = "bold"),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 12, colour = "#1F3864"))

p <- ggplot(df, aes(x = T, y = bias, colour = estimator, group = interaction(estimator, N_lbl))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_grid(param ~ N_lbl, scales = "free_y") +
  scale_colour_manual(values = c("CCE-QML (uncorrected)" = "#E63946",
                                  "HPJ (bias-corrected)" = "#1F3864"),
                      name = "Estimator") +
  scale_x_log10(breaks = c(250,500,1000,2000)) +
  labs(x = "Time series length T", y = "Bias",
       title = "Figure 3: HPJ Bias Correction — Cancels Leading O(1/T) Nickell Bias",
       subtitle = "Baseline configuration; the HPJ bias is centred on zero across all grid cells") +
  theme_paper

ggsave("figures/Figure3_HPJ_Bias_Correction.pdf", p,
       width = 13, height = 8, dpi = 300)
cat("Figure 3 saved: figures/Figure3_HPJ_Bias_Correction.pdf\n")
