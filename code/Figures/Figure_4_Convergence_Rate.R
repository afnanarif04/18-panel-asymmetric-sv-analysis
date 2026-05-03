################################################################################
#  Figure_4_Convergence_Rate.R
#  ───────────────────────────────────────────────────────────────────────────
#  Generates Figure 4 of the paper:
#  Verifies the predicted 1/sqrt(NT) convergence rate of the CCE-QML estimator
#  by plotting RMSE × sqrt(NT) — should be approximately constant if the rate
#  is correct.
#  Output: figures/Figure4_Convergence_Rate.pdf
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

# Compute scaled RMSE = sqrt(NT) * RMSE
rows <- list()
for (key in names(results)) {
  cell <- results[[key]]
  scaled_rmse <- sqrt(cell$N * cell$T) * cell$rmse_qml
  rows[[length(rows)+1]] <- data.frame(
    N = cell$N, T = cell$T, NT = cell$N * cell$T,
    param = c("phi","gamma1","sigma2_u"),
    scaled_rmse = scaled_rmse,
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)
df$param <- factor(df$param, levels = c("phi","gamma1","sigma2_u"))

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#1F3864"),
        strip.text       = element_text(colour = "white", face = "bold"),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 12, colour = "#1F3864"))

p <- ggplot(df, aes(x = NT, y = scaled_rmse, colour = factor(N))) +
  geom_line(linewidth = 0.5, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  scale_x_log10() +
  scale_colour_brewer(palette = "Set1", name = "Cross-section size N") +
  labs(x = "NT (log scale)", y = expression(sqrt(NT) %*% RMSE),
       title = "Figure 4: Verification of 1/sqrt(NT) Convergence Rate",
       subtitle = "If the predicted rate is correct, scaled RMSE should be approximately constant") +
  theme_paper

ggsave("figures/Figure4_Convergence_Rate.pdf", p,
       width = 12, height = 5.5, dpi = 300)
cat("Figure 4 saved: figures/Figure4_Convergence_Rate.pdf\n")
