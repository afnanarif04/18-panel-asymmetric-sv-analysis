################################################################################
#  Figure_2_CCE_vs_Naive.R
#  ───────────────────────────────────────────────────────────────────────────
#  Generates Figure 2 of the paper:
#  Demonstrates the dramatic improvement of CCE-QML over the Naive (no-CCE)
#  estimator. Highlights the sigma2_u channel where the difference is largest.
#  Output: figures/Figure2_CCE_vs_Naive.pdf
#
#  REQUIRES: Results/MC_Baseline.rds (run MonteCarlo/MC_01_Baseline.R first)
################################################################################

rm(list = ls())
if (requireNamespace("rstudioapi", quietly = TRUE))
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(dplyr)
library(tidyr)
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
  rows[[length(rows)+1]] <- data.frame(
    N = cell$N, T = cell$T,
    param = rep(c("phi","gamma1","sigma2_u"), 2),
    estimator = rep(c("CCE-QML (proposed)","Naive (no CCE)"), each = 3),
    rmse = c(cell$rmse_qml, cell$rmse_naive),
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)

df$param <- factor(df$param, levels = c("phi","gamma1","sigma2_u"),
                   labels = c("phi", "gamma1", "sigma2_u"))
df$N_lbl <- factor(df$N, levels = c(30,50,100,200,400),
                   labels = paste0("N=", c(30,50,100,200,400)))

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#1F3864"),
        strip.text       = element_text(colour = "white", face = "bold"),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 12, colour = "#1F3864"))

# Per-parameter panel
p <- ggplot(df, aes(x = T, y = rmse, colour = estimator, linetype = estimator,
                    group = interaction(estimator, N_lbl))) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_grid(param ~ N_lbl, scales = "free_y") +
  scale_colour_manual(values = c("CCE-QML (proposed)" = "#1F3864",
                                  "Naive (no CCE)"    = "#E63946"),
                      name = "Estimator") +
  scale_linetype_manual(values = c("solid","dashed"), guide = "none") +
  scale_x_log10(breaks = c(250,500,1000,2000)) +
  scale_y_log10() +
  labs(x = "Time series length T", y = "RMSE (log scale)",
       title = "Figure 2: CCE-QML vs Naive Estimator — RMSE Across Grid",
       subtitle = "Naive omits cross-sectional defactoring; CCE-QML applies it (Baseline configuration)") +
  theme_paper

ggsave("figures/Figure2_CCE_vs_Naive.pdf", p,
       width = 13, height = 8, dpi = 300)
cat("Figure 2 saved: figures/Figure2_CCE_vs_Naive.pdf\n")
