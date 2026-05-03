################################################################################
#  Figure_1_Baseline_Performance.R
#  ───────────────────────────────────────────────────────────────────────────
#  Generates Figure 1 of the paper:
#  Bias and RMSE of CCE-QML estimator across N×T grid for Baseline configuration
#  Three parameters: phi, gamma1, sigma2_u
#  Output: figures/Figure1_Baseline_Performance.pdf
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

# Load Baseline results
res_path <- "../MonteCarlo/Results/MC_Baseline.rds"
if (!file.exists(res_path))
  stop(paste0("Results not found at: ", res_path,
              "\nRun MonteCarlo/MC_01_Baseline.R first."))
results <- readRDS(res_path)

# Build long data frame for plotting
build_df <- function(results) {
  rows <- list()
  for (key in names(results)) {
    cell <- results[[key]]
    rows[[length(rows)+1]] <- data.frame(
      N = cell$N, T = cell$T,
      param = c("phi","gamma1","sigma2_u"),
      bias  = cell$bias_qml,
      rmse  = cell$rmse_qml,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
df <- build_df(results)

# Format
df$param <- factor(df$param, levels = c("phi","gamma1","sigma2_u"),
                   labels = c("phi (persistence)",
                              "gamma1 (leverage)",
                              "sigma2_u (transition variance)"))
df$T_lbl <- factor(df$T, levels = c(250,500,1000,2000),
                   labels = c("T=250","T=500","T=1000","T=2000"))

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#1F3864"),
        strip.text       = element_text(colour = "white", face = "bold"),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 12, colour = "#1F3864"))

# Panel A: Bias
p_bias <- ggplot(df, aes(x = N, y = bias, colour = T_lbl, group = T_lbl)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("#E63946","#F77F00","#16a34a","#1F3864"),
                      name = "Time series length") +
  scale_x_log10(breaks = c(30,50,100,200,400)) +
  labs(x = "Cross-section size N", y = "Bias",
       title = "(a) CCE-QML bias across N x T grid (Baseline)") +
  theme_paper

# Panel B: RMSE
p_rmse <- ggplot(df, aes(x = N, y = rmse, colour = T_lbl, group = T_lbl)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("#E63946","#F77F00","#16a34a","#1F3864"),
                      name = "Time series length") +
  scale_x_log10(breaks = c(30,50,100,200,400)) +
  scale_y_log10() +
  labs(x = "Cross-section size N", y = "RMSE (log scale)",
       title = "(b) CCE-QML RMSE across N x T grid (Baseline)") +
  theme_paper

fig <- p_bias / p_rmse +
  plot_layout(heights = c(1,1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures/Figure1_Baseline_Performance.pdf",
       fig, width = 12, height = 8, dpi = 300)
cat("Figure 1 saved: figures/Figure1_Baseline_Performance.pdf\n")
