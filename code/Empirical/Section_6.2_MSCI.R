################################################################################
#  APSV_App3_MSCI.R
#  ─────────────────────────────────────────────────────────────────────────────
#  APPLICATION 3: MSCI Developed Market Indices
#
#  PRIOR STUDY:  Poignard & Asai (2024), "Factor Multivariate Stochastic
#                Volatility Models of High Dimension", arXiv 2406.19033
#  MODEL USED:   fMSV (factor MSV via sparse PCA, symmetric volatility)
#                Same MSCI dataset (12/31/1998 – 03/12/2018) used in that paper
#  YOUR ADVANCE: Panel APSV replaces PCA with CCE; adds asymmetric leverage
#
#  EMPIRICAL PROBLEM SOLVED:
#    Poignard & Asai (2024) fMSV:
#      (1) Requires pre-specifying number of factors (PCA-based extraction)
#      (2) Imposes SYMMETRIC volatility — no leverage effect per country
#      (3) Models covariance matrices; individual volatility precision sacrificed
#    Panel APSV (this paper):
#      (1) CCE automatically handles factor uncertainty (no r selection needed)
#      (2) Heterogeneous γ₁ per country — leverage differs UK vs Japan vs Korea
#      (3) Individually targets log-volatility precision → better VaR/ES
#
#  DATA SOURCE NOTES:
#    Poignard & Asai (2024) use MSCI indices directly from msci.com.
#    This script uses iShares MSCI single-country ETFs as publicly-available
#    proxies — these track MSCI indices with tracking error < 0.2%.
#    For exact replication, replace with msci.com or Bloomberg data if available.
#
#  DATA: Daily ETF returns (MSCI proxies), Jan 1999 – Dec 2024
#        N = 10–14 countries (availability varies)
#        T ≈ 6,500 trading days
#        Source: Yahoo Finance (free, iShares ETFs)
################################################################################

rm(list = ls())
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
source("../00_Functions.R")
source("../00_Empirical_Utilities.R")

# ─────────────────────────────────────────────────────────────────────────────
# 1. DOWNLOAD DATA — iShares MSCI single-country ETFs
# ─────────────────────────────────────────────────────────────────────────────
tickers <- c(
  "EWA",   # MSCI Australia
  "EWC",   # MSCI Canada
  "EWD",   # MSCI Sweden
  "EWG",   # MSCI Germany
  "EWH",   # MSCI Hong Kong
  "EWJ",   # MSCI Japan
  "EWL",   # MSCI Switzerland
  "EWN",   # MSCI Netherlands
  "EWP",   # MSCI Spain
  "EWQ",   # MSCI France
  "EWU",   # MSCI United Kingdom
  "EWY",   # MSCI South Korea
  "EWT",   # MSCI Taiwan
  "EWI"    # MSCI Italy
)
names_ <- c("Australia","Canada","Sweden","Germany","HongKong",
            "Japan","Switzerland","Netherlands","Spain","France",
            "UK","SouthKorea","Taiwan","Italy")

raw <- download_returns(tickers, names_,
                        date_start = "1999-01-04",
                        date_end   = "2024-12-31",
                        min_obs    = 0.90)
saveRDS(raw, "data/app3_raw.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 2. APSV ESTIMATION
# ─────────────────────────────────────────────────────────────────────────────
res <- run_apsv(raw$ret, "App3: MSCI Developed Market Indices",
                init_phi  = 0.96,
                init_g1   = -0.10,
                init_s2u  = 0.06)
saveRDS(res, "data/app3_results.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 3. RESULTS TABLE
# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── App 3: MSCI MG Estimates ───\n")
cat(sprintf("%-14s  %10s  %10s  %10s\n", "Country","phi","gamma1","sigma2_u"))
cat(paste(rep("-",50), collapse=""), "\n")
for (i in seq_len(res$N)) {
  cat(sprintf("%-14s  %+10.5f  %+10.5f  %10.5f\n",
              rownames(res$MG_mat)[i],
              res$MG_mat[i,"phi"], res$MG_mat[i,"gamma1"],
              res$MG_mat[i,"sigma2_u"]))
}
cat(paste(rep("-",50), collapse=""), "\n")
cat(sprintf("%-14s  %+10.5f  %+10.5f  %10.5f  (HPJ pooled)\n",
            "POOLED",res$HPJ$theta[1],res$HPJ$theta[2],res$HPJ$theta[3]))
cat(sprintf("\nHalf-life: %.1f d  SNR: %.5f  Factors: %d\n",
            res$half_life, res$SNR, res$n_factors))
cat(sprintf("Hausman: chi2=%+.2f  p=%.4f [%s]\n",
            res$Hausman$stat, res$Hausman$pval,
            ifelse(res$Hausman$pval<0.01,"REJECT ***","not rejected")))

# ─────────────────────────────────────────────────────────────────────────────
# 4. FIGURES
# ─────────────────────────────────────────────────────────────────────────────
mg_df <- as.data.frame(res$MG_mat) |>
  rownames_to_column("country") |>
  filter(!is.na(gamma1)) |>
  mutate(country = factor(country, levels = country[order(gamma1)]))

# Panel (a): lollipop chart of γ₁ by country
p_a <- ggplot(mg_df, aes(x = country, y = gamma1, colour = gamma1)) +
  geom_segment(aes(xend = country, yend = 0), colour = "grey75", linewidth = 1) +
  geom_point(size = 5) +
  scale_colour_gradient2(low = "#E63946", mid = "#FFB703", high = "#4361EE",
                         midpoint = median(mg_df$gamma1, na.rm = TRUE),
                         name = expression(hat(gamma)[1])) +
  geom_hline(yintercept = res$HPJ$theta[2], linetype = "dashed",
             colour = "#E63946", linewidth = 0.8) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  labs(x = NULL, y = expression(hat(gamma)[1]),
       subtitle = "(a) Country-level leverage: MSCI developed markets") +
  coord_flip() + theme_apsv

# Panel (b): variance explained by factors
var_df <- data.frame(
  factor  = paste0("F", 1:min(6, length(res$var_exp))),
  marg    = c(res$var_exp[1],
              diff(res$var_exp[1:min(6, length(res$var_exp))])) * 100
)
p_b <- ggplot(var_df, aes(x = factor, y = marg)) +
  geom_col(fill = "#1F3864", colour = "grey25", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", marg), y = marg + 0.4),
            size = 3.2, fontface = "bold") +
  labs(x = "Common factor", y = "Variance explained (%)",
       subtitle = "(b) Variance explained by CCE-selected factors") +
  theme_apsv

# Panel (c): persistence by country
mg_df2 <- mg_df |>
  mutate(country = factor(country, levels = country[order(phi)]))
p_c <- ggplot(mg_df2, aes(x = country, y = phi)) +
  geom_col(fill = "#4361EE", colour = "grey25", linewidth = 0.3, width = 0.65) +
  geom_hline(yintercept = res$HPJ$theta[1], linetype = "dashed",
             colour = "#E63946", linewidth = 0.8) +
  labs(x = NULL, y = expression(hat(phi) ~ "(persistence)"),
       subtitle = expression("(c) Log-volatility persistence " * hat(phi))) +
  coord_flip() + theme_apsv

fig <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(
    title   = "Figure: App 3 — MSCI Developed Markets (APSV with CCE-QML)",
    subtitle = sprintf("N = %d countries, T = %d daily obs. Jan 1999–Dec 2024. Prior: Poignard & Asai (2024) fMSV.",
                       res$N, res$T),
    caption = "Notes: iShares MSCI single-country ETFs used as MSCI index proxies. Extends Poignard & Asai (2024) by adding leverage per country."
  )

ggsave("figures/App3_MSCI.pdf", fig, width = 12, height = 9, dpi = 300)
cat("\n✓ Figure saved: figures/App3_MSCI.pdf\n")

write.csv(
  data.frame(
    Country    = rownames(res$MG_mat),
    phi_MG     = round(res$MG_mat[,"phi"],     5),
    gamma1_MG  = round(res$MG_mat[,"gamma1"],  5),
    sigma2u_MG = round(res$MG_mat[,"sigma2_u"],5),
    phi_HPJ    = round(res$HPJ$theta[1], 5),
    gamma1_HPJ = round(res$HPJ$theta[2], 5),
    sigma2u_HPJ= round(res$HPJ$theta[3], 5)
  ),
  "results/App3_MSCI_estimates.csv", row.names = FALSE
)
cat("✓ Estimates saved: results/App3_MSCI_estimates.csv\n")
