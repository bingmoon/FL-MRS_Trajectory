# ==============================================================================
# Script 3: Advanced Metrics (Calibration, DCA, Bootstrap) for Top Journals
# ==============================================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(pROC)
  library(dcurves)    
  library(patchwork)  
})

cat("\n🚀 启动高级评估流水线 (Calibration & DCA & Bootstrap)...\n")

if (!exists("test_prob") | !exists("test_data")) {
  stop("❌ 错误：未找到测试集数据！请先运行 01_Main_Pipeline.R。")
}

adv_df <- data.frame(
  Prob = test_prob,
  Outcome = ifelse(test_data$Group == "Colorectal Cancer", 1, 0)
)

# --- 模块 A: Calibration Curve 数据准备 ---
adv_df <- adv_df %>% mutate(Bin = ntile(Prob, 5)) 
calib_data <- adv_df %>%
  group_by(Bin) %>%
  summarise(
    Mean_Pred = mean(Prob),
    Observed_Rate = mean(Outcome),
    SE = sqrt(Observed_Rate * (1 - Observed_Rate) / n()),
    .groups = "drop"
  )

# --- 模块 B: DCA 数据准备 ---
dca_res <- dca(Outcome ~ Prob, data = adv_df, thresholds = seq(0, 1, by = 0.01))
dca_plot_data <- as_tibble(dca_res)

# --- 顶刊级联合绘图 ---
p_calib <- ggplot(calib_data, aes(x = Mean_Pred, y = Observed_Rate)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 1.2) +
  geom_errorbar(aes(ymin = Observed_Rate - SE, ymax = Observed_Rate + SE), width = 0.05, color = "#2C3E50", linewidth = 1.2) +
  geom_point(shape = 21, size = 5.5, fill = "#E74C3C", color = "black", stroke = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#2C3E50", linetype = "dotted", linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.02, 0.02)) +
  theme_classic(base_size = 15) +
  labs(title = "A. Calibration Curve", x = "Predicted Malignancy Probability", y = "Observed Proportion of CRC") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", color = "black"),
        axis.text = element_text(color = "black", size = 13),
        axis.line = element_line(color = "black", linewidth = 1.2),
        axis.ticks = element_line(color = "black", linewidth = 1.2))

p_dca <- ggplot(dca_plot_data, aes(x = threshold, y = net_benefit, color = label)) +
  geom_line(linewidth = 1.8) +
  coord_cartesian(ylim = c(-0.05, max(dca_plot_data$net_benefit, na.rm=TRUE) * 1.15)) +
  scale_color_manual(values = c("Prob" = "#D32F2F", "Treat All" = "grey50", "Treat None" = "black"),
                     labels = c("Treat None", "Treat All", "Fecal Lipidomic Risk Score")) +
  theme_classic(base_size = 15) +
  labs(title = "B. Decision Curve Analysis", x = "Threshold Probability", y = "Clinical Net Benefit", color = "") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", color = "black"),
        axis.text = element_text(color = "black", size = 13),
        axis.line = element_line(color = "black", linewidth = 1.2),
        axis.ticks = element_line(color = "black", linewidth = 1.2),
        legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8), color = "black", linewidth = 0.5),
        legend.text = element_text(size = 13, face = "bold"))

final_adv_plot <- p_calib + p_dca + plot_layout(ncol = 2)
ggsave("Figure_S1_Advanced_Metrics.pdf", final_adv_plot, width = 11, height = 5.5, dpi = 300, device = cairo_pdf)
cat("✅ 顶刊级图表已完美输出为: Figure_S1_Advanced_Metrics.pdf\n")

# --- 模块 C: Bootstrap 1000次 ---
set.seed(2026)
n_boot <- 1000
boot_cutoffs <- numeric(n_boot)
boot_aucs <- numeric(n_boot)
cases <- which(adv_df$Outcome == 1)
controls <- which(adv_df$Outcome == 0)

for (i in 1:n_boot) {
  boot_idx <- c(sample(cases, replace = TRUE), sample(controls, replace = TRUE))
  boot_sample <- adv_df[boot_idx, ]
  suppressMessages({
    boot_roc <- roc(boot_sample$Outcome, boot_sample$Prob, direction = "<", quiet = TRUE)
    boot_aucs[i] <- auc(boot_roc)
    boot_cutoffs[i] <- coords(boot_roc, "best", ret="threshold", transpose = FALSE)$threshold[1]
  })
}
cat(sprintf("✅ Bootstrap Median Cutoff: %.3f (95%% CI: %.3f - %.3f)\n", median(boot_cutoffs), quantile(boot_cutoffs, 0.025), quantile(boot_cutoffs, 0.975)))
cat(sprintf("✅ Bootstrap Median AUC: %.3f\n", median(boot_aucs)))
