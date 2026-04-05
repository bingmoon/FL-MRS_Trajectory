# ==============================================================================
# Script 04: Feature Extraction (Table S4) & Stratified Bootstrap AUC CI
# Description: Generates values for manuscript blanks and exports Table S4.
# INSTRUCTIONS: Run this AFTER executing 'Main_Pipeline_ST003798.R'
# ==============================================================================

suppressMessages({
  library(glmnet)
  library(randomForest)
  library(pROC)
  library(dplyr)
})

cat("\n======================================================\n")
cat("🚀 正在提取手稿填空数据与生成 Table S4...\n")

# ------------------------------------------------------------------------------
# 任务 1: LASSO 特征数量提取
# ------------------------------------------------------------------------------
# selected_features 已在 Main_Pipeline 的 Step 3 中完美定义
feature_count <- length(selected_features)
cat("\n>>> [手稿填空 1] LASSO (lambda.1se) 筛选出的核心脂质数量为:", feature_count, "\n")

# ------------------------------------------------------------------------------
# 任务 2: 组装 Supplementary Table S4 (包含系数与 Gini 重要性)
# ------------------------------------------------------------------------------
# 1. 提取 LASSO 具体系数值
all_coefs <- as.matrix(predict(lasso_cv, type = "coef", s = best_lambda))
feature_coefs <- all_coefs[selected_features, 1]

# 2. 提取 Random Forest 的 MeanDecreaseGini (VIP重要性)
# rf_model 已在 Main_Pipeline 的 Step 4 中设定 importance = TRUE
rf_imp <- randomForest::importance(rf_model)
gini_scores <- rf_imp[selected_features, "MeanDecreaseGini"]

# 3. 合并为数据框并按照重要性排序
Table_S4 <- data.frame(
  Lipid_Feature = selected_features,
  LASSO_Coefficient = round(feature_coefs, 4),
  RF_Importance_Gini = round(gini_scores, 4)
) %>%
  arrange(desc(RF_Importance_Gini))

# 4. 导出为 CSV
write.csv(Table_S4, "Table_S4_Feature_Coefficient_VIP.csv", row.names = FALSE)
cat(">>> [成功] Supplementary Table S4 已生成并保存为: Table_S4_Feature_Coefficient_VIP.csv\n")

# ------------------------------------------------------------------------------
# 任务 3: 测试集 AUC 的 1000次分层 Bootstrap 95% 置信区间
# ------------------------------------------------------------------------------
cat("\n⏳ 正在运行 1000次分层 Bootstrap 计算 AUC 置信区间 (约需10-20秒)...\n")

set.seed(2026) # 保持与主流程相同的种子，保证绝对可复现
# roc_obj 已在 Main_Pipeline 的 Step 4 中完美定义
auc_ci_boot <- ci.auc(roc_obj, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)

auc_lower <- round(auc_ci_boot[1], 3)
auc_median <- round(auc_ci_boot[2], 3)
auc_upper <- round(auc_ci_boot[3], 3)

cat(sprintf(">>> [手稿填空 2] 1000次 Bootstrap AUC 结果为: %.3f (95%% CI: %.3f - %.3f)\n", 
            auc_median, auc_lower, auc_upper))
cat("======================================================\n")
