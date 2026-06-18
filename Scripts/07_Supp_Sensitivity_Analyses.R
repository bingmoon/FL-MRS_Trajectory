###############################################################################
# 补充分析完整脚本
# 用途：投稿前的敏感性分析与稳健性验证
# 前提：必须先运行过主分析脚本，确保 rf_model, selected_features,
#        data_clean_3798, df_adenoma_full, meta_2787, pca_res, curve_df 等
#        对象已存在于 R 环境中。
# 作者：[BING]
# 日期：2026-06-18
###############################################################################

# 加载必要包
library(randomForest)
library(glmnet)

# ===========================================================================
# 任务1：11个特征在外部队列 ST002787 中的检出情况与重叠表
# ===========================================================================
cat("\n===== 任务1：特征重叠表 =====\n")

# 读取外部队列数据
pos_2787 <- read.table("st002787_positive_clean.txt", header = TRUE,
                       sep = "\t", row.names = 1, check.names = FALSE)
neg_2787 <- read.table("st002787_negative_clean.txt", header = TRUE,
                       sep = "\t", row.names = 1, check.names = FALSE)

all_2787_features <- unique(c(rownames(pos_2787), rownames(neg_2787)))
matched <- selected_features %in% all_2787_features

overlap_table <- data.frame(
  Feature = selected_features,
  In_ST003798 = rep("Yes", 11),
  In_ST002787 = ifelse(matched, "Yes", "No"),
  stringsAsFactors = FALSE
)

print(overlap_table)
write.csv(overlap_table, "Table_S6_Feature_Overlap.csv",
          row.names = FALSE)

cat(sprintf("特征检出: %d / 11\n", sum(matched)))

# ===========================================================================
# 任务2：OOB 误差收敛图
# ===========================================================================
cat("\n===== 任务2：OOB误差收敛 =====\n")

pdf("Figure_S3.pdf", width = 7, height = 5)
plot(rf_model, main = "Out-of-Bag Error Convergence")
legend("topright", legend = colnames(rf_model$err.rate),
       lty = 1, col = 1:ncol(rf_model$err.rate), cex = 0.8)
dev.off()

oob_errors <- rf_model$err.rate[, "OOB"]
min_oob <- min(oob_errors)
min_tree <- which.min(oob_errors)
cat(sprintf("最小OOB误差: %.4f (在第 %d 棵树时达到)\n", min_oob, min_tree))

# ===========================================================================
# 任务3：缺失值插补敏感性分析（LOD/2 替代方案）
# ===========================================================================
cat("\n===== 任务3：缺失值插补敏感性分析 =====\n")

# 从 raw_df 重建原始未插补矩阵
refmet_col <- which(colnames(raw_df) == "RefMet_name")
if(length(refmet_col) == 0) {
  refmet_col <- which(colnames(raw_df) == "Metabolite_name")
}
sample_cols <- (refmet_col + 1):ncol(raw_df)

raw_matrix <- as.matrix(raw_df[, sample_cols])
raw_matrix <- t(raw_matrix)
raw_matrix <- apply(raw_matrix, 2, as.numeric)
raw_matrix[raw_matrix == 0] <- NA
colnames(raw_matrix) <- raw_df[, refmet_col]
rownames(raw_matrix) <- colnames(raw_df)[sample_cols]

cat(sprintf("原始矩阵: %d 样本 x %d 代谢物\n", nrow(raw_matrix), ncol(raw_matrix)))
cat(sprintf("缺失率: %.2f%%\n",
            100 * sum(is.na(raw_matrix)) / (nrow(raw_matrix) * ncol(raw_matrix))))

# LOD/2 插补函数
lod_impute <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  if (is.finite(min_val) && min_val > 0) {
    x[is.na(x)] <- min_val / 2
  } else {
    x[is.na(x)] <- 1e-6
  }
  return(x)
}

raw_lod <- apply(raw_matrix, 2, lod_impute)
raw_lod_log2 <- log2(raw_lod + 1)

# 匹配样本
sample_ids <- data_clean_3798$SampleID
group_labels <- data_clean_3798$Group
matched_samples <- intersect(rownames(raw_lod_log2), sample_ids)

lod_matrix <- raw_lod_log2[matched_samples, ]
matched_groups <- group_labels[match(matched_samples, sample_ids)]

# 极端表型子集
extreme_idx <- which(matched_groups %in% c("Healthy Control", "Colorectal Cancer"))
lod_extreme <- lod_matrix[extreme_idx, ]
lod_groups <- droplevels(as.factor(matched_groups[extreme_idx]))

cat("极端表型分组:\n")
print(table(lod_groups))

# LASSO
set.seed(42)
lod_cv <- cv.glmnet(as.matrix(lod_extreme), lod_groups,
                    family = "binomial", nfolds = 10, alpha = 1)
lod_coefs <- coef(lod_cv, s = "lambda.1se")
lod_selected <- rownames(lod_coefs)[which(lod_coefs[,1] != 0)][-1]

cat(sprintf("LOD/2 插补后 LASSO 选出特征数: %d\n", length(lod_selected)))

# 标准化特征名并比较重叠
standardize_name <- function(x) {
  x <- gsub("[()]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- gsub(":", " ", x)
  x <- gsub(";", " ", x)
  return(x)
}

orig_std <- standardize_name(selected_features)
lod_std <- standardize_name(lod_selected)

common_std <- intersect(orig_std, lod_std)
jaccard_std <- length(common_std) / length(union(orig_std, lod_std))

cat(sprintf("标准化后重叠特征数: %d / 11\n", length(common_std)))
cat("重叠特征名:", common_std, "\n")
cat(sprintf("Jaccard 相似度: %.3f\n", jaccard_std))

sens_result <- data.frame(
  Imputation_Method = c("Original 30%_KNN_k10", "LOD_half_min"),
  Feature_Count = c(length(selected_features), length(lod_selected)),
  Overlap_with_Original = c(length(selected_features), length(common_std)),
  stringsAsFactors = FALSE
)
print(sens_result)
write.csv(sens_result, "Table_S7_Sensitivity_Analysis.csv", row.names = FALSE)

# ===========================================================================
# 任务4：伪时间根节点翻转验证
# ===========================================================================
cat("\n===== 任务4：伪时间根节点翻转验证 =====\n")

# 补充分组信息
pt_row_idx <- as.numeric(rownames(curve_df))
curve_df$SampleID <- data_clean_3798$SampleID[pt_row_idx]
curve_df$Group <- data_clean_3798$Group[pt_row_idx]

# 提取腺瘤并添加高风险标签
adenoma_pt <- curve_df[curve_df$Group == "Adenoma", ]
adenoma_row_idx <- which(data_clean_3798$Group == "Adenoma")

if (identical(as.character(adenoma_row_idx), rownames(df_adenoma_full))) {
  adenoma_pt$Subtype <- df_adenoma_full$Subtype
} else {
  stop("腺瘤样本顺序不一致，请手动排查！")
}
adenoma_pt$HighRisk <- adenoma_pt$Subtype == "High-Risk Adenoma"

cat("高风险分布:\n")
print(table(adenoma_pt$HighRisk))

# 构造反向伪时间
crc_max <- max(curve_df$lambda[curve_df$Group == "Colorectal Cancer"])
adenoma_pt$lambda_rev <- crc_max - adenoma_pt$lambda

# 统计检验
wilcox_orig <- wilcox.test(lambda ~ HighRisk, data = adenoma_pt)
wilcox_rev <- wilcox.test(lambda_rev ~ HighRisk, data = adenoma_pt)

# 绘制对比箱线图
pdf("Figure_S4.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
boxplot(lambda ~ HighRisk, data = adenoma_pt,
        main = "Original Pseudotime (Root = Healthy)",
        names = c("Low-Risk", "High-Risk"),
        ylab = "Pseudotime (lambda)", col = c("#2196F3", "#FF9800"))
boxplot(lambda_rev ~ HighRisk, data = adenoma_pt,
        main = "Reversed Pseudotime (Root = CRC)",
        names = c("Low-Risk", "High-Risk"),
        ylab = "Reversed Pseudotime", col = c("#2196F3", "#FF9800"))
par(mfrow = c(1, 1))
dev.off()

cat(sprintf("正向伪时间(根=健康): 高风险均值=%.2f, 低风险均值=%.2f, p=%.4f\n",
            mean(adenoma_pt$lambda[adenoma_pt$HighRisk]),
            mean(adenoma_pt$lambda[!adenoma_pt$HighRisk]),
            wilcox_orig$p.value))
cat(sprintf("反向伪时间(根=CRC): 高风险均值=%.2f, 低风险均值=%.2f, p=%.4f\n",
            mean(adenoma_pt$lambda_rev[adenoma_pt$HighRisk]),
            mean(adenoma_pt$lambda_rev[!adenoma_pt$HighRisk]),
            wilcox_rev$p.value))

cat("\n===== 所有补充分析完成 =====\n")
cat("生成文件清单:\n")
cat("1. Table_SX_Feature_Overlap_Between_Cohorts.csv\n")
cat("2. Figure_SX_OOB_Convergence.pdf\n")
cat("3. Table_SX_Sensitivity_Analysis.csv\n")
cat("4. Figure_SX_Pseudotime_Root_Sensitivity.pdf\n")
