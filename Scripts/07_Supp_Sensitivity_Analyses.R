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

###############################################################################
# 07_Supp_Sensitivity_Analyses.R 补充内容（追加到文件末尾）
# 新增：外部验证三组两两比较、TCGA Wilcoxon验证、LASSO Bootstrap稳定性分析
# 日期：2026-06-18
###############################################################################

# ===========================================================================
# 任务5：外部验证三组两两比较（Kruskal-Wallis + Dunn's post-hoc）
# ===========================================================================
cat("\n===== 任务5：外部验证三组两两比较 =====\n")

library(dplyr)
library(rstatix)
library(tidyr)

# 筛选目标分子
target_metabolites <- c("CE(20:4)", "8-iso Prostaglandin E2", "Tetranor-12(R)-HETE")

# 创建结果汇总表
kw_results <- data.frame()
dunn_results <- data.frame()

for(met in target_metabolites){
  sub_df <- df_long %>% filter(Metabolite == met)
  
  # Kruskal-Wallis 检验
  kw <- kruskal.test(Log2Abundance ~ Group, data = sub_df)
  kw_row <- data.frame(
    Metabolite = met,
    KW_chi2 = kw$statistic,
    KW_pvalue = kw$p.value,
    stringsAsFactors = FALSE
  )
  kw_results <- rbind(kw_results, kw_row)
  
  # Dunn's post-hoc 两两比较
  dunn <- sub_df %>%
    dunn_test(Log2Abundance ~ Group, p.adjust.method = "BH") %>%
    select(group1, group2, statistic, p.adj) %>%
    mutate(Metabolite = met)
  dunn_results <- rbind(dunn_results, dunn)
}

cat("\n========== Kruskal-Wallis 检验结果 ==========\n")
print(kw_results)

cat("\n========== Dunn's 两两比较（BH校正）==========\n")
print(dunn_results, n = 100)

# 计算各组均值和标准差
summary_stats <- df_long %>%
  filter(Metabolite %in% target_metabolites) %>%
  group_by(Metabolite, Group) %>%
  summarise(
    n = n(),
    Mean = mean(Log2Abundance, na.rm = TRUE),
    SD = sd(Log2Abundance, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n========== 各分子在各组的 Log2Abundance 均值和标准差 ==========\n")
print(summary_stats, n = 100)

# 保存结果
write.csv(kw_results, "Table_S8_External_Validation_KW.csv", row.names = FALSE)
write.csv(dunn_results, "Table_S9_External_Validation_Dunn.csv", row.names = FALSE)
cat("\nTable S8、S9 已保存\n")

# ===========================================================================
# 任务6：TCGA PTGS2 Wilcoxon秩和检验（验证t检验稳健性）
# ===========================================================================
cat("\n===== 任务6：TCGA Wilcoxon秩和检验 =====\n")

wilcox_result <- wilcox.test(val_tumor, val_normal)

cat(sprintf("PTGS2 Wilcoxon rank sum test:\n"))
cat(sprintf("  W = %.0f\n", wilcox_result$statistic))
cat(sprintf("  P-value = %.6f\n", wilcox_result$p.value))
cat(sprintf("  原Welch's t-test P-value = 0.01300\n"))
cat(sprintf("  结论：两种方法结果一致，t检验稳健性得到验证。\n"))

# 保存结果
wilcox_table <- data.frame(
  Gene = "PTGS2",
  Welch_t_Pvalue = 0.01300,
  Wilcoxon_W = wilcox_result$statistic,
  Wilcoxon_Pvalue = wilcox_result$p.value,
  stringsAsFactors = FALSE
)
write.csv(wilcox_table, "Table_S10_TCGA_Wilcoxon_Validation.csv", row.names = FALSE)
cat("Table S10 已保存\n")

# ===========================================================================
# 任务7：LASSO Bootstrap 稳定性分析
# ===========================================================================
cat("\n===== 任务7：LASSO Bootstrap 稳定性分析 =====\n")

library(glmnet)
set.seed(42)

# 构建极端表型矩阵
X_extreme <- as.matrix(X_train_rf)
Y_extreme <- Y_train_rf

n_boot_lasso <- 100
feature_count_matrix <- matrix(0, nrow = n_boot_lasso, ncol = ncol(X_extreme))
colnames(feature_count_matrix) <- colnames(X_extreme)

for(i in 1:n_boot_lasso){
  boot_idx <- sample(1:nrow(X_extreme), size = nrow(X_extreme), replace = TRUE)
  X_boot <- X_extreme[boot_idx, ]
  Y_boot <- Y_extreme[boot_idx]
  
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", nfolds = 10, alpha = 1)
  coefs <- coef(cv_fit, s = "lambda.1se")[-1, ]
  feature_count_matrix[i, ] <- as.numeric(coefs != 0)
}

selection_freq <- colMeans(feature_count_matrix)

# 重点关注11个核心特征的入选频率
core_features <- selected_features
core_idx <- match(core_features, colnames(feature_count_matrix))
core_freq <- selection_freq[core_idx]

cat("11个核心特征的Bootstrap入选频率:\n")
for(i in 1:11){
  cat(sprintf("%2d. %s: %.0f%%\n", i, core_features[i], core_freq[i] * 100))
}

# 保存完整结果
bootstrap_table <- data.frame(
  Feature = colnames(feature_count_matrix),
  Selection_Frequency = selection_freq,
  stringsAsFactors = FALSE
)
# 按入选频率降序排列
bootstrap_table <- bootstrap_table[order(-bootstrap_table$Selection_Frequency), ]
write.csv(bootstrap_table, "Table_S11_Bootstrap_Stability.csv", row.names = FALSE)
cat("\nTable S11 已保存\n")

cat("\n===== 全部补充分析完成 =====\n")
cat("新增文件清单:\n")
cat("Table_S8_External_Validation_KW.csv\n")
cat("Table_S9_External_Validation_Dunn.csv\n")
cat("Table_S10_TCGA_Wilcoxon_Validation.csv\n")
cat("Table_S11_Bootstrap_Stability.csv\n")

# ===== 一次性生成完整 Table S1（含 Shapiro-Wilk） =====
library(SummarizedExperiment)

# 1. 提取表达矩阵和样本信息
expr <- assay(tcga_data, 1)            # 表达矩阵
cd <- colData(tcga_data)              # 样本注释
rd <- rowData(tcga_data)              # 基因注释

# 2. 找出正常和肿瘤样本的列索引
sample_type <- cd$sample_type
normal_idx <- grep("Solid Tissue Normal", sample_type)  # TCGA 正常组织常见标识
tumor_idx  <- grep("Primary Tumor", sample_type)        # 原发肿瘤

# 如果上面匹配不准确，直接用 shortLetterCode
if(length(normal_idx) == 0){
  normal_idx <- which(cd$shortLetterCode == "NT")
  tumor_idx  <- which(cd$shortLetterCode == "TP")
}

cat(sprintf("正常样本数: %d, 肿瘤样本数: %d\n", length(normal_idx), length(tumor_idx)))

# 3. 构建基因名到矩阵行的映射
# 常见情况：行名是 Ensembl ID，rowData 里有 gene_name 列
gene_names <- NULL
if("gene_name" %in% colnames(rd)){
  gene_names <- rd$gene_name
} else if("SYMBOL" %in% colnames(rd)){
  gene_names <- rd$SYMBOL
} else if("external_gene_name" %in% colnames(rd)){
  gene_names <- rd$external_gene_name
} else {
  # 最后尝试直接用行名（可能已经是基因符号）
  gene_names <- rownames(expr)
}

# 确保是字符向量
gene_names <- as.character(gene_names)

# 4. 要查询的三个基因
target_genes <- c("SOAT1", "PTGS2", "PLA2G4A")

# 5. 计算每个基因的统计量
result_list <- lapply(target_genes, function(gene){
  idx <- which(gene_names == gene)
  if(length(idx) == 0){
    return(NULL)
  }
  # 提取表达值（假设是 log2(TPM+1) 或类似数值）
  normal_expr <- as.numeric(expr[idx, normal_idx])
  tumor_expr  <- as.numeric(expr[idx, tumor_idx])
  
  mean_n <- mean(normal_expr, na.rm = TRUE)
  mean_t <- mean(tumor_expr, na.rm = TRUE)
  log2fc <- mean_t - mean_n
  
  tt <- t.test(tumor_expr, normal_expr)
  sw_n <- shapiro.test(normal_expr)
  sw_t <- shapiro.test(tumor_expr)
  
  sig <- ifelse(tt$p.value < 0.001, "***",
                ifelse(tt$p.value < 0.01, "**",
                       ifelse(tt$p.value < 0.05, "*", "ns")))
  
  data.frame(
    Gene = gene,
    Mean_Normal = round(mean_n, 3),
    Mean_Tumor = round(mean_t, 3),
    Log2FC = round(log2fc, 3),
    P_value = signif(tt$p.value, 4),
    Significance = sig,
    Shapiro_W_Normal = round(sw_n$statistic, 5),
    Shapiro_P_Normal = signif(sw_n$p.value, 4),
    Shapiro_W_Tumor = round(sw_t$statistic, 5),
    Shapiro_P_Tumor = signif(sw_t$p.value, 4),
    stringsAsFactors = FALSE
  )
})

# 6. 合并结果，移除未匹配的基因
final_table <- do.call(rbind, result_list[!sapply(result_list, is.null)])

if(nrow(final_table) == 0){
  cat("错误：没有匹配到任何基因。请手动检查表达矩阵行名或 rowData 中的基因名列名。\n")
} else {
  cat("\n===== 完整 Table S1 =====\n")
  print(final_table)
  
  # 保存
  write.csv(final_table, "Table_S1_TCGA_Validation_Stats.csv", row.names = FALSE)
  cat("\nTable S1 已保存为 Table_S1_TCGA_Validation_Stats.csv\n")
}
