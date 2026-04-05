# ==============================================================================
# Script 06: Final Essential Statistics & TCGA Covariate Adjustment
# Description: Calculates Hartigan's dip test, External Cohort AUCs, and 
#              performs Multivariate Regression on TCGA-COAD data.
# ==============================================================================

# 1. 强制加载所有必需的基础包
suppressMessages({
  library(diptest)  # 任务 B.1: Hartigan's dip test
  library(pROC)     # 任务 B.2: 外部独立 AUC
  library(dplyr)    # 任务 C: TCGA 数据处理
  library(tidyr)
  library(tibble)
  library(stats)
})

cat("\n====================================================================\n")
cat("🚀 启动 [任务 B & C]：最终关键统计量补充与临床混杂因素校正...\n")

# ------------------------------------------------------------------------------
# 任务 B.1: 计算 3.2 节的 Hartigan's dip test (腺瘤异质性)
# ------------------------------------------------------------------------------
cat("\n>>> [任务 B.1] 计算腺瘤预测概率分布的双模态检验 (Hartigan's dip test)...\n")

if (exists("adenoma_prob")) {
  set.seed(2026)
  dip_res <- dip.test(adenoma_prob)
  
  dip_d <- round(dip_res$statistic, 4)
  dip_p <- signif(dip_res$p.value, 3)
  
  cat(sprintf("   - D statistic: %.4f\n", dip_d))
  cat(sprintf("   - Exact P-value: %.3e\n", dip_p))
  cat("💡 提示：将这两个数字填入 3.2 节的 (Fig. 3) 描述中。\n")
} else {
  cat("❌ 错误：未找到 adenoma_prob，请确保主脚本 Step 4 已运行。\n")
}

# ------------------------------------------------------------------------------
# 任务 B.2: 计算 3.6 节外部队列 (ST002787) 中 CE(20:4) 的独立诊断效能
# ------------------------------------------------------------------------------
cat("\n>>> [任务 B.2] 评估外部验证队列中单一标志物 CE(20:4) 的独立 AUC...\n")

if (exists("merged_data")) {
  ce_col <- grep("CE\\(20:4\\)", colnames(merged_data), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(ce_col)) {
    ce204_data <- merged_data %>% 
      select(SampleID, Group, Abundance = all_of(ce_col)) %>% 
      mutate(Abundance = as.numeric(Abundance)) %>%
      filter(!is.na(Abundance) & Abundance > 0)
    
    # a. CRC vs Health (修复了换行错误)
    data_crc_hc <- ce204_data %>% filter(Group %in% c("Healthy Control", "Colorectal Cancer"))
    roc_crc <- roc(data_crc_hc$Group, data_crc_hc$Abundance, 
                   levels = c("Healthy Control", "Colorectal Cancer"), direction = "<", quiet = TRUE)
    auc_crc <- round(auc(roc_crc), 3)
    
    # b. Adenoma vs Health (修复了换行错误)
    data_ade_hc <- ce204_data %>% filter(Group %in% c("Healthy Control", "Adenoma"))
    roc_ade <- roc(data_ade_hc$Group, data_ade_hc$Abundance, 
                   levels = c("Healthy Control", "Adenoma"), direction = "<", quiet = TRUE)
    auc_ade <- round(auc(roc_ade), 3)
    
    cat(sprintf("   - CE(20:4) 独立区分 CRC vs Health AUC: %.3f\n", auc_crc))
    cat(sprintf("   - CE(20:4) 独立区分 Adenoma vs Health AUC: %.3f\n", auc_ade))
    cat("💡 提示：您可以将这些 AUC 值补充到 3.6 节中，证明其作为独立 Biomarker 的潜力。\n")
  } else {
    cat("❌ 未在 merged_data 中找到 CE(20:4) 列，请检查拼写。\n")
  }
} else {
  cat("❌ 错误：未找到 merged_data，请确保主脚本 Step 9 已运行。\n")
}

# ------------------------------------------------------------------------------
# 任务 C: TCGA 差异表达的混杂因素校正 (智能匹配 Age, Sex 列)
# ------------------------------------------------------------------------------
cat("\n>>> [任务 C] 对 TCGA 核心基因进行多变量线性回归 (仅校正 Age, Sex)...\n")

if (exists("tcga_data") && exists("target_tpm")) {
  
  # 1. 尝试用最底层的 S4 提取方法获取临床数据
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    clin_raw <- as.data.frame(SummarizedExperiment::colData(tcga_data))
  } else {
    clin_raw <- as.data.frame(tcga_data@colData)
  }
  
  # 2. 智能寻找 Age 和 Sex 列
  age_col <- grep("age_at_index|age_at_diagnosis|age", colnames(clin_raw), value=TRUE, ignore.case=TRUE)[1]
  sex_col <- grep("gender|sex", colnames(clin_raw), value=TRUE, ignore.case=TRUE)[1]
  
  if (is.na(age_col) | is.na(sex_col)) {
    cat("❌ 错误：无法在 TCGA 临床数据中找到 Age 或 Sex 列！\n")
  } else {
    
    clin_data <- data.frame(
      SampleID = clin_raw$barcode,
      Age = as.numeric(clin_raw[[age_col]]),
      Sex = as.character(clin_raw[[sex_col]])
    )
    
    # 3. 从 target_tpm 重建长矩阵并合并
    df_expr_clean <- as.data.frame(t(target_tpm)) %>% 
      rownames_to_column("SampleID") %>%
      pivot_longer(cols = -SampleID, names_to = "Gene", values_to = "TPM") %>%
      mutate(Log2TPM = log2(TPM + 1))
    
    tcga_glm_df <- df_expr_clean %>% 
      left_join(clin_data, by = "SampleID") %>%
      left_join(data.frame(SampleID = colnames(target_tpm), 
                           Type = ifelse(tcga_data$sample_type == "Solid Tissue Normal", "Normal", "Tumor")), 
                by = "SampleID") %>%
      filter(!is.na(Age)) %>%
      mutate(Type = factor(Type, levels = c("Normal", "Tumor")), 
             Sex = as.factor(Sex))
    
    # 4. 执行 GLM
    target_genes_glm <- unique(tcga_glm_df$Gene)
    glm_results <- data.frame(
      Gene = character(),
      Tumor_Effect_Estimate = numeric(),
      Tumor_Effect_Pvalue = numeric(),
      Significance_After_Adjustment = character(),
      stringsAsFactors = FALSE
    )
    
    for (g in target_genes_glm) {
      gene_df <- tcga_glm_df %>% filter(Gene == g)
      fit_glm <- glm(Log2TPM ~ Type + Age + Sex, data = gene_df, family = gaussian())
      
      coef_summary <- summary(fit_glm)$coefficients
      row_idx <- grep("TypeTumor", rownames(coef_summary))
      
      if (length(row_idx) > 0) {
        est <- coef_summary[row_idx, "Estimate"]
        pval <- coef_summary[row_idx, "Pr(>|t|)"]
        sig <- ifelse(pval < 0.05, "Yes (P<0.05)", "No (P>=0.05)")
        
        glm_results <- rbind(glm_results, data.frame(
          Gene = g,
          Tumor_Effect_Estimate = round(est, 4),
          Tumor_Effect_Pvalue = signif(pval, 3),
          Significance_After_Adjustment = sig
        ))
      }
    }
    
    cat("\n🏆 TCGA 多变量回归校正结果 (Log2TPM ~ Tumor_Status + Age + Sex)：\n")
    print(glm_results)
    write.csv(glm_results, "Table_S6_TCGA_Multivariate_Adjusted.csv", row.names = FALSE)
    cat("\n✅ [终极任务全部完成] 所有统计学漏洞均已被完美填补！\n")
    cat("====================================================================\n")
  }
} else {
  cat("❌ 错误：环境中未找到 tcga_data 或 target_tpm，请确保主脚本 Step 8 已运行。\n")
}
