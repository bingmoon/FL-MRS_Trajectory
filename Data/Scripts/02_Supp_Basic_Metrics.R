# ==============================================================================
# Supplementary Script: Detailed Diagnostic Metrics for the Independent Test Set
# (In response to Reviewer comments requesting Sensitivity, Specificity, PPV, NPV)
# 
# INSTRUCTIONS: 
# Please run this script AFTER executing the main pipeline (Step 1 to 4), 
# as it relies on the 'test_prob', 'test_data', and 'best_cutoff' objects 
# generated in the Global Environment.
# ==============================================================================

# 1. Load required library
suppressMessages(library(caret))

# 2. Safety Check: Ensure main pipeline has been executed
if (!exists("test_prob") | !exists("test_data") | !exists("best_cutoff")) {
  stop("ERROR: Required objects not found. Please run the main RF pipeline (Step 1-4) first.")
}

cat("\n======================================================\n")
cat("📊 Calculating Comprehensive Diagnostic Metrics (Test Set)...\n")

# 3. Convert probabilities to binary class predictions using the optimal Cutoff
# Note: Reordering levels so that 'Colorectal Cancer' is treated as the Positive class
pred_labels <- ifelse(test_prob >= best_cutoff, "Colorectal Cancer", "Healthy Control")
pred_labels <- factor(pred_labels, levels = c("Colorectal Cancer", "Healthy Control")) 

actual_labels <- factor(test_data$Group, levels = c("Colorectal Cancer", "Healthy Control"))

# 4. Generate Confusion Matrix
conf_matrix <- confusionMatrix(pred_labels, actual_labels, positive = "Colorectal Cancer")

# 5. Extract specific metrics
metrics_summary <- data.frame(
  Metric = c("Optimal Cutoff", "Accuracy", "Sensitivity (TPR)", "Specificity (TNR)", 
             "Positive Predictive Value (PPV)", "Negative Predictive Value (NPV)"),
  Value = c(
    round(best_cutoff, 3),
    round(conf_matrix$overall["Accuracy"] * 100, 1),
    round(conf_matrix$byClass["Sensitivity"] * 100, 1),
    round(conf_matrix$byClass["Specificity"] * 100, 1),
    round(conf_matrix$byClass["Pos Pred Value"] * 100, 1),
    round(conf_matrix$byClass["Neg Pred Value"] * 100, 1)
  ),
  Unit = c("Score", "%", "%", "%", "%", "%")
)

# 6. Print to Console
cat(sprintf(" - Optimal Cutoff: %.3f\n", metrics_summary$Value[1]))
cat(sprintf(" - Accuracy: %.1f%%\n", metrics_summary$Value[2]))
cat(sprintf(" - Sensitivity: %.1f%%\n", metrics_summary$Value[3]))
cat(sprintf(" - Specificity: %.1f%%\n", metrics_summary$Value[4]))
cat(sprintf(" - PPV: %.1f%%\n", metrics_summary$Value[5]))
cat(sprintf(" - NPV: %.1f%%\n", metrics_summary$Value[6]))
cat("======================================================\n")

# 7. Export to CSV for manuscript Supplementary Table
output_file <- "Table_S3_Detailed_Diagnostic_Metrics.csv"
write.csv(metrics_summary, file = output_file, row.names = FALSE)
cat(sprintf("✅ Success! Metrics exported to: %s\n", output_file))
