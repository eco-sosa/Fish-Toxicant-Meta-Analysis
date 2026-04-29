# =============================================================================
# Publication Bias Analysis for All Meta-Analyses
# Larval Success: Borenstein-corrected (Inclusive Approach)
# Maternal Characteristics: Correlation-based (Fisher's Z)
# Includes: Egger's test, Funnel plots, Trim-and-fill correction
# =============================================================================

library(metafor)
library(ggplot2)
library(readxl)

# Set working directory to the repo root before sourcing this script.
# In RStudio: Session > Set Working Directory > To Project Directory.
# Otherwise run: setwd("path/to/Hg_PCBs_Meta_Repo")

# =============================================================================
# Load Data
# =============================================================================

# Borenstein-corrected CSVs for larval success (Hedges' d)
hatching_data <- read.csv("data/borenstein_corrected_hatching.csv", stringsAsFactors = FALSE)
mortality_data <- read.csv("data/borenstein_corrected_mortality.csv", stringsAsFactors = FALSE)
growth_data <- read.csv("data/borenstein_corrected_growth.csv", stringsAsFactors = FALSE)

# Load Size and Muscle from ICC file (Fisher's Z - no Borenstein correction needed)
icc_data <- read_excel("data/ICC List for Meta.xlsx", sheet = "Sheet1")
names(icc_data) <- c("meta_analysis", "study_id", "toxicant", "yi", "vi",
                     "treated_conc", "control_n", "treated_n")

# Extract Size and Muscle data
size_data <- icc_data[icc_data$meta_analysis == "Size", ]
size_data$yi_corrected <- size_data$yi
size_data$vi_corrected <- size_data$vi

muscle_data <- icc_data[icc_data$meta_analysis == "Muscle", ]
muscle_data$yi_corrected <- muscle_data$yi
muscle_data$vi_corrected <- muscle_data$vi

# =============================================================================
# FUNCTION: Run complete publication bias analysis
# =============================================================================

run_pub_bias_analysis <- function(data, analysis_name, toxicant = NULL,
                                   effect_type = "Hedges' d") {

  # Subset by toxicant if specified
  if (!is.null(toxicant)) {
    data <- data[data$toxicant == toxicant, ]
    label <- paste(analysis_name, "-", toxicant)
  } else {
    label <- paste(analysis_name, "- Overall")
  }

  # Extract effect sizes and variances
  yi <- data$yi_corrected
  vi <- data$vi_corrected
  sei <- sqrt(vi)
  k <- length(yi)

  cat("\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat("PUBLICATION BIAS ANALYSIS:", label, "\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat("Number of studies (k):", k, "\n\n")

  # Check minimum sample size for meaningful analysis
  if (k < 3) {
    cat("WARNING: Too few studies (k < 3) for publication bias analysis.\n")
    cat("Skipping this analysis.\n")
    return(NULL)
  }

  # ---------------------------------------------
  # 1. Fit fixed-effect model (to match Borenstein correction approach)
  # ---------------------------------------------
  cat("--- Fixed-Effect Model ---\n")

  fe_model <- rma(yi = yi, vi = vi, method = "FE")

  cat("Effect size (", effect_type, "):", round(fe_model$b[1], 4), "\n", sep = "")
  cat("95% CI: [", round(fe_model$ci.lb, 4), ", ", round(fe_model$ci.ub, 4), "]\n", sep = "")
  cat("SE:", round(fe_model$se, 4), "\n")
  cat("z =", round(fe_model$zval, 3), ", p =", format.pval(fe_model$pval, digits = 4), "\n")
  cat("Q =", round(fe_model$QE, 2), ", df =", fe_model$k - 1,
      ", p =", format.pval(fe_model$QEp, digits = 4), "\n")
  cat("I^2 =", round(fe_model$I2, 1), "%\n\n")

  # Use fe_model for downstream analyses
  re_model <- fe_model

  # ---------------------------------------------
  # 2. Egger's Regression Test
  # ---------------------------------------------
  cat("--- Egger's Regression Test for Funnel Plot Asymmetry ---\n")

  if (k >= 10) {
    egger_test <- regtest(re_model, model = "lm")
    cat("Egger's test statistic (intercept):", round(egger_test$zval, 3), "\n")
    cat("p-value:", format.pval(egger_test$pval, digits = 4), "\n")

    if (egger_test$pval < 0.05) {
      cat("INTERPRETATION: Significant asymmetry detected (p < 0.05).\n")
      cat("  This suggests potential publication bias.\n\n")
    } else if (egger_test$pval < 0.10) {
      cat("INTERPRETATION: Marginal asymmetry (0.05 < p < 0.10).\n")
      cat("  Results should be interpreted with caution.\n\n")
    } else {
      cat("INTERPRETATION: No significant asymmetry detected (p >= 0.10).\n")
      cat("  No strong evidence of publication bias.\n\n")
    }
  } else {
    cat("WARNING: Egger's test has low power with k < 10 studies.\n")
    cat("Results should be interpreted with caution.\n")

    egger_test <- regtest(re_model, model = "lm")
    cat("Egger's test statistic (intercept):", round(egger_test$zval, 3), "\n")
    cat("p-value:", format.pval(egger_test$pval, digits = 4), "\n\n")
  }

  # ---------------------------------------------
  # 3. Trim-and-Fill Analysis
  # ---------------------------------------------
  cat("--- Trim-and-Fill Analysis ---\n")

  # Run trim-and-fill (default: right side, estimator L0)
  tf_model <- trimfill(re_model, side = "right")

  # Also try left side
  tf_model_left <- trimfill(re_model, side = "left")

  # Determine which side has imputed studies
  k0_right <- tf_model$k0
  k0_left <- tf_model_left$k0

  cat("Studies imputed (right side):", k0_right, "\n")
  cat("Studies imputed (left side):", k0_left, "\n\n")

  # Use the side with more imputed studies for adjustment
  if (k0_right >= k0_left) {
    tf_final <- tf_model
    side_used <- "right"
  } else {
    tf_final <- tf_model_left
    side_used <- "left"
  }

  k0_final <- tf_final$k0

  if (k0_final > 0) {
    cat("Trim-and-fill identified", k0_final, "missing studies on the", side_used, "side.\n\n")

    cat("ADJUSTED EFFECT SIZE (after trim-and-fill):\n")
    cat("Effect size (", effect_type, "):", round(tf_final$b[1], 4), "\n", sep = "")
    cat("95% CI: [", round(tf_final$ci.lb, 4), ", ", round(tf_final$ci.ub, 4), "]\n", sep = "")
    cat("SE:", round(tf_final$se, 4), "\n")
    cat("I^2 (adjusted):", round(tf_final$I2, 1), "%\n\n")

    # Calculate change
    original_effect <- re_model$b[1]
    adjusted_effect <- tf_final$b[1]
    change <- adjusted_effect - original_effect
    pct_change <- (change / abs(original_effect)) * 100

    cat("COMPARISON:\n")
    cat("  Original effect:", round(original_effect, 4), "\n")
    cat("  Adjusted effect:", round(adjusted_effect, 4), "\n")
    cat("  Change:", round(change, 4), "(", round(pct_change, 1), "%)\n\n")

    # Interpretation
    if (abs(pct_change) < 10) {
      cat("INTERPRETATION: Minimal change (<10%). Original results appear robust.\n")
    } else if (abs(pct_change) < 25) {
      cat("INTERPRETATION: Moderate change (10-25%). Some publication bias may exist.\n")
    } else {
      cat("INTERPRETATION: Substantial change (>25%). Publication bias likely affects results.\n")
    }

    # Check if significance changes
    orig_sig <- re_model$pval < 0.05
    adj_sig <- tf_final$pval < 0.05

    if (orig_sig && !adj_sig) {
      cat("  NOTE: Effect becomes non-significant after adjustment.\n")
    } else if (!orig_sig && adj_sig) {
      cat("  NOTE: Effect becomes significant after adjustment.\n")
    } else {
      cat("  Statistical significance unchanged after adjustment.\n")
    }

  } else {
    cat("No missing studies identified. No adjustment needed.\n")
    cat("Original results appear unaffected by publication bias.\n")
    tf_final <- re_model
  }

  cat("\n")

  # ---------------------------------------------
  # Return results
  # ---------------------------------------------
  results <- list(
    label = label,
    k = k,
    effect_type = effect_type,
    original_model = re_model,
    egger_test = egger_test,
    trimfill_model = tf_final,
    k0_imputed = k0_final,
    side_used = side_used,
    original_effect = as.numeric(re_model$b[1]),
    original_ci_lb = re_model$ci.lb,
    original_ci_ub = re_model$ci.ub,
    original_pval = re_model$pval,
    original_I2 = re_model$I2,
    adjusted_effect = as.numeric(tf_final$b[1]),
    adjusted_ci_lb = tf_final$ci.lb,
    adjusted_ci_ub = tf_final$ci.ub,
    adjusted_pval = tf_final$pval,
    adjusted_I2 = tf_final$I2,
    egger_pval = egger_test$pval
  )

  return(results)
}

# =============================================================================
# FUNCTION: Create funnel plot with trim-and-fill (kept for individual plots)
# =============================================================================

create_funnel_plot <- function(results, filename = NULL) {

  if (is.null(results)) return(NULL)

  # Determine x-axis label based on effect type
  xlab_text <- paste0("Effect Size (", results$effect_type, ")")

  # Create plot
  if (!is.null(filename)) {
    png(filename, width = 6, height = 5, units = "in", res = 300)
  }

  # Funnel plot with trim-and-fill using precision (1/SE)
  funnel(results$trimfill_model,
         yaxis = "seinv",
         main = results$label,
         xlab = xlab_text,
         ylab = "Precision (1/SE)")

  if (!is.null(filename)) {
    dev.off()
    cat("Funnel plot saved to:", filename, "\n")
  }
}

# =============================================================================
# RUN ANALYSES
# =============================================================================

cat("\n")
cat(strrep("*", 70), "\n", sep = "")
cat("PUBLICATION BIAS ANALYSIS - ALL META-ANALYSES\n")
cat("Larval Success: Borenstein-corrected (Inclusive Approach)\n")
cat("Maternal Characteristics: Correlation-based (Fisher's Z)\n")
cat(strrep("*", 70), "\n")

# Store all results
all_results <- list()
summary_table <- data.frame()

# ---------------------------------------------
# HATCHING SUCCESS (by toxicant only)
# ---------------------------------------------

hatching_hg <- run_pub_bias_analysis(hatching_data, "Hatching Success", toxicant = "Hg")
if (!is.null(hatching_hg)) {
  all_results[["Hatching_Hg"]] <- hatching_hg
  create_funnel_plot(hatching_hg, "figures/funnel_hatching_hg.png")
}

hatching_pcb <- run_pub_bias_analysis(hatching_data, "Hatching Success", toxicant = "PCBs")
if (!is.null(hatching_pcb)) {
  all_results[["Hatching_PCBs"]] <- hatching_pcb
  create_funnel_plot(hatching_pcb, "figures/funnel_hatching_pcbs.png")
}

# ---------------------------------------------
# MORTALITY (by toxicant only)
# ---------------------------------------------

mortality_hg <- run_pub_bias_analysis(mortality_data, "Mortality", toxicant = "Hg")
if (!is.null(mortality_hg)) {
  all_results[["Mortality_Hg"]] <- mortality_hg
  create_funnel_plot(mortality_hg, "figures/funnel_mortality_hg.png")
}

mortality_pcb <- run_pub_bias_analysis(mortality_data, "Mortality", toxicant = "PCBs")
if (!is.null(mortality_pcb)) {
  all_results[["Mortality_PCBs"]] <- mortality_pcb
  create_funnel_plot(mortality_pcb, "figures/funnel_mortality_pcbs.png")
}

# ---------------------------------------------
# GROWTH (by toxicant only)
# ---------------------------------------------

growth_hg <- run_pub_bias_analysis(growth_data, "Growth", toxicant = "Hg")
if (!is.null(growth_hg)) {
  all_results[["Growth_Hg"]] <- growth_hg
  create_funnel_plot(growth_hg, "figures/funnel_growth_hg.png")
}

growth_pcb <- run_pub_bias_analysis(growth_data, "Growth", toxicant = "PCBs")
if (!is.null(growth_pcb)) {
  all_results[["Growth_PCBs"]] <- growth_pcb
  create_funnel_plot(growth_pcb, "figures/funnel_growth_pcbs.png")
}

# ---------------------------------------------
# SIZE (Maternal - Fisher's Z, by toxicant only)
# ---------------------------------------------

size_hg <- run_pub_bias_analysis(size_data, "Maternal Size", toxicant = "Hg", effect_type = "Fisher's Z")
if (!is.null(size_hg)) {
  all_results[["Size_Hg"]] <- size_hg
  create_funnel_plot(size_hg, "figures/funnel_size_hg.png")
}

size_pcb <- run_pub_bias_analysis(size_data, "Maternal Size", toxicant = "PCBs", effect_type = "Fisher's Z")
if (!is.null(size_pcb)) {
  all_results[["Size_PCBs"]] <- size_pcb
  create_funnel_plot(size_pcb, "figures/funnel_size_pcbs.png")
}

# ---------------------------------------------
# MUSCLE (Maternal - Fisher's Z, by toxicant only)
# ---------------------------------------------

muscle_hg <- run_pub_bias_analysis(muscle_data, "Maternal Muscle", toxicant = "Hg", effect_type = "Fisher's Z")
if (!is.null(muscle_hg)) {
  all_results[["Muscle_Hg"]] <- muscle_hg
  create_funnel_plot(muscle_hg, "figures/funnel_muscle_hg.png")
}

muscle_pcb <- run_pub_bias_analysis(muscle_data, "Maternal Muscle", toxicant = "PCBs", effect_type = "Fisher's Z")
if (!is.null(muscle_pcb)) {
  all_results[["Muscle_PCBs"]] <- muscle_pcb
  create_funnel_plot(muscle_pcb, "figures/funnel_muscle_pcbs.png")
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n", sep = "")
cat("SUMMARY TABLE: PUBLICATION BIAS RESULTS\n")
cat(strrep("=", 70), "\n\n")

for (name in names(all_results)) {
  res <- all_results[[name]]
  if (!is.null(res)) {
    summary_table <- rbind(summary_table, data.frame(
      Analysis = res$label,
      k = res$k,
      Original_Effect = round(res$original_effect, 4),
      Original_CI = paste0("[", round(res$original_ci_lb, 4), ", ",
                           round(res$original_ci_ub, 4), "]"),
      Original_I2 = paste0(round(res$original_I2, 1), "%"),
      Egger_p = round(res$egger_pval, 4),
      k0_Imputed = res$k0_imputed,
      Adjusted_Effect = round(res$adjusted_effect, 4),
      Adjusted_CI = paste0("[", round(res$adjusted_ci_lb, 4), ", ",
                           round(res$adjusted_ci_ub, 4), "]"),
      Adjusted_I2 = paste0(round(res$adjusted_I2, 1), "%"),
      Change_Pct = round((res$adjusted_effect - res$original_effect) /
                           abs(res$original_effect) * 100, 1),
      stringsAsFactors = FALSE
    ))
  }
}

print(summary_table, row.names = FALSE)

# Save summary table
write.csv(summary_table, "publication_bias_summary.csv", row.names = FALSE)
cat("\n\nSummary table saved to: publication_bias_summary.csv\n")

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n", sep = "")
cat("INTERPRETATION GUIDE\n")
cat(strrep("=", 70), "\n\n")

cat("Egger's Test:\n")
cat("  - p < 0.05: Significant funnel plot asymmetry (potential bias)\n")
cat("  - p < 0.10: Marginal asymmetry (interpret with caution)\n")
cat("  - p >= 0.10: No significant asymmetry\n")
cat("  - Note: Low power with k < 10 studies\n\n")

cat("Trim-and-Fill (k0):\n")
cat("  - k0 = 0: No missing studies imputed\n")
cat("  - k0 > 0: Number of hypothetically missing studies\n\n")

cat("Effect Size Change:\n")
cat("  - < 10%: Minimal change, results robust\n")
cat("  - 10-25%: Moderate change, some bias possible\n")
cat("  - > 25%: Substantial change, bias likely\n\n")

# =============================================================================
# MULTI-PANEL FUNNEL PLOTS BY TOXICANT - Hg
# =============================================================================

cat("Creating multi-panel funnel plots by toxicant...\n\n")

# Helper function to get k0 value (returns 0 if NULL)
get_k0 <- function(model) {
  if (is.null(model$k0)) return(0)
  return(model$k0)
}

Leng_hg_tf <- all_results[["Size_Hg"]]$trimfill_model
Mus_hg_tf <- all_results[["Muscle_Hg"]]$trimfill_model
H_hg_tf <- all_results[["Hatching_Hg"]]$trimfill_model
Mor_hg_tf <- all_results[["Mortality_Hg"]]$trimfill_model
G_hg_tf <- all_results[["Growth_Hg"]]$trimfill_model

png("figures/trimfill_funnel_plots_Hg.png", width = 10, height = 12, units = "in", res = 900)

par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))

# Size Hg
funnel(Leng_hg_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Fisher's Z)",
       main = paste0("A) Maternal Size - Hg (k0 = ", get_k0(Leng_hg_tf), ")"))

# Muscle Hg
funnel(Mus_hg_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Fisher's Z)",
       main = paste0("B) Maternal Muscle - Hg (k0 = ", get_k0(Mus_hg_tf), ")"))

# Hatching Hg
funnel(H_hg_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Hedges' d)",
       main = paste0("C) Hatching - Hg (k0 = ", get_k0(H_hg_tf), ")"))

# Mortality Hg
funnel(Mor_hg_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Hedges' d)",
       main = paste0("D) Mortality - Hg (k0 = ", get_k0(Mor_hg_tf), ")"))

# Growth Hg
funnel(G_hg_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Hedges' d)",
       main = paste0("E) Growth - Hg (k0 = ", get_k0(G_hg_tf), ")"))

dev.off()
cat("Multi-panel funnel plot saved to: figures/trimfill_funnel_plots_Hg.png\n")

# =============================================================================
# MULTI-PANEL FUNNEL PLOTS BY TOXICANT - PCBs
# =============================================================================

Leng_pcb_tf <- all_results[["Size_PCBs"]]$trimfill_model
Mus_pcb_tf <- all_results[["Muscle_PCBs"]]$trimfill_model
H_pcb_tf <- all_results[["Hatching_PCBs"]]$trimfill_model
Mor_pcb_tf <- all_results[["Mortality_PCBs"]]$trimfill_model
G_pcb_tf <- all_results[["Growth_PCBs"]]$trimfill_model

png("figures/trimfill_funnel_plots_PCBs.png", width = 10, height = 12, units = "in", res = 900)

par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))

# Size PCBs
funnel(Leng_pcb_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Fisher's Z)",
       main = paste0("A) Maternal Size - PCBs (k0 = ", get_k0(Leng_pcb_tf), ")"))

# Muscle PCBs
funnel(Mus_pcb_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Fisher's Z)",
       main = paste0("B) Maternal Muscle - PCBs (k0 = ", get_k0(Mus_pcb_tf), ")"))

# Hatching PCBs
funnel(H_pcb_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Hedges' d)",
       main = paste0("C) Hatching - PCBs (k0 = ", get_k0(H_pcb_tf), ")"))

# Mortality PCBs
funnel(Mor_pcb_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Hedges' d)",
       main = paste0("D) Mortality - PCBs (k0 = ", get_k0(Mor_pcb_tf), ")"))

# Growth PCBs
funnel(G_pcb_tf, yaxis = "seinv",
       ylab = "Precision (1/SE)", xlab = "Effect Size (Hedges' d)",
       main = paste0("E) Growth - PCBs (k0 = ", get_k0(G_pcb_tf), ")"))

dev.off()
cat("Multi-panel funnel plot saved to: figures/trimfill_funnel_plots_PCBs.png\n")

# Reset par
par(mfrow = c(1, 1))

cat("\nAnalysis complete!\n")
