# =============================================================================
# Borenstein Correction for Shared Control Groups
# Based on Song et al. (2022) Equations 1-2
# For meta-analysis of Hg and PCB effects on fish reproductive success
# =============================================================================

library(readxl)
library(dplyr)
library(MASS)  # for ginv (generalized inverse)

# -----------------------------------------------------------------------------
# FUNCTION: Apply Borenstein correction to aggregate effect sizes within studies
# -----------------------------------------------------------------------------
#
# When multiple treatment groups share one control group, effect sizes are
# correlated. This function:
# 1. Builds a variance-covariance matrix using assumed correlation r
# 2. Calculates aggregated (pooled) effect size per study
# 3. Calculates pooled variance accounting for covariance
#
# Borenstein et al. (2009) Equations:
#   Eq 1: Covariance = r * sqrt(vi) * sqrt(vj)
#   Eq 2: Pooled variance = (1/k^2) * sum(V), where V is the var-cov matrix
#         Pooled effect = (1/k) * sum(yi)  [unweighted mean]

borenstein_aggregate <- function(yi, vi, nc = NULL, nt = NULL, r_default = 0.5) {
  # yi: vector of effect sizes from same study
  # vi: vector of corresponding sampling variances
  # nc: control group sample size (scalar, same for all effect sizes in study)
  # nt: vector of treatment group sample sizes
  # r_default: fallback correlation if sample sizes unavailable (default 0.5)

  k <- length(yi)

  if (k == 1) {
    # Single effect size - no correction needed
    return(list(yi_pooled = yi, vi_pooled = vi, r_used = NA, r_method = "none"))
  }

  # Build variance-covariance matrix (Equation 1)
  V <- matrix(0, nrow = k, ncol = k)

  # Determine if we can calculate r from sample sizes
  can_calc_r <- !is.null(nc) && !is.null(nt) &&
                !is.na(nc) && all(!is.na(nt)) &&
                nc > 0 && all(nt > 0)

  r_method <- ifelse(can_calc_r, "calculated", "assumed")
  r_values <- matrix(NA, nrow = k, ncol = k)

  for (i in 1:k) {
    for (j in 1:k) {
      if (i == j) {
        V[i, j] <- vi[i]  # Diagonal: variance
      } else {
        # Calculate correlation between effect sizes i and j
        if (can_calc_r) {
          # Formula: r_ij = nc / sqrt((nc + nt_i) * (nc + nt_j))
          # From Gleser & Olkin (2009), Borenstein et al. (2009)
          r_ij <- nc / sqrt((nc + nt[i]) * (nc + nt[j]))
        } else {
          r_ij <- r_default
        }
        r_values[i, j] <- r_ij
        V[i, j] <- r_ij * sqrt(vi[i]) * sqrt(vi[j])  # Off-diagonal: covariance
      }
    }
  }

  # Calculate pooled effect size (unweighted mean - Borenstein et al. 2009)
  yi_pooled <- mean(yi)

  # Calculate pooled variance (Equation 2)
  # vi_pooled = (1/k^2) * sum of all elements in V
  vi_pooled <- sum(V) / (k^2)

  # Mean r used (for reporting)
  r_used <- mean(r_values[upper.tri(r_values)], na.rm = TRUE)

  return(list(yi_pooled = yi_pooled, vi_pooled = vi_pooled,
              r_used = r_used, r_method = r_method))
}

# -----------------------------------------------------------------------------
# FUNCTION: Apply correction to entire dataset
# -----------------------------------------------------------------------------

apply_borenstein_correction <- function(data, study_col = "study_id",
                                         yi_col = "yi", vi_col = "vi",
                                         toxicant_col = "toxicant",
                                         nc_col = "control_n", nt_col = "treated_n",
                                         r_default = 0.5) {
  # Group by study and apply correction
  studies <- unique(data[[study_col]])

  results <- data.frame(
    study_id = character(),
    toxicant = character(),
    yi_original = character(),
    vi_original = character(),
    k = integer(),  # number of effect sizes in study
    yi_corrected = numeric(),
    vi_corrected = numeric(),
    r_used = numeric(),
    r_method = character(),
    needs_correction = logical(),
    stringsAsFactors = FALSE
  )

  for (study in studies) {
    subset_data <- data[data[[study_col]] == study, ]
    yi_vec <- subset_data[[yi_col]]
    vi_vec <- subset_data[[vi_col]]

    # Get sample sizes
    nc <- unique(subset_data[[nc_col]])[1]  # Control n (same for all in study)
    nt_vec <- subset_data[[nt_col]]          # Treatment n (may vary)

    # Get toxicant (should be same for all rows in a study)
    toxicant <- unique(subset_data[[toxicant_col]])[1]
    k <- length(yi_vec)

    correction <- borenstein_aggregate(yi_vec, vi_vec, nc = nc, nt = nt_vec,
                                        r_default = r_default)

    results <- rbind(results, data.frame(
      study_id = study,
      toxicant = toxicant,
      yi_original = paste(round(yi_vec, 3), collapse = "; "),
      vi_original = paste(round(vi_vec, 4), collapse = "; "),
      k = k,
      yi_corrected = correction$yi_pooled,
      vi_corrected = correction$vi_pooled,
      r_used = ifelse(is.na(correction$r_used), NA, round(correction$r_used, 3)),
      r_method = correction$r_method,
      needs_correction = (k > 1),
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

# -----------------------------------------------------------------------------
# FUNCTION: Calculate meta-analytic weighted mean and CI
# -----------------------------------------------------------------------------

calculate_meta_results <- function(yi, vi, label = "") {
  # Inverse-variance weights
  wi <- 1 / vi

  # Weighted mean
  weighted_mean <- sum(wi * yi) / sum(wi)

  # Variance of weighted mean
  var_weighted_mean <- 1 / sum(wi)

  # Standard error
  se <- sqrt(var_weighted_mean)

  # 95% CI
  ci_lower <- weighted_mean - 1.96 * se
  ci_upper <- weighted_mean + 1.96 * se

  # Heterogeneity (Q statistic)
  Q <- sum(wi * (yi - weighted_mean)^2)
  df <- length(yi) - 1
  p_het <- 1 - pchisq(Q, df)

  # I-squared
  I2 <- max(0, (Q - df) / Q * 100)

  cat("\n", strrep("=", 60), "\n", sep="")
  cat(label, "\n")
  cat(strrep("=", 60), "\n", sep="")
  cat("Number of studies (k):", length(yi), "\n")
  cat("Weighted mean effect (Hedge's d):", round(weighted_mean, 4), "\n")
  cat("Standard error:", round(se, 4), "\n")
  cat("95% CI: [", round(ci_lower, 4), ", ", round(ci_upper, 4), "]\n", sep="")
  cat("Q statistic:", round(Q, 2), "(df =", df, ", p =", round(p_het, 4), ")\n")
  cat("I-squared:", round(I2, 1), "%\n")

  return(list(
    k = length(yi),
    weighted_mean = weighted_mean,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    Q = Q,
    p_het = p_het,
    I2 = I2
  ))
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

# Set working directory to the repo root before sourcing this script.
# In RStudio: Session > Set Working Directory > To Project Directory.
# Otherwise run: setwd("path/to/Hg_PCBs_Meta_Repo")

# Read data
df <- read_excel("data/ICC List for Meta.xlsx", sheet = "Sheet1")
names(df) <- c("meta_analysis", "study_id", "toxicant", "yi", "vi",
               "treated_conc", "control_n", "treated_n")

# Define which analyses to run (only experimental studies with shared controls)
# Size and Muscle are correlation-based (Fisher's Z) - no shared control issue
target_analyses <- c("Hatching", "Mortality", "Growth")

# Correlation for shared control (standard assumption)
r_assumed <- 0.5

cat("\n")
cat(strrep("*", 70), "\n", sep="")
cat("BORENSTEIN CORRECTION FOR SHARED CONTROL GROUPS\n")
cat("Assumed correlation (r) =", r_assumed, "\n")
cat("Based on Song et al. (2022) Equations 1-2\n")
cat(strrep("*", 70), "\n", sep="")

# Store all results
all_results <- list()
comparison_table <- data.frame()

for (analysis in target_analyses) {

  cat("\n\n")
  cat(strrep("#", 70), "\n", sep="")
  cat("ANALYSIS:", analysis, "\n")
  cat(strrep("#", 70), "\n", sep="")

  # Subset data for this analysis
  analysis_data <- df[df$meta_analysis == analysis, ]

  cat("\nData summary:\n")
  cat("  Total effect sizes:", nrow(analysis_data), "\n")
  cat("  Unique studies:", length(unique(analysis_data$study_id)), "\n")

  # Count studies needing correction
  study_counts <- table(analysis_data$study_id)
  multi_effect_studies <- sum(study_counts > 1)
  cat("  Studies with >1 effect size (need correction):", multi_effect_studies, "\n")

  # ---------------------------------------------
  # UNCORRECTED ANALYSIS (for comparison)
  # ---------------------------------------------
  cat("\n--- UNCORRECTED (naive) analysis ---\n")
  uncorrected <- calculate_meta_results(
    analysis_data$yi,
    analysis_data$vi,
    label = paste(analysis, "- UNCORRECTED (treating all effect sizes as independent)")
  )

  # ---------------------------------------------
  # APPLY BORENSTEIN CORRECTION
  # ---------------------------------------------
  cat("\n--- Applying Borenstein correction ---\n")

  corrected_data <- apply_borenstein_correction(
    analysis_data,
    study_col = "study_id",
    yi_col = "yi",
    vi_col = "vi",
    toxicant_col = "toxicant",
    nc_col = "control_n",
    nt_col = "treated_n",
    r_default = r_assumed
  )

  # Show which studies were corrected
  cat("\nStudies requiring correction (k > 1):\n")
  corrected_studies <- corrected_data[corrected_data$needs_correction, ]
  if (nrow(corrected_studies) > 0) {
    for (i in 1:nrow(corrected_studies)) {
      cat("  ", corrected_studies$study_id[i], "\n", sep="")
      cat("    Original yi: ", corrected_studies$yi_original[i], "\n", sep="")
      cat("    Original vi: ", corrected_studies$vi_original[i], "\n", sep="")
      cat("    k =", corrected_studies$k[i], "\n")
      cat("    r =", corrected_studies$r_used[i],
          "(", corrected_studies$r_method[i], ")\n", sep="")
      cat("    Pooled yi:", round(corrected_studies$yi_corrected[i], 4), "\n")
      cat("    Pooled vi:", round(corrected_studies$vi_corrected[i], 4), "\n\n")
    }
  }

  # Summary of r methods used
  calc_count <- sum(corrected_data$r_method == "calculated", na.rm = TRUE)
  assumed_count <- sum(corrected_data$r_method == "assumed", na.rm = TRUE)
  cat("Correlation (r) summary:\n")
  cat("  Studies with r calculated from sample sizes:", calc_count, "\n")
  cat("  Studies with r assumed (", r_assumed, "):", assumed_count, "\n", sep="")

  # ---------------------------------------------
  # CORRECTED META-ANALYSIS (OVERALL)
  # ---------------------------------------------
  corrected_results <- calculate_meta_results(
    corrected_data$yi_corrected,
    corrected_data$vi_corrected,
    label = paste(analysis, "- CORRECTED OVERALL (Borenstein, r =", r_assumed, ")")
  )

  # ---------------------------------------------
  # RESULTS BY TOXICANT (UNCORRECTED AND CORRECTED)
  # ---------------------------------------------
  cat("\n--- Results by Toxicant ---\n")

  toxicants <- unique(analysis_data$toxicant)
  toxicant_results <- list()

  for (tox in toxicants) {
    # UNCORRECTED for this toxicant
    tox_raw_data <- analysis_data[analysis_data$toxicant == tox, ]
    if (nrow(tox_raw_data) > 0) {
      tox_uncorrected <- calculate_meta_results(
        tox_raw_data$yi,
        tox_raw_data$vi,
        label = paste(analysis, "-", tox, "UNCORRECTED")
      )

      # Add uncorrected to comparison table
      comparison_table <- rbind(comparison_table, data.frame(
        Analysis = analysis,
        Toxicant = tox,
        Method = "Uncorrected",
        k = tox_uncorrected$k,
        Effect = round(tox_uncorrected$weighted_mean, 4),
        SE = round(tox_uncorrected$se, 4),
        CI_lower = round(tox_uncorrected$ci_lower, 4),
        CI_upper = round(tox_uncorrected$ci_upper, 4),
        I2 = round(tox_uncorrected$I2, 1)
      ))
    }

    # CORRECTED for this toxicant
    tox_data <- corrected_data[corrected_data$toxicant == tox, ]
    if (nrow(tox_data) > 0) {
      tox_results <- calculate_meta_results(
        tox_data$yi_corrected,
        tox_data$vi_corrected,
        label = paste(analysis, "-", tox, "CORRECTED")
      )
      toxicant_results[[tox]] <- tox_results

      # Add corrected to comparison table
      comparison_table <- rbind(comparison_table, data.frame(
        Analysis = analysis,
        Toxicant = tox,
        Method = "Corrected",
        k = tox_results$k,
        Effect = round(tox_results$weighted_mean, 4),
        SE = round(tox_results$se, 4),
        CI_lower = round(tox_results$ci_lower, 4),
        CI_upper = round(tox_results$ci_upper, 4),
        I2 = round(tox_results$I2, 1)
      ))
    }
  }

  # Store results
  all_results[[analysis]] <- list(
    uncorrected = uncorrected,
    corrected = corrected_results,
    corrected_data = corrected_data,
    by_toxicant = toxicant_results
  )

  # Add overall results to comparison table
  comparison_table <- rbind(comparison_table, data.frame(
    Analysis = analysis,
    Toxicant = "Overall",
    Method = "Uncorrected",
    k = uncorrected$k,
    Effect = round(uncorrected$weighted_mean, 4),
    SE = round(uncorrected$se, 4),
    CI_lower = round(uncorrected$ci_lower, 4),
    CI_upper = round(uncorrected$ci_upper, 4),
    I2 = round(uncorrected$I2, 1)
  ))

  comparison_table <- rbind(comparison_table, data.frame(
    Analysis = analysis,
    Toxicant = "Overall",
    Method = "Corrected",
    k = corrected_results$k,
    Effect = round(corrected_results$weighted_mean, 4),
    SE = round(corrected_results$se, 4),
    CI_lower = round(corrected_results$ci_lower, 4),
    CI_upper = round(corrected_results$ci_upper, 4),
    I2 = round(corrected_results$I2, 1)
  ))

  # Save corrected data to CSV
  write.csv(corrected_data,
            paste0("data/borenstein_corrected_", tolower(analysis), ".csv"),
            row.names = FALSE)
  cat("\nCorrected data saved to: data/borenstein_corrected_", tolower(analysis), ".csv\n", sep="")
}

# =============================================================================
# SUMMARY COMPARISON TABLE
# =============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n", sep="")
cat("SUMMARY: COMPARISON OF UNCORRECTED vs BORENSTEIN-CORRECTED RESULTS\n")
cat(strrep("=", 70), "\n\n", sep="")

print(comparison_table, row.names = FALSE)

cat("\n\nNotes:\n")
cat("- k = number of units (uncorrected: effect sizes; corrected: studies)\n")
cat("- Effect = weighted mean Hedge's d\n")
cat("- SE = standard error of the weighted mean\n")
cat("- CI = 95% confidence interval\n")
cat("- I2 = heterogeneity (percentage of variance due to true differences)\n")
cat("- Positive effects indicate increased hatching success/mortality/growth\n")
cat("  (check your coding direction for interpretation)\n")

# Save comparison table
write.csv(comparison_table, "borenstein_comparison_summary.csv", row.names = FALSE)
cat("\nComparison table saved to: borenstein_comparison_summary.csv\n")

cat("\n\nAnalysis complete!\n")
