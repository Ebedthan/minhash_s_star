# Analysis 3: Idealized binomial-sketch calibration
#
# Purpose: Confirm that IF x ~ Binomial(s*, J) holds exactly (i.e. under
#          the idealized sampling model with no k-mer non-independence,
#          no hashing artifacts, no sequence-composition effects), the
#          closed-form minimum sketch size s* actually delivers the target
#          relative CI half-width rho* and the target CP coverage 1-alpha.

library(tidyverse)

set.seed(42)

RESULTS_DIR <- "results"
RDS_DIR <- "rds"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Core formulas (Fofanov k, Mash J, closed-form s*)

k_fofanov <- function(L) ceiling(log(100 * L) / log(4))

J_from_f <- function(f, k) f^k / (2 - f^k)

# Unrounded, for direct numerical comparison against s*_exact only
# (Analysis 2) -- rounding either side would distort that comparison.
s_star_wald_real <- function(J, rho_star = 0.10, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  z^2 * (1 - J) / (rho_star^2 * J)
}

# Integer design quantity used everywhere a sketch is actually
# constructed or simulated (Analyses 3, 4, 6, and the closed-form
# formula as reported/recommended in the main text).
s_star_wald <- function(J, rho_star = 0.10, alpha = 0.05) {
  ceiling(s_star_wald_real(J, rho_star, alpha))
}

# 2. Vectorized Clopper-Pearson relative half-width and coverage check

#' Compute realized rho and CI-coverage indicator for a vector of draws
#'
#' rho (relative CI half-width) is left NA for degenerate draws (x = 0 or
#' x = s), since it is genuinely undefined there (J_hat = 0 makes the
#' relative half-width blow up or the interval collapse asymmetrically).
#'
#' coverage is computed for ALL draws, including degenerate ones, using
#' the standard one-sided Clopper-Pearson convention:
#'   x = 0  ->  J_lo = 0, J_hi = qbeta(1 - alpha/2, 1, s)
#'   x = s  ->  J_lo = qbeta(alpha/2, s, 1), J_hi = 1
#' This is required for the unconditional coverage to match Clopper-Pearson
#' theory's >= 1 - alpha guarantee; excluding degenerate draws computes
#' coverage conditional on being away from the boundary instead, which has
#' no such guarantee and can read well below nominal purely as an artifact
#' of the exclusion.
#'
#' @param x vector of observed shared-hash counts (integer, 0 <= x <= s)
#' @param s sketch size (scalar, same for all draws)
#' @param J_true true Jaccard index used to generate x (for coverage check)
#' @param alpha significance level
#' @return tibble with columns x, J_hat, J_lo, J_hi, rho, degenerate, covered
cp_diagnostics <- function(x, s, J_true, alpha = 0.05) {
  n <- length(x)
  J_hat <- x / s

  is_zero <- x == 0
  is_full <- x == s
  degenerate <- is_zero | is_full
  interior <- !degenerate

  J_lo <- rep(NA_real_, n)
  J_hi <- rep(NA_real_, n)

  # Interior draws: standard two-sided CP interval
  J_lo[interior] <- qbeta(alpha / 2, x[interior], s - x[interior] + 1)
  J_hi[interior] <- qbeta(1 - alpha / 2, x[interior] + 1, s - x[interior])

  # x = 0: one-sided upper bound only, J_lo = 0
  if (any(is_zero)) {
    J_lo[is_zero] <- 0
    J_hi[is_zero] <- qbeta(1 - alpha / 2, 1, s)
  }

  # x = s: one-sided lower bound only, J_hi = 1
  if (any(is_full)) {
    J_lo[is_full] <- qbeta(alpha / 2, s, 1)
    J_hi[is_full] <- 1
  }

  # rho stays NA for degenerate draws -- undefined relative half-width
  rho <- ifelse(interior, (J_hi - J_lo) / (2 * J_hat), NA_real_)

  # coverage is defined for every draw, degenerate or not
  covered <- (J_lo <= J_true) & (J_true <= J_hi)

  tibble::tibble(
    x = x, J_hat = J_hat, J_lo = J_lo, J_hi = J_hi,
    rho = rho, degenerate = degenerate, covered = covered
  )
}

# 3. Parameter grid

L_range <- 10^seq(5, 10, length.out = 30)
f_values <- c(0.999, 0.99, 0.95, 0.90, 0.80)

n_draws <- 10000 # per cell, matching the design target

grid <- expand_grid(L = L_range, f_min = f_values) |>
  mutate(
    k = k_fofanov(L),
    J = J_from_f(f_min, k),
    s_star = s_star_wald(J, rho_star = 0.10, alpha = 0.05)
  ) |>
  filter(J > 0, s_star >= 2) # guard against degenerate cells

# 4. Simulate: draw x ~ Binomial(s*, J) n_draws times per cell

raw_rho_store <- list() # populated as a side effect, keyed by cell id

simulate_cell <- function(cell_id, L, f_min, k, J, s_star,
                          rho_star = 0.10, alpha = 0.05,
                          store_raw = FALSE) {
  x <- rbinom(n_draws, size = s_star_wald(J, rho_star, alpha), prob = J)
  diag <- cp_diagnostics(x, s_star, J_true = J, alpha = alpha)

  if (store_raw) {
    raw_rho_store[[as.character(cell_id)]] <<- diag$rho
  }

  # For a continuous statistic, median_rho <= rho_star would imply
  # prop_rho_le_target >= 0.5 by definition of the median, and vice versa.
  # rho here is discrete (only as many distinct values as distinct x), so
  # exact 50% is not guaranteed, but a *directional* mismatch (median on
  # one side of rho_star while the proportion strongly implies the other)
  # signals something worth inspecting rather than dismissing as
  # discreteness noise.
  med_rho <- median(diag$rho, na.rm = TRUE)
  prop_le <- mean(diag$rho <= rho_star, na.rm = TRUE)

  median_says_under <- med_rho <= rho_star
  prop_says_majority_under <- prop_le >= 0.5
  consistency_flag <- median_says_under != prop_says_majority_under

  tibble::tibble(
    cell_id = cell_id,
    L = L, f_min = f_min, k = k, J = J, s_star = s_star,
    rho_star = rho_star, alpha = alpha,
    n_draws = n_draws,
    n_degenerate = sum(diag$degenerate),
    prop_rho_le_target = prop_le,
    median_rho = med_rho,
    mean_rho = mean(diag$rho, na.rm = TRUE),
    coverage = mean(diag$covered),
    median_prop_consistency_flag = consistency_flag
  )
}

grid <- grid |> dplyr::mutate(cell_id = row_number())

results <- purrr::pmap_dfr(
  list(grid$cell_id, grid$L, grid$f_min, grid$k, grid$J, grid$s_star),
  simulate_cell
)

# 5. Summary diagnostics

cat("[02] === Overall calibration summary (rho* = 0.10, alpha = 0.05) ===\n")
cat(sprintf("[02] Cells evaluated: %d\n", nrow(results)))
cat(sprintf(
  "[02] Median coverage across cells: %.4f (target: %.4f)\n",
  median(results$coverage), 1 - 0.05
))
cat(sprintf(
  "[02] Range of coverage: [%.4f, %.4f]\n",
  min(results$coverage), max(results$coverage)
))
cat(sprintf(
  "[02] Median realized rho across cells: %.4f (target rho* = 0.10)\n",
  median(results$median_rho)
))
cat(sprintf(
  "[02] Median proportion of draws with rho <= rho*: %.4f\n",
  median(results$prop_rho_le_target)
))

n_deg_total <- sum(results$n_degenerate)
if (n_deg_total > 0) {
  message(sprintf(
    "[02] %d degenerate draws (x = 0 or x = s) out of %d total across all cells (now included in coverage, excluded from rho).",
    n_deg_total, nrow(results) * n_draws
  ))
}

alpha <- 0.05
mc_tol <- 3 * sqrt(0.95 * 0.05 / n_draws) # ~3 SE under nominal binomial MC error

coverage_below_flags <- results |>
  dplyr::filter(coverage < (1 - alpha) - mc_tol) |>
  dplyr::select(cell_id, L, f_min, k, J, s_star, coverage)

if (nrow(coverage_below_flags) > 0) {
  cat("[02] === Cells with coverage BELOW nominal 1-alpha beyond MC tolerance (unexpected under CP theory) ===\n")
  print(coverage_below_flags, n = Inf)
} else {
  cat("[02] No cells show coverage below nominal beyond Monte Carlo tolerance -- consistent with Clopper-Pearson theory.\n")
}

coverage_above_flags <- results |>
  dplyr::filter(coverage > (1 - alpha) + 0.02) |>
  dplyr::select(cell_id, L, f_min, k, J, s_star, coverage, n_degenerate)

cat(sprintf(
  "[02] %d cells show coverage >2pp above nominal (expected CP conservatism, not a concern).\n",
  nrow(coverage_above_flags)
))

n_consistency_flags <- sum(results$median_prop_consistency_flag)
cat(sprintf(
  "[02] === Median/proportion consistency check: %d of %d cells flagged ===\n",
  n_consistency_flags, nrow(results)
))
if (n_consistency_flags > 0) {
  print(
    results |>
      dplyr::filter(median_prop_consistency_flag) |>
      dplyr::select(cell_id, L, f_min, s_star, J, median_rho, prop_rho_le_target, rho_star),
    n = Inf
  )
}

# 6. Sweep rho_star and alpha as well, at a fixed representative L, f_min

rep_L <- 1e6 # 1 Mb, matching the paper's anchor case
rep_f_min <- 0.90

rep_k <- k_fofanov(rep_L)
rep_J <- J_from_f(rep_f_min, rep_k)

param_grid <- expand_grid(
  rho_star = c(0.05, 0.10, 0.30),
  alpha = c(0.01, 0.05, 0.10)
) |>
  dplyr::mutate(
    s_star = s_star_wald(rep_J, rho_star, alpha),
    cell_id = paste0("param_", row_number())
  )

param_results <- pmap_dfr(
  list(param_grid$cell_id, param_grid$rho_star, param_grid$alpha, param_grid$s_star),
  function(cell_id, rho_star, alpha, s_star) {
    simulate_cell(cell_id, rep_L, rep_f_min, rep_k, rep_J, s_star,
      rho_star, alpha,
      store_raw = TRUE
    )
  }
)

cat("[02] === Calibration across (rho*, alpha) at L = 1 Mb, f_min = 0.90 ===\n")
print(
  param_results |>
    dplyr::select(
      cell_id, rho_star, alpha, s_star, coverage, median_rho,
      prop_rho_le_target, median_prop_consistency_flag
    ),
  n = Inf
)


# 7. Gap between rho* and achieved rho distribution (representative cell)

target_cell <- param_grid |>
  dplyr::filter(rho_star == 0.10, alpha == 0.05) |>
  dplyr::pull(cell_id)
rho_sample <- raw_rho_store[[target_cell]]

target_rho_star <- 0.10

cat("[02] === Gap between rho* and achieved rho distribution (representative cell) ===\n")
cat(sprintf("[02] rho* target: %.5f\n", target_rho_star))
cat(sprintf(
  "[02] Observed median rho: %.5f (%.2f%% above target)\n",
  median(rho_sample, na.rm = TRUE),
  100 * (median(rho_sample, na.rm = TRUE) / target_rho_star - 1)
))
cat(sprintf("[02] P(rho <= rho*): %.4f\n", mean(rho_sample <= target_rho_star, na.rm = TRUE)))

quantiles_check <- quantile(rho_sample, probs = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), na.rm = TRUE)
cat("[02] Empirical quantiles of realized rho:\n")
print(quantiles_check)

results <- results |>
  dplyr::mutate(median_rho_overshoot_pct = 100 * (median_rho / rho_star - 1))

cat("[02] === Median rho overshoot vs rho_star, across full (L, f_min) grid ===\n")
cat(sprintf("[02] Median overshoot: %.2f%%\n", median(results$median_rho_overshoot_pct, na.rm = TRUE)))
cat(sprintf(
  "[02] Range: [%.2f%%, %.2f%%]\n",
  min(results$median_rho_overshoot_pct, na.rm = TRUE),
  max(results$median_rho_overshoot_pct, na.rm = TRUE)
))

# 8. Save outputs

write_csv(results, file.path(RESULTS_DIR, "analysis3_L_fmin_calibration.csv"))
write_csv(
  param_results |> dplyr::select(-median_prop_consistency_flag) |>
    dplyr::bind_cols(param_results |> select(median_prop_consistency_flag)),
  file.path(RESULTS_DIR, "analysis3_rho_alpha_calibration.csv")
)

saveRDS(
  list(
    results = results,
    param_results = param_results,
    rho_sample_representative_cell = rho_sample,
    representative_cell_params = list(L = rep_L, f_min = rep_f_min, rho_star = 0.10, alpha = 0.05),
    quantiles_check = quantiles_check
  ),
  file.path(RDS_DIR, "03_idealized_binom.rds")
)

cat("\n[02] === Analysis 2 complete ===\n")
cat("[02] Outputs: results/analysis3_L_fmin_calibration.csv, results/analysis3_rho_alpha_calibration.csv, rds/03_idealized_binom.rds\n")
