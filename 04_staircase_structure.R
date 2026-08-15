# Analysis 5: Staircase structure validation
#
# Purpose: Formalize and validate the staircase structure of s*(L) at fixed
#          f_min. Because k*(L) = ceil(log4(100L)) is a ceiling function,
#          s*(L) is piecewise constant within each Fofanov k-band and jumps
#          discontinuously at each band edge, this analysis (a) resolves
#          those edges precisely on a fine L grid, (b) checks them against
#          the closed-form prediction L = 4^k / 100, and (c) checks the
#          step-height ratio at each jump against the small-J asymptotic
#          prediction 1/f_min.

library(tidyverse)

options(tibble.print_max = Inf, tibble.width = Inf)
print_full <- function(x) print(dplyr::as_tibble(x))

RESULTS_DIR <- "results"
RDS_DIR <- "rds"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)


# 1. Core formulas

k_fofanov <- function(L) ceiling(log(100 * L) / log(4))
J_from_f <- function(f, k) f^k / (2 - f^k)
s_star_wald <- function(J, rho_star = 0.10, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  z^2 * (1 - J) / (rho_star^2 * J)
}

# 2. Fine L grid

L_min <- 1e5
L_max <- 1e10
n_grid_points <- 2000

L_grid <- 10^seq(log10(L_min), log10(L_max), length.out = n_grid_points)

f_values <- c(0.999, 0.99, 0.95, 0.90, 0.80)
rho_star <- 0.10
alpha <- 0.05

grid <- expand_grid(L = L_grid, f_min = f_values) |>
  mutate(
    k      = k_fofanov(L),
    J      = J_from_f(f_min, k),
    s_star = s_star_wald(J, rho_star, alpha)
  )

# 3. Step edge detection

detect_steps <- function(df_one_fmin) {
  df_one_fmin <- df_one_fmin |> dplyr::arrange(L)
  k_vec <- df_one_fmin$k
  L_vec <- df_one_fmin$L
  s_vec <- df_one_fmin$s_star

  step_idx <- which(diff(k_vec) != 0)

  if (length(step_idx) == 0) {
    return(tibble::tibble())
  }

  tibble(
    k_before = k_vec[step_idx],
    k_after = k_vec[step_idx + 1],
    L_before = L_vec[step_idx],
    L_after = L_vec[step_idx + 1],
    s_star_before = s_vec[step_idx],
    s_star_after = s_vec[step_idx + 1],
    L_predicted = 4^k_before / 100
  ) |>
    mutate(
      bracket_ok = L_predicted >= L_before & L_predicted <= L_after,
      bracket_relative_width_pct = 100 * (L_after - L_before) / L_predicted,
      observed_ratio = s_star_after / s_star_before,
      predicted_ratio = 1 / unique(df_one_fmin$f_min)[1]
    )
}

steps_by_fmin <- grid |>
  group_by(f_min) |>
  group_modify(~ detect_steps(.x), .keep = TRUE) |>
  ungroup() |>
  mutate(
    ratio_relative_error_pct = 100 * (observed_ratio / predicted_ratio - 1)
  )


# 4. Report: crossing-point validation and step-ratio validation

cat("[04] === Step edge crossing-point validation ===\n")
cat(sprintf("[04] Total steps detected across all f_min: %d\n", nrow(steps_by_fmin)))
cat(sprintf(
  "[04] Steps where predicted crossing falls inside the grid bracket: %d / %d\n",
  sum(steps_by_fmin$bracket_ok), nrow(steps_by_fmin)
))
cat(sprintf(
  "[04] Max bracket width (grid resolution) as %% of crossing point: %.4f%%\n",
  max(steps_by_fmin$bracket_relative_width_pct)
))

if (!all(steps_by_fmin$bracket_ok)) {
  cat("[04] WARNING: some predicted crossings fall OUTSIDE the grid bracket -- check grid density.\n")
  print_full(steps_by_fmin |> filter(!bracket_ok))
}

cat("[04] === Step-height ratio validation (observed s* ratio vs predicted 1/f_min) ===\n")
print_full(
  steps_by_fmin |>
    select(
      f_min, k_before, k_after, L_predicted, observed_ratio,
      predicted_ratio, ratio_relative_error_pct
    )
)

cat("[04] === Ratio error summary by f_min ===\n")
ratio_error_summary <- steps_by_fmin |>
  group_by(f_min) |>
  summarise(
    n_steps = n(),
    median_ratio_error_pct = median(ratio_relative_error_pct),
    max_ratio_error_pct = max(abs(ratio_relative_error_pct)),
    .groups = "drop"
  )
print_full(ratio_error_summary)

# 5. Save outputs

write_csv(grid, file.path(RESULTS_DIR, "analysis5_staircase_grid.csv"))
write_csv(steps_by_fmin, file.path(RESULTS_DIR, "analysis5_step_validation.csv"))

saveRDS(
  list(
    grid = grid,
    steps_by_fmin = steps_by_fmin,
    ratio_error_summary = ratio_error_summary,
    L_max = L_max
  ),
  file.path(RDS_DIR, "05_staircase_structure.rds")
)

cat("\n[04] === Analysis 4 complete ===\n")
cat("[04] Outputs written:\n")
cat("[04]   - results/analysis5_staircase_grid.csv\n")
cat("[04]   - results/analysis5_step_validation.csv\n")
cat("[04]   - rds/05_staircase_structure.rds\n")
