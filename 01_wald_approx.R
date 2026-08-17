# Analysis 2: Wald approximation vs. exact Clopper-Pearson
# inversion for minimum sketch size s*

library(tidyverse)
library(furrr)
library(future)

options(tibble.print_max = Inf, tibble.width = Inf)
print_full <- function(x) print(dplyr::as_tibble(x))

RESULTS_DIR <- "results"
RDS_DIR <- "rds"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Closed-form Wald estimator
s_star_wald <- function(J, rho_star = 0.10, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  z^2 * (1 - J) / (rho_star^2 * J)
}


# 2. Vectorized exact relative half-width rho(s) for a VECTOR of sketch
#    sizes s, at fixed J and alpha. Degenerate draws (x=0 or x=s) return NA.
rho_exact_vec <- function(s_vec, J, alpha = 0.05) {
  x <- round(J * s_vec)
  degenerate <- (x <= 0) | (x >= s_vec)

  J_lo <- rep(NA_real_, length(s_vec))
  J_hi <- rep(NA_real_, length(s_vec))

  ok <- !degenerate
  J_lo[ok] <- qbeta(alpha / 2, x[ok], s_vec[ok] - x[ok] + 1)
  J_hi[ok] <- qbeta(1 - alpha / 2, x[ok] + 1, s_vec[ok] - x[ok])

  J_hat <- x / s_vec
  rho <- ifelse(ok, (J_hi - J_lo) / (2 * J_hat), NA_real_)
  rho
}


# 3. Exact first-crossing minimum sketch size, s*_exact, via full
#    vectorized scan
#
#    s*_exact  = min{s : rho_exact(s,J,alpha) <= rho_star} (TRUE
#                minimum as formally defined, no monotonicity assumed)
#    s*_stable = 1 + (largest s at which the criterion FAILS) (the
#                smallest s after which it holds for every larger s
#                within the scanned range)
#
#    By construction from a single scan, s*_stable >= s*_exact always,
#    with equality iff the criterion, once first satisfied, is never
#    again violated at any larger scanned s (i.e. genuinely monotone
#    behavior in the relevant sense).
#
#    Scan range starts at a generous multiple of s_wald (cheap, and
#    almost certainly sufficient given how closely s*_exact tracks
#    s*_wald elsewhere in this analysis) and doubles up to s_max if no
#    success is found.

s_star_exact_and_stable <- function(
  J,
  rho_star = 0.10,
  alpha = 0.05,
  s_min = 1,
  s_max = 5e6,
  initial_multiple = 5
) {
  s_wald <- s_star_wald(J, rho_star, alpha)
  ceiling_try <- min(s_max, max(100, ceiling(s_wald * initial_multiple)))

  found <- FALSE
  s_exact <- NA_integer_
  s_stable <- NA_integer_
  any_nonmonotonic <- NA

  while (!found) {
    s_vec <- s_min:ceiling_try
    rho_vec <- rho_exact_vec(s_vec, J, alpha)
    success <- !is.na(rho_vec) & (rho_vec <= rho_star)

    if (any(success)) {
      found <- TRUE
      first_success_idx <- which(success)[1]
      s_exact <- s_vec[first_success_idx]

      failure_idx <- which(!success)
      # only failures AFTER the first success are relevant to stability
      failures_after <- failure_idx[failure_idx > first_success_idx]

      if (length(failures_after) == 0) {
        s_stable <- s_exact
        any_nonmonotonic <- FALSE
      } else {
        last_failure_idx <- max(failures_after)
        s_stable <- s_vec[last_failure_idx] + 1
        any_nonmonotonic <- TRUE
        # If the scan window ended exactly at a failure (last element
        # failed), s_stable's stability is unconfirmed within this
        # window
        if (last_failure_idx == length(s_vec) && ceiling_try < s_max) {
          ceiling_try <- min(s_max, ceiling_try * 2)
          found <- FALSE # rescan with extended window
        }
      }
    } else {
      if (ceiling_try >= s_max) {
        found <- TRUE # give up: infeasible within s_max
        s_exact <- NA_integer_
        s_stable <- NA_integer_
        any_nonmonotonic <- NA
      } else {
        ceiling_try <- min(s_max, ceiling_try * 2)
      }
    }
  }

  tibble::tibble(
    s_exact = s_exact,
    s_stable = s_stable,
    nonmonotonic = any_nonmonotonic,
    scan_ceiling_used = ceiling_try
  )
}

# 4. Parallel grid evaluation
n_workers <- parallel::detectCores()
plan(multisession, workers = n_workers)
cat(sprintf("[01] Using %d parallel workers (all detected cores)\n", n_workers))

grid <- expand_grid(
  J = 10^seq(log10(1e-4), log10(0.9), length.out = 40),
  rho_star = c(0.05, 0.10, 0.30),
  alpha = c(0.01, 0.05, 0.10)
) |>
  mutate(cell_id = row_number())

cat(sprintf("[01] Evaluating %d grid cells...\n", nrow(grid)))
t0 <- Sys.time()

exact_results <- future_pmap_dfr(
  list(grid$J, grid$rho_star, grid$alpha),
  ~ s_star_exact_and_stable(..1, ..2, ..3),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

cat(sprintf(
  "[01] Grid evaluation: %.1fs\n\n",
  as.numeric(Sys.time() - t0, units = "secs")
))

results <- bind_cols(grid, exact_results) |>
  mutate(
    s_wald = s_star_wald(J, rho_star, alpha),
    ratio_exact_wald = s_exact / s_wald,
    rel_error_exact_pct = 100 * (ratio_exact_wald - 1),
    stable_exceeds_exact = s_stable > s_exact,
    stable_exact_gap = s_stable - s_exact
  )

write_csv(results, file.path(RESULTS_DIR, "analysis2_wald_vs_exact_grid.csv"))

# 5. Summary: Wald vs. exact (s*_exact, the true minimum)

cat("[01] === s*_exact (true first crossing) vs s*_wald ===\n")
summary_table <- results |>
  filter(!is.na(s_exact)) |>
  group_by(rho_star, alpha) |>
  summarise(
    n_cells = dplyr::n(),
    median_ratio = median(ratio_exact_wald),
    max_ratio = max(ratio_exact_wald),
    min_ratio = min(ratio_exact_wald),
    median_relerr = median(rel_error_exact_pct),
    max_relerr = max(rel_error_exact_pct),
    .groups = "drop"
  ) |>
  arrange(rho_star, alpha)
print_full(summary_table)
write_csv(summary_table, file.path(RESULTS_DIR, "analysis2_summary_table.csv"))

n_na <- sum(is.na(results$s_exact))
if (n_na > 0) {
  message(sprintf(
    "[01] %d of %d grid cells returned NA for s*_exact (infeasible within s_max = 5e6).",
    n_na,
    nrow(results)
  ))
}

# 6. Non-monotonicity diagnostic
cat("\n[01] === s*_exact vs s*_stable: how often do they differ? ===\n")
valid <- results |> filter(!is.na(s_exact))

n_differ <- sum(valid$stable_exceeds_exact)
pct_coincide <- 100 * (1 - n_differ / nrow(valid))

cat(sprintf(
  "[01] s*_exact == s*_stable in %d / %d cells (%.1f%% coincide)\n",
  nrow(valid) - n_differ,
  nrow(valid),
  pct_coincide
))

nonmono_detail <- NULL
if (n_differ > 0) {
  cat(sprintf(
    "[01] Where they differ: gap ranges from %d to %d integer steps (median %d)\n",
    min(valid$stable_exact_gap[valid$stable_exceeds_exact]),
    max(valid$stable_exact_gap[valid$stable_exceeds_exact]),
    median(valid$stable_exact_gap[valid$stable_exceeds_exact])
  ))

  nonmono_detail <- valid |>
    filter(stable_exceeds_exact) |>
    select(J, rho_star, alpha, s_exact, s_stable, stable_exact_gap) |>
    arrange(desc(stable_exact_gap))
  print_full(nonmono_detail)
}

# 7. Checkpoint
saveRDS(
  list(
    results = results,
    summary_table = summary_table,
    nonmono_detail = nonmono_detail,
    pct_coincide = pct_coincide
  ),
  file.path(RDS_DIR, "02_wald_approx.rds")
)

cat("\n[01] === Analysis 1 complete ===\n")
cat(
  "[01] Outputs: results/analysis2_wald_vs_exact_grid.csv, results/analysis2_summary_table.csv, rds/02_wald_approx.rds\n"
)
