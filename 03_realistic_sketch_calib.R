# Realistic-sequence sketch calibration
#
# Purpose:
# Test whether the idealized assumption behind s* (x ~ Binomial(s,J)
# exactly) survives contact with REAL overlapping k-mers and REAL
# hashing, rather than an idealized binomial draw (Analysis 3).

library(tidyverse)
library(digest)
library(furrr)
library(future)

set.seed(42)

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
  ceiling(z^2 * (1 - J) / (rho_star^2 * J))
}

# 2. Clopper-Pearson diagnostics (with the coverage fix from Analysis 3:
#    degenerate draws use the one-sided CP convention rather than being
#    excluded from coverage)

cp_diagnostics <- function(x, s, J_true, alpha = 0.05) {
  n <- length(x)
  J_hat <- x / s
  is_zero <- x == 0
  is_full <- x == s
  degenerate <- is_zero | is_full
  interior <- !degenerate

  J_lo <- rep(NA_real_, n)
  J_hi <- rep(NA_real_, n)

  if (any(interior)) {
    J_lo[interior] <- qbeta(alpha / 2, x[interior], s - x[interior] + 1)
    J_hi[interior] <- qbeta(1 - alpha / 2, x[interior] + 1, s - x[interior])
  }
  if (any(is_zero)) {
    J_lo[is_zero] <- 0
    J_hi[is_zero] <- qbeta(1 - alpha / 2, 1, s)
  }
  if (any(is_full)) {
    J_lo[is_full] <- qbeta(alpha / 2, s, 1)
    J_hi[is_full] <- 1
  }

  rho <- ifelse(interior, (J_hi - J_lo) / (2 * J_hat), NA_real_)
  covered <- (J_lo <= J_true) & (J_true <= J_hi)

  dplyr::tibble(
    x = x, J_hat = J_hat, J_lo = J_lo, J_hi = J_hi,
    rho = rho, degenerate = degenerate, covered = covered
  )
}

# 3. Sequence generation and mutation

generate_iid_sequence <- function(N) {
  sample.int(4L, N, replace = TRUE) - 1L
}

mutate_sequence <- function(seq_int, p) {
  N <- length(seq_int)
  mutate_mask <- runif(N) < p
  n_mut <- sum(mutate_mask)
  if (n_mut > 0) {
    offset <- sample.int(3L, n_mut, replace = TRUE)
    seq_int[mutate_mask] <- (seq_int[mutate_mask] + offset) %% 4L
  }
  seq_int
}

# 4. K-mer encoding and hashing
#    hash_codes(): MurmurHash3-based, via digest::digest2int() (a tested,
#    vectorized MurmurHash32 binding. the same hash family Mash and
#    sourmash use for k-mer hashing). Two independently seeded 32-bit
#    hashes are combined into a 64-bit-equivalent space to keep collision
#    probability negligible at the L used here.

kmer_codes <- function(seq_int, k) {
  N <- length(seq_int)
  n_kmers <- N - k + 1
  code <- numeric(n_kmers)
  for (j in 0:(k - 1)) {
    code <- code + seq_int[(1 + j):(n_kmers + j)] * 4^(k - 1 - j)
  }
  code
}

hash_codes <- function(codes) {
  codes_chr <- as.character(as.integer(codes))
  h_hi <- digest2int(codes_chr, seed = 1L)
  h_lo <- digest2int(codes_chr, seed = 2L)
  h_hi_u <- bitwAnd(h_hi, 0x7FFFFFFF)
  h_lo_u <- bitwAnd(h_lo, 0x7FFFFFFF)
  h_hi_u * 2^31 + h_lo_u
}


# 5. Bottom-s MinHash sketch construction and Mash-style pairwise comparison

build_sketch_from_hashes <- function(uh, s) {
  if (length(uh) < s) {
    return(list(sketch = sort(uh), full = FALSE))
  }
  list(sketch = sort(uh)[1:s], full = TRUE)
}

mash_compare <- function(sketchA, sketchB, s_target) {
  merged <- sort(unique(c(sketchA, sketchB)))
  s_used <- min(s_target, length(merged))
  merged_s <- merged[seq_len(s_used)]
  x <- sum(merged_s %in% sketchA & merged_s %in% sketchB)
  list(x = x, s_used = s_used)
}


# 6. STAGE 1, hash_codes() correctness tests. Hard stop if either fails.

cat("[03] === Stage 1: hash_codes() correctness tests ===\n")

test_hash_not_monotonic <- function() {
  codes <- 0:2000
  h <- hash_codes(codes)
  is_monotonic <- all(diff(h) > 0) || all(diff(h) < 0)
  list(pass = !is_monotonic, detail = sprintf(
    "Monotonic increasing/decreasing: %s (should be FALSE)", is_monotonic
  ))
}

test_hash_no_collisions_realistic <- function(n = 1e6) {
  codes <- sample.int(4^14, n) # realistic-scale, largest k used (k=14)
  h <- hash_codes(codes)
  n_unique_codes <- length(unique(codes))
  n_unique_hashes <- length(unique(h))
  list(
    pass = n_unique_hashes == n_unique_codes,
    detail = sprintf("Unique codes: %d, unique hashes: %d", n_unique_codes, n_unique_hashes)
  )
}

test_hash_deterministic <- function() {
  codes <- c(1L, 100L, 100000L)
  h1 <- hash_codes(codes)
  h2 <- hash_codes(codes)
  list(pass = identical(h1, h2), detail = "Same input -> same output, checked twice")
}

hash_tests <- list(
  "not monotonic" = test_hash_not_monotonic(),
  "no collisions (1e6)" = test_hash_no_collisions_realistic(),
  "deterministic" = test_hash_deterministic()
)

all_passed <- TRUE
for (name in names(hash_tests)) {
  t <- hash_tests[[name]]
  status <- if (t$pass) "PASS" else "FAIL"
  cat(sprintf("[03]   [%s] %s -- %s\n", status, name, t$detail))
  if (!t$pass) all_passed <- FALSE
}

if (!all_passed) {
  stop(
    "[03] hash_codes() failed one or more correctness tests -- fix before proceeding. ",
    "Do not run the simulation grid against a hash that hasn't passed these checks."
  )
}
cat("[03] All hash_codes() tests passed.\n\n")

# 7. STAGE 2:.groups replicate-level simulation function and cell-level wrapper

simulate_one_replicate <- function(L, k, p, s_star) {
  N <- L + k - 1
  base_seq <- generate_iid_sequence(N)
  mut_seq <- mutate_sequence(base_seq, p)

  hA <- hash_codes(kmer_codes(base_seq, k))
  hB <- hash_codes(kmer_codes(mut_seq, k))
  uhA <- unique(hA)
  uhB <- unique(hB)

  J_exact <- length(intersect(uhA, uhB)) / length(union(uhA, uhB))

  skA <- build_sketch_from_hashes(uhA, s_star)$sketch
  skB <- build_sketch_from_hashes(uhB, s_star)$sketch

  cmp <- mash_compare(skA, skB, s_star)

  dplyr::tibble(x = cmp$x, s_used = cmp$s_used, J_exact = J_exact)
}

simulate_real_cell <- function(cell_id, L, f_min, k, J, s_star,
                               n_reps, alpha = 0.05, rho_star = 0.10,
                               verbose = TRUE) {
  p <- 1 - f_min
  t0 <- Sys.time()

  reps <- future_map_dfr(
    seq_len(n_reps),
    ~ simulate_one_replicate(L, k, p, s_star),
    .options = furrr_options(seed = TRUE),
    .progress = verbose
  )

  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  if (verbose) {
    cat(sprintf(
      "[03]   [cell %s] %d replicates in %.1fs (%.4fs/rep)\n",
      cell_id, n_reps, elapsed, elapsed / n_reps
    ))
  }

  s_used_mode <- as.integer(names(sort(table(reps$s_used), decreasing = TRUE))[1])
  if (s_used_mode != s_star) {
    warning(sprintf(
      "[03] Cell %s: modal achieved sketch size (%d) != s_star (%d)",
      cell_id, s_used_mode, s_star
    ))
  }

  diag_vs_theory <- cp_diagnostics(reps$x, s_used_mode, J_true = J, alpha = alpha)

  covered_vs_exact <- logical(n_reps)
  rho_vs_exact <- rep(NA_real_, n_reps)
  for (r in seq_len(n_reps)) {
    d <- cp_diagnostics(reps$x[r], s_used_mode, J_true = reps$J_exact[r], alpha = alpha)
    covered_vs_exact[r] <- d$covered
    rho_vs_exact[r] <- d$rho
  }

  x_ideal <- rbinom(n_reps, size = s_star, prob = J)
  diag_ideal <- cp_diagnostics(x_ideal, s_star, J_true = J, alpha = alpha)

  dplyr::tibble(
    cell_id = cell_id, L = L, f_min = f_min, k = k, J = J,
    s_star = s_star, s_used_real = s_used_mode, n_reps = n_reps,
    elapsed_sec = elapsed,
    mean_J_exact = mean(reps$J_exact),
    J_formula_bias_pct = 100 * (mean(reps$J_exact) / J - 1),
    coverage_vs_theory = mean(diag_vs_theory$covered),
    coverage_vs_exact = mean(covered_vs_exact),
    coverage_ideal = mean(diag_ideal$covered),
    median_rho_vs_theory = median(diag_vs_theory$rho, na.rm = TRUE),
    median_rho_vs_exact = median(rho_vs_exact, na.rm = TRUE),
    median_rho_ideal = median(diag_ideal$rho, na.rm = TRUE)
  )
}

# 8. Parallel backend

n_workers <- max(1, parallel::detectCores() - 2)
plan(multisession, workers = n_workers)
cat(sprintf("[03] Using %d parallel workers\n\n", n_workers))

# 9. Grid definition

L_range_real <- c(1e4, 1e5, 1e6)
f_values <- c(0.999, 0.99, 0.95, 0.90, 0.80)
n_reps_default <- 10000

grid_real <- expand_grid(L = L_range_real, f_min = f_values) |>
  mutate(
    k = k_fofanov(L),
    J = J_from_f(f_min, k),
    s_star = s_star_wald(J, rho_star = 0.10, alpha = 0.05),
    cell_id = row_number()
  ) |>
  filter(J > 0, s_star >= 2)

# 10. STAGE 3: runtime probe. Always runs, always printed, before any
#     commitment to the full grid. Uses the LARGEST L in the grid (the
#     dominant cost driver) so the estimate is not optimistic.

cat("[03] === Stage 3: runtime probe (largest-L cell, 20 replicates) ===\n")

probe_row <- grid_real |>
  filter(L == max(L)) |>
  slice(1)
t0 <- Sys.time()
probe_result <- simulate_real_cell(
  probe_row$cell_id, probe_row$L, probe_row$f_min, probe_row$k,
  probe_row$J, probe_row$s_star,
  n_reps = 20, verbose = TRUE
)
probe_elapsed <- as.numeric(Sys.time() - t0, units = "secs")

est_full_grid_min <- (probe_elapsed / 20) * n_reps_default * nrow(grid_real) / 60

cat(sprintf(
  "[03] Probe: 20 replicates at L=%.0f took %.1fs (%.3fs/rep)\n",
  probe_row$L, probe_elapsed, probe_elapsed / 20
))
cat(sprintf(
  "[03] Estimated FULL grid time (%d cells x %d reps): %.1f minutes\n\n",
  nrow(grid_real), n_reps_default, est_full_grid_min
))

write_csv(probe_result, file.path(RESULTS_DIR, "analysis4_runtime_probe.csv"))
cat("[03] Probe result saved to results/analysis4_runtime_probe.csv\n\n")

# 11. STAGE 4: full grid, with per-cell checkpointing to disk

run_full_grid <- toupper(Sys.getenv("RUN_FULL_GRID", "TRUE")) != "FALSE"

if (run_full_grid) {
  cat(sprintf(
    "[03] === Stage 4: running full grid: %d cells x %d replicates ===\n",
    nrow(grid_real), n_reps_default
  ))

  partial_path <- file.path(RESULTS_DIR, "analysis4_real_vs_idealized_partial.csv")
  results_list <- vector("list", nrow(grid_real))

  for (i in seq_len(nrow(grid_real))) {
    row <- grid_real[i, ]
    cat(sprintf(
      "[03] Cell %d/%d: L=%.0f, f_min=%.3f, k=%d, s*=%d\n",
      i, nrow(grid_real), row$L, row$f_min, row$k, row$s_star
    ))

    results_list[[i]] <- simulate_real_cell(
      row$cell_id, row$L, row$f_min, row$k, row$J, row$s_star,
      n_reps = n_reps_default
    )

    # per-cell checkpoint, internal resumability aid, distinct from the
    # top-level orchestrator's single skip-check on the final rds output
    dplyr::bind_rows(results_list) |> write_csv(partial_path)
  }

  results_real <- dplyr::bind_rows(results_list)
  write_csv(results_real, file.path(RESULTS_DIR, "analysis4_real_vs_idealized.csv"))

  # 12. Reporting

  cat("[03] === Three-way coverage comparison ===\n")
  tbl_main <- results_real |>
    dplyr::select(
      cell_id, L, f_min, s_star, mean_J_exact, J_formula_bias_pct,
      coverage_vs_theory, coverage_vs_exact, coverage_ideal,
      median_rho_vs_theory, median_rho_vs_exact, median_rho_ideal
    )
  print_full(tbl_main)

  sketching_flags <- results_real |>
    dplyr::mutate(coverage_gap_vs_ideal = coverage_vs_exact - coverage_ideal) |>
    dplyr::filter(abs(coverage_gap_vs_ideal) > 0.02 |
      abs(median_rho_vs_exact / median_rho_ideal - 1) > 0.05)

  cat(sprintf(
    "[03] %d of %d cells show >2pp coverage gap or >5%% median-rho gap between real sketching (vs exact J) and idealized binomial.\n",
    nrow(sketching_flags), nrow(results_real)
  ))
  if (nrow(sketching_flags) > 0) {
    print_full(
      sketching_flags |>
        dplyr::select(
          cell_id, L, f_min, s_star, coverage_vs_exact, coverage_ideal,
          median_rho_vs_exact, median_rho_ideal, coverage_gap_vs_ideal
        )
    )
  }

  cat("[03] === Mash/Ondov formula bias (mean J_exact vs J_theory), by f_min ===\n")
  tbl_bias <- results_real |>
    dplyr::group_by(f_min) |>
    dplyr::summarise(mean_bias_pct = mean(J_formula_bias_pct), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(f_min))
  print_full(tbl_bias)

  # 13. Checkpoint

  saveRDS(
    list(
      results_real = results_real,
      probe_result = probe_result,
      sketching_flags = sketching_flags,
      tbl_bias = tbl_bias
    ),
    file.path(RDS_DIR, "04_realistic_sketch_calib.rds")
  )

  cat("\n[03] === Analysis 3 complete ===\n")
  cat("[03] Outputs written:\n")
  cat("[03]   - results/analysis4_runtime_probe.csv\n")
  cat("[03]   - results/analysis4_real_vs_idealized_partial.csv (checkpoint, safe to delete)\n")
  cat("[03]   - results/analysis4_real_vs_idealized.csv\n")
  cat("[03]   - rds/04_realistic_sketch_calib.rds\n")
} else {
  cat("\n[03] === RUN_FULL_GRID=FALSE: probe only, full grid NOT run ===\n")
  cat("[03] No rds/04_realistic_sketch_calib.rds written -- orchestrator will re-run this step next time.\n")
}
