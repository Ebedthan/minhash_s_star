#!/bin/bash
# Analysis 6: real genome benchmark

set -euo pipefail

FORCE_RERUN_06="${FORCE_RERUN_06:-FALSE}"

DATA_DIR="data"
GENOME_DIR="${DATA_DIR}/genomes"
RESULTS_DIR="results"
RDS_DIR="rds"
mkdir -p "$DATA_DIR" "$GENOME_DIR" "$RESULTS_DIR" "$RDS_DIR"

SELECTION_TSV="${DATA_DIR}/genome_selection.tsv"
ACCESSIONS_TXT="${DATA_DIR}/accessions.txt"
FASTANI_TSV="${DATA_DIR}/fastani_results.tsv"

log() { echo "[05] $1"; }

# PART A (shell): xgt-based taxonomically anchored genome selection

if [[ -f "$SELECTION_TSV" && "$FORCE_RERUN_06" != "TRUE" ]]; then
  log "Part A: $SELECTION_TSV already exists -- skipping genome selection."
else
  log "Part A: selecting genome pairs via xgt..."

  N_PER_RANK=3
  RANKS=("genus" "family" "order")
  RANK_PREFIX=(g f o)

  ANCHOR_SPECIES=(
    "s__Escherichia coli"
    "s__Staphylococcus aureus"
    "s__Pseudomonas aeruginosa"
    "s__Salmonella enterica"
    "s__Bacillus subtilis"
    "s__Klebsiella pneumoniae"
    "s__Mycobacterium tuberculosis"
    "s__Streptococcus pneumoniae"
  )

  echo -e "pair_id\tanchor_species\taccession_A\taccession_B\tshared_rank" > "$SELECTION_TSV"

  get_accessions() {
    local query="$1" flags="$2"
    xgt taxon "$query" -g $flags -O json 2>/dev/null | jq -r '.[]' 2>/dev/null || true
  }

  get_lineage() {
    xgt search "$1" -O json 2>/dev/null | jq -r '.[0].gtdbTaxonomy // empty' 2>/dev/null || true
  }

  pair_counter=0

  for anchor in "${ANCHOR_SPECIES[@]}"; do
    log "Part A: anchor $anchor"

    rep_acc=$(get_accessions "$anchor" "-r" | head -n 1)
    if [[ -z "$rep_acc" ]]; then
      echo "[05] WARNING: no representative genome for $anchor, skipping..." >&2
      continue
    fi
    log "  Representative: $rep_acc"

    lineage=$(get_lineage "$rep_acc")
    if [[ -z "$lineage" ]]; then
      echo "[05] WARNING: could not fetch lineage for $rep_acc, skipping..." >&2
      continue
    fi

    declare -A rank_name
    IFS=';' read -ra parts <<< "$lineage"
    for part in "${parts[@]}"; do
      trimmed=$(echo "$part" | sed 's/^ *//;s/ *$//')
      prefix="${trimmed:0:1}"
      name="${trimmed:3}"
      rank_name["$prefix"]="$name"
    done

    species_accs=$(get_accessions "$anchor" "")
    species_candidates=$(echo "$species_accs" | grep -v "^${rep_acc}$" || true)
    sampled=$(echo "$species_candidates" | shuf -n "$N_PER_RANK" --random-source=<(yes "$anchor") 2>/dev/null || true)
    while read -r acc; do
      [[ -z "$acc" ]] && continue
      pair_counter=$((pair_counter + 1))
      echo -e "pair_${pair_counter}\t${anchor}\t${rep_acc}\t${acc}\tspecies" >> "$SELECTION_TSV"
    done <<< "$sampled"

    prev_accs="$species_accs"
    for i in "${!RANKS[@]}"; do
      rank="${RANKS[$i]}"
      prefix="${RANK_PREFIX[$i]}"
      taxon_name="${rank_name[$prefix]:-}"

      if [[ -z "$taxon_name" ]]; then
        echo "[05] WARNING: no $rank name in lineage for $anchor, skipping this tier" >&2
        continue
      fi

      query="${prefix}__${taxon_name}"
      log "  Fetching $rank-level genomes: $query"
      rank_accs=$(get_accessions "$query" "")

      candidates=$(comm -23 <(echo "$rank_accs" | sort -u) <(echo "$prev_accs" | sort -u))

      if [[ -z "$candidates" ]]; then
        echo "[05] WARNING: no candidates at $rank level, skipping..." >&2
      else
        sampled=$(echo "$candidates" | shuf -n "$N_PER_RANK" --random-source=<(yes "${anchor}_${rank}") 2>/dev/null || true)
        while read -r acc; do
          [[ -z "$acc" ]] && continue
          pair_counter=$((pair_counter + 1))
          echo -e "pair_${pair_counter}\t${anchor}\t${rep_acc}\t${acc}\t${rank}" >> "$SELECTION_TSV"
        done <<< "$sampled"
      fi

      prev_accs="$rank_accs"
    done

    unset rank_name
  done

  n_pairs=$(tail -n +2 "$SELECTION_TSV" | wc -l)
  log "Part A complete: $n_pairs pairs written to $SELECTION_TSV"
  if [[ "$n_pairs" -eq 0 ]]; then
    echo "[05] ERROR: zero pairs produced by Part A." >&2
    exit 1
  fi
fi

# PART B (shell): NCBI `datasets` genome FASTA download

tail -n +2 "$SELECTION_TSV" | cut -f3,4 | tr '\t' '\n' | sort -u > "$ACCESSIONS_TXT"
n_total=$(wc -l < "$ACCESSIONS_TXT")
n_have=$(find "$GENOME_DIR" -maxdepth 1 -name "*.fna" 2>/dev/null | wc -l)

if [[ "$n_have" -ge "$n_total" && "$FORCE_RERUN_06" != "TRUE" ]]; then
  log "Part B: $n_have/$n_total genome FASTAs already on disk -- skipping download."
else
  log "Part B: downloading genome FASTAs ($n_total accessions)..."
  n_done=0
  while read -r acc; do
    n_done=$((n_done + 1))
    fasta_path="${GENOME_DIR}/${acc}.fna"

    if [[ -f "$fasta_path" ]]; then
      log "  [$n_done/$n_total] $acc already downloaded, skipping"
      continue
    fi

    log "  [$n_done/$n_total] Downloading $acc"
    datasets download genome accession "$acc" --include genome --filename "${DATA_DIR}/${acc}.zip"
    unzip -q -o "${DATA_DIR}/${acc}.zip" -d "${DATA_DIR}/${acc}_extracted"

    found_fasta=$(find "${DATA_DIR}/${acc}_extracted" -name "*.fna" | head -n 1)
    if [[ -z "$found_fasta" ]]; then
      echo "[05] WARNING: no .fna found for $acc -- check ${DATA_DIR}/${acc}_extracted manually" >&2
      continue
    fi
    mv "$found_fasta" "$fasta_path"
    rm -rf "${DATA_DIR}/${acc}.zip" "${DATA_DIR}/${acc}_extracted"
  done < "$ACCESSIONS_TXT"

  log "Part B complete: $(find "$GENOME_DIR" -maxdepth 1 -name "*.fna" | wc -l) / $n_total genomes on disk"
fi

# PART C (shell): FastANI ground-truth ANI computation

if [[ -f "$FASTANI_TSV" && "$FORCE_RERUN_06" != "TRUE" ]]; then
  log "Part C: $FASTANI_TSV already exists -- skipping FastANI."
else
  log "Part C: computing ground-truth ANI via FastANI..."

  echo -e "pair_id\tanchor_species\taccession_A\taccession_B\tshared_rank\tani\taligned_fraction" > "$FASTANI_TSV"

  n_lines=$(tail -n +2 "$SELECTION_TSV" | wc -l)
  i=0

  tail -n +2 "$SELECTION_TSV" | while IFS=$'\t' read -r pair_id anchor accA accB rank; do
    i=$((i + 1))
    fa="${GENOME_DIR}/${accA}.fna"
    fb="${GENOME_DIR}/${accB}.fna"

    if [[ ! -f "$fa" || ! -f "$fb" ]]; then
      echo "[05] [$i/$n_lines] SKIP $pair_id -- missing FASTA for $accA or $accB" >&2
      continue
    fi

    result=$(fastANI -q "$fa" -r "$fb" -o /dev/stdout 2>/dev/null | head -n 1)

    if [[ -z "$result" ]]; then
      echo "[05] [$i/$n_lines] $pair_id: fastANI returned no result (likely ANI too low, <~75-80%)" >&2
      ani="NA"; aligned_frac="NA"
    else
      ani=$(echo "$result" | cut -f3)
      matched=$(echo "$result" | cut -f4)
      total=$(echo "$result" | cut -f5)
      aligned_frac=$(echo "scale=4; $matched / $total" | bc)
    fi

    echo -e "${pair_id}\t${anchor}\t${accA}\t${accB}\t${rank}\t${ani}\t${aligned_frac}" >> "$FASTANI_TSV"
    log "  [$i/$n_lines] $pair_id ($accA vs $accB): ANI=$ani"
  done

  log "Part C complete. Results in $FASTANI_TSV"
  log "Note: pairs below ~75-80% ANI may return no result (known FastANI limitation)."
fi

# PARTS D + ADDENDUM + E (R): written out below via heredoc, then executed.

R_SCRIPT_PATH="06_real_genomes_analysis.R"

log "Writing embedded R analysis to $R_SCRIPT_PATH ..."
cat > "$R_SCRIPT_PATH" << 'RSCRIPT_EOF'
# Analysis 6, R portion: merged Part D (retrospective benchmark) + Part D
# addendum (bias diagnostic + taxonomic rank cross-check) + Part E
# (prospective validation at fixed a priori targets).

library(tidyverse)
library(digest)
library(Biostrings)
library(furrr)
library(future)

# Biostrings pulls in S4Vectors/IRanges/BiocGenerics, several of which mask
# dplyr verbs (rename, first, collapse, desc, slice, combine) AND base
# verbs (intersect, union, setdiff) with incompatible versions.
conflicted::conflicts_prefer(
  dplyr::rename, dplyr::filter, dplyr::first, dplyr::collapse,
  dplyr::desc, dplyr::slice, dplyr::combine,
  base::intersect, base::union, base::setdiff,
  .quiet = TRUE
)

options(tibble.print_max = Inf, tibble.width = Inf)
print_full <- function(x) print(dplyr::as_tibble(x))

DATA_DIR <- "data"
GENOME_DIR  <- file.path(DATA_DIR, "genomes")
RESULTS_DIR <- "results"
RDS_DIR <- "rds"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

# core formulas
k_fofanov  <- function(L) ceiling(log(100 * L) / log(4))
s_star_wald_real <- function(J, rho_star = 0.10, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  z^2 * (1 - J) / (rho_star^2 * J)
}
s_star_wald <- function(J, rho_star = 0.10, alpha = 0.05) {
  ceiling(s_star_wald_real(J, rho_star, alpha))
}

# fast genome loading (utf8ToInt, not strsplit)
load_genome_as_int <- function(fasta_path) {
  seqs <- readDNAStringSet(fasta_path)
  full <- paste(as.character(seqs), collapse = paste(rep("N", 50), collapse = ""))
  codes <- utf8ToInt(full)
  map <- rep(NA_integer_, 128L)
  map[utf8ToInt("A") + 1L] <- 0L
  map[utf8ToInt("C") + 1L] <- 1L
  map[utf8ToInt("G") + 1L] <- 2L
  map[utf8ToInt("T") + 1L] <- 3L
  mapped <- map[codes + 1L]
  mapped[!is.na(mapped)]
}
genome_length_only <- function(fasta_path) length(load_genome_as_int(fasta_path))

# k-mer encoding, canonicalization, hashing
kmer_codes <- function(seq_int, k) {
  N <- length(seq_int)
  n_kmers <- N - k + 1
  if (n_kmers < 1) return(integer(0))
  code <- numeric(n_kmers)
  for (j in 0:(k - 1)) {
    code <- code + seq_int[(1 + j):(n_kmers + j)] * 4^(k - 1 - j)
  }
  code
}
revcomp_int <- function(seq_int) rev(3L - seq_int)

#' Canonical k-mer codes: minimum of forward and reverse-complement code
#' at each position, making k-mer identity independent of the arbitrary
#' strand/contig orientation of independently assembled genomes.
kmer_codes_canonical <- function(seq_int, k) {
  fwd <- kmer_codes(seq_int, k)
  rc_aligned <- rev(kmer_codes(revcomp_int(seq_int), k))
  pmin(fwd, rc_aligned)
}
hash_codes <- function(codes) {
  codes_chr <- as.character(as.integer(codes))
  h_hi <- digest2int(codes_chr, seed = 1L)
  h_lo <- digest2int(codes_chr, seed = 2L)
  bitwAnd(h_hi, 0x7FFFFFFF) * 2^31 + bitwAnd(h_lo, 0x7FFFFFFF)
}
genome_hashes_at_k <- function(fasta_path, k) {
  seq_int <- load_genome_as_int(fasta_path)
  unique(hash_codes(kmer_codes_canonical(seq_int, k)))
}

# sketch comparison
mash_compare <- function(sketchA, sketchB, s_target) {
  merged <- sort(base::union(sketchA, sketchB))
  s_used <- min(s_target, length(merged))
  merged_s <- merged[seq_len(s_used)]
  x <- sum(merged_s %in% sketchA & merged_s %in% sketchB)
  list(x = x, s_used = s_used)
}
cp_diagnostics <- function(x, s, J_true, alpha = 0.05) {
  if (x <= 0 || x >= s) {
    if (x <= 0) { J_lo <- 0; J_hi <- qbeta(1 - alpha/2, 1, s) }
    else        { J_lo <- qbeta(alpha/2, s, 1); J_hi <- 1 }
    return(tibble(x = x, J_hat = x/s, J_lo = J_lo, J_hi = J_hi,
                  rho = NA_real_, degenerate = TRUE,
                  covered = (J_lo <= J_true) & (J_true <= J_hi)))
  }
  J_hat <- x / s
  J_lo <- qbeta(alpha / 2,     x,     s - x + 1)
  J_hi <- qbeta(1 - alpha / 2, x + 1, s - x)
  rho <- (J_hi - J_lo) / (2 * J_hat)
  tibble(x = x, J_hat = J_hat, J_lo = J_lo, J_hi = J_hi, rho = rho,
         degenerate = FALSE, covered = (J_lo <= J_true) & (J_true <= J_hi))
}

# PART D: retrospective real genome benchmark

alpha <- 0.05

fastani <- read_tsv(file.path(DATA_DIR, "fastani_results.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(ani)) %>%
  mutate(f_min_fastani = ani / 100) %>%
  filter(
    file.exists(file.path(GENOME_DIR, paste0(accession_A, ".fna"))),
    file.exists(file.path(GENOME_DIR, paste0(accession_B, ".fna")))
  )

cat(sprintf("[05] === %d pairs with both genomes present on disk ===\n", nrow(fastani)))

n_workers <- max(1, parallel::detectCores() - 2)
plan(multisession, workers = n_workers)
cat(sprintf("[05] Using %d parallel workers\n", n_workers))

# Pass 1: genome lengths for EVERY distinct accession appearing as either
# A or B -- needed because k depends on max(L_A, L_B) per pair.

all_accessions <- unique(c(fastani$accession_A, fastani$accession_B))
cat(sprintf("[05] Computing genome length for %d distinct accessions (A and B combined)...\n",
            length(all_accessions)))
t0 <- Sys.time()
all_lengths <- future_map_dfr(
  all_accessions,
  function(acc) {
    L <- genome_length_only(file.path(GENOME_DIR, paste0(acc, ".fna")))
    tibble::tibble(accession = acc, L_est = L)
  },
  .options = furrr_options(seed = TRUE), .progress = TRUE
)
cat(sprintf("[05] Pass 1 (lengths): %.1fs\n\n", as.numeric(Sys.time() - t0, units = "secs")))

fastani <- fastani %>%
  dplyr::left_join(all_lengths %>% dplyr::rename(accession_A = accession, L_A = L_est), by = "accession_A") %>%
  dplyr::left_join(all_lengths %>% dplyr::rename(accession_B = accession, L_B = L_est), by = "accession_B") %>%
  dplyr::mutate(
    L_pair = pmax(L_A, L_B),
    k = k_fofanov(L_pair)
  )

n_would_differ <- sum(k_fofanov(fastani$L_A) != fastani$k)
cat(sprintf("[05] Pairs where k(max(L_A,L_B)) differs from k(L_A) alone: %d / %d\n\n",
            n_would_differ, nrow(fastani)))

write_csv(fastani, file.path(RESULTS_DIR, "analysis6_pair_table_enriched.csv"))

# Build the deduplicated (accession, k) work list.
work_list <- bind_rows(
  fastani %>% transmute(accession = accession_A, k = k),
  fastani %>% transmute(accession = accession_B, k = k)
) %>%
  distinct(accession, k)

cat(sprintf("[05] === %d distinct (accession, k) hash sets needed (vs %d pairs x 2 = %d without dedup) ===\n",
            nrow(work_list), nrow(fastani), nrow(fastani) * 2))

# Pass 2: hash every distinct (accession, k) once, in parallel.
# Checkpointed: if rds/06_hash_results.rds already exists AND covers every
# (accession, k) key currently needed, re-use it instead of re-hashing,
# hashing genomes is by far the most expensive step in this analysis.

hash_checkpoint_path <- file.path(RDS_DIR, "06_hash_results.rds")
needed_keys <- paste(work_list$accession, work_list$k, sep = "|")

hash_results <- NULL
if (file.exists(hash_checkpoint_path)) {
  cached <- readRDS(hash_checkpoint_path)
  if (all(needed_keys %in% names(cached))) {
    cat("[05] Reusing cached hash_results from rds/06_hash_results.rds (all needed keys present).\n\n")
    hash_results <- cached[needed_keys]
  } else {
    cat("[05] Cached hash_results is missing some needed (accession, k) keys -- re-hashing.\n\n")
  }
}

if (is.null(hash_results)) {
  cat("[05] Hashing all distinct (accession, k) genome/k combinations...\n")
  t0 <- Sys.time()
  hash_results <- future_map(
    seq_len(nrow(work_list)),
    function(i) {
      acc <- work_list$accession[i]
      k   <- work_list$k[i]
      fa  <- file.path(GENOME_DIR, paste0(acc, ".fna"))
      genome_hashes_at_k(fa, k)
    },
    .options = furrr_options(seed = TRUE), .progress = TRUE
  )
  names(hash_results) <- needed_keys
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  saveRDS(hash_results, hash_checkpoint_path)
  cat(sprintf("[05] Pass 2 (hashing): %.1fs total, %.2fs per (accession,k)\n\n",
              elapsed, elapsed / nrow(work_list)))
}


# Pass 3: pair-level assembly from the precomputed hash-set cache
# (RETROSPECTIVE design: s* from each pair's OWN FastANI ANI)

process_pair_fast <- function(pair_id, accA, accB, f_min, k) {
  uhA <- hash_results[[paste(accA, k, sep = "|")]]
  uhB <- hash_results[[paste(accB, k, sep = "|")]]
  if (is.null(uhA) || is.null(uhB)) return(NULL)

  J_exact  <- length(base::intersect(uhA, uhB)) / length(base::union(uhA, uhB))
  J_theory <- f_min^k / (2 - f_min^k)
  s_star   <- s_star_wald(J_theory)

  run_at_s <- function(s_target, label) {
    skA <- if (length(uhA) < s_target) sort(uhA) else sort(uhA)[1:s_target]
    skB <- if (length(uhB) < s_target) sort(uhB) else sort(uhB)[1:s_target]
    cmp <- mash_compare(skA, skB, s_target)
    diag <- cp_diagnostics(cmp$x, cmp$s_used, J_true = J_exact, alpha = alpha)
    tibble(s_label = label, s_target = s_target, s_used = cmp$s_used,
           x = cmp$x, J_hat = diag$J_hat, rho = diag$rho, covered = diag$covered,
           tier = case_when(
             is.na(diag$rho) ~ "unreliable (degenerate)",
             diag$rho <= 0.10 ~ "reliable",
             diag$rho <= 0.30 ~ "borderline",
             TRUE ~ "unreliable"
           ))
  }

  bind_rows(run_at_s(1000, "default_s1000"), run_at_s(s_star, "formula_s_star")) %>%
    mutate(pair_id = pair_id, accA = accA, accB = accB,
           k = k, f_min_fastani = f_min, J_exact = J_exact,
           J_theory = J_theory, s_star = s_star)
}

t0 <- Sys.time()
results_pairs <- pmap_dfr(
  list(fastani$pair_id, fastani$accession_A, fastani$accession_B,
       fastani$f_min_fastani, fastani$k),
  process_pair_fast
)
cat(sprintf("[05] Pass 3 (pair assembly): %.1fs\n\n", as.numeric(Sys.time() - t0, units = "secs")))

write_csv(results_pairs, file.path(RESULTS_DIR, "analysis6_real_genome_comparison.csv"))

cat("[05] === Reliability tier distribution: default s=1000 vs formula s* (RETROSPECTIVE) ===\n")
tier_distribution <- results_pairs %>%
  count(s_label, tier) %>%
  pivot_wider(names_from = s_label, values_from = n, values_fill = 0)
print_full(tier_distribution)

reclass <- results_pairs %>%
  select(pair_id, f_min_fastani, s_label, tier) %>%
  pivot_wider(names_from = s_label, values_from = tier) %>%
  filter(default_s1000 != formula_s_star)

cat(sprintf("[05] %d of %d pairs (%.1f%%) reclassified between s=1000 and s*.\n",
            nrow(reclass), n_distinct(results_pairs$pair_id),
            100 * nrow(reclass) / n_distinct(results_pairs$pair_id)))

cat("\n[05] === Part D complete ===\n\n")


# PART D ADDENDUM: J_theory vs J_exact bias diagnostic + taxonomic rank
# cross-check. Reuses results_pairs directly from Part D above (same
# session); falls back to reading the CSV if run standalone.

if (!exists("results_pairs")) {
  results_pairs <- read_csv(file.path(RESULTS_DIR, "analysis6_real_genome_comparison.csv"),
                             show_col_types = FALSE)
}

pair_level <- results_pairs %>%
  distinct(pair_id, accA, accB, k, f_min_fastani, J_exact, J_theory) %>%
  mutate(
    J_bias_pct = 100 * (J_exact / J_theory - 1),
    J_bias_ratio = J_exact / J_theory
  )

cat(sprintf("[05] === J_theory vs J_exact bias: %d real genome pairs ===\n", nrow(pair_level)))

bias_overall_summary <- pair_level %>%
  summarise(
    n_pairs = n(),
    median_bias_pct = median(J_bias_pct),
    mean_bias_pct   = mean(J_bias_pct),
    min_bias_pct    = min(J_bias_pct),
    max_bias_pct    = max(J_bias_pct),
    n_negative_bias = sum(J_bias_pct < 0),
    n_over_50pct    = sum(J_bias_pct > 50)
  )
cat("[05] === Overall bias summary ===\n")
print_full(bias_overall_summary)

f_bands <- c(-Inf, 0.80, 0.85, 0.90, 0.95, 0.99, Inf)
f_labels <- c("<0.80", "0.80-0.85", "0.85-0.90", "0.90-0.95", "0.95-0.99", ">0.99")

pair_level <- pair_level %>%
  mutate(f_min_band = cut(f_min_fastani, breaks = f_bands, labels = f_labels))

cat("[05] === Bias by f_min band ===\n")
bias_by_band <- pair_level %>%
  group_by(f_min_band) %>%
  summarise(
    n_pairs = n(),
    median_bias_pct = median(J_bias_pct),
    mean_bias_pct   = mean(J_bias_pct),
    max_bias_pct    = max(J_bias_pct),
    .groups = "drop"
  )
print_full(bias_by_band)

analysis4_bias_reference <- tribble(
  ~f_min, ~analysis4_bias_pct,
  0.999,  0.0087,
  0.99,   0.0721,
  0.95,   0.519,
  0.90,   1.55,
  0.80,   8.05
)

cat("[05] === Real-genome bias vs Analysis 4's simulated-sequence bias (nearest f_min match) ===\n")
comparison_tbl <- pair_level %>%
  mutate(f_min_nearest = case_when(
    f_min_fastani >= 0.995 ~ 0.999,
    f_min_fastani >= 0.97  ~ 0.99,
    f_min_fastani >= 0.925 ~ 0.95,
    f_min_fastani >= 0.85  ~ 0.90,
    TRUE ~ 0.80
  )) %>%
  group_by(f_min_nearest) %>%
  summarise(
    n_pairs = n(),
    real_median_bias_pct = median(J_bias_pct),
    real_max_bias_pct = max(J_bias_pct),
    .groups = "drop"
  ) %>%
  left_join(analysis4_bias_reference, by = c("f_min_nearest" = "f_min")) %>%
  mutate(fold_excess_vs_simulation = real_median_bias_pct / analysis4_bias_pct)
print_full(comparison_tbl)

cat("[05] === Pairs with J_bias_pct > 100% (J_exact more than double J_theory) ===\n")
extreme_bias <- pair_level %>%
  filter(J_bias_pct > 100) %>%
  arrange(desc(J_bias_pct)) %>%
  select(pair_id, accA, accB, f_min_fastani, J_theory, J_exact, J_bias_pct)
print_full(extreme_bias)
cat(sprintf("[05] %d of %d pairs (%.1f%%) show >100%% bias.\n",
            nrow(extreme_bias), nrow(pair_level), 100 * nrow(extreme_bias) / nrow(pair_level)))

write_csv(pair_level, file.path(RESULTS_DIR, "analysis6_bias_diagnostic.csv"))

# --- taxonomic rank cross-check ---
selection_path <- file.path(DATA_DIR, "genome_selection.tsv")
rank_bias_by_rank <- NULL
rank_bias_by_anchor <- NULL
rank_bias_extreme_by_rank <- NULL
rank_bias_by_rank_anchor <- NULL

if (file.exists(selection_path)) {
  selection <- read_tsv(selection_path, show_col_types = FALSE)

  pair_level_ranked <- pair_level %>%
    left_join(selection %>% select(pair_id, shared_rank, anchor_species), by = "pair_id")

  n_unmatched <- sum(is.na(pair_level_ranked$shared_rank))
  if (n_unmatched > 0) {
    cat(sprintf("[05] WARNING: %d of %d pairs could not be matched to a shared_rank.\n",
                n_unmatched, nrow(pair_level_ranked)))
  }

  rank_order <- c("species", "genus", "family", "order")

  cat("[05] === J bias by shared taxonomic rank ===\n")
  rank_bias_by_rank <- pair_level_ranked %>%
    filter(!is.na(shared_rank)) %>%
    mutate(shared_rank = factor(shared_rank, levels = rank_order)) %>%
    group_by(shared_rank) %>%
    summarise(
      n_pairs = n(), median_f_min = median(f_min_fastani),
      median_bias_pct = median(J_bias_pct), mean_bias_pct = mean(J_bias_pct),
      max_bias_pct = max(J_bias_pct), .groups = "drop"
    ) %>%
    arrange(shared_rank)
  print_full(rank_bias_by_rank)

  cat("[05] === J bias by anchor species ===\n")
  rank_bias_by_anchor <- pair_level_ranked %>%
    filter(!is.na(shared_rank)) %>%
    group_by(anchor_species) %>%
    summarise(n_pairs = n(), median_bias_pct = median(J_bias_pct),
              max_bias_pct = max(J_bias_pct), .groups = "drop") %>%
    arrange(desc(median_bias_pct))
  print_full(rank_bias_by_anchor)

  cat("[05] === shared_rank distribution among pairs with J_bias_pct > 100% ===\n")
  rank_bias_extreme_by_rank <- pair_level_ranked %>%
    filter(!is.na(shared_rank), J_bias_pct > 100) %>%
    count(shared_rank, name = "n_extreme_bias_pairs") %>%
    arrange(desc(n_extreme_bias_pairs))
  print_full(rank_bias_extreme_by_rank)

  cat("[05] === Median bias: rank x anchor_species ===\n")
  rank_bias_by_rank_anchor <- pair_level_ranked %>%
    filter(!is.na(shared_rank)) %>%
    mutate(shared_rank = factor(shared_rank, levels = rank_order)) %>%
    group_by(anchor_species, shared_rank) %>%
    summarise(median_bias_pct = median(J_bias_pct), n = n(), .groups = "drop") %>%
    arrange(shared_rank, anchor_species)
  print_full(rank_bias_by_rank_anchor)
} else {
  cat(sprintf("[05] NOTE: %s not found -- skipping taxonomic rank cross-check.\n", selection_path))
}

cat("\n[05] === Part D addendum complete ===\n\n")

# PART E: PROSPECTIVE validation at fixed, a priori design targets
#
# s* here is computed from a FIXED target f_min_target (0.95, then 0.90),
# identical for every pair regardless of that pair's true divergence --
# mimicking a user who specifies a design target with no prior knowledge
# of any pair's actual ANI. k still varies per pair (via max(L_A, L_B)),
# since k selection is a property of genome length, not the target
# identity. Reuses hash_results and fastani directly from Part D above.

process_pair_prospective <- function(pair_id, accA, accB, k, f_min_target, alpha = 0.05) {
  uhA <- hash_results[[paste(accA, k, sep = "|")]]
  uhB <- hash_results[[paste(accB, k, sep = "|")]]
  if (is.null(uhA) || is.null(uhB)) return(NULL)

  J_exact  <- length(base::intersect(uhA, uhB)) / length(base::union(uhA, uhB))
  J_target <- f_min_target^k / (2 - f_min_target^k)   # fixed target, NOT this pair's ANI
  s_star   <- s_star_wald(J_target, rho_star = 0.10, alpha = alpha)

  skA <- if (length(uhA) < s_star) sort(uhA) else sort(uhA)[1:s_star]
  skB <- if (length(uhB) < s_star) sort(uhB) else sort(uhB)[1:s_star]
  cmp <- mash_compare(skA, skB, s_star)

  diag <- cp_diagnostics(cmp$x, cmp$s_used, J_true = J_exact, alpha = alpha)

  tibble(
    pair_id = pair_id, accA = accA, accB = accB, k = k,
    f_min_target = f_min_target, J_target = J_target, s_star = s_star,
    s_used = cmp$s_used, x = cmp$x, J_hat = diag$J_hat, rho = diag$rho,
    covered = diag$covered,
    tier = case_when(
      is.na(diag$rho) ~ "unreliable (degenerate)",
      diag$rho <= 0.10 ~ "reliable",
      diag$rho <= 0.30 ~ "borderline",
      TRUE ~ "unreliable"
    )
  )
}

stopifnot(
  "fastani is missing the k column" = "k" %in% names(fastani),
  "fastani is missing f_min_fastani" = "f_min_fastani" %in% names(fastani),
  "hash_results is empty or missing" = exists("hash_results") && length(hash_results) > 0
)

targets <- c(0.95, 0.90)

t0 <- Sys.time()
results_prospective <- map_dfr(targets, function(ft) {
  pmap_dfr(
    list(fastani$pair_id, fastani$accession_A, fastani$accession_B, fastani$k),
    ~ process_pair_prospective(..1, ..2, ..3, ..4, f_min_target = ft)
  )
})
cat(sprintf("[05] Part E pair assembly: %.1fs\n\n", as.numeric(Sys.time() - t0, units = "secs")))

results_prospective <- results_prospective %>%
  dplyr::left_join(fastani %>% dplyr::select(pair_id, f_min_fastani), by = "pair_id") %>%
  dplyr::mutate(
    ani_gap_from_target = f_min_fastani - f_min_target,
    ani_band = case_when(
      f_min_fastani >= f_min_target ~ "at or above target (easier than assumed)",
      f_min_fastani >= f_min_target - 0.05 ~ "slightly below target",
      TRUE ~ "well below target (harder than assumed)"
    )
  )

write_csv(results_prospective, file.path(RESULTS_DIR, "analysis6_prospective_validation.csv"))

cat("[05] === Prospective validation: coverage and precision by fixed target f_min ===\n")
prospective_summary <- results_prospective %>%
  group_by(f_min_target) %>%
  summarise(
    n_pairs = n(),
    s_star_range = paste0(min(s_star), "-", max(s_star)),
    coverage = mean(covered),
    median_rho = median(rho, na.rm = TRUE),
    n_reliable = sum(tier == "reliable"),
    n_borderline = sum(tier == "borderline"),
    n_unreliable = sum(tier == "unreliable"),
    n_degenerate = sum(tier == "unreliable (degenerate)"),
    .groups = "drop"
  )
print_full(prospective_summary)

cat("[05] === Coverage stratified by true-ANI relative to the fixed target ===\n")
prospective_by_gap_band <- results_prospective %>%
  group_by(f_min_target, ani_band) %>%
  summarise(n = n(), coverage = mean(covered), median_rho = median(rho, na.rm = TRUE), .groups = "drop")
print_full(prospective_by_gap_band)

cat("\n[05] === Part E complete ===\n")

# FINAL CHECKPOINT: bundle everything Analysis 6 produced for 07_figures.R

saveRDS(
  list(
    fastani = fastani,
    results_pairs = results_pairs,
    tier_distribution = tier_distribution,
    reclass = reclass,
    pair_level = pair_level,
    bias_overall_summary = bias_overall_summary,
    bias_by_band = bias_by_band,
    analysis4_bias_reference = analysis4_bias_reference,
    comparison_tbl = comparison_tbl,
    extreme_bias = extreme_bias,
    rank_bias_by_rank = rank_bias_by_rank,
    rank_bias_by_anchor = rank_bias_by_anchor,
    rank_bias_extreme_by_rank = rank_bias_extreme_by_rank,
    rank_bias_by_rank_anchor = rank_bias_by_rank_anchor,
    results_prospective = results_prospective,
    prospective_summary = prospective_summary,
    prospective_by_gap_band = prospective_by_gap_band
  ),
  file.path(RDS_DIR, "06_real_genomes.rds")
)

cat("\n[05] === Analysis 5 (R portion: Part D + addendum + Part E) complete ===\n")
cat("[05] Outputs written:\n")
cat("[05]   - results/analysis6_pair_table_enriched.csv\n")
cat("[05]   - results/analysis6_real_genome_comparison.csv\n")
cat("[05]   - results/analysis6_bias_diagnostic.csv\n")
cat("[05]   - results/analysis6_prospective_validation.csv\n")
cat("[05]   - rds/06_hash_results.rds\n")
cat("[05]   - rds/06_real_genomes.rds\n")

RSCRIPT_EOF

log "Running R analysis (Part D + addendum + Part E)..."
Rscript "$R_SCRIPT_PATH"

log "=== Analysis 6 (fully merged: Parts A-E) complete ==="
