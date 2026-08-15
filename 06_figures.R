# 07_figures.R: centralized figure generation for the entire analysis
# pipeline.


library(tidyverse)

RDS_DIR <- "rds"
FIGURES_DIR <- "figures"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

require_rds <- function(name) {
  path <- file.path(RDS_DIR, name)
  if (!file.exists(path)) {
    stop(sprintf("[06] Required checkpoint %s not found -- run the corresponding analysis step first.", path))
  }
  readRDS(path)
}

# Shared helper: consistent log-scale treatment with visible minor
# gridlines, applied throughout this script instead of plain
# scale_*_log10(labels = scales::label_number()).
log_scale_x <- function(...) {
  scale_x_log10(
    labels = scales::label_number(),
    breaks = scales::breaks_log(),
    minor_breaks = scales::minor_breaks_log(10),
    ...
  )
}
log_scale_y <- function(...) {
  scale_y_log10(
    labels = scales::label_number(),
    breaks = scales::breaks_log(),
    minor_breaks = scales::minor_breaks_log(10),
    ...
  )
}
minor_grid_theme <- theme(panel.grid.minor = element_line(linewidth = 0.15, colour = "grey85"))

cat("[06] === Generating all figures from rds/ checkpoints ===\n\n")

# Figures 01-02: Analysis 2 (Wald vs. exact interval inversion)

cat("[06] Loading rds/02_wald_approx.rds ...\n")
d02 <- require_rds("02_wald_approx.rds")
results02 <- d02$results

p01_heatmap <- results02 |>
  filter(!is.na(ratio_exact_wald)) |>
  ggplot(aes(x = J, y = factor(rho_star), fill = ratio_exact_wald)) +
  geom_tile() +
  log_scale_x() +
  scale_fill_distiller(
    palette = "RdYlBu",
    direction = -1,
    name = expression(s[exact]^"*" / s[Wald]^"*"),
    limits = c(1, NA)
  ) +
  facet_wrap(~ paste0("\u03b1 = ", alpha)) +
  labs(
    x = "True Jaccard index J (log scale)",
    y = expression(paste("Target relative half-width ", rho^"*")),
    title = expression(paste("Exact (true first crossing) / Wald ratio for minimum sketch size ", s^"*"))
  ) +
  theme_bw(base_size = 11) +
  minor_grid_theme

ggsave(file.path(FIGURES_DIR, "01_wald_vs_exact_heatmap.pdf"), p01_heatmap, width = 9, height = 4)
cat("[06]   wrote figures/01_wald_vs_exact_heatmap.pdf\n")

p02_lines <- results02 |>
  filter(!is.na(ratio_exact_wald)) |>
  ggplot(aes(x = J, y = ratio_exact_wald, colour = factor(alpha))) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  log_scale_x() +
  scale_colour_brewer(palette = "Dark2", name = expression(alpha)) +
  facet_wrap(~ paste0("\u03c1* = ", rho_star)) +
  labs(
    x = "True Jaccard index J (log scale)",
    y = expression(s[exact]^"*" / s[Wald]^"*"),
    title = "Wald approximation accuracy across the sampling range (true first-crossing minimum)"
  ) +
  theme_bw(base_size = 11) +
  minor_grid_theme

ggsave(file.path(FIGURES_DIR, "02_wald_vs_exact_lines.pdf"), p02_lines, width = 9, height = 4)
cat("[06]   wrote figures/02_wald_vs_exact_lines.pdf\n\n")

rm(d02, results02)

# Figures 03-05: Analysis 3 (idealized binomial calibration)

cat("[06] Loading rds/03_idealized_binom.rds ...\n")
d03 <- require_rds("03_idealized_binom.rds")
results03 <- d03$results
rho_sample <- d03$rho_sample_representative_cell
rep_params <- d03$representative_cell_params

p03_hist <- tibble(rho = rho_sample) |>
  filter(!is.na(rho)) |>
  ggplot(aes(x = rho)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "white") +
  geom_vline(xintercept = rep_params$rho_star, colour = "red", linetype = "dashed") +
  geom_vline(
    xintercept = median(rho_sample, na.rm = TRUE),
    colour = "blue", linetype = "dotted"
  ) +
  labs(
    x = expression(rho), y = "Count",
    title = sprintf(
      "Raw rho distribution: rho* = %.2f, alpha = %.2f, L = %.0e, f_min = %.2f",
      rep_params$rho_star, rep_params$alpha, rep_params$L, rep_params$f_min
    ),
    subtitle = "Red dashed = rho* target; blue dotted = observed median"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(FIGURES_DIR, "03_idealized_rho_histogram.pdf"), p03_hist, width = 7, height = 5)
cat("[06]   wrote figures/03_idealized_rho_histogram.pdf\n")

p04_coverage <- results03 |>
  ggplot(aes(x = f_min, y = L, fill = coverage)) +
  geom_tile() +
  log_scale_y() +
  scale_fill_gradient2(
    low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
    midpoint = 1 - 0.05,
    name = "Coverage"
  ) +
  labs(
    x = expression(f[min]), y = "Genome length L (log scale)",
    title = "Empirical CP coverage under idealized Binomial(s*, J) sampling",
    subtitle = expression(paste("Target: 1 - ", alpha, " = 0.95 (degenerate draws included)"))
  ) +
  theme_bw(base_size = 11) +
  minor_grid_theme

ggsave(file.path(FIGURES_DIR, "04_idealized_coverage_heatmap.pdf"), p04_coverage, width = 7, height = 5)
cat("[06]   wrote figures/04_idealized_coverage_heatmap.pdf\n")

p05_rho <- results03 |>
  ggplot(aes(x = f_min, y = L, fill = median_rho)) +
  geom_tile() +
  log_scale_y() +
  scale_fill_distiller(
    palette = "RdYlBu", direction = -1, name = expression(paste("Median ", rho))
  ) +
  labs(
    x = expression(f[min]), y = "Genome length L (log scale)",
    title = expression(paste("Realized median ", rho, " at formula-recommended s*")),
    subtitle = expression(paste("Target: ", rho, "* = 0.10"))
  ) +
  theme_bw(base_size = 11) +
  minor_grid_theme

ggsave(file.path(FIGURES_DIR, "05_idealized_rho_heatmap.pdf"), p05_rho, width = 7, height = 5)
cat("[06]   wrote figures/05_idealized_rho_heatmap.pdf\n\n")

rm(d03, results03, rho_sample, rep_params)


# Figures 06-07: Analysis 4 (realistic-sequence calibration)

cat("[06] Loading rds/04_realistic_sketch_calib.rds ...\n")
d04 <- require_rds("04_realistic_sketch_calib.rds")
results_real <- d04$results_real

p06_compare <- results_real |>
  dplyr::select(cell_id, L, f_min, coverage_vs_theory, coverage_vs_exact, coverage_ideal) |>
  tidyr::pivot_longer(c(coverage_vs_theory, coverage_vs_exact, coverage_ideal),
    names_to = "source", values_to = "coverage"
  ) |>
  dplyr::mutate(source = dplyr::recode(source,
    coverage_vs_theory = "Real sketch vs formula J",
    coverage_vs_exact = "Real sketch vs exact J",
    coverage_ideal = "Idealized binomial"
  )) |>
  ggplot(aes(x = factor(L), y = coverage, colour = source, shape = source)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  facet_wrap(~f_min, labeller = label_both) +
  labs(
    x = "Genome length L", y = "Coverage",
    title = "Real MinHash sketches (MurmurHash3-based): coverage vs formula J, exact J, idealized binomial",
    colour = NULL, shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "06_realistic_coverage_comparison.pdf"), p06_compare, width = 9, height = 6)
cat("[06]   wrote figures/06_realistic_coverage_comparison.pdf\n")

p07_bias <- results_real |>
  ggplot(aes(x = f_min, y = J_formula_bias_pct, colour = factor(L))) +
  geom_point(size = 2) +
  geom_line(aes(group = L)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(
    x = expression(f[min]), y = "Formula bias (%)",
    title = "Mash/Ondov J-formula bias vs exact realized J", colour = "L"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(FIGURES_DIR, "07_realistic_formula_bias.pdf"), p07_bias, width = 7, height = 5)
cat("[06]   wrote figures/07_realistic_formula_bias.pdf\n\n")

rm(d04, results_real)

# Figures 08-09: Analysis 5 (staircase structure)

cat("[06] Loading rds/05_staircase_structure.rds ...\n")
d05 <- require_rds("05_staircase_structure.rds")
grid05 <- d05$grid
steps05 <- d05$steps_by_fmin
L_max <- d05$L_max

edge_L_values <- steps05 |>
  distinct(L_predicted) |>
  pull(L_predicted)

p08_staircase <- grid05 |>
  mutate(L_Mb = L / 1e6) |>
  ggplot(aes(x = L_Mb, y = s_star, colour = factor(f_min))) +
  geom_vline(
    xintercept = edge_L_values / 1e6, linetype = "dotted",
    colour = "grey70", linewidth = 0.3
  ) +
  geom_point(size = 0.6) +
  geom_hline(yintercept = 1000, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 10000, linetype = "dashed", colour = "grey40") +
  annotate("text",
    x = L_max / 1e6 * 0.4, y = 1300, label = "s = 1000 (Mash default)",
    size = 3, colour = "grey30"
  ) +
  annotate("text",
    x = L_max / 1e6 * 0.4, y = 13000, label = "s = 10,000",
    size = 3, colour = "grey30"
  ) +
  log_scale_x() +
  log_scale_y() +
  scale_colour_brewer(palette = "RdYlBu", name = "Min. identity") +
  labs(
    x = "Genome length (Mb, log scale)",
    y = expression(paste("Minimum sketch size ", s^"*", " (log scale)"))
  ) +
  theme_bw(base_size = 11) +
  minor_grid_theme

ggsave(file.path(FIGURES_DIR, "08_staircase_full.pdf"), p08_staircase, width = 8, height = 5.5)
cat("[06]   wrote figures/08_staircase_full.pdf\n")

p09_zoom <- grid05 |>
  filter(L >= 3e5, L <= 3e6) |>
  mutate(L_Mb = L / 1e6) |>
  ggplot(aes(x = L_Mb, y = s_star, colour = factor(f_min))) +
  geom_vline(
    xintercept = edge_L_values[edge_L_values >= 3e5 & edge_L_values <= 3e6] / 1e6,
    linetype = "dotted", colour = "grey60"
  ) +
  geom_point(size = 1.2) +
  log_scale_x() +
  log_scale_y() +
  scale_colour_brewer(palette = "RdYlBu", name = "Min. identity") +
  labs(
    x = "Genome length (Mb, log scale)", y = expression(s^"*"),
    title = "Zoom near L = 1 Mb: a single step transition"
  ) +
  theme_bw(base_size = 11) +
  minor_grid_theme

ggsave(file.path(FIGURES_DIR, "09_staircase_zoom.pdf"), p09_zoom, width = 7, height = 5)
cat("[06]   wrote figures/09_staircase_zoom.pdf\n\n")

rm(d05, grid05, steps05, edge_L_values, L_max)

# Figure 10: Analysis 6 addendum (real-genome bias vs simulated-sequence bias)


cat("[06] Loading rds/06_real_genomes.rds ...\n")
d06 <- require_rds("06_real_genomes.rds")
pair_level <- d06$pair_level
analysis4_bias_reference <- d06$analysis4_bias_reference

p10_bias_compare <- ggplot() +
  geom_point(
    data = pair_level, aes(x = f_min_fastani, y = J_bias_pct),
    colour = "#D73027", alpha = 0.6, size = 2
  ) +
  geom_line(
    data = analysis4_bias_reference, aes(x = f_min, y = analysis4_bias_pct),
    colour = "#4575B4", linewidth = 0.8
  ) +
  geom_point(
    data = analysis4_bias_reference, aes(x = f_min, y = analysis4_bias_pct),
    colour = "#4575B4", size = 2.5, shape = 17
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(-10, -1, 0, 1, 10, 100),
    labels = scales::label_number()
  ) +
  labs(
    x = "True ANI / f_min",
    y = "J bias (%, signed pseudo-log scale)"
  ) +
  theme_bw(base_size = 11) +
  minor_grid_theme

ggsave(file.path(FIGURES_DIR, "10_real_genome_bias_vs_simulated.pdf"), p10_bias_compare, width = 8, height = 5.5)
cat("[06]   wrote figures/10_real_genome_bias_vs_simulated.pdf\n\n")

rm(d06, pair_level, analysis4_bias_reference)

cat("[06] === All figures generated ===\n")
cat(sprintf(
  "[06] %d PDF files written to %s/\n",
  length(list.files(FIGURES_DIR, pattern = "\\\\.pdf$")), FIGURES_DIR
))
