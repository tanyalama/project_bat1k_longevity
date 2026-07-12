#!/usr/bin/env Rscript
# ==============================================================================
# Western blot quantification — mouse vs. Myotis bat
#
# Densitometry of mitophagy / ubiquitin-pathway proteins in primary fibroblasts
# from 3 mouse and 9 Myotis individuals. Each protein's band intensity is
# normalized to a loading control (alpha-tubulin or GAPDH, as noted per target).
# For each protein we draw a boxplot (mouse vs. bat) and run an unpaired
# Welch two-sample t-test on species.
#
# Outputs:
#   western_blots_panel.png        - multi-protein boxplot panel
#   western_blots_ttests.csv       - tidy table of per-protein t-test results
#
# Reproduces the manuscript western-blot statistics.
# ==============================================================================

library(tidyverse)
library(viridis)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. Sample metadata
# ------------------------------------------------------------------------------
# 3 mouse (M*) + 9 Myotis bat (MM*) individuals, in fixed lane order.
individuals <- c(
  "M1A", "M4D", "M5E",
  "MM6309A", "MM6309B", "MM6309C",
  "MM6313A", "MM6313B", "MM6313C",
  "MMY6678A", "MMY6678B", "MMY6678E"
)
# Bat IDs start with "MM"; the three mouse IDs start with a single "M".
species_of <- ifelse(str_starts(individuals, "MM"), "Bat", "Mouse")

# ------------------------------------------------------------------------------
# 2. Densitometry values (normalized band intensities), one vector per protein
# ------------------------------------------------------------------------------
# Group A: normalized per alpha-tubulin (or GAPDH), low dynamic range (0-1).
protein_values <- list(
  NDP52    = c(0, 0, 0,
               0.2810266698, 0.367919397, 0.38742969,
               0.4710629646, 0.3551021042, 0.6406251231,
               0.4790105313, 0.2809685024, 0.5876667772),
  LC3B     = c(0.006984163017, 0, 0.006671043193,
               0.1996997688, 0.2549926058, 0.1348086993,
               0.378599972, 0.2083777865, 0.2342872004,
               0.3295423931, 0.4388619656, 0.394206074),
  P62      = c(0.04140616062, 0.2702235794, 0.1428580547,
               0.2019369078, 0.1671138313, 0.2003487468,
               0.3754961026, 0.04235991797, 0.02834543071,
               0.1689547306, 0.4395205426, 0.3839037965),
  OPTN     = c(0.09827979359, 0.3114585526, 0.1451053267,
               0.1125326842, 0.2387565038, 0.2068318554,
               0.3116408486, 0.05821120887, 0.01310877658,
               0.2131217949, 0.270703505, 0.3745616096),
  PRKN     = c(0.01173437803, 0.02583821384, 0.0191009673,
               0.02179153741, 0.05676406801, 0.04271370912,
               0.1025403754, 0.04816504654, 0.1477388017,
               0.07798231563, 0.03991493781, 0.04329136292),
  PINK1    = c(0.5835102167, 0.2724059173, 0.2820678357,
               0.161527961, 0.08563573417, 0.05539745294,
               0.1919096914, 0.07671585338, 0,
               0.1950472365, 0.5571220082, 0.1676661272),
  PINK1_v2 = c(1.179759451, 0.4793777719, 0.5576390259,
               0.3716869987, 0.1582324624, 0.1250568114,
               0.3241460607, 0.01140496819, 0,
               0.334582339, 0.9213947605, 0.03769879527),
  BNIP3    = c(0.3317771066, 0.7015505402, 0.2300220222,
               0.4640523983, 0.6433225826, 0.3201925797,
               0.2554232671, 0.5048260952, 1.95013516,
               0.3468700065, 0.4295199564, 0.2998846455)
)

# Group B: ubiquitin smears, high dynamic range (0-20).
big_protein_values <- list(
  `Phospho-Ub (per tubulin)`  = c(8.339775337, 11.6431863, 11.40283576,
                                  5.989895478, 6.520194058, 9.82912605,
                                  5.113486211, 6.322033978, 19.70682628,
                                  6.223376666, 6.792291249, 4.936742529),
  `K48-linked Ub (per GAPDH)` = c(8.647616475, 9.77801167, 9.714880671,
                                  7.850559465, 9.87239565, 8.033906558,
                                  8.898997844, 9.662638215, 9.565385772,
                                  10.9261841, 7.239464451, 13.68624336),
  `K63-linked Ub (per GAPDH)` = c(4.885333272, 4.523965796, 4.49811918,
                                  6.346429522, 6.190410517, 5.886588512,
                                  6.562982778, 6.718964511, 6.689992826,
                                  7.708040999, 6.786724831, 8.354939324),
  `Total Ub (per GAPDH)`      = c(8.87292998, 7.670514311, 8.267310593,
                                  10.80780405, 10.92114776, 10.67978054,
                                  10.79098217, 11.2329235, 11.50513434,
                                  12.15899197, 8.356434122, 14.00562606)
)

# ------------------------------------------------------------------------------
# 3. Assemble one long-format data frame
# ------------------------------------------------------------------------------
make_protein_df <- function(values, protein_name) {
  tibble(
    individual = individuals,
    species    = species_of,
    protein    = protein_name,
    value      = values
  )
}

low_df <- imap(protein_values, ~ make_protein_df(.x, .y)) |> bind_rows()
big_df <- imap(big_protein_values, ~ make_protein_df(.x, .y)) |> bind_rows()
all_df <- bind_rows(low_df, big_df) |>
  mutate(species = factor(species, levels = c("Mouse", "Bat")))

# ------------------------------------------------------------------------------
# 4. Boxplots
# ------------------------------------------------------------------------------
species_cols <- c(Mouse = "black", Bat = "#21908CFF")

make_boxplot <- function(df, ymax) {
  ggplot(df, aes(x = species, y = value, color = species)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.5) +
    geom_jitter(width = 0.15, size = 1.1, alpha = 0.9) +
    scale_color_manual(values = species_cols) +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(x = NULL, y = NULL, title = unique(df$protein)) +
    theme_classic(base_size = 10.5) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
}

low_panels <- protein_values |>
  names() |>
  set_names() |>
  map(~ make_boxplot(filter(low_df, protein == .x), ymax = 1.0))

big_panels <- big_protein_values |>
  names() |>
  set_names() |>
  map(~ make_boxplot(filter(big_df, protein == .x), ymax = 20.0))

panel <- wrap_plots(c(low_panels, big_panels), ncol = 4) +
  plot_annotation(title = "Western blot densitometry: mouse vs. Myotis bat")

ggsave("western_blots_panel.png", plot = panel,
       width = 11, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Unpaired (Welch) two-sample t-tests, species as the grouping variable
# ------------------------------------------------------------------------------
star <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

ttest_one <- function(df) {
  tt <- t.test(value ~ species, data = df)   # Welch, unpaired
  means <- df |> group_by(species) |> summarise(m = mean(value), .groups = "drop")
  m_mouse <- means$m[means$species == "Mouse"]
  m_bat   <- means$m[means$species == "Bat"]
  tibble(
    protein     = unique(df$protein),
    mean_mouse  = m_mouse,
    mean_bat    = m_bat,
    t           = unname(tt$statistic),
    df          = unname(tt$parameter),
    p_value     = tt$p.value,
    signif      = star(tt$p.value),
    higher_in   = ifelse(m_bat > m_mouse, "Bat", "Mouse")
  )
}

ttests <- all_df |>
  group_split(protein) |>
  map(ttest_one) |>
  bind_rows() |>
  arrange(p_value)

write_csv(ttests, "western_blots_ttests.csv")

print(ttests, n = Inf, width = Inf)
