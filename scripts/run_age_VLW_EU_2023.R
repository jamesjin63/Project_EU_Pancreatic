################################################################################
# EU Pancreatic Cancer — Age-Stratified VLW Analysis, 2023
# GBD 2023 | EU-28 countries | Standard 5-year age groups
#
# Method:
#   VSL_i  = VSL_US × (GDP_pc_i / GDP_pc_US)^IE
#   VSLY   = VSL_i / (HALE_allages / 2)   [all-ages HALE, sex-specific]
#   VLW    = VSLY × DALY_age / 1e9        [billion USD]
#   VLW/GDP = VLW × 1e9 / GDP_total × 100 [%]
#
# NOTE: All-ages HALE used (not age-specific) for consistency with main analysis.
#       Age variation comes entirely from the age-stratified DALYs.
#
# Outputs (results_VLW/age_analysis_EU_IE1/):
#   CSVs: T_age_EU_2023.csv, T_age_subregion_2023.csv, T_age_country_2023.csv
#   PDFs: F_age1_EU_ab.pdf        — EU-28 total: VLW & VLW/GDP by age × sex
#         F_age2_DALY_ab.pdf      — EU-28 total: DALYs by age × sex
#         F_age3_subregion_VLW.pdf  — VLW by age × sex, facet by subregion (3-panel)
#         F_age4_subregion_GDP.pdf  — VLW/GDP by age × sex, facet by subregion
################################################################################

library(tidyverse)
library(patchwork)

# ===========================================================================
# PARAMETERS
# ===========================================================================
VSL_peak_USA  <- 13.2e6
GDP_pc_USA    <- 82304.62
discount_rate <- 0.03

args   <- commandArgs(trailingOnly = TRUE)
IE     <- if (length(args) >= 1) as.numeric(args[1]) else 1.0
ie_tag <- if (IE == 1.0) "IE1" else paste0("IE", IE)

daly_file <- "data/merged.csv"
hale_file <- "data/HALE.csv"
gdp_file  <- "data/gdp.csv"

out_dir <- file.path("results_VLW", sprintf("age_analysis_EU_IE%.1f", IE))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# COUNTRY / REGION DEFINITIONS
# ===========================================================================
eu28_ids <- c(
  75L,76L,45L,46L,77L,47L,78L,58L,79L,80L,81L,82L,
  48L,84L,86L,59L,60L,87L,88L,89L,51L,91L,52L,54L,55L,92L,93L,95L
)

subregion_map <- c(
  "75"="Western Europe","76"="Western Europe","45"="Eastern Europe",
  "46"="Central Europe","77"="Western Europe","47"="Central Europe",
  "78"="Western Europe","58"="Eastern Europe","79"="Western Europe",
  "80"="Western Europe","81"="Western Europe","82"="Western Europe",
  "48"="Central Europe","84"="Western Europe","86"="Western Europe",
  "59"="Eastern Europe","60"="Eastern Europe","87"="Western Europe",
  "88"="Western Europe","89"="Western Europe","51"="Central Europe",
  "91"="Western Europe","52"="Central Europe","54"="Central Europe",
  "55"="Central Europe","92"="Western Europe","93"="Western Europe",
  "95"="Western Europe"
)
subregion_levels <- c("Western Europe", "Central Europe", "Eastern Europe")
subregion_colors <- c(
  "Western Europe" = "#1B4F72",
  "Central Europe" = "#117A65",
  "Eastern Europe" = "#C0392B"
)

gbd_to_wb_eu <- c("Slovakia" = "Slovak Republic")

# Standard 5-year age groups (20 groups, <5 to 95+)
age_5yr <- c(
  "<5 years","5-9 years","10-14 years","15-19 years",
  "20-24 years","25-29 years","30-34 years","35-39 years",
  "40-44 years","45-49 years","50-54 years","55-59 years",
  "60-64 years","65-69 years","70-74 years","75-79 years",
  "80-84 years","85-89 years","90-94 years","95+ years"
)
age_labels_short <- c(
  "<5","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
  "40-44","45-49","50-54","55-59","60-64","65-69","70-74",
  "75-79","80-84","85-89","90-94","95+"
)

sex_colors <- c("Male" = "#2E86C1", "Female" = "#C0392B", "Both" = "#5D6D7E")

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================
norm_name <- function(x) {
  x %>%
    str_replace_all("\u2019", "'") %>%
    str_replace_all("\u00a0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish()
}

fmt_ui <- function(val, lower, upper, digits = 2, scale = 1) {
  v <- val / scale; l <- lower / scale; u <- upper / scale
  paste0(
    formatC(v, format = "f", digits = digits, big.mark = ","),
    " (", formatC(l, format = "f", digits = digits, big.mark = ","),
    "\u2013", formatC(u, format = "f", digits = digits, big.mark = ","), ")"
  )
}

theme_nm <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title       = element_text(size = base_size, face = "bold", colour = "grey20"),
      axis.text        = element_text(size = base_size - 1, colour = "grey30"),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1.5),
      legend.key.size  = unit(0.35, "cm"),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = base_size, face = "bold"),
      plot.margin      = margin(4, 4, 4, 4, "pt")
    )
}

cat("══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  EU Pancreatic Cancer Age-Stratified VLW — 2023 (IE = %.1f)\n", IE))
cat("══════════════════════════════════════════════════════════════════\n\n")

# ===========================================================================
# 1. READ DATA
# ===========================================================================
cat("[1/4] Reading data...\n")

norm_name <- function(x) {
  x %>%
    str_replace_all("\u2019", "'") %>%
    str_replace_all("\u00a0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish()
}

# GDP 2023
df_gdp_raw <- read_csv(gdp_file, show_col_types = FALSE) %>%
  filter(year == 2023, !is.na(NY.GDP.PCAP.PP.CD), !is.na(NY.GDP.MKTP.PP.CD)) %>%
  select(country,
         GDP_pc_PPP    = NY.GDP.PCAP.PP.CD,
         GDP_PPP_total = NY.GDP.MKTP.PP.CD) %>%
  mutate(country_norm = norm_name(country))

df_eu_meta <- tibble(
  location_id   = eu28_ids,
  location_name = c("Austria","Belgium","Bulgaria","Croatia","Cyprus",
                    "Czechia","Denmark","Estonia","Finland","France",
                    "Germany","Greece","Hungary","Ireland","Italy",
                    "Latvia","Lithuania","Luxembourg","Malta","Netherlands",
                    "Poland","Portugal","Romania","Slovakia","Slovenia",
                    "Spain","Sweden","United Kingdom")
) %>%
  mutate(
    subregion  = factor(subregion_map[as.character(location_id)], levels = subregion_levels),
    wb_country = case_when(
      location_name %in% names(gbd_to_wb_eu) ~ gbd_to_wb_eu[location_name],
      TRUE ~ location_name
    ),
    wb_norm = norm_name(wb_country)
  )

df_gdp <- df_eu_meta %>%
  left_join(df_gdp_raw, by = c("wb_norm" = "country_norm")) %>%
  filter(!is.na(GDP_pc_PPP)) %>%
  select(location_id, location_name, subregion, GDP_pc_PPP, GDP_PPP_total)

cat("   EU countries with GDP data:", nrow(df_gdp), "/ 28\n")

# All-ages HALE (2023, sex-specific)
df_hale <- read_csv(hale_file, show_col_types = FALSE) %>%
  filter(metric_name == "Years", year == 2023, age_name == "All ages",
         location_id %in% eu28_ids) %>%
  select(location_id, sex_id, HALE = val)

cat("   HALE rows:", nrow(df_hale), "(expect 84 = 28 × 3)\n")

# Age-stratified DALYs — 2023, pancreatic cancer, standard 5-year groups
cat("   Loading age-stratified DALYs from merged.csv...\n")
df_daly_age <- read_csv(daly_file, show_col_types = FALSE) %>%
  filter(
    cause_name   == "Pancreatic cancer",
    measure_name == "DALYs (Disability-Adjusted Life Years)",
    metric_name  == "Number",
    year         == 2023,
    age_name     %in% age_5yr,
    location_id  %in% eu28_ids
  ) %>%
  select(location_id, location_name, sex_id, sex_name,
         age_name, DALY = val, lower, upper)

cat("   Age DALYs rows:", nrow(df_daly_age), "\n")
cat("   Age groups found:", n_distinct(df_daly_age$age_name),
    "(expect 20)\n")
cat("   All data loaded.\n\n")

# ===========================================================================
# 2. COMPUTE AGE-SPECIFIC VLW
# ===========================================================================
cat("[2/4] Computing age-specific VLW...\n")

df_vlw_age <- df_daly_age %>%
  left_join(df_gdp,  by = c("location_id", "location_name")) %>%
  left_join(df_hale, by = c("location_id", "sex_id")) %>%
  filter(!is.na(GDP_pc_PPP), !is.na(HALE)) %>%
  mutate(
    age_name  = factor(age_name,  levels = age_5yr),
    sex_name  = factor(sex_name,  levels = c("Male", "Female", "Both")),
    subregion = factor(subregion, levels = subregion_levels),
    # VSL → VSLY (all-ages HALE, consistent with main analysis)
    VSL_i    = VSL_peak_USA * (GDP_pc_PPP / GDP_pc_USA)^IE,
    VSLY     = VSL_i / (HALE / 2),
    # VLW (billion USD)
    VLW           = VSLY * DALY  / 1e9,
    VLW_lower     = VSLY * lower / 1e9,
    VLW_upper     = VSLY * upper / 1e9,
    # VLW/GDP (%)
    VLW_GDP_pct       = VLW       * 1e9 / GDP_PPP_total * 100,
    VLW_GDP_pct_lower = VLW_lower * 1e9 / GDP_PPP_total * 100,
    VLW_GDP_pct_upper = VLW_upper * 1e9 / GDP_PPP_total * 100
  )

cat("   VLW rows:", nrow(df_vlw_age), "\n")

# ── Aggregated EU-28 total: age × sex ──────────────────────────────────────
agg_eu <- df_vlw_age %>%
  group_by(age_name, sex_name) %>%
  summarise(
    VLW_total       = sum(VLW),       VLW_lower       = sum(VLW_lower),      VLW_upper       = sum(VLW_upper),
    DALY_total      = sum(DALY),      DALY_lower      = sum(lower),          DALY_upper      = sum(upper),
    GDP_total       = sum(unique(GDP_PPP_total)),
    .groups = "drop"
  ) %>%
  mutate(
    VLW_GDP_pct       = VLW_total * 1e9 / GDP_total * 100,
    VLW_GDP_pct_lower = VLW_lower * 1e9 / GDP_total * 100,
    VLW_GDP_pct_upper = VLW_upper * 1e9 / GDP_total * 100
  )

# ── Subregion aggregate: age × sex × subregion ────────────────────────────
agg_sub <- df_vlw_age %>%
  group_by(age_name, sex_name, subregion) %>%
  summarise(
    VLW_total       = sum(VLW),       VLW_lower       = sum(VLW_lower),      VLW_upper       = sum(VLW_upper),
    DALY_total      = sum(DALY),      DALY_lower      = sum(lower),          DALY_upper      = sum(upper),
    GDP_total       = sum(unique(GDP_PPP_total)),
    .groups = "drop"
  ) %>%
  mutate(
    VLW_GDP_pct       = VLW_total * 1e9 / GDP_total * 100,
    VLW_GDP_pct_lower = VLW_lower * 1e9 / GDP_total * 100,
    VLW_GDP_pct_upper = VLW_upper * 1e9 / GDP_total * 100
  )

cat("   Done.\n\n")

# ===========================================================================
# 3. TABLES
# ===========================================================================
cat("[3/4] Saving tables...\n")

# ── T_age1: EU-28 aggregate, age × sex ────────────────────────────────────
tbl_age_eu <- agg_eu %>%
  mutate(
    `DALYs`             = fmt_ui(DALY_total,  DALY_lower,        DALY_upper,        0, 1),
    `VLW (billion USD)` = fmt_ui(VLW_total,   VLW_lower,         VLW_upper,         3, 1),
    `VLW/GDP (%)`       = fmt_ui(VLW_GDP_pct, VLW_GDP_pct_lower, VLW_GDP_pct_upper, 5, 1)
  ) %>%
  arrange(sex_name, age_name) %>%
  select(`Age Group` = age_name, Sex = sex_name,
         `DALYs`, `VLW (billion USD)`, `VLW/GDP (%)`)

write_csv(tbl_age_eu, file.path(out_dir, "T_age_EU_2023.csv"))
cat("   T_age1 saved (EU-28 total, age × sex).\n")

# ── T_age2: Subregion, age × sex ──────────────────────────────────────────
tbl_age_sub <- agg_sub %>%
  mutate(
    `DALYs`             = fmt_ui(DALY_total,  DALY_lower,        DALY_upper,        0, 1),
    `VLW (billion USD)` = fmt_ui(VLW_total,   VLW_lower,         VLW_upper,         3, 1),
    `VLW/GDP (%)`       = fmt_ui(VLW_GDP_pct, VLW_GDP_pct_lower, VLW_GDP_pct_upper, 5, 1)
  ) %>%
  arrange(subregion, sex_name, age_name) %>%
  select(Subregion = subregion, `Age Group` = age_name, Sex = sex_name,
         `DALYs`, `VLW (billion USD)`, `VLW/GDP (%)`)

write_csv(tbl_age_sub, file.path(out_dir, "T_age_subregion_2023.csv"))
cat("   T_age2 saved (subregion, age × sex).\n")

# ── T_age3: Country-level, age × sex (2023, numeric) ──────────────────────
tbl_age_country <- df_vlw_age %>%
  filter(sex_name != "Both") %>%
  arrange(location_name, sex_name, age_name) %>%
  transmute(
    Country     = location_name,
    Subregion   = subregion,
    Sex         = sex_name,
    `Age Group` = age_name,
    DALY        = round(DALY, 1),
    VLW_bn      = round(VLW, 4),
    VLW_GDP_pct = round(VLW_GDP_pct, 5)
  )

write_csv(tbl_age_country, file.path(out_dir, "T_age_country_2023.csv"))
cat("   T_age3 saved (country level, age × sex).\n\n")

# ===========================================================================
# 4. FIGURES
# ===========================================================================
cat("[4/4] Generating figures...\n")

# Shared x-axis scale (short labels)
scale_x_age <- scale_x_discrete(
  labels = setNames(age_labels_short, age_5yr)
)

ax_theme <- theme(
  axis.text.x       = element_text(angle = 45, hjust = 1, size = 6),
  axis.ticks        = element_line(colour = "grey30", linewidth = 0.3),
  axis.ticks.length = unit(0.12, "cm"),
  panel.border      = element_blank(),
  axis.line         = element_blank()
)

# ── F_age1: EU-28 aggregate VLW by age — a: VLW abs, b: VLW/GDP ────────────
eu_mf <- agg_eu %>% filter(sex_name %in% c("Male","Female","Both"))

pa1 <- ggplot(eu_mf,
              aes(x = age_name, y = VLW_total, colour = sex_name,
                  group = sex_name)) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VLW_lower, ymax = VLW_upper, fill = sex_name),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = sex_colors, name = "Sex") +
  scale_fill_manual(values   = sex_colors, guide = "none") +
  scale_x_age +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x)
                       ifelse(x >= 1, paste0(round(x, 1)), formatC(x, format = "f", digits = 3))) +
  labs(x = NULL, y = "VLW (billion USD, 2023)", tag = "a") +
  theme_nm(base_size = 9) +
  theme(legend.position = "top", plot.tag = element_text(size = 10, face = "bold")) +
  ax_theme

pb1 <- ggplot(eu_mf,
              aes(x = age_name, y = VLW_GDP_pct, colour = sex_name,
                  group = sex_name)) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VLW_GDP_pct_lower, ymax = VLW_GDP_pct_upper,
                  fill = sex_name),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = sex_colors, name = "Sex") +
  scale_fill_manual(values   = sex_colors, guide = "none") +
  scale_x_age +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(formatC(x, format="f", digits=4), "%")) +
  labs(x = "Age group", y = "VLW / GDP (%, 2023)", tag = "b") +
  theme_nm(base_size = 9) +
  theme(legend.position = "none", plot.tag = element_text(size = 10, face = "bold")) +
  ax_theme

p_f_age1 <- ((pa1 / pb1) + plot_layout(heights = c(1, 1)))

pdf(file.path(out_dir, "F_age1_EU_ab.pdf"), width = 9, height = 7)
print(p_f_age1)
dev.off()
cat("   F_age1 saved (EU-28 total VLW by age: a=abs, b=VLW/GDP).\n")

# ── F_age2: EU-28 DALYs by age — a: Male/Female/Both lines, b: Male vs Female bar ──
pa2 <- ggplot(eu_mf,
              aes(x = age_name, y = DALY_total, colour = sex_name, group = sex_name)) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_ribbon(aes(ymin = DALY_lower, ymax = DALY_upper, fill = sex_name),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = sex_colors, name = "Sex") +
  scale_fill_manual(values   = sex_colors, guide = "none") +
  scale_x_age +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x)
                       ifelse(x >= 1000, paste0(round(x/1000, 0), "K"), round(x, 0))) +
  labs(x = NULL, y = "DALYs (2023)", tag = "a") +
  theme_nm(base_size = 9) +
  theme(legend.position = "top", plot.tag = element_text(size = 10, face = "bold")) +
  ax_theme

# Panel b: Male vs Female stacked bar showing DALYs composition by age
eu_mf_bar <- agg_eu %>%
  filter(sex_name %in% c("Male", "Female")) %>%
  mutate(sex_name = factor(sex_name, levels = c("Male", "Female")))

pb2 <- ggplot(eu_mf_bar,
              aes(x = age_name, y = DALY_total, fill = sex_name)) +
  geom_col(position = "stack", width = 0.75, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("Male" = "#2E86C1", "Female" = "#C0392B"), name = "Sex") +
  scale_x_age +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06)),
                     labels = function(x)
                       ifelse(x >= 1000, paste0(round(x/1000, 0), "K"), round(x, 0))) +
  labs(x = "Age group", y = "DALYs — Male + Female (2023)", tag = "b") +
  theme_nm(base_size = 9) +
  theme(legend.position = "none", plot.tag = element_text(size = 10, face = "bold")) +
  ax_theme

p_f_age2 <- ((pa2 / pb2) + plot_layout(heights = c(1, 1)))

pdf(file.path(out_dir, "F_age2_DALY_ab.pdf"), width = 9, height = 7)
print(p_f_age2)
dev.off()
cat("   F_age2 saved (EU-28 DALYs by age: a=line, b=stacked bar).\n")

# ── F_age3: VLW by age × sex — 3-panel facet by subregion ──────────────────
sub_mf <- agg_sub %>% filter(sex_name %in% c("Male", "Female"))

p_f_age3 <- ggplot(sub_mf,
                   aes(x = age_name, y = VLW_total, colour = sex_name,
                       group = sex_name)) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VLW_lower, ymax = VLW_upper, fill = sex_name),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = sex_colors[c("Male","Female")], name = "Sex") +
  scale_fill_manual(values   = sex_colors[c("Male","Female")], guide = "none") +
  scale_x_age +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x)
                       ifelse(x >= 1, round(x, 1), formatC(x, format = "f", digits = 3))) +
  facet_wrap(~ subregion, scales = "free_y", ncol = 3) +
  labs(x = "Age group", y = "VLW (billion USD, 2023)") +
  theme_nm(base_size = 8) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 5.5),
    axis.ticks        = element_line(colour = "grey30", linewidth = 0.3),
    axis.ticks.length = unit(0.12, "cm"),
    legend.position   = "top",
    panel.border      = element_blank(),
    axis.line         = element_blank()
  )

pdf(file.path(out_dir, "F_age3_subregion_VLW.pdf"), width = 11, height = 4.5)
print(p_f_age3)
dev.off()
cat("   F_age3 saved (VLW by age × sex, 3-panel subregion).\n")

# ── F_age4: VLW/GDP by age × sex — 3-panel facet by subregion ──────────────
p_f_age4 <- ggplot(sub_mf,
                   aes(x = age_name, y = VLW_GDP_pct, colour = sex_name,
                       group = sex_name)) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VLW_GDP_pct_lower, ymax = VLW_GDP_pct_upper,
                  fill = sex_name),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = sex_colors[c("Male","Female")], name = "Sex") +
  scale_fill_manual(values   = sex_colors[c("Male","Female")], guide = "none") +
  scale_x_age +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x)
                       paste0(formatC(x, format = "f", digits = 4), "%")) +
  facet_wrap(~ subregion, scales = "free_y", ncol = 3) +
  labs(x = "Age group", y = "VLW / GDP (%, 2023)") +
  theme_nm(base_size = 8) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 5.5),
    axis.ticks        = element_line(colour = "grey30", linewidth = 0.3),
    axis.ticks.length = unit(0.12, "cm"),
    legend.position   = "top",
    panel.border      = element_blank(),
    axis.line         = element_blank()
  )

pdf(file.path(out_dir, "F_age4_subregion_GDP.pdf"), width = 11, height = 4.5)
print(p_f_age4)
dev.off()
cat("   F_age4 saved (VLW/GDP by age × sex, 3-panel subregion).\n\n")

# ===========================================================================
# 5. SUMMARY
# ===========================================================================
cat("[5/5] Summary...\n")

peak_both <- agg_eu %>%
  filter(sex_name == "Both") %>%
  arrange(desc(VLW_total)) %>%
  slice(1)

peak_male <- agg_eu %>%
  filter(sex_name == "Male") %>%
  arrange(desc(VLW_total)) %>%
  slice(1)

peak_female <- agg_eu %>%
  filter(sex_name == "Female") %>%
  arrange(desc(VLW_total)) %>%
  slice(1)

total_both <- agg_eu %>%
  filter(sex_name == "Both") %>%
  summarise(total_vlw = sum(VLW_total)) %>%
  pull(total_vlw)

cat(sprintf("\n   Peak VLW age group (Both):   %s  (%.3f bn USD)\n",
            as.character(peak_both$age_name), peak_both$VLW_total))
cat(sprintf("   Peak VLW age group (Male):   %s  (%.3f bn USD)\n",
            as.character(peak_male$age_name), peak_male$VLW_total))
cat(sprintf("   Peak VLW age group (Female): %s  (%.3f bn USD)\n",
            as.character(peak_female$age_name), peak_female$VLW_total))
cat(sprintf("   EU-28 total VLW, all ages summed (Both): %.2f bn USD\n", total_both))
cat(sprintf("\n   Output : %s\n", out_dir))
cat(sprintf("   Tables : T_age_EU_2023.csv  T_age_subregion_2023.csv  T_age_country_2023.csv\n"))
cat(sprintf("   Figures: F_age1–F_age4 (4 PDF)\n\n"))
cat("══ Age analysis complete ══\n")
