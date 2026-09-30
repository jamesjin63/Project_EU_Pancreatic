################################################################################
# EU Pancreatic Cancer — Value of Lost Welfare (VLW) Analysis
# GBD 2023 | EU-28 Countries | 1990–2023
#
# Method  : VSL income-elasticity transfer → VSLY → VLW
# Disease : Pancreatic cancer
# Outputs (results_VLW/pancreatic_EU_IE1/):
#   CSVs  : T1–T4
#   PDFs  : F1–F5 (9 PDF files)
################################################################################

library(tidyverse)
library(sf)
library(patchwork)
library(scales)

# ===========================================================================
# PARAMETERS
# ===========================================================================
VSL_peak_USA  <- 13.2e6      # USD 2023 PPP — US EPA (2023)
GDP_pc_USA    <- 82304.62    # USD PPP 2023 — World Bank WDI
discount_rate <- 0.03        # 3 % annual discount rate

# IE overrideable: Rscript run_pancreatic_EU_VLW.R 0.5
args   <- commandArgs(trailingOnly = TRUE)
IE     <- if (length(args) >= 1) as.numeric(args[1]) else 1.0
ie_tag <- if (IE == 1.0) "IE1" else paste0("IE", IE)

# ── Paths ──────────────────────────────────────────────────────────────────
daly_file  <- "data/merged.csv"
hale_file  <- "data/HALE.csv"
gdp_file   <- "data/gdp.csv"
world_file <- "data/df_world2.geojson"

out_dir <- file.path("results_VLW", paste0("pancreatic_EU_", ie_tag))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# COUNTRY / REGION DEFINITIONS
# ===========================================================================

# EU-28 location IDs (GBD 2023)
eu28_ids <- c(
  75L,  # Austria
  76L,  # Belgium
  45L,  # Bulgaria
  46L,  # Croatia
  77L,  # Cyprus
  47L,  # Czechia
  78L,  # Denmark
  58L,  # Estonia
  79L,  # Finland
  80L,  # France
  81L,  # Germany
  82L,  # Greece
  48L,  # Hungary
  84L,  # Ireland
  86L,  # Italy
  59L,  # Latvia
  60L,  # Lithuania
  87L,  # Luxembourg
  88L,  # Malta
  89L,  # Netherlands
  51L,  # Poland
  91L,  # Portugal
  52L,  # Romania
  54L,  # Slovakia
  55L,  # Slovenia
  92L,  # Spain
  93L,  # Sweden
  95L   # United Kingdom
)

# GBD subregion assignment (follows GBD 2023 European subregions)
subregion_map <- c(
  "75"  = "Western Europe",  # Austria
  "76"  = "Western Europe",  # Belgium
  "45"  = "Eastern Europe",  # Bulgaria
  "46"  = "Central Europe",  # Croatia
  "77"  = "Western Europe",  # Cyprus
  "47"  = "Central Europe",  # Czechia
  "78"  = "Western Europe",  # Denmark
  "58"  = "Eastern Europe",  # Estonia
  "79"  = "Western Europe",  # Finland
  "80"  = "Western Europe",  # France
  "81"  = "Western Europe",  # Germany
  "82"  = "Western Europe",  # Greece
  "48"  = "Central Europe",  # Hungary
  "84"  = "Western Europe",  # Ireland
  "86"  = "Western Europe",  # Italy
  "59"  = "Eastern Europe",  # Latvia
  "60"  = "Eastern Europe",  # Lithuania
  "87"  = "Western Europe",  # Luxembourg
  "88"  = "Western Europe",  # Malta
  "89"  = "Western Europe",  # Netherlands
  "51"  = "Central Europe",  # Poland
  "91"  = "Western Europe",  # Portugal
  "52"  = "Central Europe",  # Romania
  "54"  = "Central Europe",  # Slovakia
  "55"  = "Central Europe",  # Slovenia
  "92"  = "Western Europe",  # Spain
  "93"  = "Western Europe",  # Sweden
  "95"  = "Western Europe"   # United Kingdom
)

subregion_levels <- c("Western Europe", "Central Europe", "Eastern Europe")
subregion_colors <- c(
  "Western Europe" = "#1B4F72",
  "Central Europe" = "#117A65",
  "Eastern Europe" = "#C0392B"
)

# GBD name → World Bank name mapping (EU-specific)
gbd_to_wb_eu <- c(
  "Slovakia" = "Slovak Republic"
)

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================
norm_name <- function(x) {
  x %>%
    str_replace_all("\u2019", "'") %>%
    str_replace_all("\u2018", "'") %>%
    str_replace_all("\u00a0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish()
}

fmt_ui <- function(val, lower, upper, digits = 1, scale = 1) {
  v <- val / scale; l <- lower / scale; u <- upper / scale
  paste0(
    formatC(v, format = "f", digits = digits, big.mark = ","),
    " (", formatC(l, format = "f", digits = digits, big.mark = ","),
    "\u2013", formatC(u, format = "f", digits = digits, big.mark = ","), ")"
  )
}

calc_discount <- function(remaining_years, r = 0.03) {
  ifelse(remaining_years <= 0, 0,
         (1 - (1 / (1 + r))^remaining_years) / (r * remaining_years))
}

theme_nm <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_blank(),
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

theme_nm_map <- function(base_size = 8) {
  theme_void(base_size = base_size) +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = base_size - 1, face = "bold"),
      legend.text     = element_text(size = base_size - 1.5),
      legend.key.size = unit(0.35, "cm"),
      plot.margin     = margin(2, 2, 2, 2, "pt")
    )
}

# ===========================================================================
cat("══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  EU Pancreatic Cancer VLW Analysis — IE = %.1f\n", IE))
cat("══════════════════════════════════════════════════════════════════\n\n")

# ===========================================================================
# 1. READ DATA
# ===========================================================================
cat("[1/5] Reading data...\n")

# ── GDP (World Bank WDI, 2023) ──────────────────────────────────────────────
df_gdp_raw <- read_csv(gdp_file, show_col_types = FALSE) %>%
  filter(year == 2023, !is.na(NY.GDP.PCAP.PP.CD), !is.na(NY.GDP.MKTP.PP.CD)) %>%
  select(country,
         GDP_pc_PPP    = NY.GDP.PCAP.PP.CD,
         GDP_PPP_total = NY.GDP.MKTP.PP.CD) %>%
  mutate(country_norm = norm_name(country))

# Build EU-28 country table with WB name mapping
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
    subregion = subregion_map[as.character(location_id)],
    subregion = factor(subregion, levels = subregion_levels),
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

# ── DALYs — Pancreatic cancer, All ages, Number, 1990–2023 ─────────────────
cat("   Loading DALYs from merged.csv...\n")
df_daly <- read_csv(daly_file, show_col_types = FALSE) %>%
  filter(
    cause_name   == "Pancreatic cancer",
    measure_name == "DALYs (Disability-Adjusted Life Years)",
    metric_name  == "Number",
    age_name     == "All ages",
    location_id  %in% eu28_ids
  ) %>%
  select(location_id, location_name, sex_id, sex_name,
         year, DALY = val, lower, upper)

cat("   DALYs rows:", nrow(df_daly), "\n")

# ── EU aggregate epidemiology (loc_id=4743) for overview figures ────────────
df_eu_epi <- read_csv(daly_file, show_col_types = FALSE) %>%
  filter(
    cause_name  == "Pancreatic cancer",
    location_id == 4743L,
    age_name    == "All ages",
    sex_name    == "Both"
  ) %>%
  select(measure_name, metric_name, year, val, upper, lower)

# ── HALE — All ages, 2023 ───────────────────────────────────────────────────
df_hale <- read_csv(hale_file, show_col_types = FALSE) %>%
  filter(
    metric_name == "Years",
    year        == 2023,
    age_name    == "All ages",
    location_id %in% eu28_ids
  ) %>%
  select(location_id, sex_id, HALE = val)

cat("   HALE rows:", nrow(df_hale), " (expect 84 = 28 × 3 sexes)\n")

# ── World map geometry ──────────────────────────────────────────────────────
df_world <- st_read(world_file, quiet = TRUE)
df_world$location_id <- as.integer(round(as.numeric(df_world$location_id)))

cat("   All data loaded.\n\n")

# ===========================================================================
# 2. COMPUTE VLW
# ===========================================================================
cat(sprintf("[2/5] Computing VLW (IE = %.1f)...\n", IE))

df_vlw <- df_daly %>%
  left_join(df_gdp,  by = c("location_id", "location_name")) %>%
  left_join(df_hale, by = c("location_id", "sex_id")) %>%
  filter(!is.na(GDP_pc_PPP), !is.na(HALE)) %>%
  mutate(
    subregion = factor(subregion, levels = subregion_levels),
    # VSL transfer
    VSL_i         = VSL_peak_USA * (GDP_pc_PPP / GDP_pc_USA)^IE,
    VSLY          = VSL_i / (HALE / 2),
    # Discount factor (annuity over HALE × 40% remaining years)
    remaining_yrs = HALE * 0.4,
    disc_factor   = calc_discount(remaining_yrs, discount_rate),
    # VLW undiscounted (billion USD)
    VLW           = VSLY * DALY  / 1e9,
    VLW_lower     = VSLY * lower / 1e9,
    VLW_upper     = VSLY * upper / 1e9,
    # VLW discounted
    VLW_disc       = VLW       * disc_factor,
    VLW_disc_lower = VLW_lower * disc_factor,
    VLW_disc_upper = VLW_upper * disc_factor,
    # VLW / GDP (%)
    VLW_GDP_pct       = VLW       * 1e9 / GDP_PPP_total * 100,
    VLW_GDP_pct_lower = VLW_lower * 1e9 / GDP_PPP_total * 100,
    VLW_GDP_pct_upper = VLW_upper * 1e9 / GDP_PPP_total * 100,
    # Discounted VLW / GDP (%)
    VLW_disc_GDP_pct       = VLW_disc       * 1e9 / GDP_PPP_total * 100,
    VLW_disc_GDP_pct_lower = VLW_disc_lower * 1e9 / GDP_PPP_total * 100,
    VLW_disc_GDP_pct_upper = VLW_disc_upper * 1e9 / GDP_PPP_total * 100
  )

n_countries <- n_distinct(df_vlw$location_id)
n_years     <- n_distinct(df_vlw$year)
cat(sprintf("   VLW computed: %d countries × %d years × 3 sexes\n", n_countries, n_years))

# ── EU-28 aggregate trend (sum across countries, Both sexes) ───────────────
eu_trend <- df_vlw %>%
  filter(sex_name == "Both") %>%
  group_by(year) %>%
  summarise(
    VLW_total       = sum(VLW),       VLW_lower       = sum(VLW_lower),      VLW_upper       = sum(VLW_upper),
    VLW_disc_total  = sum(VLW_disc),  VLW_disc_lower  = sum(VLW_disc_lower), VLW_disc_upper  = sum(VLW_disc_upper),
    DALY_total      = sum(DALY),      DALY_lower      = sum(lower),          DALY_upper      = sum(upper),
    GDP_total       = sum(GDP_PPP_total),
    .groups = "drop"
  ) %>%
  mutate(
    VLW_GDP_pct       = VLW_total      * 1e9 / GDP_total * 100,
    VLW_GDP_pct_lower = VLW_lower      * 1e9 / GDP_total * 100,
    VLW_GDP_pct_upper = VLW_upper      * 1e9 / GDP_total * 100,
    VLW_disc_GDP_pct  = VLW_disc_total * 1e9 / GDP_total * 100
  )

# ── Subregion trend ────────────────────────────────────────────────────────
subregion_trend <- df_vlw %>%
  filter(sex_name == "Both") %>%
  group_by(year, subregion) %>%
  summarise(
    VLW_total      = sum(VLW),      VLW_lower      = sum(VLW_lower),      VLW_upper      = sum(VLW_upper),
    VLW_disc_total = sum(VLW_disc), VLW_disc_lower = sum(VLW_disc_lower), VLW_disc_upper = sum(VLW_disc_upper),
    DALY_total     = sum(DALY),     DALY_lower     = sum(lower),          DALY_upper     = sum(upper),
    GDP_total      = sum(GDP_PPP_total),
    .groups = "drop"
  ) %>%
  mutate(
    VLW_GDP_pct       = VLW_total * 1e9 / GDP_total * 100,
    VLW_GDP_pct_lower = VLW_lower * 1e9 / GDP_total * 100,
    VLW_GDP_pct_upper = VLW_upper * 1e9 / GDP_total * 100
  )

# 2023 slice (Both sexes)
df_2023 <- df_vlw %>% filter(year == 2023, sex_name == "Both")

cat("   Done.\n\n")

# ===========================================================================
# 3. TABLES (T1–T4)
# ===========================================================================
cat("[3/5] Saving tables...\n")

# ── T1: EU-28 aggregate annual VLW trend (1990–2023, Both sexes) ───────────
tbl1 <- eu_trend %>%
  mutate(
    `DALYs`                        = fmt_ui(DALY_total,     DALY_lower,      DALY_upper,      0,   1),
    `VLW (billion USD)`            = fmt_ui(VLW_total,      VLW_lower,       VLW_upper,       3,   1),
    `VLW discounted (billion USD)` = fmt_ui(VLW_disc_total, VLW_disc_lower,  VLW_disc_upper,  3,   1),
    `VLW/GDP (%)`                  = fmt_ui(VLW_GDP_pct,    VLW_GDP_pct_lower, VLW_GDP_pct_upper, 4, 1)
  ) %>%
  select(Year = year, `DALYs`, `VLW (billion USD)`,
         `VLW discounted (billion USD)`, `VLW/GDP (%)`,
         GDP_total_bn = GDP_total) %>%
  mutate(GDP_total_bn = round(GDP_total_bn / 1e9, 1))

write_csv(tbl1, file.path(out_dir, "T1_EU_VLW_trend.csv"))
cat("   T1 saved (EU-28 aggregate annual trend, 1990–2023).\n")

# ── T2: EU-28 country-level VLW (Both sexes, 2023) ────────────────────────
tbl2 <- df_2023 %>%
  arrange(desc(VLW_GDP_pct)) %>%
  mutate(
    VLW_rank  = row_number(),
    `DALYs`                        = fmt_ui(DALY,           lower,              upper,              0, 1),
    `VLW (billion USD)`            = fmt_ui(VLW,            VLW_lower,          VLW_upper,          3, 1),
    `VLW discounted (billion USD)` = fmt_ui(VLW_disc,       VLW_disc_lower,     VLW_disc_upper,     3, 1),
    `VLW/GDP (%)`                  = fmt_ui(VLW_GDP_pct,    VLW_GDP_pct_lower,  VLW_GDP_pct_upper,  4, 1),
    `VLW disc/GDP (%)`             = fmt_ui(VLW_disc_GDP_pct, VLW_disc_GDP_pct_lower, VLW_disc_GDP_pct_upper, 4, 1)
  ) %>%
  transmute(
    Rank         = VLW_rank,
    Country      = location_name,
    Subregion    = subregion,
    GDP_pc_PPP   = round(GDP_pc_PPP),
    GDP_total_bn = round(GDP_PPP_total / 1e9, 1),
    `DALYs`, `VLW (billion USD)`, `VLW discounted (billion USD)`,
    `VLW/GDP (%)`, `VLW disc/GDP (%)`
  )

write_csv(tbl2, file.path(out_dir, "T2_country_VLW_2023.csv"))
cat("   T2 saved (EU-28 country VLW, 2023, ranked by VLW/GDP).\n")

# ── T3: Subregion VLW trend (Both sexes, 1990–2023) ────────────────────────
tbl3 <- subregion_trend %>%
  mutate(
    `DALYs`             = fmt_ui(DALY_total,  DALY_lower,        DALY_upper,        0, 1),
    `VLW (billion USD)` = fmt_ui(VLW_total,   VLW_lower,         VLW_upper,         3, 1),
    `VLW/GDP (%)`       = fmt_ui(VLW_GDP_pct, VLW_GDP_pct_lower, VLW_GDP_pct_upper, 4, 1)
  ) %>%
  select(Year = year, Subregion = subregion, `DALYs`,
         `VLW (billion USD)`, `VLW/GDP (%)`)

write_csv(tbl3, file.path(out_dir, "T3_subregion_VLW_trend.csv"))
cat("   T3 saved (subregion VLW trend, 1990–2023).\n")

# ── T4: Country ranking summary (2023) ────────────────────────────────────
tbl4 <- df_2023 %>%
  transmute(
    Country      = location_name,
    Subregion    = subregion,
    GDP_pc_PPP   = round(GDP_pc_PPP),
    DALY         = round(DALY),
    VLW_bn       = round(VLW, 3),
    VLW_disc_bn  = round(VLW_disc, 3),
    VLW_GDP_pct  = round(VLW_GDP_pct, 4),
    VLW_disc_GDP = round(VLW_disc_GDP_pct, 4)
  ) %>%
  arrange(desc(VLW_GDP_pct)) %>%
  mutate(Rank_VLW_GDP = row_number()) %>%
  arrange(desc(VLW_bn)) %>%
  mutate(Rank_VLW_abs = row_number()) %>%
  arrange(Rank_VLW_GDP) %>%
  rename(
    `GDP pc (PPP USD)` = GDP_pc_PPP,
    `DALYs` = DALY,
    `VLW (bn USD)` = VLW_bn,
    `VLW disc (bn USD)` = VLW_disc_bn,
    `VLW/GDP (%)` = VLW_GDP_pct,
    `VLW disc/GDP (%)` = VLW_disc_GDP,
    `Rank (VLW/GDP)` = Rank_VLW_GDP,
    `Rank (VLW abs)` = Rank_VLW_abs
  )

write_csv(tbl4, file.path(out_dir, "T4_country_ranked_2023.csv"))
cat("   T4 saved (EU-28 country rankings, 2023).\n\n")

# ===========================================================================
# 4. FIGURES (F1–F5)
# ===========================================================================
cat("[4/5] Generating figures...\n")

# ── F1: EU-28 aggregate burden overview (4-panel) ──────────────────────────
# Panel a/b: Epidemiology (DALYs Number, Deaths Number from EU aggregate)
# Panel c/d: Economic burden (VLW absolute, VLW/GDP)

eu_daly_num <- df_eu_epi %>%
  filter(measure_name == "DALYs (Disability-Adjusted Life Years)", metric_name == "Number")
eu_death_num <- df_eu_epi %>%
  filter(measure_name == "Deaths", metric_name == "Number")

p_f1a <- ggplot(eu_daly_num, aes(x = year, y = val / 1e3)) +
  geom_ribbon(aes(ymin = lower / 1e3, ymax = upper / 1e3), fill = "#1B4F72", alpha = 0.18) +
  geom_line(colour = "#1B4F72", linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Year", y = "DALYs (thousands)",
       tag = "a") +
  theme_nm() +
  theme(axis.line = element_blank(), panel.border = element_blank(),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.3))

p_f1b <- ggplot(eu_death_num, aes(x = year, y = val / 1e3)) +
  geom_ribbon(aes(ymin = lower / 1e3, ymax = upper / 1e3), fill = "#C0392B", alpha = 0.18) +
  geom_line(colour = "#C0392B", linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Year", y = "Deaths (thousands)",
       tag = "b") +
  theme_nm() +
  theme(axis.line = element_blank(), panel.border = element_blank(),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.3))

p_f1c <- ggplot(eu_trend, aes(x = year, y = VLW_total)) +
  geom_ribbon(aes(ymin = VLW_lower, ymax = VLW_upper), fill = "#117A65", alpha = 0.18) +
  geom_line(colour = "#117A65", linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Year", y = "VLW (billion USD)",
       tag = "c") +
  theme_nm() +
  theme(axis.line = element_blank(), panel.border = element_blank(),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.3))

p_f1d <- ggplot(eu_trend, aes(x = year, y = VLW_GDP_pct)) +
  geom_ribbon(aes(ymin = VLW_GDP_pct_lower, ymax = VLW_GDP_pct_upper),
              fill = "#6C3483", alpha = 0.18) +
  geom_line(colour = "#6C3483", linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Year", y = "VLW / GDP (%)",
       tag = "d") +
  theme_nm() +
  theme(axis.line = element_blank(), panel.border = element_blank(),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.3))

p_f1 <- (p_f1a | p_f1b) / (p_f1c | p_f1d)

pdf(file.path(out_dir, "F1_EU_overview_4panel.pdf"), width = 7, height = 5.5)
print(p_f1)
dev.off()
cat("   F1 saved (EU overview 4-panel: DALYs / Deaths / VLW / VLW/GDP).\n")

# ── F2: 2023 EU-28 country ranking (a = VLW abs, b = VLW/GDP) ──────────────
# Panel a: countries ordered by VLW absolute value
df_rank_a <- df_2023 %>%
  arrange(VLW) %>%
  mutate(country_fct = factor(location_name, levels = location_name))

p_f2a <- ggplot(df_rank_a,
                aes(x = VLW, xmin = VLW_lower, xmax = VLW_upper,
                    y = country_fct, colour = subregion)) +
  geom_errorbarh(height = 0.3, linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = subregion_colors, name = "Subregion") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "VLW (billion USD, 2023)", y = NULL, tag = "a") +
  theme_nm(base_size = 7) +
  theme(axis.text.y   = element_text(size = 6),
        axis.line.x   = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.x  = element_line(colour = "grey30", linewidth = 0.3),
        legend.position = "top")

# Panel b: countries ordered by VLW/GDP
df_rank_b <- df_2023 %>%
  arrange(VLW_GDP_pct) %>%
  mutate(country_fct = factor(location_name, levels = location_name))

p_f2b <- ggplot(df_rank_b,
                aes(x = VLW_GDP_pct,
                    xmin = VLW_GDP_pct_lower, xmax = VLW_GDP_pct_upper,
                    y = country_fct, colour = subregion)) +
  geom_errorbarh(height = 0.3, linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = subregion_colors, name = "Subregion") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "VLW / GDP (%, 2023)", y = NULL, tag = "b") +
  theme_nm(base_size = 7) +
  theme(axis.text.y   = element_text(size = 6),
        axis.line.x   = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.x  = element_line(colour = "grey30", linewidth = 0.3),
        legend.position = "top")

p_f2 <- p_f2a | p_f2b

fig_h_f2 <- max(6, 28 * 0.22 + 1.5)
pdf(file.path(out_dir, "F2_country_ranking_ab.pdf"), width = 8, height = fig_h_f2)
print(p_f2)
dev.off()
cat("   F2 saved (country ranking: a=VLW abs, b=VLW/GDP).\n")

# ── F3: Subregion VLW trend comparison (a = VLW abs, b = VLW/GDP) ──────────
p_f3a <- ggplot(subregion_trend,
                aes(x = year, y = VLW_total, colour = subregion, fill = subregion)) +
  geom_ribbon(aes(ymin = VLW_lower, ymax = VLW_upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_colour_manual(values = subregion_colors, name = NULL) +
  scale_fill_manual(values   = subregion_colors, name = NULL) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Year", y = "VLW (billion USD)", tag = "a") +
  theme_nm() +
  theme(axis.line    = element_blank(),
        panel.border = element_blank(),
        axis.ticks   = element_line(colour = "grey30", linewidth = 0.3),
        legend.position = "top")

p_f3b <- ggplot(subregion_trend,
                aes(x = year, y = VLW_GDP_pct, colour = subregion, fill = subregion)) +
  geom_ribbon(aes(ymin = VLW_GDP_pct_lower, ymax = VLW_GDP_pct_upper),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_colour_manual(values = subregion_colors, name = NULL) +
  scale_fill_manual(values   = subregion_colors, name = NULL) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Year", y = "VLW / GDP (%)", tag = "b") +
  theme_nm() +
  theme(axis.line    = element_blank(),
        panel.border = element_blank(),
        axis.ticks   = element_line(colour = "grey30", linewidth = 0.3),
        legend.position = "top")

p_f3 <- p_f3a | p_f3b

pdf(file.path(out_dir, "F3_subregion_trend_ab.pdf"), width = 8, height = 4)
print(p_f3)
dev.off()
cat("   F3 saved (subregion VLW trend: a=abs, b=GDP%).\n")

# ── F4: Sex-stratified EU-28 aggregate trend (a = Male/Female VLW, b = VLW/GDP) ──
sex_trend <- df_vlw %>%
  filter(sex_name != "Both") %>%
  group_by(year, sex_name) %>%
  summarise(
    VLW_total       = sum(VLW),       VLW_lower       = sum(VLW_lower),      VLW_upper       = sum(VLW_upper),
    DALY_total      = sum(DALY),
    GDP_total       = sum(GDP_PPP_total),
    .groups = "drop"
  ) %>%
  mutate(
    VLW_GDP_pct       = VLW_total * 1e9 / GDP_total * 100,
    VLW_GDP_pct_lower = VLW_lower * 1e9 / GDP_total * 100,
    VLW_GDP_pct_upper = VLW_upper * 1e9 / GDP_total * 100,
    sex_name          = factor(sex_name, levels = c("Male", "Female"))
  )

sex_colors <- c("Male" = "#1B4F72", "Female" = "#C0392B")

p_f4a <- ggplot(sex_trend, aes(x = year, y = VLW_total, colour = sex_name, fill = sex_name)) +
  geom_ribbon(aes(ymin = VLW_lower, ymax = VLW_upper), alpha = 0.14, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_colour_manual(values = sex_colors, name = NULL) +
  scale_fill_manual(values   = sex_colors, name = NULL) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Year", y = "VLW (billion USD)", tag = "a") +
  theme_nm() +
  theme(axis.line = element_blank(), panel.border = element_blank(),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
        legend.position = "top")

p_f4b <- ggplot(sex_trend,
                aes(x = year, y = VLW_GDP_pct, colour = sex_name, fill = sex_name)) +
  geom_ribbon(aes(ymin = VLW_GDP_pct_lower, ymax = VLW_GDP_pct_upper),
              alpha = 0.14, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.4) +
  scale_colour_manual(values = sex_colors, name = NULL) +
  scale_fill_manual(values   = sex_colors, name = NULL) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2023)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Year", y = "VLW / GDP (%)", tag = "b") +
  theme_nm() +
  theme(axis.line = element_blank(), panel.border = element_blank(),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
        legend.position = "top")

p_f4 <- p_f4a | p_f4b

pdf(file.path(out_dir, "F4_sex_trend_ab.pdf"), width = 8, height = 4)
print(p_f4)
dev.off()
cat("   F4 saved (sex-stratified trend: a=VLW abs, b=VLW/GDP).\n")

# ── F5: Europe choropleth map (a = VLW abs, b = VLW/GDP) — left | right ────
sf::sf_use_s2(FALSE)
world_valid <- sf::st_make_valid(df_world)
world_eu    <- world_valid %>%
  filter(location_id %in% eu28_ids) %>%
  left_join(
    df_2023 %>% select(location_id, VLW, VLW_GDP_pct),
    by = "location_id"
  )

# Label points — point_on_surface is more robust for non-convex shapes
world_eu_pts <- st_point_on_surface(world_eu)
sf::sf_use_s2(TRUE)

# Short 2–3 letter labels for tiny countries to reduce clutter
label_abbrev <- c(
  "Luxembourg"     = "LUX",
  "Malta"          = "MLT",
  "Cyprus"         = "CYP",
  "Slovenia"       = "SVN",
  "Estonia"        = "EST",
  "Latvia"         = "LVA",
  "Lithuania"      = "LTU",
  "Slovakia"       = "SVK",
  "Croatia"        = "HRV",
  "Ireland"        = "IRL",
  "Denmark"        = "DNK",
  "Finland"        = "FIN",
  "Sweden"         = "SWE",
  "Netherlands"    = "NLD",
  "Belgium"        = "BEL",
  "Austria"        = "AUT",
  "Czech Republic" = "CZE",
  "Hungary"        = "HUN",
  "Romania"        = "ROU",
  "Bulgaria"       = "BGR",
  "Greece"         = "GRC",
  "Portugal"       = "PRT",
  "Spain"          = "ESP",
  "France"         = "FRA",
  "Germany"        = "DEU",
  "Italy"          = "ITA",
  "Poland"         = "POL",
  "United Kingdom" = "GBR"
)

world_eu_pts <- world_eu_pts %>%
  mutate(label = label_abbrev[name_long])

# Background (non-EU) layer
bg_map <- world_valid %>% filter(!location_id %in% eu28_ids)

p_f5a <- ggplot(world_eu) +
  geom_sf(data = bg_map, fill = "grey93", colour = "grey80", linewidth = 0.15) +
  geom_sf(aes(fill = VLW), colour = "white", linewidth = 0.3) +
  geom_sf_text(data = world_eu_pts, aes(label = label),
               size = 1.9, colour = "grey15", fontface = "bold",
               check_overlap = FALSE) +
  scale_fill_distiller(
    palette  = "Blues", direction = 1, na.value = "grey85",
    name     = "VLW\n(bn USD)",
    labels   = function(x) round(x, 1)
  ) +
  coord_sf(xlim = c(-11, 32), ylim = c(34, 71), expand = FALSE) +
  labs(tag = "a") +
  theme_nm_map()

p_f5b <- ggplot(world_eu) +
  geom_sf(data = bg_map, fill = "grey93", colour = "grey80", linewidth = 0.15) +
  geom_sf(aes(fill = VLW_GDP_pct), colour = "white", linewidth = 0.3) +
  geom_sf_text(data = world_eu_pts, aes(label = label),
               size = 1.9, colour = "grey15", fontface = "bold",
               check_overlap = FALSE) +
  scale_fill_distiller(
    palette  = "Oranges", direction = 1, na.value = "grey85",
    name     = "VLW/GDP\n(%)",
    labels   = function(x) paste0(round(x, 3), "%")
  ) +
  coord_sf(xlim = c(-11, 32), ylim = c(34, 71), expand = FALSE) +
  labs(tag = "b") +
  theme_nm_map()

# Left | right layout
p_f5 <- p_f5a | p_f5b

pdf(file.path(out_dir, "F5_EU_map_ab.pdf"), width = 12, height = 5.5)
print(p_f5)
dev.off()
cat("   F5 saved (Europe choropleth left|right: a=VLW abs, b=VLW/GDP).\n\n")

# ===========================================================================
# 5. SUMMARY
# ===========================================================================
cat("[5/5] Summary statistics...\n")

vlw_2023_total <- eu_trend %>% filter(year == 2023) %>% pull(VLW_total)
vlw_1990_total <- eu_trend %>% filter(year == 1990) %>% pull(VLW_total)
gdp_pct_2023   <- eu_trend %>% filter(year == 2023) %>% pull(VLW_GDP_pct)
gdp_pct_1990   <- eu_trend %>% filter(year == 1990) %>% pull(VLW_GDP_pct)
top_country     <- df_2023 %>% arrange(desc(VLW_GDP_pct)) %>% slice(1)

cat(sprintf("\n   ┌─────────────────────────────────────────────────────┐\n"))
cat(sprintf("   │  EU-28 Pancreatic Cancer VLW  (IE = %.1f, 2023 PPP)  │\n", IE))
cat(sprintf("   ├─────────────────────────────────────────────────────┤\n"))
cat(sprintf("   │  VLW (2023): %.2f billion USD                       │\n", vlw_2023_total))
cat(sprintf("   │  VLW (1990): %.2f billion USD                       │\n", vlw_1990_total))
cat(sprintf("   │  Change:     +%.1f%%                                  │\n",
            (vlw_2023_total - vlw_1990_total) / vlw_1990_total * 100))
cat(sprintf("   │  VLW/GDP (2023): %.4f%%                              │\n", gdp_pct_2023))
cat(sprintf("   │  VLW/GDP (1990): %.4f%%                              │\n", gdp_pct_1990))
cat(sprintf("   │  Highest VLW/GDP: %s (%.4f%%)          │\n",
            top_country$location_name, top_country$VLW_GDP_pct))
cat(sprintf("   └─────────────────────────────────────────────────────┘\n\n"))

cat(sprintf("   Output directory : %s\n", out_dir))
cat(sprintf("   Tables  : T1–T4 (4 CSV)\n"))
cat(sprintf("   Figures : F1–F5 (5 PDF)\n\n"))
cat("══ Analysis complete ══\n")
