################################################################################
# EU Pancreatic Cancer — VLW Forecast to 2050 (ETS)
# GBD 2023 | EU-28 countries | IE = 1.0 (overrideable)
#
# Forecasts:
#   • EU-28 aggregate
#   • 3 subregions (Western / Central / Eastern Europe)
#   • 28 individual EU-28 countries
#
# Method: ETS (exponential smoothing, auto AIC) on historical VLW 1990–2023
#   VSL_i  = VSL_US × (GDP_pc_i / GDP_pc_US)^IE
#   VSLY   = VSL_i / (HALE_allages / 2)
#   VLW    = VSLY × DALY / 1e9  [billion USD]
#   VLW/GDP uses 2023 GDP (fixed denominator)
#
# Outputs (results_VLW/forecast_EU_IE1/):
#   CSVs: Forecast_EU_raw.csv  Forecast_subregion_raw.csv
#         Forecast_country_raw.csv  Forecast_summary_key_years.csv
#         ETS_model_summary.csv
#   PDFs: Forecast_EU_ab.pdf          — EU total VLW & VLW/GDP
#         Forecast_subregion_ab.pdf   — 3 subregions overlay (ab)
#         Forecast_sub3panel_VLW.pdf  — 3-panel VLW by subregion
#         Forecast_sub3panel_GDP.pdf  — 3-panel VLW/GDP by subregion
#         Forecast_country_28panel.pdf — 28-panel all countries (VLW)
#         Forecast_country_top10.pdf  — top 10 countries ab (VLW & VLW/GDP)
################################################################################

library(tidyverse)
library(forecast)
library(patchwork)
library(scales)

options(warn = -1)

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

out_dir   <- file.path("results_VLW", sprintf("forecast_EU_IE%.1f", IE))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

HORIZON  <- 27L   # 2024–2050
CI_LVL   <- 95

# ===========================================================================
# COUNTRY / REGION DEFINITIONS
# ===========================================================================
eu28_ids <- c(
  75L,76L,45L,46L,77L,47L,78L,58L,79L,80L,81L,82L,
  48L,84L,86L,59L,60L,87L,88L,89L,51L,91L,52L,54L,55L,92L,93L,95L
)
eu28_names <- c(
  "Austria","Belgium","Bulgaria","Croatia","Cyprus","Czechia","Denmark",
  "Estonia","Finland","France","Germany","Greece","Hungary","Ireland",
  "Italy","Latvia","Lithuania","Luxembourg","Malta","Netherlands",
  "Poland","Portugal","Romania","Slovakia","Slovenia","Spain","Sweden",
  "United Kingdom"
)

subregion_map <- c(
  "75"="Western Europe","76"="Western Europe","45"="Eastern Europe",
  "46"="Central Europe","47"="Central Europe","77"="Western Europe",
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

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================
norm_name <- function(x) {
  x |> str_replace_all("\u2019","'") |> str_replace_all("\u00a0"," ") |>
    iconv(from="UTF-8", to="ASCII//TRANSLIT") |> str_squish()
}

calc_discount <- function(n, r = 0.03) {
  ifelse(n <= 0, 0, (1 - (1/(1+r))^n) / (r * n))
}

fmt_ui <- function(v, lo, hi, digits = 2) {
  sprintf("%s (%s\u2013%s)",
          formatC(v,  "f", digits, big.mark = ","),
          formatC(lo, "f", digits, big.mark = ","),
          formatC(hi, "f", digits, big.mark = ","))
}

# ETS wrapper with tryCatch for robustness on small series
run_ets <- function(vlw_vec, h = HORIZON, level = CI_LVL) {
  tryCatch({
    ts_obj <- ts(vlw_vec, start = 1990, frequency = 1)
    fit    <- ets(ts_obj)
    fc     <- forecast(fit, h = h, level = level)
    list(
      fc_df = tibble(
        year   = seq(2024, 2023 + h),
        VLW    = as.numeric(fc$mean),
        VLW_lo = as.numeric(fc$lower[, 1]),
        VLW_hi = as.numeric(fc$upper[, 1])
      ),
      model_str = sprintf("%s (AIC = %.1f)", fit$method, fit$aic),
      ok = TRUE
    )
  }, error = function(e) {
    list(
      fc_df = tibble(year = seq(2024, 2023 + h), VLW = NA_real_,
                     VLW_lo = NA_real_, VLW_hi = NA_real_),
      model_str = paste("ERROR:", conditionMessage(e)),
      ok = FALSE
    )
  })
}

theme_nm <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title        = element_text(size = base_size, face = "bold", colour = "grey20"),
      axis.text         = element_text(size = base_size - 1, colour = "grey30"),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 1.5),
      legend.key.size   = unit(0.35, "cm"),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold"),
      plot.margin       = margin(4, 4, 4, 4, "pt")
    )
}

SPLIT_YEAR <- 2023.5

# Panel builder — single series
make_panel <- function(hist_df, fc_df, col, y_col, ylo_col, yhi_col,
                       y_label, title_str, pct = FALSE, tag = NULL) {
  y_fmt <- if (pct) {
    function(x) paste0(formatC(x, "f", 3), "%")
  } else {
    function(x) ifelse(x >= 1000, paste0(round(x/1000, 1), "K"),
                       ifelse(x >= 1, round(x, 1), formatC(x, "f", 3)))
  }

  p <- ggplot() +
    geom_ribbon(data = hist_df,
                aes(x = year, ymin = .data[[ylo_col]], ymax = .data[[yhi_col]]),
                fill = col, alpha = 0.22, colour = NA) +
    geom_ribbon(data = fc_df,
                aes(x = year, ymin = .data[[ylo_col]], ymax = .data[[yhi_col]]),
                fill = col, alpha = 0.12, colour = NA) +
    geom_line(data = hist_df, aes(x = year, y = .data[[y_col]]),
              colour = col, linewidth = 0.75) +
    geom_line(data = fc_df, aes(x = year, y = .data[[y_col]]),
              colour = col, linewidth = 0.75, linetype = "dashed") +
    geom_vline(xintercept = SPLIT_YEAR, colour = "grey45", linewidth = 0.35,
               linetype = "dotted") +
    geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.35) +
    geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.35) +
    scale_x_continuous(limits = c(1990, 2050),
                       breaks = c(1990, 2000, 2010, 2023, 2040, 2050),
                       labels = c("1990","2000","2010","2023","2040","2050")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = y_fmt) +
    labs(x = "Year", y = y_label, title = title_str) +
    theme_nm(base_size = 8) +
    theme(
      plot.title        = element_text(size = 8, face = "bold", colour = col,
                                       margin = margin(b = 2)),
      axis.text.x       = element_text(angle = 30, hjust = 1, size = 6),
      axis.ticks        = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks.length = unit(0.12, "cm"),
      panel.border      = element_blank(),
      axis.line         = element_blank()
    )
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}

cat("══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  EU Pancreatic Cancer VLW Forecast 2024–2050 (IE = %.1f)\n", IE))
cat("══════════════════════════════════════════════════════════════════\n\n")

# ===========================================================================
# 1. LOAD & COMPUTE HISTORICAL VLW
# ===========================================================================
cat("[1/5] Loading data and computing historical VLW (1990–2023)...\n")

df_gdp_raw <- read_csv(gdp_file, show_col_types = FALSE) |>
  filter(year == 2023, !is.na(NY.GDP.PCAP.PP.CD), !is.na(NY.GDP.MKTP.PP.CD)) |>
  transmute(country_norm = norm_name(country),
            GDP_pc_PPP    = NY.GDP.PCAP.PP.CD,
            GDP_PPP_total = NY.GDP.MKTP.PP.CD)

df_eu_meta <- tibble(location_id = eu28_ids, location_name = eu28_names) |>
  mutate(
    subregion  = factor(subregion_map[as.character(location_id)], levels = subregion_levels),
    wb_country = if_else(location_name %in% names(gbd_to_wb_eu),
                         gbd_to_wb_eu[location_name], location_name),
    wb_norm    = norm_name(wb_country)
  )

df_gdp <- df_eu_meta |>
  left_join(df_gdp_raw, by = c("wb_norm" = "country_norm")) |>
  filter(!is.na(GDP_pc_PPP)) |>
  select(location_id, location_name, subregion, GDP_pc_PPP, GDP_PPP_total)

df_hale <- read_csv(hale_file, show_col_types = FALSE) |>
  filter(metric_name == "Years", year == 2023, age_name == "All ages",
         location_id %in% eu28_ids, sex_name == "Both") |>
  select(location_id, HALE = val)

df_daly <- read_csv(daly_file, show_col_types = FALSE) |>
  filter(cause_name == "Pancreatic cancer",
         measure_name == "DALYs (Disability-Adjusted Life Years)",
         metric_name  == "Number",
         age_name     == "All ages",
         sex_name     == "Both",
         location_id  %in% eu28_ids) |>
  select(location_id, location_name, year, DALY = val, lower, upper)

df_vlw_hist <- df_daly |>
  left_join(df_gdp,  by = c("location_id","location_name")) |>
  left_join(df_hale, by = "location_id") |>
  filter(!is.na(GDP_pc_PPP), !is.na(HALE)) |>
  mutate(
    VSL_i     = VSL_peak_USA * (GDP_pc_PPP / GDP_pc_USA)^IE,
    VSLY      = VSL_i / (HALE / 2),
    VLW       = VSLY * DALY  / 1e9,
    VLW_lo    = VSLY * lower / 1e9,
    VLW_hi    = VSLY * upper / 1e9
  )

cat(sprintf("   Countries: %d  |  Years: %d\n",
            n_distinct(df_vlw_hist$location_id), n_distinct(df_vlw_hist$year)))

# ── EU-28 total GDP (fixed 2023 denominator) ──────────────────────────────
total_gdp <- sum(df_gdp$GDP_PPP_total)
cat(sprintf("   EU-28 total GDP (2023 PPP): %.2f trillion USD\n\n",
            total_gdp / 1e12))

# ── Build aggregated series ────────────────────────────────────────────────
# EU total
hist_eu <- df_vlw_hist |>
  group_by(year) |>
  summarise(VLW = sum(VLW), VLW_lo = sum(VLW_lo), VLW_hi = sum(VLW_hi),
            .groups = "drop") |>
  mutate(VLW_GDP    = VLW    / (total_gdp / 1e9) * 100,
         VLW_GDP_lo = VLW_lo / (total_gdp / 1e9) * 100,
         VLW_GDP_hi = VLW_hi / (total_gdp / 1e9) * 100,
         series = "EU-28 Total")

# Subregion totals
sub_gdp <- df_gdp |>
  group_by(subregion) |>
  summarise(gdp_total = sum(GDP_PPP_total), .groups = "drop")

hist_sub <- df_vlw_hist |>
  group_by(year, subregion) |>
  summarise(VLW = sum(VLW), VLW_lo = sum(VLW_lo), VLW_hi = sum(VLW_hi),
            .groups = "drop") |>
  left_join(sub_gdp, by = "subregion") |>
  mutate(VLW_GDP    = VLW    / (gdp_total / 1e9) * 100,
         VLW_GDP_lo = VLW_lo / (gdp_total / 1e9) * 100,
         VLW_GDP_hi = VLW_hi / (gdp_total / 1e9) * 100,
         series = as.character(subregion))

# Country-level — GDP_PPP_total already present from initial join
hist_country <- df_vlw_hist |>
  mutate(VLW_GDP    = VLW    / (GDP_PPP_total / 1e9) * 100,
         VLW_GDP_lo = VLW_lo / (GDP_PPP_total / 1e9) * 100,
         VLW_GDP_hi = VLW_hi / (GDP_PPP_total / 1e9) * 100,
         series = location_name)

cat("[1/5] Historical VLW series built.\n\n")

# ===========================================================================
# 2. ETS FORECASTING
# ===========================================================================
cat("[2/5] Fitting ETS models...\n")

model_log <- list()
fc_eu <- fc_sub <- fc_country <- list()

# ── EU total ──────────────────────────────────────────────────────────────
res <- run_ets(hist_eu |> arrange(year) |> pull(VLW))
fc_eu[["EU-28 Total"]] <- res$fc_df |>
  mutate(series = "EU-28 Total",
         VLW_GDP    = VLW    / (total_gdp / 1e9) * 100,
         VLW_GDP_lo = VLW_lo / (total_gdp / 1e9) * 100,
         VLW_GDP_hi = VLW_hi / (total_gdp / 1e9) * 100)
model_log[["EU-28 Total"]] <- tibble(Series = "EU-28 Total",
                                      Level  = "EU Aggregate",
                                      ETS_model = res$model_str)
cat(sprintf("   EU-28 Total:  %s\n", res$model_str))

# ── Subregions ────────────────────────────────────────────────────────────
for (s in subregion_levels) {
  gdp_s <- sub_gdp |> filter(subregion == s) |> pull(gdp_total)
  vec   <- hist_sub |> filter(subregion == s) |> arrange(year) |> pull(VLW)
  res   <- run_ets(vec)
  fc_sub[[s]] <- res$fc_df |>
    mutate(series = s, subregion = factor(s, levels = subregion_levels),
           VLW_GDP    = VLW    / (gdp_s / 1e9) * 100,
           VLW_GDP_lo = VLW_lo / (gdp_s / 1e9) * 100,
           VLW_GDP_hi = VLW_hi / (gdp_s / 1e9) * 100)
  model_log[[s]] <- tibble(Series = s, Level = "Subregion",
                            ETS_model = res$model_str)
  cat(sprintf("   %-20s  %s\n", s, res$model_str))
}

# ── 28 countries ──────────────────────────────────────────────────────────
cat("   Country-level ETS:\n")
for (cty in eu28_names) {
  gdp_c <- df_gdp |> filter(location_name == cty) |> pull(GDP_PPP_total)
  vec   <- hist_country |> filter(location_name == cty) |>
    arrange(year) |> pull(VLW)
  res   <- run_ets(vec)
  fc_country[[cty]] <- res$fc_df |>
    mutate(series = cty, location_name = cty,
           VLW_GDP    = VLW    / (gdp_c / 1e9) * 100,
           VLW_GDP_lo = VLW_lo / (gdp_c / 1e9) * 100,
           VLW_GDP_hi = VLW_hi / (gdp_c / 1e9) * 100)
  model_log[[cty]] <- tibble(Series = cty, Level = "Country",
                               ETS_model = res$model_str)
  cat(sprintf("     %-20s  %s\n", cty, res$model_str))
}
cat("\n")

# Combine model log
df_model_log <- bind_rows(model_log)
write_csv(df_model_log, file.path(out_dir, "ETS_model_summary.csv"))

# ===========================================================================
# 3. BUILD FULL SERIES & SAVE CSV
# ===========================================================================
cat("[3/5] Building full series and saving CSVs...\n")

build_full <- function(hist_df, fc_list, key_col = "series") {
  bind_rows(
    hist_df |> transmute(year, series = .data[[key_col]], VLW, VLW_lo, VLW_hi,
                         VLW_GDP, VLW_GDP_lo, VLW_GDP_hi, period = "Historical"),
    bind_rows(fc_list) |>
      transmute(year, series, VLW, VLW_lo, VLW_hi,
                VLW_GDP, VLW_GDP_lo, VLW_GDP_hi, period = "Forecast")
  ) |> mutate(period = factor(period, c("Historical","Forecast"))) |>
    arrange(series, year)
}

full_eu      <- build_full(hist_eu,      fc_eu)
full_sub     <- build_full(hist_sub,     fc_sub)
full_country <- build_full(hist_country, fc_country, key_col = "location_name")

write_csv(full_eu,      file.path(out_dir, "Forecast_EU_raw.csv"))
write_csv(full_sub,     file.path(out_dir, "Forecast_subregion_raw.csv"))
write_csv(full_country, file.path(out_dir, "Forecast_country_raw.csv"))

# Key-year summary
key_years <- c(2023, 2030, 2040, 2050)
fmt_series <- function(df) {
  df |> filter(year %in% key_years) |>
    mutate(
      `VLW (billion USD)` = fmt_ui(VLW,     VLW_lo,     VLW_hi,     2),
      `VLW/GDP (%)`       = fmt_ui(VLW_GDP, VLW_GDP_lo, VLW_GDP_hi, 4)
    ) |>
    select(Series = series, Year = year, Period = period,
           `VLW (billion USD)`, `VLW/GDP (%)`)
}

tbl_key <- bind_rows(
  fmt_series(full_eu)      |> mutate(Level = "EU Aggregate"),
  fmt_series(full_sub)     |> mutate(Level = "Subregion"),
  fmt_series(full_country) |> mutate(Level = "Country")
) |> arrange(Level, Series, Year)

write_csv(tbl_key, file.path(out_dir, "Forecast_summary_key_years.csv"))
cat("   CSVs saved (EU + subregion + 28 countries).\n\n")

# ===========================================================================
# 4. FIGURES
# ===========================================================================
cat("[4/5] Generating figures...\n")

# ── F1: EU-28 total a/b (VLW | VLW/GDP) ───────────────────────────────────
eu_h <- full_eu |> filter(period == "Historical")
eu_f <- full_eu |> filter(period == "Forecast")

p_eu_a <- make_panel(eu_h, eu_f, "#1B4F72",
                     "VLW", "VLW_lo", "VLW_hi",
                     "VLW (billion USD)", "EU-28 Total", pct = FALSE, tag = "a")
p_eu_b <- make_panel(eu_h, eu_f, "#1B4F72",
                     "VLW_GDP", "VLW_GDP_lo", "VLW_GDP_hi",
                     "VLW / GDP (%)", "EU-28 Total", pct = TRUE, tag = "b")

p_f1 <- ((p_eu_a / p_eu_b) + plot_layout(heights = c(1,1)))

pdf(file.path(out_dir, "Forecast_EU_ab.pdf"), width = 6, height = 7)
print(p_f1)
dev.off()
cat("   F1 saved: Forecast_EU_ab.pdf\n")

# ── F2: Subregion overlay on one plot (a = VLW, b = VLW/GDP) ──────────────
sub_h <- full_sub |> filter(period == "Historical") |>
  mutate(subregion = factor(series, levels = subregion_levels))
sub_f <- full_sub |> filter(period == "Forecast") |>
  mutate(subregion = factor(series, levels = subregion_levels))

p_sub_a <- ggplot() +
  geom_ribbon(data = sub_h,
              aes(x = year, ymin = VLW_lo, ymax = VLW_hi, fill = subregion),
              alpha = 0.15, colour = NA) +
  geom_ribbon(data = sub_f,
              aes(x = year, ymin = VLW_lo, ymax = VLW_hi, fill = subregion),
              alpha = 0.08, colour = NA) +
  geom_line(data = sub_h, aes(x = year, y = VLW, colour = subregion), linewidth = 0.8) +
  geom_line(data = sub_f, aes(x = year, y = VLW, colour = subregion),
            linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = SPLIT_YEAR, colour = "grey45", linewidth = 0.35,
             linetype = "dotted") +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.35) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.35) +
  scale_colour_manual(values = subregion_colors, name = NULL) +
  scale_fill_manual(values   = subregion_colors, name = NULL) +
  scale_x_continuous(limits = c(1990,2050),
                     breaks = c(1990,2000,2010,2023,2040,2050),
                     labels = c("1990","2000","2010","2023","2040","2050")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) round(x, 1)) +
  labs(x = NULL, y = "VLW (billion USD)", tag = "a") +
  theme_nm(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        axis.ticks  = element_line(colour = "grey30", linewidth = 0.3),
        panel.border = element_blank(), axis.line = element_blank(),
        legend.position = "top",
        plot.tag = element_text(size = 9, face = "bold"))

p_sub_b <- ggplot() +
  geom_ribbon(data = sub_h,
              aes(x = year, ymin = VLW_GDP_lo, ymax = VLW_GDP_hi, fill = subregion),
              alpha = 0.15, colour = NA) +
  geom_ribbon(data = sub_f,
              aes(x = year, ymin = VLW_GDP_lo, ymax = VLW_GDP_hi, fill = subregion),
              alpha = 0.08, colour = NA) +
  geom_line(data = sub_h, aes(x = year, y = VLW_GDP, colour = subregion),
            linewidth = 0.8) +
  geom_line(data = sub_f, aes(x = year, y = VLW_GDP, colour = subregion),
            linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = SPLIT_YEAR, colour = "grey45", linewidth = 0.35,
             linetype = "dotted") +
  geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.35) +
  geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.35) +
  scale_colour_manual(values = subregion_colors, name = NULL) +
  scale_fill_manual(values   = subregion_colors, name = NULL) +
  scale_x_continuous(limits = c(1990,2050),
                     breaks = c(1990,2000,2010,2023,2040,2050),
                     labels = c("1990","2000","2010","2023","2040","2050")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(formatC(x,"f",3),"%")) +
  labs(x = "Year", y = "VLW / GDP (%)", tag = "b") +
  theme_nm(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        axis.ticks  = element_line(colour = "grey30", linewidth = 0.3),
        panel.border = element_blank(), axis.line = element_blank(),
        legend.position = "none",
        plot.tag = element_text(size = 9, face = "bold"))

p_f2 <- ((p_sub_a / p_sub_b) + plot_layout(heights = c(1,1)))

pdf(file.path(out_dir, "Forecast_subregion_ab.pdf"), width = 7, height = 7)
print(p_f2)
dev.off()
cat("   F2 saved: Forecast_subregion_ab.pdf\n")

# ── F3: Subregion 3-panel (VLW) ────────────────────────────────────────────
sub_panels_vlw <- lapply(subregion_levels, function(s) {
  col <- subregion_colors[[s]]
  h   <- full_sub |> filter(series == s, period == "Historical")
  f   <- full_sub |> filter(series == s, period == "Forecast")
  make_panel(h, f, col, "VLW","VLW_lo","VLW_hi", "VLW (billion USD)", s)
})

pdf(file.path(out_dir, "Forecast_sub3panel_VLW.pdf"), width = 10, height = 4)
print(wrap_plots(sub_panels_vlw, ncol = 3))
dev.off()
cat("   F3 saved: Forecast_sub3panel_VLW.pdf\n")

# ── F4: Subregion 3-panel (VLW/GDP) ───────────────────────────────────────
sub_panels_gdp <- lapply(subregion_levels, function(s) {
  col <- subregion_colors[[s]]
  h   <- full_sub |> filter(series == s, period == "Historical")
  f   <- full_sub |> filter(series == s, period == "Forecast")
  make_panel(h, f, col, "VLW_GDP","VLW_GDP_lo","VLW_GDP_hi",
             "VLW / GDP (%)", s, pct = TRUE)
})

pdf(file.path(out_dir, "Forecast_sub3panel_GDP.pdf"), width = 10, height = 4)
print(wrap_plots(sub_panels_gdp, ncol = 3))
dev.off()
cat("   F4 saved: Forecast_sub3panel_GDP.pdf\n")

# ── F5: 28-country VLW forecast — 4×7 panel ───────────────────────────────
# Order countries by 2023 VLW descending
order_by_vlw <- full_country |>
  filter(year == 2023, period == "Historical") |>
  arrange(desc(VLW)) |>
  pull(series)

country_panels <- lapply(order_by_vlw, function(cty) {
  col_idx <- which(eu28_names == cty)
  col     <- colorRampPalette(c("#1B4F72","#AED6F1"))(28)[col_idx]
  h <- full_country |> filter(series == cty, period == "Historical")
  f <- full_country |> filter(series == cty, period == "Forecast")
  make_panel(h, f, "#2E86C1", "VLW","VLW_lo","VLW_hi",
             "VLW (bn USD)", cty)
})

pdf(file.path(out_dir, "Forecast_country_28panel.pdf"), width = 14, height = 10)
print(wrap_plots(country_panels, ncol = 7))
dev.off()
cat("   F5 saved: Forecast_country_28panel.pdf\n")

# ── F6: Top-10 countries ab (VLW / VLW/GDP) ───────────────────────────────
top10 <- full_country |>
  filter(year == 2023, period == "Historical") |>
  arrange(desc(VLW)) |>
  slice_head(n = 10) |>
  pull(series)

top10_colors <- setNames(
  colorRampPalette(c("#1B4F72","#117A65","#C0392B","#E67E22","#6C3483"))(10),
  top10
)

make_top10_panel <- function(y_col, ylo_col, yhi_col, y_lab, pct) {
  hist_t <- full_country |> filter(series %in% top10, period == "Historical") |>
    mutate(series = factor(series, levels = top10))
  fc_t   <- full_country |> filter(series %in% top10, period == "Forecast") |>
    mutate(series = factor(series, levels = top10))

  y_fmt <- if (pct) function(x) paste0(formatC(x,"f",3),"%") else
    function(x) ifelse(x >= 1000, paste0(round(x/1000,1),"K"), round(x,1))

  ggplot() +
    geom_ribbon(data = hist_t,
                aes(x = year, ymin = .data[[ylo_col]], ymax = .data[[yhi_col]],
                    fill = series), alpha = 0.14, colour = NA) +
    geom_ribbon(data = fc_t,
                aes(x = year, ymin = .data[[ylo_col]], ymax = .data[[yhi_col]],
                    fill = series), alpha = 0.07, colour = NA) +
    geom_line(data = hist_t, aes(x = year, y = .data[[y_col]], colour = series),
              linewidth = 0.75) +
    geom_line(data = fc_t, aes(x = year, y = .data[[y_col]], colour = series),
              linewidth = 0.75, linetype = "dashed") +
    geom_vline(xintercept = SPLIT_YEAR, colour = "grey45", linewidth = 0.35,
               linetype = "dotted") +
    geom_hline(yintercept = -Inf, colour = "grey30", linewidth = 0.35) +
    geom_vline(xintercept = -Inf, colour = "grey30", linewidth = 0.35) +
    scale_colour_manual(values = top10_colors, name = NULL) +
    scale_fill_manual(values   = top10_colors, name = NULL) +
    scale_x_continuous(limits = c(1990,2050),
                       breaks = c(1990,2000,2010,2023,2040,2050),
                       labels = c("1990","2000","2010","2023","2040","2050")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = y_fmt) +
    labs(x = "Year", y = y_lab) +
    theme_nm(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
          axis.ticks  = element_line(colour = "grey30", linewidth = 0.3),
          panel.border = element_blank(), axis.line = element_blank(),
          legend.position = "right")
}

p_t10a <- make_top10_panel("VLW",     "VLW_lo",     "VLW_hi",
                            "VLW (billion USD)", FALSE) + labs(tag = "a")
p_t10b <- make_top10_panel("VLW_GDP", "VLW_GDP_lo", "VLW_GDP_hi",
                            "VLW / GDP (%)",     TRUE)  + labs(tag = "b")

p_f6 <- ((p_t10a / p_t10b) + plot_layout(heights = c(1,1)))

pdf(file.path(out_dir, "Forecast_country_top10.pdf"), width = 8, height = 7)
print(p_f6)
dev.off()
cat("   F6 saved: Forecast_country_top10.pdf\n\n")

# ===========================================================================
# 5. SUMMARY
# ===========================================================================
cat("[5/5] Summary...\n")

eu_2050 <- full_eu |> filter(year == 2050)
eu_2023 <- full_eu |> filter(year == 2023)

cat(sprintf("\n   ┌──────────────────────────────────────────────────────┐\n"))
cat(sprintf("   │  EU-28 Pancreatic Cancer VLW Forecast (IE = %.1f)      │\n", IE))
cat(sprintf("   ├──────────────────────────────────────────────────────┤\n"))
cat(sprintf("   │  VLW 2023 (observed):  %.2f bn USD  (%.4f%% GDP)    │\n",
            eu_2023$VLW, eu_2023$VLW_GDP))
cat(sprintf("   │  VLW 2050 (forecast):  %.2f bn USD  (%.4f%% GDP)    │\n",
            eu_2050$VLW, eu_2050$VLW_GDP))
cat(sprintf("   │  Projected change:     %+.1f%%                         │\n",
            (eu_2050$VLW - eu_2023$VLW) / eu_2023$VLW * 100))
cat(sprintf("   └──────────────────────────────────────────────────────┘\n\n"))

cat(sprintf("   Output: %s\n", out_dir))
cat(sprintf("   CSVs:   Forecast_EU_raw | Forecast_subregion_raw | Forecast_country_raw\n"))
cat(sprintf("           Forecast_summary_key_years | ETS_model_summary\n"))
cat(sprintf("   PDFs:   F1 EU_ab | F2 subregion_ab | F3 sub3panel_VLW\n"))
cat(sprintf("           F4 sub3panel_GDP | F5 country_28panel | F6 country_top10\n\n"))
cat("══ Forecast complete ══\n")
