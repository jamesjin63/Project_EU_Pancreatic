# Mathematical definitions used in the revised EU-28 pancreatic cancer article.
# Currency inputs are 2023 PPP international dollars; VLW output is billions.

reference_vsl <- 13.2e6
reference_gdppc <- 82304.62

vsl <- function(gdppc, income_elasticity = 1) {
  reference_vsl * (gdppc / reference_gdppc)^income_elasticity
}

vsly <- function(gdppc, hale, income_elasticity = 1,
                 annualization_fraction = 0.5) {
  stopifnot(all(hale > 0), annualization_fraction > 0)
  vsl(gdppc, income_elasticity) / (annualization_fraction * hale)
}

vlw_billions <- function(dalys, gdppc, hale, income_elasticity = 1,
                         annualization_fraction = 0.5) {
  dalys * vsly(gdppc, hale, income_elasticity, annualization_fraction) / 1e9
}

vlw_gdp_percent <- function(vlw_billions, gdp_2023) {
  stopifnot(all(gdp_2023 > 0))
  100 * vlw_billions * 1e9 / gdp_2023
}

discount_multiplier <- function(hale, horizon_fraction = 0.4,
                                annual_rate = 0.03) {
  n <- horizon_fraction * hale
  ifelse(annual_rate == 0, 1,
         (1 - (1 + annual_rate)^(-n)) / (annual_rate * n))
}

annuity_denominator <- function(hale, annual_rate = 0.03,
                                horizon_fraction = 0.5) {
  n <- horizon_fraction * hale
  ifelse(annual_rate == 0, n,
         (1 - (1 + annual_rate)^(-n)) / annual_rate)
}

component_share_percent <- function(component, dalys) {
  ifelse(dalys == 0, NA_real_, 100 * component / dalys)
}
