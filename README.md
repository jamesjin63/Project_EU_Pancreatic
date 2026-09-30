# EU-28 Pancreatic Cancer Welfare Burden

Code for the revised manuscript *Pancreatic Cancer in the Original 28 European Union Members: Welfare Burden Beyond Conventional Metrics*. The analyses calculate pancreatic-cancer value of lost welfare (VLW) for 28 study countries from 1990 to 2023, project trends to 2050, and implement the revision's YLL/YLD and valuation sensitivity analyses. The study cohort includes the United Kingdom.

**This repository contains code only.** Source data, generated CSV tables, simulation draws, maps, and rendered figures are not distributed here.

## Code structure

```text
scripts/
├── run_pancreatic_EU_VLW.R     Main VLW calculations; country, region and sex plots
├── run_age_VLW_EU_2023.R      Age and sex calculations and plots
├── run_forecast_EU_2050.R    ETS forecasting and plots
├── formulas.R                 Valuation, discounting and ratio functions
├── prepare_source_data.py     Build local study inputs from separately obtained data
├── build_main_tables.py       Revised Tables 1–3
└── build_vsl_sensitivity.py   Reference-VSL sensitivity table
analysis/
├── R1_1_YLL_YLD/
│   ├── run_R1_1_YLL_YLD.R      DALY = YLL + YLD and YLL-based VLW analysis
│   └── plot_R1_1_YLL_YLD.R     Decomposition and sensitivity plots
└── Review_completion/
    ├── run_remaining_analyses.R  Age valuation, discounting, income-ratio checks, PSA
    └── plot_remaining_analyses.R Corrected Figure 1 and additional revision plots
```

## Data sources and availability

- **Disease burden and HALE:** GBD estimates from the [IHME GBD Results Tool](https://ghdx.healthdata.org/gbd-results-tool). Obtain the version and selections matching the study: pancreatic cancer; deaths, DALYs, YLLs and YLDs; Number metric; 1990–2023; study countries; age and sex strata. HALE inputs use 2023 estimates. These estimates can be queried and downloaded through IHME, subject to its [data-use agreement](https://www.ihmeclientservices.org/noncommercial.html); public access does not imply unrestricted redistribution or commercial use.
- **Economic inputs:** 2023 GDP per capita, PPP (`NY.GDP.PCAP.PP.CD`), and GDP, PPP (`NY.GDP.MKTP.PP.CD`), from the [World Bank World Development Indicators](https://datacatalog.worldbank.org/search/dataset/0037712/world-development-indicators). WDI is publicly available under CC BY 4.0 with attribution.
- **Map geometry:** Provide a country-boundary GeoJSON with names compatible with the plotting script. The geometry file used for the original analysis is not distributed here.

The source project contained merged GBD-format CSVs, not the pre-merge download files or their merge script. Reproducing the exact published values requires the same GBD release, query settings and input versions. The code repository does not make those source files public or grant rights to redistribute them.

## Local requirements

Use R with `tidyverse`, `sf`, `patchwork`, `scales`, `forecast`, `data.table`, `ggplot2`, and `ragg`; Python 3 with `pandas`; and a UTF-8 locale. Supply the full source files locally, then run `scripts/prepare_source_data.py --help` for its required paths. It creates the ignored `data/` inputs and `EU28_location_list.csv`. Run the main R scripts from the repository root before the revision scripts; generated files are written to the ignored `results_VLW/` and `analysis/` output paths. Income elasticity is passed to the three main R scripts as `1` (main analysis), `0.5`, or `1.5`.

## License

The code in this repository is licensed under the [Apache License 2.0](LICENSE). External datasets retain their own terms of use; this license does not apply to data obtained from IHME, the World Bank, or map providers.
