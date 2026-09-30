"""Regenerate revised Tables 1–3 from R1 country/year/sex calculations."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "analysis/Review_completion/R2_8_PSA/Input_country_year_sex_summary_bounds.csv"
data = pd.read_csv(source).query("sex_id == 3")

t1 = data.groupby("year", as_index=False).agg(
    DALYs=("D", "sum"),
    VLW_billion_2023_international_dollars=("VLW", "sum"),
    GDP_2023_PPP=("GDP", "sum"),
)
t2 = data.query("year == 2023")[[
    "location_id", "location_name", "subregion", "D", "VLW", "GDP"
]].rename(columns={
    "D": "DALYs", "VLW": "VLW_billion_2023_international_dollars",
    "GDP": "GDP_2023_PPP",
}).sort_values("location_name")
t3 = data.groupby(["subregion", "year"], as_index=False).agg(
    DALYs=("D", "sum"),
    VLW_billion_2023_international_dollars=("VLW", "sum"),
    GDP_2023_PPP=("GDP", "sum"),
)

for number, table in enumerate((t1, t2, t3), 1):
    table["VLW_relative_to_fixed_2023_GDP_pct"] = (
        100e9 * table.VLW_billion_2023_international_dollars / table.GDP_2023_PPP
    )
    table.to_csv(ROOT / f"analysis/Table_{number}_revised_points.csv", index=False)

assert abs(t1.loc[t1.year == 2023, "DALYs"].iloc[0] - 2311589.764641) < 1e-5
assert abs(t1.loc[t1.year == 2023, "VLW_billion_2023_international_dollars"].iloc[0] - 634.135197) < 1e-6
