"""Rebuild Table S7 from the three income-elasticity forecast archives."""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
rows = []
for ie in (0.5, 1.0, 1.5):
    source = ROOT / f"results_VLW/forecast_EU_IE{ie:.1f}/Forecast_EU_raw.csv"
    with source.open(newline="") as handle:
        base = next(row for row in csv.DictReader(handle) if row["year"] == "2023")
    for reference_millions in (7, 10, 13.2, 14):
        scale = reference_millions / 13.2
        rows.append({
            "VSL_reference_million_2023_int_dollars": reference_millions,
            "IE": ie,
            "VLW_billion": float(base["VLW"]) * scale,
            "VLW_GDP_percent": float(base["VLW_GDP"]) * scale,
        })
destination = ROOT / "analysis/Response_support/Table_R1_reference_VSL_sensitivity.csv"
with destination.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0])
    writer.writeheader()
    writer.writerows(rows)
