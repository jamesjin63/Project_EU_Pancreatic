"""Create local study inputs from separately obtained source data files.

No source data or derived CSV files are included in this code repository.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

STUDY_COUNTRIES = (
    (75, "Austria"), (76, "Belgium"), (45, "Bulgaria"),
    (46, "Croatia"), (77, "Cyprus"), (47, "Czechia"),
    (78, "Denmark"), (58, "Estonia"), (79, "Finland"),
    (80, "France"), (81, "Germany"), (82, "Greece"),
    (48, "Hungary"), (84, "Ireland"), (86, "Italy"),
    (59, "Latvia"), (60, "Lithuania"), (87, "Luxembourg"),
    (88, "Malta"), (89, "Netherlands"), (51, "Poland"),
    (91, "Portugal"), (52, "Romania"), (54, "Slovakia"),
    (55, "Slovenia"), (92, "Spain"), (93, "Sweden"),
    (95, "United Kingdom"),
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract(source, target, keep):
    count = 0
    with source.open(newline="", encoding="utf-8-sig") as input_file, \
            target.open("w", newline="", encoding="utf-8") as output_file:
        reader = csv.DictReader(input_file)
        writer = csv.DictWriter(output_file, fieldnames=reader.fieldnames)
        writer.writeheader()
        for record in reader:
            if keep(record):
                writer.writerow(record)
                count += 1
    return {
        "input": target.name,
        "rows": count,
        "source_sha256": sha256(source),
        "subset_sha256": sha256(target),
        "subset_bytes": target.stat().st_size,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("merged", "yll_yld", "hale", "gdp", "geojson"):
        parser.add_argument(f"--{name.replace('_', '-')}", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    destination = root / "data"
    destination.mkdir(exist_ok=True)
    with (root / "EU28_location_list.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("location_id", "location_name", "group_type"))
        writer.writerow((4743, "European Union", "EU-28 aggregate"))
        writer.writerows((identifier, name, "EU-28 member")
                         for identifier, name in STUDY_COUNTRIES)
    ids = {identifier for identifier, _ in STUDY_COUNTRIES}
    names = {name for _, name in STUDY_COUNTRIES} | {"Slovak Republic"}
    manifest = [
        extract(args.merged, destination / "merged.csv",
                lambda row: int(row["cause_id"]) == 456
                and int(row["metric_id"]) == 1
                and int(row["measure_id"]) in (1, 2)
                and int(row["location_id"]) in ids | {4743}),
        extract(args.yll_yld,
                destination / "GBD2023_Pancreatic_cancer_32locations_YLLs_YLDs_1990_2023.csv",
                lambda row: int(row["cause_id"]) == 456
                and int(row["metric_id"]) == 1
                and int(row["measure_id"]) in (3, 4)
                and int(row["location_id"]) in ids),
        extract(args.hale, destination / "HALE.csv",
                lambda row: int(row["year"]) == 2023
                and int(row["location_id"]) in ids
                and row["metric_name"] == "Years"),
        extract(args.gdp, destination / "gdp.csv",
                lambda row: int(row["year"]) == 2023
                and row["country"] in names),
    ]
    shutil.copy2(args.geojson, destination / "df_world2.geojson")
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
