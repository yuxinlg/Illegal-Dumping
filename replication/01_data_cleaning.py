"""
01_data_cleaning.py
-------------------
Load, merge, and clean all PPR 311 illegal-dumping records (2018–2025).
Outputs:
    output/illegal_dumping_full.geojson   -- main event dataset used by 03_analysis.py

Data sources required (not tracked in git, see README.md for download instructions):
    data/311/PPR3112016_to_2019/PPR311DumpingDataforDDDI_CY2018.csv
    data/311/PPR3112016_to_2019/PPR311DumpingDataforDDDI_CY2019.csv
    data/311/PPR3112020_to_2022/PPR311DumpingDataforDDDI_CY2020.csv
    data/311/PPR3112020_to_2022/PPR311DumpingDataforDDDI_CY2021.csv
    data/311/PPR3112020_to_2022/PPR311DumpingDataforDDDI_CY2022.csv
    data/311/PPR3112023_to_2025/PPR311DumpingDataforDDDI_CY2023.csv
    data/311/PPR3112023_to_2025/PPR311DumpingDataforDDDI_CY2024.csv
    data/311/PPR3112023_to_2025/PPR311DumpingDataforDDDI_CY2025.csv
    data/311/PPR311DumpingDataforDDDI_Aug2025toDec2025.csv

Additional input files (tracked in output/ as they are derived artifacts):
    output/picked_yes.csv        -- manually confirmed "Illegal Dumping" cases
    output/dfv_CRT.csv           -- pending cases confirmed through case review
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd

os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")

# ── Column rename map ─────────────────────────────────────────────────────────

RENAME_MAP = {
    "Date/Time Opened":                        "start_time",
    "Date/Time Closed":                        "close_time",
    "Service Request Number":                  "request_id",
    "Service Request Type":                    "type",
    "Department":                              "department",
    "Address/Intersection":                    "address",
    "ZipCode":                                 "zipcode",
    "Status":                                  "status",
    "Problem Category":                        "problem_category",
    "Problem":                                 "problem",
    "Location Type":                           "location_type",
    "Park District":                           "park_district",
    "Large Crew&Heavy Machinery To Clean Up?": "large_crew",
    "Type of Materials":                       "material_type",
    "Condition of Materials":                  "material_condition",
    "Description":                             "description",
    "Centerline (Latitude)":                   "lat",
    "Centerline (Longitude)":                  "lon",
    "Case ID":                                 "case_id",
    "Parent Case ID":                          "parent_case_id",
    "Mobile App Photos":                       "media_url",
}

FILE_PATHS = {
    "ppr311_18": "data/311/PPR3112016_to_2019/PPR311DumpingDataforDDDI_CY2018.csv",
    "ppr311_19": "data/311/PPR3112016_to_2019/PPR311DumpingDataforDDDI_CY2019.csv",
    "ppr311_20": "data/311/PPR3112020_to_2022/PPR311DumpingDataforDDDI_CY2020.csv",
    "ppr311_21": "data/311/PPR3112020_to_2022/PPR311DumpingDataforDDDI_CY2021.csv",
    "ppr311_22": "data/311/PPR3112020_to_2022/PPR311DumpingDataforDDDI_CY2022.csv",
    "ppr311_23": "data/311/PPR3112023_to_2025/PPR311DumpingDataforDDDI_CY2023.csv",
    "ppr311_24": "data/311/PPR3112023_to_2025/PPR311DumpingDataforDDDI_CY2024.csv",
    "ppr311_25": "data/311/PPR3112023_to_2025/PPR311DumpingDataforDDDI_CY2025.csv",
    "ppr311_25_aug_dec": "data/311/PPR311DumpingDataforDDDI_Aug2025toDec2025.csv",
}


# ── Step 1: Load and clean each year ─────────────────────────────────────────

def load_year(var_name: str, path: str) -> gpd.GeoDataFrame:
    """Read one CSV file, standardize columns, clean coordinates, return GeoDataFrame."""
    with open(path, encoding="utf-8", errors="replace") as f:
        df = pd.read_csv(f).rename(columns=RENAME_MAP)

    year = int(var_name.split("_")[1]) + 2000
    df["year"] = year

    for col in ("start_time", "close_time"):
        if col in df.columns:
            df[col] = pd.to_datetime(df[col], errors="coerce")
            df[col] = df[col].dt.strftime("%Y/%m/%d %H:%M")

    df[["lon", "lat"]] = df[["lon", "lat"]].apply(pd.to_numeric, errors="coerce")
    df.loc[(df["lon"] == 0) & (df["lat"] == 0), ["lon", "lat"]] = np.nan
    df = df[df["lon"].between(-76, -74) & df["lat"].between(39, 41)].dropna(
        subset=["lon", "lat"]
    )

    return gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"]),
        crs=4326,
    ).to_crs(26918)


print("Loading 311 CSVs …")
yearly = {name: load_year(name, path) for name, path in FILE_PATHS.items()}


# ── Step 2: Deduplicate overlapping 2025 files ───────────────────────────────

CUTOFF_2025 = "2025/08/01 00:00"
yearly["ppr311_25"] = yearly["ppr311_25"][
    yearly["ppr311_25"]["start_time"] < CUTOFF_2025
].copy()

ppr311_25_full = pd.concat(
    [yearly["ppr311_25"], yearly["ppr311_25_aug_dec"]],
    ignore_index=True,
).sort_values("start_time").reset_index(drop=True)


# ── Step 3: Concatenate analysis years (2021–2025) ────────────────────────────

illegal_dumping = pd.concat(
    [yearly["ppr311_21"], yearly["ppr311_22"], yearly["ppr311_23"],
     yearly["ppr311_24"], ppr311_25_full],
    ignore_index=True,
)
print(f"Total records (2021–2025): {illegal_dumping.shape[0]:,}")


# ── Step 4: Filter to confirmed illegal-dumping cases ────────────────────────

# 4a. Cases flagged as 'Illegal Dumping' by service type
illegal_dumping_ppr = illegal_dumping[illegal_dumping["type"] == "Illegal Dumping"]

# 4b. Manually labelled confirmed cases
picked_yes = pd.read_csv("output/picked_yes.csv", header=None, names=["request_id"])
picked_yes["request_id"] = picked_yes["request_id"].astype(str)
illegal_dumping["request_id"] = illegal_dumping["request_id"].astype(str)
illegal_dumping_yes = illegal_dumping.merge(picked_yes, on="request_id", how="right")

# 4c. Pending cases confirmed through manual case review
pending_confirmed = pd.read_csv("output/dfv_CRT.csv")
pending_confirmed = pending_confirmed[pending_confirmed["confirm"] == 1][["request_id"]]
pending_confirmed["request_id"] = pending_confirmed["request_id"].astype(str)
illegal_dumping_confirmed = illegal_dumping.merge(
    pending_confirmed, on="request_id", how="right"
)

# 4d. Merge all confirmed cases
illegal_dumping_full = pd.concat(
    [illegal_dumping_ppr, illegal_dumping_yes, illegal_dumping_confirmed],
    ignore_index=True,
).drop_duplicates(subset=["request_id"])

print(f"Confirmed illegal-dumping events: {illegal_dumping_full.shape[0]:,}")


# ── Step 5: Export ────────────────────────────────────────────────────────────

os.makedirs("output", exist_ok=True)
illegal_dumping_full.to_file("output/illegal_dumping_full.geojson", driver="GeoJSON")
print("Saved → output/illegal_dumping_full.geojson")
