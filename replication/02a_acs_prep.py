"""
02a_acs_prep.py
===============
Replication script for ACS covariate construction.

Workflow
--------
  1. Load ACS data  — from Census API (optional) or static raw files in replication/data/acs/
  2. Attach CBG geometry and build GEOID
  3. Merge 4-year ACS into a single cross-sectional panel (2021–2024 averages)
  4. Add area and population density
  5. Export final GeoJSON to replication/data/

Required datasets (place in replication/data/ before running):
  - acs/raw_acs_21.csv … acs/raw_acs_24.csv       (ACS raw pulls; or re-download via API)
  - tiger/tl_2023_42_bg/tl_2023_42_bg.shp         (TIGER CBG shapefile for PA)

Output
------
  replication/data/philly_cbg_acs.geojson
"""

# ---------------------------------------------------------------------------
# 0. Imports & paths
# ---------------------------------------------------------------------------
import os
import pandas as pd
import geopandas as gpd
import numpy as np

REPL_DIR  = os.path.dirname(os.path.abspath(__file__))
REPL_DATA = os.path.join(REPL_DIR, "data")
ACS_DATA  = os.path.join(REPL_DATA, "acs")
TIGER_DIR = os.path.join(REPL_DATA, "tiger", "tl_2023_42_bg")

os.makedirs(ACS_DATA, exist_ok=True)

YEARS = [2021, 2022, 2023, 2024]


# ---------------------------------------------------------------------------
# 1. Load ACS data
#    Option A (default): read static CSVs saved from the Census API pull.
#    Option B:           re-download from the Census API by setting
#                        USE_API = True and supplying your API key below.
# ---------------------------------------------------------------------------
USE_API = False          # set True to re-download from Census API
API_KEY = "YOUR_CENSUS_API_KEY_HERE"

# ACS variables pulled at the census block-group level (Philadelphia, PA 42101)
ACS_FIELDS = [
    "NAME", "B01003_001E", "B19013_001E",
    # Limited English Proficiency — ages 5-17
    "B16004_007E", "B16004_008E",   # Spanish
    "B16004_012E", "B16004_013E",   # Other Indo-European
    "B16004_017E", "B16004_018E",   # Asian/PI
    "B16004_022E", "B16004_023E",   # Other
    # LEP ages 18-64
    "B16004_029E", "B16004_030E",
    "B16004_034E", "B16004_035E",
    "B16004_039E", "B16004_040E",
    "B16004_044E", "B16004_045E",
    # LEP ages 65+
    "B16004_051E", "B16004_052E",
    "B16004_056E", "B16004_057E",
    "B16004_061E", "B16004_062E",
    "B16004_066E", "B16004_067E",
    # Education (B15003 universe = population 25+)
    "B15003_001E",                  # total
    "B15003_002E",  "B15003_003E",  "B15003_004E",  "B15003_005E",
    "B15003_006E",  "B15003_007E",  "B15003_008E",  "B15003_009E",
    "B15003_010E",  "B15003_011E",  "B15003_012E",  "B15003_013E",
    "B15003_014E",  "B15003_015E",  "B15003_016E",  "B15003_017E",
    "B15003_018E",  "B15003_019E",  "B15003_020E",  "B15003_021E",
    "B15003_022E",  "B15003_023E",  "B15003_024E",  "B15003_025E",
    # Median age
    "B01002_001E",
]

ACS_GEO = {"for": "block group:*", "in": "state:42 county:101"}


def download_acs_year(year: int) -> pd.DataFrame:
    """Pull one year of ACS 5-year estimates from the Census API."""
    from census import Census
    c = Census(API_KEY)
    data = c.acs5.get(ACS_FIELDS, ACS_GEO, year=year)
    return pd.DataFrame(data)


def load_acs_year(year: int) -> pd.DataFrame:
    """Load ACS data for one year — from API or saved CSV."""
    yr2 = str(year)[-2:]
    csv_path = os.path.join(ACS_DATA, f"raw_acs_{yr2}.csv")

    if USE_API:
        print(f"  Downloading ACS {year} from Census API...")
        df = download_acs_year(year)
        df.to_csv(csv_path, index=False)
        print(f"  Saved to {csv_path}")
    else:
        print(f"  Loading ACS {year} from {csv_path}")
        df = pd.read_csv(csv_path, dtype=str)
        numeric_cols = [c for c in df.columns
                        if c not in ("NAME", "GEO_ID", "state", "county", "tract", "block group")]
        df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")

    return df


def clean_acs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename key columns, derive covariates, and handle suppressed/invalid values.

    Derived columns
    ---------------
    lep    : count of limited-English-proficient residents
             (density per km² is computed later in add_area_and_density)
    edu_hs : share of population 25+ with HS diploma or higher
    """
    df = df.copy()

    df.rename(columns={
        "B01003_001E": "pop",
        "B19013_001E": "med_inc",
        "B01002_001E": "med_age",
        "B15003_001E": "edu_total",
    }, inplace=True)

    # Median income: negative = suppressed/unreliable → impute with column mean
    df.loc[df["med_inc"] < 0, "med_inc"] = np.nan
    df["med_inc"] = df["med_inc"].fillna(df["med_inc"].mean())

    # Median age: implausibly low values → impute with column mean
    df.loc[df["med_age"] < 5, "med_age"] = np.nan
    df["med_age"] = df["med_age"].fillna(df["med_age"].mean())

    # Limited English Proficiency: sum across all language groups and age brackets
    lep_vars = [
        "B16004_007E", "B16004_008E", "B16004_012E", "B16004_013E",
        "B16004_017E", "B16004_018E", "B16004_022E", "B16004_023E",
        "B16004_029E", "B16004_030E", "B16004_034E", "B16004_035E",
        "B16004_039E", "B16004_040E", "B16004_044E", "B16004_045E",
        "B16004_051E", "B16004_052E", "B16004_056E", "B16004_057E",
        "B16004_061E", "B16004_062E", "B16004_066E", "B16004_067E",
    ]
    df["lep"] = df[lep_vars].sum(axis=1)
    df.loc[df["lep"] < 0, "lep"] = np.nan
    df["lep"] = df["lep"].fillna(df["lep"].mean())

    # Education: share with HS diploma or higher (GED through Doctorate)
    edu_hs_vars = [f"B15003_{i:03d}E" for i in range(17, 26)]  # 017–025
    df["edu_hs"] = df[edu_hs_vars].sum(axis=1) / df["edu_total"]
    df["edu_hs"] = df["edu_hs"].fillna(df["edu_hs"].mean())

    return df


# ---------------------------------------------------------------------------
# 2. Attach CBG geometry and build GEOID
# ---------------------------------------------------------------------------

def attach_geometry(df: pd.DataFrame, philly_shp: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Merge ACS tabular data onto the Philadelphia CBG shapefile.
    Keeps geometry, FIPS identifiers, population, and derived covariates.
    Adds a 12-digit GEOID (state + county + tract + block group).
    Returns a GeoDataFrame in EPSG:26918 (NAD83 UTM Zone 18N).
    """
    gdf = philly_shp.merge(
        df,
        left_on=["TRACTCE", "BLKGRPCE"],
        right_on=["tract", "block group"],
        how="inner",
    )
    gdf = gdf.to_crs("EPSG:26918")

    keep_cols = ["geometry", "state", "county", "tract", "block group",
                 "pop", "med_inc", "med_age", "lep", "edu_hs"]
    gdf = gdf[keep_cols].copy()

    gdf["GEOID"] = (
        gdf["state"].astype(str).str.zfill(2)
        + gdf["county"].astype(str).str.zfill(3)
        + gdf["tract"].astype(str).str.zfill(6)
        + gdf["block group"].astype(str)
    )
    return gdf


# ---------------------------------------------------------------------------
# 3. Merge 4-year ACS into cross-sectional averages
# ---------------------------------------------------------------------------

def merge_4year_acs(acs_gdfs: dict) -> pd.DataFrame:
    """
    Wide-merge acs_21 … acs_24 GeoDataFrames on FIPS keys, then
    compute 4-year averages for each covariate.

    Parameters
    ----------
    acs_gdfs : dict  {year (int): GeoDataFrame}

    Returns
    -------
    DataFrame with columns: state, county, tract, block group,
        pop_avg, med_inc_avg, med_age_avg, lep_avg, edu_hs_avg
    """
    merge_keys = ["state", "county", "tract", "block group"]
    drop_cols  = ["geometry", "GEOID", "NAME"]
    covariates = ["pop", "med_inc", "med_age", "lep", "edu_hs"]

    def _suffix(df, yr):
        yr2 = str(yr)[-2:]
        df  = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")
        return df.rename(columns={c: f"{c}_{yr2}" for c in covariates if c in df.columns})

    base = _suffix(acs_gdfs[2021], 2021)
    for yr in [2022, 2023, 2024]:
        base = base.merge(_suffix(acs_gdfs[yr], yr), on=merge_keys, how="inner")

    for cov in covariates:
        cols = [f"{cov}_{str(yr)[-2:]}" for yr in YEARS]
        base[f"{cov}_avg"] = base[cols].mean(axis=1)

    avg_cols = merge_keys + [f"{c}_avg" for c in covariates]
    panel = base[avg_cols].copy()

    for col in [f"{c}_avg" for c in covariates]:
        panel[col] = panel[col].fillna(panel[col].mean())

    return panel


# ---------------------------------------------------------------------------
# 4. Add area and population density
# ---------------------------------------------------------------------------

def add_area_and_density(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Compute CBG area (m² and km²), population density, and LEP density.
    Requires EPSG:26918 (metres).

    Derived columns
    ---------------
    pop_density : residents per km²
    lep_density : LEP residents per km²
    """
    gdf = gdf.copy()
    gdf["area_m2"]     = gdf.geometry.area
    gdf["area_km2"]    = gdf["area_m2"] / 1_000_000
    gdf["pop_density"] = gdf["pop_avg"] / gdf["area_km2"]
    gdf["lep_density"] = gdf["lep_avg"] / gdf["area_km2"]
    return gdf


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("ACS Covariates — Replication Pipeline")
    print("=" * 60)

    # --- 1. Load ACS raw data ---
    print("\n[1/4] Loading ACS data...")
    raw_acs = {}
    for yr in YEARS:
        print(f"  Year {yr}:")
        raw_acs[yr] = clean_acs(load_acs_year(yr))

    # --- 2. Attach CBG geometry ---
    print("\n[2/4] Attaching CBG geometry (TIGER 2023)...")
    pa_shp     = gpd.read_file(os.path.join(TIGER_DIR, "tl_2023_42_bg.shp"))
    philly_shp = pa_shp[pa_shp["COUNTYFP"] == "101"].copy()

    acs_geo = {}
    for yr in YEARS:
        acs_geo[yr] = attach_geometry(raw_acs[yr], philly_shp)
        print(f"  acs_{str(yr)[-2:]} with geometry: {len(acs_geo[yr])} CBGs")

    # --- 3. Merge 4-year ACS averages ---
    print("\n[3/4] Merging 4-year ACS into cross-sectional averages...")
    acs_4year = merge_4year_acs(acs_geo)
    print(f"  Panel shape: {acs_4year.shape}")

    acs_4year_geo = philly_shp.to_crs("EPSG:26918").merge(
        acs_4year,
        left_on=["TRACTCE", "BLKGRPCE"],
        right_on=["tract", "block group"],
        how="inner",
    )
    acs_4year_geo = acs_4year_geo[
        ["geometry", "state", "county", "tract", "block group",
         "pop_avg", "med_inc_avg", "med_age_avg", "lep_avg", "edu_hs_avg"]
    ].copy()

    acs_4year_geo["GEOID"] = (
        acs_4year_geo["state"].astype(str).str.zfill(2)
        + acs_4year_geo["county"].astype(str).str.zfill(3)
        + acs_4year_geo["tract"].astype(str).str.zfill(6)
        + acs_4year_geo["block group"].astype(str)
    )

    # --- 4. Area and population density ---
    print("\n[4/4] Computing CBG area and population density...")
    final = add_area_and_density(acs_4year_geo)
    print(f"  Final dataset: {len(final)} CBGs, {len(final.columns)} columns")
    print(f"  Columns: {list(final.columns)}")

    # --- Export ---
    repl_out = os.path.join(REPL_DIR, "output")
    os.makedirs(repl_out, exist_ok=True)
    out_path = os.path.join(repl_out, "philly_cbg_acs.geojson")
    final.to_file(out_path, driver="GeoJSON")
    print(f"\nSaved: {out_path}")
    print("Done.")


if __name__ == "__main__":
    main()
