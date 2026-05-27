"""
acs_draft.py
============
Replication script for ACS covariates and 311 reporting rate construction.

Workflow
--------
  1. Load ACS data  — from Census API (optional) or static raw files in replication/data/
  2. Attach CBG geometry and build GEOID
  3. Load & clean 311 data — filter out Illegal Dumping calls
  4. Spatial join: count 311 reports per CBG per year
  5. Impute reporting rate for zero-population CBGs via spatial neighbor median
  6. Merge 4-year ACS into a single cross-sectional panel (2021–2024 averages)
  7. Add area and population density
  8. Merge 4-year average 311 reporting rate
  9. Export final GeoJSON to replication/data/

Required datasets (place in replication/data/ before running):
  - acs/raw_acs_21.csv … acs/raw_acs_24.csv       (ACS raw pulls; or re-download via API)
  - tiger/tl_2023_42_bg/tl_2023_42_bg.shp        (TIGER CBG shapefile for PA)
  - 311/{2021..2024}/public_cases_fc.shp          (Philly 311 shapefiles, one per year)

Output
------
  replication/data/philly_cbg_acs_311.geojson
"""

# ---------------------------------------------------------------------------
# 0. Imports & paths
# ---------------------------------------------------------------------------
import os
import pandas as pd
import geopandas as gpd
import numpy as np
from census import Census

# All paths are relative to the project root (Illegal-Dumping/).
# Run this script from the project root:  python replication/acs_draft.py
REPL_DIR     = os.path.dirname(os.path.abspath(__file__))
REPL_DATA    = os.path.join(REPL_DIR, "data")
ACS_DATA     = os.path.join(REPL_DATA, "acs")
TIGER_DIR    = os.path.join(REPL_DATA, "tiger", "tl_2023_42_bg")
PC311_DIR    = os.path.join(REPL_DATA, "311")

os.makedirs(ACS_DATA, exist_ok=True)

YEARS = [2021, 2022, 2023, 2024]


# ---------------------------------------------------------------------------
# 1. Load ACS data
#    Option A (default): read static CSVs saved from the Census API pull.
#    Option B:           re-download from the Census API by setting
#                        USE_API = True and supplying your API key below.
# ---------------------------------------------------------------------------
USE_API  = False          # set True to re-download from Census API
API_KEY  = "YOUR_CENSUS_API_KEY_HERE"

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
        # Cast numeric columns back from string
        numeric_cols = [c for c in df.columns if c not in ("NAME", "GEO_ID", "state", "county", "tract", "block group")]
        df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")

    return df


def clean_acs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename key columns, derive covariates, and handle suppressed/invalid values.

    Derived columns
    ---------------
    lep           : count of limited-English speakers (all age groups)
    lep_rate      : lep / pop  (0 for zero-pop CBGs)
    edu_hs        : share of population 25+ with HS diploma or higher
    """
    df = df.copy()

    # Rename core variables
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
        # 5-17
        "B16004_007E", "B16004_008E", "B16004_012E", "B16004_013E",
        "B16004_017E", "B16004_018E", "B16004_022E", "B16004_023E",
        # 18-64
        "B16004_029E", "B16004_030E", "B16004_034E", "B16004_035E",
        "B16004_039E", "B16004_040E", "B16004_044E", "B16004_045E",
        # 65+
        "B16004_051E", "B16004_052E", "B16004_056E", "B16004_057E",
        "B16004_061E", "B16004_062E", "B16004_066E", "B16004_067E",
    ]
    df["lep"] = df[lep_vars].sum(axis=1)
    df["lep_rate"] = df["lep"] / df["pop"]
    df.loc[df["pop"] == 0, "lep_rate"] = 0
    df["lep_rate"] = df["lep_rate"].fillna(df["lep_rate"].mean())

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
    Keeps only geometry, FIPS identifiers, population, and derived covariates.
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

    keep_cols = ["geometry", "state", "county", "tract", "block group", "pop",
                 "med_inc", "med_age", "lep_rate", "edu_hs"]
    gdf = gdf[keep_cols].copy()

    gdf["GEOID"] = (
        gdf["state"].astype(str).str.zfill(2)
        + gdf["county"].astype(str).str.zfill(3)
        + gdf["tract"].astype(str).str.zfill(6)
        + gdf["block group"].astype(str)
    )
    return gdf


# ---------------------------------------------------------------------------
# 3. Load & clean 311 data
# ---------------------------------------------------------------------------

def load_311_year(year: int) -> gpd.GeoDataFrame:
    """
    Load Philly 311 shapefile for one year.
    Removes 'Illegal Dumping' service requests (those are modeled separately)
    and reprojects to EPSG:26918.
    """
    shp_path = os.path.join(PC311_DIR, str(year), "public_cases_fc.shp")
    print(f"  Loading 311 data: {shp_path}")
    gdf = gpd.read_file(shp_path)

    n_before = len(gdf)
    gdf = gdf[gdf["service_na"] != "Illegal Dumping"].copy()
    print(f"  {year}: removed {n_before - len(gdf)} 'Illegal Dumping' rows "
          f"({len(gdf)} remaining)")

    gdf["year"] = year
    gdf = gdf.to_crs("EPSG:26918")
    return gdf


# ---------------------------------------------------------------------------
# 4. Spatial join: count 311 reports per CBG and compute reporting rate
# ---------------------------------------------------------------------------

def compute_reporting_rate(
    acs_gdf: gpd.GeoDataFrame,
    pc311_gdf: gpd.GeoDataFrame,
    year: int,
) -> gpd.GeoDataFrame:
    """
    Assign each 311 report to a CBG (spatial join), count reports per CBG,
    and compute reporting_rate_per_1000 = (reports / pop) * 1000.
    Zero-population CBGs receive NaN (imputed later).
    """
    # Align CRS
    if pc311_gdf.crs != acs_gdf.crs:
        pc311_gdf = pc311_gdf.to_crs(acs_gdf.crs)

    # Assign each report to the CBG it falls within
    reports_in_cbg = gpd.sjoin(
        pc311_gdf,
        acs_gdf[["geometry", "GEOID"]],
        how="left",
        predicate="within",
    )

    assigned   = reports_in_cbg["GEOID"].notna().sum()
    unassigned = reports_in_cbg["GEOID"].isna().sum()
    print(f"  {year}: {len(pc311_gdf)} total reports | "
          f"{assigned} assigned to CBGs | {unassigned} unassigned")

    # Count per CBG
    report_counts = (
        reports_in_cbg.groupby("GEOID")
        .size()
        .reset_index(name="total_reports")
    )

    # Merge counts back to ACS GeoDataFrame
    result = acs_gdf.merge(report_counts, on="GEOID", how="left")
    result["total_reports"] = result["total_reports"].fillna(0)

    # Rate per 1,000 residents; NaN for zero-pop CBGs (imputed in next step)
    result["reporting_rate_per_1000"] = result.apply(
        lambda row: (row["total_reports"] / row["pop"] * 1000)
        if row["pop"] > 0 else float("nan"),
        axis=1,
    )

    zero_pop_with_reports = (
        (result["pop"] == 0) & (result["total_reports"] > 0)
    ).sum()
    if zero_pop_with_reports:
        print(f"  Warning: {zero_pop_with_reports} CBGs have reports but zero population "
              "(will be imputed from neighbors)")

    return result


# ---------------------------------------------------------------------------
# 5. Impute reporting rate for zero-population CBGs
# ---------------------------------------------------------------------------

def _build_neighbor_map(gdf: gpd.GeoDataFrame, buffer_m: float = 1.0) -> dict:
    """Return {row_index: [neighbor_row_indices]} using a 1-metre buffer."""
    gdf = gdf.reset_index(drop=True).copy()
    gdf["_i"] = gdf.index

    buffered = gdf.copy()
    buffered["geometry"] = buffered.geometry.buffer(buffer_m)

    joined = gpd.sjoin(
        buffered[["geometry", "_i"]],
        gdf[["geometry", "_i"]],
        how="left",
        predicate="intersects",
    )
    joined = joined[joined["_i_left"] != joined["_i_right"]]
    return joined.groupby("_i_left")["_i_right"].apply(list).to_dict()


def impute_zero_pop_rates(
    gdf: gpd.GeoDataFrame,
    rate_col: str = "reporting_rate_per_1000",
) -> gpd.GeoDataFrame:
    """
    Fill NaN reporting rates (zero-pop CBGs) with the median rate of
    spatially contiguous neighbors.  Falls back to the global median for
    isolated zero-pop clusters (rare).
    """
    gdf = gdf.copy().reset_index(drop=True)
    nan_mask = gdf[rate_col].isna()

    if nan_mask.sum() == 0:
        return gdf

    neighbor_map  = _build_neighbor_map(gdf)
    global_median = gdf[rate_col].median()
    imputed = fallback = 0

    for idx in gdf.index[nan_mask]:
        nbrs = neighbor_map.get(idx, [])
        nbr_rates = gdf.loc[nbrs, rate_col].dropna() if nbrs else pd.Series(dtype=float)
        if len(nbr_rates) > 0:
            gdf.loc[idx, rate_col] = nbr_rates.median()
            imputed += 1
        else:
            gdf.loc[idx, rate_col] = global_median
            fallback += 1

    print(f"  Imputed {imputed} CBGs from neighbors, "
          f"{fallback} fallback to global median; "
          f"remaining NaNs: {gdf[rate_col].isna().sum()}")
    return gdf


# ---------------------------------------------------------------------------
# 6. Merge 4-year ACS into cross-sectional averages
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
        pop_avg, med_inc_avg, med_age_avg, lep_rate_avg, edu_hs_avg
    """
    merge_keys   = ["state", "county", "tract", "block group"]
    drop_cols    = ["geometry", "GEOID", "NAME"]
    covariates   = ["pop", "med_inc", "med_age", "lep_rate", "edu_hs"]

    def _suffix(df, yr):
        yr2 = str(yr)[-2:]
        df  = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")
        return df.rename(columns={c: f"{c}_{yr2}" for c in covariates if c in df.columns})

    # Start from 2021 (carries geometry for the base merge)
    base = _suffix(acs_gdfs[2021], 2021)
    for yr in [2022, 2023, 2024]:
        base = base.merge(_suffix(acs_gdfs[yr], yr), on=merge_keys, how="inner")

    for cov in covariates:
        cols = [f"{cov}_{str(yr)[-2:]}" for yr in YEARS]
        base[f"{cov}_avg"] = base[cols].mean(axis=1)

    avg_cols = merge_keys + [f"{c}_avg" for c in covariates]
    panel = base[avg_cols].copy()

    # Fill any remaining NaNs in averages with column mean
    for col in [f"{c}_avg" for c in covariates]:
        panel[col] = panel[col].fillna(panel[col].mean())

    # lep_rate_avg expressed as percentage (0-100)
    panel["lep_rate_avg"] = panel["lep_rate_avg"] * 100

    return panel


# ---------------------------------------------------------------------------
# 7. Add area and population density
# ---------------------------------------------------------------------------

def add_area_and_density(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Compute CBG area (m² and km²) and population density (residents / km²).
    Requires EPSG:26918 (metres).
    """
    gdf = gdf.copy()
    gdf["area_m2"]      = gdf.geometry.area
    gdf["area_km2"]     = gdf["area_m2"] / 1_000_000
    gdf["pop_density"]  = gdf["pop_avg"] / gdf["area_km2"]
    return gdf


# ---------------------------------------------------------------------------
# 8. Average 311 reporting rate across 4 years
# ---------------------------------------------------------------------------

def merge_avg_reporting_rate(
    base_gdf: gpd.GeoDataFrame,
    rate_gdfs: dict,
) -> gpd.GeoDataFrame:
    """
    Build a 4-year average reporting rate per CBG and merge it onto the
    base GeoDataFrame.

    Parameters
    ----------
    base_gdf  : GeoDataFrame with GEOID (the merged ACS panel)
    rate_gdfs : dict  {year (int): GeoDataFrame with GEOID + reporting_rate_per_1000}
    """
    rate_df = rate_gdfs[2021][["GEOID", "reporting_rate_per_1000"]].rename(
        columns={"reporting_rate_per_1000": "rate_21"}
    )
    for yr, suffix in [(2022, "rate_22"), (2023, "rate_23"), (2024, "rate_24")]:
        tmp = rate_gdfs[yr][["GEOID", "reporting_rate_per_1000"]].rename(
            columns={"reporting_rate_per_1000": suffix}
        )
        rate_df = rate_df.merge(tmp, on="GEOID", how="outer")

    rate_df["reporting_rate_avg"] = rate_df[["rate_21", "rate_22", "rate_23", "rate_24"]].mean(axis=1)

    return base_gdf.merge(rate_df[["GEOID", "reporting_rate_avg"]], on="GEOID", how="left")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("ACS + 311 Reporting Rate — Replication Pipeline")
    print("=" * 60)

    # --- 1. Load ACS raw data ---
    print("\n[1/8] Loading ACS data...")
    raw_acs = {}
    for yr in YEARS:
        print(f"  Year {yr}:")
        raw_acs[yr] = clean_acs(load_acs_year(yr))

    # --- 2. Attach CBG geometry ---
    print("\n[2/8] Attaching CBG geometry (TIGER 2023)...")
    pa_shp      = gpd.read_file(os.path.join(TIGER_DIR, "tl_2023_42_bg.shp"))
    philly_shp  = pa_shp[pa_shp["COUNTYFP"] == "101"].copy()

    acs_geo = {}
    for yr in YEARS:
        acs_geo[yr] = attach_geometry(raw_acs[yr], philly_shp)
        print(f"  acs_{str(yr)[-2:]} with geometry: {len(acs_geo[yr])} CBGs")

    # --- 3. Load & clean 311 data ---
    print("\n[3/8] Loading 311 service request data...")
    pc311 = {}
    for yr in YEARS:
        pc311[yr] = load_311_year(yr)

    # --- 4. Spatial join: 311 reports → CBG ---
    print("\n[4/8] Spatial join: assigning 311 reports to CBGs...")
    rate_gdfs = {}
    for yr in YEARS:
        print(f"  Year {yr}:")
        rate_gdfs[yr] = compute_reporting_rate(acs_geo[yr], pc311[yr], yr)

    # --- 5. Impute zero-population CBG rates ---
    print("\n[5/8] Imputing reporting rates for zero-population CBGs...")
    for yr in YEARS:
        print(f"  Year {yr}:")
        rate_gdfs[yr] = impute_zero_pop_rates(rate_gdfs[yr])

    # --- 6. Merge 4-year ACS averages ---
    print("\n[6/8] Merging 4-year ACS into cross-sectional averages...")
    acs_4year = merge_4year_acs(acs_geo)
    print(f"  Panel shape: {acs_4year.shape}")

    # Attach geometry from 2021 base (uses TIGER 2023 polygons for all years)
    acs_4year_geo = philly_shp.to_crs("EPSG:26918").merge(
        acs_4year,
        left_on=["TRACTCE", "BLKGRPCE"],
        right_on=["tract", "block group"],
        how="inner",
    )
    acs_4year_geo = acs_4year_geo[
        ["geometry", "state", "county", "tract", "block group",
         "pop_avg", "med_inc_avg", "med_age_avg", "lep_rate_avg", "edu_hs_avg"]
    ].copy()

    # Build GEOID on the merged panel
    acs_4year_geo["GEOID"] = (
        acs_4year_geo["state"].astype(str).str.zfill(2)
        + acs_4year_geo["county"].astype(str).str.zfill(3)
        + acs_4year_geo["tract"].astype(str).str.zfill(6)
        + acs_4year_geo["block group"].astype(str)
    )

    # --- 7. Area and population density ---
    print("\n[7/8] Computing CBG area and population density...")
    acs_4year_geo = add_area_and_density(acs_4year_geo)

    # --- 8. Merge 4-year average reporting rate ---
    print("\n[8/8] Merging 4-year average 311 reporting rate...")
    final = merge_avg_reporting_rate(acs_4year_geo, rate_gdfs)
    print(f"  Final dataset: {len(final)} CBGs, {len(final.columns)} columns")
    print(f"  Columns: {list(final.columns)}")

    # --- Export ---
    out_path = os.path.join(REPL_DATA, "philly_cbg_acs_311.geojson")
    final.to_crs("EPSG:4326").to_file(out_path, driver="GeoJSON")
    print(f"\nSaved: {out_path}")
    print("Done.")


if __name__ == "__main__":
    main()
