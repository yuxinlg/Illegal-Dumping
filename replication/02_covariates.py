"""
02_covariates.py
================
Build the CBG-level covariate panel used as spatial background-intensity
inputs for the Cox-Hawkes model.

Run order
---------
  1. 02a_acs_prep.py               ← produces output/philly_cbg_acs.geojson
  2. 02b_network_connectivity.py  ← produces output/philadelphia_walk_network_connectivity.geojson
                                   (run ONCE; computationally expensive ~2-6 h)
  3. THIS SCRIPT (02_covariates.py)

Required input files (not tracked in git — see README.md for download links)
-----------------------------------------------------------------------------
  data/tiger/tl_2023_42_bg/tl_2023_42_bg.shp          Census TIGER CBG 2023
  data/tiger/tl_2018_42_bg/tl_2018_42_bg.shp          Census TIGER CBG 2018  (foot traffic)
  data/Land_Care/phs_landcare.geojson                  PHS LandCare parcels
  data/Land_Use/Land_use.geojson                       City of Philadelphia land use
  data/Vacant_Block_Percent_Land.geojson               Vacant block percent land
  data/NDVI/HLS.S30.T18TVK.<date>.v2.0.B8A/B04.tif    HLS Sentinel-2 raster tiles (2021-2025)
  data/311/2021/public_cases_fc.shp  ...  2024/        Open Data Philly 311 records
  output/philly_cbg_acs.geojson                        ACS panel (from 02a_acs_prep.py)
  output/philadelphia_walk_network_connectivity.geojson Walk network edge metrics
                                                        (from 02b_network_connectivity.py,
                                                         OR upload the pre-computed file)

  *** FOOT TRAFFIC [PROPRIETARY — see Section 6 below] ***
  data/foot_traffic/*.csv.gz   SafeGraph Neighborhood Patterns monthly files
                               (licensed via Dewey Data; users must obtain access
                                independently before running Section 6)

Output
------
  output/cov_cbg.geojson   CBG-level covariate GeoDataFrame (EPSG 26918)

Covariate inventory
-------------------
  From ACS (2021-2024 5-yr averages):
    pop_density       population per km²
    med_inc_avg       median household income (USD)
    lep_density       limited-English-proficient persons per km²
    edu_hs_avg        share of population 25+ with ≥ high-school diploma

  Land care:
    landcare_area     total PHS LandCare area within CBG (m²)

  Land use (area in m²; 'c_dig2' classification):
    RLD  RMD  RHD        residential low / medium / high density
    CC   CBP  CMR        commercial consumer / business / mixed residential
    I                    industrial
    Civic  Cult  AR      civic-institution / cultural / active recreation
    Park   Cemetery      park/open space / cemetery
    Water                water bodies
    VL   Other           vacant land / other
    use_res  use_com  use_civ   aggregated residential / commercial / civic area

  Vacant parcels:
    vac_area          total estimated vacant land area within CBG (m²)

  NDVI:
    ndvi_mean_4yr     5-year average NDVI (2021-2025 summer scenes; column kept as ndvi_mean_4yr for model compatibility)

  Foot traffic [requires SafeGraph access]:
    alloc_avg_d_cnt       avg monthly unique device count (2021-2025, 2023 CBG)
    alloc_avg_s_cnt       avg monthly stop count         (2021-2025, 2023 CBG)
    unique_device_ratio_aw  unique devices / stops (area-weighted avg)

  311 reporting rate [derived from foot traffic + 311 data]:
    reporting_rate    monthly 311 reports per RAW_STOP_COUNTS (avg 2021-2025)

  Walk-network connectivity (from 02b_network_connectivity.py):
    total_length_w        total walk-street length in CBG (m)
    betweenness_avg_w     length-weighted avg betweenness_1000 (1 km radius)
    pagerank_avg_w        mean PageRank of walk edges in CBG
"""

# ============================================================
# SECTION 0: Imports and global paths
# ============================================================
import os
import glob
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterstats import zonal_stats
from shapely.validation import make_valid
from shapely.geometry import Point
from tqdm import tqdm

warnings.filterwarnings("ignore")

REPL_DIR = os.path.dirname(os.path.abspath(__file__))   # replication/
DATA     = os.path.join(REPL_DIR, "data")               # replication/data/
OUT      = os.path.join(REPL_DIR, "output")             # replication/output/
os.makedirs(OUT, exist_ok=True)

print("=" * 60)
print("02_covariates.py — CBG covariate panel")
print("=" * 60)


# ============================================================
# CONFIGURATION  (edit these flags before running)
# ============================================================

# ── Walk-network connectivity ─────────────────────────────────
# Network metrics are always included.  Two ways to supply the file:
#   Option A (recommended): place the pre-computed file at
#     output/philadelphia_walk_network_connectivity.geojson
#     → this script loads it directly, no extra work needed.
#   Option B: run 02b_network_connectivity.py first (~2-6 h on a
#     laptop); it writes the file to the same output/ path.
# If the file is absent for any reason, network columns are skipped.
NETWORK_PATH = os.path.join(OUT, "philadelphia_walk_network_connectivity.geojson")

# ── Foot traffic (SafeGraph Neighborhood Patterns) ────────────
# True  → process raw .csv.gz files in data/foot_traffic/
#          Requires a licensed Dewey Data / SafeGraph subscription.
#          Set False for public replication packages.
# False → skip foot-traffic processing; columns will be absent
#         from cov_cbg.geojson.
USE_FOOT_TRAFFIC = True


# ============================================================
# SECTION 1: Base CBG geometry
# Load 2023 Philadelphia Census block-group boundaries (TIGER/Line).
# All spatial operations are in EPSG:26918 (NAD83 UTM 18N, metres).
# ============================================================
print("\n[1] Loading base CBG geometry (TIGER 2023) ...")

tiger_2023_path = os.path.join(DATA, "tiger", "tl_2023_42_bg", "tl_2023_42_bg.shp")
pa_2023  = gpd.read_file(tiger_2023_path)
philly_2023   = pa_2023[pa_2023["COUNTYFP"] == "101"].copy()
philly_26918  = philly_2023.to_crs(epsg=26918)

print(f"  {len(philly_26918)} Philadelphia CBGs loaded.")


# ============================================================
# SECTION 2: ACS covariates
# Loaded from the pre-processed panel produced by 02a_acs_prep.py.
# Make sure to run 02a_acs_prep.py first.
# Output file: output/philly_cbg_acs.geojson
#
# Columns used:
#   GEOID, pop_density, med_inc_avg, lep_density, edu_hs_avg
# ============================================================
print("\n[2] Loading ACS panel (output of 02a_acs_prep.py) ...")

acs_path = os.path.join(OUT, "philly_cbg_acs.geojson")   # output of 02a_acs_prep.py
acs_21_24 = gpd.read_file(acs_path)
acs_21_24["GEOID"] = acs_21_24["GEOID"].astype(str)

acs_cols = ["GEOID", "pop_density", "med_inc_avg", "lep_density", "edu_hs_avg"]
# keep only available ACS columns (robustness in case column names differ slightly)
acs_cols = [c for c in acs_cols if c in acs_21_24.columns]
acs_base = acs_21_24[acs_cols + ["geometry"]].copy()

print(f"  {len(acs_base)} CBGs, columns: {acs_cols}")


# ============================================================
# SECTION 3: Land Care
# PHS LandCare parcel layer (2020–2025).
# Metric: total LandCare area (m²) within each CBG.
#
# Source:  https://data-phl.opendata.arcgis.com/datasets/phl::phs-landcare
# Downloaded: July 17, 2025
# ============================================================
print("\n[3] Processing Land Care area by CBG ...")

landcare_path = os.path.join(DATA, "Land_Care", "phs_landcare.geojson")
land_care = gpd.read_file(landcare_path)

# filter to study years (2020 included for slight buffer)
land_care["year"] = pd.to_numeric(land_care["year"], errors="coerce")
land_care_filt = land_care[
    (land_care["year"] >= 2020) & (land_care["year"] <= 2025)
].copy()
land_care_filt = land_care_filt.to_crs(epsg=26918)

# overlay: clip LandCare polygons to CBG boundaries
lc_overlay = gpd.overlay(
    land_care_filt, philly_26918[["GEOID", "geometry"]], how="intersection"
)
lc_overlay["landcare_area"] = lc_overlay.geometry.area

# aggregate total LandCare area per CBG
land_care_cbg = lc_overlay.groupby("GEOID")["landcare_area"].sum().reset_index()
print(f"  LandCare area computed for {len(land_care_cbg)} CBGs.")


# ============================================================
# SECTION 4: Land Use
# Current land-use classification from City of Philadelphia.
# Metric: area (m²) per land-use category (c_dig2) within each CBG.
#         Aggregated into residential / commercial / civic composites.
#
# Source:  https://opendataphilly.org/datasets/land-use/
# Downloaded: March 28, 2025
# ============================================================
print("\n[4] Processing Land Use composition by CBG ...")

landuse_path = os.path.join(DATA, "Land_Use", "Land_use.geojson")
land_use = gpd.read_file(landuse_path)
land_use_26918 = land_use.to_crs(epsg=26918)

# fix invalid geometries before overlay
land_use_26918["geometry"] = land_use_26918["geometry"].buffer(0)

# clip land-use polygons to CBG boundaries
lu_overlay = gpd.overlay(
    land_use_26918, philly_26918[["GEOID", "geometry"]], how="intersection"
)
lu_overlay["intersect_area"] = lu_overlay.geometry.area

# group by CBG and land-use category code
cbg_lu_area = (
    lu_overlay.groupby(["GEOID", "c_dig2"])["intersect_area"].sum().reset_index()
)

# pivot to wide format: one column per category
cbg_lu_wide = (
    cbg_lu_area.pivot(index="GEOID", columns="c_dig2", values="intersect_area")
    .fillna(0)
    .reset_index()
)

# rename numeric codes to readable labels
land_use_dict = {
    11: "RLD",  12: "RMD",   13: "RHD",
    21: "CC",   22: "CBP",   23: "CMR",
    31: "I",    41: "Civic", 51: "T",
    52: "GW",   61: "Cult",  62: "AR",
    71: "Park", 72: "Cemetery",
    81: "Water", 91: "VL",   92: "Other",
}
cbg_lu_wide = cbg_lu_wide.rename(columns=land_use_dict)

# composite land-use areas
resid_cols = [c for c in ["RLD", "RMD", "RHD"]       if c in cbg_lu_wide.columns]
com_cols   = [c for c in ["CC",  "CBP", "CMR"]        if c in cbg_lu_wide.columns]
civ_cols   = [c for c in ["Civic","Cult","AR","Park","Cemetery"] if c in cbg_lu_wide.columns]
cbg_lu_wide["use_res"] = cbg_lu_wide[resid_cols].sum(axis=1)
cbg_lu_wide["use_com"] = cbg_lu_wide[com_cols].sum(axis=1)
cbg_lu_wide["use_civ"] = cbg_lu_wide[civ_cols].sum(axis=1)

print(f"  Land-use composition computed for {len(cbg_lu_wide)} CBGs.")


# ============================================================
# SECTION 5: Vacant Block Percent Land
# Block-level vacancy estimates from the City of Philadelphia.
# Metric: estimated total vacant area (m²) within each CBG.
#         vac_area = LANDVACPERCENTAGE × parcel_area, summed over parcels
#         that spatially fall within each CBG.
#
# Source:  https://data-phl.opendata.arcgis.com/datasets/phl::vacant-block-percent-land
# Downloaded: June 11, 2025
# ============================================================
print("\n[5] Processing Vacant Block area by CBG ...")

vacant_path = os.path.join(DATA, "Vacant_Block_Percent_Land.geojson")
vacant_block = gpd.read_file(vacant_path)
vacant_block = vacant_block.to_crs(philly_26918.crs)

# compute vacant area for each parcel
vacant_block["parcel_area"] = vacant_block.geometry.area
vacant_block["vac_area"]    = (
    vacant_block["LANDVACPERCENTAGE"] * vacant_block["parcel_area"]
)

# spatial join: assign each parcel to its CBG, then sum
vac_joined = gpd.sjoin(
    vacant_block, philly_26918[["GEOID", "geometry"]],
    how="inner", predicate="intersects"
)
cbg_vac_area = vac_joined.groupby("GEOID")["vac_area"].sum().reset_index()
print(f"  Vacant area computed for {len(cbg_vac_area)} CBGs.")


# ============================================================
# SECTION 6: NDVI
# 5-year average Normalized Difference Vegetation Index
# from HLS Sentinel-2 (HLS.S30) summer scenes.
#
# Band files (place in data/NDVI/):
#   B8A (NIR) and B04 (Red) for 2021-06-29, 2022-06-29,
#   2023-06-19, 2024-06-28, 2025-07-05.
#
# Source:  NASA Earthdata  https://www.earthdata.nasa.gov/data/catalog/lpcloud-hlss30-2.0
# Downloaded: July 15, 2025 (2021-2024); June 16, 2026 (2025)
# ============================================================
print("\n[6] Computing NDVI by CBG (5-year average) ...")

# reproject CBGs to match HLS raster CRS (UTM 18N, EPSG 32618)
shp_32618 = philly_2023.to_crs("32618")

# dictionary: year-suffix → (B8A path, B04 path)
ndvi_files = {
    "21": (
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2021180T154809.v2.0.B8A.tif"),
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2021180T154809.v2.0.B04.tif"),
    ),
    "22": (
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2022180T154821.v2.0.B8A.tif"),
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2022180T154821.v2.0.B04.tif"),
    ),
    "23": (
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2023170T154819.v2.0.B8A.tif"),
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2023170T154819.v2.0.B04.tif"),
    ),
    "24": (
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2024180T154941.v2.0.B8A.tif"),
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2024180T154941.v2.0.B04.tif"),
    ),
    "25": (
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2025186T155711.v2.0.B8A.tif"),
        os.path.join(DATA, "NDVI", "HLS.S30.T18TVK.2025186T155711.v2.0.B04.tif"),
    ),
}

cbg_ndvi_dict = {}
for yr, (b8a_path, b04_path) in ndvi_files.items():
    # load NIR and Red bands
    with rasterio.open(b8a_path) as nir_src:
        nir      = nir_src.read(1).astype("float32")
        nir_meta = nir_src.meta.copy()
    with rasterio.open(b04_path) as red_src:
        red = red_src.read(1).astype("float32")

    # NDVI = (NIR - Red) / (NIR + Red)
    ndvi = np.where((nir + red) == 0, 0, (nir - red) / (nir + red))
    nir_meta.update(count=1, dtype="float32")

    # save NDVI raster temporarily for masking
    tmp_path = f"/tmp/_ndvi_{yr}.tif"
    with rasterio.open(tmp_path, "w", **nir_meta) as dst:
        dst.write(ndvi, 1)

    # mask to Philadelphia extent and rewrite
    with rasterio.open(tmp_path) as src:
        masked_ndvi, masked_transform = mask(src, shp_32618.geometry, crop=True)
        masked_ndvi = masked_ndvi[0]

    masked_meta = nir_meta.copy()
    masked_meta.update({
        "height":    masked_ndvi.shape[0],
        "width":     masked_ndvi.shape[1],
        "transform": masked_transform,
    })
    masked_path = f"/tmp/_masked_ndvi_{yr}.tif"
    with rasterio.open(masked_path, "w", **masked_meta) as dst:
        dst.write(masked_ndvi, 1)

    # zonal statistics: mean NDVI per CBG
    stats = zonal_stats(shp_32618, masked_path, stats=["mean"], geojson_out=True)
    cbg_ndvi = gpd.GeoDataFrame.from_features(stats).rename(
        columns={"mean": f"ndvi_mean_{yr}"}
    )
    cbg_ndvi_dict[yr] = cbg_ndvi
    print(f"  NDVI {yr} done.")

# merge yearly NDVI by GEOID
ndvi_merged = cbg_ndvi_dict["21"][["GEOID", "ndvi_mean_21"]].copy()
for yr in ["22", "23", "24", "25"]:
    ndvi_merged = ndvi_merged.merge(
        cbg_ndvi_dict[yr][["GEOID", f"ndvi_mean_{yr}"]], on="GEOID", how="left"
    )

# 2025 sanity check — print city-wide mean and range before averaging
ndvi_25_vals = ndvi_merged["ndvi_mean_25"].dropna()
print(f"  2025 NDVI sanity check — n={len(ndvi_25_vals)} CBGs | "
      f"mean={ndvi_25_vals.mean():.3f} | "
      f"min={ndvi_25_vals.min():.3f} | "
      f"max={ndvi_25_vals.max():.3f}")

# 5-year average (skip NaNs)
ndvi_merged["ndvi_mean_5yr"] = ndvi_merged[
    ["ndvi_mean_21", "ndvi_mean_22", "ndvi_mean_23", "ndvi_mean_24", "ndvi_mean_25"]
].mean(axis=1, skipna=True)

ndvi_cbg = ndvi_merged[["GEOID", "ndvi_mean_5yr"]].rename(
    columns={"ndvi_mean_5yr": "ndvi_mean_4yr"}   # keep column name for model compatibility
).copy()
print(f"  5-year average NDVI computed for {len(ndvi_cbg)} CBGs.")


# ============================================================
# SECTION 7: Foot Traffic  [PROPRIETARY DATA — SafeGraph]
# ============================================================
# NOTE ─────────────────────────────────────────────────────
# The foot-traffic data used in this study is SafeGraph
# Neighborhood Patterns (monthly CBG-level visit counts),
# accessed via the Dewey Data platform:
#   https://app.deweydata.io/data/advan/neighborhood-patterns
#
# This is a LICENSED product.  Researchers must obtain their
# own access credentials before running this section.
# The raw files should be downloaded as .csv.gz and placed in:
#   data/foot_traffic/
#
# Once you have the files, this script processes them to:
#   7a. Compute a 5-year (2021-2025) monthly average per CBG,
#       area-weighted to 2023 Census block-group boundaries.
#   7b. Compute a monthly 311-to-foot-traffic reporting rate
#       (311 non-dumping reports / RAW_STOP_COUNTS) and average
#       to a single cross-sectional value per CBG.
# ─────────────────────────────────────────────────────────────
print("\n[7] Processing Foot Traffic ...")
print(
    "  NOTE: SafeGraph Neighborhood Patterns is a licensed dataset.\n"
    "  Researchers must obtain their own access via Dewey Data:\n"
    "  https://app.deweydata.io/data/advan/neighborhood-patterns\n"
    "  Raw .csv.gz files must be placed in data/foot_traffic/\n"
)

FT_DIR = os.path.join(DATA, "foot_traffic")
ft_files = glob.glob(os.path.join(FT_DIR, "*.csv.gz"))

if not USE_FOOT_TRAFFIC:
    print(
        "  [SKIPPED] USE_FOOT_TRAFFIC = False.\n"
        "  Columns alloc_avg_d_cnt, alloc_avg_s_cnt, unique_device_ratio_aw,\n"
        "  and reporting_rate will be absent from cov_cbg.geojson."
    )
    result_avg       = None
    ft_all_month_agg = None
elif not ft_files:
    print(
        "  [SKIPPED] No foot-traffic files found in data/foot_traffic/.\n"
        "  Columns alloc_avg_d_cnt, alloc_avg_s_cnt, unique_device_ratio_aw,\n"
        "  and reporting_rate will be absent from cov_cbg.geojson."
    )
    result_avg       = None
    ft_all_month_agg = None

else:
    # ------------------------------------------------------------------
    # 7a. 5-year average foot traffic (area-weighted to 2023 CBGs)
    # ------------------------------------------------------------------
    print("  [7a] Loading foot-traffic files ...")

    # Tiger 2018: SafeGraph uses 2018 CBG boundaries
    tiger_2018_path = os.path.join(DATA, "tiger", "tl_2018_42_bg", "tl_2018_42_bg.shp")
    # → replication/data/tiger/tl_2018_42_bg/tl_2018_42_bg.shp
    pa_2018   = gpd.read_file(tiger_2018_path)
    philly_2018 = pa_2018[pa_2018["COUNTYFP"] == "101"].copy()
    philly_2018["GEOID"] = philly_2018["GEOID"].astype(str)

    # Column names changed between schema versions (2021-2024 vs 2025+).
    # Normalise to the 2021-2024 names so the rest of the code is uniform.
    _COL_RENAMES = {
        "DEVICE_COUNTS": "RAW_DEVICE_COUNTS",
        "STOP_COUNTS":   "RAW_STOP_COUNTS",
        "YEAR":          "Y",
        "MONTH":         "M",
    }

    df_list = []
    for fp in tqdm(ft_files, desc="  reading files"):
        df = pd.read_csv(fp)
        df = df.rename(columns={k: v for k, v in _COL_RENAMES.items() if k in df.columns})
        df = df.loc[df["AREA"].astype(str).str[:5] == "42101"]  # Philadelphia county
        df_list.append(df)

    df_all = pd.concat(df_list, ignore_index=True)

    # keep relevant columns (all years 2021-2025 included)
    df_visit = df_all[
        ["AREA", "DATE_RANGE_START", "DATE_RANGE_END", "RAW_DEVICE_COUNTS",
         "RAW_STOP_COUNTS", "Y"]
    ].copy()
    df_visit["AREA"] = df_visit["AREA"].astype(str)
    df_visit = df_visit[df_visit["AREA"].isin(philly_2018["GEOID"])]

    # 5-year aggregate: mean monthly counts per CBG (2018-boundary GEOIDs)
    df_avg = (
        df_visit
        .groupby("AREA", as_index=False)
        .agg(avg_d_cnt=("RAW_DEVICE_COUNTS", "mean"),
             avg_s_cnt=("RAW_STOP_COUNTS",   "mean"))
        .rename(columns={"AREA": "GEOID"})
    )
    df_avg["unique_device_ratio"] = (df_avg["avg_d_cnt"] / df_avg["avg_s_cnt"]).round(2)

    # merge onto 2018 CBG geometries
    df_avg_geo = philly_2018.merge(df_avg, on="GEOID", how="left")
    df_avg_geo = gpd.GeoDataFrame(df_avg_geo, geometry="geometry", crs=philly_2018.crs)

    # -------------------------------------------------------------------
    # Area-weighted spatial interpolation: 2018 CBG → 2023 CBG boundaries
    # Extensive quantities (counts): split proportionally by area overlap.
    # Intensive quantity (ratio):    area-weighted average.
    # -------------------------------------------------------------------
    print("  [7a] Area-weighted interpolation: 2018 → 2023 CBGs ...")

    A = df_avg_geo.copy().to_crs(26918)
    B = philly_26918.copy()
    A["geometry"] = A.geometry.apply(make_valid)
    B["geometry"] = B.geometry.apply(make_valid)

    extensive_cols = ["avg_d_cnt", "avg_s_cnt"]
    intensive_col  = "unique_device_ratio"

    A = A[extensive_cols + [intensive_col, "geometry"]].copy()
    B = B[["GEOID", "geometry"]].copy()

    A["aid"]    = np.arange(len(A))
    A["area_A"] = A.area
    B["area_B"] = B.area
    B_area = B[["GEOID", "area_B"]].copy()

    AB = gpd.overlay(
        A[["aid", "area_A"] + extensive_cols + [intensive_col, "geometry"]],
        B[["GEOID", "geometry"]],
        how="intersection",
        keep_geom_type=False,
    )
    AB = AB[AB.geometry.notna() & ~AB.geometry.is_empty].copy()
    AB["area_int"] = AB.area
    AB = AB[AB["area_int"] > 0].copy()
    AB["share_A"] = AB["area_int"] / AB["area_A"]

    # extensive: proportional split
    for col in extensive_cols:
        AB[f"alloc_{col}"] = AB[col] * AB["share_A"]
    ext_alloc = AB.groupby("GEOID", as_index=False)[
        [f"alloc_{c}" for c in extensive_cols]
    ].sum()

    # intensive: area-weighted average
    AB = AB.merge(B_area, on="GEOID", how="left")
    AB["weighted_rate"] = AB[intensive_col] * AB["area_int"]
    rate = AB.groupby("GEOID", as_index=False).agg(
        weighted_sum=("weighted_rate", "sum"),
        area_B=("area_B", "first"),
    )
    rate["unique_device_ratio_aw"] = rate["weighted_sum"] / rate["area_B"]

    result_avg = (
        B
        .merge(ext_alloc, on="GEOID", how="left")
        .merge(rate[["GEOID", "unique_device_ratio_aw"]], on="GEOID", how="left")
        .fillna({f"alloc_{c}": 0 for c in extensive_cols})
    )
    print(f"  Foot-traffic 5-yr avg computed for {len(result_avg)} CBGs (2023 boundaries).")

    # ------------------------------------------------------------------
    # 7b. Monthly 311 reporting rate
    # 311 reports / RAW_STOP_COUNTS per CBG-month, then cross-sectional avg.
    # ------------------------------------------------------------------
    print("  [7b] Computing monthly 311 reporting rate ...")

    # ── Monthly foot-traffic: area-weight each monthly CSV to 2023 CBGs ──
    print("  [7b] Interpolating monthly foot-traffic to 2023 CBGs ...")
    interpolated_results = []
    for df_monthly in tqdm(df_list, desc="  monthly FT interpolation"):
        cols_needed = [
            "AREA", "DATE_RANGE_START", "DATE_RANGE_END",
            "RAW_DEVICE_COUNTS", "RAW_STOP_COUNTS", "Y",
        ]
        df_m = df_monthly[cols_needed].copy()
        df_m["AREA"] = df_m["AREA"].astype(str)
        df_m = df_m[df_m["AREA"].isin(philly_2018["GEOID"])]

        if df_m.empty:
            continue

        # get month from this batch
        date_str   = df_m["DATE_RANGE_START"].dropna().iloc[0]
        month_value = int(pd.to_datetime(date_str).strftime("%m"))

        # merge onto 2018 CBG geometries
        A_m = philly_2018.merge(df_m, left_on="GEOID", right_on="AREA", how="inner")
        A_m = gpd.GeoDataFrame(A_m, geometry="geometry", crs=philly_2018.crs)
        A_m = A_m.to_crs(26918)
        B_m = philly_26918[["GEOID", "geometry"]].copy()

        A_m["geometry"] = A_m.geometry.apply(make_valid)
        B_m["geometry"] = B_m.geometry.apply(make_valid)

        ext_m = ["RAW_DEVICE_COUNTS", "RAW_STOP_COUNTS"]
        int_m = "Y"
        A_m = A_m[ext_m + [int_m, "geometry"]].copy()

        A_m["aid"]    = np.arange(len(A_m))
        A_m["area_A"] = A_m.area
        B_m["area_B"] = B_m.area
        Barea_m = B_m[["GEOID", "area_B"]].copy()

        AB_m = gpd.overlay(
            A_m[["aid", "area_A"] + ext_m + [int_m, "geometry"]],
            B_m[["GEOID", "geometry"]],
            how="intersection",
            keep_geom_type=False,
        )
        AB_m = AB_m[AB_m.geometry.notna() & ~AB_m.geometry.is_empty].copy()
        AB_m["area_int"] = AB_m.area
        AB_m = AB_m[AB_m["area_int"] > 0].copy()
        AB_m["share_A"] = AB_m["area_int"] / AB_m["area_A"]

        for col in ext_m:
            AB_m[f"alloc_{col}"] = AB_m[col] * AB_m["share_A"]
        ext_alloc_m = AB_m.groupby("GEOID", as_index=False)[
            [f"alloc_{c}" for c in ext_m]
        ].sum()

        AB_m = AB_m.merge(Barea_m, on="GEOID", how="left")
        AB_m["weighted_Y"] = AB_m[int_m] * AB_m["area_int"]
        rate_m = AB_m.groupby("GEOID", as_index=False).agg(
            weighted_sum=("weighted_Y", "sum"),
            area_B=("area_B", "first"),
        )
        rate_m["Y_aw"] = rate_m["weighted_sum"] / rate_m["area_B"]

        result_m = (
            B_m
            .merge(ext_alloc_m, on="GEOID", how="left")
            .merge(rate_m[["GEOID", "Y_aw"]], on="GEOID", how="left")
            .fillna({f"alloc_{c}": 0 for c in ext_m})
        )
        result_m["month"] = month_value
        interpolated_results.append(result_m)

    ft_all = gpd.GeoDataFrame(pd.concat(interpolated_results, ignore_index=True))
    ft_all = ft_all[(ft_all["Y_aw"] >= 2021) & (ft_all["Y_aw"] <= 2025)].copy()
    ft_all["Y_aw"] = ft_all["Y_aw"].astype(str).str[:4].astype(int)
    ft_all["year_month"] = ft_all["Y_aw"] * 100 + ft_all["month"]

    # ── Load 311 data 2021–2025, remove illegal-dumping records ──────────
    print("  [7b] Loading 311 records (2021–2025) ...")
    pc311_list = []
    for yr in [2021, 2022, 2023, 2024, 2025]:
        shp_path = os.path.join(DATA, "311", str(yr), "public_cases_fc.shp")
        # → replication/data/311/{yr}/public_cases_fc.shp  (already present)
        df_yr = gpd.read_file(shp_path)
        # exclude illegal dumping reports (those are the outcome variable)
        df_yr = df_yr[df_yr["service_na"] != "Illegal Dumping"].copy()
        df_yr = df_yr.to_crs(epsg=26918)
        pc311_list.append(df_yr)

    all_311 = gpd.GeoDataFrame(pd.concat(pc311_list, ignore_index=True))
    all_311["requested_"] = pd.to_datetime(all_311["requested_"])
    all_311["year_month"]  = (
        all_311["requested_"].dt.year * 100 + all_311["requested_"].dt.month
    )
    all_311 = all_311[all_311["geometry"].notna()][["year_month", "objectid", "geometry"]]

    # assign 311 reports to CBGs via spatial join
    all_311 = all_311.sjoin(philly_26918[["GEOID", "geometry"]], how="left", predicate="within")
    all_311 = all_311[["year_month", "objectid", "GEOID"]].dropna(subset=["GEOID"])

    # count 311 reports per CBG-month
    all_311_month = (
        all_311.groupby(["year_month", "GEOID"]).size().reset_index(name="count")
    )

    # merge monthly 311 counts into monthly foot-traffic
    ft_month_merged = ft_all[["GEOID", "year_month", "alloc_RAW_STOP_COUNTS"]].merge(
        all_311_month, on=["year_month", "GEOID"], how="left"
    )
    ft_month_merged["count"] = ft_month_merged["count"].fillna(0)
    ft_month_merged["reporting_rate"] = (
        ft_month_merged["count"] / ft_month_merged["alloc_RAW_STOP_COUNTS"]
    )

    # cross-sectional average reporting rate per CBG
    ft_all_month_agg = (
        ft_month_merged.groupby("GEOID")["reporting_rate"].mean().reset_index()
    )
    print(f"  Reporting rate computed for {len(ft_all_month_agg)} CBGs.")


# ============================================================
# SECTION 8: Walk-network connectivity  [PRE-COMPUTED]
# ============================================================
# The walking-network edge-level connectivity metrics are
# computationally expensive to produce.  They were computed
# ONCE by running 02b_network_connectivity.py and are loaded
# here from the resulting file.
#
# If you need to regenerate this file:
#   python replication/02b_network_connectivity.py
#
# The user can also upload the pre-computed file directly to:
#   output/philadelphia_walk_network_connectivity.geojson
#
# Metrics used (edge level, aggregated to CBG):
#   betweenness_1000  — local (1 km radius) betweenness centrality
#   pagerank          — global PageRank (damping=0.85, length-weighted)
#   length            — street segment length (m)
# ============================================================
print("\n[8] Loading walk-network connectivity ...")

if not os.path.exists(NETWORK_PATH):
    print(
        f"  [SKIPPED] File not found:\n  {NETWORK_PATH}\n"
        "  Either place the pre-computed file there, or run\n"
        "  02b_network_connectivity.py to generate it (~2-6 h).\n"
        "  Network metrics will be absent from cov_cbg.geojson."
    )
    cbg_walk_metrics = None

else:
    philly_walk = gpd.read_file(NETWORK_PATH)
    philly_walk = philly_walk.to_crs(epsg=26918)
    philly_walk["GEOID"] = None   # will be filled via spatial split below

    print(f"  Loaded {len(philly_walk):,} walk edges.")

    # ── Aggregate edge metrics to CBG level ────────────────────────────
    # Each edge is split at CBG boundaries; metrics are allocated
    # proportionally by length.  PageRank is averaged (not weighted).
    print("  Splitting edges to CBG boundaries and aggregating ...")

    cbg_26918 = philly_26918[["GEOID", "geometry"]].copy()

    all_segments = []
    for _, row in tqdm(
        philly_walk.iterrows(), total=len(philly_walk), desc="  splitting walk edges"
    ):
        edge_geom  = row.geometry
        orig_len   = row.get("length", edge_geom.length)
        if orig_len == 0:
            continue

        intersecting = cbg_26918[cbg_26918.intersects(edge_geom)]
        for _, cbg_row in intersecting.iterrows():
            seg_geom = edge_geom.intersection(cbg_row.geometry)
            if seg_geom.is_empty or seg_geom.geom_type not in (
                "LineString", "MultiLineString"
            ):
                continue

            seg_len  = (
                sum(ln.length for ln in seg_geom.geoms)
                if seg_geom.geom_type == "MultiLineString"
                else seg_geom.length
            )
            prop = seg_len / orig_len

            all_segments.append({
                "GEOID":            cbg_row["GEOID"],
                "seg_len":          seg_len,
                "betweenness_1000": row.get("betweenness_1000", 0) * prop,
                "pagerank":         row.get("pagerank", 0),    # not length-weighted
                "pagerank_count":   1,
            })

    seg_df = pd.DataFrame(all_segments)

    cbg_walk_agg = seg_df.groupby("GEOID", as_index=False).agg(
        total_length_w=("seg_len",          "sum"),
        betweenness_sum_w=("betweenness_1000", "sum"),
        pagerank_sum_w=("pagerank",         "sum"),
        pagerank_count_w=("pagerank_count", "sum"),
    )
    cbg_walk_agg["betweenness_avg_w"] = (
        cbg_walk_agg["betweenness_sum_w"] / cbg_walk_agg["total_length_w"]
    )
    cbg_walk_agg["pagerank_avg_w"] = (
        cbg_walk_agg["pagerank_sum_w"] / cbg_walk_agg["pagerank_count_w"]
    )
    cbg_walk_metrics = cbg_walk_agg[[
        "GEOID", "total_length_w", "betweenness_avg_w", "pagerank_avg_w"
    ]].copy()

    print(f"  Walk-network metrics aggregated for {len(cbg_walk_metrics)} CBGs.")


# ============================================================
# SECTION 9: Merge all covariates
# Build the final CBG-level covariate panel.
# Start from the ACS base (with geometry) and left-join each
# covariate table on GEOID.
# ============================================================
print("\n[9] Merging all covariates ...")

# start from ACS base (GeoDataFrame, EPSG 26918)
cov = acs_base.copy()
cov = gpd.GeoDataFrame(cov, geometry="geometry")

# merge GEOID from CBG shapefile (to ensure all CBGs present even if ACS has gaps)
all_geoids = philly_26918[["GEOID", "geometry"]].copy()
all_geoids["GEOID"] = all_geoids["GEOID"].astype(str)
cov["GEOID"] = cov["GEOID"].astype(str)

# re-base on the official CBG list so no CBG is accidentally dropped
cov = all_geoids.merge(
    cov.drop(columns="geometry"), on="GEOID", how="left"
)
cov = gpd.GeoDataFrame(cov, geometry="geometry", crs=philly_26918.crs)

# vacant area
cbg_vac_area["GEOID"] = cbg_vac_area["GEOID"].astype(str)
cov = cov.merge(cbg_vac_area, on="GEOID", how="left")

# NDVI
ndvi_cbg["GEOID"] = ndvi_cbg["GEOID"].astype(str)
cov = cov.merge(ndvi_cbg, on="GEOID", how="left")

# land use
cbg_lu_wide["GEOID"] = cbg_lu_wide["GEOID"].astype(str)
cov = cov.merge(cbg_lu_wide, on="GEOID", how="left")

# land care
land_care_cbg["GEOID"] = land_care_cbg["GEOID"].astype(str)
cov = cov.merge(land_care_cbg, on="GEOID", how="left").fillna({"landcare_area": 0})

# foot traffic (5-yr avg)  — only if data was available
if result_avg is not None:
    ft_cols = ["GEOID", "alloc_avg_d_cnt", "alloc_avg_s_cnt", "unique_device_ratio_aw"]
    ft_merge = result_avg[ft_cols].drop_duplicates(subset="GEOID").copy()
    ft_merge["GEOID"] = ft_merge["GEOID"].astype(str)
    cov = cov.merge(ft_merge, on="GEOID", how="left")

# 311 reporting rate — only if data was available
if ft_all_month_agg is not None:
    ft_all_month_agg["GEOID"] = ft_all_month_agg["GEOID"].astype(str)
    cov = cov.merge(ft_all_month_agg, on="GEOID", how="left")

# walk-network metrics — only if file was available
if cbg_walk_metrics is not None:
    cbg_walk_metrics["GEOID"] = cbg_walk_metrics["GEOID"].astype(str)
    cov = cov.merge(cbg_walk_metrics, on="GEOID", how="left")

print(f"  Final covariate table: {len(cov)} CBGs × {len(cov.columns)} columns.")
print(f"  Columns: {list(cov.columns)}")


# ============================================================
# SECTION 10: Export
# Output is written in EPSG:26918 (NAD83 UTM 18N, metres).
# ============================================================
print("\n[10] Exporting cov_cbg.geojson ...")

out_path = os.path.join(OUT, "cov_cbg.geojson")
cov.to_file(out_path, driver="GeoJSON")

print(f"  Saved → {out_path}")
print("\nDone.")
