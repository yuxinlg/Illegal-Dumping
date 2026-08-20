"""
03_analysis.py
==============
Fit a Cox-Hawkes spatiotemporal point process model (via BSTPP) for the
selected study area and produce all diagnostic figures.

Run order
---------
  01_data_cleaning.py  →  02a/02b/02c_covariates.py  →  THIS SCRIPT

Required inputs  (produced by earlier scripts or provided directly)
-------------------------------------------------------------------
  replication/output/illegal_dumping_full.geojson  from 01_data_cleaning.py
  replication/output/all_parks_gdf.geojson         from 01_data_cleaning.py
  replication/output/all_boxes_gdf.geojson         from 01_data_cleaning.py
  replication/output/cov_cbg.geojson               from 02_covariates.py

Optional inputs
---------------
  replication/output/litter_df_for_model.csv       from 01_data_cleaning.py
  replication/output/cleanup_df_for_model.csv      from 01_data_cleaning.py

Auxiliary data for spatial/context plots
-----------------------------------------
  data/tiger/tl_2023_42_bg/tl_2023_42_bg.shp
  data/Land_Care/phs_landcare.geojson
  data/Land_Use/Land_use.geojson
  data/Vacant_Block_Percent_Land.geojson
  data/City_Limits/City_Limits.shp

Outputs  (written to replication/output/figures/{prefix}/)
----------------------------------------------------------
  {prefix}_cov_coef.png              covariate posterior coefficients
  {prefix}_prop_excitation.png       proportion of Hawkes excitation
  {prefix}_trigger_posterior.png     trigger-kernel posterior
  {prefix}_spatial_trigger.png       3-panel spatial trigger + land-use
  {prefix}_uncertainty.png           mean & std-dev of f_s surface
  {prefix}_road_network.png          road network coloured by mean f_s
  {prefix}_trigger_time_decay.png    trigger time-decay curve
  {prefix}_freq_dist.png             Δt histogram between repeat requests
  {prefix}_trigger_time_dist.png     observed vs simulated Δt (trimmed)
  {prefix}_temporal_components.png   GP temporal intensity + histogram
  {prefix}_seasonal_components.png   GP seasonal intensity + histogram

    Where {prefix} is:
    USE_PARK = True   →  PARK_NAME as-is, e.g. "mifflin" or "fairmount_aoi"
    USE_PARK = False  →  "{lat:.4f}_{lon:.4f}", e.g. "39.9526_-75.1652"
"""

# =============================================================================
# SECTION 0 – Imports and global paths
# =============================================================================

import os
import warnings
import calendar
from datetime import datetime

os.environ["JAX_PLATFORM_NAME"] = "cpu"

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import contextily as ctx
from matplotlib.colors import ListedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import Point, box as shapely_box
from shapely.ops import unary_union

import numpyro.distributions as dist

import bstpp
import bstpp.trigger

# ── bstpp environment verification ────────────────────────────────────────────
# This analysis requires a CUSTOMIZED version of bstpp that adds the
# cox_hawkes_shared module (Cox_Hawkes_Shared, attach_litter_survey_args,
# attach_dumping_report_args).  The standard PyPI package (pip install bstpp)
# does NOT include this module and will fail here.
#
# To verify you are using the correct version, check the path printed below.
# It should point to your local customized bstpp, NOT the site-packages copy.
#
# How to set up the customized bstpp:
#   Option A — editable install from the local source:
#       pip install -e /path/to/your/customized/bstpp
#   Option B — add the parent directory of your customized bstpp to sys.path:
#       import sys; sys.path.insert(0, "/path/to/customized/")
#
print(f"  bstpp loaded from: {bstpp.__file__}")

_REQUIRED_CUSTOM = ["cox_hawkes_shared"]
_REQUIRED_SYMBOLS = {
    "bstpp.cox_hawkes_shared": ["Cox_Hawkes_Shared",
                                 "attach_litter_survey_args",
                                 "attach_dumping_report_args"],
}
import importlib as _il
for _mod_name, _symbols in _REQUIRED_SYMBOLS.items():
    try:
        _mod = _il.import_module(_mod_name)
    except ImportError:
        raise ImportError(
            f"\n\n{'='*60}\n"
            f"WRONG BSTPP VERSION\n"
            f"{'='*60}\n"
            f"Module '{_mod_name}' not found in the bstpp loaded from:\n"
            f"  {bstpp.__file__}\n\n"
            f"This script requires the CUSTOMIZED bstpp that includes\n"
            f"cox_hawkes_shared.  The standard PyPI package does not have it.\n"
            f"See the comment above this block for setup instructions.\n"
        )
    for _sym in _symbols:
        if not hasattr(_mod, _sym):
            raise AttributeError(
                f"\n\nCustom symbol '{_sym}' missing from {_mod_name}.\n"
                f"bstpp loaded from: {bstpp.__file__}\n"
                f"Make sure you are using the correct customized version.\n"
            )
print(f"  bstpp custom modules verified (cox_hawkes_shared OK).")
del _il, _REQUIRED_CUSTOM, _REQUIRED_SYMBOLS

import bstpp.cox_hawkes_shared
from bstpp.trigger import Temporal_Power_Law
from bstpp.cox_hawkes_shared import (
    Cox_Hawkes_Shared,
    attach_litter_survey_args,
    attach_dumping_report_args,
)

warnings.filterwarnings("ignore")

# ── Coordinate parser ─────────────────────────────────────────────────────────

def parse_coord(value) -> float:
    """
    Parse a latitude or longitude value that may be supplied in either:
      • Decimal degrees  (float or numeric string):   39.9526  /  -75.1652
      • DMS string with degree/minute/second symbols:
            "39°57'54.1\"N"   "75°9'54.7\"W"
            "39°57'54.1N"     "75°9'54.7W"
            "39 57 54.1 N"    "75 9 54.7 W"
        The hemisphere letter (N/S/E/W) is optional; a leading minus sign on
        the numeric part works for W/S as well.

    Returns
    -------
    float  — signed decimal degrees (negative for South / West)
    """
    if isinstance(value, (int, float)):
        return float(value)

    s = str(value).strip()

    # Try plain decimal first
    try:
        return float(s)
    except ValueError:
        pass

    # Extract hemisphere letter if present (N/S/E/W at start or end)
    import re
    hemi_match = re.search(r"([NSEWnsew])$", s.strip())
    hemisphere = ""
    if hemi_match:
        hemisphere = hemi_match.group(1).upper()
        s = s[: hemi_match.start()].strip()
    else:
        hemi_match = re.match(r"^([NSEWnsew])\s*", s)
        if hemi_match:
            hemisphere = hemi_match.group(1).upper()
            s = s[hemi_match.end():].strip()

    # Tokenise: split on °, ', ", whitespace
    tokens = re.split(r'[°\'"″′\s]+', s.strip(" °'\""))
    tokens = [t for t in tokens if t]
    if not tokens:
        raise ValueError(f"Cannot parse coordinate: {value!r}")

    degrees  = float(tokens[0]) if len(tokens) > 0 else 0.0
    minutes  = float(tokens[1]) if len(tokens) > 1 else 0.0
    seconds  = float(tokens[2]) if len(tokens) > 2 else 0.0

    decimal = abs(degrees) + minutes / 60.0 + seconds / 3600.0
    # Preserve sign of degrees (handles "-75°…" style)
    if degrees < 0:
        decimal = -decimal

    if hemisphere in ("S", "W"):
        decimal = -abs(decimal)
    elif hemisphere in ("N", "E"):
        decimal = abs(decimal)

    return decimal


REPL_DIR = os.path.dirname(os.path.abspath(__file__))   # replication/
DATA     = os.path.join(REPL_DIR, "data")               # replication/data/
OUT      = os.path.join(REPL_DIR, "output")             # replication/output/
FIG_BASE = os.path.join(OUT, "figures")                 # replication/output/figures/
# Per-run subfolder (output/figures/{prefix}/) is created in __main__ after
# the prefix is determined from PARK_NAME or the custom coordinate.

print("=" * 60)
print("03_analysis.py — Cox-Hawkes model fitting + figures")
print("=" * 60)


# =============================================================================
# SECTION 1 – Configuration  (edit these values to change the analysis)
# =============================================================================
#
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  FOUR MODES — choose ONE                                                │
# │                                                                         │
# │  MODE A · Named park  (USE_PARK = True,  USE_CITY = USE_DISTRICT=False) │
# │    Study area = buffered box around a pre-built park polygon.           │
# │    Available parks: cobbscreek, mifflin, tacony, fairmount,             │
# │                     fairmount_aoi, cobbs_aoi, tacony_aoi, cobbs_aoi_focus│
# │    → Set PARK_NAME (no "_box" suffix needed).                           │
# │                                                                         │
# │  MODE B · Custom coordinate  (USE_PARK = False, USE_CITY/DISTRICT=False)│
# │    Study area = square centred on any WGS84 lat/lon in Philadelphia.    │
# │    → Set CUSTOM_LAT, CUSTOM_LON, CUSTOM_SIDE.                          │
# │    Accepts decimal degrees or DMS strings, e.g. "39°57'54.1N".         │
# │                                                                         │
# │  MODE C · Whole city  (USE_CITY = True)                                 │
# │    Study area = Philadelphia city boundary (City_Limits.shp).           │
# │    Overrides USE_PARK and USE_DISTRICT.                                 │
# │                                                                         │
# │  MODE D · Planning district  (USE_DISTRICT = True, USE_CITY = False)    │
# │    Study area = the named district's own polygon, as-is — no buffer,    │
# │    no bounding box.  Overrides USE_PARK.                                │
# │    → Set DISTRICT_NAME (see Step 2C for available names).               │
# └─────────────────────────────────────────────────────────────────────────┘

# ── Step 1: choose mode ───────────────────────────────────────────────────────
USE_CITY     = False    # True = whole Philadelphia (Mode C) — overrides USE_PARK/USE_DISTRICT
USE_DISTRICT = True   # True = named planning district (Mode D) — overrides USE_PARK; ignored if USE_CITY=True
USE_PARK     = False    # True = named park (Mode A) | False = custom point (Mode B)

# ── Step 2A: named-park settings  (only used in Mode A) ──────────────────────
PARK_NAME = "cobbscreek"       # e.g. "mifflin", "cobbscreek", "tacony", "fairmount"
BOX_NAME  = f"{PARK_NAME}_box"   # resolved automatically — do not edit

# ── Step 2B: custom-coordinate settings  (only used in Mode B) ───────────────
# Either decimal degrees OR DMS strings are accepted, e.g.:
#   CUSTOM_LAT = 39.9526   or   "39°57'54.1N"   or   "39 57 54.1 N"
#   CUSTOM_LON = -75.1652  or   "75°9'54.7W"
CUSTOM_LAT  = "39 57 54.1N"
CUSTOM_LON  = "75 12 19.2W"
CUSTOM_SIDE =  3000.0       # side length of the study square in metres

# ── Step 2C: planning-district settings  (only used in Mode D) ───────────────
# Study area = the district's own polygon as-is (no buffer, no bounding box).
# Available districts (from Planning_Districts.geojson, column "dist_name"):
#   River Wards, North Delaware, Lower Far Northeast, Central,
#   University Southwest, Upper Northwest, Upper North, South, North,
#   Lower South, Lower Northeast, Central Northeast, West,
#   Upper Far Northeast, Lower Southwest, West Park, Lower North,
#   Lower Northwest
DISTRICT_NAME = "University Southwest"   # exact dist_name value, used as-is for the figure prefix

# ── Step 3: covariates ────────────────────────────────────────────────────────
# Columns to pull from cov_cbg.geojson as spatial background-intensity inputs.
# Any column listed here that is absent from cov_cbg.geojson will be silently
# dropped with a warning at runtime — safe to leave extras in the list.
COV_NAMES = [
    # ACS demographics (from 02a_acs_prep.py)
    "pop_density", "med_inc_avg", "lep_density", "edu_hs_avg",
    # environment
    "ndvi_mean_4yr",
    # land use (area in m² per CBG; from 02_covariates.py Section 4)
    "RLD", "RMD", "RHD", "use_com", "use_civ", "I", "T", "GW",
    # vacancy & green space
    "vac_area", "landcare_area",
    # foot traffic — proprietary SafeGraph data; remove if USE_FOOT_TRAFFIC=False
    "alloc_avg_d_cnt", "alloc_avg_s_cnt", "unique_device_ratio_aw", "reporting_rate",
    # walk-network connectivity (from 02b_network_connectivity.py / pre-computed file)
    "betweenness_avg_w", "pagerank_avg_w",
]

# ── Study years ───────────────────────────────────────────────────────────────
# List the calendar years to include.  Events outside these years are dropped.
# TOTAL_DAYS is derived automatically as (Jan 1 of last year+1) − (Jan 1 of
# first year), so the model horizon always matches the selected window exactly.
# Examples:
#   STUDY_YEARS = [2021, 2022, 2023, 2024, 2025]   # full 5-year window
#   STUDY_YEARS = [2021, 2022]                      # 2-year subset
#   STUDY_YEARS = [2023]                            # single year
# STUDY_YEARS     = [2021, 2022, 2023, 2024, 2025]
STUDY_YEARS     = [2021, 2022, 2023, 2024, 2025]
TOTAL_DAYS      = (
    pd.Timestamp(f"{max(STUDY_YEARS) + 1}-01-01")
    - pd.Timestamp(f"{min(STUDY_YEARS)}-01-01")
).days
OFFSET_SEASONAL = 0         # seasonal cycle shift in days (0 = Jan 1 = day 0)

print(f"  Study years: {STUDY_YEARS}  →  TOTAL_DAYS = {TOTAL_DAYS}")

# ── Covariate filtering ───────────────────────────────────────────────────────
# True  → events are restricted to CBGs that have full covariate coverage.
#         Requires cov_cbg.geojson (from 02_covariates.py).
# False → no covariate filter; use all events inside the study box.
FILTER_TO_COV = True

# ── Model window parameters ───────────────────────────────────────────────────
# These are passed to model.set_window() and directly shape the model's
# internal representation of space and time.
#
# TEMPORAL_WINDOW  (int)
#   Number of temporal grid cells used for the background intensity GP.
#   Larger values → finer temporal resolution, slower fitting.
#   Typical range: 50–200.  Default: 100.
#
# SPATIAL_WINDOW  (float)
#   Bandwidth of the spatial background GP (unitless, relative to study box).
#   Smaller → sharper spatial features; larger → smoother surface.
#   Typical range: 0.05–0.3.  Default: 0.1.
TEMPORAL_WINDOW = 0.8       # temporal grid cells for background GP
SPATIAL_WINDOW  = 0.001       # spatial bandwidth for background GP

# ── SVI optimiser parameters ──────────────────────────────────────────────────
SVI_LR        = 0.01        # learning rate for Adam optimiser
SVI_NUM_STEPS = 20_000      # total SVI gradient steps

# ── Road-network plot ─────────────────────────────────────────────────────────
# True  → downloads the OSMnx walk network live and colours edges by mean f_s.
#         Can take 30 s – several minutes depending on internet speed.
# False → skip this figure.
PLOT_ROAD_NETWORK = True

# ── Figure settings ───────────────────────────────────────────────────────────
FIG_DPI             = 300
TEMPORAL_START_YEAR = 2021  # Jan 1 of this year = day 0 on the temporal axis
SEASONAL_REF_YEAR   = 2021  # reference year used to draw month gridlines
TRIGGER_XLIM_MAX    = 200   # x-axis cap (days) for the trimmed trigger-time plot
N_SIM_TRIGGER       = 500   # Monte-Carlo draws for simulated Δt distribution

# ── File paths (no edits needed below this line) ──────────────────────────────
ILLEGAL_DUMPING_PATH = os.path.join(OUT,  "illegal_dumping_full.geojson")
ALL_PARKS_PATH       = os.path.join(OUT,  "all_parks_gdf.geojson")
ALL_BOXES_PATH       = os.path.join(OUT,  "all_boxes_gdf.geojson")
COV_CBG_PATH         = os.path.join(OUT,  "cov_cbg.geojson")
LITTER_DF_PATH       = os.path.join(OUT,  "litter_df_for_model.csv")
CLEANUP_DF_PATH      = os.path.join(OUT,  "cleanup_df_for_model.csv")
CITY_LIMITS_PATH     = os.path.join(DATA, "City_Limits", "City_Limits.shp")
PLANNING_DISTRICTS_PATH = os.path.join(DATA, "Planning_Districts.geojson")
TIGER_BG_PATH        = os.path.join(DATA, "tiger", "tl_2023_42_bg", "tl_2023_42_bg.shp")
LAND_CARE_PATH       = os.path.join(DATA, "Land_Care", "phs_landcare.geojson")
LAND_USE_PATH        = os.path.join(DATA, "Land_Use", "Land_use.geojson")
VACANT_BLOCK_PATH    = os.path.join(DATA, "Vacant_Block_Percent_Land.geojson")


# =============================================================================
# SECTION 2 – Functions
# =============================================================================

# ─── 2-A: Figure-name prefix ──────────────────────────────────────────────────

def make_fig_prefix(
    use_park: bool,
    park_name: str,
    lat: float,
    lon: float,
) -> str:
    """
    Return the filename prefix used for every figure in this run.

    Parameters
    ----------
    use_park : bool
        True  → prefix is derived from park_name (first underscore token).
        False → prefix encodes the lat/lon coordinates.
    park_name : str
        Box label such as "mifflin_box".
    lat, lon : float
        WGS84 coordinates (used when use_park=False).

    Returns
    -------
    str  — e.g. "mifflin", "fairmount_aoi"  or  "39.9526_-75.1652"
    """
    if use_park:
        # Strip only the trailing "_box" so AOI variants keep their full name:
        #   "mifflin_box"       → "mifflin"
        #   "fairmount_aoi_box" → "fairmount_aoi"
        return park_name.removesuffix("_box")
    return f"{lat:.4f}_{lon:.4f}"


# ─── 2-B: Load supporting shapefiles ─────────────────────────────────────────

def load_supporting_layers(
    tiger_bg_path: str,
    land_care_path: str,
    land_use_path: str,
    vacant_block_path: str,
    city_limits_path: str,
) -> dict:
    """
    Load auxiliary geographic layers used in the spatial-trigger and
    road-network figures.

    Parameters
    ----------
    tiger_bg_path : str    TIGER/Line 2023 block-group shapefile.
    land_care_path : str   PHS land-care GeoJSON.
    land_use_path : str    Philadelphia land-use GeoJSON.
    vacant_block_path : str  Vacant-block GeoJSON.
    city_limits_path : str Philadelphia city-limits shapefile.

    Returns
    -------
    dict with keys:
        philly_cbg    – GeoDataFrame, EPSG:26918  (or None)
        land_care     – GeoDataFrame, EPSG:26918, 2021–2025  (or None)
        land_use      – GeoDataFrame, EPSG:3857 with 'type' column  (or None)
        vacant_block  – GeoDataFrame, EPSG:26918  (or None)
        philly_border – GeoDataFrame, EPSG:26918  (or None)
        use_color     – dict {hex_color: land_use_label}
    """
    layers = {}

    # ── Census block groups ───────────────────────────────────────────────────
    if os.path.exists(tiger_bg_path):
        pa = gpd.read_file(tiger_bg_path)
        layers["philly_cbg"] = pa[pa["COUNTYFP"] == "101"].to_crs(26918)
        print(f"  CBG layer: {len(layers['philly_cbg'])} block groups")
    else:
        print(f"  [WARN] CBG shapefile not found: {tiger_bg_path}")
        layers["philly_cbg"] = None

    # ── Land-care ─────────────────────────────────────────────────────────────
    if os.path.exists(land_care_path):
        lc = gpd.read_file(land_care_path)
        lc["year"] = pd.to_numeric(lc["year"], errors="coerce")
        layers["land_care"] = (
            lc[(lc["year"] >= 2021) & (lc["year"] <= 2025)].to_crs(26918)
        )
        print(f"  Land-care layer: {len(layers['land_care'])} records (2021–2025)")
    else:
        print(f"  [WARN] Land-care file not found: {land_care_path}")
        layers["land_care"] = None

    # ── Land-use colour palette: {hex: label}  (mirrors cbg-park-litter.ipynb) ──
    # Key = hex color, value = display label.  Order here defines both the
    # ListedColormap order and the pd.Categorical category order, so colors
    # and labels always align in the legend.
    _use_lst = [
        "Res Low Den", "Res Med Den", "Res High Den", "Commercial",
        "Civic", "Industry", "Transportation", "Green ROW",
        "Water", "Vacant", "Other",
    ]
    _colors = [
        "#543f32", "#db4d6d", "#69821b", "#e6b422", "#5654a2",
        "#0095d9", "#a73836", "#a8bf93", "#164a84", "#b3ada0", "#d3cfd9",
    ]
    layers["use_color"] = dict(zip(_colors, _use_lst))

    # ── Land-use ──────────────────────────────────────────────────────────────
    if os.path.exists(land_use_path):
        lu = gpd.read_file(land_use_path)
        lu["type"] = pd.Series(dtype="object")
        lu.loc[lu["c_dig2"] == 11, "type"] = "Res Low Den"
        lu.loc[lu["c_dig2"] == 12, "type"] = "Res Med Den"
        lu.loc[lu["c_dig2"] == 13, "type"] = "Res High Den"
        lu.loc[lu["c_dig2"].isin([21, 22, 23]), "type"]         = "Commercial"
        lu.loc[lu["c_dig2"].isin([41, 61, 62, 71, 72]), "type"] = "Civic"
        lu.loc[lu["c_dig2"] == 31, "type"] = "Industry"
        lu.loc[lu["c_dig2"] == 51, "type"] = "Green ROW"
        lu.loc[lu["c_dig2"] == 81, "type"] = "Water"
        lu.loc[lu["c_dig2"] == 91, "type"] = "Vacant"
        lu.loc[lu["c_dig2"].isin([52, 92]), "type"] = "Other"
        lu["type"] = lu["type"].fillna("Other")
        # enforce category order so ListedColormap colors always match labels
        lu["type"] = pd.Categorical(
            lu["type"],
            categories=list(layers["use_color"].values()),  # values = labels
            ordered=True,
        )
        layers["land_use"] = lu.to_crs(3857)
        print(f"  Land-use layer: {len(lu)} parcels")
    else:
        print(f"  [WARN] Land-use file not found: {land_use_path}")
        layers["land_use"] = None

    # ── Vacant block ──────────────────────────────────────────────────────────
    if os.path.exists(vacant_block_path):
        layers["vacant_block"] = gpd.read_file(vacant_block_path).to_crs(26918)
        print(f"  Vacant-block layer: {len(layers['vacant_block'])} records")
    else:
        print(f"  [WARN] Vacant-block file not found: {vacant_block_path}")
        layers["vacant_block"] = None

    # ── City limits ───────────────────────────────────────────────────────────
    if os.path.exists(city_limits_path):
        layers["philly_border"] = gpd.read_file(city_limits_path).to_crs(26918)
        print(f"  City-limits layer loaded")
    else:
        print(f"  [WARN] City-limits file not found: {city_limits_path}")
        layers["philly_border"] = None

    return layers


# ─── 2-C: Resolve study-area box ─────────────────────────────────────────────

def resolve_study_box(
    use_park: bool,
    all_boxes_gdf: gpd.GeoDataFrame,
    park_name: str,
    lat: float,
    lon: float,
    side_meters: float,
    city_limits_path: str,
    use_city: bool = False,
    use_district: bool = False,
    district_name: str | None = None,
    planning_districts_path: str | None = None,
    district_id_col: str = "dist_name",
) -> gpd.GeoDataFrame:
    """
    Return the single-row GeoDataFrame (EPSG:26918) for the study-area box.

    Mode C (use_city=True):     uses the Philadelphia city-limits polygon.
    Mode D (use_district=True): uses the named planning district's own
                                 polygon as-is — no buffer, no bounding box.
    Mode A (use_park=True):     looks up the named box from all_boxes_gdf.
    Mode B (use_park=False):    builds a square centred on the WGS84 coordinate.

    Parameters
    ----------
    use_city : bool  If True, use city limits as study area (overrides use_park/use_district).
    use_district : bool  If True, use the named planning district (overrides use_park).
    use_park : bool
    all_boxes_gdf : GeoDataFrame
    park_name : str   PARKNAME value (Mode A).
    lat, lon : float  WGS84 (Mode B).
    side_meters : float  Side length in metres (Mode B).
    city_limits_path : str
    district_name : str | None   dist_name value (Mode D).
    planning_districts_path : str | None  Path to Planning_Districts.geojson (Mode D).
    district_id_col : str  Column in the planning-districts file holding the name (Mode D).

    Returns
    -------
    gpd.GeoDataFrame  (single row, EPSG:26918, column PARKNAME)
    """
    if use_city:
        if not os.path.exists(city_limits_path):
            raise FileNotFoundError(
                f"City-limits file not found: {city_limits_path}\n"
                "Needed for USE_CITY = True."
            )
        city_gdf = gpd.read_file(city_limits_path).to_crs(26918)
        city_union = unary_union(city_gdf.geometry)
        city_box_geom = shapely_box(*city_union.bounds)
        result = gpd.GeoDataFrame(
            {"PARKNAME": ["philadelphia"]},
            geometry=[city_box_geom],
            crs=26918,
        )
        print(f"  Study box: whole Philadelphia "
              f"({city_box_geom.area / 1e6:.1f} km²)")
        return result

    if use_district:
        if not planning_districts_path or not os.path.exists(planning_districts_path):
            raise FileNotFoundError(
                f"Planning-districts file not found: {planning_districts_path}\n"
                "Needed for USE_DISTRICT = True."
            )
        districts_gdf = gpd.read_file(planning_districts_path).to_crs(26918)
        district_row = districts_gdf[districts_gdf[district_id_col] == district_name]
        if district_row.empty:
            raise ValueError(
                f"District '{district_name}' not found in {planning_districts_path}. "
                f"Check DISTRICT_NAME in Section 1 — available districts: "
                f"{districts_gdf[district_id_col].tolist()}"
            )
        # The district's own polygon is used as-is: no buffer, no bounding box.
        district_geom = district_row.geometry.iloc[0]
        result = gpd.GeoDataFrame(
            {"PARKNAME": [district_name]}, geometry=[district_geom], crs=26918
        )
        print(f"  Study box: planning district '{district_name}' "
              f"({district_geom.area / 1e6:.1f} km², own polygon, no buffer)")
        return result

    if use_park:
        box_gdf = all_boxes_gdf[all_boxes_gdf["PARKNAME"] == park_name].copy()
        if box_gdf.empty:
            raise ValueError(
                f"Box '{park_name}' not found in all_boxes_gdf. "
                f"Check PARK_NAME in Section 1 — available boxes: "
                f"{all_boxes_gdf['PARKNAME'].tolist()}"
            )
        print(f"  Study box: {park_name}")
        return box_gdf.reset_index(drop=True)

    # ── Custom coordinate mode ─────────────────────────────────────────────────
    point_26918 = gpd.GeoDataFrame(
        geometry=[Point(lon, lat)], crs=4326
    ).to_crs(26918)
    px = point_26918.geometry.iloc[0].x
    py = point_26918.geometry.iloc[0].y

    # Optional city-limits check
    if os.path.exists(city_limits_path):
        philly_union = unary_union(
            gpd.read_file(city_limits_path).to_crs(4326).geometry
        )
        if not philly_union.contains(Point(lon, lat)):
            warnings.warn(
                f"Coordinate ({lat}, {lon}) appears to be outside the "
                "Philadelphia city boundary.",
                stacklevel=2,
            )

    half = side_meters / 2.0
    square_geom = shapely_box(px - half, py - half, px + half, py + half)
    label = f"{lat:.4f}_{lon:.4f}_box"
    print(f"  Study box: custom ({lat}, {lon}), side={side_meters} m → {label}")
    return gpd.GeoDataFrame(
        {"PARKNAME": [label]}, geometry=[square_geom], crs=26918
    )


# ─── 2-D: Filter events and build locs_s ─────────────────────────────────────

def filter_events_and_build_locs_s(
    illegal_dumping: gpd.GeoDataFrame,
    study_box: gpd.GeoDataFrame,
    cov_gdf: gpd.GeoDataFrame | None,
    filter_to_cov: bool,
    study_years: list | None = None,
    offset_seasonal: int = 0,
) -> tuple:
    """
    Spatially restrict illegal-dumping events to the study box, optionally
    further filter to CBGs with covariate coverage, then build the
    model-ready event table (locs_s).

    Parameters
    ----------
    illegal_dumping : GeoDataFrame  (EPSG:26918)
    study_box : GeoDataFrame  single-row, EPSG:26918
    cov_gdf : GeoDataFrame | None  CBG covariate panel
    filter_to_cov : bool  restrict to CBGs with covariate coverage
    study_years : list | None
        Calendar years to include, e.g. [2021, 2022].  None = all years.
    offset_seasonal : int  seasonal cycle shift in days

    Returns
    -------
    (park_points, locs_s, baseline_time)
    park_points   – raw event GeoDataFrame for diagnostic plots
    locs_s        – pd.DataFrame with columns X, Y, T, A
    baseline_time – pd.Timestamp (Jan 1 of first study year)
    """
    print("\n[D] Filtering events and building locs_s ...")

    park_name = study_box["PARKNAME"].iloc[0]

    park_points = illegal_dumping.sjoin(study_box, predicate="within")
    if "index_right" in park_points.columns:
        park_points = park_points.drop(columns=["index_right"])
    park_points = park_points[park_points["PARKNAME"] == park_name].copy()
    print(f"  Events inside study box: {len(park_points):,}")

    if filter_to_cov and cov_gdf is not None:
        park_points = park_points.sjoin(
            cov_gdf[["geometry"]], predicate="within", how="inner"
        )
        if "index_right" in park_points.columns:
            park_points = park_points.drop(columns=["index_right"])
        print(f"  Events after covariate-coverage filter: {len(park_points):,}")

    park_points = park_points.reset_index(drop=True)

    # Build locs_s (replicates 01_data_cleaning.build_locs_s)
    pts = park_points.copy()
    pts["start_time"] = pd.to_datetime(pts["start_time"], errors="coerce")
    pts = pts.dropna(subset=["start_time"])

    # ── Year filter ────────────────────────────────────────────────────────────
    if study_years is not None:
        before = len(pts)
        pts = pts[pts["start_time"].dt.year.isin(study_years)].copy()
        print(f"  Year filter {study_years}: {before:,} → {len(pts):,} events")
        if pts.empty:
            raise ValueError(
                f"No events remain after filtering to STUDY_YEARS={study_years}."
            )
        baseline_time = pd.Timestamp(f"{min(study_years)}-01-01")
    else:
        baseline_time = pd.Timestamp(f"{pts['start_time'].min().year}-01-01")

    time_days     = (pts["start_time"] - baseline_time).dt.total_seconds() / 86400.0
    seasonal      = (time_days + offset_seasonal) % 365

    locs_s = pd.DataFrame({
        "X": pts.geometry.centroid.x.astype(float).values,
        "Y": pts.geometry.centroid.y.astype(float).values,
        "T": time_days.astype(float).values,
        "A": seasonal.astype(float).values,
    })

    print(f"  locs_s shape: {locs_s.shape}")
    print(f"  Baseline time: {baseline_time.date()}")
    print(f"  T range: [{locs_s['T'].min():.1f}, {locs_s['T'].max():.1f}] days")
    return park_points, locs_s, baseline_time


# ─── 2-E: Validate covariate names ───────────────────────────────────────────

def validate_cov_names(
    cov_names: list,
    cov_gdf: gpd.GeoDataFrame,
) -> list:
    """
    Return the subset of cov_names that exist as columns in cov_gdf,
    printing a warning for any that are absent.

    Parameters
    ----------
    cov_names : list  requested covariate column names
    cov_gdf : GeoDataFrame

    Returns
    -------
    list  — validated column names
    """
    available = set(cov_gdf.columns)
    missing = [c for c in cov_names if c not in available]
    if missing:
        warnings.warn(
            f"COV_NAMES columns missing from cov_cbg.geojson and dropped: "
            f"{missing}",
            stacklevel=2,
        )
    active = [c for c in cov_names if c in available]
    print(f"  Covariates: {len(active)} of {len(cov_names)} active")
    return active


# ─── 2-F: Set up and fit the Cox-Hawkes model ────────────────────────────────

def setup_and_fit_model(
    locs_s: pd.DataFrame,
    study_box: gpd.GeoDataFrame,
    total_days: int,
    cov_gdf: gpd.GeoDataFrame,
    active_cov_names: list,
    litter_df: pd.DataFrame | None,
    cleanup_df: pd.DataFrame | None,
    temporal_window: int = 100,
    spatial_window: float = 0.1,
    svi_lr: float = 0.01,
    svi_num_steps: int = 20_000,
) -> Cox_Hawkes_Shared:
    """
    Construct, configure, and fit a Cox_Hawkes_Shared model via SVI.

    Parameters
    ----------
    locs_s : pd.DataFrame
        Model-ready event table (X, Y, T, A).
    study_box : GeoDataFrame
        Single-row study-area box (EPSG:26918).
    total_days : int
        Study-window length in days.
    cov_gdf : GeoDataFrame
        CBG covariate panel.
    active_cov_names : list
        Validated covariate column names.
    litter_df : pd.DataFrame | None
        Litter-survey data; None skips attachment.
    cleanup_df : pd.DataFrame | None
        PPR cleanup data; None skips attachment.
    temporal_window : int
        Number of temporal grid cells for the background intensity GP.
    spatial_window : float
        Spatial bandwidth of the background GP (relative to study box).
    svi_lr : float
        Learning rate for the Adam optimiser.
    svi_num_steps : int
        Total SVI gradient steps.

    Returns
    -------
    Cox_Hawkes_Shared  — fitted model with posterior samples
    """
    print("\n[F] Setting up Cox-Hawkes model ...")
    model = Cox_Hawkes_Shared(
        locs_s,
        study_box,
        total_days,
        True,
        spatial_cov=cov_gdf,
        cov_names=active_cov_names,
        standardize_cov=True,
        a_0=dist.Normal(0, 0.5),
        alpha=dist.Beta(2, 4),
        temporal_trig=Temporal_Power_Law,
    )
    model.args["model"] = "cox_hawkes"

    if litter_df is not None:
        print("  Attaching litter-survey args ...")
        model.args = attach_litter_survey_args(model, litter_df, "score")

    if cleanup_df is not None:
        print("  Attaching cleanup-report args ...")
        model.args = attach_dumping_report_args(model, cleanup_df, "log_cleanup_cost")

    model.set_window(window=temporal_window, spatial_window=spatial_window)
    print(f"  Temporal window: {temporal_window},  spatial window: {spatial_window}")

    print(f"\n[F] Fitting model (lr={svi_lr}, num_steps={svi_num_steps}) ...")
    model.run_svi(lr=svi_lr, num_steps=svi_num_steps)
    print("  Fitting complete.")
    return model


# ─── 2-G: Check covariate fitting ────────────────────────────────────────────

def plot_covariate_check(
    model: Cox_Hawkes_Shared,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
) -> None:
    """
    Save a posterior covariate-coefficient forest plot
    (model.cov_weight_post_summary).

    Parameters
    ----------
    model : Cox_Hawkes_Shared
    prefix : str   figure filename prefix
    fig_out : str  output directory
    dpi : int
    """
    print("\n[G] Plotting covariate coefficient check ...")
    _show = plt.show
    plt.show = lambda: None
    model.cov_weight_post_summary(trace=True)
    fig = plt.gcf()
    out_path = os.path.join(fig_out, f"{prefix}_cov_coef.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.show = _show
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-H: Proportion of excitation ───────────────────────────────────────────

def plot_prop_excitation(
    model: Cox_Hawkes_Shared,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
) -> None:
    """Save model.plot_prop_excitation() figure."""
    print("\n[H] Plotting proportion excitation ...")
    _show = plt.show
    plt.show = lambda: None
    model.plot_prop_excitation()
    fig = plt.gcf()
    out_path = os.path.join(fig_out, f"{prefix}_prop_excitation.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.show = _show
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-I: Trigger posterior ──────────────────────────────────────────────────

def plot_trigger_posterior(
    model: Cox_Hawkes_Shared,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
) -> None:
    """Save model.plot_trigger_posterior(trace=False) figure."""
    print("\n[I] Plotting trigger posterior ...")
    _show = plt.show
    plt.show = lambda: None
    model.plot_trigger_posterior(trace=False)
    fig = plt.gcf()
    out_path = os.path.join(fig_out, f"{prefix}_trigger_posterior.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.show = _show
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-J: Spatial trigger 3-panel ────────────────────────────────────────────

def _fetch_major_roads(study_box: gpd.GeoDataFrame):
    """
    Fetch arterial-only OSM roads (motorway/trunk/primary/secondary, incl.
    link ramps) covering study_box, reprojected to study_box.crs.

    Returns None on failure (osmnx not installed, or a network/Overpass
    error) so callers can skip the overlay without failing the whole figure.
    """
    try:
        import osmnx as ox
    except ImportError:
        print("  [SKIP] osmnx is not installed. `pip install osmnx`")
        return None
    polygon_ll = (
        gpd.GeoSeries([study_box.geometry.iloc[0]], crs=study_box.crs)
        .to_crs(4326)
        .iloc[0]
    )
    custom_filter = (
        '["highway"~"motorway|motorway_link|trunk|trunk_link|'
        'primary|primary_link|secondary|secondary_link"]'
    )
    try:
        G = ox.graph_from_polygon(
            polygon_ll, custom_filter=custom_filter, simplify=True, retain_all=True
        )
        return ox.graph_to_gdfs(G, nodes=False).to_crs(study_box.crs)
    except Exception as e:
        print(f"  [SKIP] major-road fetch failed: {e}")
        return None


def plot_spatial_trigger_3panel(
    model: Cox_Hawkes_Shared,
    study_box: gpd.GeoDataFrame,
    all_parks_gdf: gpd.GeoDataFrame,
    layers: dict,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
    add_major_roads: bool = True,
) -> None:
    """
    Three-panel spatial trigger figure.

    Left/centre panels show the model's posterior spatial intensity surface
    (produced by model.plot_spatial); the right panel adds a land-use map
    with an OSM basemap for geographic context.

    Parameters
    ----------
    model : Cox_Hawkes_Shared
    study_box : GeoDataFrame  (EPSG:26918)
    all_parks_gdf : GeoDataFrame  (EPSG:26918)
    layers : dict  from load_supporting_layers()
    prefix, fig_out, dpi
    """
    print("\n[J] Plotting spatial trigger 3-panel ...")
    _show = plt.show
    plt.show = lambda: None

    model.plot_spatial(include_cov=True)

    fig  = plt.gcf()
    axs  = fig.get_axes()

    # ── Overlay park boundary and city outline on model panels ────────────────
    # all_parks_gdf here is pre-filtered to the current park polygon only.
    xlim = axs[0].get_xlim()
    ylim = axs[0].get_ylim()

    all_parks_gdf.plot(ax=axs[0], edgecolor="white", facecolor="none", linewidth=2)
    if layers["philly_border"] is not None:
        layers["philly_border"].plot(
            ax=axs[0], edgecolor="darkgray", facecolor="none", linewidth=2
        )
    axs[0].set_xlim(xlim)
    axs[0].set_ylim(ylim)
    axs[0].get_xaxis().set_visible(False)
    axs[0].get_yaxis().set_visible(False)

    if layers["philly_border"] is not None:
        layers["philly_border"].plot(
            ax=axs[1], edgecolor="darkgray", facecolor="none", linewidth=2
        )
    all_parks_gdf.plot(ax=axs[1], edgecolor="white", facecolor="none", linewidth=2)
    axs[1].set_xlim(xlim)
    axs[1].set_ylim(ylim)
    axs[1].get_xaxis().set_visible(False)
    axs[1].get_yaxis().set_visible(False)

    if add_major_roads:
        roads = _fetch_major_roads(study_box)
        if roads is not None and not roads.empty:
            roads.plot(ax=axs[1], color="white", linewidth=0.6, alpha=0.7, zorder=5)

    # ── Right panel: land-use + OSM basemap ───────────────────────────────────
    osm_ax  = fig.add_axes([1.0, 0.1, 0.3, 0.8])
    osm_box = study_box.to_crs(3857)

    if layers["philly_cbg"] is not None:
        layers["philly_cbg"].to_crs(3857).plot(
            ax=osm_ax, edgecolor="#3f312b", facecolor="none", linewidth=1, zorder=3
        )
    # Show the park/AOI boundary on the land-use panel; omit the buffer-zone box.
    all_parks_gdf.to_crs(3857).plot(
        ax=osm_ax, color="none", edgecolor="black", linewidth=2, zorder=2
    )

    if layers["land_use"] is not None:
        lu = layers["land_use"].copy()
        # use_color = {hex: label}; keys are colors in the same order as _use_lst
        cmap = ListedColormap(list(layers["use_color"].keys()))
        lu_plot = lu[lu["type"].notna()].copy()
        lu_plot.plot(
            ax=osm_ax,
            column="type",
            cmap=cmap,
            alpha=0.8,
            legend=True,
            legend_kwds={
                "loc": "center left",
                "bbox_to_anchor": (1, 0.5),
                "title": "Land Use",
            },
        )

    box_bounds = osm_box.total_bounds
    osm_ax.set_xlim(box_bounds[[0, 2]])
    osm_ax.set_ylim(box_bounds[[1, 3]])
    ctx.add_basemap(osm_ax, source=ctx.providers.CartoDB.Positron)

    scalebar = ScaleBar(
        dx=1, units="m", dimension="si-length", length_fraction=0.25,
        location="lower right", scale_loc="bottom",
        box_color="white", box_alpha=0.8, border_pad=0.6, color="black",
        font_properties={"size": 10},
    )
    osm_ax.add_artist(scalebar)
    osm_ax.set_title("Land Use")
    osm_ax.set_axis_off()

    out_path = os.path.join(fig_out, f"{prefix}_spatial_trigger.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.show = _show
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-K: Uncertainty surface (non-transparent) ──────────────────────────────

def plot_uncertainty(
    model: Cox_Hawkes_Shared,
    all_parks_gdf: gpd.GeoDataFrame,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
) -> None:
    """
    Two-panel figure: mean posterior f_s surface (left) and posterior
    standard-deviation (right).  Colour is fully opaque (non-transparent).

    Parameters
    ----------
    model : Cox_Hawkes_Shared
    all_parks_gdf : GeoDataFrame  (EPSG:26918)
        Pre-filtered to the current park/AOI polygon only (not all parks).
    prefix, fig_out, dpi
    """
    print("\n[K] Plotting uncertainty surface ...")

    f_xy_samples = model.samples["f_xy"]
    f_xy_mean    = np.mean(f_xy_samples, axis=0)
    f_xy_std     = np.std(f_xy_samples,  axis=0)

    bounds = model.A[["geometry"]].total_bounds
    minx, miny, maxx, maxy = bounds
    grid_size = int(np.sqrt(len(f_xy_mean)))
    print(f"  Grid size: {grid_size}×{grid_size}")

    f_xy_mean_grid = f_xy_mean.reshape(grid_size, grid_size)
    f_xy_std_grid  = f_xy_std.reshape(grid_size, grid_size)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    im1 = ax1.imshow(
        f_xy_mean_grid, origin="lower", cmap="viridis",
        extent=[minx, maxx, miny, maxy],
    )
    model.points.plot(ax=ax1, color="red", marker="x", alpha=0.1)
    all_parks_gdf.plot(ax=ax1, edgecolor="white", facecolor="none", linewidth=2)
    ax1.set_title("Mean Posterior $f_s$ With Events")
    ax1.set_xlim(minx, maxx)
    ax1.set_ylim(miny, maxy)
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    plt.colorbar(im1, ax=ax1, shrink=0.7)

    im2 = ax2.imshow(
        f_xy_std_grid, origin="lower", cmap="RdYlGn_r",
        extent=[minx, maxx, miny, maxy], vmin=0, vmax=0.5,
    )
    all_parks_gdf.plot(ax=ax2, edgecolor="gray", facecolor="none", linewidth=2)
    ax2.set_title("$f_s$ Uncertainty (Std Dev)")
    ax2.set_xlim(minx, maxx)
    ax2.set_ylim(miny, maxy)
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    plt.colorbar(im2, ax=ax2, shrink=0.7)

    plt.tight_layout()
    out_path = os.path.join(fig_out, f"{prefix}_uncertainty.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved → {out_path}")
    print(f"  Mean f_xy range: [{f_xy_mean.min():.3f}, {f_xy_mean.max():.3f}]")
    print(f"  Std  f_xy range: [{f_xy_std.min():.3f}, {f_xy_std.max():.3f}]")


# ─── 2-L: Colored road network ───────────────────────────────────────────────

def plot_road_network(
    model: Cox_Hawkes_Shared,
    study_box: gpd.GeoDataFrame,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
) -> None:
    """
    Download the OSM road network for the study area, colour each edge by the
    posterior mean spatial intensity f_s at that location, and save the figure.

    Requires osmnx and a live internet connection.

    Parameters
    ----------
    model : Cox_Hawkes_Shared
    study_box : GeoDataFrame  (EPSG:26918)
    prefix, fig_out, dpi
    """
    print("\n[L] Plotting colored road network ...")
    try:
        import osmnx as ox
    except ImportError:
        print("  [SKIP] osmnx is not installed. `pip install osmnx`")
        return

    # Build a 25×25 grid over the study box to interpolate f_s onto roads
    park_polygon = study_box.geometry.iloc[0]
    minx, miny, maxx, maxy = park_polygon.bounds
    x_coords = np.linspace(minx, maxx, 26)
    y_coords = np.linspace(miny, maxy, 26)
    grid_cells = [
        shapely_box(x_coords[i], y_coords[j], x_coords[i + 1], y_coords[j + 1])
        for i in range(25)
        for j in range(25)
    ]
    grid_gdf = gpd.GeoDataFrame({"geometry": grid_cells}, crs=study_box.crs)

    # Assign posterior mean f_s to each grid cell
    f_s_mean = np.array(model.samples["f_xy"]).mean(axis=0)
    model.comp_grid["mean_f_s"] = f_s_mean

    joined = gpd.sjoin(
        grid_gdf[["geometry"]],
        model.comp_grid[["geometry", "mean_f_s"]],
        how="left",
        predicate="intersects",
    )
    grid_gdf["mean_f_s"] = joined.groupby(joined.index)["mean_f_s"].mean()

    # Fetch road network (converts to WGS84 for OSMnx, then reprojects)
    park_polygon_ll = (
        gpd.GeoSeries([park_polygon], crs=study_box.crs)
        .to_crs(4326)
        .iloc[0]
    )
    G = ox.graph_from_polygon(park_polygon_ll, network_type="all")
    edges_gdf = ox.graph_to_gdfs(G, nodes=False).to_crs(grid_gdf.crs)

    edges_joined = gpd.sjoin(
        edges_gdf,
        grid_gdf[["geometry", "mean_f_s"]],
        how="left",
        predicate="intersects",
    )
    edges_gdf["f_s_value"] = edges_joined.groupby(edges_joined.index)["mean_f_s"].mean()

    fig, ax = plt.subplots(figsize=(9, 9))
    edges_gdf.plot(
        ax=ax, column="f_s_value", cmap="viridis", linewidth=2, legend=True,
        missing_kwds={"color": "lightgray", "linewidth": 1},
    )
    gpd.GeoDataFrame(
        [{"geometry": park_polygon}], crs=edges_gdf.crs
    ).plot(ax=ax, facecolor="none", edgecolor="white", linewidth=2)

    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect("equal")
    plt.title("Road Network colored by Mean Posterior $f_s$")
    plt.tight_layout()

    out_path = os.path.join(fig_out, f"{prefix}_road_network.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-M: Trigger time decay ──────────────────────────────────────────────────

def plot_trigger_time_decay(
    model: Cox_Hawkes_Shared,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
) -> None:
    """Save model.plot_trigger_time_decay() figure."""
    print("\n[M] Plotting trigger time decay ...")
    _show = plt.show
    plt.show = lambda: None
    model.plot_trigger_time_decay()
    fig = plt.gcf()
    out_path = os.path.join(fig_out, f"{prefix}_trigger_time_decay.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.show = _show
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-N: Frequency of Δt between repeat requests ─────────────────────────────

def plot_freq_time_diff(
    park_points: gpd.GeoDataFrame,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
) -> None:
    """
    Histogram of time differences (days) between consecutive illegal-dumping
    reports filed at the same address.

    Parameters
    ----------
    park_points : GeoDataFrame  events within the study box
    prefix, fig_out, dpi
    """
    print("\n[N] Plotting frequency of Δt between requests ...")
    df = park_points.copy()
    df["start_time"] = pd.to_datetime(df["start_time"])
    df.sort_values(by=["address", "start_time"], inplace=True)
    df["next_start_time"] = df.groupby("address")["start_time"].shift(-1)
    df["delta_t_days"] = (
        (df["next_start_time"] - df["start_time"]).dt.total_seconds() / 86400.0
    )
    df.dropna(subset=["delta_t_days"], inplace=True)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.hist(df["delta_t_days"], bins=50, edgecolor="black", alpha=0.7, color="gray")
    fig.suptitle("Frequency Distribution of Time Difference Between Requests")
    ax.set_xlim(-50, 1500)
    ax.set_xlabel("Time Difference (days)")
    ax.set_ylabel("Frequency")
    ax.grid(axis="y", alpha=0.75)

    out_path = os.path.join(fig_out, f"{prefix}_freq_dist.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-O: Trigger time distribution (trimmed, observed vs simulated) ──────────

def plot_trigger_time_distribution(
    model: Cox_Hawkes_Shared,
    locs_s: pd.DataFrame,
    park_points: gpd.GeoDataFrame,
    total_days: int,
    prefix: str,
    fig_out: str,
    dpi: int = 300,
    n_sim: int = 500,
    xlim_max: int = 200,
) -> None:
    """
    Observed vs simulated inter-event Δt distribution, trimmed to
    [0, xlim_max] days.

    Posterior mean trigger parameters (beta, gamma from Temporal_Power_Law)
    are extracted from model.samples.  The simulated distribution is the
    average of n_sim Monte-Carlo draws, each with as many samples as events
    in locs_s.

    Parameters
    ----------
    model : Cox_Hawkes_Shared
    locs_s : pd.DataFrame
    park_points : GeoDataFrame
    total_days : int
    prefix, fig_out, dpi
    n_sim : int    Monte-Carlo iterations
    xlim_max : int x-axis cap in days
    """
    print("\n[O] Plotting trigger time distribution (trimmed) ...")

    # ── Observed Δt ────────────────────────────────────────────────────────────
    df = park_points.copy()
    df["start_time"] = pd.to_datetime(df["start_time"])
    df.sort_values(by=["address", "start_time"], inplace=True)
    df["next_start_time"] = df.groupby("address")["start_time"].shift(-1)
    df["delta_t_days"] = (
        (df["next_start_time"] - df["start_time"]).dt.total_seconds() / 86400.0
    )
    df.dropna(subset=["delta_t_days"], inplace=True)

    min_val, max_val = df["delta_t_days"].min(), df["delta_t_days"].max()
    bins_obs = np.linspace(min_val, max_val, 1200)
    bin_centers_obs = (bins_obs[:-1] + bins_obs[1:]) / 2
    counts, _ = np.histogram(df["delta_t_days"], bins=bins_obs)

    # ── Simulated Δt ───────────────────────────────────────────────────────────
    par_names   = model.args["t_trig"].get_par_names()
    params_mean = {
        name: float(np.array(model.samples[name]).mean())
        for name in par_names
    }
    tpl   = Temporal_Power_Law(locs_s)
    scale = total_days / 50.0

    bin_width    = (max_val - min_val) / 1200.0
    max_time_sim = 365 * 4
    n_bins_sim   = int(max_time_sim / bin_width)
    bin_edges_sim  = np.linspace(0, max_time_sim, n_bins_sim + 1)
    bin_centers_sim = (bin_edges_sim[:-1] + bin_edges_sim[1:]) / 2

    histogram_matrix = []
    for i in range(n_sim):
        sim_samples = [
            float(tpl.simulate_trigger(params_mean)) * scale
            for _ in range(locs_s.shape[0])
        ]
        h, _ = np.histogram(sim_samples, bins=bin_edges_sim)
        histogram_matrix.append(h)
        if (i + 1) % 100 == 0:
            print(f"    Simulation {i + 1}/{n_sim}")

    col_avg = np.mean(np.array(histogram_matrix), axis=0)

    # ── Trimmed plot ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.plot(
        bin_centers_obs, counts,
        marker="o", linestyle="-", color="#69821b",
        linewidth=1.5, markersize=2, label="Observed", alpha=0.7,
    )
    ax.plot(
        bin_centers_sim, col_avg,
        marker="s", linestyle="-", color="#4d5aaf",
        linewidth=1.5, markersize=2, label="Simulated",
    )
    ax.set_title(
        "Frequency Distribution of Time Difference Between Requests", pad=30
    )
    ax.set_xlabel("Time Difference (days)")
    ax.set_ylabel("Frequency / Count")
    ax.set_xlim(-10, xlim_max)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    out_path = os.path.join(fig_out, f"{prefix}_trigger_time_dist.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-P: Temporal components ─────────────────────────────────────────────────

def plot_temporal_components(
    model: Cox_Hawkes_Shared,
    prefix: str,
    fig_out: str,
    start_year: int = 2021,
    label_every_n_months: int = 3,
    dpi: int = 300,
) -> None:
    """
    Two-row figure: GP temporal intensity (top) and event histogram (bottom),
    with a calendar-date x-axis derived from start_year.

    Parameters
    ----------
    model : Cox_Hawkes_Shared
    prefix, fig_out, dpi
    start_year : int  calendar year of model baseline (Jan 1 = day 0)
    label_every_n_months : int  tick-label spacing
    """
    print("\n[P] Plotting temporal components ...")
    start_date = datetime(start_year, 1, 1)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(15, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    model.plot_temporal_components(
        rescale=True, start_date=None, show_gp=True, show_hist=False, ax=ax_top
    )
    ax_top.tick_params(labelbottom=False)
    ax_top.set_xlabel("")

    model.plot_temporal_components(
        rescale=True, start_date=None, show_gp=False, show_hist=True, ax=ax_bot
    )

    # Calendar x-axis labels
    x0, x1 = ax_bot.get_xlim()
    end_date    = start_date + pd.Timedelta(days=int(np.ceil(x1)))
    month_starts = pd.date_range(start=start_date, end=end_date, freq="MS")
    labeled      = month_starts[::label_every_n_months]

    xticks  = [(d.to_pydatetime() - start_date).days for d in labeled]
    xlabels = [d.strftime("%Y-%m") for d in labeled]
    ax_bot.set_xticks(xticks)
    ax_bot.set_xticklabels(xlabels)
    ax_bot.set_xlabel("Date")

    # Month gridlines on both panels
    for ax in (ax_top, ax_bot):
        for d in labeled:
            ax.axvline(
                (d.to_pydatetime() - start_date).days,
                linestyle=":", linewidth=0.6, alpha=0.25,
            )

    plt.tight_layout()
    out_path = os.path.join(fig_out, f"{prefix}_temporal_components.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved → {out_path}")


# ─── 2-Q: Seasonal components ─────────────────────────────────────────────────

def plot_seasonal_components(
    model: Cox_Hawkes_Shared,
    prefix: str,
    fig_out: str,
    ref_year: int = 2021,
    dpi: int = 300,
) -> None:
    """
    Two-row figure: GP within-year seasonal intensity (top) and event
    histogram (bottom), with calendar-month gridlines.

    Parameters
    ----------
    model : Cox_Hawkes_Shared
    prefix, fig_out, dpi
    ref_year : int  reference year used to compute month-boundary positions
    """
    print("\n[Q] Plotting seasonal components ...")
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(15, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    model.plot_seasonal_components(
        rescale=True, ref_year=ref_year, show_gp=True, show_hist=False, ax=ax_top
    )
    ax_top.tick_params(labelbottom=False)
    ax_top.set_xlabel("")

    model.plot_seasonal_components(
        rescale=True, ref_year=ref_year, show_gp=False, show_hist=True, ax=ax_bot
    )

    # Month-boundary gridlines in seasonal day units
    offset = model.args.get("offset_seasonal", 0)
    S      = getattr(model, "S", 365)
    month_days = np.cumsum(
        [0] + [calendar.monthrange(ref_year, m)[1] for m in range(1, 12)]
    )
    month_boundaries = np.sort((month_days + offset) % S)

    for ax in (ax_top, ax_bot):
        for x_pos in month_boundaries:
            ax.axvline(x_pos, linestyle=":", linewidth=0.6, alpha=0.25)

    plt.tight_layout()
    out_path = os.path.join(fig_out, f"{prefix}_seasonal_components.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved → {out_path}")


# =============================================================================
# SECTION 3 – Main execution  (runs when script is called directly)
# =============================================================================

if __name__ == "__main__":

    # ── [1] Load upstream outputs ──────────────────────────────────────────────
    print("\n[1] Loading upstream outputs ...")
    illegal_dumping = gpd.read_file(ILLEGAL_DUMPING_PATH).to_crs(26918)
    all_parks_gdf   = gpd.read_file(ALL_PARKS_PATH).to_crs(26918)
    all_boxes_gdf   = gpd.read_file(ALL_BOXES_PATH).to_crs(26918)
    cov_gdf         = gpd.read_file(COV_CBG_PATH).to_crs(26918)
    print(f"  Illegal-dumping events:   {len(illegal_dumping):,}")
    print(f"  Covariate CBGs:           {len(cov_gdf):,}")
    print(f"  Park boxes available:     {all_boxes_gdf['PARKNAME'].tolist()}")

    # Parse CUSTOM_LAT/LON once — they may be decimal floats or DMS strings.
    _lat = parse_coord(CUSTOM_LAT)
    _lon = parse_coord(CUSTOM_LON)
    if not USE_PARK:
        print(f"  Custom coord (parsed): lat={_lat:.6f}, lon={_lon:.6f}")

    litter_df  = pd.read_csv(LITTER_DF_PATH)  if os.path.exists(LITTER_DF_PATH)  else None
    cleanup_df = pd.read_csv(CLEANUP_DF_PATH) if os.path.exists(CLEANUP_DF_PATH) else None
    if litter_df  is not None: print(f"  Litter-survey rows:  {len(litter_df):,}")
    if cleanup_df is not None: print(f"  Cleanup-report rows: {len(cleanup_df):,}")

    # ── [2] Load supporting spatial layers ─────────────────────────────────────
    print("\n[2] Loading supporting layers ...")
    layers = load_supporting_layers(
        tiger_bg_path     = TIGER_BG_PATH,
        land_care_path    = LAND_CARE_PATH,
        land_use_path     = LAND_USE_PATH,
        vacant_block_path = VACANT_BLOCK_PATH,
        city_limits_path  = CITY_LIMITS_PATH,
    )

    # ── [3] Resolve study-area box ─────────────────────────────────────────────
    print("\n[3] Resolving study-area box ...")
    study_box = resolve_study_box(
        use_city                = USE_CITY,
        use_park                = USE_PARK,
        all_boxes_gdf           = all_boxes_gdf,
        park_name               = BOX_NAME,
        lat                     = _lat,
        lon                     = _lon,
        side_meters             = CUSTOM_SIDE,
        city_limits_path        = CITY_LIMITS_PATH,
        use_district            = USE_DISTRICT,
        district_name           = DISTRICT_NAME,
        planning_districts_path = PLANNING_DISTRICTS_PATH,
    )
    active_park_name = study_box["PARKNAME"].iloc[0]
    # For city mode the prefix is always "philadelphia"; for district mode the
    # prefix is the district name as-is (already free of any "_box" suffix).
    if USE_CITY:
        prefix = "philadelphia"
    elif USE_DISTRICT:
        prefix = active_park_name
    else:
        prefix = make_fig_prefix(USE_PARK, active_park_name, _lat, _lon)
    print(f"  Figure prefix: '{prefix}'")

    # Park polygon for the active study area only — used in all plot overlays.
    # PARKNAME in all_parks_gdf uses the park key (e.g. "cobbs_aoi"), same as
    # active_park_name without its "_box" suffix.
    if USE_DISTRICT:
        # study_box already holds the district's own polygon (no buffer) —
        # use it directly as the overlay instead of looking it up in
        # all_parks_gdf, which only contains the named park polygons.
        current_park_gdf = study_box.copy()
    else:
        park_key = active_park_name.removesuffix("_box")
        current_park_gdf = all_parks_gdf[all_parks_gdf["PARKNAME"] == park_key].copy()
        if current_park_gdf.empty:
            # Fallback for stale output files: spatial intersection
            print(f"  [WARN] PARKNAME '{park_key}' not found in all_parks_gdf; "
                  "using spatial intersection fallback. Re-run 01_data_cleaning.py "
                  "to fix.")
            current_park_gdf = all_parks_gdf[
                all_parks_gdf.intersects(study_box.geometry.iloc[0])
            ].copy()
    print(f"  Current park polygons: {len(current_park_gdf)}")

    # ── Per-run figure folder  output/figures/{prefix}/ ───────────────────────
    FIG_OUT = os.path.join(FIG_BASE, prefix)
    os.makedirs(FIG_OUT, exist_ok=True)
    print(f"  Figures → {FIG_OUT}")

    # ── [4] Filter events and build locs_s ─────────────────────────────────────
    park_points, locs_s, baseline_time = filter_events_and_build_locs_s(
        illegal_dumping = illegal_dumping,
        study_box       = study_box,
        cov_gdf         = cov_gdf if FILTER_TO_COV else None,
        filter_to_cov   = FILTER_TO_COV,
        study_years     = STUDY_YEARS,
        offset_seasonal = OFFSET_SEASONAL,
    )

    # ── [5] Validate covariates ────────────────────────────────────────────────
    print("\n[5] Validating covariate names ...")
    active_cov_names = validate_cov_names(COV_NAMES, cov_gdf)

    # ── [6] Set up and fit model ───────────────────────────────────────────────
    model = setup_and_fit_model(
        locs_s           = locs_s,
        study_box        = study_box,
        total_days       = TOTAL_DAYS,
        cov_gdf          = cov_gdf,
        active_cov_names = active_cov_names,
        litter_df        = litter_df,
        cleanup_df       = cleanup_df,
        temporal_window  = TEMPORAL_WINDOW,
        spatial_window   = SPATIAL_WINDOW,
        svi_lr           = SVI_LR,
        svi_num_steps    = SVI_NUM_STEPS,
    )

    # ── [7] Generate all figures ───────────────────────────────────────────────
    print("\n[7] Generating figures ...")

    plot_covariate_check(model, prefix, FIG_OUT, FIG_DPI)

    plot_prop_excitation(model, prefix, FIG_OUT, FIG_DPI)

    plot_trigger_posterior(model, prefix, FIG_OUT, FIG_DPI)

    plot_spatial_trigger_3panel(
        model, study_box, current_park_gdf, layers, prefix, FIG_OUT, FIG_DPI
    )

    plot_uncertainty(model, current_park_gdf, prefix, FIG_OUT, FIG_DPI)

    if PLOT_ROAD_NETWORK:
        plot_road_network(model, study_box, prefix, FIG_OUT, FIG_DPI)
    else:
        print("\n[L] Road-network plot skipped (PLOT_ROAD_NETWORK = False).")

    plot_trigger_time_decay(model, prefix, FIG_OUT, FIG_DPI)

    plot_freq_time_diff(park_points, prefix, FIG_OUT, FIG_DPI)

    plot_trigger_time_distribution(
        model      = model,
        locs_s     = locs_s,
        park_points = park_points,
        total_days = TOTAL_DAYS,
        prefix     = prefix,
        fig_out    = FIG_OUT,
        dpi        = FIG_DPI,
        n_sim      = N_SIM_TRIGGER,
        xlim_max   = TRIGGER_XLIM_MAX,
    )

    plot_temporal_components(
        model, prefix, FIG_OUT,
        start_year = TEMPORAL_START_YEAR,
        dpi        = FIG_DPI,
    )

    plot_seasonal_components(
        model, prefix, FIG_OUT,
        ref_year = SEASONAL_REF_YEAR,
        dpi      = FIG_DPI,
    )

    # ── Summary ────────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"  Mode:             {'city' if USE_CITY else 'district' if USE_DISTRICT else 'park' if USE_PARK else 'custom coordinate'}")
    print(f"  Study area:       {active_park_name}")
    print(f"  Figure prefix:    {prefix}")
    print(f"  Baseline date:    {baseline_time.date()}")
    print(f"  Events in locs_s: {len(locs_s):,}")
    print(f"  Covariates used:  {len(active_cov_names)}")
    print(f"  Litter survey:    {'yes' if litter_df  is not None else 'not found'}")
    print(f"  Cleanup data:     {'yes' if cleanup_df is not None else 'not found'}")
    print(f"  Road network:     {'plotted' if PLOT_ROAD_NETWORK else 'skipped'}")
    print(f"\nFigures written to: {FIG_OUT}/")  
    print("Done.")
