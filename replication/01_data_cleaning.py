"""
01_data_cleaning.py
===================
Prepare all event data and park geometries for the Cox-Hawkes
spatiotemporal illegal-dumping model.

Run order
---------
  THIS SCRIPT  →  02_covariates.py  →  03_analysis.py

Required input files  (place in the paths shown; see README for download links)
-------------------------------------------------------------------------------
  data/311/{year}/public_cases_fc.shp   OpenDataPhilly 311 service records
      year = 2021, 2022, 2023, 2024, 2025
      https://opendataphilly.org/datasets/311-service-and-information-requests/

  data/parks/PPR_Properties.geojson     PPR park boundary polygons
      https://opendataphilly.org/datasets/ppr-properties/

  data/parks/park_AOI_from_google.kml   Custom AOI polygons (Google My Maps,
      Feb 3 2026; provided by research team)

  data/City_Limits/City_Limits.shp      Philadelphia city boundary
      https://opendataphilly.org/datasets/city-limits/
      (only needed for make_centroid_box() boundary check)

  output/cov_cbg.geojson                CBG covariate panel produced by
      02_covariates.py.  Used here only to restrict events to CBGs with full
      covariate coverage.  Run 02_covariates.py first, OR set
      FILTER_TO_COV = False to skip this step.

Optional model-fitting inputs (only needed if running auxiliary model components):
  data/CSVFiles/litter_surveys_raw_20250730.csv
  data/PPR_BaseFiles_IllegalDumpingModel.gdb  (layer: PPR_Ops_Dumping_Cleanup_Complete)

Outputs  (written to replication/output/)
-----------------------------------------
  illegal_dumping_full.geojson   full cleaned illegal-dumping event dataset
  all_parks_gdf.geojson          park boundary polygons
  all_boxes_gdf.geojson          square study-area boxes around each park
  locs_s.csv                     model-ready event table: X, Y, T, A
  litter_df_for_model.csv        litter-survey rows aligned to model time axis
                                   (written only when litter data is available)
  cleanup_df_for_model.csv       PPR cleanup rows aligned to model time axis
                                   (written only when GDB data is available)
"""

# =============================================================================
# SECTION 0 – Imports and global paths
# =============================================================================

import os
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
import shapely.ops
from shapely.geometry import Point, box
from shapely.ops import unary_union

warnings.filterwarnings("ignore")

REPL_DIR = os.path.dirname(os.path.abspath(__file__))   # replication/
DATA     = os.path.join(REPL_DIR, "data")               # replication/data/
OUT      = os.path.join(REPL_DIR, "output")             # replication/output/
os.makedirs(OUT, exist_ok=True)

print("=" * 60)
print("01_data_cleaning.py — event data + park geometry prep")
print("=" * 60)


# =============================================================================
# SECTION 1 – Configuration  (edit these values to change the analysis)
# =============================================================================

# ── Study window ──────────────────────────────────────────────────────────────
STUDY_YEARS = [2021, 2022, 2023, 2024, 2025]
STUDY_END   = "2025-12-31"
TOTAL_DAYS  = 365 * 5          # 1825 days; used as the model time horizon

# ── Park selection ────────────────────────────────────────────────────────────
# ObjectIDs from PPR_Properties.geojson
# To browse available IDs: ppr = gpd.read_file(PPR_PROPERTIES_PATH); ppr[['OBJECTID','PARENT_NAME']].drop_duplicates()
PPR_PARK_IDS = {
    "cobbscreek": [432],
    "mifflin":    [208],
    "tacony":     [10, 64],
    "fairmount":  [7, 9],
}

# 'Name' values inside park_AOI_from_google.kml
# Dict key  = park key used throughout the code (must be snake_case / lowercase)
# Dict value = exact 'Name' attribute in the KML file
AOI_PARK_NAMES = {
    "fairmount_aoi":  "fairmount_AOI",
    "cobbs_aoi":      "cobbs_AOI",
    "tacony_aoi":     "taconyAOI",       # KML name kept as-is; key normalised
    "cobbs_aoi_focus": "cobbs_AOI_focus",
}

# ── Active park (the park whose events are used for locs_s) ──────────────────
# Set to the park key (no "_box" suffix); the suffix is appended automatically.
# Valid park keys:   "mifflin", "cobbscreek", "tacony", "fairmount"
# Valid AOI keys:    "fairmount_aoi", "cobbs_aoi", "tacony_aoi", "cobbs_aoi_focus"
PARK_NAME = "mifflin"

# ── Time jittering (optional) ─────────────────────────────────────────────────
# Public 311 records carry a datetime stamp so ties are rare.
# Setting JITTER_SECONDS > 0 still adds a small Uniform(0, N) random offset
# to ensure no two events share an identical timestamp.
# Set to 0 to disable (default).
JITTER_SECONDS = 0          # e.g.  3600  for a 1-hour jitter
JITTER_SEED    = 999

# ── Seasonal offset ───────────────────────────────────────────────────────────
# Shifts the start of the seasonal cycle; 0 = no shift (Jan 1 = day 0).
OFFSET_SEASONAL = 0

# ── Covariate filtering ───────────────────────────────────────────────────────
# When True, events are additionally restricted to lie inside a CBG that has
# full covariate coverage (requires output/cov_cbg.geojson from 02_covariates.py).
# Set to False if you have not yet run 02_covariates.py.
FILTER_TO_COV = True

# ── File paths ────────────────────────────────────────────────────────────────
PPR_PROPERTIES_PATH = os.path.join(DATA, "parks", "PPR_Properties.geojson")
PARK_AOI_PATH       = os.path.join(DATA, "parks", "park_AOI_from_google.kml")
CITY_LIMITS_PATH    = os.path.join(DATA, "City_Limits", "City_Limits.shp")
COV_CBG_PATH        = os.path.join(OUT,  "cov_cbg.geojson")

LITTER_CSV_PATH  = os.path.join(DATA, "CSVFiles",
                                 "litter_surveys_raw_20250730.csv")
CLEANUP_GDB_PATH = os.path.join(DATA, "PPR_BaseFiles_IllegalDumpingModel.gdb")


# =============================================================================
# SECTION 2 – Functions
# =============================================================================

# ─── 2-A: Public 311 data loading ────────────────────────────────────────────

def load_public_311_year(year: int, data_dir: str) -> gpd.GeoDataFrame:
    """
    Load one year of OpenDataPhilly 311 service records, filter to
    "Illegal Dumping" cases, standardise the timestamp column, and return
    a GeoDataFrame in EPSG:26918.

    Parameters
    ----------
    year : int
        Calendar year (e.g. 2021).
    data_dir : str
        Path to the data/ folder containing  data/311/{year}/public_cases_fc.shp.

    Returns
    -------
    gpd.GeoDataFrame  (EPSG:26918) with columns including
        start_time  (datetime64[ns])
        objectid    unique case identifier
        geometry    point
    """
    shp_path = os.path.join(data_dir, "311", str(year), "public_cases_fc.shp")
    if not os.path.exists(shp_path):
        print(f"  [SKIP] {shp_path} not found.")
        return gpd.GeoDataFrame()

    gdf = gpd.read_file(shp_path)
    gdf = gdf[gdf["service_na"] == "Illegal Dumping"].copy()

    # Standardise timestamp column name
    if "requested_" in gdf.columns:
        gdf = gdf.rename(columns={"requested_": "start_time"})
    gdf["start_time"] = pd.to_datetime(gdf["start_time"], errors="coerce")
    gdf["year"] = year

    return gdf.to_crs(epsg=26918)


def load_and_merge_public_311(
    study_years: list,
    data_dir: str,
    study_end: str,
) -> gpd.GeoDataFrame:
    """
    Load all study-year public 311 shapefiles, filter to "Illegal Dumping"
    service type, concatenate, and apply the date ceiling.

    Parameters
    ----------
    study_years : list[int]
        Years to load (e.g. [2021, 2022, 2023, 2024, 2025]).
    data_dir : str
        Path to the replication data/ folder.
    study_end : str
        ISO-date string; events after this date are excluded.

    Returns
    -------
    gpd.GeoDataFrame  (EPSG:26918)  — deduplicated illegal-dumping events
    """
    print("\n[A] Loading public 311 shapefiles ...")
    parts = []
    for yr in study_years:
        gdf = load_public_311_year(yr, data_dir)
        if not gdf.empty:
            parts.append(gdf)
            print(f"  {yr}: {len(gdf):,} illegal-dumping records")

    if not parts:
        raise FileNotFoundError(
            "No public 311 shapefiles found. "
            "Place data/311/{year}/public_cases_fc.shp in replication/data/."
        )

    illegal_dumping = gpd.GeoDataFrame(
        pd.concat(parts, ignore_index=True), geometry="geometry", crs=26918
    )
    illegal_dumping = illegal_dumping[
        illegal_dumping["start_time"] <= study_end
    ].dropna(subset=["start_time", "geometry"]).reset_index(drop=True)

    print(f"  Total illegal-dumping events "
          f"({min(study_years)}–{max(study_years)}): {len(illegal_dumping):,}")
    return illegal_dumping


# ─── 2-B: Park boundary loading ──────────────────────────────────────────────

def load_park_geometries(
    ppr_path: str,
    aoi_path: str,
    ppr_park_ids: dict,
    aoi_park_names: dict,
) -> dict:
    """
    Load park boundary polygons from the PPR Properties layer and a Google
    My Maps KML file.

    Parameters
    ----------
    ppr_path : str
        Path to PPR_Properties.geojson.
    aoi_path : str
        Path to park_AOI_from_google.kml.
    ppr_park_ids : dict
        {park_key: [OBJECTID, ...]} — which PPR objects to include per park.
    aoi_park_names : dict
        {park_key: 'Name value in KML'} — AOI polygons from the KML file.

    Returns
    -------
    dict  {park_key: GeoDataFrame (EPSG:26918)}
    """
    print("\n[B] Loading park geometries ...")
    parks = {}

    if os.path.exists(ppr_path):
        ppr = gpd.read_file(ppr_path).to_crs(26918)
        for key, ids in ppr_park_ids.items():
            subset = ppr[ppr["OBJECTID"].isin(ids)].copy()
            if subset.empty:
                print(f"  [WARN] No PPR features for '{key}' (OBJECTID in {ids})")
            else:
                parks[key] = subset
                print(f"  {key}: {len(subset)} PPR polygon(s)")
    else:
        print(f"  [WARN] PPR Properties not found: {ppr_path}")

    if os.path.exists(aoi_path):
        park_aoi = gpd.read_file(aoi_path, driver="KML").to_crs(26918)
        for key, name_val in aoi_park_names.items():
            subset = park_aoi[park_aoi["Name"] == name_val].copy()
            if subset.empty:
                print(f"  [WARN] No AOI features named '{name_val}' in KML")
            else:
                parks[key] = subset
                print(f"  {key}: {len(subset)} AOI polygon(s)")
    else:
        print(f"  [WARN] Park AOI KML not found: {aoi_path}")

    return parks


# ─── 2-C: Park-box creation ───────────────────────────────────────────────────

# Buffer distances applied to the actual park polygon boundary.
# AOI keys contain "aoi" (case-insensitive); all other keys are PPR parks.
_BUFFER_PARK = 500   # metres — PPR park keys (mifflin, cobbscreek, tacony, fairmount)
_BUFFER_AOI  = 500   # metres — AOI keys  (fairmount_aoi, cobbs_aoi, tacony_aoi, cobbs_aoi_focus …)


def create_parks_and_boxes(parks: dict) -> tuple:
    """
    Build two GeoDataFrames:
      all_parks_gdf  – the raw park polygon(s) with a PARKNAME column
      all_boxes_gdf  – a buffered study-area polygon per park

    Study-area construction
    -----------------------
    The bounding box of each park's polygon(s) is expanded by a fixed margin
    on all four sides.  This always produces a sharp-cornered rectangle
    (no rounded edges).

      • PPR park keys (mifflin, cobbscreek, tacony, fairmount …) → +300 m
      • AOI keys (names containing "aoi", case-insensitive)       → +500 m

    AOI keys produce their own separate polygon so they can be selected
    independently in downstream analysis (e.g. "fairmount" vs "fairmount_aoi").

    Special case: 'cobbs_aoi_focus' keeps its original polygon with no buffer.

    Parameters
    ----------
    parks : dict
        {park_key: GeoDataFrame}  — output of load_park_geometries().

    Returns
    -------
    (all_parks_gdf, all_boxes_gdf) : (GeoDataFrame, GeoDataFrame)
        Both in EPSG:26918 with columns ['PARKNAME', 'geometry'].
    """
    all_parks_gdf = gpd.GeoDataFrame(columns=["PARKNAME", "geometry"])
    all_boxes_gdf = gpd.GeoDataFrame(columns=["PARKNAME", "geometry"])

    for name, gdf in parks.items():
        if gdf is None or gdf.empty:
            continue

        tmp = gdf.copy()
        if tmp.crs is None:
            tmp = tmp.set_crs(4326).to_crs(26918)
        elif tmp.crs.to_epsg() != 26918:
            tmp = tmp.to_crs(26918)

        # Normalise PARKNAME to the park key so downstream filtering is reliable.
        # (KML AOI parks carry a raw "Name" like "cobbs_AOI"; we overwrite that
        # with the consistent snake_case key, e.g. "cobbs_aoi".)
        tmp = tmp.copy()
        tmp["PARKNAME"] = name
        tmp = tmp[["PARKNAME", "geometry"]]
        all_parks_gdf = pd.concat([all_parks_gdf, tmp], ignore_index=True)

        # Special case: cobbs_aoi_focus keeps its polygon as-is (no buffer)
        if name == "cobbs_aoi_focus":
            focus = tmp.copy()
            focus["geometry"] = focus["geometry"].apply(
                lambda g: shapely.force_2d(g) if g is not None else g
            )
            focus["PARKNAME"] = name
            all_boxes_gdf = pd.concat(
                [all_boxes_gdf, focus[["PARKNAME", "geometry"]]],
                ignore_index=True,
            )
            continue

        # Union all parts of this park into one polygon, then find the four
        # extreme points (N/S/E/W edges) and extend each outward by buf_dist.
        # This always produces a sharp-cornered rectangle — no rounded edges.
        park_union = shapely.ops.unary_union(tmp.geometry)
        buf_dist   = _BUFFER_AOI if "aoi" in name.lower() else _BUFFER_PARK
        minx, miny, maxx, maxy = park_union.bounds  # W, S, E, N extreme points
        buffered = box(
            minx - buf_dist, miny - buf_dist,
            maxx + buf_dist, maxy + buf_dist,
        )

        area_km2 = buffered.area / 1_000_000
        w_km = (maxx - minx + 2 * buf_dist) / 1000
        h_km = (maxy - miny + 2 * buf_dist) / 1000
        print(
            f"  {name}_box: park bounds "
            f"W={minx:.0f} S={miny:.0f} E={maxx:.0f} N={maxy:.0f} "
            f"+ {buf_dist} m → {w_km:.2f}×{h_km:.2f} km, "
            f"area≈{area_km2:.2f} km²"
        )

        buf_gdf = gpd.GeoDataFrame(
            {"PARKNAME": [f"{name}_box"]},
            geometry=[buffered],
            crs=tmp.crs,
        )
        all_boxes_gdf = pd.concat(
            [all_boxes_gdf, buf_gdf[["PARKNAME", "geometry"]]],
            ignore_index=True,
        )

    all_parks_gdf = gpd.GeoDataFrame(all_parks_gdf, geometry="geometry", crs=26918)
    all_boxes_gdf = gpd.GeoDataFrame(all_boxes_gdf, geometry="geometry", crs=26918)
    return all_parks_gdf, all_boxes_gdf


# ─── 2-D: Custom centroid box (new) ───────────────────────────────────────────

def make_centroid_box(
    lat: float,
    lon: float,
    city_limits_path: str,
    side_meters: float = 1000.0,
    park_name: str = "custom",
) -> gpd.GeoDataFrame | None:
    """
    Create a square study-area box centred on a user-supplied coordinate.

    The input must be a WGS84 (lat, lon) pair.  The function projects the
    point to EPSG:26918 (NAD83 / UTM 18N, metres), draws a square of the
    requested side length, and checks that the centre falls inside the
    Philadelphia city boundary.  A warning is printed and None is returned
    if the coordinate lies outside Philadelphia.

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees (WGS84).
    lon : float
        Longitude in decimal degrees (WGS84).
    city_limits_path : str
        Path to the Philadelphia City Limits shapefile.
    side_meters : float
        Side length of the square in metres (default 1 000 m = 1 km²).
    park_name : str
        Label stored in the 'PARKNAME' column of the output GeoDataFrame.

    Returns
    -------
    gpd.GeoDataFrame  (EPSG:26918, PARKNAME = park_name)  or  None
    """
    # Project input point to EPSG:26918
    point_26918 = gpd.GeoDataFrame(
        geometry=[Point(lon, lat)], crs=4326
    ).to_crs(26918)
    px = point_26918.geometry.iloc[0].x
    py = point_26918.geometry.iloc[0].y

    # Check containment within Philadelphia boundary
    if not os.path.exists(city_limits_path):
        warnings.warn(
            f"City limits file not found at {city_limits_path}. "
            "Boundary check skipped.",
            stacklevel=2,
        )
    else:
        philly = gpd.read_file(city_limits_path).to_crs(4326)
        philly_union = unary_union(philly.geometry)
        if not philly_union.contains(Point(lon, lat)):
            warnings.warn(
                f"Coordinate ({lat}, {lon}) is OUTSIDE the Philadelphia "
                "city boundary. No box was created.",
                stacklevel=2,
            )
            return None

    half = side_meters / 2.0
    square_geom = box(px - half, py - half, px + half, py + half)
    return gpd.GeoDataFrame(
        {"PARKNAME": [park_name]},
        geometry=[square_geom],
        crs=26918,
    )


# ─── 2-E: Spatial filter events to park box ───────────────────────────────────

def filter_events_to_park(
    illegal_dumping: gpd.GeoDataFrame,
    all_boxes_gdf: gpd.GeoDataFrame,
    park_name: str,
    cov_gdf: gpd.GeoDataFrame | None = None,
) -> gpd.GeoDataFrame:
    """
    Retain only illegal-dumping events that fall within the selected park box.
    Optionally further restrict to events whose CBG has full covariate coverage.

    Parameters
    ----------
    illegal_dumping : GeoDataFrame  (EPSG:26918)
        Full cleaned illegal-dumping event dataset.
    all_boxes_gdf : GeoDataFrame  (EPSG:26918)
        Park bounding boxes from create_parks_and_boxes().
    park_name : str
        PARKNAME value of the target box (e.g. "mifflin_box").
    cov_gdf : GeoDataFrame | None
        CBG covariate panel (from 02_covariates.py).  When provided, events
        are additionally filtered to lie inside a CBG present in cov_gdf.

    Returns
    -------
    GeoDataFrame  — events within the park box (EPSG:26918)
    """
    print(f"\n[E] Filtering events to park box: {park_name} ...")

    park_points = illegal_dumping.sjoin(all_boxes_gdf, predicate="within")
    park_points = park_points.drop(columns=["index_right"])
    park_points = park_points[park_points["PARKNAME"] == park_name].copy()
    print(f"  Events inside box: {len(park_points):,}")

    if cov_gdf is not None:
        park_points = park_points.sjoin(
            cov_gdf[["geometry"]], predicate="within", how="inner"
        ).drop(columns=["index_right"])
        print(f"  Events after covariate-coverage filter: {len(park_points):,}")

    return park_points.reset_index(drop=True)


# ─── 2-F: Build model-ready locs_s ───────────────────────────────────────────

def build_locs_s(
    park_points: gpd.GeoDataFrame,
    offset_seasonal: int = 0,
    jitter_seconds: int = 0,
    seed: int = 999,
) -> tuple:
    """
    Construct the model-ready event table ``locs_s`` and return the baseline
    time used for the temporal axis.

    Columns produced
    ----------------
    X : float   easting  (EPSG:26918, metres)
    Y : float   northing (EPSG:26918, metres)
    T : float   days since Jan 1 of the earliest event year
    A : float   seasonal position  =  (T + offset_seasonal) % 365

    Parameters
    ----------
    park_points : GeoDataFrame  (EPSG:26918)
        Filtered events from filter_events_to_park().
    offset_seasonal : int
        Shifts the origin of the seasonal cycle (days).
    jitter_seconds : int
        If > 0, adds Uniform(0, jitter_seconds) random noise to each event's
        timestamp before computing T.
    seed : int
        Random seed for reproducibility when jitter_seconds > 0.

    Returns
    -------
    (locs_s, baseline_time) : (pd.DataFrame, pd.Timestamp)
    """
    print("\n[F] Building locs_s ...")

    pts = park_points.copy()
    pts["start_time"] = pd.to_datetime(pts["start_time"], errors="coerce")
    pts = pts.dropna(subset=["start_time"])

    # Optional temporal jitter
    if jitter_seconds > 0:
        rng = np.random.default_rng(seed)
        offsets = rng.uniform(0, jitter_seconds, size=len(pts))
        pts["start_time"] = pts["start_time"] + pd.to_timedelta(offsets, unit="s")
        print(f"  Jitter applied: Uniform(0, {jitter_seconds} s) per event")

    min_year  = pts["start_time"].min()
    baseline  = pd.Timestamp(f"{min_year.year}-01-01")
    time_days = (pts["start_time"] - baseline).dt.total_seconds() / (24 * 60 * 60)
    seasonal  = (time_days + offset_seasonal) % 365

    locs_s = pd.DataFrame({
        "X": pts.geometry.centroid.x.astype(float).values,
        "Y": pts.geometry.centroid.y.astype(float).values,
        "T": time_days.astype(float).values,
        "A": seasonal.astype(float).values,
    })

    print(f"  locs_s shape: {locs_s.shape}")
    print(f"  Baseline time: {baseline.date()}")
    print(f"  T range: [{locs_s['T'].min():.1f}, {locs_s['T'].max():.1f}] days")
    return locs_s, baseline


# ─── 2-G: Optional — litter survey data prep ─────────────────────────────────

def prep_litter_survey(
    litter_path: str,
    baseline_time: pd.Timestamp,
    all_boxes_gdf: gpd.GeoDataFrame,
    park_name: str,
    locs_s: pd.DataFrame,
    offset_seasonal: int = 0,
) -> pd.DataFrame | None:
    """
    Load litter survey CSV, project to EPSG:26918, align time axis to the
    model's baseline, and restrict to the park box.

    Returns a DataFrame with columns  [X, Y, T, A, score]  ready to pass to
    ``bstpp.cox_hawkes_shared.attach_litter_survey_args()``, or None if the
    file is missing.

    Parameters
    ----------
    litter_path : str
        Path to litter_surveys_raw_20250730.csv.
    baseline_time : pd.Timestamp
        Jan 1 of the earliest event year (from build_locs_s).
    all_boxes_gdf : GeoDataFrame  (EPSG:26918)
    park_name : str
        PARKNAME value of the park box.
    locs_s : pd.DataFrame
        Model event table; used to determine the T ceiling.
    offset_seasonal : int
        Must match the value used in build_locs_s.
    """
    if not os.path.exists(litter_path):
        print("  [SKIP] Litter survey file not found.")
        return None

    print("\n[G] Preparing litter survey data ...")
    df = pd.read_csv(litter_path)
    df["survey_start"] = pd.to_datetime(df["survey_start"], errors="coerce")
    df = df.dropna(subset=["survey_start", "x", "y"])

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["x"], df["y"]),
        crs=4326,
    ).to_crs(26918)

    gdf["X"] = gdf.geometry.x.astype(float)
    gdf["Y"] = gdf.geometry.y.astype(float)
    gdf["score"] = gdf["Litter Rating"].astype(float).astype(int) - 1
    gdf["T"] = (
        (gdf["survey_start"] - baseline_time).dt.total_seconds() / (24 * 60 * 60)
    )
    gdf["A"] = (gdf["T"] + offset_seasonal) % 365

    park_boundary = all_boxes_gdf[all_boxes_gdf["PARKNAME"] == park_name].to_crs(26918)
    litter_in_park = gpd.sjoin(gdf, park_boundary, how="inner", predicate="within")

    t_max  = locs_s["T"].max()
    result = litter_in_park[["X", "Y", "T", "A", "score"]].copy()
    result = result[(result["T"] > 0) & (result["T"] < t_max)].reset_index(drop=True)

    print(f"  Litter rows within park & time window: {len(result):,}")
    return result


# ─── 2-H: Optional — PPR cleanup data prep ───────────────────────────────────

def prep_cleanup_data(
    gdb_path: str,
    baseline_time: pd.Timestamp,
    all_boxes_gdf: gpd.GeoDataFrame,
    park_name: str,
    locs_s: pd.DataFrame,
    offset_seasonal: int = 0,
) -> pd.DataFrame | None:
    """
    Load the PPR cleanup operations layer from the GDB, project, align the
    time axis, and restrict to the park box.

    Returns a DataFrame with columns  [X, Y, T, A, log_cleanup_cost]  ready
    for ``bstpp.cox_hawkes_shared.attach_dumping_report_args()``, or None if
    the GDB is missing.

    Parameters
    ----------
    gdb_path : str
        Path to PPR_BaseFiles_IllegalDumpingModel.gdb.
    baseline_time : pd.Timestamp
        Jan 1 of the earliest event year (from build_locs_s).
    all_boxes_gdf : GeoDataFrame  (EPSG:26918)
    park_name : str
        PARKNAME value of the park box.
    locs_s : pd.DataFrame
        Model event table; used to determine the T ceiling.
    offset_seasonal : int
        Must match the value used in build_locs_s.
    """
    if not os.path.exists(gdb_path):
        print("  [SKIP] PPR GDB not found.")
        return None

    print("\n[H] Preparing PPR cleanup data ...")
    cleanup_gdf = gpd.read_file(
        gdb_path, layer="PPR_Ops_Dumping_Cleanup_Complete"
    ).to_crs(26918)

    # GDB timestamps carry timezone; align to UTC baseline
    baseline_utc = pd.Timestamp(
        f"{baseline_time.year}-01-01 17:00:00", tz="UTC"
    )
    cleanup_gdf["T"] = (
        (cleanup_gdf["CLEANUP_DATE"] - baseline_utc).dt.total_seconds()
        / (24 * 60 * 60)
    )
    cleanup_gdf["A"] = (cleanup_gdf["T"] + offset_seasonal) % 365

    park_boundary = all_boxes_gdf[all_boxes_gdf["PARKNAME"] == park_name].to_crs(26918)
    within_park = gpd.sjoin(cleanup_gdf, park_boundary, how="inner", predicate="within")
    within_park = within_park[
        within_park.geometry.within(park_boundary.iloc[0]["geometry"])
    ].copy()

    within_park["X"] = within_park.geometry.x.astype(float)
    within_park["Y"] = within_park.geometry.y.astype(float)

    t_max  = locs_s["T"].max()
    result = within_park[["X", "Y", "T", "A"]].copy()
    result["log_cleanup_cost"] = np.log(
        within_park["TOTAL_CLEANUP_COST"].replace(0, np.nan)
    )
    result = result[
        (result["T"] > 0) & (result["T"] < t_max)
    ].dropna(subset=["log_cleanup_cost"]).reset_index(drop=True)

    print(f"  Cleanup rows within park & time window: {len(result):,}")
    return result


# =============================================================================
# SECTION 3 – Main execution  (runs when script is called directly)
# =============================================================================

if __name__ == "__main__":

    # ── [1] Load public 311 data ──────────────────────────────────────────────
    print("\n[1] Loading public 311 shapefiles ...")
    illegal_dumping_full = load_and_merge_public_311(
        study_years = STUDY_YEARS,
        data_dir    = DATA,
        study_end   = STUDY_END,
    )

    out_dump = os.path.join(OUT, "illegal_dumping_full.geojson")
    illegal_dumping_full.to_crs(4326).to_file(out_dump, driver="GeoJSON")
    print(f"  Saved → {out_dump}")

    # ── [2] Load park geometries ──────────────────────────────────────────────
    print("\n[2] Loading park geometries ...")
    parks = load_park_geometries(
        ppr_path       = PPR_PROPERTIES_PATH,
        aoi_path       = PARK_AOI_PATH,
        ppr_park_ids   = PPR_PARK_IDS,
        aoi_park_names = AOI_PARK_NAMES,
    )

    # ── [3] Build park bounding boxes ─────────────────────────────────────────
    print("\n[3] Building park bounding boxes ...")
    all_parks_gdf, all_boxes_gdf = create_parks_and_boxes(parks)
    print(f"  Parks: {len(all_parks_gdf)}, Boxes: {len(all_boxes_gdf)}")
    print(f"  Box names: {all_boxes_gdf['PARKNAME'].tolist()}")

    out_parks = os.path.join(OUT, "all_parks_gdf.geojson")
    out_boxes = os.path.join(OUT, "all_boxes_gdf.geojson")
    all_parks_gdf.to_crs(4326).to_file(out_parks, driver="GeoJSON")
    all_boxes_gdf.to_crs(4326).to_file(out_boxes, driver="GeoJSON")
    print(f"  Saved → {out_parks}")
    print(f"  Saved → {out_boxes}")

    # ── [Optional] Custom centroid box ───────────────────────────────────────
    # To study a custom area centred on a lat/lon coordinate, call:
    #
    #   custom_box = make_centroid_box(
    #       lat=39.9526, lon=-75.1652,
    #       city_limits_path=CITY_LIMITS_PATH,
    #       side_meters=1000,
    #       park_name="custom_area",
    #   )
    #   if custom_box is not None:
    #       all_boxes_gdf = pd.concat([all_boxes_gdf, custom_box],
    #                                  ignore_index=True)
    #       all_boxes_gdf = gpd.GeoDataFrame(all_boxes_gdf,
    #                                         geometry="geometry", crs=26918)
    #
    # The function raises a warning and returns None if the point is
    # outside the Philadelphia city boundary.

    # ── [4] Load covariates (optional) ───────────────────────────────────────
    cov_gdf = None
    if FILTER_TO_COV:
        if os.path.exists(COV_CBG_PATH):
            cov_gdf = gpd.read_file(COV_CBG_PATH)
            print(f"\n[4] Covariate panel loaded: {len(cov_gdf)} CBGs")
        else:
            print(
                "\n[4] cov_cbg.geojson not found — skipping covariate filter.\n"
                "    Run 02_covariates.py first to enable this filter."
            )

    # ── [5] Filter events to selected park box ────────────────────────────────
    BOX_NAME = f"{PARK_NAME}_box"
    park_points = filter_events_to_park(
        illegal_dumping = illegal_dumping_full,
        all_boxes_gdf   = all_boxes_gdf,
        park_name       = BOX_NAME,
        cov_gdf         = cov_gdf,
    )

    # ── [6] Build locs_s ──────────────────────────────────────────────────────
    locs_s, baseline_time = build_locs_s(
        park_points     = park_points,
        offset_seasonal = OFFSET_SEASONAL,
        jitter_seconds  = JITTER_SECONDS,
        seed            = JITTER_SEED,
    )

    out_locs = os.path.join(OUT, "locs_s.csv")
    locs_s.to_csv(out_locs, index=False)
    print(f"  Saved → {out_locs}")

    # ── [7] Optional auxiliary model data ────────────────────────────────────
    litter_df = prep_litter_survey(
        litter_path     = LITTER_CSV_PATH,
        baseline_time   = baseline_time,
        all_boxes_gdf   = all_boxes_gdf,
        park_name       = BOX_NAME,
        locs_s          = locs_s,
        offset_seasonal = OFFSET_SEASONAL,
    )
    if litter_df is not None:
        out_litter = os.path.join(OUT, "litter_df_for_model.csv")
        litter_df.to_csv(out_litter, index=False)
        print(f"  Saved → {out_litter}")

    cleanup_df = prep_cleanup_data(
        gdb_path        = CLEANUP_GDB_PATH,
        baseline_time   = baseline_time,
        all_boxes_gdf   = all_boxes_gdf,
        park_name       = BOX_NAME,
        locs_s          = locs_s,
        offset_seasonal = OFFSET_SEASONAL,
    )
    if cleanup_df is not None:
        out_cleanup = os.path.join(OUT, "cleanup_df_for_model.csv")
        cleanup_df.to_csv(out_cleanup, index=False)
        print(f"  Saved → {out_cleanup}")

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"  Park:               {PARK_NAME} ({BOX_NAME})")
    print(f"  Events in locs_s:   {len(locs_s):,}")
    print(f"  Baseline date:      {baseline_time.date()}")
    print(f"  Study window:       {TOTAL_DAYS} days")
    print(f"  Jitter:             {JITTER_SECONDS} seconds")
    print(f"  Seasonal offset:    {OFFSET_SEASONAL} days")
    print(f"  Litter survey:      {'yes' if litter_df is not None else 'not found'}")
    print(f"  Cleanup data:       {'yes' if cleanup_df is not None else 'not found'}")
    print("\nOutputs written to:", OUT)
    print("Done.")
