"""
03a_analysis_by_division.py
============================
Spatial-hotspot views of illegal-dumping reports across two official
Philadelphia geographies, instead of the single study area used by
03_analysis.py.

Two methods, chosen by DIVISION below:

  DIVISION = "planning_district"   (18 units)
      Re-fits the full Cox-Hawkes model independently per planning district
      (reusing filter_events_and_build_locs_s / validate_cov_names /
      setup_and_fit_model from 03_analysis.py) and maps a handful of scalar
      posterior summaries (self-excitation share, Hawkes alpha, background
      level, covariate effects) across the 18 districts.

  DIVISION = "neighborhood"   (~150 units)
      Too many units to safely re-fit a full spatiotemporal model per unit,
      so this is raw event density (count / area) with no model fitting --
      a standard descriptive hotspot map.

This script only *reads* functions and config constants already defined in
03_analysis.py (loaded as a module below). It does not modify that file.

Required inputs
----------------
  replication/output/illegal_dumping_full.geojson   from 01_data_cleaning.py
  replication/output/cov_cbg.geojson                from 02_covariates.py
  replication/data/Planning_Districts.geojson       (planning_district mode)
  data/philadelphia-neighborhoods.geojson           (neighborhood mode; repo root)

Outputs
-------
  planning_district mode:
    output/planning_district_fit_summary.csv / .geojson
    output/figures/planning_district_maps/pd_map_{stat}.png

  neighborhood mode:
    output/neighborhood_density.csv / .geojson
    output/figures/neighborhood_density_map.png
"""

# =============================================================================
# SECTION 0 - Imports
# =============================================================================

import os
import time
import traceback
import warnings
import importlib.util

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

REPL_DIR = os.path.dirname(os.path.abspath(__file__))   # replication/
ROOT_DIR = os.path.dirname(REPL_DIR)                     # repo root
OUT      = os.path.join(REPL_DIR, "output")
FIG_BASE = os.path.join(OUT, "figures")

# ── Load 03_analysis.py as a module ───────────────────────────────────────────
# Its filename starts with a digit, so it can't be `import`-ed by name; load it
# by path instead. This reuses its functions/config without copying or editing
# that file.
_spec = importlib.util.spec_from_file_location(
    "analysis_base", os.path.join(REPL_DIR, "03_analysis.py")
)
analysis_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(analysis_base)

filter_events_and_build_locs_s = analysis_base.filter_events_and_build_locs_s
validate_cov_names              = analysis_base.validate_cov_names
setup_and_fit_model              = analysis_base.setup_and_fit_model

COV_NAMES       = analysis_base.COV_NAMES
OFFSET_SEASONAL = analysis_base.OFFSET_SEASONAL
SVI_LR          = analysis_base.SVI_LR
SVI_NUM_STEPS   = analysis_base.SVI_NUM_STEPS
FILTER_TO_COV   = analysis_base.FILTER_TO_COV

ILLEGAL_DUMPING_PATH = analysis_base.ILLEGAL_DUMPING_PATH
COV_CBG_PATH         = analysis_base.COV_CBG_PATH

print("=" * 60)
print("03a_analysis_by_division.py — hotspot maps by spatial division")
print("=" * 60)

# =============================================================================
# SECTION 1 - Configuration  (edit these values to change the analysis)
# =============================================================================

DIVISION = "planning_district"   # "planning_district" | "neighborhood"

# ── Study years ────────────────────────────────────────────────────────────────
# Years to include in this batch run, independent of 03_analysis.py's own
# STUDY_YEARS config. TOTAL_DAYS is derived automatically from this range.
STUDY_YEARS = [2021, 2022, 2023, 2024, 2025]              # e.g. [2021, 2022, 2023, 2024, 2025]
TOTAL_DAYS  = (
    pd.Timestamp(f"{max(STUDY_YEARS) + 1}-01-01")
    - pd.Timestamp(f"{min(STUDY_YEARS)}-01-01")
).days
print(f"  Study years: {STUDY_YEARS} → TOTAL_DAYS = {TOTAL_DAYS}")

# ── Model window parameters ──────────────────────────────────────────────────
# Independent of 03_analysis.py's own TEMPORAL_WINDOW/SPATIAL_WINDOW config --
# see 03_analysis.py Section 1 for what these control.
TEMPORAL_WINDOW = 0.8       # temporal grid cells for background GP
SPATIAL_WINDOW  = 0.025     # spatial bandwidth for background GP

PLANNING_DISTRICTS_PATH  = os.path.join(REPL_DIR, "data", "Planning_Districts.geojson")
PLANNING_DISTRICT_ID_COL = "dist_name"

NEIGHBORHOODS_PATH  = os.path.join(ROOT_DIR, "data", "philadelphia-neighborhoods.geojson")
NEIGHBORHOOD_ID_COL = "NAME"

# Optional: restrict to a handful of unit ids for a quick smoke test before
# running the full batch (18 districts / ~150 neighborhoods).
UNIT_SUBSET = None   # e.g. ["CENTRAL", "SOUTH"] or None for all units

FIG_DPI = 300

DIV_FIG_OUT = os.path.join(FIG_BASE, f"{DIVISION}_maps")
os.makedirs(DIV_FIG_OUT, exist_ok=True)


# =============================================================================
# SECTION 2 - Functions
# =============================================================================

def build_unit_study_box(geometry, unit_id: str, crs) -> gpd.GeoDataFrame:
    """
    Wrap a single boundary-layer unit's own polygon (no bounding box, no
    buffer) as the single-row study-area GeoDataFrame expected by
    filter_events_and_build_locs_s / setup_and_fit_model. The customized
    bstpp model already handles irregular polygons correctly (area-fraction
    correction + grid-cell intersection), so the unit's real shape is used
    directly, unlike the buffered rectangular boxes used for parks.
    """
    return gpd.GeoDataFrame({"PARKNAME": [unit_id]}, geometry=[geometry], crs=crs)


def extract_fit_stats(model, active_cov_names: list) -> dict:
    """
    Pull a small set of scalar posterior summaries out of a fitted
    Cox_Hawkes_Shared model. Walks model.samples generically rather than
    hardcoding trigger-kernel parameter names, since those vary by the
    configured temporal-trigger kernel.
    """
    samples = model.samples
    stats = {}

    if "Itot_excite" in samples and "Itot_txy" in samples:
        ratio = np.asarray(samples["Itot_excite"]) / np.asarray(samples["Itot_txy"])
        stats["prop_excitation_mean"] = float(np.mean(ratio))

    if "alpha" in samples:
        alpha = np.asarray(samples["alpha"])
        stats["alpha_mean"]  = float(np.mean(alpha))
        stats["alpha_p_pos"] = float(np.mean(alpha > 0))

    if "a_0" in samples:
        stats["a_0_mean"] = float(np.mean(np.asarray(samples["a_0"])))

    if "w" in samples:
        w_mean = np.asarray(samples["w"]).mean(axis=0).reshape(-1)
        for name, val in zip(active_cov_names, w_mean):
            stats[f"w_mean_{name}"] = float(val)

    # Any other scalar-per-posterior-draw parameter not already captured
    # above (e.g. trigger time-decay rate) -- whatever the configured
    # Temporal_* kernel exposes, picked up generically by shape.
    n_samples = len(np.asarray(samples.get("alpha", samples.get("a_0", []))))
    handled = {"Itot_excite", "Itot_txy", "alpha", "a_0", "w"}
    for key, val in samples.items():
        if key in handled:
            continue
        arr = np.asarray(val)
        if arr.ndim == 1 and n_samples and arr.shape[0] == n_samples:
            stats[f"{key}_mean"] = float(np.mean(arr))

    return stats


def plot_choropleth(
    gdf: gpd.GeoDataFrame,
    column: str,
    title: str,
    out_path: str,
    dpi: int = 300,
    cmap: str = "OrRd",
    boundary_overlay: gpd.GeoDataFrame | None = None,
) -> None:
    """
    Save a single choropleth PNG of `column` over `gdf`'s geometry.
    `boundary_overlay`, if given, is drawn as outline-only on top (e.g. the
    18 planning-district boundaries for spatial reference).
    """
    if column not in gdf.columns or gdf[column].dropna().empty:
        print(f"  [skip] '{column}' has no data to map.")
        return
    fig, ax = plt.subplots(figsize=(9, 9))
    gdf.plot(
        column=column, ax=ax, legend=True, cmap=cmap,
        edgecolor="white", linewidth=0.2,
        missing_kwds={"color": "lightgray", "label": "no data"},
    )
    if boundary_overlay is not None:
        boundary_overlay.boundary.plot(ax=ax, color="black", linewidth=0.8)
    ax.set_title(title)
    ax.set_axis_off()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")


def extract_cbg_fxy(model, district_geometry, cov_gdf: gpd.GeoDataFrame, crs) -> gpd.GeoDataFrame:
    """
    Aggregate the fitted spatial GP surface (model.samples["f_xy"]) from the
    model's own grid (model.comp_grid -- the same real-world-coordinate grid
    cells used by plot_uncertainty/plot_road_network in 03_analysis.py) onto
    CBG polygons, the same granularity used for the covariates.

    Values are then min-max normalized *within this district* (0-1), since
    each district is fit independently and isn't on a shared absolute scale
    -- see the f_xy mosaic discussion: this answers "where within this
    district is the model surface highest/lowest", not absolute cross-
    district comparison (the scalar-stat maps already cover that).

    Returns a GeoDataFrame with columns GEOID, f_xy_norm, geometry.
    """
    f_xy_mean = np.asarray(model.samples["f_xy"]).mean(axis=0)
    grid = model.comp_grid[["geometry"]].copy()
    grid["f_xy_mean"] = f_xy_mean
    if grid.crs is None:
        grid = grid.set_crs(crs)

    district_cbgs = cov_gdf[cov_gdf.intersects(district_geometry)][["GEOID", "geometry"]].copy()
    if district_cbgs.empty:
        return gpd.GeoDataFrame(columns=["GEOID", "f_xy_norm", "geometry"], geometry="geometry", crs=crs)

    joined = gpd.sjoin(grid, district_cbgs, predicate="intersects")
    cbg_means = joined.groupby("GEOID")["f_xy_mean"].mean()

    vmin, vmax = cbg_means.min(), cbg_means.max()
    normalized = (cbg_means - vmin) / (vmax - vmin) if vmax > vmin else cbg_means * 0 + 0.5

    out = district_cbgs.set_index("GEOID").loc[normalized.index].copy()
    out["f_xy_norm"] = normalized
    return out.reset_index()[["GEOID", "f_xy_norm", "geometry"]]


def extract_grid_fxy(model, district_geometry, crs) -> gpd.GeoDataFrame:
    """
    Same within-district min-max normalization as extract_cbg_fxy, but kept
    at the model's own native grid resolution (model.comp_grid, the 25x25
    cells underlying model.samples["f_xy"]) instead of aggregating onto CBG
    polygons -- finer detail, at the cost of cells not aligning to any
    administrative boundary. Grid cells are clipped to the district's actual
    polygon so the rectangular bounding-box grid doesn't spill into
    neighboring districts.
    """
    f_xy_mean = np.asarray(model.samples["f_xy"]).mean(axis=0)
    grid = model.comp_grid[["geometry"]].copy()
    grid["f_xy_mean"] = f_xy_mean
    if grid.crs is None:
        grid = grid.set_crs(crs)

    vmin, vmax = grid["f_xy_mean"].min(), grid["f_xy_mean"].max()
    grid["f_xy_norm"] = (
        (grid["f_xy_mean"] - vmin) / (vmax - vmin) if vmax > vmin else grid["f_xy_mean"] * 0 + 0.5
    )

    clipped = grid.clip(district_geometry)
    return clipped[["f_xy_norm", "geometry"]]


# =============================================================================
# SECTION 3 - Planning-district mode (model fit per district)
# =============================================================================

def run_planning_district_mode() -> None:
    print("\n[1] Loading inputs (planning_district mode) ...")
    illegal_dumping = gpd.read_file(ILLEGAL_DUMPING_PATH).to_crs(26918)
    cov_gdf         = gpd.read_file(COV_CBG_PATH).to_crs(26918)
    districts       = gpd.read_file(PLANNING_DISTRICTS_PATH).to_crs(26918)

    if UNIT_SUBSET is not None:
        districts = districts[districts[PLANNING_DISTRICT_ID_COL].isin(UNIT_SUBSET)].copy()
    print(f"  Districts to fit: {len(districts)}")

    rows = []
    cbg_fxy_rows = []
    grid_fxy_rows = []
    for _, district in districts.iterrows():
        unit_id = district[PLANNING_DISTRICT_ID_COL]
        print(f"\n--- District: {unit_id} ---")
        t0 = time.time()
        row = {PLANNING_DISTRICT_ID_COL: unit_id, "geometry": district["geometry"]}
        try:
            study_box = build_unit_study_box(district["geometry"], unit_id, districts.crs)

            _, locs_s, _ = filter_events_and_build_locs_s(
                illegal_dumping = illegal_dumping,
                study_box       = study_box,
                cov_gdf         = cov_gdf if FILTER_TO_COV else None,
                filter_to_cov   = FILTER_TO_COV,
                study_years     = STUDY_YEARS,
                offset_seasonal = OFFSET_SEASONAL,
            )
            row["n_events"] = len(locs_s)

            active_cov_names = validate_cov_names(COV_NAMES, cov_gdf)

            model = setup_and_fit_model(
                locs_s           = locs_s,
                study_box        = study_box,
                total_days       = TOTAL_DAYS,
                cov_gdf          = cov_gdf,
                active_cov_names = active_cov_names,
                litter_df        = None,
                cleanup_df       = None,
                temporal_window  = TEMPORAL_WINDOW,
                spatial_window   = SPATIAL_WINDOW,
                svi_lr           = SVI_LR,
                svi_num_steps    = SVI_NUM_STEPS,
            )
            row.update(extract_fit_stats(model, active_cov_names))
            cbg_fxy_rows.append(extract_cbg_fxy(model, district["geometry"], cov_gdf, districts.crs))
            grid_fxy_rows.append(extract_grid_fxy(model, district["geometry"], districts.crs))
            row["status"] = "fit"
            row["error"]  = None
        except Exception as exc:
            row["status"] = "failed"
            row["error"]  = str(exc)
            print(f"  [FAILED] {unit_id}: {exc}")
            traceback.print_exc()
        row["runtime_sec"] = round(time.time() - t0, 1)
        rows.append(row)

    summary = gpd.GeoDataFrame(pd.DataFrame(rows), geometry="geometry", crs=districts.crs)

    out_geojson = os.path.join(OUT, "planning_district_fit_summary.geojson")
    out_csv     = os.path.join(OUT, "planning_district_fit_summary.csv")
    summary.to_crs(4326).to_file(out_geojson, driver="GeoJSON")
    summary.drop(columns="geometry").to_csv(out_csv, index=False)
    print(f"\nSaved summary → {out_geojson}, {out_csv}")

    print("\n[2] Rendering choropleth maps ...")
    for col in ["prop_excitation_mean", "alpha_mean", "a_0_mean"]:
        plot_choropleth(
            summary, col,
            title=f"{col} by planning district",
            out_path=os.path.join(DIV_FIG_OUT, f"pd_map_{col}.png"),
            dpi=FIG_DPI,
        )

    print("\n[3] Building CBG-level spatial-surface mosaic ...")
    if cbg_fxy_rows:
        cbg_mosaic = pd.concat(cbg_fxy_rows, ignore_index=True)
        cbg_mosaic = gpd.GeoDataFrame(cbg_mosaic, geometry="geometry", crs=districts.crs)
        cbg_mosaic = cbg_mosaic.drop_duplicates(subset="GEOID", keep="first")

        out_geojson = os.path.join(OUT, "planning_district_cbg_fxy.geojson")
        cbg_mosaic.to_crs(4326).to_file(out_geojson, driver="GeoJSON")
        print(f"  Saved → {out_geojson}")

        plot_choropleth(
            cbg_mosaic, "f_xy_norm",
            title="Spatial intensity surface (f_xy, normalized within district) by CBG",
            out_path=os.path.join(DIV_FIG_OUT, "pd_map_f_xy_norm_cbg.png"),
            dpi=FIG_DPI,
            cmap="viridis",
            boundary_overlay=districts,
        )
    else:
        print("  [skip] No districts fit successfully; nothing to mosaic.")

    print("\n[4] Building native-grid (25×25 per district) spatial-surface mosaic ...")
    if grid_fxy_rows:
        grid_mosaic = pd.concat(grid_fxy_rows, ignore_index=True)
        grid_mosaic = gpd.GeoDataFrame(grid_mosaic, geometry="geometry", crs=districts.crs)

        out_geojson = os.path.join(OUT, "planning_district_grid_fxy.geojson")
        grid_mosaic.to_crs(4326).to_file(out_geojson, driver="GeoJSON")
        print(f"  Saved → {out_geojson}")

        plot_choropleth(
            grid_mosaic, "f_xy_norm",
            title="Spatial intensity surface (f_xy, normalized within district) — native 25×25 grid",
            out_path=os.path.join(DIV_FIG_OUT, "pd_map_f_xy_norm_grid.png"),
            dpi=FIG_DPI,
            cmap="viridis",
            boundary_overlay=districts,
        )
    else:
        print("  [skip] No districts fit successfully; nothing to mosaic.")

    n_fit = int((summary["status"] == "fit").sum())
    print(f"\nDone. {n_fit}/{len(summary)} districts fit successfully.")


# =============================================================================
# SECTION 4 - Neighborhood mode (density only, no model fit)
# =============================================================================

def run_neighborhood_mode() -> None:
    print("\n[1] Loading inputs (neighborhood mode) ...")
    illegal_dumping = gpd.read_file(ILLEGAL_DUMPING_PATH).to_crs(26918)
    neighborhoods   = gpd.read_file(NEIGHBORHOODS_PATH).to_crs(26918)

    if UNIT_SUBSET is not None:
        neighborhoods = neighborhoods[neighborhoods[NEIGHBORHOOD_ID_COL].isin(UNIT_SUBSET)].copy()
    print(f"  Neighborhoods: {len(neighborhoods)}")

    events = illegal_dumping.copy()
    if STUDY_YEARS is not None:
        events["start_time"] = pd.to_datetime(events["start_time"], errors="coerce")
        events = events.dropna(subset=["start_time"])
        events = events[events["start_time"].dt.year.isin(STUDY_YEARS)]

    joined = events.sjoin(
        neighborhoods[[NEIGHBORHOOD_ID_COL, "geometry"]], predicate="within"
    )
    counts = joined.groupby(NEIGHBORHOOD_ID_COL).size()

    neighborhoods = neighborhoods.set_index(NEIGHBORHOOD_ID_COL)
    neighborhoods["n_events"]     = counts.reindex(neighborhoods.index).fillna(0).astype(int)
    neighborhoods["area_km2"]     = neighborhoods.geometry.area / 1e6
    neighborhoods["rate_per_km2"] = neighborhoods["n_events"] / neighborhoods["area_km2"]
    neighborhoods = neighborhoods.reset_index()

    out_geojson = os.path.join(OUT, "neighborhood_density.geojson")
    out_csv     = os.path.join(OUT, "neighborhood_density.csv")
    neighborhoods.to_crs(4326).to_file(out_geojson, driver="GeoJSON")
    neighborhoods.drop(columns="geometry").to_csv(out_csv, index=False)
    print(f"  Total events mapped: {int(neighborhoods['n_events'].sum()):,} "
          f"(of {len(events):,} filtered citywide events)")
    print(f"  Saved → {out_geojson}, {out_csv}")

    print("\n[2] Rendering density choropleth ...")
    plot_choropleth(
        neighborhoods, "rate_per_km2",
        title="Illegal-dumping density by neighborhood (events / km²)",
        out_path=os.path.join(FIG_BASE, "neighborhood_density_map.png"),
        dpi=FIG_DPI,
    )
    print("\nDone.")


# =============================================================================
# SECTION 5 - Main execution
# =============================================================================

if __name__ == "__main__":
    if DIVISION == "planning_district":
        run_planning_district_mode()
    elif DIVISION == "neighborhood":
        run_neighborhood_mode()
    else:
        raise ValueError(
            f"Unknown DIVISION: {DIVISION!r}. Use 'planning_district' or 'neighborhood'."
        )
