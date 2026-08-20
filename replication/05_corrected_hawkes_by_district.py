"""
05_corrected_hawkes_by_district.py
===================================
Reporting-corrected Cox-Hawkes fits per planning district (all 18 by
default). Injects the log reporting-thinning offset

    offset_i = beta * w_R_i          [log p(s), per CBG]

into the background intensity of the per-district Cox-Hawkes model
(03a_analysis_by_division.py's planning-district mode). The offset comes
from 05a_streetlight_decomposition.py's single CITYWIDE decomposition
(NT = Street Light Outage, poles as exposure): one citywide beta and a
per-CBG w_R field that is directly comparable across the whole city, so no
district-level anchor term is needed — a constant level shift is absorbed
by each district's a_0. See 05_corrected_hawkes.md for derivation.

Three specs per district:
  baseline      — no offset (re-fit: unlike the old 03a run, `reporting_rate`
                  is dropped from the covariates so all three specs share
                  the same covariate set)
  offset_fixed  — offset coefficient pinned at 1 (main spec: true thinning)
  offset_free   — offset coefficient ~ Normal(1, 0.5) (robustness: gamma ~ 1
                  validates the correction from the point-process side)

For every district, the full 03-style diagnostic figure set (temporal +
seasonal components, self-excitation / trigger figures, covariate forest,
spatial surface, uncertainty, road network, Δt distributions) is rendered
for the baseline and offset_fixed fits under
output/figures/{RUN_NAME}/districts/{district}/{spec}/, plus per-district
Δ f_xy_norm maps (offset_fixed − baseline, CBG + native grid) at each
district folder's top level.

The offset-aware model lives in cox_hawkes_offset.py; this script re-uses
03/03a functions by loading them as modules.

Required inputs
----------------
  replication/output/illegal_dumping_full.geojson         01_data_cleaning.py
  replication/output/cov_cbg.geojson                      02_covariates.py
  replication/output/corrected_hawkes_streetlight/
      reporting_decomp_cbg_sl.geojson                     05a
      reporting_decomp_city_summary_sl.csv                05a
  replication/data/Planning_Districts.geojson

Outputs (per run, RUN_NAME = "corrected_hawkes" by default)
-------
  output/{RUN_NAME}/fit_summary.csv / .geojson            district x spec
  output/{RUN_NAME}/cbg_fxy.geojson                       CBG mosaic, all specs
  output/{RUN_NAME}/grid_fxy.geojson                      native-grid mosaic
  output/{RUN_NAME}/offset_cbg.geojson                    the offset itself
  output/figures/{RUN_NAME}/*.png                         citywide maps
  output/figures/{RUN_NAME}/districts/*/                  per-district figures
"""

# =============================================================================
# SECTION 0 - Imports
# =============================================================================

import os
import sys
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
import numpyro.distributions as dist

warnings.filterwarnings("ignore")

REPL_DIR = os.path.dirname(os.path.abspath(__file__))   # replication/
OUT      = os.path.join(REPL_DIR, "output")

# Each run writes its data files into output/{RUN_NAME}/ and its figures into
# output/figures/{RUN_NAME}/.
RUN_NAME = "corrected_hawkes"
RUN_OUT  = os.path.join(OUT, RUN_NAME)
FIG_OUT  = os.path.join(OUT, "figures", RUN_NAME)

# ── Load 03a as a module (which itself loads 03_analysis.py) ─────────────────
_spec = importlib.util.spec_from_file_location(
    "analysis_by_division", os.path.join(REPL_DIR, "03a_analysis_by_division.py")
)
div = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(div)
base = div.analysis_base            # the loaded 03_analysis.py module

# reused pieces — nothing in 03/03a is modified
filter_events_and_build_locs_s = base.filter_events_and_build_locs_s
validate_cov_names             = base.validate_cov_names
build_unit_study_box = div.build_unit_study_box
extract_fit_stats    = div.extract_fit_stats
extract_cbg_fxy      = div.extract_cbg_fxy
extract_grid_fxy     = div.extract_grid_fxy
plot_choropleth      = div.plot_choropleth

sys.path.insert(0, REPL_DIR)
from cox_hawkes_offset import Cox_Hawkes_Shared_Offset

print("=" * 60)
print("05_corrected_hawkes_by_district.py — reporting-corrected Cox-Hawkes")
print("=" * 60)

# =============================================================================
# SECTION 1 - Configuration
# =============================================================================

# Same study window / model settings as the 03a planning-district run
STUDY_YEARS = [2021, 2022, 2023, 2024, 2025]
TOTAL_DAYS  = (
    pd.Timestamp(f"{max(STUDY_YEARS) + 1}-01-01")
    - pd.Timestamp(f"{min(STUDY_YEARS)}-01-01")
).days
TEMPORAL_WINDOW = 0.8
SPATIAL_WINDOW  = 0.025

PLANNING_DISTRICTS_PATH  = os.path.join(REPL_DIR, "data", "Planning_Districts.geojson")
PLANNING_DISTRICT_ID_COL = "dist_name"

# 05a citywide streetlight-decomposition outputs feeding the offset.
DECOMP_DIR        = os.path.join(OUT, "corrected_hawkes_streetlight")
DECOMP_CBG_PATH   = os.path.join(DECOMP_DIR, "reporting_decomp_cbg_sl.geojson")
CITY_SUMMARY_PATH = os.path.join(DECOMP_DIR, "reporting_decomp_city_summary_sl.csv")
DECOMP_SPEC       = "sl"                  # "sl_pre24" = PSIP-robust variant
DECOMP_WR_COL     = {"sl": "w_R_mean_sl", "sl_pre24": "w_R_mean_pre24"}[DECOMP_SPEC]

OFFSET_COL       = "log_p_report"         # column added to the covariate gdf
FREE_COEF_PRIOR  = dist.Normal(1.0, 0.5)  # offset_free spec

# `reporting_rate` is dropped in ALL specs: the offset replaces it, and the
# baseline must share the covariate set for a clean comparison.
COV_NAMES = [c for c in base.COV_NAMES if c != "reporting_rate"]

SPECS = ["baseline", "offset_fixed", "offset_free"]

# Quick smoke test before the full batch, e.g. ["Central"]. None = all 18.
UNIT_SUBSET = None

# Fit one district with an all-zero offset and compare with baseline: the
# posteriors must match (same PRNG seed), proving both injection points are
# consistent. Runs before the main batch on the first district.
RUN_EQUIVALENCE_CHECK = False

FIG_DPI = 300

# ── Per-district 03-style diagnostic figures ─────────────────────────────────
MAKE_DISTRICT_DIAGNOSTICS = True
DIAG_SPECS        = ["baseline", "offset_fixed"]   # specs that get figure sets
PLOT_ROAD_NETWORK = True        # OSMnx download per district (disk-cached)
N_SIM_TRIGGER     = 500         # MC draws for the Δt-distribution figure
                                # (dominant runtime cost; 100 = fast mode)
MAKE_GRID_FXY_MAP = True        # citywide native-grid f_xy mosaic maps


def slug(name: str) -> str:
    """Filename-safe district id: 'University Southwest' → 'university_southwest'."""
    return name.lower().replace(" ", "_")


# Disk caches so the road-network / basemap figures download each district's
# data once, not once per spec. Kept out of output/ — these are transient
# download caches, not analysis outputs.
CACHE_DIR = os.path.join(REPL_DIR, ".cache")
try:
    import osmnx as ox
    ox.settings.use_cache = True
    ox.settings.cache_folder = os.path.join(CACHE_DIR, "osmnx")
except ImportError:
    pass
if hasattr(base, "ctx"):
    try:
        base.ctx.set_cache_dir(os.path.join(CACHE_DIR, "contextily"))
    except Exception:
        pass


# =============================================================================
# SECTION 2 - Offset construction (citywide streetlight decomposition, 05a)
# =============================================================================

def build_offset(cov_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Compute offset_i = beta * w_R_i from 05a's citywide streetlight
    decomposition and merge it into cov_gdf as OFFSET_COL.

    The decomposition is a single citywide fit, so w_R is directly comparable
    across all CBGs and one beta applies everywhere. There is no district
    level term: a citywide-constant shift is absorbed identically by each
    district's a_0, leaving the fitted a_0's comparable across districts.

    CBGs without a decomposition estimate get 0 (the citywide mean).
    """
    print("\n[offset] Building citywide reporting offset  beta * w_R ...")
    decomp = gpd.read_file(DECOMP_CBG_PATH)[
        ["GEOID", "dist_name", DECOMP_WR_COL, "T_count"]
    ]
    decomp["GEOID"] = decomp["GEOID"].astype(str)
    city = pd.read_csv(CITY_SUMMARY_PATH)
    beta = float(city.loc[city["spec"] == DECOMP_SPEC, "beta_mean"].iloc[0])
    print(f"  beta (citywide, spec '{DECOMP_SPEC}') = {beta:.4f}")

    decomp[OFFSET_COL] = beta * decomp[DECOMP_WR_COL]
    cov_gdf = cov_gdf.merge(
        decomp[["GEOID", "dist_name", OFFSET_COL, "T_count"]],
        on="GEOID", how="left",
    )

    n_missing = int(cov_gdf[OFFSET_COL].isna().sum())
    cov_gdf[OFFSET_COL] = cov_gdf[OFFSET_COL].fillna(0.0)
    cov_gdf["T_count"] = cov_gdf["T_count"].fillna(0.0)
    print(f"  Offset merged: {len(cov_gdf) - n_missing}/{len(cov_gdf)} CBGs "
          f"from 05a; {n_missing} filled with 0 (the citywide mean)")
    print(f"  Offset range: [{cov_gdf[OFFSET_COL].min():.3f}, "
          f"{cov_gdf[OFFSET_COL].max():.3f}], "
          f"mean {cov_gdf[OFFSET_COL].mean():.3f}")

    # Human-readable unit: how many of a CBG's observed reports are
    # attributable to reporting bias rather than true dumping, on an
    # absolute-count scale (see 05_corrected_hawkes.md for the offset's
    # log-rate derivation).
    cov_gdf["excess_reports"] = (
        cov_gdf["T_count"] - cov_gdf["T_count"] / np.exp(cov_gdf[OFFSET_COL])
    )
    print(f"  Excess reports (bias-attributable) range: "
          f"[{cov_gdf['excess_reports'].min():.1f}, "
          f"{cov_gdf['excess_reports'].max():.1f}], "
          f"total {cov_gdf['excess_reports'].sum():.1f} of "
          f"{cov_gdf['T_count'].sum():,.0f} citywide reports")

    return cov_gdf


# =============================================================================
# SECTION 3 - Model setup (clone of 03's setup_and_fit_model, offset-aware)
# =============================================================================

def setup_and_fit_model_offset(
    locs_s, study_box, total_days, cov_gdf, active_cov_names,
    temporal_window, spatial_window, svi_lr, svi_num_steps,
    spatial_offset=None, offset_coef_prior=None,
):
    """
    Same construction as 03_analysis.py's setup_and_fit_model (same priors,
    trigger kernel, windows, SVI settings; litter/cleanup attachments not
    used here), but builds Cox_Hawkes_Shared_Offset so a per-CBG background
    offset can be supplied. spatial_offset=None reproduces the base model.
    """
    model = Cox_Hawkes_Shared_Offset(
        locs_s,
        study_box,
        total_days,
        True,
        spatial_cov=cov_gdf,
        cov_names=active_cov_names,
        standardize_cov=True,
        a_0=dist.Normal(0, 0.5),
        alpha=dist.Beta(2, 4),
        temporal_trig=base.Temporal_Power_Law,
        spatial_offset=spatial_offset,
        offset_coef_prior=offset_coef_prior,
    )
    model.args["model"] = "cox_hawkes"
    model.set_window(window=temporal_window, spatial_window=spatial_window)
    model.run_svi(lr=svi_lr, num_steps=svi_num_steps)
    return model


SPEC_KWARGS = {
    "baseline":     dict(spatial_offset=None,       offset_coef_prior=None),
    "offset_fixed": dict(spatial_offset=OFFSET_COL, offset_coef_prior=None),
    "offset_free":  dict(spatial_offset=OFFSET_COL, offset_coef_prior=FREE_COEF_PRIOR),
}


def offset_coef_stats(model) -> dict:
    """Posterior mean and 94% CI of the free offset coefficient."""
    if "offset_coef" not in model.samples:
        return {}
    g = np.asarray(model.samples["offset_coef"]).reshape(-1)
    return {
        "offset_coef_mean": float(np.mean(g)),
        "offset_coef_lo":   float(np.quantile(g, 0.03)),
        "offset_coef_hi":   float(np.quantile(g, 0.97)),
    }


def a_0_stats(model) -> dict:
    a0 = np.asarray(model.samples["a_0"]).reshape(-1)
    return {
        "a_0_lo": float(np.quantile(a0, 0.03)),
        "a_0_hi": float(np.quantile(a0, 0.97)),
    }


# =============================================================================
# SECTION 4 - Equivalence check (offset ≡ 0 must reproduce baseline)
# =============================================================================

def run_equivalence_check(illegal_dumping, cov_gdf, districts, active_cov_names):
    """
    Fit the first district twice — baseline vs offset_fixed with the offset
    column zeroed — and compare key posterior means. bstpp seeds SVI with a
    fixed PRNGKey, so a missed injection point would show up as a mismatch.
    """
    district = districts.iloc[0]
    unit_id  = district[PLANNING_DISTRICT_ID_COL]
    print(f"\n[equivalence check] district = {unit_id}")

    study_box = build_unit_study_box(district["geometry"], unit_id, districts.crs)
    _, locs_s, _ = filter_events_and_build_locs_s(
        illegal_dumping=illegal_dumping, study_box=study_box,
        cov_gdf=cov_gdf if base.FILTER_TO_COV else None,
        filter_to_cov=base.FILTER_TO_COV,
        study_years=STUDY_YEARS, offset_seasonal=base.OFFSET_SEASONAL,
    )

    cov_zero = cov_gdf.copy()
    cov_zero[OFFSET_COL] = 0.0

    results = {}
    for name, cg, kw in [
        ("baseline",    cov_gdf,  SPEC_KWARGS["baseline"]),
        ("zero_offset", cov_zero, SPEC_KWARGS["offset_fixed"]),
    ]:
        model = setup_and_fit_model_offset(
            locs_s, study_box, TOTAL_DAYS, cg, active_cov_names,
            TEMPORAL_WINDOW, SPATIAL_WINDOW, base.SVI_LR, base.SVI_NUM_STEPS,
            **kw,
        )
        results[name] = extract_fit_stats(model, active_cov_names)

    print(f"  {'stat':<28}{'baseline':>12}{'zero_offset':>12}")
    worst = 0.0
    for k in sorted(set(results["baseline"]) & set(results["zero_offset"])):
        b, z = results["baseline"][k], results["zero_offset"][k]
        diff = abs(b - z)
        worst = max(worst, diff)
        print(f"  {k:<28}{b:>12.5f}{z:>12.5f}   |diff|={diff:.2e}")
    print(f"  worst |diff| = {worst:.2e} "
          f"({'OK — specs equivalent' if worst < 1e-4 else 'MISMATCH — check injection points!'})")
    return worst


# =============================================================================
# SECTION 4b - Per-district 03-style diagnostic figures
# =============================================================================

def make_district_diagnostics(model, spec, study_box, park_points, locs_s,
                              layers, dist_dir, dist_slug):
    """
    Render the full 03_analysis.py diagnostic figure set for one fitted
    (district, spec) model under {dist_dir}/{spec}/, prefixed
    {dist_slug}_{spec}. Each figure is wrapped in its own try/except so a
    plotting failure (e.g. a network outage for the basemap / OSMnx figures)
    never marks the fit itself as failed.

    study_box doubles as the boundary-overlay GeoDataFrame, exactly as 03's
    own district mode does (03_analysis.py builds current_park_gdf from
    study_box when USE_DISTRICT is on).
    """
    spec_dir = os.path.join(dist_dir, spec)
    os.makedirs(spec_dir, exist_ok=True)
    prefix = f"{dist_slug}_{spec}"
    print(f"  [{spec}] diagnostics → {spec_dir}")

    def _try(fn, *args, **kwargs):
        try:
            fn(*args, **kwargs)
        except Exception:
            traceback.print_exc()

    _try(base.plot_covariate_check, model, prefix, spec_dir, FIG_DPI)
    _try(base.plot_prop_excitation, model, prefix, spec_dir, FIG_DPI)
    _try(base.plot_trigger_posterior, model, prefix, spec_dir, FIG_DPI)
    _try(base.plot_spatial_trigger_3panel,
         model, study_box, study_box, layers, prefix, spec_dir, FIG_DPI)
    _try(base.plot_uncertainty, model, study_box, prefix, spec_dir, FIG_DPI)
    if PLOT_ROAD_NETWORK:
        _try(base.plot_road_network, model, study_box, prefix, spec_dir, FIG_DPI)
    _try(base.plot_trigger_time_decay, model, prefix, spec_dir, FIG_DPI)
    _try(base.plot_trigger_time_distribution,
         model=model, locs_s=locs_s, park_points=park_points,
         total_days=TOTAL_DAYS, prefix=prefix, fig_out=spec_dir,
         dpi=FIG_DPI, n_sim=N_SIM_TRIGGER, xlim_max=base.TRIGGER_XLIM_MAX)
    _try(base.plot_temporal_components, model, prefix, spec_dir,
         start_year=base.TEMPORAL_START_YEAR, dpi=FIG_DPI)
    _try(base.plot_seasonal_components, model, prefix, spec_dir,
         ref_year=base.SEASONAL_REF_YEAR, dpi=FIG_DPI)


# =============================================================================
# SECTION 5 - Main batch: districts × specs
# =============================================================================

def main():
    os.makedirs(RUN_OUT, exist_ok=True)
    os.makedirs(FIG_OUT, exist_ok=True)
    print(f"\n[1] Loading inputs ... (outputs → {RUN_OUT})")
    illegal_dumping = gpd.read_file(base.ILLEGAL_DUMPING_PATH).to_crs(26918)
    cov_gdf         = gpd.read_file(base.COV_CBG_PATH).to_crs(26918)
    districts       = gpd.read_file(PLANNING_DISTRICTS_PATH).to_crs(26918)

    if UNIT_SUBSET is not None:
        districts = districts[districts[PLANNING_DISTRICT_ID_COL].isin(UNIT_SUBSET)].copy()

    cov_gdf = build_offset(cov_gdf)
    active_cov_names = validate_cov_names(COV_NAMES, cov_gdf)
    print(f"  Districts: {len(districts)}, specs: {SPECS}, "
          f"covariates: {len(active_cov_names)} (reporting_rate dropped)")
    if MAKE_DISTRICT_DIAGNOSTICS and len(districts) > 2:
        print(f"  [note] per-district diagnostics ON for {len(districts)} districts "
              f"× {DIAG_SPECS} — expect a multi-hour run "
              f"(N_SIM_TRIGGER={N_SIM_TRIGGER}, road network={PLOT_ROAD_NETWORK}).")

    layers = None
    if MAKE_DISTRICT_DIAGNOSTICS:
        layers = base.load_supporting_layers(
            base.TIGER_BG_PATH, base.LAND_CARE_PATH, base.LAND_USE_PATH,
            base.VACANT_BLOCK_PATH, base.CITY_LIMITS_PATH,
        )

    # save the offset itself (map + record of what was fed to the model)
    off_path = os.path.join(RUN_OUT, "offset_cbg.geojson")
    cov_gdf[["GEOID", "dist_name", OFFSET_COL, "geometry"]].to_crs(4326).to_file(
        off_path, driver="GeoJSON")
    print(f"  Saved offset → {off_path}")
    plot_choropleth(
        cov_gdf, OFFSET_COL,
        title="Log reporting-thinning offset  β·ŵ_R  (citywide streetlight decomposition)",
        out_path=os.path.join(FIG_OUT, "ch_map_offset_cbg.png"),
        dpi=FIG_DPI, cmap="RdBu_r", boundary_overlay=districts,
    )
    plot_delta_map(
        cov_gdf, "excess_reports",
        title="Citywide: Excess 311 Reports Attributable to Reporting Bias\n"
              "(observed T_count − expected true count), CBG",
        out_path=os.path.join(FIG_OUT, "ch_map_excess_reports_cbg.png"),
        boundary=districts, dpi=FIG_DPI,
    )

    if RUN_EQUIVALENCE_CHECK:
        run_equivalence_check(illegal_dumping, cov_gdf, districts, active_cov_names)

    print("\n[2] Fitting districts × specs ...")
    rows = []
    cbg_fxy  = {s: [] for s in SPECS}
    grid_fxy = {s: [] for s in DIAG_SPECS}
    for _, district in districts.iterrows():
        unit_id = district[PLANNING_DISTRICT_ID_COL]
        print(f"\n--- District: {unit_id} ---")

        try:
            study_box = build_unit_study_box(district["geometry"], unit_id, districts.crs)
            park_points, locs_s, _ = filter_events_and_build_locs_s(
                illegal_dumping=illegal_dumping, study_box=study_box,
                cov_gdf=cov_gdf if base.FILTER_TO_COV else None,
                filter_to_cov=base.FILTER_TO_COV,
                study_years=STUDY_YEARS, offset_seasonal=base.OFFSET_SEASONAL,
            )
        except Exception as exc:
            for spec in SPECS:
                rows.append({PLANNING_DISTRICT_ID_COL: unit_id, "spec": spec,
                             "geometry": district["geometry"],
                             "status": "failed", "error": f"event prep: {exc}"})
            traceback.print_exc()
            continue

        dist_slug = slug(unit_id)
        dist_dir  = os.path.join(FIG_OUT, "districts", dist_slug)
        if MAKE_DISTRICT_DIAGNOSTICS:
            os.makedirs(dist_dir, exist_ok=True)
            try:  # data-only figure — identical across specs, made once
                base.plot_freq_time_diff(park_points, dist_slug, dist_dir, FIG_DPI)
            except Exception:
                traceback.print_exc()

        for spec in SPECS:
            print(f"  [{spec}]")
            t0 = time.time()
            row = {PLANNING_DISTRICT_ID_COL: unit_id, "spec": spec,
                   "geometry": district["geometry"], "n_events": len(locs_s)}
            try:
                model = setup_and_fit_model_offset(
                    locs_s, study_box, TOTAL_DAYS, cov_gdf, active_cov_names,
                    TEMPORAL_WINDOW, SPATIAL_WINDOW,
                    base.SVI_LR, base.SVI_NUM_STEPS,
                    **SPEC_KWARGS[spec],
                )
                row.update(extract_fit_stats(model, active_cov_names))
                row.update(a_0_stats(model))
                row.update(offset_coef_stats(model))
                fxy = extract_cbg_fxy(model, district["geometry"], cov_gdf, districts.crs)
                cbg_fxy[spec].append(
                    fxy.assign(spec=spec, **{PLANNING_DISTRICT_ID_COL: unit_id}))
                if MAKE_GRID_FXY_MAP and spec in grid_fxy:
                    gfx = extract_grid_fxy(model, district["geometry"], districts.crs)
                    grid_fxy[spec].append(
                        gfx.assign(spec=spec, **{PLANNING_DISTRICT_ID_COL: unit_id}))
                row["status"], row["error"] = "fit", None
                if MAKE_DISTRICT_DIAGNOSTICS and spec in DIAG_SPECS:
                    make_district_diagnostics(
                        model, spec, study_box, park_points, locs_s,
                        layers, dist_dir, dist_slug,
                    )
            except Exception as exc:
                row["status"], row["error"] = "failed", str(exc)
                traceback.print_exc()
            row["runtime_sec"] = round(time.time() - t0, 1)
            rows.append(row)
            plt.close("all")

    summary = gpd.GeoDataFrame(pd.DataFrame(rows), geometry="geometry",
                               crs=districts.crs)

    out_geojson = os.path.join(RUN_OUT, "fit_summary.geojson")
    out_csv     = os.path.join(RUN_OUT, "fit_summary.csv")
    summary.to_crs(4326).to_file(out_geojson, driver="GeoJSON")
    summary.drop(columns="geometry").to_csv(out_csv, index=False)
    print(f"\nSaved summary → {out_geojson}, {out_csv}")

    make_outputs(summary, cbg_fxy, districts, grid_fxy)

    n_fit = int((summary["status"] == "fit").sum())
    print(f"\nDone. {n_fit}/{len(summary)} fits successful.")


# =============================================================================
# SECTION 6 - Mosaics and comparison figures
# =============================================================================

def plot_delta_map(gdf, column, title, out_path, boundary, dpi):
    """
    Diverging choropleth for a delta column, color scale symmetric about 0
    (plot_choropleth doesn't center its colormap, which is misleading for
    signed deltas on a single district).
    """
    vals = gdf[column].dropna()
    if vals.empty:
        print(f"  [skip] '{column}' has no data to map.")
        return
    a = float(np.abs(vals).max()) or 1e-9
    fig, ax = plt.subplots(figsize=(9, 9))
    gdf.plot(
        column=column, ax=ax, legend=True, cmap="RdBu_r", vmin=-a, vmax=a,
        edgecolor="white", linewidth=0.2,
        missing_kwds={"color": "lightgray", "label": "no data"},
    )
    if boundary is not None:
        boundary.boundary.plot(ax=ax, color="black", linewidth=0.8)
    ax.set_title(title)
    ax.set_axis_off()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")


def make_district_delta_figures(cbg_by_spec, grid_by_spec, districts):
    """
    Per-district Δ f_xy_norm maps (offset_fixed − baseline, both normalized
    0–1 within the district) saved into each district's own figure folder:

        districts/{slug}/{slug}_f_xy_norm_delta_cbg.png     CBG aggregation
        districts/{slug}/{slug}_f_xy_norm_delta_grid.png    native model grid

    cbg_by_spec / grid_by_spec are {spec: GeoDataFrame} with a
    PLANNING_DISTRICT_ID_COL column, PRE-dedup so boundary CBGs shared by
    two districts keep each district's own within-district normalization.
    CBG rows merge on GEOID; grid cells have no id, so they merge on
    geometry (WKB) — baseline and offset_fixed share the district's
    comp_grid, so the clipped cells are identical across specs.
    """
    have_cbg  = all(s in cbg_by_spec for s in ("baseline", "offset_fixed"))
    have_grid = all(s in grid_by_spec for s in ("baseline", "offset_fixed"))
    if not (have_cbg or have_grid):
        return
    print("  Per-district Δ f_xy_norm maps ...")
    for _, district in districts.iterrows():
        unit_id   = district[PLANNING_DISTRICT_ID_COL]
        dist_slug = slug(unit_id)
        dist_dir  = os.path.join(FIG_OUT, "districts", dist_slug)
        os.makedirs(dist_dir, exist_ok=True)
        boundary = gpd.GeoDataFrame(geometry=[district["geometry"]],
                                    crs=districts.crs)

        if have_cbg:
            b = cbg_by_spec["baseline"]
            c = cbg_by_spec["offset_fixed"]
            b = b[b[PLANNING_DISTRICT_ID_COL] == unit_id]
            c = c[c[PLANNING_DISTRICT_ID_COL] == unit_id]
            if len(b) and len(c):
                d = b[["GEOID", "f_xy_norm", "geometry"]].merge(
                    c[["GEOID", "f_xy_norm"]],
                    on="GEOID", suffixes=("_base", "_corr"))
                d["f_xy_norm_delta"] = d["f_xy_norm_corr"] - d["f_xy_norm_base"]
                plot_delta_map(
                    gpd.GeoDataFrame(d, geometry="geometry", crs=districts.crs),
                    "f_xy_norm_delta",
                    title=f"{unit_id}: Δ f_xy_norm (corrected − baseline), CBG\n"
                          f"+ = hotter after reporting correction",
                    out_path=os.path.join(
                        dist_dir, f"{dist_slug}_f_xy_norm_delta_cbg.png"),
                    boundary=boundary, dpi=FIG_DPI,
                )

        if have_grid:
            b = grid_by_spec["baseline"]
            c = grid_by_spec["offset_fixed"]
            b = b[b[PLANNING_DISTRICT_ID_COL] == unit_id].copy()
            c = c[c[PLANNING_DISTRICT_ID_COL] == unit_id].copy()
            if len(b) and len(c):
                b["_cell"] = b.geometry.apply(lambda g: g.wkb)
                c["_cell"] = c.geometry.apply(lambda g: g.wkb)
                d = b[["_cell", "f_xy_norm", "geometry"]].merge(
                    c[["_cell", "f_xy_norm"]],
                    on="_cell", suffixes=("_base", "_corr")).drop(columns="_cell")
                d["f_xy_norm_delta"] = d["f_xy_norm_corr"] - d["f_xy_norm_base"]
                plot_delta_map(
                    gpd.GeoDataFrame(d, geometry="geometry", crs=districts.crs),
                    "f_xy_norm_delta",
                    title=f"{unit_id}: Δ f_xy_norm (corrected − baseline), "
                          f"native grid\n+ = hotter after reporting correction",
                    out_path=os.path.join(
                        dist_dir, f"{dist_slug}_f_xy_norm_delta_grid.png"),
                    boundary=boundary, dpi=FIG_DPI,
                )


def make_outputs(summary, cbg_fxy, districts, grid_fxy=None):
    print("\n[3] Mosaics and comparison figures ...")

    # ── CBG f_xy mosaics, one geojson with a spec column ──────────────────────
    mosaics = {}
    cbg_full = {}   # pre-dedup, keeps every district's own copy of shared CBGs
    frames = []
    for spec in SPECS:
        if not cbg_fxy[spec]:
            continue
        m = pd.concat(cbg_fxy[spec], ignore_index=True)
        m = gpd.GeoDataFrame(m, geometry="geometry", crs=districts.crs)
        cbg_full[spec] = m
        m = m.drop_duplicates(subset="GEOID", keep="first")
        mosaics[spec] = m
        frames.append(m)
    if frames:
        all_m = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True),
                                 geometry="geometry", crs=districts.crs)
        path = os.path.join(RUN_OUT, "cbg_fxy.geojson")
        all_m.to_crs(4326).to_file(path, driver="GeoJSON")
        print(f"  Saved → {path}")

    for spec in ["baseline", "offset_fixed"]:
        if spec in mosaics:
            plot_choropleth(
                mosaics[spec], "f_xy_norm",
                title=f"Spatial surface f_xy ({spec}, normalized within district)",
                out_path=os.path.join(FIG_OUT, f"ch_map_f_xy_norm_{spec}.png"),
                dpi=FIG_DPI, cmap="viridis", boundary_overlay=districts,
            )

    # ── native-grid mosaics (smoother sub-CBG texture, same normalization) ───
    grid_full = {}
    if grid_fxy:
        grid_frames = []
        for spec, parts in grid_fxy.items():
            if not parts:
                continue
            g = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True),
                                 geometry="geometry", crs=districts.crs)
            grid_full[spec] = g
            grid_frames.append(g)
            plot_choropleth(
                g, "f_xy_norm",
                title=f"Spatial surface f_xy ({spec}, native grid, "
                      f"normalized within district)",
                out_path=os.path.join(FIG_OUT, f"ch_map_f_xy_norm_{spec}_grid.png"),
                dpi=FIG_DPI, cmap="viridis", boundary_overlay=districts,
            )
        if grid_frames:
            all_g = gpd.GeoDataFrame(pd.concat(grid_frames, ignore_index=True),
                                     geometry="geometry", crs=districts.crs)
            path = os.path.join(RUN_OUT, "grid_fxy.geojson")
            all_g.to_crs(4326).to_file(path, driver="GeoJSON")
            print(f"  Saved → {path}")

    # ── (a) delta map: corrected − baseline (both within-district 0–1) ────────
    if "baseline" in mosaics and "offset_fixed" in mosaics:
        d = mosaics["baseline"][["GEOID", "f_xy_norm", "geometry"]].merge(
            mosaics["offset_fixed"][["GEOID", "f_xy_norm"]],
            on="GEOID", suffixes=("_base", "_corr"),
        )
        d["f_xy_norm_delta"] = d["f_xy_norm_corr"] - d["f_xy_norm_base"]
        d = gpd.GeoDataFrame(d, geometry="geometry", crs=districts.crs)
        plot_choropleth(
            d, "f_xy_norm_delta",
            title="Δ f_xy_norm (corrected − baseline); + = hotter after correction",
            out_path=os.path.join(FIG_OUT, "ch_map_f_xy_norm_delta.png"),
            dpi=FIG_DPI, cmap="RdBu_r", boundary_overlay=districts,
        )

    # ── (a2) per-district delta maps, into each district's own folder ─────────
    make_district_delta_figures(cbg_full, grid_full, districts)

    ok = summary[summary["status"] == "fit"]

    # ── (b) corrected a_0 district map ────────────────────────────────────────
    fixed = ok[ok["spec"] == "offset_fixed"]
    if len(fixed):
        plot_choropleth(
            gpd.GeoDataFrame(fixed, geometry="geometry", crs=districts.crs),
            "a_0_mean",
            title="Corrected background level a_0 (offset_fixed) — citywide comparable",
            out_path=os.path.join(FIG_OUT, "ch_map_a_0_offset_fixed.png"),
            dpi=FIG_DPI,
        )

    # ── (c) alpha / excitation-share shift ────────────────────────────────────
    piv = ok.pivot_table(index="dist_name", columns="spec",
                         values=["alpha_mean", "prop_excitation_mean", "a_0_mean"])
    if ("alpha_mean", "baseline") in piv.columns and ("alpha_mean", "offset_fixed") in piv.columns:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        for ax, stat in zip(axes, ["alpha_mean", "prop_excitation_mean"]):
            x = piv[(stat, "baseline")]; y = piv[(stat, "offset_fixed")]
            ax.scatter(x, y)
            lims = [min(x.min(), y.min()), max(x.max(), y.max())]
            ax.plot(lims, lims, "k--", lw=1)
            for name in piv.index:
                ax.annotate(name, (x[name], y[name]), fontsize=6, alpha=0.8)
            ax.set_xlabel(f"{stat} (baseline)"); ax.set_ylabel(f"{stat} (offset_fixed)")
            ax.set_title(stat)
        fig.suptitle("Self-excitation before vs after reporting correction")
        fig.savefig(os.path.join(FIG_OUT, "ch_scatter_excitation_shift.png"),
                    dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)

    # ── (d) covariate-weight shifts (mean across districts) ──────────────────
    w_cols = [c for c in ok.columns if c.startswith("w_mean_")]
    if w_cols:
        wb = ok[ok["spec"] == "baseline"][w_cols].mean()
        wc = ok[ok["spec"] == "offset_fixed"][w_cols].mean()
        names = [c.replace("w_mean_", "") for c in w_cols]
        ypos = np.arange(len(names))
        fig, ax = plt.subplots(figsize=(8, 0.35 * len(names) + 2))
        ax.scatter(wb.values, ypos, label="baseline", zorder=3)
        ax.scatter(wc.values, ypos, label="offset_fixed", zorder=3)
        for i in ypos:
            ax.plot([wb.values[i], wc.values[i]], [i, i], color="gray", lw=1, zorder=2)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(ypos); ax.set_yticklabels(names, fontsize=7)
        ax.set_xlabel("posterior-mean weight (avg over districts)")
        ax.set_title("Covariate weights: baseline vs corrected")
        ax.legend()
        fig.savefig(os.path.join(FIG_OUT, "ch_dumbbell_cov_weights.png"),
                    dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)

    # ── (e) offset_coef forest plot (headline validation) ─────────────────────
    free = ok[(ok["spec"] == "offset_free")].dropna(subset=["offset_coef_mean"])
    if len(free):
        free = free.sort_values("offset_coef_mean")
        ypos = np.arange(len(free))
        fig, ax = plt.subplots(figsize=(7, 0.35 * len(free) + 2))
        ax.errorbar(
            free["offset_coef_mean"], ypos,
            xerr=[free["offset_coef_mean"] - free["offset_coef_lo"],
                  free["offset_coef_hi"] - free["offset_coef_mean"]],
            fmt="o", capsize=3,
        )
        ax.axvline(1.0, color="red", ls="--", lw=1, label="γ = 1 (pure thinning)")
        ax.set_yticks(ypos); ax.set_yticklabels(free["dist_name"], fontsize=7)
        ax.set_xlabel("offset coefficient γ (posterior mean, 94% CI)")
        ax.set_title("Free offset coefficient by district — γ≈1 validates the correction")
        ax.legend()
        fig.savefig(os.path.join(FIG_OUT, "ch_forest_offset_coef.png"),
                    dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  offset_coef: mean of district means = "
              f"{free['offset_coef_mean'].mean():.3f}")


# =============================================================================
# SECTION 7 - Main execution
# =============================================================================

if __name__ == "__main__":
    main()
