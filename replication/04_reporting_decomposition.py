"""
04_reporting_decomposition.py
=============================
Separate the 311 *reporting propensity* from the *true illegal-dumping rate*
using a two-equation shared-latent-factor count model, fit independently per
planning district (same analysis geography as 03a_analysis_by_division.py)
on CBG-level counts.

The idea (whiteboard: W_R -> {NT, T},  W_T -> {T only})
-------------------------------------------------------
Observed dumping reports = (true dumping) x (willingness to report).  A second
report stream that shares the willingness component but carries no dumping
signal identifies the split: non-trash (NT) 311 reports are filed by the same
residents through the same channel, so they pin down the latent reporting
propensity W_R; whatever variation is left in trash (T) reports after
crediting W_R is attributed to W_T — the reporting-corrected dumping rate.

Model (per district d, CBG i, counts pooled over STUDY_YEARS)
--------------------------------------------------------------
  NT_i ~ Poisson( E_i * exp( a_NT + w_R,i ) )            eq.1  NT pins down W_R
  T_i  ~ Poisson( E_i * exp( a_T + beta*w_R,i + w_T,i )) eq.2  leftover = W_T

  w_R,i ~ Normal(0, sigma_R^2),  w_T,i ~ Normal(0, sigma_T^2)
  beta  ~ LogNormal(0, 0.5),     sigma_R, sigma_T ~ HalfNormal(1)

  E_i = exposure offset (SafeGraph avg monthly stop counts, coefficient
        fixed at 1) — converts counts to rates.

Key identifying assumption: cross-CBG variation in the NT categories reflects
reporting behaviour, not variation in the true rate of those problems.  The
main spec therefore uses infrastructure-failure categories (street lights
burn out roughly uniformly everywhere); the all-non-dumping spec is the
robustness check on that assumption.

Cross-district scalar: under the assumption that true NT incidence per unit
exposure is uniform citywide, (a_NT - a_ref) is the district's log propensity
deviation from the citywide NT rate a_ref = log(sum NT / sum E), so
  corrected_level = a_T - beta * (a_NT - a_ref)
(computed per posterior sample) is the log reporting-corrected dumping rate
per unit exposure, comparable across districts.  The a_ref centring matters:
beta differs by district, so without it the scalar would be driven
mechanically by beta rather than by dumping.

Downstream hook (not implemented here): beta * w_R,i is the estimated log
reporting-thinning log p(s); it can enter the per-district Cox-Hawkes fits
(03a) as a fixed offset in the background intensity, making each district's
f_s a reporting-corrected dumping surface.

Required inputs
---------------
  data/311/{year}/public_cases_fc.shp      OpenDataPhilly 311 records
  data/tiger/tl_2023_42_bg/tl_2023_42_bg.shp  TIGER 2023 CBG boundaries
  data/Planning_Districts.geojson          18 planning districts (dist_name)
  output/cov_cbg.geojson                   from 02_covariates.py
                                           (exposure alloc_avg_s_cnt +
                                            reporting_rate validation column)

Optional validation inputs (skipped gracefully if absent)
---------------------------------------------------------
  data/CSVFiles/litter_surveys_raw_20250730.csv
  data/PPR_BaseFiles_IllegalDumpingModel.gdb (layer: PPR_Ops_Dumping_Cleanup_Complete)

Outputs  (written to replication/output/)
-----------------------------------------
  reporting_decomp_cbg.geojson              per CBG: w_R / w_T posterior mean
                                            & sd + within-district min-max
                                            normalised values, per NT spec
  reporting_decomp_district_summary.csv     per district x spec: a_NT, a_T,
  reporting_decomp_district_summary.geojson beta, sigmas, corrected_level,
                                            raw level, rankings, fit status
  figures/reporting_decomposition/
    rd_map_w_T_norm_cbg.png       HEADLINE: whole-Philly CBG mosaic of W_T,
                                  min-max normalised within district (03a style)
    rd_map_w_R_norm_cbg.png       same mosaic for reporting propensity W_R
    rd_map_{stat}.png             district choropleths (corrected level, raw
                                  level, propensity level, beta)
    rd_rank_raw_vs_corrected.png  district ranking shift, raw vs corrected
    rd_spec_robustness.png        spec A vs spec B w_T scatter + correlation
    rd_beta_forest.png            beta posterior by district
    rd_ppc_fit.png                observed vs posterior-mean fitted counts
    rd_validation_reporting_rate.png  w_R vs existing reporting_rate covariate
"""

# =============================================================================
# SECTION 0 – Imports and global paths
# =============================================================================

import os
import time
import warnings

os.environ["JAX_PLATFORM_NAME"] = "cpu"

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS

warnings.filterwarnings("ignore")

REPL_DIR = os.path.dirname(os.path.abspath(__file__))   # replication/
DATA     = os.path.join(REPL_DIR, "data")               # replication/data/
OUT      = os.path.join(REPL_DIR, "output")             # replication/output/
FIG_OUT  = os.path.join(OUT, "figures", "reporting_decomposition")
os.makedirs(FIG_OUT, exist_ok=True)

print("=" * 60)
print("04_reporting_decomposition.py — W_R / W_T decomposition")
print("=" * 60)


# =============================================================================
# SECTION 1 – Configuration
# =============================================================================

STUDY_YEARS = [2021, 2022, 2023, 2024, 2025]

# ── NT (non-trash) report specifications ─────────────────────────────────────
# Spec A (main): infrastructure-failure categories whose *true* incidence is
# plausibly uniform across neighbourhoods, so cross-CBG variation in their
# counts reflects willingness to report.
NT_CURATED = [
    "Street Light Outage",
    "Alley Light Outage",
    "Traffic Signal Emergency",
    "Stop Sign Repair",
    "Dead Animal in Street",
]

# Spec B (robustness): all non-dumping 311, excluding phone-based information
# requests (not geolocated incidents) and trash-adjacent categories that
# overlap the outcome.
NT_ALL_EXCLUDE = [
    "Illegal Dumping",
    "Information Request",
    "Rubbish/Recyclable Material Collection",
    "Sanitation Violation",
    "Dumpster Violation",
    "Sanitation / Dumpster Violation",
]

# {spec key: column name of the NT count}
NT_SPECS = {"a": "NT_curated", "b": "NT_all"}
MAIN_SPEC = "a"     # spec used for headline maps / district choropleths

# ── Analysis geography ────────────────────────────────────────────────────────
PLANNING_DISTRICTS_PATH  = os.path.join(DATA, "Planning_Districts.geojson")
PLANNING_DISTRICT_ID_COL = "dist_name"
UNIT_SUBSET = None          # e.g. ["Central", "South"] to fit a subset; None = all 18

# ── Exposure ──────────────────────────────────────────────────────────────────
EXPOSURE_COL = "alloc_avg_s_cnt"    # SafeGraph avg monthly stop counts (cov_cbg)
# Fallback if EXPOSURE_COL is absent from cov_cbg.geojson: population
# (pop_density [per km2] x CBG area).  One exposure definition is used for
# ALL CBGs — mixing definitions across CBGs would distort the within-district
# comparison.

# ── MCMC settings ─────────────────────────────────────────────────────────────
MCMC_WARMUP  = 1000
MCMC_SAMPLES = 1000
MCMC_CHAINS  = 1
MCMC_SEED    = 0

FIG_DPI = 300

# ── File paths ────────────────────────────────────────────────────────────────
TIGER_BG_PATH  = os.path.join(DATA, "tiger", "tl_2023_42_bg", "tl_2023_42_bg.shp")
COV_CBG_PATH   = os.path.join(OUT, "cov_cbg.geojson")
LITTER_CSV_PATH  = os.path.join(DATA, "CSVFiles", "litter_surveys_raw_20250730.csv")
CLEANUP_GDB_PATH = os.path.join(DATA, "PPR_BaseFiles_IllegalDumpingModel.gdb")


# =============================================================================
# SECTION 2 – Build CBG-level count table
# =============================================================================

def build_cbg_counts() -> gpd.GeoDataFrame:
    """
    Count T (Illegal Dumping), NT_curated and NT_all 311 reports per 2023 CBG,
    pooled over STUDY_YEARS.  Returns the Philadelphia CBG GeoDataFrame
    (EPSG:26918) with count columns T, NT_curated, NT_all.
    """
    print("\n[2] Counting 311 reports per CBG ...")

    pa_2023      = gpd.read_file(TIGER_BG_PATH)
    philly_26918 = pa_2023[pa_2023["COUNTYFP"] == "101"].to_crs(26918)
    philly_26918["GEOID"] = philly_26918["GEOID"].astype(str)
    print(f"  {len(philly_26918)} Philadelphia CBGs loaded.")

    parts = []
    for yr in STUDY_YEARS:
        shp_path = os.path.join(DATA, "311", str(yr), "public_cases_fc.shp")
        if not os.path.exists(shp_path):
            print(f"  [SKIP] {shp_path} not found.")
            continue
        gdf = gpd.read_file(shp_path, columns=["service_na"])
        gdf = gdf[gdf["geometry"].notna()].to_crs(26918)
        print(f"  {yr}: {len(gdf):,} geolocated 311 records")
        parts.append(gdf)

    if not parts:
        raise FileNotFoundError(
            "No public 311 shapefiles found under data/311/{year}/."
        )
    all_311 = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=26918)

    # classify each record into the count groups (curated is a subset of all)
    svc = all_311["service_na"]
    all_311["is_T"]      = svc == "Illegal Dumping"
    all_311["is_NT_cur"] = svc.isin(NT_CURATED)
    all_311["is_NT_all"] = ~svc.isin(NT_ALL_EXCLUDE)

    joined = all_311.sjoin(
        philly_26918[["GEOID", "geometry"]], how="inner", predicate="within"
    )
    counts = joined.groupby("GEOID")[["is_T", "is_NT_cur", "is_NT_all"]].sum()
    counts.columns = ["T_count", "NT_curated", "NT_all"]

    out = philly_26918[["GEOID", "geometry"]].merge(
        counts.reset_index(), on="GEOID", how="left"
    )
    for c in ["T_count", "NT_curated", "NT_all"]:
        out[c] = out[c].fillna(0).astype(int)

    print(f"  Totals — T: {out['T_count'].sum():,}  "
          f"NT_curated: {out['NT_curated'].sum():,}  "
          f"NT_all: {out['NT_all'].sum():,}")
    return gpd.GeoDataFrame(out, geometry="geometry", crs=26918)


# =============================================================================
# SECTION 3 – Exposure and district assignment
# =============================================================================

def attach_exposure_and_district(cbg: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Merge the exposure column from cov_cbg.geojson (population fallback when
    the foot-traffic column is entirely absent) and assign each CBG to the
    planning district containing its representative point.
    """
    print("\n[3] Attaching exposure and planning district ...")

    cov = gpd.read_file(COV_CBG_PATH)
    cov["GEOID"] = cov["GEOID"].astype(str)

    if EXPOSURE_COL in cov.columns:
        expo = cov[["GEOID", EXPOSURE_COL]].rename(columns={EXPOSURE_COL: "exposure"})
        print(f"  Exposure: {EXPOSURE_COL} (SafeGraph foot traffic)")
    else:
        # population fallback: pop_density [per km2] x CBG area [km2]
        cov_pop = cov[["GEOID", "pop_density"]].merge(
            cbg[["GEOID"]].assign(area_km2=cbg.geometry.area / 1e6),
            on="GEOID", how="right",
        )
        cov_pop["exposure"] = cov_pop["pop_density"] * cov_pop["area_km2"]
        expo = cov_pop[["GEOID", "exposure"]]
        print(f"  Exposure: population (fallback — '{EXPOSURE_COL}' absent from cov_cbg)")

    cbg = cbg.merge(expo, on="GEOID", how="left")

    # keep the existing reporting_rate covariate for validation, if present
    if "reporting_rate" in cov.columns:
        cbg = cbg.merge(cov[["GEOID", "reporting_rate"]], on="GEOID", how="left")

    districts = gpd.read_file(PLANNING_DISTRICTS_PATH).to_crs(26918)
    rep_pts = cbg.copy()
    rep_pts["geometry"] = cbg.geometry.representative_point()
    assign = gpd.sjoin(
        rep_pts[["GEOID", "geometry"]],
        districts[[PLANNING_DISTRICT_ID_COL, "geometry"]],
        how="left", predicate="within",
    )[["GEOID", PLANNING_DISTRICT_ID_COL]]
    cbg = cbg.merge(assign, on="GEOID", how="left")

    n_no_dist = cbg[PLANNING_DISTRICT_ID_COL].isna().sum()
    n_no_expo = (~(cbg["exposure"] > 0)).sum()
    print(f"  CBGs without district: {n_no_dist}   "
          f"CBGs with missing/zero exposure (excluded from fits): {n_no_expo}")
    return cbg, districts


# =============================================================================
# SECTION 4 – Two-equation model and per-district fitting
# =============================================================================

def two_eq_model(log_E, nt_obs=None, t_obs=None, a_ref=0.0):
    """
    Joint NT/T Poisson model with a shared latent reporting-propensity field.
    Non-centred parameterisation for NUTS stability.

    a_ref : citywide log NT rate per unit exposure, log(sum NT / sum E) over
        all fitted CBGs — the reference point for the cross-district
        corrected_level scalar (see module docstring).
    """
    n = log_E.shape[0]

    a_nt    = numpyro.sample("a_NT",    dist.Normal(0.0, 5.0))
    a_t     = numpyro.sample("a_T",     dist.Normal(0.0, 5.0))
    beta    = numpyro.sample("beta",    dist.LogNormal(0.0, 0.5))
    sigma_r = numpyro.sample("sigma_R", dist.HalfNormal(1.0))
    sigma_t = numpyro.sample("sigma_T", dist.HalfNormal(1.0))

    with numpyro.plate("cbg", n):
        z_r = numpyro.sample("z_R", dist.Normal(0.0, 1.0))
        z_t = numpyro.sample("z_T", dist.Normal(0.0, 1.0))

    w_r = numpyro.deterministic("w_R", z_r * sigma_r)
    w_t = numpyro.deterministic("w_T", z_t * sigma_t)

    # cross-district scalar: log corrected dumping rate per unit exposure.
    # (a_nt - a_ref) = district log propensity deviation from the citywide
    # NT rate; centring on a_ref keeps the scalar from being driven
    # mechanically by cross-district differences in beta.
    numpyro.deterministic("corrected_level", a_t - beta * (a_nt - a_ref))

    with numpyro.plate("obs", n):
        numpyro.sample("NT", dist.Poisson(jnp.exp(log_E + a_nt + w_r)),
                       obs=nt_obs)
        numpyro.sample("T",  dist.Poisson(jnp.exp(log_E + a_t + beta * w_r + w_t)),
                       obs=t_obs)


def fit_district(sub: pd.DataFrame, nt_col: str, seed: int,
                 a_ref: float = 0.0) -> dict:
    """
    Fit the two-equation model for one district's CBGs.

    Parameters
    ----------
    sub : DataFrame with columns GEOID, T_count, {nt_col}, exposure (all > 0)
    nt_col : NT count column for this spec
    seed : PRNG seed
    a_ref : citywide log NT rate per unit exposure (see two_eq_model)

    Returns
    -------
    dict with keys:
      cbg     – DataFrame (GEOID, w_R_mean/sd, w_T_mean/sd, T_fitted, NT_fitted)
      scalars – dict of district-level posterior summaries
    """
    log_E = jnp.asarray(np.log(sub["exposure"].to_numpy(float)))
    nt    = jnp.asarray(sub[nt_col].to_numpy(int))
    t     = jnp.asarray(sub["T_count"].to_numpy(int))

    kernel = NUTS(two_eq_model, target_accept_prob=0.9)
    mcmc = MCMC(
        kernel, num_warmup=MCMC_WARMUP, num_samples=MCMC_SAMPLES,
        num_chains=MCMC_CHAINS, progress_bar=False,
    )
    mcmc.run(jax.random.PRNGKey(seed), log_E=log_E, nt_obs=nt, t_obs=t,
             a_ref=a_ref)
    s = mcmc.get_samples()

    w_r = np.asarray(s["w_R"])          # (draws, n)
    w_t = np.asarray(s["w_T"])
    a_nt, a_t = np.asarray(s["a_NT"]), np.asarray(s["a_T"])
    beta = np.asarray(s["beta"])

    # posterior-mean fitted counts (for the PPC figure)
    log_E_np = np.asarray(log_E)
    nt_fit = np.exp(log_E_np[None, :] + a_nt[:, None] + w_r).mean(axis=0)
    t_fit  = np.exp(log_E_np[None, :] + a_t[:, None]
                    + beta[:, None] * w_r + w_t).mean(axis=0)

    cbg = pd.DataFrame({
        "GEOID":     sub["GEOID"].to_numpy(),
        "w_R_mean":  w_r.mean(axis=0), "w_R_sd": w_r.std(axis=0),
        "w_T_mean":  w_t.mean(axis=0), "w_T_sd": w_t.std(axis=0),
        "NT_fitted": nt_fit, "T_fitted": t_fit,
    })

    corr = np.asarray(s["corrected_level"])
    scalars = {
        "a_NT_mean": a_nt.mean(),  "a_T_mean": a_t.mean(),
        "beta_mean": beta.mean(),
        "beta_lo":   np.quantile(beta, 0.03), "beta_hi": np.quantile(beta, 0.97),
        "sigma_R_mean": np.asarray(s["sigma_R"]).mean(),
        "sigma_T_mean": np.asarray(s["sigma_T"]).mean(),
        "corrected_level_mean": corr.mean(),
        "corrected_level_sd":   corr.std(),
    }
    return {"cbg": cbg, "scalars": scalars}


def run_all_districts(cbg: gpd.GeoDataFrame, districts: gpd.GeoDataFrame):
    """
    Loop over districts x NT specs, fit, and assemble:
      cbg_results  – cbg GeoDataFrame with w_R/w_T columns per spec (suffix _a/_b)
      summary      – district-level GeoDataFrame, one row per district x spec
    """
    print("\n[4] Fitting the two-equation model per district ...")

    dist_names = districts[PLANNING_DISTRICT_ID_COL].tolist()
    if UNIT_SUBSET is not None:
        dist_names = [d for d in dist_names if d in UNIT_SUBSET]

    fit_mask = (cbg["exposure"] > 0) & cbg[PLANNING_DISTRICT_ID_COL].notna()

    # citywide reference NT rate per spec: a_ref = log(sum NT / sum E) over
    # all fitted CBGs — anchors the cross-district corrected_level scalar
    a_refs = {
        spec: float(np.log(
            cbg.loc[fit_mask, nt_col].sum() / cbg.loc[fit_mask, "exposure"].sum()
        ))
        for spec, nt_col in NT_SPECS.items()
    }
    for spec, nt_col in NT_SPECS.items():
        print(f"  Citywide reference a_ref (spec {spec}, {nt_col}): "
              f"{a_refs[spec]:.3f}")

    summary_rows = []
    cbg_results = cbg.copy()
    for spec, nt_col in NT_SPECS.items():
        for col in ["w_R_mean", "w_R_sd", "w_T_mean", "w_T_sd",
                    "NT_fitted", "T_fitted"]:
            cbg_results[f"{col}_{spec}"] = np.nan

    for di, dname in enumerate(dist_names):
        sub = cbg[fit_mask & (cbg[PLANNING_DISTRICT_ID_COL] == dname)]
        geom = districts.loc[
            districts[PLANNING_DISTRICT_ID_COL] == dname, "geometry"
        ].iloc[0]

        for si, (spec, nt_col) in enumerate(NT_SPECS.items()):
            print(f"  [{dname} | spec {spec.upper()} ({nt_col})] "
                  f"n_cbg={len(sub)}, T={sub['T_count'].sum():,}, "
                  f"NT={sub[nt_col].sum():,} ... ", end="", flush=True)
            t0 = time.time()
            row = {
                PLANNING_DISTRICT_ID_COL: dname, "spec": spec,
                "nt_col": nt_col, "n_cbg": len(sub),
                "T_total": int(sub["T_count"].sum()),
                "NT_total": int(sub[nt_col].sum()),
                # raw (uncorrected) district log dumping rate per unit exposure
                "raw_level": float(np.log(
                    max(sub["T_count"].sum(), 0.5) / sub["exposure"].sum()
                )),
                "geometry": geom,
            }
            try:
                if len(sub) < 5:
                    raise ValueError(f"too few CBGs with exposure ({len(sub)})")
                res = fit_district(sub, nt_col, seed=MCMC_SEED + 100 * di + si,
                                   a_ref=a_refs[spec])
                row.update(res["scalars"])
                row["status"], row["error"] = "fit", None

                idx = cbg_results["GEOID"].isin(res["cbg"]["GEOID"])
                merged = res["cbg"].set_index("GEOID")
                for col in ["w_R_mean", "w_R_sd", "w_T_mean", "w_T_sd",
                            "NT_fitted", "T_fitted"]:
                    cbg_results.loc[idx, f"{col}_{spec}"] = (
                        cbg_results.loc[idx, "GEOID"].map(merged[col]).values
                    )
            except Exception as exc:
                row["status"], row["error"] = "failed", str(exc)
                print(f"[FAILED] {exc}")
            else:
                print(f"done ({time.time() - t0:.1f}s, "
                      f"beta={row.get('beta_mean', float('nan')):.2f})")
            row["runtime_sec"] = round(time.time() - t0, 1)
            summary_rows.append(row)

    summary = gpd.GeoDataFrame(
        pd.DataFrame(summary_rows), geometry="geometry", crs=districts.crs
    )

    # within-district min-max normalisation (03a style) for the mosaic maps
    for spec in NT_SPECS:
        for base in ["w_T", "w_R"]:
            col, norm_col = f"{base}_mean_{spec}", f"{base}_norm_{spec}"
            grp = cbg_results.groupby(PLANNING_DISTRICT_ID_COL)[col]
            vmin, vmax = grp.transform("min"), grp.transform("max")
            rng = vmax - vmin
            cbg_results[norm_col] = np.where(
                rng > 0, (cbg_results[col] - vmin) / rng, 0.5
            )
            cbg_results.loc[cbg_results[col].isna(), norm_col] = np.nan

    n_fit = int((summary["status"] == "fit").sum())
    print(f"\n  {n_fit}/{len(summary)} district-spec fits succeeded.")
    return cbg_results, summary


# =============================================================================
# SECTION 5 – Figures
# =============================================================================

def plot_choropleth(gdf, column, title, out_path, dpi=300,
                    cmap="viridis", boundary_overlay=None):
    """Single choropleth PNG (same helper style as 03a)."""
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


def make_figures(cbg_results, summary, districts):
    print("\n[5] Rendering figures ...")
    m = MAIN_SPEC

    # ── Headline mosaics: within-district normalised W_T and W_R ─────────────
    plot_choropleth(
        cbg_results, f"w_T_norm_{m}",
        title="Reporting-corrected dumping W_T (normalized within district) by CBG",
        out_path=os.path.join(FIG_OUT, "rd_map_w_T_norm_cbg.png"),
        dpi=FIG_DPI, boundary_overlay=districts,
    )
    plot_choropleth(
        cbg_results, f"w_R_norm_{m}",
        title="Reporting propensity W_R (normalized within district) by CBG",
        out_path=os.path.join(FIG_OUT, "rd_map_w_R_norm_cbg.png"),
        dpi=FIG_DPI, boundary_overlay=districts,
    )

    # ── District choropleths (main spec) ──────────────────────────────────────
    summ_m = summary[summary["spec"] == m].copy()
    for col, title in [
        ("corrected_level_mean", "Reporting-corrected dumping level (a_T − β·(a_NT − ā))"),
        ("raw_level",            "Raw dumping level  log(ΣT / ΣE)"),
        ("a_NT_mean",            "Reporting propensity level (a_NT)"),
        ("beta_mean",            "Shared-factor loading β"),
    ]:
        plot_choropleth(
            summ_m, col, title=f"{title} by planning district",
            out_path=os.path.join(FIG_OUT, f"rd_map_{col}.png"), dpi=FIG_DPI,
        )

    # ── Raw vs corrected district ranking ─────────────────────────────────────
    ok = summ_m[summ_m["status"] == "fit"].copy()
    if len(ok) > 1:
        ok["rank_raw"]  = ok["raw_level"].rank(ascending=False)
        ok["rank_corr"] = ok["corrected_level_mean"].rank(ascending=False)
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(ok["rank_raw"], ok["rank_corr"], color="tab:blue", zorder=3)
        for _, r in ok.iterrows():
            ax.annotate(r[PLANNING_DISTRICT_ID_COL],
                        (r["rank_raw"], r["rank_corr"]),
                        fontsize=8, xytext=(4, 4), textcoords="offset points")
        lim = [0, len(ok) + 1]
        ax.plot(lim, lim, ls="--", color="gray", zorder=1)
        ax.set_xlabel("District rank, raw dumping reports (1 = highest)")
        ax.set_ylabel("District rank, reporting-corrected (1 = highest)")
        ax.set_title("Does correcting for reporting propensity reshuffle districts?")
        out = os.path.join(FIG_OUT, "rd_rank_raw_vs_corrected.png")
        fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {out}")

    # ── Spec A vs spec B robustness ────────────────────────────────────────────
    both = cbg_results[["w_T_mean_a", "w_T_mean_b"]].dropna()
    if len(both) > 2:
        r = both["w_T_mean_a"].corr(both["w_T_mean_b"])
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(both["w_T_mean_a"], both["w_T_mean_b"], s=8, alpha=0.4)
        ax.set_xlabel("w_T — spec A (curated NT categories)")
        ax.set_ylabel("w_T — spec B (all non-dumping 311)")
        ax.set_title(f"Robustness of W_T to NT definition (r = {r:.3f}, n = {len(both):,})")
        ax.axline((0, 0), slope=1, ls="--", color="gray")
        out = os.path.join(FIG_OUT, "rd_spec_robustness.png")
        fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {out}")
        print(f"  [robustness] corr(w_T spec A, w_T spec B) = {r:.3f}")

    # ── Beta forest plot ───────────────────────────────────────────────────────
    ok_beta = summary[(summary["status"] == "fit")].copy()
    if not ok_beta.empty:
        fig, ax = plt.subplots(figsize=(7, 9))
        offsets, colors = {"a": -0.15, "b": 0.15}, {"a": "tab:blue", "b": "tab:orange"}
        names = ok_beta[PLANNING_DISTRICT_ID_COL].unique().tolist()
        for spec in NT_SPECS:
            ss = ok_beta[ok_beta["spec"] == spec]
            y = [names.index(n) + offsets[spec] for n in ss[PLANNING_DISTRICT_ID_COL]]
            ax.errorbar(
                ss["beta_mean"], y,
                xerr=[ss["beta_mean"] - ss["beta_lo"], ss["beta_hi"] - ss["beta_mean"]],
                fmt="o", color=colors[spec], label=f"spec {spec.upper()}",
                markersize=4, capsize=2, lw=1,
            )
        ax.axvline(1.0, ls="--", color="gray")
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("β  (propensity loading into trash reports, 94% CI)")
        ax.set_title("Shared-factor loading β by district")
        ax.legend()
        out = os.path.join(FIG_OUT, "rd_beta_forest.png")
        fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {out}")

    # ── Posterior-predictive fit check ────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, obs_col, fit_col, label in [
        (axes[0], "T_count",           f"T_fitted_{m}",  "Trash (T)"),
        (axes[1], NT_SPECS[m] if NT_SPECS[m] in cbg_results.columns else "NT_curated",
                  f"NT_fitted_{m}", "Non-trash (NT)"),
    ]:
        d = cbg_results[[obs_col, fit_col]].dropna()
        ax.scatter(d[obs_col] + 0.5, d[fit_col] + 0.5, s=8, alpha=0.4)
        ax.set_xscale("log"); ax.set_yscale("log")
        lims = [0.4, max((d[obs_col].max() if len(d) else 1), 1) * 2]
        ax.plot(lims, lims, ls="--", color="gray")
        ax.set_xlabel(f"Observed {label} count (+0.5, log)")
        ax.set_ylabel(f"Posterior-mean fitted (+0.5, log)")
        ax.set_title(label)
    fig.suptitle(f"Posterior-predictive fit — spec {m.upper()}")
    out = os.path.join(FIG_OUT, "rd_ppc_fit.png")
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")


# =============================================================================
# SECTION 6 – Validation checks
# =============================================================================

def run_validation(cbg_results: gpd.GeoDataFrame):
    """
    (i)  w_R should correlate positively with the existing reporting_rate
         covariate (non-dumping 311 / foot traffic from 02_covariates.py).
    (ii) w_T should track ground-truth trash measures (litter survey, PPR
         cleanups) at least as well as raw T counts — skipped when those
         proprietary files are absent.
    """
    print("\n[6] Validation checks ...")
    m = MAIN_SPEC

    # ── (i) against the existing reporting_rate covariate ────────────────────
    if "reporting_rate" in cbg_results.columns:
        d = cbg_results[[f"w_R_mean_{m}", "reporting_rate"]].dropna()
        d = d[np.isfinite(d["reporting_rate"])]
        if len(d) > 2:
            r_p = d[f"w_R_mean_{m}"].corr(d["reporting_rate"])
            r_s = d[f"w_R_mean_{m}"].corr(d["reporting_rate"], method="spearman")
            print(f"  corr(w_R, reporting_rate covariate): "
                  f"pearson={r_p:.3f}, spearman={r_s:.3f}  (n={len(d):,}) "
                  f"— expected: strongly positive")
            fig, ax = plt.subplots(figsize=(7, 7))
            ax.scatter(d["reporting_rate"], d[f"w_R_mean_{m}"], s=8, alpha=0.4)
            ax.set_xscale("log")
            ax.set_xlabel("Existing reporting_rate covariate (311 / stops, log)")
            ax.set_ylabel("Estimated reporting propensity w_R")
            ax.set_title(f"w_R vs. reporting_rate (spearman r = {r_s:.3f})")
            out = os.path.join(FIG_OUT, "rd_validation_reporting_rate.png")
            fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved → {out}")
    else:
        print("  [skip] reporting_rate column absent from cov_cbg.geojson.")

    # ── (ii) against litter survey ────────────────────────────────────────────
    if os.path.exists(LITTER_CSV_PATH):
        litter = pd.read_csv(LITTER_CSV_PATH).dropna(subset=["x", "y"])
        litter_gdf = gpd.GeoDataFrame(
            litter, geometry=gpd.points_from_xy(litter["x"], litter["y"]),
            crs=4326,
        ).to_crs(26918)
        joined = gpd.sjoin(
            litter_gdf, cbg_results[["GEOID", "geometry"]],
            how="inner", predicate="within",
        )
        litter_cbg = joined.groupby("GEOID")["Litter Rating"].mean().rename("litter_score")
        d = cbg_results.set_index("GEOID").join(litter_cbg, how="inner")
        d["raw_T_rate"] = d["T_count"] / d["exposure"]
        d = d[[f"w_T_mean_{m}", "raw_T_rate", "litter_score"]].dropna()
        if len(d) > 2:
            r_wt  = d[f"w_T_mean_{m}"].corr(d["litter_score"], method="spearman")
            r_raw = d["raw_T_rate"].corr(d["litter_score"], method="spearman")
            print(f"  spearman corr with litter score (n={len(d)}): "
                  f"w_T={r_wt:.3f}  vs  raw T rate={r_raw:.3f} "
                  f"— w_T should be ≥ raw")
    else:
        print("  [skip] Litter survey file not found — ground-truth check (i) skipped.")

    # ── (ii) against PPR cleanup operations ───────────────────────────────────
    if os.path.exists(CLEANUP_GDB_PATH):
        cleanup = gpd.read_file(
            CLEANUP_GDB_PATH, layer="PPR_Ops_Dumping_Cleanup_Complete"
        ).to_crs(26918)
        joined = gpd.sjoin(
            cleanup, cbg_results[["GEOID", "geometry"]],
            how="inner", predicate="within",
        )
        clean_cbg = joined.groupby("GEOID").size().rename("cleanup_count")
        d = cbg_results.set_index("GEOID").join(clean_cbg, how="inner")
        d["raw_T_rate"] = d["T_count"] / d["exposure"]
        d = d[[f"w_T_mean_{m}", "raw_T_rate", "cleanup_count"]].dropna()
        if len(d) > 2:
            r_wt  = d[f"w_T_mean_{m}"].corr(d["cleanup_count"], method="spearman")
            r_raw = d["raw_T_rate"].corr(d["cleanup_count"], method="spearman")
            print(f"  spearman corr with PPR cleanups (n={len(d)}): "
                  f"w_T={r_wt:.3f}  vs  raw T rate={r_raw:.3f} "
                  f"— w_T should be ≥ raw")
    else:
        print("  [skip] PPR cleanup GDB not found — ground-truth check (ii) skipped.")


# =============================================================================
# SECTION 7 – Main
# =============================================================================

if __name__ == "__main__":
    cbg_counts = build_cbg_counts()
    cbg, districts = attach_exposure_and_district(cbg_counts)

    cbg_results, summary = run_all_districts(cbg, districts)

    # ── Export ────────────────────────────────────────────────────────────────
    print("\n[7] Exporting results ...")
    out_cbg = os.path.join(OUT, "reporting_decomp_cbg.geojson")
    cbg_results.to_crs(4326).to_file(out_cbg, driver="GeoJSON")
    print(f"  Saved → {out_cbg}")

    out_geo = os.path.join(OUT, "reporting_decomp_district_summary.geojson")
    out_csv = os.path.join(OUT, "reporting_decomp_district_summary.csv")
    summary.to_crs(4326).to_file(out_geo, driver="GeoJSON")
    summary.drop(columns="geometry").to_csv(out_csv, index=False)
    print(f"  Saved → {out_geo}")
    print(f"  Saved → {out_csv}")

    make_figures(cbg_results, summary, districts)
    run_validation(cbg_results)

    # ── Summary ───────────────────────────────────────────────────────────────
    ok = summary[(summary["spec"] == MAIN_SPEC) & (summary["status"] == "fit")]
    print("\n" + "=" * 60)
    print(f"Summary (spec {MAIN_SPEC.upper()} — {NT_SPECS[MAIN_SPEC]})")
    print("=" * 60)
    if not ok.empty:
        print(f"  Districts fit: {len(ok)}/{summary['spec'].eq(MAIN_SPEC).sum()}")
        print(f"  beta range across districts: "
              f"[{ok['beta_mean'].min():.2f}, {ok['beta_mean'].max():.2f}]")
        top = ok.nlargest(5, "corrected_level_mean")
        print("  Top-5 districts by reporting-corrected dumping level:")
        for _, r in top.iterrows():
            print(f"    {r[PLANNING_DISTRICT_ID_COL]:<25s} "
                  f"corrected={r['corrected_level_mean']:.2f}  "
                  f"raw={r['raw_level']:.2f}  beta={r['beta_mean']:.2f}")
    print("\nOutputs written to:", OUT)
    print("Figures written to:", FIG_OUT)
    print("Done.")
