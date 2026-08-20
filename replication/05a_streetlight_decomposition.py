"""
05a_streetlight_decomposition.py
================================
Citywide street-light reporting decomposition — the PREREQUISITE for
05_corrected_hawkes_by_district.py. Documented in 05a_streetlight_citywide.md.

What it does:

  1. Fits ONE citywide two-equation decomposition over all Philadelphia CBGs —
     no per-district beta_d / a_NT,d, no cross-district anchor scalar. Every
     CBG's w_R / w_T is directly comparable citywide.
  2. Uses the street-pole inventory (data/Street_Poles.geojson) as the
     exposure in the NT equation: outage reports per POLE, not per
     foot-traffic unit. The identifying assumption becomes "true breakage per
     pole is uniform across the city".
  3. Validates that assumption: pole EDA, proxy-validation correlations, and
     two PSIP robustness checks (the LED conversion landed mostly inside the
     study window's tail, so true breakage fell non-uniformly during 2024-25):
       - pre-conversion time slice: refit with NT counts from 2021-2023 only;
       - LED-share sensitivity: corr(w_R, per-CBG LED / pre-2024-relit share).

Downstream, 05_corrected_hawkes_by_district.py reads this script's outputs
(reporting_decomp_cbg_sl.geojson + reporting_decomp_city_summary_sl.csv) and
builds the per-CBG reporting offset  beta * w_R_i  for its per-district
Cox-Hawkes fits. This script is self-contained: the CBG count / exposure /
district builders formerly imported from 04_reporting_decomposition.py (now
archived under archive/) are absorbed below.

Usage
-----
  python 05a_streetlight_decomposition.py

Inputs
------
  data/311/{2021..2025}/public_cases_fc.shp     geolocated 311 reports
  data/tiger/tl_2023_42_bg/tl_2023_42_bg.shp    2023 census block groups
  data/Planning_Districts.geojson               18 planning districts
  data/Street_Poles.geojson                     street-pole inventory
  output/cov_cbg.geojson                        exposure (alloc_avg_s_cnt)

Outputs
-------
  output/corrected_hawkes_streetlight/
    reporting_decomp_cbg_sl.geojson            citywide decomposition per CBG
    reporting_decomp_city_summary_sl.csv       specs sl / sl_pre24: beta etc.
    reporting_decomp_*_districtfit.*           archived pre-refactor outputs
    poles_cbg.geojson                          poles per CBG + LED/PECO shares
    pole_eda_*.csv                             EDA tables
    validation_summary.csv                     all validation correlations
  output/figures/corrected_hawkes_streetlight/*.png
"""

import os
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

REPL_DIR = os.path.dirname(os.path.abspath(__file__))
DATA     = os.path.join(REPL_DIR, "data")
OUT      = os.path.join(REPL_DIR, "output")

RUN_NAME = "corrected_hawkes_streetlight"
RUN_OUT  = os.path.join(OUT, RUN_NAME)
FIG_OUT  = os.path.join(OUT, "figures", RUN_NAME)
os.makedirs(RUN_OUT, exist_ok=True)
os.makedirs(FIG_OUT, exist_ok=True)

NT_STREETLIGHT = ["Street Light Outage"]

POLES_PATH   = os.path.join(DATA, "Street_Poles.geojson")
PSIP_CUTOFF  = pd.Timestamp("2024-01-01", tz="UTC")    # mass LED relighting starts
PRE24_YEARS  = [2021, 2022, 2023]                      # pre-conversion NT window

DECOMP_CBG_PATH   = os.path.join(RUN_OUT, "reporting_decomp_cbg_sl.geojson")
CITY_SUMMARY_PATH = os.path.join(RUN_OUT, "reporting_decomp_city_summary_sl.csv")
POLES_CBG_PATH    = os.path.join(RUN_OUT, "poles_cbg.geojson")
VALIDATION_PATH   = os.path.join(RUN_OUT, "validation_summary.csv")

# same MCMC settings as the archived 04, one citywide run instead of 18
MCMC_WARMUP  = 1000
MCMC_SAMPLES = 1000
MCMC_CHAINS  = 1
MCMC_SEED    = 0

FIG_DPI = 300

# ── constants absorbed from archive/04_reporting_decomposition.py ────────────
STUDY_YEARS = [2021, 2022, 2023, 2024, 2025]

NT_CURATED = list(NT_STREETLIGHT)   # NT = Street Light Outage only

# Spec B companion count (kept verbatim from 04; the NT_all column is written
# but not used as the headline NT here).
NT_ALL_EXCLUDE = [
    "Illegal Dumping",
    "Information Request",
    "Rubbish/Recyclable Material Collection",
    "Sanitation Violation",
    "Dumpster Violation",
    "Sanitation / Dumpster Violation",
]

PLANNING_DISTRICTS_PATH  = os.path.join(DATA, "Planning_Districts.geojson")
PLANNING_DISTRICT_ID_COL = "dist_name"

EXPOSURE_COL = "alloc_avg_s_cnt"    # SafeGraph avg monthly stop counts (cov_cbg)

TIGER_BG_PATH = os.path.join(DATA, "tiger", "tl_2023_42_bg", "tl_2023_42_bg.shp")
COV_CBG_PATH  = os.path.join(OUT, "cov_cbg.geojson")


# =============================================================================
# SECTION 0 - CBG counts / exposure / district assignment (absorbed from 04)
# =============================================================================

def build_cbg_counts(study_years=None) -> gpd.GeoDataFrame:
    """
    Count T (Illegal Dumping), NT_curated and NT_all 311 reports per 2023 CBG,
    pooled over study_years (default STUDY_YEARS).  Returns the Philadelphia
    CBG GeoDataFrame (EPSG:26918) with count columns T, NT_curated, NT_all.
    """
    if study_years is None:
        study_years = STUDY_YEARS
    print(f"\n[counts] Counting 311 reports per CBG (years {study_years}) ...")

    pa_2023      = gpd.read_file(TIGER_BG_PATH)
    philly_26918 = pa_2023[pa_2023["COUNTYFP"] == "101"].to_crs(26918)
    philly_26918["GEOID"] = philly_26918["GEOID"].astype(str)
    print(f"  {len(philly_26918)} Philadelphia CBGs loaded.")

    parts = []
    for yr in study_years:
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


def attach_exposure_and_district(cbg: gpd.GeoDataFrame):
    """
    Merge the exposure column from cov_cbg.geojson (population fallback when
    the foot-traffic column is entirely absent) and assign each CBG to the
    planning district containing its representative point.
    """
    print("\n[exposure] Attaching exposure and planning district ...")

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


# =============================================================================
# SECTION 1 - Street-pole inventory: EDA and poles-per-CBG join
# =============================================================================

def load_poles() -> gpd.GeoDataFrame:
    print("\n[poles] Loading street-pole inventory ...")
    poles = gpd.read_file(POLES_PATH).to_crs(26918)
    poles["owner"] = poles["owner"].replace("", np.nan)
    poles["light_dt"] = pd.to_datetime(poles["light_date"], errors="coerce", utc=True)
    print(f"  {len(poles):,} poles loaded.")
    return poles


def pole_eda(poles: gpd.GeoDataFrame):
    """Descriptive tables + figures. Denominator decision: ALL poles count,
    regardless of owner (residents report outages without knowing the owner)."""
    print("\n[poles] EDA tables ...")

    owner = poles["owner"].fillna("(blank)").value_counts()
    owner.rename("n_poles").to_csv(os.path.join(RUN_OUT, "pole_eda_owner.csv"))
    print(f"  Owners: Streets {owner.get('Streets', 0):,} "
          f"({owner.get('Streets', 0) / len(poles):.1%}), "
          f"PECO {owner.get('PECO', 0):,} ({owner.get('PECO', 0) / len(poles):.1%}), "
          f"{len(owner) - 2} other categories")

    pd.crosstab(poles["owner"].fillna("(blank)"), poles["psip_status"]).to_csv(
        os.path.join(RUN_OUT, "pole_eda_psip_by_owner.csv"))
    pd.crosstab(poles["owner"].fillna("(blank)"), poles["bulb_type"]).to_csv(
        os.path.join(RUN_OUT, "pole_eda_bulb_by_owner.csv"))

    peco = poles[poles["owner"] == "PECO"]
    print(f"  PECO pattern: tap_id=0 for {(peco['tap_id'] == 0).mean():.1%}, "
          f"psip_status='N/A' for {(peco['psip_status'] == 'N/A').mean():.1%}, "
          f"bulb_type unknown for {(peco['bulb_type'] == 'UNKNOWN').mean():.1%} "
          f"— PECO poles carry no PSIP/attribute data but stay in the denominator.")

    attr_cols = ["type", "nlumin", "lum_size", "height", "pole_date", "up_date",
                 "light_date", "psip_status", "bulb_type", "block", "plate"]
    miss = (poles[attr_cols].isna().mean() * 100).round(1).rename("pct_missing")
    miss.to_csv(os.path.join(RUN_OUT, "pole_eda_missingness.csv"))

    desc = pd.concat([
        poles["type"].value_counts().head(20).rename("count").to_frame()
             .assign(variable="type"),
        poles["nlumin"].value_counts(dropna=False).head(10).rename("count").to_frame()
             .assign(variable="nlumin"),
        poles["lum_size"].value_counts(dropna=False).head(10).rename("count").to_frame()
             .assign(variable="lum_size"),
    ])
    desc.to_csv(os.path.join(RUN_OUT, "pole_eda_type_nlumin_lumsize.csv"))
    n_zero_lum = int((poles["nlumin"] == 0).sum())
    print(f"  nlumin: {n_zero_lum:,} poles with 0 luminaires, "
          f"{poles['nlumin'].isna().sum():,} missing (incl. all PECO) — the "
          f"denominator is POLE count, not luminaire count.")
    print(f"  height: median {poles['height'].median():.0f} ft "
          f"(missing {poles['height'].isna().mean():.1%})")

    # light_date (PSIP relighting) year histogram — the assumption-(ii) threat
    yr = poles["light_dt"].dt.year.value_counts().sort_index()
    yr.rename("n_relit").to_csv(os.path.join(RUN_OUT, "pole_eda_light_date_years.csv"))
    n_dated = int(poles["light_dt"].notna().sum())
    n_pre24 = int((poles["light_dt"] < PSIP_CUTOFF).sum())
    print(f"  light_date: {n_dated:,} dated relights ({n_dated / len(poles):.0%}); "
          f"{n_pre24:,} before 2024 vs {n_dated - n_pre24:,} in 2024+ — the LED "
          f"conversion landed mostly inside the study window's tail.")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(yr.index.astype(int), yr.values, color="tab:blue")
    ax.axvspan(2020.5, 2025.5, color="orange", alpha=0.15, label="study window 2021–2025")
    ax.set_xlabel("light_date year (PSIP relighting)")
    ax.set_ylabel("poles relit")
    ax.set_title("PSIP LED relighting dates vs the study window")
    ax.legend()
    out = os.path.join(FIG_OUT, "pole_light_date_hist.png")
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")


def join_poles_to_cbg(poles: gpd.GeoDataFrame, cbg: gpd.GeoDataFrame,
                      districts: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Spatial join poles → CBGs; returns cbg with n_poles / led_share /
    pre24_share / peco_share columns. Also writes poles_cbg.geojson + maps."""
    print("\n[poles] Joining poles to CBGs ...")
    j = poles.sjoin(cbg[["GEOID", "geometry"]], how="inner", predicate="within")
    print(f"  {len(j):,}/{len(poles):,} poles fall inside a CBG "
          f"({len(poles) - len(j):,} outside the CBG layer).")

    j["is_led"]   = j["bulb_type"].eq("LED")
    j["is_peco"]  = j["owner"].eq("PECO")
    j["is_pre24"] = j["light_dt"] < PSIP_CUTOFF
    agg = j.groupby("GEOID").agg(
        n_poles=("pole_num", "size"),
        n_led=("is_led", "sum"),
        n_peco=("is_peco", "sum"),
        n_pre24=("is_pre24", "sum"),
    ).reset_index()

    cbg = cbg.merge(agg, on="GEOID", how="left")
    for c in ["n_poles", "n_led", "n_peco", "n_pre24"]:
        cbg[c] = cbg[c].fillna(0).astype(int)
    with np.errstate(invalid="ignore"):
        cbg["led_share"]   = np.where(cbg["n_poles"] > 0, cbg["n_led"]   / cbg["n_poles"], np.nan)
        cbg["peco_share"]  = np.where(cbg["n_poles"] > 0, cbg["n_peco"]  / cbg["n_poles"], np.nan)
        cbg["pre24_share"] = np.where(cbg["n_poles"] > 0, cbg["n_pre24"] / cbg["n_poles"], np.nan)

    q = cbg["n_poles"]
    print(f"  Poles per CBG: median {q.median():.0f} "
          f"(p5 {q.quantile(.05):.0f}, p95 {q.quantile(.95):.0f}); "
          f"{(q == 0).sum()} zero-pole CBGs (excluded from the fit), "
          f"{(q < 10).sum()} CBGs with <10 poles (kept — the prior shrinks them).")
    if (q == 0).any():
        print(f"  Zero-pole GEOIDs: {cbg.loc[q == 0, 'GEOID'].tolist()}")
    print(f"  PECO share per CBG: median {cbg['peco_share'].median():.2f}, "
          f"p90 {cbg['peco_share'].quantile(.9):.2f}, "
          f"{int((cbg['peco_share'] > 0.5).sum())} CBGs >50% PECO — PECO cannot "
          f"drive cross-CBG variation; keep all poles per design decision.")
    print(f"  LED share per CBG: p10 {cbg['led_share'].quantile(.1):.2f}, "
          f"median {cbg['led_share'].median():.2f}, "
          f"p90 {cbg['led_share'].quantile(.9):.2f}")

    keep = ["GEOID", "n_poles", "n_led", "led_share", "n_peco", "peco_share",
            "n_pre24", "pre24_share", "geometry"]
    gpd.GeoDataFrame(cbg[keep], geometry="geometry", crs=cbg.crs).to_crs(4326).to_file(
        POLES_CBG_PATH, driver="GeoJSON")
    print(f"  Saved → {POLES_CBG_PATH}")

    plot_choropleth(
        cbg, "n_poles", title="Street poles per CBG (all owners)",
        out_path=os.path.join(FIG_OUT, "pole_map_per_cbg.png"),
        dpi=FIG_DPI, boundary_overlay=districts,
    )
    plot_choropleth(
        cbg, "led_share", title="LED share of poles per CBG (bulb_type = LED)",
        out_path=os.path.join(FIG_OUT, "pole_map_led_share.png"),
        dpi=FIG_DPI, cmap="magma", boundary_overlay=districts,
    )
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(q, bins=60, color="tab:blue")
    ax.set_xlabel("poles per CBG"); ax.set_ylabel("CBGs")
    ax.set_title("Distribution of street poles per CBG")
    out = os.path.join(FIG_OUT, "pole_hist_per_cbg.png")
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")
    return cbg


# =============================================================================
# SECTION 2 - Citywide two-equation model (poles as NT exposure)
# =============================================================================

def two_eq_model_city(log_P, log_E, nt_obs=None, t_obs=None):
    """
    04's two_eq_model with (a) a single citywide fit — no district grouping,
    no corrected_level / a_ref anchor — and (b) separate exposures:

      NT_i ~ Poisson( P_i * exp( a_NT + w_R,i ) )         P_i = poles in CBG i
      T_i  ~ Poisson( E_i * exp( a_T + beta*w_R,i + w_T,i ) )   E_i = foot traffic

    Identifying assumption: true breakage per pole is uniform citywide, so
    cross-CBG variation in outages-per-pole is reporting propensity w_R.
    """
    n = log_P.shape[0]

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

    with numpyro.plate("obs", n):
        numpyro.sample("NT", dist.Poisson(jnp.exp(log_P + a_nt + w_r)),
                       obs=nt_obs)
        numpyro.sample("T",  dist.Poisson(jnp.exp(log_E + a_t + beta * w_r + w_t)),
                       obs=t_obs)


def fit_citywide(sub: pd.DataFrame, nt_col: str, seed: int) -> dict:
    """One NUTS run over all fitted CBGs. sub needs GEOID, T_count, {nt_col},
    exposure > 0, n_poles > 0."""
    log_P = jnp.asarray(np.log(sub["n_poles"].to_numpy(float)))
    log_E = jnp.asarray(np.log(sub["exposure"].to_numpy(float)))
    nt    = jnp.asarray(sub[nt_col].to_numpy(int))
    t     = jnp.asarray(sub["T_count"].to_numpy(int))

    kernel = NUTS(two_eq_model_city, target_accept_prob=0.9)
    mcmc = MCMC(kernel, num_warmup=MCMC_WARMUP, num_samples=MCMC_SAMPLES,
                num_chains=MCMC_CHAINS, progress_bar=False)
    mcmc.run(jax.random.PRNGKey(seed), log_P=log_P, log_E=log_E,
             nt_obs=nt, t_obs=t, extra_fields=("diverging",))
    s = mcmc.get_samples()
    n_div = int(np.asarray(mcmc.get_extra_fields()["diverging"]).sum())

    w_r = np.asarray(s["w_R"])
    w_t = np.asarray(s["w_T"])
    a_nt, a_t = np.asarray(s["a_NT"]), np.asarray(s["a_T"])
    beta = np.asarray(s["beta"])

    log_P_np, log_E_np = np.asarray(log_P), np.asarray(log_E)
    nt_fit = np.exp(log_P_np[None, :] + a_nt[:, None] + w_r).mean(axis=0)
    t_fit  = np.exp(log_E_np[None, :] + a_t[:, None]
                    + beta[:, None] * w_r + w_t).mean(axis=0)

    cbg = pd.DataFrame({
        "GEOID":     sub["GEOID"].to_numpy(),
        "w_R_mean":  w_r.mean(axis=0), "w_R_sd": w_r.std(axis=0),
        "w_T_mean":  w_t.mean(axis=0), "w_T_sd": w_t.std(axis=0),
        "NT_fitted": nt_fit, "T_fitted": t_fit,
    })
    scalars = {
        "a_NT_mean": a_nt.mean(), "a_T_mean": a_t.mean(),
        "beta_mean": beta.mean(),
        "beta_lo": np.quantile(beta, 0.03), "beta_hi": np.quantile(beta, 0.97),
        "sigma_R_mean": np.asarray(s["sigma_R"]).mean(),
        "sigma_T_mean": np.asarray(s["sigma_T"]).mean(),
        "n_cbg": len(sub), "n_divergences": n_div,
    }
    return {"cbg": cbg, "scalars": scalars}


# =============================================================================
# SECTION 3 - Validation and PSIP robustness
# =============================================================================

def _corr_row(name, x, y, note=""):
    d = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    return {
        "check": name, "n": len(d),
        "pearson":  d["x"].corr(d["y"]) if len(d) > 2 else np.nan,
        "spearman": d["x"].corr(d["y"], method="spearman") if len(d) > 2 else np.nan,
        "note": note,
    }


def _scatter(x, y, xlabel, ylabel, title, fname):
    d = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 3:
        return
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(d["x"], d["y"], s=8, alpha=0.4)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    r = d["x"].corr(d["y"], method="spearman")
    ax.set_title(f"{title} (spearman r = {r:.3f}, n = {len(d):,})")
    out = os.path.join(FIG_OUT, fname)
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")


def run_validation(cbg: gpd.GeoDataFrame, archive_cbg_path: str):
    """Proxy validation + PSIP LED-share sensitivity. w_R columns already
    attached to cbg (w_R_mean_sl = full window, w_R_mean_pre24 = 2021-23 NT)."""
    print("\n[validate] Proxy validation and PSIP sensitivity ...")
    rows = []
    w_full = cbg["w_R_mean_sl"]

    # (1) vs the archived district-wise fit (the refactor's effect)
    if os.path.exists(archive_cbg_path):
        old = gpd.read_file(archive_cbg_path)[["GEOID", "w_R_mean_sl"]]
        old["GEOID"] = old["GEOID"].astype(str)
        m = cbg[["GEOID", "w_R_mean_sl"]].merge(
            old, on="GEOID", suffixes=("_city", "_district"))
        rows.append(_corr_row(
            "w_R citywide vs district-wise fit",
            m["w_R_mean_sl_city"], m["w_R_mean_sl_district"],
            "old fit was within-district centred; citywide adds level variation"))
        _scatter(m["w_R_mean_sl_district"], m["w_R_mean_sl_city"],
                 "w_R — district-wise fit (archived)", "w_R — citywide fit",
                 "Refactor effect on w_R", "val_scatter_wR_city_vs_district.png")
    else:
        print("  [skip] no archived district-wise decomposition to compare against.")

    # (2) vs raw outages per pole (mechanical sanity check)
    with np.errstate(divide="ignore", invalid="ignore"):
        opp = np.log(cbg["NT_curated"].replace(0, np.nan) / cbg["n_poles"].replace(0, np.nan))
    rows.append(_corr_row("w_R vs log outages-per-pole", w_full, opp,
                          "sanity: w_R is shrunk log outage-per-pole"))
    _scatter(opp, w_full, "log(outage reports / pole), 2021–25",
             "w_R — citywide fit", "w_R vs raw outages per pole",
             "val_scatter_wR_vs_outage_per_pole.png")

    # (3) PSIP pre-conversion time slice
    if "w_R_mean_pre24" in cbg.columns:
        rows.append(_corr_row("w_R full window vs 2021–23 NT window",
                              w_full, cbg["w_R_mean_pre24"],
                              "high corr → assumption (ii) robust to PSIP timing"))
        _scatter(cbg["w_R_mean_pre24"], w_full,
                 "w_R — NT from 2021–2023 only", "w_R — NT from full window",
                 "PSIP robustness: pre-conversion time slice",
                 "val_scatter_wR_full_vs_pre24.png")

    # (4) PSIP LED-share sensitivity (no refit)
    rows.append(_corr_row("w_R vs LED share", w_full, cbg["led_share"],
                          "material corr → conversion geography contaminates proxy"))
    rows.append(_corr_row("w_R vs pre-2024 relit share", w_full, cbg["pre24_share"],
                          "material corr → conversion geography contaminates proxy"))
    _scatter(cbg["led_share"], w_full, "LED share of poles",
             "w_R — citywide fit", "PSIP sensitivity: LED share",
             "val_scatter_wR_vs_led_share.png")

    # (5) descriptive: the two denominators
    with np.errstate(divide="ignore"):
        rows.append(_corr_row("log poles vs log foot-traffic exposure",
                              np.log(cbg["n_poles"].replace(0, np.nan)),
                              np.log(cbg["exposure"].where(cbg["exposure"] > 0)),
                              "how different the NT and T denominators are"))

    val = pd.DataFrame(rows)
    val.to_csv(VALIDATION_PATH, index=False)
    print(f"  Saved → {VALIDATION_PATH}")
    print(val[["check", "n", "pearson", "spearman"]].to_string(index=False))

    led_r = val.loc[val["check"] == "w_R vs LED share", "spearman"].iloc[0]
    if abs(led_r) > 0.3:
        print(f"  [decision rule] |corr(w_R, LED share)| = {abs(led_r):.2f} > 0.3 — "
              f"PSIP geography contaminates the proxy; prefer the 2021–23 NT window "
              f"as the headline spec.")
    else:
        print(f"  [decision rule] |corr(w_R, LED share)| = {abs(led_r):.2f} ≤ 0.3 — "
              f"full-window NT stands as the headline spec.")


# =============================================================================
# SECTION 4 - STEP 1: citywide decomposition
# =============================================================================

def run_decomposition():
    print("\n" + "=" * 60)
    print("05a — citywide W_R decomposition, NT = Street Light Outage,")
    print("      poles as NT exposure (single fit, no district layer)")
    print("=" * 60)
    # pre-refactor (district-wise) outputs, kept for the validation comparison
    archive_cbg = DECOMP_CBG_PATH.replace(".geojson", "_districtfit.geojson")

    # counts: full window, then the pre-conversion NT slice (2021-2023)
    counts = build_cbg_counts()
    counts_pre = build_cbg_counts(PRE24_YEARS)[["GEOID", "NT_curated"]].rename(
        columns={"NT_curated": "NT_pre24"})
    counts = counts.merge(counts_pre, on="GEOID", how="left")
    counts["NT_pre24"] = counts["NT_pre24"].fillna(0).astype(int)

    cbg, districts = attach_exposure_and_district(counts)

    # pole inventory
    poles = load_poles()
    pole_eda(poles)
    cbg = join_poles_to_cbg(poles, cbg, districts)

    fit_mask = (cbg["exposure"] > 0) & (cbg["n_poles"] > 0)
    sub = cbg[fit_mask]
    print(f"\n[fit] Citywide fit: {len(sub)}/{len(cbg)} CBGs "
          f"(excluded: {int((~(cbg['exposure'] > 0)).sum())} without exposure, "
          f"{int((cbg['n_poles'] == 0).sum())} without poles); "
          f"T={sub['T_count'].sum():,}, NT={sub['NT_curated'].sum():,}")

    import time
    t0 = time.time()
    res = fit_citywide(sub, "NT_curated", seed=MCMC_SEED)
    sc = res["scalars"]
    print(f"  done ({time.time() - t0:.0f}s, {sc['n_divergences']} divergences)  "
          f"beta={sc['beta_mean']:.3f} [{sc['beta_lo']:.3f}, {sc['beta_hi']:.3f}]  "
          f"a_NT={sc['a_NT_mean']:.3f}  a_T={sc['a_T_mean']:.3f}  "
          f"sigma_R={sc['sigma_R_mean']:.3f}  sigma_T={sc['sigma_T_mean']:.3f}")

    # PSIP robustness refit: NT restricted to 2021-2023 (T stays full window)
    print("\n[fit] PSIP robustness refit: NT window 2021–2023 ...")
    t0 = time.time()
    res_pre = fit_citywide(sub, "NT_pre24", seed=MCMC_SEED + 1)
    sp = res_pre["scalars"]
    print(f"  done ({time.time() - t0:.0f}s, {sp['n_divergences']} divergences)  "
          f"beta={sp['beta_mean']:.3f}")

    # attach per-CBG results (suffix _sl keeps the established naming)
    for col in ["w_R_mean", "w_R_sd", "w_T_mean", "w_T_sd", "NT_fitted", "T_fitted"]:
        cbg[f"{col}_sl"] = cbg["GEOID"].map(res["cbg"].set_index("GEOID")[col])
    cbg["w_R_mean_pre24"] = cbg["GEOID"].map(res_pre["cbg"].set_index("GEOID")["w_R_mean"])

    cbg.to_crs(4326).to_file(DECOMP_CBG_PATH, driver="GeoJSON")
    print(f"  Saved → {DECOMP_CBG_PATH}")

    city = pd.DataFrame([
        dict(spec="sl", nt_window="2021-2025", **sc),
        dict(spec="sl_pre24", nt_window="2021-2023", **sp),
    ])
    city.to_csv(CITY_SUMMARY_PATH, index=False)
    print(f"  Saved → {CITY_SUMMARY_PATH}")

    # figures: citywide maps (directly comparable across all CBGs) + PPC
    plot_choropleth(
        cbg, "w_R_mean_sl",
        title="Reporting propensity w_R (citywide fit, poles as NT exposure)",
        out_path=os.path.join(FIG_OUT, "sl_map_w_R.png"),
        dpi=FIG_DPI, boundary_overlay=districts,
    )
    plot_choropleth(
        cbg, "w_T_mean_sl",
        title="Reporting-corrected dumping w_T (citywide fit)",
        out_path=os.path.join(FIG_OUT, "sl_map_w_T.png"),
        dpi=FIG_DPI, cmap="inferno", boundary_overlay=districts,
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, obs_col, fit_col, label in [
        (axes[0], "T_count", "T_fitted_sl", "Trash (T), exposure = foot traffic"),
        (axes[1], "NT_curated", "NT_fitted_sl", "Street-light outages (NT), exposure = poles"),
    ]:
        d = cbg[[obs_col, fit_col]].dropna()
        ax.scatter(d[obs_col] + 0.5, d[fit_col] + 0.5, s=8, alpha=0.4)
        ax.set_xscale("log"); ax.set_yscale("log")
        lims = [0.4, max((d[obs_col].max() if len(d) else 1), 1) * 2]
        ax.plot(lims, lims, ls="--", color="gray")
        ax.set_xlabel(f"Observed (+0.5, log)"); ax.set_ylabel("Posterior-mean fitted (+0.5, log)")
        ax.set_title(label)
    fig.suptitle("Posterior-predictive fit — citywide two-exposure model")
    out = os.path.join(FIG_OUT, "sl_ppc_fit.png")
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")

    run_validation(cbg, archive_cbg)


if __name__ == "__main__":
    run_decomposition()
