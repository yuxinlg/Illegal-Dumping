"""
05b_corrected_hawkes_streetlight.py
====================================
Street-light-only, CITYWIDE variant of the reporting-corrected Cox-Hawkes
pipeline. Documented in 05b_streetlight_citywide.md.

Compared to the earlier version of this script (which re-ran 04's per-district
decomposition with NT = Street Light Outage), this version:

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
  4. Feeds the simpler citywide offset  beta * w_R_i  into the per-district
     Cox-Hawkes fits (baseline / offset_fixed / offset_free) for Central
     (Center City) and University Southwest (University City).

No existing file is modified: 04 and 05 are loaded as modules; 04's count /
exposure builders are reused as-is; 05's district-based build_offset is
replaced in memory with the citywide version.

Usage
-----
  python 05b_corrected_hawkes_streetlight.py [decomp|hawkes|all]   (default all)

Outputs
-------
  output/corrected_hawkes_streetlight/
    reporting_decomp_cbg_sl.geojson            citywide decomposition per CBG
    reporting_decomp_city_summary_sl.csv       one row: a_NT, a_T, beta, sigmas
    reporting_decomp_*_districtfit.*           archived pre-refactor outputs
    poles_cbg.geojson                          poles per CBG + LED/PECO shares
    pole_eda_*.csv                             EDA tables
    validation_summary.csv                     all validation correlations
    offset_cbg.geojson, fit_summary.csv/.geojson, cbg_fxy.geojson  (05 outputs)
  output/figures/corrected_hawkes_streetlight/*.png
"""

import os
import sys
import shutil
import importlib.util
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
FIT_DISTRICTS  = ["Central", "University Southwest"]   # Center City, University City

POLES_PATH   = os.path.join(DATA, "Street_Poles.geojson")
PSIP_CUTOFF  = pd.Timestamp("2024-01-01", tz="UTC")    # mass LED relighting starts
PRE24_YEARS  = [2021, 2022, 2023]                      # pre-conversion NT window

DECOMP_CBG_PATH   = os.path.join(RUN_OUT, "reporting_decomp_cbg_sl.geojson")
CITY_SUMMARY_PATH = os.path.join(RUN_OUT, "reporting_decomp_city_summary_sl.csv")
POLES_CBG_PATH    = os.path.join(RUN_OUT, "poles_cbg.geojson")
VALIDATION_PATH   = os.path.join(RUN_OUT, "validation_summary.csv")

# same MCMC settings as 04, one citywide run instead of 18 district runs
MCMC_WARMUP  = 1000
MCMC_SAMPLES = 1000
MCMC_CHAINS  = 1
MCMC_SEED    = 0

FIG_DPI = 300


def _load(fname, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(REPL_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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
                      districts: gpd.GeoDataFrame, mod04) -> gpd.GeoDataFrame:
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

    mod04.plot_choropleth(
        cbg, "n_poles", title="Street poles per CBG (all owners)",
        out_path=os.path.join(FIG_OUT, "pole_map_per_cbg.png"),
        dpi=FIG_DPI, boundary_overlay=districts,
    )
    mod04.plot_choropleth(
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
    print("STEP 1 — citywide W_R decomposition, NT = Street Light Outage,")
    print("         poles as NT exposure (single fit, no district layer)")
    print("=" * 60)
    mod04 = _load("04_reporting_decomposition.py", "reporting_decomposition")
    mod04.NT_CURATED = NT_STREETLIGHT      # 04 reads this at call time

    # archive the pre-refactor (district-wise) outputs once, for comparison
    archive_cbg = DECOMP_CBG_PATH.replace(".geojson", "_districtfit.geojson")
    old_district_summary = os.path.join(RUN_OUT, "reporting_decomp_district_summary_sl.csv")
    for src, dst in [
        (DECOMP_CBG_PATH, archive_cbg),
        (old_district_summary, old_district_summary.replace(".csv", "_districtfit.csv")),
    ]:
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)
            print(f"  Archived pre-refactor output → {dst}")

    # counts: full window, then the pre-conversion NT slice (2021-2023)
    counts = mod04.build_cbg_counts()
    full_years = list(mod04.STUDY_YEARS)
    mod04.STUDY_YEARS = PRE24_YEARS
    counts_pre = mod04.build_cbg_counts()[["GEOID", "NT_curated"]].rename(
        columns={"NT_curated": "NT_pre24"})
    mod04.STUDY_YEARS = full_years
    counts = counts.merge(counts_pre, on="GEOID", how="left")
    counts["NT_pre24"] = counts["NT_pre24"].fillna(0).astype(int)

    cbg, districts = mod04.attach_exposure_and_district(counts)

    # STEP 0: pole inventory
    poles = load_poles()
    pole_eda(poles)
    cbg = join_poles_to_cbg(poles, cbg, districts, mod04)

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
    mod04.plot_choropleth(
        cbg, "w_R_mean_sl",
        title="Reporting propensity w_R (citywide fit, poles as NT exposure)",
        out_path=os.path.join(FIG_OUT, "sl_map_w_R.png"),
        dpi=FIG_DPI, boundary_overlay=districts,
    )
    mod04.plot_choropleth(
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


# =============================================================================
# SECTION 5 - STEP 2: Cox-Hawkes with the citywide offset
# =============================================================================

def build_offset_citywide(cov_gdf, ch05):
    """Replaces 05's district-based build_offset: offset_i = beta * w_R_i.
    With a single citywide fit there is no district layer — the level is
    absorbed identically everywhere by a_0."""
    print("\n[offset] Building citywide reporting offset  beta * w_R ...")
    decomp = gpd.read_file(DECOMP_CBG_PATH)[["GEOID", "dist_name", "w_R_mean_sl"]]
    decomp["GEOID"] = decomp["GEOID"].astype(str)
    city = pd.read_csv(CITY_SUMMARY_PATH)
    beta = float(city.loc[city["spec"] == "sl", "beta_mean"].iloc[0])
    print(f"  beta (citywide) = {beta:.4f}")

    decomp[ch05.OFFSET_COL] = beta * decomp["w_R_mean_sl"]
    cov_gdf = cov_gdf.merge(decomp[["GEOID", "dist_name", ch05.OFFSET_COL]],
                            on="GEOID", how="left")
    n_missing = int(cov_gdf[ch05.OFFSET_COL].isna().sum())
    cov_gdf[ch05.OFFSET_COL] = cov_gdf[ch05.OFFSET_COL].fillna(0.0)
    print(f"  Offset merged: {len(cov_gdf) - n_missing}/{len(cov_gdf)} CBGs; "
          f"{n_missing} filled with 0 (the citywide mean)")
    print(f"  Offset range: [{cov_gdf[ch05.OFFSET_COL].min():.3f}, "
          f"{cov_gdf[ch05.OFFSET_COL].max():.3f}]")

    # dsum kept API-compatible with 05's main(); corrected_level is retired —
    # NaN makes the downstream a_0-vs-corrected_level figure skip itself.
    dsum = pd.DataFrame({
        "dist_name": decomp["dist_name"].dropna().unique(),
        "corrected_level_mean": np.nan,
    })
    return cov_gdf, dsum


def run_cox_hawkes():
    print("\n" + "=" * 60)
    print(f"STEP 2 — Cox-Hawkes with citywide street-light offset: {FIT_DISTRICTS}")
    print("=" * 60)
    ch05 = _load("05_corrected_hawkes_by_district.py", "corrected_hawkes")

    ch05.RUN_NAME = RUN_NAME
    ch05.RUN_OUT  = RUN_OUT
    ch05.FIG_OUT  = FIG_OUT
    ch05.UNIT_SUBSET = FIT_DISTRICTS
    ch05.build_offset = lambda cov_gdf: build_offset_citywide(cov_gdf, ch05)

    ch05.main()


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("decomp", "all"):
        run_decomposition()
    if step in ("hawkes", "all"):
        run_cox_hawkes()
