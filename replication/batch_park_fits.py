"""
Batch Cox–Hawkes fits for all park boxes; save figures and metrics for comparison.

Used from cbg-park-seasonal-refactor-preview.ipynb after data-prep cells
(``illegal_dumping_in_park``, etc.). This preview preserves the legacy
citywide-row covariate standardization explicitly and uses real spatial units.
"""

from __future__ import annotations

import contextlib
import gc
import io
import json
import time
from datetime import datetime
from pathlib import Path

import bstpp
import bstpp.main
import jax
import matplotlib.pyplot as plt
import numpy as np
import numpyro.distributions as dist
import pandas as pd
from tqdm.auto import tqdm

from bstpp.cutoffs import days_to_internal, internal_to_days

COLUMN_NAMES = [
    "pop_density",
    "med_inc_avg",
    "lep_density",
    "edu_hs_avg",
    "ndvi_mean_4yr",
    "RLD",
    "RMD",
    "RHD",
    "use_com",
    "use_civ",
    "I",
    "T",
    "GW",
    "vac_area",
    "landcare_area",
    "alloc_avg_d_cnt",
    "alloc_avg_s_cnt",
    "unique_device_ratio_aw",
    "reporting_rate",
    "betweenness_avg_w",
    "pagerank_avg_w",
    "betweenness_avg_d",
    "pagerank_avg_d",
]

# PRE-S3 MIGRATION EDIT: spatial trigger quantities are real metric units.
# ``spatial_window`` is a per-axis square cutoff, so a value of 1,000 keeps
# pairs with |dx| <= 1,000 m and |dy| <= 1,000 m. The variance prior remains a
# provisional sensitivity choice; 50/100/250 m are the analysis candidates.
SPATIAL_WINDOW_M = 1_000.0
SPATIAL_SIGMA_SENSITIVITY_M = (50.0, 100.0, 250.0)
DEFAULT_SPATIAL_SIGMA_PRIOR_SCALE_M = 250.0

# PRE-S3 MIGRATION EDIT: four calendar years include the 2024 leap day.
ANALYSIS_START = pd.Timestamp("2021-01-01")
ANALYSIS_END = pd.Timestamp("2025-01-01")
ANALYSIS_TOTAL_DAYS = int((ANALYSIS_END - ANALYSIS_START).days)

PARK_NAMES = [
    "tacony_box",
    "cobbscreek_box",
    "mifflin_box",
    "fairmount_box",
]

PARK_LABELS = {
    "tacony_box": "Tacony Creek",
    "cobbscreek_box": "Cobbs Creek",
    "mifflin_box": "Mifflin Square",
    "fairmount_box": "Fairmount Park",
}

# Default Hawkes cutoffs (see ``resolve_temporal_window`` / ``get_park_fit_kwargs``).
# Prefer ``temporal_cutoff_days`` (real calendar days) when changing the temporal
# computational cutoff — ``window`` is the legacy *internal*-unit form
# (T_INTERNAL=50 maps the full analysis horizon). Passing both is an error.
# With the 1,461-day preview horizon, ``window=100`` ≈ 2,922 real days and does
# not thin pairs at all; use a finite ``temporal_cutoff_days`` to reduce pair count.
DEFAULT_FIT_KWARGS = {
    "window": 50.0,
    "spatial_window": SPATIAL_WINDOW_M,
    "spatial_sigma_prior_scale_m": DEFAULT_SPATIAL_SIGMA_PRIOR_SCALE_M,
    "num_samples": 750,
}

# Per-park overrides. Optional ``num_steps``: int or list[int] (else uses ``num_steps_list``).
# Temporal cutoffs may be set with ``temporal_cutoff_days`` (preferred) or ``window``.
PARK_FIT_OVERRIDES: dict[str, dict] = {
    "fairmount_box": {
        "num_samples": 500,
    },
}


def _resolve_log_level(verbose: bool | None, log_level: str | None) -> str:
    """verbose=False → quiet; verbose=True → normal; else use log_level (default normal)."""
    if verbose is False:
        return "quiet"
    if log_level is not None:
        return log_level
    return "normal"


def _ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


def _log(msg: str, *, log_level: str = "normal") -> None:
    """Print a short status line (use log_level='verbose' for detail, 'quiet' for errors only)."""
    if log_level == "quiet":
        return
    if log_level == "normal" and msg.startswith("  "):
        return
    tqdm.write(f"[{_ts()}] {msg}")


def release_jax_memory(*, log_level: str = "normal") -> None:
    """Drop unreachable Python objects and JAX compilation / staging caches.

    Safe between independent park fits: the next ``run_svi`` recompiles. Does
    **not** free arrays still referenced by live model/sample objects — delete
    those first. Does not mutate ``jax_enable_x64`` or other process-global
    precision flags.
    """
    gc.collect()
    jax.clear_caches()
    _log("Cleared JAX compilation caches", log_level=log_level)


def resolve_temporal_window(
    fit_kwargs: dict,
    total_days: int | float,
) -> tuple[float, float, str]:
    """Resolve the temporal computational cutoff for one fit.

    Accepts either ``temporal_cutoff_days`` (real days; preferred) or legacy
    ``window`` (internal units on ``[0, T_INTERNAL]``). Returns
    ``(window_internal, window_days, source)`` where ``source`` is
    ``\"temporal_cutoff_days\"`` or ``\"window\"``.
    """
    has_days = fit_kwargs.get("temporal_cutoff_days") is not None
    has_window = fit_kwargs.get("window") is not None
    if has_days and has_window:
        raise ValueError(
            "Pass temporal_cutoff_days or window, not both "
            "(real days vs legacy internal units)."
        )
    if not has_days and not has_window:
        raise ValueError(
            "A temporal cutoff is required: set temporal_cutoff_days "
            "(preferred, real days) or window (legacy internal units)."
        )

    horizon = float(total_days)
    if not (np.isfinite(horizon) and horizon > 0):
        raise ValueError(f"total_days must be finite and > 0; got {total_days}")

    if has_days:
        days = float(fit_kwargs["temporal_cutoff_days"])
        if not (np.isfinite(days) and days > 0):
            raise ValueError(
                f"temporal_cutoff_days must be finite and > 0; got {days}"
            )
        return days_to_internal(days, horizon), days, "temporal_cutoff_days"

    window = float(fit_kwargs["window"])
    if not (np.isfinite(window) and window > 0):
        raise ValueError(f"window must be finite and > 0; got {window}")
    return window, internal_to_days(window, horizon), "window"


@contextlib.contextmanager
def _suppress_bstpp_chatter(enabled: bool):
    """Hide set_window DEBUG prints so the SVI progress bar stays visible."""
    if not enabled:
        yield
        return
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        yield


@contextlib.contextmanager
def _suppress_plt_show():
    """Prevent BSTPP helpers from popping figures during save-only batch plots."""
    _show = plt.show
    plt.show = lambda *args, **kwargs: None
    try:
        yield
    finally:
        plt.show = _show


def ensure_b_0_in_samples(model) -> None:
    """
    Spatial plots use b_0 = X(s)w (log-intensity offset per CBG).

    Older run_svi only drew sampled sites (w, not deterministic b_0), which caused
    KeyError('b_0') in plot_spatial(include_cov=True).
    """
    if "spatial_cov" not in model.args:
        return
    if "b_0" in model.samples:
        return
    if "w" not in model.samples:
        raise KeyError(
            "Spatial covariates in model but posterior has neither 'b_0' nor 'w'."
        )
    w = np.asarray(model.samples["w"])
    if w.ndim == 1:
        w = w[:, None]
    spatial_cov = np.asarray(model.args["spatial_cov"])
    model.samples["b_0"] = w @ spatial_cov.T


def get_park_fit_kwargs(
    park_name: str,
    defaults: dict | None = None,
    park_overrides: dict[str, dict] | None = None,
    **overrides,
) -> dict:
    """Merge default Hawkes settings for one park.

    Cutoff keys: ``temporal_cutoff_days`` (real days; preferred) or ``window``
    (legacy internal), plus ``spatial_window``, ``spatial_sigma_prior_scale_m``,
    ``num_samples``. See ``resolve_temporal_window``.
    """
    base = defaults if defaults is not None else DEFAULT_FIT_KWARGS
    per_park = park_overrides if park_overrides is not None else PARK_FIT_OVERRIDES
    merged = {**base, **per_park.get(park_name, {}), **overrides}
    merged.pop("num_steps", None)  # batch job planner only, not passed to set_window
    # Prefer an explicit real-day cutoff over the inherited legacy default.
    if merged.get("temporal_cutoff_days") is not None:
        merged.pop("window", None)
    return merged


def legacy_citywide_standardize_covariates(cov_gdf):
    """Reproduce the old package's unweighted all-row standardization.

    The historical notebook supplied all 1,338 CBG rows to every park fit and
    ``standardize_cov=True`` calculated population moments before geographic
    clipping. The refactored package rejects that boolean, so the preview does
    the identical operation explicitly and passes ``standardize_cov=None``.
    """
    values = cov_gdf[COLUMN_NAMES].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Legacy standardization requires finite covariate values")
    mean = values.mean(axis=0)
    scale = values.var(axis=0) ** 0.5
    if np.any(scale <= 0.0):
        bad = [COLUMN_NAMES[i] for i in np.flatnonzero(scale <= 0.0)]
        raise ValueError(f"Legacy standardization has zero-scale columns: {bad}")

    standardized = cov_gdf.copy()
    standardized.loc[:, COLUMN_NAMES] = (values - mean) / scale
    record = {
        "method": "legacy_unweighted_all_supplied_rows",
        "row_count": int(values.shape[0]),
        "columns": list(COLUMN_NAMES),
        "mean": mean.tolist(),
        "scale": scale.tolist(),
    }
    return standardized, record


def _normalize_num_steps(steps: int | list[int]) -> list[int]:
    if isinstance(steps, int):
        return [steps]
    return [int(s) for s in steps]


def num_steps_for_park(
    park_name: str,
    num_steps_list: list[int],
    park_num_steps: dict[str, int | list[int]] | None = None,
    park_overrides: dict[str, dict] | None = None,
) -> list[int]:
    """
    SVI step counts for one park.

    Priority: ``park_num_steps[park]`` → ``park_fit_overrides[park]["num_steps"]``
    → ``num_steps_list``.
    """
    if park_num_steps and park_name in park_num_steps:
        return _normalize_num_steps(park_num_steps[park_name])
    per_park = park_overrides if park_overrides is not None else PARK_FIT_OVERRIDES
    steps = per_park.get(park_name, {}).get("num_steps")
    if steps is not None:
        return _normalize_num_steps(steps)
    return list(num_steps_list)


def build_batch_jobs(
    park_names: list[str],
    num_steps_list: list[int],
    park_num_steps: dict[str, int | list[int]] | None = None,
    park_overrides: dict[str, dict] | None = None,
) -> list[tuple[str, int]]:
    """(park_name, num_steps) jobs — per-park step overrides expand independently."""
    jobs: list[tuple[str, int]] = []
    for park_name in park_names:
        for num_steps in num_steps_for_park(
            park_name, num_steps_list, park_num_steps, park_overrides
        ):
            jobs.append((park_name, num_steps))
    return jobs


def _tag(
    park_name: str,
    num_steps: int,
    lr: float,
    window: float,
    num_samples: int,
) -> str:
    return (
        f"{park_name}_{num_steps}steps_lr{lr:g}"
        f"_w{window:g}_ns{num_samples}"
    )


def _timestamped_run_dir(base_dir: Path) -> Path:
    """Create a unique run folder under base_dir, e.g. output/batch_fits/20260529_153045/."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = base_dir / stamp
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def _fig_path(
    out_dir: Path,
    park_name: str,
    num_steps: int,
    lr: float,
    window: float,
    num_samples: int,
    plot_name: str,
) -> Path:
    return out_dir / f"{_tag(park_name, num_steps, lr, window, num_samples)}_{plot_name}.png"


def _artifact_path(
    out_dir: Path,
    park_name: str,
    num_steps: int,
    lr: float,
    window: float,
    num_samples: int,
    name: str,
    ext: str,
) -> Path:
    return out_dir / f"{_tag(park_name, num_steps, lr, window, num_samples)}_{name}.{ext}"


KEY_TRIGGER_PARAMS = ("alpha", "beta", "gamma", "sigmax_2")


def _format_key_trigger_means(model) -> list[str]:
    lines = ["Key trigger posterior means:"]
    for name in KEY_TRIGGER_PARAMS:
        if name in model.samples:
            lines.append(f"  {name}={float(model.samples[name].mean()):.6g}")
        else:
            lines.append(f"  {name}=— (not in posterior)")
    return lines


def _trigger_post_summary(model) -> pd.DataFrame:
    par_names = (
        model.args["t_trig"].get_par_names() + model.args["sp_trig"].get_par_names()
    )
    names = ["alpha"] + par_names
    trig_pos = np.stack([model.samples[name] for name in names]).T
    mean = trig_pos.mean(axis=0)
    std = trig_pos.var(axis=0) ** 0.5
    p_val = [(model.samples[name] > 0).mean() for name in names]
    quantiles = np.quantile(trig_pos, [0.025, 0.975], axis=0)
    return pd.DataFrame(
        {
            "Post Mean": mean,
            "Post Std": std,
            "P(w>0)": p_val,
            "[0.025": quantiles[0],
            "0.975]": quantiles[1],
        },
        index=names,
    )


def _cov_weight_post_summary(model) -> pd.DataFrame:
    w_samples = model.samples["w"]
    if w_samples.ndim == 1:
        w_samples = w_samples[:, None]
    w_samps = np.concatenate(
        (w_samples, model.samples["a_0"].reshape(-1, 1)), axis=1
    )
    mean = w_samps.mean(axis=0)
    std = w_samps.var(axis=0) ** 0.5
    p = (w_samps > 0).mean(axis=0)
    quantiles = np.quantile(w_samps, [0.025, 0.975], axis=0)
    return pd.DataFrame(
        {
            "Post Mean": mean,
            "Post Std": std,
            "P(w>0)": p,
            "[0.025": quantiles[0],
            "0.975]": quantiles[1],
        },
        index=model.cov_names + ["a_0"],
    )


def build_fit_summary_text(
    model,
    locs_s: pd.DataFrame,
    *,
    park_name: str,
    num_steps: int,
    lr: float,
    window: float,
    spatial_window: float,
    num_samples: int,
    n_pairs: int | None,
    include_covariate_weights: bool = True,
    key_trigger_means_only: bool = False,
) -> str:
    if key_trigger_means_only:
        return "\n".join(_format_key_trigger_means(model))

    label = PARK_LABELS.get(park_name, park_name)
    loss = np.asarray(model.svi_results.losses)
    prop_excitation = float(
        (model.samples["Itot_excite"] / model.samples["Itot_txy"]).mean()
    )
    loglik = float(model.log_expected_likelihood(locs_s))
    aic = float(model.expected_AIC())

    lines = [
        f"{label} | {num_steps} SVI steps | lr={lr:g}",
        f"window={window:g} | spatial_window={spatial_window:g} | "
        f"posterior draws={num_samples}",
        f"events={len(locs_s):,}"
        + (f" | trigger pairs={n_pairs:,}" if n_pairs is not None else ""),
        f"final SVI loss={float(loss[-1]):.4f}"
        + (" | NaN in loss trace" if np.isnan(loss).any() else ""),
        f"prop. self-excitation (mean)={prop_excitation:.4f}",
        f"log E[likelihood]={loglik:.4f} | E[AIC]={aic:.2f}",
        "",
        "Trigger parameters:",
        _trigger_post_summary(model).to_string(),
    ]
    if include_covariate_weights:
        lines.extend(
            [
                "",
                "Covariate weights:",
                _cov_weight_post_summary(model).to_string(),
            ]
        )
    return "\n".join(lines)


def prepare_park_locs(
    park_points: pd.DataFrame,
    analysis_start: pd.Timestamp = ANALYSIS_START,
) -> tuple[pd.DataFrame, pd.Timestamp, float]:
    """Build locs_s and season/temporal metadata for one park (copy to avoid SettingWithCopy)."""
    park_points = park_points.copy()
    coords = pd.DataFrame(
        {
            "x": park_points.geometry.centroid.x,
            "y": park_points.geometry.centroid.y,
        }
    )
    offset_seasonal = 30 * 0
    min_time = pd.Timestamp(analysis_start)
    time = (park_points["start_time"] - min_time).dt.total_seconds() / (24 * 60 * 60)
    seasonal = (time + offset_seasonal) % 365
    park_points["time_diff"] = time
    park_points["seasonal_diff"] = seasonal
    locs_s = pd.DataFrame(
        {
            "X": coords["x"].astype(float),
            "Y": coords["y"].astype(float),
            "T": time.astype(float),
            "A": seasonal.astype(float),
        }
    )
    return locs_s, min_time, offset_seasonal


def build_model(
    locs_s: pd.DataFrame,
    park_name: str,
    all_boxes_gdf,
    cov_gdf,
    total_days: int,
    offset_seasonal: float,
    spatial_sigma_prior_scale_m: float = DEFAULT_SPATIAL_SIGMA_PRIOR_SCALE_M,
):
    park_box = all_boxes_gdf[all_boxes_gdf["PARKNAME"] == park_name]
    cov_gdf_model, standardization_record = legacy_citywide_standardize_covariates(
        cov_gdf
    )
    model = bstpp.main.Hawkes_Model(
        locs_s,
        park_box,
        total_days,
        cox_background=True,
        spatial_cov=cov_gdf_model,
        cov_names=COLUMN_NAMES,
        standardize_cov=None,
        # PRE-S3 preview: preserve the old zero-valued behavior in the 3.44%
        # uncovered Tacony region, but surface/export the contract finding.
        data_contracts="report",
        # Preserve the historical bounding-rectangle excitation compensator.
        excitation_support="rectangle",
        a_0=dist.Normal(0, 2),
        alpha=dist.Beta(2, 2),
        beta=dist.Exponential(1.0),
        gamma=dist.Exponential(1.0),
        # PRE-S3 MIGRATION EDIT: sigmax_2 is a variance in square metres.
        sigmax_2=dist.HalfNormal(float(spatial_sigma_prior_scale_m) ** 2),
        offset_seasonal=offset_seasonal,
        temporal_trig=bstpp.trigger.Temporal_Power_Law,
    )
    # Preview-only provenance until results schema v1 lands at the S3 boundary.
    model.preview_standardization = standardization_record
    model.preview_spatial_sigma_prior_scale_m = float(spatial_sigma_prior_scale_m)
    return model


def _fit_title(
    park_name: str,
    num_steps: int,
    lr: float,
    plot_name: str | None = None,
) -> str:
    """Standard figure caption: park label, SVI steps, learning rate."""
    label = PARK_LABELS.get(park_name, park_name)
    base = f"{label} | {num_steps} steps | lr={lr:g}"
    return f"{plot_name} — {base}" if plot_name else base


def plot_svi_loss(
    model,
    park_name: str,
    num_steps: int,
    lr: float,
    *,
    out_path: Path | None = None,
    show: bool = False,
) -> None:
    loss = np.asarray(model.svi_results.losses)
    fig, ax = plt.subplots(figsize=(10, 4))
    start = int(0.01 * len(loss))
    if start < len(loss):
        ax.plot(np.arange(start, len(loss)), loss[start:])
    ax.set_xlabel("SVI iteration")
    ax.set_ylabel("Loss (−ELBO)")
    ax.set_title(_fit_title(park_name, num_steps, lr, "SVI loss"))
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if out_path is not None:
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def _save_current_fig(out_path: Path, title: str | None = None) -> None:
    fig = plt.gcf()
    if title:
        fig.suptitle(title, y=1.02)
    with contextlib.suppress(Exception):
        fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_all_diagnostic_figures(
    model,
    park_name: str,
    num_steps: int,
    lr: float,
    min_time: pd.Timestamp,
    out_dir: Path,
    window: float,
    num_samples: int,
    log_level: str = "normal",
) -> None:
    """Save full BSTPP diagnostic figures to disk (never displayed)."""
    def title(plot_name):
        return _fit_title(park_name, num_steps, lr, plot_name)

    def _path(name: str) -> Path:
        return _fig_path(out_dir, park_name, num_steps, lr, window, num_samples, name)

    if log_level == "verbose":
        _log(f"Saving diagnostic figures to {out_dir.resolve()}", log_level=log_level)

    with _suppress_plt_show():
        model.plot_prop_excitation()
        fig = plt.gcf()
        fig.suptitle(title("Proportion self-excitation"), y=1.02)
        plt.xlabel("Proportion of intensity due to self-excitation")
        _save_current_fig(_path("prop_excitation"))

        model.plot_trigger_posterior(trace=False)
        plt.gcf().suptitle(title("Trigger parameter posteriors"), y=1.02)
        _save_current_fig(_path("trigger_posterior"))

        model.plot_trigger_time_decay()
        fig = plt.gcf()
        fig.suptitle(title("Temporal trigger decay"), y=1.02)
        _save_current_fig(_path("trigger_time_decay"))

        model.plot_temporal(start_date=min_time)
        plt.gcf().suptitle(title("Temporal background $f_t$"), y=1.02)
        _save_current_fig(_path("temporal_ft"))

        model.plot_seasonal()
        plt.gcf().suptitle(title("Seasonal background $f_a$"), y=1.02)
        _save_current_fig(_path("seasonal_fa"))

        model.cov_weight_post_summary(trace=False)
        plt.gcf().suptitle(title("Covariate weights"), y=1.02)
        _save_current_fig(_path("covariate_weights"))
        plt.close("all")


def _plot_spatial_background_figure(
    model,
    park_name: str,
    num_steps: int,
    lr: float,
    window: float,
    num_samples: int,
    out_dir: Path | None,
    *,
    show: bool = False,
) -> None:
    """
    Spatial background using BSTPP's plot_spatial (grid × CBG overlay, 2 panels + colorbar).
    """
    out_path = (
        _fig_path(out_dir, park_name, num_steps, lr, window, num_samples, "spatial")
        if out_dir is not None
        else None
    )

    model.plot_spatial(include_cov=True, alpha=0.5)
    fig = plt.gcf()
    fig.suptitle(
        _fit_title(park_name, num_steps, lr, "Spatial background $f_{xy}$ + X(s)w"),
        y=1.02,
    )
    with contextlib.suppress(Exception):
        fig.tight_layout()
    if out_path is not None:
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def _plot_uncertainty_map(
    model,
    park_name: str,
    num_steps: int,
    lr: float,
    window: float,
    num_samples: int,
    out_dir: Path | None,
    *,
    parks_gdf=None,
    std_vmax: float | None = 0.33,
    show: bool = False,
) -> None:
    """
    Mean / std-dev map of posterior spatial field ``f_xy`` (notebook uncertainty map).

    Left: mean posterior ``f_s`` with events. Right: cell-wise posterior std.
    ``std_vmax`` matches the notebook default (0.33); pass ``None`` to auto-scale.
    """
    if "f_xy" not in model.samples:
        raise KeyError("Posterior has no 'f_xy'; cannot build uncertainty map.")

    out_path = (
        _fig_path(
            out_dir, park_name, num_steps, lr, window, num_samples, "uncertainty"
        )
        if out_dir is not None
        else None
    )

    f_xy_samples = np.asarray(model.samples["f_xy"])
    f_xy_mean = np.mean(f_xy_samples, axis=0)
    f_xy_std = np.std(f_xy_samples, axis=0)

    n_xy = int(model.args.get("n_xy", int(np.sqrt(f_xy_mean.size))))
    if n_xy * n_xy != f_xy_mean.size:
        raise ValueError(
            f"f_xy length {f_xy_mean.size} is not a square grid "
            f"(n_xy={n_xy})"
        )

    minx, miny, maxx, maxy = model.A[["geometry"]].total_bounds
    extent = [minx, maxx, miny, maxy]
    f_xy_mean_grid = f_xy_mean.reshape(n_xy, n_xy)
    f_xy_std_grid = f_xy_std.reshape(n_xy, n_xy)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    im1 = ax1.imshow(
        f_xy_mean_grid,
        origin="lower",
        cmap="viridis",
        extent=extent,
    )
    model.points.plot(ax=ax1, color="red", marker="x", alpha=0.1)
    if parks_gdf is not None:
        parks_gdf.plot(ax=ax1, edgecolor="white", facecolor="none", linewidth=2)
    ax1.set_title("Mean Posterior $f_s$ With Events")
    ax1.set_xlim(minx, maxx)
    ax1.set_ylim(miny, maxy)
    fig.colorbar(im1, ax=ax1, shrink=0.7)

    im2_kw: dict = {
        "origin": "lower",
        "cmap": "RdYlGn_r",
        "extent": extent,
        "vmin": 0,
    }
    if std_vmax is not None:
        im2_kw["vmax"] = std_vmax
    im2 = ax2.imshow(f_xy_std_grid, **im2_kw)
    if parks_gdf is not None:
        parks_gdf.plot(ax=ax2, edgecolor="gray", facecolor="none", linewidth=2)
    ax2.set_title("$f_s$ Uncertainty (Std Dev)")
    ax2.set_xlim(minx, maxx)
    ax2.set_ylim(miny, maxy)
    fig.colorbar(im2, ax=ax2, shrink=0.7)

    fig.suptitle(
        _fit_title(park_name, num_steps, lr, "Spatial field uncertainty"),
        y=1.02,
    )
    fig.tight_layout()
    if out_path is not None:
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def save_fit_outputs(
    model,
    locs_s: pd.DataFrame,
    park_name: str,
    num_steps: int,
    lr: float,
    min_time: pd.Timestamp,
    out_dir: Path,
    window: float,
    spatial_window: float,
    num_samples: int,
    n_pairs: int | None = None,
    *,
    parks_gdf=None,
    show_svi_loss: bool = False,
    show_maps: bool = False,
    show_summary: bool = False,
    log_level: str = "normal",
) -> dict:
    """Save all figures; optionally display SVI loss, spatial maps, and text summary."""
    label = PARK_LABELS.get(park_name, park_name)
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_svi_loss(
        model,
        park_name,
        num_steps,
        lr,
        out_path=_fig_path(
            out_dir, park_name, num_steps, lr, window, num_samples, "svi_loss"
        ),
        show=show_svi_loss,
    )

    save_all_diagnostic_figures(
        model,
        park_name,
        num_steps,
        lr,
        min_time,
        out_dir,
        window,
        num_samples,
        log_level=log_level,
    )

    summary_path = _artifact_path(
        out_dir, park_name, num_steps, lr, window, num_samples, "params", "txt"
    )
    summary_path.write_text(
        build_fit_summary_text(
            model,
            locs_s,
            park_name=park_name,
            num_steps=num_steps,
            lr=lr,
            window=window,
            spatial_window=spatial_window,
            num_samples=num_samples,
            n_pairs=n_pairs,
            include_covariate_weights=True,
            key_trigger_means_only=False,
        ),
        encoding="utf-8",
    )
    if show_summary and log_level != "quiet":
        print(
            build_fit_summary_text(
                model,
                locs_s,
                park_name=park_name,
                num_steps=num_steps,
                lr=lr,
                window=window,
                spatial_window=spatial_window,
                num_samples=num_samples,
                n_pairs=n_pairs,
                include_covariate_weights=False,
                key_trigger_means_only=True,
            )
        )

    _plot_spatial_background_figure(
        model,
        park_name,
        num_steps,
        lr,
        window,
        num_samples,
        out_dir,
        show=show_maps,
    )

    _plot_uncertainty_map(
        model,
        park_name,
        num_steps,
        lr,
        window,
        num_samples,
        out_dir,
        parks_gdf=parks_gdf,
        show=show_maps,
    )

    trig_df = _trigger_post_summary(model)
    prop_excitation = float(
        (model.samples["Itot_excite"] / model.samples["Itot_txy"]).mean()
    )
    loglik = float(model.log_expected_likelihood(locs_s))
    aic = float(model.expected_AIC())
    loss = np.asarray(model.svi_results.losses)

    row = {
        "park_name": park_name,
        "park_label": label,
        "num_steps": num_steps,
        "lr": lr,
        "window": window,
        "spatial_window": spatial_window,
        "spatial_sigma_prior_scale_m": model.preview_spatial_sigma_prior_scale_m,
        "num_samples": num_samples,
        "n_events": int(len(locs_s)),
        "n_pairs": n_pairs,
        "prop_excitation_mean": prop_excitation,
        "log_expected_likelihood": loglik,
        "expected_aic": aic,
        "final_loss": float(loss[-1]),
        "any_nan_loss": bool(np.isnan(loss).any()),
        "params_txt": str(summary_path.resolve()),
    }
    for name in trig_df.index:
        row[f"trigger_{name}_mean"] = float(trig_df.loc[name, "Post Mean"])
    return row


def fit_one_park(
    park_name: str,
    illegal_dumping_in_park: pd.DataFrame,
    all_boxes_gdf,
    cov_gdf,
    total_days: int,
    num_steps: int,
    lr: float,
    out_dir: Path,
    fit_kwargs: dict | None = None,
    defaults: dict | None = None,
    park_overrides: dict[str, dict] | None = None,
    parks_gdf=None,
    job_index: int | None = None,
    job_total: int | None = None,
    log_level: str = "normal",
    show_svi_loss: bool = True,
    show_maps: bool = True,
    show_summary: bool = True,
) -> dict:
    label = PARK_LABELS.get(park_name, park_name)
    job_hdr = ""
    if job_index is not None and job_total is not None:
        job_hdr = f"[{job_index}/{job_total}] "

    t0 = time.perf_counter()
    park_points = illegal_dumping_in_park[
        illegal_dumping_in_park["PARKNAME"] == park_name
    ].copy()
    n_events = len(park_points)
    if n_events < 10:
        raise ValueError(f"{park_name}: only {n_events} events after filtering")

    kw = get_park_fit_kwargs(
        park_name, defaults=defaults, park_overrides=park_overrides, **(fit_kwargs or {})
    )
    window, window_days, window_source = resolve_temporal_window(kw, total_days)
    spatial_window = float(kw["spatial_window"])
    spatial_sigma_prior_scale_m = float(kw["spatial_sigma_prior_scale_m"])
    num_samples = int(kw["num_samples"])

    _log(
        f"{job_hdr}{label} | {num_steps} steps | n={n_events:,} | "
        f"temporal={window_days:g} d "
        f"(internal={window:g} via {window_source}) | "
        f"spatial_window={spatial_window:g} m | "
        f"sigma-prior-scale={spatial_sigma_prior_scale_m:g} m | ns={num_samples}",
        log_level=log_level,
    )
    locs_s, min_time, offset_seasonal = prepare_park_locs(park_points)

    model = build_model(
        locs_s,
        park_name,
        all_boxes_gdf,
        cov_gdf,
        total_days,
        offset_seasonal,
        spatial_sigma_prior_scale_m=spatial_sigma_prior_scale_m,
    )

    with _suppress_bstpp_chatter(log_level != "verbose"):
        model.set_window(window=window, spatial_window=spatial_window)
    n_pairs = int(model.args["coords"].shape[0])
    if log_level == "verbose":
        _log(f"  Trigger pairs: {n_pairs:,}", log_level=log_level)

    _log(f"{job_hdr}SVI {num_steps} steps + {num_samples} draws…", log_level=log_level)
    t_svi = time.perf_counter()
    model.run_svi(
        lr=lr,
        num_steps=num_steps,
        num_samples=num_samples,
        plot_loss=False,
    )
    ensure_b_0_in_samples(model)
    svi_sec = time.perf_counter() - t_svi

    loss = np.asarray(model.svi_results.losses)
    final_loss = float(loss[-1])
    nan_loss = bool(np.isnan(loss).any())
    _log(
        f"{job_hdr}SVI done {svi_sec / 60:.1f} min | loss={final_loss:.1f}"
        + (f" | pairs={n_pairs:,}" if log_level == "normal" else "")
        + (" | NaN in loss" if nan_loss else ""),
        log_level=log_level,
    )

    if log_level == "normal":
        _log(f"{job_hdr}Saving figures…", log_level=log_level)

    row = save_fit_outputs(
        model,
        locs_s,
        park_name,
        num_steps,
        lr,
        min_time,
        out_dir,
        window,
        spatial_window,
        num_samples,
        n_pairs=n_pairs,
        parks_gdf=parks_gdf,
        show_svi_loss=show_svi_loss,
        show_maps=show_maps,
        show_summary=show_summary,
        log_level=log_level,
    )

    elapsed = time.perf_counter() - t0
    _log(
        f"{job_hdr}Done {label} ({num_steps} steps) {elapsed / 60:.1f} min | "
        f"AIC={row['expected_aic']:.0f}",
        log_level=log_level,
    )

    # Free model state before the next fit (caller may also clear JAX caches)
    del model
    gc.collect()

    return row


def run_all_park_fits(
    illegal_dumping_in_park,
    all_boxes_gdf,
    cov_gdf,
    total_days: int,
    output_dir: str | Path = "output/batch_fits",
    timestamped: bool = True,
    park_names: list[str] | None = None,
    num_steps_list: list[int] | None = None,
    park_num_steps: dict[str, int | list[int]] | None = None,
    lr: float = 0.01,
    default_fit_kwargs: dict | None = None,
    park_fit_overrides: dict[str, dict] | None = None,
    temporal_cutoff_days: float | None = None,
    parks_gdf=None,
    log_level: str = "normal",
    verbose: bool | None = None,
    show_svi_loss: bool = True,
    show_maps: bool = True,
    show_summary: bool = True,
    clear_jax_cache: bool = True,
) -> pd.DataFrame:
    """
    Fit all parks × step counts; save outputs under output_dir.

    Per fit, all diagnostic PNGs are saved (including the ``f_xy`` uncertainty map);
    only SVI loss, spatial maps, and the text parameter summary are shown inline
    (plus ``fit_summary.csv`` for the batch).

    Optional ``parks_gdf`` (e.g. ``all_parks_gdf``) is overlaid on the uncertainty
    map; omit it to skip park outlines.

    When timestamped=True (default), each call writes to a new subfolder
    ``output_dir/YYYYMMDD_HHMMSS/`` so runs do not overwrite prior outputs.

    Global default SVI steps: ``num_steps_list`` (default ``[6000, 10000]``) for parks
    not listed in ``park_num_steps``.

    Per-park steps (preferred): pass ``park_num_steps`` with park name → int or list,
    e.g. ``{"mifflin_box": 6000, "tacony_box": [6000, 10000]}``.

    Temporal computational cutoff: pass ``temporal_cutoff_days`` (real days) as a
    batch-wide default, and/or set ``temporal_cutoff_days`` / legacy ``window`` in
    ``default_fit_kwargs`` or ``park_fit_overrides``. Real-day values win over the
    inherited legacy ``window`` default when both would otherwise be present.

    Other per-park Hawkes settings (``spatial_window``, ``num_samples``, …) use
    ``park_fit_overrides``.

    When ``clear_jax_cache`` is True (default), call ``jax.clear_caches`` at batch
    start and after each fit so compile caches do not accumulate across parks.
    """
    log_level = _resolve_log_level(verbose, log_level)

    defaults = {**DEFAULT_FIT_KWARGS, **(default_fit_kwargs or {})}
    if temporal_cutoff_days is not None:
        defaults["temporal_cutoff_days"] = float(temporal_cutoff_days)
        defaults.pop("window", None)
    park_overrides = {**PARK_FIT_OVERRIDES, **(park_fit_overrides or {})}
    park_num_steps = park_num_steps or {}

    base_dir = Path(output_dir)
    base_dir.mkdir(parents=True, exist_ok=True)
    out_dir = _timestamped_run_dir(base_dir) if timestamped else base_dir
    run_started_at = datetime.now().isoformat(timespec="seconds")
    park_names = park_names or PARK_NAMES
    num_steps_list = num_steps_list or [6000, 10000]

    jobs = build_batch_jobs(
        park_names, num_steps_list, park_num_steps, park_overrides
    )
    resolved_park_num_steps = {
        p: num_steps_for_park(p, num_steps_list, park_num_steps, park_overrides)
        for p in park_names
    }

    _, standardization_record = legacy_citywide_standardize_covariates(cov_gdf)
    config = {
        "output_base_dir": str(base_dir.resolve()),
        "run_dir": str(out_dir.resolve()),
        "run_started_at": run_started_at,
        "timestamped": timestamped,
        "park_names": park_names,
        "num_steps_list": num_steps_list,
        "park_num_steps": resolved_park_num_steps,
        "park_num_steps_arg": park_num_steps,
        "jobs": [{"park_name": p, "num_steps": s} for p, s in jobs],
        "lr": lr,
        "default_fit_kwargs": defaults,
        "park_fit_overrides": park_overrides,
        "temporal_cutoff_days_arg": temporal_cutoff_days,
        "clear_jax_cache": clear_jax_cache,
        "total_days": total_days,
        "standardization": standardization_record,
        "covariate_crs": str(cov_gdf.crs),
    }
    with open(out_dir / "run_config.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    rows = []
    errors = []
    n_jobs = len(jobs)

    _log(
        f"Batch: {n_jobs} fits → {out_dir.resolve()}",
        log_level=log_level,
    )
    if clear_jax_cache:
        release_jax_memory(log_level=log_level)

    batch_t0 = time.perf_counter()
    for job_i, (park_name, num_steps) in enumerate(jobs, start=1):
        try:
            row = fit_one_park(
                park_name,
                illegal_dumping_in_park,
                all_boxes_gdf,
                cov_gdf,
                total_days,
                num_steps,
                lr,
                out_dir,
                defaults=defaults,
                park_overrides=park_overrides,
                parks_gdf=parks_gdf,
                job_index=job_i,
                job_total=n_jobs,
                log_level=log_level,
                show_svi_loss=show_svi_loss,
                show_maps=show_maps,
                show_summary=show_summary,
            )
            rows.append(row)
            if clear_jax_cache:
                release_jax_memory(log_level="quiet")
        except Exception as exc:
            errors.append(
                {
                    "park_name": park_name,
                    "num_steps": num_steps,
                    "lr": lr,
                    "error": repr(exc),
                }
            )
            _log(f"FAILED [{job_i}/{n_jobs}] {park_name} {num_steps} steps: {exc}", log_level="normal")
            if clear_jax_cache:
                release_jax_memory(log_level="quiet")
            else:
                gc.collect()

    summary = pd.DataFrame(rows)
    if len(summary):
        summary.insert(0, "run_dir", str(out_dir.resolve()))
        summary.insert(1, "run_started_at", run_started_at)
    summary_path = out_dir / "fit_summary.csv"
    summary.to_csv(summary_path, index=False)

    batch_min = (time.perf_counter() - batch_t0) / 60
    if errors:
        err_path = out_dir / "fit_errors.csv"
        pd.DataFrame(errors).to_csv(err_path, index=False)
        _log(
            f"Finished {batch_min:.1f} min — {len(rows)}/{n_jobs} OK, {len(errors)} failed → {err_path}",
            log_level=log_level,
        )
    else:
        _log(f"Finished {batch_min:.1f} min — all {n_jobs} OK → {summary_path}", log_level=log_level)

    if log_level == "verbose" and len(summary):
        cols = ["park_name", "num_steps", "n_events", "final_loss", "expected_aic"]
        _log("Completed runs:\n" + summary[cols].to_string(index=False), log_level=log_level)
    return summary
