"""Focused checks for the provisional cbg-park-seasonal analysis adapter."""

import os
import sys
from pathlib import Path

_REPL = Path(__file__).resolve().parent
_ANALYSIS = _REPL.parent
_PREVIEW = Path(
    os.environ.get("BSTPP_PREVIEW_ROOT", _ANALYSIS.parent / "BSTPP_preview")
).resolve()
if str(_PREVIEW) not in sys.path:
    sys.path.insert(0, str(_PREVIEW))
if str(_REPL) not in sys.path:
    sys.path.insert(0, str(_REPL))

import numpy as np
import pandas as pd
import pytest

from batch_park_fits import (
    ANALYSIS_END,
    ANALYSIS_START,
    ANALYSIS_TOTAL_DAYS,
    COLUMN_NAMES,
    DEFAULT_SPATIAL_SIGMA_PRIOR_SCALE_M,
    SPATIAL_SIGMA_SENSITIVITY_M,
    SPATIAL_WINDOW_M,
    get_park_fit_kwargs,
    legacy_citywide_standardize_covariates,
    resolve_temporal_window,
)
from bstpp.cutoffs import days_to_internal
from bstpp.preparation import T_INTERNAL


def test_legacy_standardization_uses_all_supplied_rows_and_population_scale():
    values = np.arange(5 * len(COLUMN_NAMES), dtype=float).reshape(5, -1)
    covariates = pd.DataFrame(values, columns=COLUMN_NAMES)

    standardized, record = legacy_citywide_standardize_covariates(covariates)

    expected_mean = values.mean(axis=0)
    expected_scale = values.var(axis=0) ** 0.5
    np.testing.assert_allclose(
        standardized[COLUMN_NAMES].to_numpy(),
        (values - expected_mean) / expected_scale,
    )
    np.testing.assert_allclose(record["mean"], expected_mean)
    np.testing.assert_allclose(record["scale"], expected_scale)
    assert record["row_count"] == len(covariates)
    assert record["columns"] == COLUMN_NAMES


def test_preview_defaults_record_real_unit_spatial_choices():
    settings = get_park_fit_kwargs("tacony_box")

    assert settings["spatial_window"] == SPATIAL_WINDOW_M == 1_000.0
    assert settings["spatial_sigma_prior_scale_m"] == (
        DEFAULT_SPATIAL_SIGMA_PRIOR_SCALE_M
    ) == 250.0
    assert SPATIAL_SIGMA_SENSITIVITY_M == (50.0, 100.0, 250.0)


def test_preview_calendar_horizon_includes_the_2024_leap_day():
    assert ANALYSIS_START == pd.Timestamp("2021-01-01")
    assert ANALYSIS_END == pd.Timestamp("2025-01-01")
    assert ANALYSIS_TOTAL_DAYS == 1_461


def test_resolve_temporal_window_from_real_days():
    window, days, source = resolve_temporal_window(
        {"temporal_cutoff_days": 30.0}, ANALYSIS_TOTAL_DAYS
    )
    assert source == "temporal_cutoff_days"
    assert days == 30.0
    np.testing.assert_allclose(
        window, days_to_internal(30.0, ANALYSIS_TOTAL_DAYS)
    )
    assert window < T_INTERNAL  # finite cutoff inside the horizon


def test_resolve_temporal_window_from_legacy_internal():
    window, days, source = resolve_temporal_window(
        {"window": 100.0}, ANALYSIS_TOTAL_DAYS
    )
    assert source == "window"
    assert window == 100.0
    np.testing.assert_allclose(
        days, 100.0 * (ANALYSIS_TOTAL_DAYS / T_INTERNAL)
    )


def test_resolve_temporal_window_rejects_both_forms():
    with pytest.raises(ValueError, match="temporal_cutoff_days or window"):
        resolve_temporal_window(
            {"window": 100.0, "temporal_cutoff_days": 30.0},
            ANALYSIS_TOTAL_DAYS,
        )


def test_get_park_fit_kwargs_prefers_temporal_cutoff_days_over_default_window():
    settings = get_park_fit_kwargs(
        "tacony_box",
        defaults={
            "window": 100.0,
            "temporal_cutoff_days": 45.0,
            "spatial_window": SPATIAL_WINDOW_M,
            "spatial_sigma_prior_scale_m": 250.0,
            "num_samples": 750,
        },
    )
    assert "window" not in settings
    assert settings["temporal_cutoff_days"] == 45.0


def test_get_park_fit_kwargs_park_days_override_beats_default_window():
    settings = get_park_fit_kwargs(
        "mifflin_box",
        park_overrides={"mifflin_box": {"temporal_cutoff_days": 14.0}},
    )
    assert "window" not in settings
    assert settings["temporal_cutoff_days"] == 14.0
