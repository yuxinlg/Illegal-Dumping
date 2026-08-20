"""Construct and optionally fit the PRE-S3 cbg-park-seasonal preview model."""

# ruff: noqa: E402 - load BSTPP preview and this replication dir before imports.

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("JAX_PLATFORM_NAME", "cpu")

REPL_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = REPL_DIR.parent
PREVIEW_ROOT = Path(
    os.environ.get(
        "BSTPP_PREVIEW_ROOT",
        ANALYSIS_ROOT.parent / "BSTPP_preview",
    )
).resolve()
if str(PREVIEW_ROOT) not in sys.path:
    sys.path.insert(0, str(PREVIEW_ROOT))
if str(REPL_DIR) not in sys.path:
    sys.path.insert(0, str(REPL_DIR))

import geopandas as gpd
import numpy as np
import pandas as pd

import bstpp
from batch_park_fits import (
    ANALYSIS_END,
    ANALYSIS_START,
    ANALYSIS_TOTAL_DAYS,
    DEFAULT_SPATIAL_SIGMA_PRIOR_SCALE_M,
    SPATIAL_WINDOW_M,
    build_model,
    prepare_park_locs,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=ANALYSIS_ROOT,
        help="Illegal-Dumping checkout containing output/illegal_dumping_full.geojson, "
        "output/all_boxes_gdf.geojson, and output/cov_cbg.geojson.",
    )
    parser.add_argument("--park", default="tacony_box")
    parser.add_argument("--total-days", type=int, default=ANALYSIS_TOTAL_DAYS)
    parser.add_argument("--temporal-window", type=float, default=100.0)
    parser.add_argument(
        "--spatial-window-m", type=float, default=SPATIAL_WINDOW_M
    )
    parser.add_argument(
        "--spatial-sigma-prior-scale-m",
        type=float,
        default=DEFAULT_SPATIAL_SIGMA_PRIOR_SCALE_M,
    )
    parser.add_argument("--svi-steps", type=int, default=0)
    parser.add_argument("--num-samples", type=int, default=20)
    parser.add_argument("--lr", type=float, default=0.01)
    return parser.parse_args()


def _load_analysis_inputs(analysis_root: Path):
    output = analysis_root.resolve() / "output"
    illegal_dumping = gpd.read_file(output / "illegal_dumping_full.geojson")
    illegal_dumping["start_time"] = pd.to_datetime(illegal_dumping["start_time"])
    illegal_dumping = illegal_dumping[
        (illegal_dumping["start_time"] >= ANALYSIS_START)
        & (illegal_dumping["start_time"] < ANALYSIS_END)
    ].copy()

    boxes = gpd.read_file(output / "all_boxes_gdf.geojson").to_crs(26918)
    covariates = gpd.read_file(output / "cov_cbg.geojson").to_crs(26918)
    illegal_dumping = illegal_dumping.to_crs(boxes.crs)

    in_park = illegal_dumping.sjoin(boxes, predicate="within")
    in_park = in_park.drop(columns=["index_right"])
    in_park = in_park.sjoin(covariates, predicate="within", how="inner")
    return in_park, boxes, covariates


def main() -> None:
    args = _parse_args()
    in_park, boxes, covariates = _load_analysis_inputs(args.analysis_root)
    park_points = in_park[in_park["PARKNAME"] == args.park].copy()
    if len(park_points) < 10:
        raise ValueError(f"{args.park}: only {len(park_points)} events after filtering")

    locs_s, _min_time, offset_seasonal = prepare_park_locs(park_points)
    model = build_model(
        locs_s,
        args.park,
        boxes,
        covariates,
        args.total_days,
        offset_seasonal,
        spatial_sigma_prior_scale_m=args.spatial_sigma_prior_scale_m,
    )
    model.set_window(
        window=args.temporal_window,
        spatial_window=args.spatial_window_m,
    )

    standardization = model.preview_standardization
    print(f"bstpp_file={bstpp.__file__}")
    print(f"park={args.park} events={len(locs_s)} pairs={model.args['coords'].shape[0]}")
    print(
        "standardization="
        f"{standardization['method']} rows={standardization['row_count']}"
    )
    print(
        f"spatial_window_m={args.spatial_window_m:g} "
        f"spatial_sigma_prior_scale_m={args.spatial_sigma_prior_scale_m:g}"
    )
    assert standardization["row_count"] == len(covariates)
    assert np.isclose(float(model.args["spatial_window"]), args.spatial_window_m)

    if args.svi_steps > 0:
        model.run_svi(
            num_steps=args.svi_steps,
            lr=args.lr,
            num_samples=args.num_samples,
            plot_loss=False,
        )
        losses = np.asarray(model.svi_results.losses)
        if not np.isfinite(losses).all():
            raise RuntimeError("SVI produced nonfinite losses")
        print(
            f"svi_steps={args.svi_steps} samples={args.num_samples} "
            f"final_loss={float(losses[-1]):.6g}"
        )


if __name__ == "__main__":
    main()
