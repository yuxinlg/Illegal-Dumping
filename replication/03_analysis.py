"""
03_analysis.py
--------------
Fit a Cox-Hawkes spatiotemporal point process model (via BSTPP) for each of the
four study parks (Cobbs Creek, Mifflin Square, Tacony Creek, Fairmount Park)
and produce all diagnostic figures saved to output/ and output_litteridx_emb/.

Inputs:
    output/illegal_dumping_full.geojson   -- from 01_data_cleaning.py
    output/cov_cbg.geojson                -- from 02_covariates.py
    parks/PPR/PPR_Properties.geojson      -- PPR park boundaries (see README.md)
    data/City_Limits/City_Limits.shp      -- Philadelphia boundary
    data/tiger/tl_2023_42_bg/             -- CBG boundaries
    data/Land_Use/Land_use.geojson
    data/Land_Care/phs_landcare.geojson
    data/Vacant_Block_Percent_Land.geojson
    data/PPR_BaseFiles_IllegalDumpingModel.gdb  -- PPR dumping hotspots layer
    output/camera_gdf.geojson             -- surveillance camera locations

Outputs  (written to output_litteridx_emb/<park>_*.png):
    seasonal trend, temporal intensity, spatial trigger, covariate coefficients,
    proportion of excitation, trigger time decay, uncertainty surface, etc.

Dependencies:
    bstpp  -- see README.md for installation of the modified version
"""
