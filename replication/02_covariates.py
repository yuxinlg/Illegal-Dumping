"""
02_covariates.py
----------------
Build the CBG-level covariate panel used as spatial background intensity inputs
for the Hawkes model.

Inputs  (raw data, not tracked in git -- see README.md):
    data/tiger/tl_2023_42_bg/          -- Census block group boundaries
    data/Land_Care/phs_landcare.geojson
    data/Land_Use/Land_use.geojson
    data/Vacant_Block_Percent_Land.geojson
    data/NDVI/                         -- HLS Sentinel-2 raster tiles
    data/foot_traffic/                 -- SafeGraph monthly foot-traffic CSVs
    output/philly_census_data_21_23.geojson  -- ACS 5-yr estimates (from acs.ipynb)
    output/all_parks_gdf.geojson       -- park polygons (from 03_analysis setup)

Output:
    output/cov_cbg.geojson             -- CBG-level covariate table with geometry
"""
