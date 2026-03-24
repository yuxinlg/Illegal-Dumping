# Replication Package

This folder contains the replication code for:

> **[Paper title]**  
> [Authors] · [Journal, Year]

---

## Overview

The analysis models illegal-dumping reports in four Philadelphia parks using a
Cox-Hawkes spatiotemporal point process fitted with the BSTPP package.
Three scripts must be run in order:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_data_cleaning.py` | Load and merge raw PPR 311 records (2018–2025), apply case-type filters and manual labels, export `output/illegal_dumping_full.geojson` |
| 2 | `02_covariates.py` | Build CBG-level covariate panel (ACS demographics, land use, NDVI, foot traffic, street-network metrics), export `output/cov_cbg.geojson` |
| 3 | `03_analysis.py` | Fit Cox-Hawkes model per park, generate all diagnostic figures to `output_litteridx_emb/` |

---

## Setup

### 1. Python environment

```bash
pip install -r requirements.txt
```

### 2. BSTPP package

This study uses a modified version of the BSTPP package. Install it separately
before running the scripts (instructions to be added once the package repo is
finalised):

```bash
# pip install git+https://github.com/<your-org>/BSTPP.git
```

---

## Data

Raw data are **not tracked in this repository**. Obtain the following and place
them in the paths listed:

| Dataset | Source | Local path |
|---------|--------|------------|
| PPR 311 illegal-dumping records (2018–2025) | Provided by Philadelphia Parks & Recreation | `data/311/` |
| Census block group boundaries (TIGER/Line 2023) | [census.gov](https://www.census.gov/cgi-bin/geo/shapefiles/index.php) | `data/tiger/tl_2023_42_bg/` |
| Philadelphia City Limits | [OpenDataPhilly](https://www.opendataphilly.org/) | `data/City_Limits/` |
| Land Use (current) | [OpenDataPhilly](https://www.opendataphilly.org/) | `data/Land_Use/` |
| PHS LandCare parcels | [OpenDataPhilly](https://www.opendataphilly.org/) | `data/Land_Care/` |
| Vacant block % land | [OpenDataPhilly](https://www.opendataphilly.org/) | `data/Vacant_Block_Percent_Land.geojson` |
| HLS Sentinel-2 NDVI tiles | [NASA Earthdata](https://search.earthdata.nasa.gov/) | `data/NDVI/` |
| SafeGraph foot-traffic (2021–2024) | SafeGraph (licensed) | `data/foot_traffic/` |
| PPR base GDB (hotspots layer) | Provided by Philadelphia Parks & Recreation | `data/PPR_BaseFiles_IllegalDumpingModel.gdb` |
| PPR property boundaries | [OpenDataPhilly](https://www.opendataphilly.org/) | `parks/PPR/PPR_Properties.geojson` |
