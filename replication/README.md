# Replication Package

This folder contains the replication code for:

> **[Paper title]**  
> [Authors] · [Journal, Year]

---

## Overview

The analysis models illegal-dumping reports in four Philadelphia parks using a
Cox-Hawkes spatiotemporal point process fitted with the BSTPP package.
Scripts must be run in the following order:


| Step | Script                          | Description                                                                                                                                         |
| ---- | ------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1    | `01_data_cleaning.py`           | Load public 311 shapefiles (2021–2025), filter to "Illegal Dumping", build park boxes, export `output/illegal_dumping_full.geojson`, `locs_s.csv`, `all_parks_gdf.geojson`, `all_boxes_gdf.geojson` |
| 2a   | `02a_acs_prep.py`               | Build ACS 5-year covariate panel (2021–2024), export `output/philly_cbg_acs.geojson`                                                               |
| 2b   | `02b_network_connectivity.py`   | **One-time, expensive (2–6 h).** Download Philadelphia walking network (OSMnx), compute centrality metrics, export `output/philadelphia_walk_network_connectivity.geojson`. Skip if uploading the pre-computed file directly. |
| 2c   | `02_covariates.py`              | Assemble full CBG-level covariate panel (ACS, land care, land use, vacant parcels, NDVI, foot traffic\*, walk-network), export `output/cov_cbg.geojson` |
| 3    | `03_analysis.py`                | Fit Cox-Hawkes model per park, generate all diagnostic figures                                                                                      |

\* **Foot-traffic data is proprietary** (SafeGraph Neighborhood Patterns via Dewey Data). Researchers must obtain their own licensed access. `02_covariates.py` will skip foot-traffic and reporting-rate sections gracefully if `data/foot_traffic/` is empty.


---

## Setup

### 1. Python environment

```bash
pip install -r requirements.txt
```

### 2. BSTPP package

This study uses the customized BSTPP preview checkout (`BSTPP_preview`),
not the stock PyPI package. `03_analysis.py` prepends that checkout to
`sys.path` and constructs Cox–Hawkes as
`Hawkes_Model(..., cox_background=True)`.

The old `cox_hawkes_shared` module is **not** on GitHub and is not used.

Default location (sibling of this analysis repo):

```text
C:\Users\Terhi\dev\BSTPP_preview
```

Override with the environment variable `BSTPP_PREVIEW_ROOT`, or install
editable from that checkout:

```bash
pip install -e C:\Users\Terhi\dev\BSTPP_preview
```

---

## Data

Raw data are **not tracked in this repository**. Obtain the following and place
them in the paths listed:


| Dataset | Source | Local path | Access Time |
| ------- | ------ | ---------- |-------------|
| Public Philly 311 Records (2021–2025) | [https://opendataphilly.org/datasets/311-service-and-information-requests/](https://opendataphilly.org/datasets/311-service-and-information-requests/) | `data/311/` | - |
| American Community Survey ()
| Census block group boundaries (TIGER/Line 2023) | [census.gov](https://www.census.gov/cgi-bin/geo/shapefiles/index.php) | `data/tiger/tl_2023_42_bg/` | - |
| Philadelphia City Limits | [OpenDataPhilly](https://opendataphilly.org/datasets/city-limits/) | `data/City_Limits/` | - |
| Land Care | [OpenDataPhilly](https://data-phl.opendata.arcgis.com/datasets/phl::phs-landcare/explore?location=39.972982%2C-75.153989%2C13.29) | - |July 17, 2025 |
| Land Use (current) | [OpenDataPhilly](https://opendataphilly.org/datasets/land-use/) | `data/Land_Use/` |March 28, 2025|
| Vacant Block| [OpenDataPhilly](https://data-phl.opendata.arcgis.com/datasets/phl::vacant-block-percent-land/) | `data/Vacant_Block_Percent_Land.geojson` |June 11, 2025|
| HLS Sentinel-2 NDVI tiles | [NASA Earthdata](https://www.earthdata.nasa.gov/data/catalog/lpcloud-hlss30-2.0) | `data/NDVI/` |July 15, 2025|
| PPR Property Boundaries | [OpenDataPhilly](https://opendataphilly.org/datasets/ppr-properties/) | `data/parks/PPR_Properties.geojson` |-|
| Park AOI polygons | Google My Maps (provided by research team) | `data/parks/park_AOI_from_google.kml` | Feb 3, 2026 |
| Network Connectivity (walk) | OpenStreetMap (OSM), accessed via OSMnx | `output/philadelphia_walk_network_connectivity.geojson` (pre-computed, upload directly) or regenerate with `02b_network_connectivity.py` | October 15, 2025 |

<!-- | SafeGraph foot-traffic (2021–2024) | SafeGraph (licensed) | `data/foot_traffic/` |-| -->
<!-- | PPR base GDB (hotspots layer) | Provided by Philadelphia Parks & Recreation | `data/PPR_BaseFiles_IllegalDumpingModel.gdb` |-| -->