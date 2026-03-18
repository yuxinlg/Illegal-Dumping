from dill.temp import b
import os
os.environ["JAX_PLATFORM_NAME"] = "cpu"

import jax
print("JAX path:", jax.__file__)

## load module
import sys

bstpp_module_dir = os.path.abspath('/Users/ctwomey/Research/WPF Illegal Dumping/project/Elena/BSTPP')
#bstpp_module_dir = os.path.abspath('/Users/ctwomey/Research/WPF Illegal Dumping/project/BSTPP')
sys.path.insert(0, bstpp_module_dir)
import bstpp

import numpy as np
import numpyro
import numpyro.distributions as dist
import pandas as pd
import geopandas as gpd
from datetime import timedelta
import matplotlib.pyplot as plt
import seaborn as sns
#!pip install torch
import torch
from shapely.geometry import box
np.random.seed(999)


from bstpp.main import Hawkes_Model


## set working directory
import os
os.chdir("/Users/ctwomey/Research/WPF Illegal Dumping/project/Elena/Illegal-Dumping-main")



## Pre settings
# - Cobbs Creek: #e6b422
# - Mifflin Square: #F80E07
# - Tacony Creek: #19448e
# - Fairmount: #543f32

## 311 data preprocessing

## load 2024 data
litter_24 = gpd.read_file('data/311/2024/public_cases_fc.shp')
## filter data
illegal_dumping_24 = litter_24[litter_24['service_na'] == "Illegal Dumping"]

## load 2023 data
litter_23 = gpd.read_file('data/311/2023/public_cases_fc.shp')
## filter data
illegal_dumping_23 = litter_23[litter_23['service_na'] == "Illegal Dumping"]

## load 2022 data
litter_22 = gpd.read_file('data/311/2022/public_cases_fc.shp')
## filter data
illegal_dumping_22 = litter_22[litter_22['service_na'] == "Illegal Dumping"]

## load 2021 data
litter_21 = gpd.read_file('data/311/2021/public_cases_fc.shp')
## filter data
illegal_dumping_21 = litter_21[litter_21['service_na'] == "Illegal Dumping"]

# NOTE: something wrong with this file, seems to be a dup of 2022
## load 2020 data
litter_20 = gpd.read_file('data/311/2020/public_cases_fc.shp')
## filter data
illegal_dumping_20 = litter_20[litter_20['service_na'] == "Illegal Dumping"]


## concact 2022, 2023, and 2024 data
illegal_dumping = pd.concat([
    illegal_dumping_20,
    illegal_dumping_21,
    illegal_dumping_22,
    illegal_dumping_23,
    illegal_dumping_24
])
#illegal_dumping = illegal_dumping_23
total_days = 365 * 5



### Loadding park boundaries

## load philly's city map for mapping
philly_border = gpd.read_file('data/City_Limits/City_Limits.shp')

## load PPR data
ppr = gpd.read_file('parks/PPR/PPR_Properties.geojson')
ppr = ppr.to_crs(4326)

## Tacony
tacony = ppr[ppr["OBJECTID"].isin([10, 64])]

## current park
park_name = "Tacony"
park = tacony



## plot park

fig, ax = plt.subplots(figsize=(10, 10))

philly_border.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=1)
park.plot(ax=ax, color='#19448e', label=park_name)

plt.legend()
plt.title('Selected Park in Philadelphia')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()




# set buffer value
buffer = 0.005 #0.003

park_gdf = park
box_gdf  = gpd.GeoDataFrame()

# compute total bounds
minx, miny, maxx, maxy = park_gdf.total_bounds
minx -= buffer
maxx += buffer
miny -= buffer
maxy += buffer

# compute the side length of the square
side = max(maxx - minx, maxy - miny)
center_x = (minx + maxx) / 2
center_y = (miny + maxy) / 2

# create square box
square_minx = center_x - side / 2
square_maxx = center_x + side / 2
square_miny = center_y - side / 2
square_maxy = center_y + side / 2
square_geom = box(square_minx, square_miny, square_maxx, square_maxy)

# create GeoDataFrame for box
square_gdf = gpd.GeoDataFrame(
    {'PARKNAME': [f"{park_name}_box"]},
    geometry=[square_geom],
    crs=park_gdf.crs
)

# add box to collection
box_gdf = square_gdf


ax = tacony.plot(edgecolor='blue', alpha=0.5)
box_gdf[box_gdf['PARKNAME']==f"{park_name}_box"].plot(ax=ax, edgecolor='red', facecolor='none', linewidth=2)
box_gdf.plot()

plt.show()



### Filtering points in each park

## filter points in parks
illegal_dumping_in_park = illegal_dumping.sjoin(box_gdf, predicate='within')

park_points = illegal_dumping_in_park[illegal_dumping_in_park['PARKNAME'] == f"{park_name}_box"]
park_points.shape


# jittered times
random_offsets = np.random.uniform(0, 24 * 60 * 60, size=len(park_points))
illegal_dumping_in_park['jittered_time'] = illegal_dumping_in_park['requested_'] + pd.to_timedelta(random_offsets, unit='s')


# get coordinates from the geometry
coords = pd.DataFrame({
    'x': illegal_dumping_in_park.geometry.centroid.x,
    'y': illegal_dumping_in_park.geometry.centroid.y
})

# calculate time difference in days from the minimum date
offset_seasonal = 30 * 0
min_year        = pd.to_datetime(illegal_dumping_in_park['jittered_time'].min())
min_time        = pd.Timestamp(f"{min_year.year}-01-01")
time            = (illegal_dumping_in_park['jittered_time'] - min_time).dt.total_seconds() / (24 * 60 * 60)
seasonal        = (time + offset_seasonal) % 365
illegal_dumping_in_park['time_diff']     = time
illegal_dumping_in_park['seasonal_diff'] = seasonal

# Create locs_s DataFrame
locs_s = pd.DataFrame({
    'X': coords['x'].astype(float),
    'Y': coords['y'].astype(float),
    'T': time.astype(float),
    'A': seasonal.astype(float)
})



## Model fitting

## Cox-Hawkes model
# a_0 = baseline log-intensity
# alpha = magnitute of excitation
# bate = temporal trigger param
# sigmax_2 = spatial trigger param
model = bstpp.main.Hawkes_Model(
    locs_s,     # spatiotemporal points
    box_gdf[box_gdf['PARKNAME']==f"{park_name}_box"], # philly boundaries
    total_days, # Time frame
    True,       # use Cox as background
    a_0      = dist.Normal(0,2),
    alpha    = dist.Beta(2,2),
    beta     = dist.HalfNormal(1.0),
    sigmax_2 = dist.HalfNormal(0.25),
    offset_seasonal = offset_seasonal
)

#model.set_window(window = 12, spatial_window= 0.25)
#model.set_window(window = 30, spatial_window= 0.25)
#model.set_window(window = 100, spatial_window= 10.0)
model.set_window(window = 100, spatial_window= 0.1)

# less interations and bigger learning rate
model.run_svi(lr=0.01,num_steps=10000)




## Results

model.log_expected_likelihood(locs_s)

model.expected_AIC()


model.plot_prop_excitation()
plt.show()


model.plot_trigger_posterior(trace=False)
plt.show()


model.plot_trigger_time_decay()
plt.gcf()


model.plot_spatial(include_cov=False)
plt.show()


model.plot_temporal(start_date=min_time)
plt.show()


model.plot_seasonal()
plt.show()


# Create histogram
plt.hist(locs_s['A'], bins=12*5, edgecolor='black')

# Add labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of Column A')

# Show the plot
plt.show()




# Get all axes in the figure
fig = plt.gcf()
axs = fig.get_axes()

# Match limits from left subplot
xlim = axs[0].get_xlim()
ylim = axs[0].get_ylim()

# Try plotting on the second subplot (with events)
all_boxes_gdf.plot(ax=axs[1], edgecolor='blue', facecolor='none', linewidth=2)
philly_border.plot(ax=axs[1], edgecolor='red', facecolor='none', linewidth=2)
all_parks_gdf.plot(ax=axs[1], edgecolor='green', facecolor='none', linewidth=2)


# Set the same limits on the right subplot
axs[1].set_xlim(xlim)
axs[1].set_ylim(ylim)

# Save and show
plt.savefig('output/spatial_trigger_tcn.png', dpi=450, bbox_inches='tight')
plt.show()
