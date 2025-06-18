setwd("~/Dropbox/DDDI/Illegal-Dumping")
#install.packages("/Users/yuxin/Downloads/stelfi_1.0.1.tar.gz", repos = NULL, type = "source")
library(stelfi)
library(dplyr)
library(mapview)
library(sf)
library(INLA) 
library(lubridate)
library(ggplot2)
library(ggmap)
library(future)
plan(multisession)



# load 311 data
shp_24 <- read_sf("data/311/2024/public_cases_fc.shp")
shp_24_litter <- shp_24 %>%
  filter(service_na == "Illegal Dumping")
head(shp_24_litter)

sum(!is.na(shp_24_litter$media_url))
7322/22173


shp_23 <- read_sf("data/311/2023/public_cases_fc.shp")
shp_23_litter <- shp_23 %>%
  filter(service_na == "Illegal Dumping")
sum(!is.na(shp_23_litter$media_url))


head(shp_23_litter$media_url,10)
shp_litter_combined <- bind_rows(shp_24_litter, shp_23_litter)

# check time distribution
ggplot(shp_24_litter, aes(x = as.Date(requested_))) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(title = "Daily Event Distribution", x = "Date", y = "Count") +
  theme_minimal()

#########################################
# stelfi way
#########################################
#### 311 data ###
philly <- st_read("data/City_Limits/City_Limits.shp")
mapview(philly)

shp_24_litter_clean <- shp_24_litter %>%
  filter(!is.na(requested_), !st_is_empty(geometry))

locs_311 <- sf::st_coordinates(shp_24_litter_clean) |>
  as.data.frame()
names(locs_311) <- c("x", "y")

smesh_311 <- fmesher::fm_mesh_2d(loc = locs_311,
                                 max.edge = 0.01, cutoff = 0.005,
                                 crs = st_crs(philly))

plot(smesh_311)

fit_311 <- fit_lgcp(locs = locs_311, sf = philly, smesh = smesh_311,
                parameters = c(beta = 0, log_tau = log(1),
                               log_kappa = log(1)))
st_crs(philly)
get_coefs(fit_311) 

get_fields(fit_311, smesh_311) |>
  show_field(smesh = smesh_311, sf = philly, clip = TRUE) + ggplot2::theme_classic()
lgcp_plot <- show_lambda(fit_311, smesh = smesh_311, sf = philly, clip = TRUE) + ggplot2::theme_classic()
lgcp_plot

#### stelfi
philly <- st_read("data/City_Limits/City_Limits.shp")
mapview(philly)

shp_24_litter_clean <- shp_24_litter %>%
  filter(!is.na(requested_), !st_is_empty(geometry))

# time jitter
set.seed(999)
shp_24_litter_clean <- shp_24_litter_clean %>%
  mutate(
    jittered_time = as.POSIXct(requested_) + runif(n(), 0, 86400),
    date_day = as.Date(requested_)
  )

## locs stelfi
locs_s <- sf::st_coordinates(shp_24_litter_clean) %>%
  as.data.frame() %>%
  mutate(t = difftime(shp_24_litter_clean$jittered_time, 
                      min(shp_24_litter_clean$jittered_time), units = "days")) %>%
  mutate(t = as.numeric(t)) %>%
  rename(., c("x" = "X", "y" = "Y"))

locs_s <- locs_s[!duplicated(locs_s[, c('t')]), ]
locs_s <- locs_s[order(locs_s$t), ]
times <- locs_s$t
locs_s <- data.frame(x = locs_s$x, y = locs_s$y)


sf_use_s2(FALSE)

## Delauney triangluation of domain
#smesh_s <- fmesher::fm_mesh_2d(loc = locs_s[, 1:2], max.edge = 0.015, cutoff = 0.01)

smesh_s <- inla.mesh.2d(
  loc = locs_s,
  max.edge = c(0.03, 0.03),  # tune these based on your CRS
  cutoff = 0.01,
  offset = c(0.05, 0.1)  # ← increase to expand beyond boundary
)


plot(smesh_s)

st_crs(philly) <- NA

# set parameters
param <- list( mu = 10.604, alpha = 0.4, beta = 0.5, # ~1.4 days
               kappa = exp(4.0449), tau = exp(-6.1167), 
               xsigma = 0.005, ysigma= 0.005, rho = 0) # 500m

fit_stelfi_311 <- fit_stelfi(
  times = times,
  locs = locs_s, 
  sf = philly,
  smesh = smesh_s,
  parameters = param,
  GMRF = F,
)

get_coefs(fit_stelfi_311)

# mu: background intensity
# alpha : trigger intensity -- size of the spike caused by a past event
# beta: temporal decay -- how quickly excitation fades over time (1/beta days decay window)
# xsigma, ysigma, rho: spatial Guassian kernel param
# kappa, tau: latent field GP params -- spatial field smoothness
# coefs: coefficient on covariates

names(fit_stelfi_311)

## visualization
# set up param
mu     <- 2.298189641
alpha  <- 0.003926338
beta   <- 0.073410180
xsigma <- 0.019466328
ysigma <- 0.013326678

# define the target time
t0 <- 365

# extract event data
events <- data.frame(
  x = locs_s$x,
  y = locs_s$y,
  t = times
)


# build a prediction grid from my mesh
grid_pts <- smesh_s$loc[,1:2]
colnames(grid_pts) <- c("x", "y")

######################## philly grid
# 2. Create 150m x 150m-ish grid in degrees (approximate!)
grid_150m <- st_make_grid(
  philly,
  cellsize = c(0.00176, 0.00135),  # lon, lat
  square = TRUE
) %>%
  st_as_sf() %>%
  st_filter(philly)  # Clip to philly

# 3. Compute centroids (for evaluating lambda at each cell center)
grid_centroids <- st_centroid(grid_150m)

# 4. Convert to coords
grid_pts <- st_coordinates(grid_centroids)
colnames(grid_pts) <- c("x", "y")
########################

# define Guassian spatial kernel
gaussian_kernel <- function(x0, y0, x1, y1, xsigma, ysigma) {
  exp(-((x0 - x1)^2) / (2 * xsigma^2) - ((y0 - y1)^2) / (2 * ysigma^2))
}

# calculate lamdas (s, t0)
lambda_hat <- apply(grid_pts, 1, function(pt) {
  xg <- pt[1]
  yg <- pt[2]
  lam <- mu +
    sum(
      ifelse(events$t < t0, 
             alpha * exp(-beta * (t0 - events$t)) *
               gaussian_kernel(xg, yg, events$x, events$y, xsigma, ysigma),
             0
      )
    )
  lam
})

# df for plotting
lambda_df <- as.data.frame(cbind(grid_pts, lambda_hat))

# plot
## base map
register_stadiamaps(key = "314edc79-5661-4f96-8818-ce1de3666a6b")
philly <- st_transform(philly, 4326)
bbox <- st_bbox(philly)
bbox

philly_map <- get_stadiamap(
  bbox = c(left = -75.28031, bottom = 39.86747,
           right = -74.95575, top = 40.13793 ),
  zoom = 12,
  maptype = "alidade_smooth"  # or "stamen-terrain", "stamen-watercolor"
)

ggmap(philly_map) +
  geom_tile(data = lambda_df, aes(x = x, y = y, fill = lambda_hat), alpha = 0.8) +
  scale_fill_viridis_c(name = "λ(s, t=10)") +
  labs(
    title = "Estimated Intensity with Stadia Map Base",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

#########################################
# Manually fit the model 
#########################################
# Filter for valid geometry and event time
df <- shp_24_litter %>%
  filter(!is.na(requested_), !st_is_empty(geometry)) %>%
  mutate(
    time = as.numeric(as.Date(requested_) - min(as.Date(requested_))),  # Days since start
    lon = st_coordinates(geometry)[,1],
    lat = st_coordinates(geometry)[,2]
  )


# Discretize time and space
df$time_bin <- as.integer(cut(df$time, breaks = 20))
df$space_bin <- as.integer(cut(df$lon + df$lat, breaks = 30))  # crude 1D spatial bin

# Build model formula: additive time + space GP with log-link
df$y <- 1  # every row is a "point"
formula <- y ~ f(time_bin, model = "rw2") + f(space_bin, model = "rw2")

# Fit background only (log-Gaussian Cox)
fit_inla <- inla(formula, family = "poisson", data = df,
                 control.predictor = list(compute = TRUE),
                 control.compute = list(dic = TRUE))

# Predict background intensity
df$log_mu_hat <- fit_inla$summary.linear.predictor$mean
df$mu_hat <- exp(df$log_mu_hat)

## fit hawkes kernel on residuals
# Calculate triggering component manually
df <- df %>% arrange(time)
alpha <- 0.5
beta <- 7
sigma_x <- 0.0005
sigma_y <- 0.0005

df$triggered <- 0
for (i in 2:nrow(df)) {
  past <- df[1:(i - 1), ]
  delta_t <- df$time[i] - past$time
  dx <- df$lon[i] - past$lon
  dy <- df$lat[i] - past$lat
  kern <- alpha * beta * exp(-beta * delta_t) *
    (1 / (2 * pi * sigma_x * sigma_y)) *
    exp(- (dx^2 / (2 * sigma_x^2) + dy^2 / (2 * sigma_y^2)))
  df$triggered[i] <- sum(kern)
  
  if (i %% 500 == 0) {
    cat("Finished row", i, "of", nrow(df), "\n")
  }
}

df$lambda_hat <- df$mu_hat + df$triggered


# fial fit
loglik <- sum(log(df$lambda_hat)) - sum(df$lambda_hat)

ggplot(df, aes(x = lon, y = lat, color = log(lambda_hat))) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "C") +
  theme_minimal() +
  labs(title = "Estimated Background Intensity (LGCP)",
       subtitle = expression(mu(t, x, y) == exp(f[s](x,y) + f[t](t))),
       x = "Longitude", y = "Latitude", color = "lambda_hat")

##### smooth heatmap
# load philly shapefile in 2023
penn_24 <- read_sf("data/tl_2024_42_bg/tl_2024_42_bg.shp")
philly_24 <- penn_24 %>%
  filter(COUNTYFP == "101")
mapview(philly_24)

#Spatial join: assign each point to a CBG
df <- st_transform(df, st_crs(philly_24))
df_joined <- st_join(df, philly_24, join = st_within)

#Aggregate lambda_hat to each CBG (mean or sum depending on interpretation)
cbg_agg <- df_joined %>%
  st_drop_geometry() %>%
  group_by(GEOID) %>%
  summarise(
    lambda_avg = mean(lambda_hat, na.rm = TRUE),
    n_events = n()
  ) %>%
  left_join(philly_24, by = "GEOID") %>%
  st_as_sf()

#Plot
ggplot(cbg_agg) +
  geom_sf(aes(fill = log(lambda_avg)), color = NA) +
  scale_fill_viridis_c(name = "Avg logged λ in CBG", option = "C") +
  theme_minimal() +
  labs(
    title = "Heatmap of Illegal Dumping Risk by Census Block Group",
    subtitle = "Based on LGCP + Hawkes model fit over 1 year",
    fill = "λ (risk)"
  )

# Get bounding box
bbox <- st_bbox(philly_24)

# Get map image (Google Maps requires API key)
ggmap::register_stadiamaps(key = "314edc79-5661-4f96-8818-ce1de3666a6b")
map <- get_stadiamap(
  bbox = c(left = bbox["xmin"], bottom = bbox["ymin"],
           right = bbox["xmax"], top = bbox["ymax"]),
  zoom = 10,
  maptype = "stamen_toner_lite"
)

ggmap(map) +
  geom_sf(data = philly_24, inherit.aes = FALSE, fill = NA, color = "black")


bbox <- c(left = -75.3, bottom = 39.85, right = -74.9, top = 40.1)

# Download Stamen basemap (this should work)
map <- get_stadiamap(
  bbox = bbox,
  zoom = 10,
  maptype = "stamen_toner_lite"
)

ggmap(map)


#########################################
# Functions
#########################################

# model fitting function
fit_cox_hawkes_model <- function(df_raw) {
  # Preprocess: drop missing and empty geometry
  df <- df_raw %>%
    filter(!is.na(requested_), !st_is_empty(geometry)) %>%
    st_transform(4326) %>%  # Use WGS 84 (or switch to 26918 for better grid alignment)
    mutate(
      time = as.numeric(as.Date(requested_) - min(as.Date(requested_))),
      lon = st_coordinates(geometry)[, 1],
      lat = st_coordinates(geometry)[, 2]
    )
  
  # -----------------------------
  # TIME BINNING (~weekly bins)
  # -----------------------------
  time_breaks <- seq(min(df$time), max(df$time), length.out = 53)  # 52 weeks
  df$time_bin <- as.integer(cut(df$time, breaks = time_breaks, include.lowest = TRUE))
  
  # -----------------------------
  # SPACE BINNING using 2D grid
  # -----------------------------
  grid <- st_make_grid(df, cellsize = 0.003, square = TRUE)  # ~300m grid
  grid_sf <- st_sf(grid_id = seq_along(grid), geometry = grid)
  
  # Join each point to its grid cell
  df <- st_join(df, grid_sf, join = st_within)
  df$space_bin <- df$grid_id
  
  # -----------------------------
  # LGCP Background Model (INLA)
  # -----------------------------
  df$y <- 1
  formula <- y ~ f(time_bin, model = "rw2") + f(space_bin, model = "rw2")
#  formula <- y ~ f(space_bin, model = "rw2")
  
  fit_inla <- inla(
    formula,
    family = "poisson",
    data = df,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
  )
  
  df$log_mu_hat <- fit_inla$summary.linear.predictor$mean
  df$mu_hat <- exp(df$log_mu_hat)
  
  # -----------------------------
  # Hawkes Triggering Component
  # -----------------------------
  df <- df %>% arrange(time)
  alpha <- 1
  beta <- 0.5
  sigma_x <- sigma_y <- 0.001
  
  df$triggered <- 0
  for (i in 2:nrow(df)) {
    past <- df[1:(i - 1), ]
    delta_t <- df$time[i] - past$time
    dx <- df$lon[i] - past$lon
    dy <- df$lat[i] - past$lat
    
    kern <- alpha * beta * exp(-beta * delta_t) *
      (1 / (2 * pi * sigma_x * sigma_y)) *
      exp(- (dx^2 / (2 * sigma_x^2) + dy^2 / (2 * sigma_y^2)))
    
    df$triggered[i] <- sum(kern)
    
    if (i %% 1000 == 0) cat("Finished", i, "of", nrow(df), "observations\n")
  }
  
  df$lambda_hat <- df$mu_hat + df$triggered
  
  # -----------------------------
  # Return everything needed
  # -----------------------------
  return(list(
    df = df,
    log_mu_hat = log_mu_hat,
    fit_inla = fit_inla,
    alpha = alpha,
    beta = beta,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    min_time = min(as.Date(df_raw$requested_)),
    time_breaks = time_breaks,
    spatial_grid = grid_sf  # keep spatial bin structure for prediction
  ))
}



## prediction function
predict_lambda_hat <- function(new_points_df, model) {
  cat("DEBUG: input has", nrow(new_points_df), "rows\n")
  cat("Starting prediction for", nrow(new_points_df), "points\n")
  
  # Convert to sf and match CRS
  new_points_sf <- new_points_df %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform(crs = st_crs(model$spatial_grid))
  
  # Assign time bin using model's breaks
  new_points_sf$time <- as.numeric(as.Date(new_points_sf$date) - model$min_time)
  new_points_sf$time_bin <- as.integer(cut(new_points_sf$time, breaks = model$time_breaks, include.lowest = TRUE))
  
  # Assign space bin using spatial join
  new_points_sf <- st_join(new_points_sf, model$spatial_grid, join = st_within)
  new_points_sf$space_bin <- new_points_sf$grid_id
  
  # Get mu_hat lookup
  mu_lookup <- model$df %>%
    st_drop_geometry() %>% 
    group_by(time_bin, space_bin) %>%
    summarise(mu_hat = mean(mu_hat, na.rm = TRUE), .groups = "drop")
  
#  mu_lookup <- model$df %>%
#    st_drop_geometry() %>% 
#    group_by(space_bin) %>%
#    summarise(mu_hat = mean(mu_hat, na.rm = TRUE), .groups = "drop")
  
  new_points_sf <- new_points_sf %>%
    left_join(mu_lookup, by = c("time_bin", "space_bin")) %>%
    mutate(mu_hat = ifelse(is.na(mu_hat), mean(model$df$mu_hat, na.rm = TRUE), mu_hat))
  
#  new_points_sf <- new_points_sf %>%
#    left_join(mu_lookup, by = c("space_bin")) %>%
#    mutate(mu_hat = ifelse(is.na(mu_hat), mean(model$df$mu_hat, na.rm = TRUE), mu_hat))
  
  # Triggered component
  df_past <- model$df
  triggered_vec <- numeric(nrow(new_points_sf))
  
  for (i in seq_len(nrow(new_points_sf))) {
    point <- new_points_sf[i, ]
    valid_past <- df_past %>% filter(time < point$time)
    
    if (nrow(valid_past) > 0) {
      delta_t <- point$time - valid_past$time
      dx <- point$lon - valid_past$lon
      dy <- point$lat - valid_past$lat
      
      kern <- model$alpha * model$beta * exp(-model$beta * delta_t) *
        (1 / (2 * pi * model$sigma_x * model$sigma_y)) *
        exp(- (dx^2 / (2 * model$sigma_x^2) + dy^2 / (2 * model$sigma_y^2)))
      
      triggered_vec[i] <- sum(kern)
    }
    
    if (i %% 500 == 0) {
      cat("Finished", i, "of", nrow(new_points_sf), "points\n")
    }
  }
  
#  new_points_sf$lambda_hat <- new_points_sf$mu_hat + triggered_vec
  new_points_sf$lambda_hat <- new_points_sf$mu_hat
  return(new_points_sf)
}


# fit model
model <- fit_cox_hawkes_model(shp_24_litter)
   
## gridding philadelphia: points
philly <- st_read("data/City_Limits/City_Limits.shp")

# Create grid with spacing ~0.002 degrees (~200 meters, rough)
grid_points <- st_make_grid(philly, cellsize = 0.003, what = "centers") %>%
  st_as_sf(crs = 4326) %>%
  st_filter(philly)

# Extract lon/lat
grid_df <- grid_points %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

# add time
grid_df <- grid_df %>%
  mutate(date = as.Date("2025-01-01"))

### predict model
predicted <- predict_lambda_hat(grid_df, model)

# visualize
predicted_sf <- st_as_sf(predicted, coords = c("lon", "lat"), crs = 4269, remove = FALSE)
philly <- st_transform(philly, crs = 4269)
predicted_sf <- st_filter(predicted_sf, philly)

predicted_df <- st_coordinates(predicted_sf) %>%
  as.data.frame() %>%
  rename(x = X, y = Y) %>%
  bind_cols(predicted_sf)

ggplot(predicted_df, aes(x = lon, y = lat, color = mu_hat)) +
  geom_point(size = 2, alpha = 0.8, shape = 15) +
  scale_color_viridis_c(
    option = "A",
    limits = c(0.96, 1),
    breaks = seq(0.96, 1, by = 0.01),
#    name = expression(log[1+p](lambda[hat]))
  ) +
  theme_minimal() +
  labs(
    title = "Predicted Litter Dumping Index (2024-12-31)",
#    subtitle = expression(mu(t, x, y) == exp(f[s](x,y) + f[t](t))),
    x = "Longitude", y = "Latitude"
  )

summary(model$df$mu_hat)
summary(shp_24_litter$requested_)
head(predicted_df)

table(shp_24_litter)
range(shp_24_litter$requested_)
summary(shp_24_litter)

