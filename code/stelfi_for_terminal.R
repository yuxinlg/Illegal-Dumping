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


# load 311 data
shp_24 <- read_sf("data/311/2024/public_cases_fc.shp")
shp_24_litter <- shp_24 %>%
  filter(service_na == "Illegal Dumping")
head(shp_24_litter)
#### stelfi
philly <- st_read("data/City_Limits/City_Limits.shp")

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

smesh_s <- inla.mesh.2d(
  loc = locs_s,
  max.edge = c(0.05, 0.05),  # tune these based on your CRS
  cutoff = 0.02,
  offset = c(0.05, 0.1)  # ← increase to expand beyond boundary
)


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
# Save result
saveRDS(fit_stelfi_311, file = "fit_stelfi_311.rds")

coefs <- get_coefs(fit_stelfi_311)
write.csv(coefs, file = "fit_stelfi_311_coefs.csv", row.names = FALSE)

cat("✔️ Model finished fitting at", Sys.time(), "\n")
writeLines("DONE", "fit_status.txt")