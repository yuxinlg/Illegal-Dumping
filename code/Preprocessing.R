setwd("~/Dropbox/DDDI/Illegal-Dumping")
library(dplyr)
library(tidyr)
library(stringr)
library(sf)
library(mapview)
library(tidycensus)

################################### Geom ref ###################################
# Pennsylvania state FIPS: 42
# Philadelphia county FIPS: 101
# GEOID: 42101

############### import 311 litter dumping survey data: 2017-2019 ###############
data <- read.csv("data/Litter_Index_Surveys/Litter_Index_Surveys.csv")
head(data,20)
## get litter types
# split the strings by comma
split_litter <- str_split(data$LITTER_TYPES, ",")

# unlist to flattern the list into a vector
flat_litter <- unlist(split_litter)

# get the uniqur litter types
unique_litter_types <- unique(flat_litter)
unique_litter_types

################################ import 311 report #############################
# link here: https://opendataphilly.org/datasets/311-service-and-information-requests/

# set base directory
base_dir <- "data/311"

# rename map for shp
rename_map <- c(
  "request_id" = "service_re" ,
  "status_notes" = "status_not",
  "service_name" = "service_na",
  "service_code" = "service_co",
  "agency_responsible" = "agency_res",
  "service_notice" = "service_no",
  "requested_datetime" = "requested_",
  "updated_datetime" = "updated_da",
  "expected_datetime" = "expected_d",
  "closed_datetime" = "closed_dat"
)

# loop through years to get shp files
for (year in 2014:2024) {
  year_short <- substr(as.character(year), 3, 4)
  folder_path <- file.path(base_dir, as.character(year))
  
  # find shp files
  shp_file <- list.files(folder_path, pattern = "\\.shp$", full.names = TRUE)[1]
  
  # load files and rename
  shp <- st_read(shp_file, quiet = TRUE)
  
  # rename shp cols
  for (new_name in names(rename_map)) {
    old_name <- rename_map[[new_name]]
    if (old_name %in% names(shp)) {
      names(shp)[names(shp) == old_name] <- new_name
    }
  }
  
  assign(paste0("shp_", year_short), shp)
  rm(shp)
}

# loop for filtering "Illegal Dumping" service
for (year in 2014:2024) {
  year_short <- substr(as.character(year), 3, 4)
  shp_name <- paste0("shp_", year_short)
  shp_var <- get(shp_name)
  
  # Filter for service_name == "Illegal Dumping"
  shp_illegal <- shp_var[shp_var$service_name == "Illegal Dumping", ]
  
  # Assign new variable
  assign(paste0("shp_illegal_", year_short), shp_illegal)
  
  # Clean up
  rm(shp_illegal)
  rm(shp_var)
  
  # Remove original shp_XX if year is between 2014 and 2021
  if (as.numeric(year_short) >= 14 && as.numeric(year_short) <= 24) {
    if (exists(shp_name)) rm(list = shp_name)
  }
}

# loop for count NA in `media_url`
na_counts <- data.frame(name = character(), 
                        count_media = integer(), 
                        observation = integer(),
                        pct = character(),
                        stringsAsFactors = FALSE)

for (year in 2014:2024) {
  year_short <- substr(as.character(year),3,4)
  shp_name <- paste0("shp_illegal_", year_short)
  shp_var <- get(shp_name)
  
  # count NA in media
  count_na <- sum(is.na(shp_var$media_url))
  
  # get the number of media
  count_media <- nrow(shp_var) - count_na
  
  # get number of observations
  obs = nrow(shp_var)
  
  # get percentage
  pct = paste0(round((count_media/obs),2)*100, "%")
  
  # append result
  na_counts <- rbind(na_counts, data.frame(name = shp_name, 
                                           count_media = count_media,
                                           observation = obs,
                                           pct = pct))
}

na_counts

shp_illegal_16 %>%
  filter(!is.na(media_url))




############################# Get Philly ACS data ##############################
# manually download ACS shapefile
# link here: https://www2.census.gov/geo/tiger/TIGER2017/BG/tl_2017_42_bg.zip
penn_shp <- read_sf("data/tl_2017_42_bg/tl_2017_42_bg.shp")
philly_shp <- penn_shp %>%
  filter(COUNTYFP == "101")
mapview(philly_shp)

# get ACS survey data
my_vars <- c(
  #population
  total_pop = "B01003_001",
  #gender
  male_total = "B01001_002",
  female_total = "B01001_026",
  # Male
  m_under_5       = "B01001_003",
  m_5_9           = "B01001_004",
  m_10_14         = "B01001_005",
  m_15_17         = "B01001_006",
  m_18_19         = "B01001_007",
  m_20            = "B01001_008",
  m_21            = "B01001_009",
  m_22_24         = "B01001_010",
  m_25_29         = "B01001_011",
  m_30_34         = "B01001_012",
  m_35_39         = "B01001_013",
  m_40_44         = "B01001_014",
  m_45_49         = "B01001_015",
  m_50_54         = "B01001_016",
  m_55_59         = "B01001_017",
  m_60_61         = "B01001_018",
  m_62_64         = "B01001_019",
  m_65_66         = "B01001_020",
  m_67_69         = "B01001_021",
  m_70_74         = "B01001_022",
  m_75_79         = "B01001_023",
  m_80_84         = "B01001_024",
  m_85_over       = "B01001_025",
  # Female
  f_under_5       = "B01001_027",
  f_5_9           = "B01001_028",
  f_10_14         = "B01001_029",
  f_15_17         = "B01001_030",
  f_18_19         = "B01001_031",
  f_20            = "B01001_032",
  f_21            = "B01001_033",
  f_22_24         = "B01001_034",
  f_25_29         = "B01001_035",
  f_30_34         = "B01001_036",
  f_35_39         = "B01001_037",
  f_40_44         = "B01001_038",
  f_45_49         = "B01001_039",
  f_50_54         = "B01001_040",
  f_55_59         = "B01001_041",
  f_60_61         = "B01001_042",
  f_62_64         = "B01001_043",
  f_65_66         = "B01001_044",
  f_67_69         = "B01001_045",
  f_70_74         = "B01001_046",
  f_75_79         = "B01001_047",
  f_80_84         = "B01001_048",
  f_85_over       = "B01001_049",
  #commute
  commute_total = "B08301_001",
  commute_car = "B08301_002",
  commute_public = "B08301_010",
  commute_walked = "B08301_018",
  commute_bike = "B08301_019",
  commute_home = "B08301_020",
  #tenure
  tenure_total = "B25003_001", #Total occupied housing units
  tenure_owner = "B25003_002", #	Owner-occupied units
  tenure_renter = "B25003_003",
  #housing & vacancy
  owner_1_1_5 = "B25014_005",
  owner_1_5_plus = "B25014_006",
  owner_2_plus = "B25014_007",
  renter_1_1_5 = "B25014_011",
  renter_1_5_plus = "B25014_012",
  renter_2_plus = "B25014_013",
  occ_total = "B25002_001", 
  occ_occupied = "B25002_002",
  occ_vacant = "B25002_003",
  median_built = "B25035_001",
  #age
  median_age = "B01002_001",
  #race
  race_white = "B02001_002",
  race_black = "B02001_003",
  race_native = "B02001_004",
  race_asian = "B02001_005",
  race_islander = "B02001_006",
  race_other = "B02001_007",
  race_two = "B02001_008",
  #ethnic
  hispanic = "B03003_003",
  #income
  median_income = "B19013_001",
  #education 
  edu_total = "B15003_001",
  edu_no_school = "B15003_002",
  edu_4th = "B15003_003",
  edu_56th = "B15003_004",
  edu_78th = "B15003_005",
  edu_9th = "B15003_006",
  edu_10th = "B15003_007", 
  edu_11th = "B15003_008",
  edu_12th = "B15003_009", 
  edu_highschool = "B15003_010",
  edu_college_less1 = "B15003_011",
  edu_somecollege = "B15003_012", 
  edu_associate = "B15003_013",
  edu_bachelor = "B15003_014",
  edu_master = "B15003_015",
  edu_professional = "B15003_016",
  edu_doctorate = "B15003_017",
  #employment
  employ_total = "B23025_001",   # total 16+
  employ_labor_force = "B23025_002",
  employ_unemployed = "B23025_005",
  #Limited English Proficiency: speaks English < "very well"
  sp_lep1 = "B16004_007", # age group 5-17
  sp_lep2 = "B16004_008", # age group 18-64
  sp_lep3 = "B16004_009", # age group 65+
  ie_lep1 = "B16004_015",
  ie_lep2 = "B16004_016",
  ie_lep3 = "B16004_017",
  ap_lep1 = "B16004_023",
  ap_lep2 = "B16004_024",
  ap_lep3 = "B16004_025",
  other_lep1 = "B16004_031",
  other_lep2 = "B16004_032",
  other_lep3 = "B16004_033",
  #kitchen
  kitchen_total = "B25051_001",
  kitchen_complate = "B25051_002",
  kitchen_lack = "B25051_003",
  #plumbing
  plumbing_total = "B25049_001",
  plumbing_complete = "B25049_002",
  plumbing_lack = "B25049_003",
  #poverty
  poverty_total = "B17001_001",
  poverty_below = "B17001_002",
  #HH structure
  hh_total = "B11016_001E",
  hh_fam = "B11016_002",
  hh_fam_5 = "B11016_006",
  hh_fam_6 = "B11016_007",
  hh_fam_7 = "B11016_008",
  hh_non_fam = "B11016_009"
)

acs_17 <- get_acs(
  geography = "block group",
  variables = my_vars,
  state = "PA",
  county = "Philadelphia",
  year = 2017,
  geometry = FALSE
)

# reshape philly_17: long to wide
philly_17 <- acs_17 %>%
  select(GEOID, variable, estimate) %>%
  pivot_wider(names_from = variable, values_from = estimate)

# variable aggregation
philly_17 <- philly_17 %>%
  mutate(
    # age groups: 5 level
    age_under_18 = rowSums(select(., m_under_5, m_5_9, m_10_14, m_15_17,
                                  f_under_5, f_5_9, f_10_14, f_15_17), na.rm = TRUE),
    age_18_24 = rowSums(select(., m_18_19, m_20, m_21, m_22_24,
                               f_18_19, f_20, f_21, f_22_24), na.rm = TRUE),
    age_25_44 = rowSums(select(., m_25_29, m_30_34, m_35_39, m_40_44,
                               f_25_29, f_30_34, f_35_39, f_40_44), na.rm = TRUE),
    age_45_64 = rowSums(select(., m_45_49, m_50_54, m_55_59, m_60_61, m_62_64,
                               f_45_49, f_50_54, f_55_59, f_60_61, f_62_64), na.rm = TRUE),
    age_older_65 = rowSums(select(., m_65_66, m_67_69, m_70_74, m_75_79, m_80_84, m_85_over,
                                  f_65_66, f_67_69, f_70_74, f_75_79, f_80_84, f_85_over), na.rm = TRUE),
    
    # education: 5 levels
    edu_none = rowSums(select(., edu_no_school, edu_4th, edu_56th, edu_78th, edu_9th), na.rm = TRUE),
    edu_somehs = rowSums(select(., edu_10th, edu_11th, edu_12th), na.rm = TRUE),
    edu_hsgrad = edu_highschool,
    edu_somecollege = rowSums(select(., edu_college_less1, edu_somecollege, edu_associate), na.rm = TRUE),
    edu_bachplus = rowSums(select(., edu_bachelor, edu_master, edu_professional, edu_doctorate), na.rm = TRUE),
    
    # race: recode other
    race_other_comb = rowSums(select(., race_native, race_islander, race_other,race_two), na.rm = TRUE),
    
    # limited english proficiency
    lep_total = rowSums(select(., sp_lep1, sp_lep2, sp_lep3,ie_lep1, ie_lep2, ie_lep3,
                               ap_lep1, ap_lep2, ap_lep3,other_lep1, other_lep2, other_lep3), na.rm = TRUE),
    
    # overcrowded
    owner_overcrowded = rowSums(select(., owner_1_5_plus, owner_2_plus), na.rm = TRUE),
    renter_overcrowded = rowSums(select(., renter_1_5_plus, renter_2_plus), na.rm = TRUE),
    total_overcrowded = owner_overcrowded + renter_overcrowded,
    
    #
  )

head(philly_17)

########################### Get Philly land use data ###########################
land_use_shp <- read_sf("data/Land_Use/Land_Use.shp")
landuse_13 <- land_use_shp %>%
  filter(C_DIG2 == 13)
mapview(landuse_23)

head(land_use_shp$YEAR)


