
#### Script for studying magnetic orientation of raptors, here REKI 
#### Script written by Filibert Heim, filibert.heim@posteo.de, following instructions from Steffen Oppel and hints from Will Schneider 
#### This script is a trail to extend the previous script to improve visualisation and presicion of output 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load required packages 
library(oce)
library(tmap)
library(basemaps)
library(sf)
library(terra)
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
library(tidyverse)
library(lubridate)
select <- dplyr::select
filer <- dplyr::filter 
rename <- dplyr::rename

# read in data 
data <- read.csv(file = 'data/REKI_sample_locations10.csv') 

# check data structure 
head(data)
str(data)

# improve data structure 
data <- data %>% 
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), 
         migration = as.factor(if_else(MIGRATION == 'migrating', 'migratory', 'stationary')),
         long = as.numeric(sub("\\|.*", "", geometry)), 
         lat = as.numeric(sub(".*\\|", "", geometry))) %>% # also as.character
  select(-MIGRATION, -geometry, -inclination, -declination, -intensity)

# create column which tells whether the current migration was the first autumn migration 
data <- data %>%
  group_by(id) %>%
  mutate(
    first_migration = cumsum(migration == 'migratory' & lag(migration, default = 'stationary') != 'migratory'),
    first_migration = ifelse(migration == 'migratory', first_migration, NA)) %>%
  ungroup() 

data <- data %>% 
  group_by(id) %>% 
  filter(migration == 'migratory') %>% 
  mutate(
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)), 1, first_migration), 
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)) + 1 & timestamp < as.POSIXct(paste0(year(min(timestamp, na.rm = TRUE)) + 1, '-07-15 00:00:00'), format = "%Y-%m-%d %H:%M:%S"), 2, first_migration)) %>%
  filter(first_migration %in% c(1,2)) %>% select(id, timestamp, first_migration) %>% 
  right_join(data %>% select(-first_migration), by = join_by(id, timestamp))

data <- data %>% mutate(first_migration = factor(case_when(first_migration == 1 ~ 'autumn', 
                                             first_migration == 2 ~ 'spring', 
                                            TRUE ~ as.character(first_migration))))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Create a raster for the migration route of REKI's  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create binding box 
bbox <- terra::ext(range(data$long), range(data$lat))

# create raster with geographic CRS for Europe where grid cell size varies with location
rast <- rast(ext = bbox, resolution = 0.5, crs = "EPSG:4326") # needs to be 'EPSG:4326' because that is the way in which data are stored

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Get magnetic field values for first spring migration in raster  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# filter spring migration data to 1 location per day
days_spring <- data %>% filter(first_migration == 'spring') %>%
  mutate(date = date(timestamp)) %>% 
  group_by(id, date) %>% 
  slice(1) %>% 
  select(id, date) # possibly add timestamp if needed for magneticField() extraction 


# create df to extract magnetic field values for spring migration for every bird 

birds <- unique(data$id) # get ids from all birds 

df_rast <- as.data.frame(xyFromCell(rast, 1:ncell(rast))) # get the coordinates from raster
df_rast <- df_rast %>% rename(long = x, lat = y) # rename raster

df_rast <- data.frame(long = rep(df_rast$long, times = length(birds)), 
                      lat = rep(df_rast$lat, times = length(birds)), 
                      id = rep(birds, each = nrow(df_rast)),
                      inclination = NA, 
                      declination = NA, 
                      intensity = NA)

df_rast <- df_rast %>% full_join(days_spring, by = join_by(id), relationship = 'many-to-many') # add the mean_spring migration dates 

# get magnetic field values for mean date of each birds first spring migration
df_rast$inclination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$date)$inclination
df_rast$declination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$date)$declination
df_rast$intensity <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$date)$intensity


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Get magnetic field values for first autumn migration for all locations of each bird --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# filter for autumn migration
data_autumn <- data %>% 
  filter(first_migration == 'autumn') %>%
  select(-migration, -first_migration, -distance_to_next, -speed)

# get magnetic field values for all locations in first autumn migration for the exact timestamp
data_autumn$declination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$declination
data_autumn$inclination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$inclination
data_autumn$intensity <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$intensity

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Get information of grid cells where autumn migration values lied within spring migration  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# round values which represents birds ability to sense magnetic field values, ATTENTION: only autumn migration magnetic field values are rounded, not spring migration
# information on ability to sense magnetic field values used form Schneider et al. (2023) https://www.nature.com/articles/s42003-023-04530-w
# used the narrower/more optimistic uncertainties, which in the end lead to narrower bands of locations the birds already experienced - leads to more conservative outcomes 
data_autumn$declination_max <- data_autumn$declination + 0.5 # declination sensitivity 0.5
data_autumn$declination_min <- data_autumn$declination - 0.5 
data_autumn$inclination_max <- data_autumn$inclination + 0.5 # inclination sensitivity 0.5
data_autumn$inclination_min <- data_autumn$inclination - 0.5 
data_autumn$intensity_max <- data_autumn$intensity + 200 # intensity sensitivity 200
data_autumn$intensity_min <- data_autumn$intensity - 200


# PAIR EVERY GRID CELL IN SPRING PER DAY WITH EVERY LOCATION IN AUTUMN TO SELECT CELLS WITH 'experienced' MAGNETIC FIELD
birds <- unique(data$id) # get ids from all birds 
overlap_list <- list() # create a list to store all the output

# start the loop over every bird individual (birds), otherwise my R collapses
for(i in 1:length(birds)){
  # first filter for one individual and conduct calculation
  overlap_id <- df_rast %>% filter(id == birds[i]) %>%
    full_join(data_autumn %>% select(-long, -lat), by = 'id', suffix = c('_spring', '_autumn'), relationship = "many-to-many") %>% 
    mutate(match_inclination = inclination_spring >= inclination_min & inclination_spring <= inclination_max, # create a column that checks whether inclination from spring lies between max and min of autumn migration
           match_declination = declination_spring >= declination_min & declination_spring <= declination_max,
           match_intensity = intensity_spring >= intensity_min & intensity_spring <= intensity_max, 
           match_all = match_inclination & match_declination & match_intensity) %>% # creates a TRUE if all other columns above are TRUE
    #select(-match_inclination, -match_declination, -match_intensity) %>%  # this is needed when only the overview map should be extracted
    filter(match_all == TRUE | match_inclination == TRUE | match_declination == TRUE | match_intensity == TRUE) # remove all grid cells that are not needed
  # save result of one bird to the output file
  overlap_list[[i]] <- overlap_id
}

# summarise output list in one df
overlap <- bind_rows(overlap_list)
overlap_relevant <- overlap %>% # filter for relevant grid cells where at least one 
  group_by(id, lat, long, date) %>% 
  slice(1)

# remove all data which is duplicate because we did the above check for overlap for daily data 
overlap <- overlap %>% select(-date) %>% group_by(id, lat, long, match_inclination, match_declination, match_intensity, match_all) %>% slice(1)

dim(overlap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Plot data --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a plot 
ggplot() + 
  geom_point(data = overlap %>% filter(match_all == TRUE), 
             mapping = aes(x = long, y = lat, color = 'Experienced Magnetic Field in Autumn'), size = 1.5) +
  geom_path(data = data %>% filter(first_migration == 'spring'), 
            mapping = aes(x = long, y = lat, color = 'Spring Migration')) +
  coord_equal() + 
  facet_wrap(~id) +
  labs(title = "Experienced Magnetic Field During Autumn Migration and Spring Migration route for ten Red Kites",
    x = "Longitude",
    y = "Latitude") +
  scale_color_manual(name = 'Legend', values = c('Experienced Magnetic Field in Autumn' = 'grey80', 'Spring Migration' = 'red')) +
  theme_bw() + 
  theme(legend.position = c(0.75,0.12), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/REKI_magnetic_field_autumn_spring.png', height = 9, width = 12)



# create a data frame in long format for easier plotting
overlap_plot <- overlap %>% 
  mutate(mag_field = case_when( # be careful with a sensitive order of recoding due to potential to overwrite cases from earlier assigned labels 
    match_all == TRUE ~ 'Overlap Inclination, Declination & Intensity',
    match_inclination == TRUE & match_declination == TRUE ~ 'Overlap Inclination & Declination', 
    match_inclination == TRUE & match_intensity == TRUE ~ 'Overlap Inclination & Intensity', 
    match_declination == TRUE & match_intensity == TRUE ~ 'Overlap Declination & Intensity', 
    match_inclination == TRUE ~ 'Inclination', 
    match_declination == TRUE ~ 'Declination', 
    match_intensity == TRUE ~ 'Intensity', 
    TRUE ~ NA_character_))
  
# create plot where all combinations of matching magnetic field values are colored distinctly 
ggplot() + 
  geom_point(data = overlap_plot , 
             mapping = aes(x = long, y = lat, color = mag_field), size = 2.5, alpha = .8) +
  geom_path(data = data %>% filter(first_migration == 'spring'), 
            mapping = aes(x = long, y = lat, color = 'Spring Migration'), linewidth = 2) +
  coord_equal() + 
  facet_wrap(~id) +
  labs(title = "Experienced Magnetic Field During Autumn Migration and Spring Migration route for ten Red Kites",
       x = "Longitude", y = "Latitude", color = 'Legend') +
  scale_color_manual(values = c(
    "Inclination" = "#1f77b4",     
    "Declination" = "#ff7f0e",
    "Intensity" = "#2ca02c",  
    "Overlap Inclination & Declination" = "#9467bd", 
    "Overlap Inclination & Intensity" = "#8c564b",   
    "Overlap Declination & Intensity" = "#e377c2",   
    "Overlap Inclination, Declination & Intensity" = "darkgrey", 
    "Spring Migration" = "#d62728")) +
  theme_bw() + 
  theme(legend.position = c(0.75,0.12), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/REKI_magnetic_field_autumn_spring.png', height = 9, width = 12)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Plot data on a map background --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## convert data to spatial objects - only needed if tmap is used
str(overlap)
magField_sf<-overlap %>%
  st_as_sf(coords = c("long", "lat"), crs=4326) %>%
  st_transform(3857)


springMig_sf<-data %>% filter(first_migration == 'spring') %>%
  st_as_sf( coords = c("long", "lat"), crs=4326) %>%
  group_by(id) %>%
  #arrange(timestamp) %>%
  summarize() %>%
  st_cast("MULTILINESTRING")%>%
  st_transform(3857)

data_sf<-data %>%
  st_as_sf( coords = c("long", "lat"), crs=4326) %>%
  st_transform(3857)

# extract background map tiles - check the provider argument, there are many options
basemap <- basemap_ggplot(ext=st_bbox(data_sf), map_service = "esri", map_type = "world_dark_gray_base")



# create a plot 
basemap +
  geom_sf(data = magField_sf, aes(color = 'Experienced Magnetic Field in Autumn'), size=1.5) +
  geom_sf(data = springMig_sf,aes(color = 'Spring Migration')) +
  coord_sf(default_crs = sf::st_crs(3857)) +
  facet_wrap(~id) +
  labs(title = "Experienced Magnetic Field During Autumn Migration and Spring Migration route for ten Red Kites",
       x = "Longitude",
       y = "Latitude") +
  scale_color_manual(name = 'Legend', values = c('Experienced Magnetic Field in Autumn' = 'forestgreen', 'Spring Migration' = 'firebrick')) +
  theme(legend.position = c(0.75,0.12),  
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/REKI_magnetic_field_autumn_spring.png', height = 9, width = 12)



