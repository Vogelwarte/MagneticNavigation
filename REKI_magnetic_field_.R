
#### Script for studying magnetic orientation of raptors, here REKI 
#### Script written by Filibert Heim, filibert.heim@posteo.de, following instructions from Steffen Oppel and hints from Will Scheider 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load required packages 
library(oce)
library(tmap)
library(sf)
library(terra)
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
library(tidyverse)
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

# calculate mean spring migration date a
data <- data %>%
  group_by(id) %>%
  mutate(mean_spring = if_else(first_migration == 'spring', 
                               mean(timestamp[first_migration == 'spring'], na.rm = TRUE), 
                               as.POSIXct(NA)))

# create df to extract magnetic field values for spring migration for every bird 

birds <- unique(data$id) # get ids from all birds 
mean_spring <- data %>% filter(first_migration == 'spring') %>% group_by(id) %>% summarise(mean_spring = mean(timestamp)) # get mean spring migration dates 

df_rast <- as.data.frame(xyFromCell(rast, 1:ncell(rast))) # get the coordinates from raster
df_rast <- df_rast %>% rename(long = x, lat = y) # rename raster

df_rast <- data.frame(long = rep(df_rast$long, times = length(birds)), 
                      lat = rep(df_rast$lat, times = length(birds)), 
                      id = rep(birds, each = nrow(df_rast)),
                      inclination = NA, 
                      declination = NA, 
                      intensity = NA)

df_rast <- df_rast %>% left_join(mean_spring, by = join_by(id)) # add the mean_spring migration dates 

# get magnetic field values for mean date of each birds first spring migration
df_rast$inclination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$inclination
df_rast$declination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$declination
df_rast$intensity <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$intensity




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Get magnetic field values for first autumn migration for all locations of each bird --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# filter for autumn migration
data_autumn <- data %>% 
  filter(first_migration == 'autumn') %>%
  select(-migration, -first_migration, -distance_to_next, -speed,-mean_spring)

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

# PAIR EVERY GRID CELL IN SPRING WITH EVERY LOCATION IN AUTUMN TO SELECT CELLS WITH 'experienced' MAGNETIC FIELD
overlap <- df_rast %>% select(-mean_spring) %>% 
  full_join(data_autumn %>% select(-long, -lat), by = 'id', suffix = c('_spring', '_autumn'), relationship = "many-to-many") %>% 
  mutate(match_inclination = inclination_spring >= inclination_min & inclination_spring <= inclination_max, # create a column that checks whether inclination from spring lies between max and min of autumn migration
         match_declination = declination_spring >= declination_min & declination_spring <= declination_max,
         match_intensity = intensity_spring >= intensity_min & intensity_spring <= intensity_max, 
         match_all = match_inclination & match_declination & match_intensity) %>% # creates a TRUE if all other columns above are TRUE
  select(-match_inclination, -match_declination, -match_intensity) %>% 
  filter(match_all == TRUE) # remove all grid cells that are not needed

dim(overlap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Plot data --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a plot 
ggplot() + 
  geom_point(data = overlap, 
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




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Plot data on a map background --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## convert data to spatial objects
str(overlap)
magField_sf<-overlap %>%
  st_as_sf(coords = c("long", "lat"), crs=4326)


springMig_sf<-data %>% filter(first_migration == 'spring') %>%
  st_as_sf( coords = c("long", "lat"), crs=4326) %>% 
  group_by(id) %>%
  #arrange(timestamp) %>%
  summarize() %>%
  st_cast("MULTILINESTRING") 


# extract background map tiles - check the provider argument, there are many options
basemap <- maptiles::get_tiles(x = bbox, 
                               zoom = 3,
                               crop = TRUE, provider = "OpenTopoMap")

# plot map
tmap_mode("plot")
tm_shape(basemap)+
  tm_rgb()+
  tm_shape(magField_sf)  +
  tm_symbols(col = "green", size = 2, alpha=0.5) +
  tm_shape(springMig_sf)  +
  tm_lines(col = "red", size = 0.2) +
  tm_facets(by="id", free.scales=FALSE)



# create a plot 
ggmap(basemap) + 
  geom_point(data = overlap, 
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



