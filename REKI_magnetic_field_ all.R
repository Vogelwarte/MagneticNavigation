
#### Script for studying magnetic orientation of raptors, here REKI 
#### Script written by Filibert Heim, filibert.heim@posteo.de, following instructions from Steffen Oppel and hints from Will Schneider 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
gc()
# load required packages 
library(oce)
library(tmap)
library(basemaps)
library(sf)
library(terra)
library(dplyr)
library(dtplyr)
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
library(tidyverse)
select <- dplyr::select
filter <- dplyr::filter 
rename <- dplyr::rename

# read in data 
data <- readRDS('data/REKI_sample_tracks100.rds') # this is data which is filtered for 1h intervals and comprises 100 individuals, including REKI_928 and REKI_515 - updated 23 Dec 2024 from 22 to 100 birds

# check data structure 
head(data)
str(data)

# improve data structure 
data <- data %>% 
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), 
         migration = as.factor(if_else(MIGRATION == 'migrating', 'migratory', 'stationary'))) %>% 
  cbind(st_coordinates(.)) %>%
  rename(long = X, lat = Y) %>% 
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
  dplyr::filter(migration == 'migratory') %>% # if all tracks should be included, remove this filter 
  mutate(
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)), 1, first_migration), 
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)) + 1 & timestamp < as.POSIXct(paste0(year(min(timestamp, na.rm = TRUE)) + 1, '-07-15 00:00:00'), format = "%Y-%m-%d %H:%M:%S"), 2, first_migration)) %>%
  filter(first_migration %in% c(1,2)) %>%
  select(id, timestamp, first_migration) %>% 
  right_join(as.data.frame(data) %>%
               select(-first_migration, -geometry), by = join_by(id, timestamp)) %>%
  st_as_sf() # troubleshooting issues with right join and sf object - temporally convert to df and back to sf afterwards

data <- data %>% mutate(first_migration = factor(case_when(first_migration == 1 ~ 'autumn', 
                                             first_migration == 2 ~ 'spring', 
                                            TRUE ~ as.character(first_migration))))


### NEED TO EXTRACT THE FIRST AUTUMN AND FIRST SPRING MIGRATION OF ALL BIRDS ####
## somehow the script only retains 10 spring migrations...



data %>%
  mutate(month=month(timestamp)) %>%
  filter(month %in% c(1,2,3,4,5,6)) %>%
  dplyr::filter(migration == 'migratory')







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Create a raster for the migration route of REKI's  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create binding box 
bbox <- terra::ext(range(data$long, na.rm = T), range(data$lat, na.rm = T))

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
mean_spring <- data %>%
  filter(first_migration == 'spring') %>%
  group_by(id) %>%
  summarise(mean_spring = mean(timestamp)) # get mean spring migration dates 

df_rast <- as.data.frame(xyFromCell(rast, 1:ncell(rast))) # get the coordinates from raster
df_rast <- df_rast %>%
  rename(long = x, lat = y) # rename raster

df_rast <- data.frame(long = rep(df_rast$long, times = length(birds)), 
                      lat = rep(df_rast$lat, times = length(birds)), 
                      id = rep(birds, each = nrow(df_rast)),
                      inclination = NA, 
                      declination = NA, 
                      intensity = NA)

df_rast <- df_rast %>%
  left_join(mean_spring, by = join_by(id)) # add the mean_spring migration dates 

df_rast <- df_rast %>%
  drop_na(mean_spring) # drop all NAs to avoid error in magneticField()

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
head(overlap)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. SUMMARIZE OVERLAP of what proportion of red kites migrated fully within the autumn magnetic field --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

overlap


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Plot data --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

length(unique(df_rast$id))

# create a plot 
ggplot() + 
  geom_point(data = overlap, 
             mapping = aes(x = long, y = lat, color = 'Experienced Magnetic Field in Autumn'), size = 1.5) +
  geom_path(data = data %>% filter(first_migration == 'spring'), 
            mapping = aes(x = long, y = lat, color = 'Spring Migration')) +
  coord_equal() + 
  facet_wrap(~id) +
  labs(title = "Experienced Magnetic Field During Autumn Migration and Spring Migration route in Red Kites",
    x = "Longitude",
    y = "Latitude") +
  scale_color_manual(name = 'Legend', values = c('Experienced Magnetic Field in Autumn' = 'grey80', 'Spring Migration' = 'red')) +
  theme_bw() + 
  theme(legend.position = c(0.75,0.12), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/REKI_magnetic_field_autumn_spring_2.png', height = 9, width = 12)

# what happened with REKI_515?
data %>% filter(id == 'REKI_515') %>% slice(5400) %>% filter(timestamp == max(timestamp)) # this is the last signal from REKI_515
ggplot() +
  geom_path(data %>% filter(id == "REKI_515"), mapping = aes(x = long, y = lat), color = "blue", size = 1) +  # Add path with blue color
  geom_point(data %>% filter(id == 'REKI_515') %>% slice(5400) %>% filter(timestamp == max(timestamp)), mapping = aes(x = long, y = lat), color = "red", size = 5)
  # geom_point(color = "red", size = 2) +  # Optional: Add points for better visibility
  labs(title = "Path of REKI_515", x = "Longitude", y = "Latitude") +
  theme_bw()

# looks like we lost signal at wintering site - no spring migration available

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Plot data on a map background --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract background map tiles - check the provider argument, there are many options
basemap <- basemap_ggplot(ext=st_bbox(data_sf), map_service = "esri", map_type = "world_dark_gray_base")

  
  
# create a plot with tracklines for 25 birds at a time
  
allbirds<-unique(data$id)  
plotgroup<-rep(seq(1:9), each=20)

for (p in 1:max(plotgroup)) {
  selbirds<-allbirds[which(plotgroup==p)]
  selbirds<-selbirds[-which(is.na(selbirds))]  ### remove NA when group is smaller than 20
  
  ## convert data to spatial objects
  magField_sf<-overlap %>%
    filter(id %in% selbirds) %>%
    st_as_sf(coords = c("long", "lat"), crs=4326) %>%
    st_transform(3857)
  
  
  data_sf<-data %>%
    filter(id %in% selbirds) %>%
    filter(!is.na(first_migration)) %>%
    st_as_sf( coords = c("long", "lat"), crs=4326) %>%
    group_by(id, first_migration) %>%
    #arrange(timestamp) %>%
    summarize() %>%
    st_cast("MULTILINESTRING")%>%
    st_transform(3857)
  
  
  basemap +
    geom_sf(data = magField_sf, color = 'lightgreen', size=1.5) +
    geom_sf(data = data_sf,aes(color = first_migration), linewidth=0.7) +
    
    facet_wrap(~id) +
    labs(title = "Experienced magnetic field during first autumn migration",
         x = "Longitude",
         y = "Latitude") +
    scale_color_manual(name = 'First migration', values = c('autumn' = 'navyblue', 'spring' = 'firebrick')) +
    theme(legend.position.inside = c(0.75,0.12),  
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12), 
          plot.title = element_text(size = 15), 
          axis.title = element_text(size = 12))
  ggsave(filename = sprintf('output/REKI_magnetic_field_migration_group%s.png',p), height = 10, width = 12)
  
}

######## COMPILE ALL GRAPHS INTO ONE FILE #######

library(png)
library(grid)
library(gridExtra)
plots<-list.files("output/.", pattern=".png")
plots

pdf("output/REKI_migration_mag_field_all.pdf")
for (i in 9:16) {
  p1 <- readPNG(paste0("output/",plots[i]), native = FALSE)
  grid.raster(p1, width=unit(1, "npc"), height= unit(1, "npc"))
  plot.new()
}
dev.off()