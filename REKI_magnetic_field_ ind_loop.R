##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## EXTRACTING RED KITE DATA FOR NAVIGATION EXPLORATION #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### facilitate collaboration with Richard Holland on bird navigation
### Steffen Oppel developed a quantitative migration index - in 'REKI_weekly_migration_index.r'

#### Script for studying magnetic orientation of raptors, here REKI 
#### Script written by Filibert Heim, filibert.heim@posteo.de, following instructions from Steffen Oppel and hints from Will Schneider

### goal is to quantify what proportion of red kites migrate back within magnetic field experienced on first autumn migration

## revised by Steffen Oppel on 3 Jan 2025 to process birds individually (for easier troubleshooting)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
gc()
# load required packages 
library(oce)
library(tmap)
library(basemaps)  ### loads basemap in EPSG:3857 and therefore all other sf objects must be st_transform(3857)
library(sf)
library(terra)
library(dplyr)
library(dtplyr)
library("rnaturalearth")
library("rnaturalearthdata")
library(gridExtra)
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
library(tidyverse)
library(units)
select <- dplyr::select
filter <- dplyr::filter 
rename <- dplyr::rename

# read in data 
data <- readRDS('data/REKI_sample_tracks100.rds') # this is data which is filtered for 1h intervals and comprises 100 individuals, including REKI_928 and REKI_515 - updated 23 Dec 2024 from 22 to 100 birds
mindist<-25
units(mindist)<-"m"

# improve data structure 
data <- data %>% 
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), 
         migration = as.factor(if_else(MIGRATION == 'migrating', 'migratory', 'stationary'))) %>% 
  mutate(long = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>% 
  filter(distance_to_next>mindist) %>%  ## eliminate duplicate coordinates
  st_transform(3857) %>%
  select(-MIGRATION, -distance_to_next,-speed)

dim(data)
length(unique(data$id))
head(data)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Create a raster for the migration route of REKI's  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create grid across Europe bounded by data extent
europe <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent=="Europe") %>%
  st_transform(3857) %>%
  st_crop(y=st_bbox(data))

grid_EU <- europe %>%
    st_make_grid(cellsize = 50000, what = "polygons",
                 square = FALSE) %>% # This statements leads to hexagons
  st_sf() %>% ## turn it into a data frame to add ID column
  mutate(poly_id=seq_along(geometry))

# extract background map tiles - check the provider argument, there are many options
basemap <- basemap_ggplot(ext=st_bbox(data), map_service = "esri", map_type = "world_dark_gray_base")
basemap +
  geom_sf(data = grid_EU, color = 'lightgreen', size=1.5, fill=NA)


# create a dataframe with lat and long of grid cells
df_rast <- st_centroid(grid_EU) %>%
  st_sf() %>%
  st_transform(4326) %>%
  mutate(long = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
  
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. EXTRACT FIRST AUTUMN AND SPRING MIGRATION FROM ALL JUVENILES --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### NEED TO EXTRACT THE FIRST AUTUMN AND FIRST SPRING MIGRATION OF ALL BIRDS ####

nav_out_summary<-data %>%
  st_drop_geometry() %>%
  filter(migration=="migratory") %>%
  mutate(n=1, month=month(timestamp), mig_year=year(timestamp)) %>%
  mutate(season=ifelse(month<8,"spring","autumn")) %>%
  mutate(mig_year=ifelse(season=="spring", mig_year-1,mig_year)) %>%
  group_by(id, mig_year, season) %>%
  summarise(locs=sum(n), start=min(timestamp),end=max(timestamp))

nav_out_spring<-nav_out_summary %>%
  filter(season=="spring") %>%
  select(id,mig_year,locs,start,end) %>%
  rename(spring_locs=locs,spring_start=start,spring_end=end) 
nav_out_autumn<-nav_out_summary %>%
  filter(season=="autumn") %>%
  select(id,mig_year,locs,start,end) %>%
  rename(autumn_locs=locs,autumn_start=start,autumn_end=end)

 
### FIND ANIMALS THAT HAVE BOTH AUTUMN AND SPRING MIGRATION AND PROCESS DATA FOR THOSE
nav_out_summary<-nav_out_summary %>% select (id, mig_year, season, locs) %>%
  spread(key=season, value=locs) %>%
  ungroup() %>%
  mutate(elim=autumn+spring) %>%
  filter(!is.na(elim)) %>%
  select(-elim,-autumn,-spring) %>%
  arrange(id, mig_year) %>%
  group_by(id) %>%
  slice_min(mig_year) %>%
  left_join(nav_out_autumn, by=c("id","mig_year")) %>%
  left_join(nav_out_spring, by=c("id","mig_year")) %>%
  mutate(mean_spring=spring_start+((spring_end-spring_start)/2))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. LOOP OVER ALL INDIVIDUALS TO EXTRACT FIRST AUTUMN AND SPRING MIGRATION  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OUTPUT<-data.frame()
for (i in nav_out_summary$id) {
  
  mi<-nav_out_summary %>% filter(id==i)
  if(i=="REKI_442"){mi$autumn_end<-mi$spring_start}
  
  # get magnetic field values for mean date of each birds first spring migration across entire migration flyway raster
  df_rast$inclination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = mi$mean_spring)$inclination
  df_rast$declination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = mi$mean_spring)$declination
  df_rast$intensity <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = mi$mean_spring)$intensity
  
  
  
  # select locations in first autumn migration with magnetic field values for the exact timestamp
  data_autumn <- data %>%
    #st_drop_geometry() %>%
    filter(id==i) %>% 
    filter(timestamp>=mi$autumn_start) %>% 
    filter(timestamp<=mi$autumn_end) %>%
    mutate(long = st_coordinates(.)[,1], lat = st_coordinates(.)[,2])
  
  # 4. OVERLAP grid cells where autumn migration values are within spring migration values --------
  
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
  experienced_mag_field <- df_rast %>% mutate(id=i) %>%
    full_join(data_autumn, by = 'id', suffix = c('_spring', '_autumn'), relationship = "many-to-many") %>% 
    mutate(match_inclination = inclination_spring >= inclination_min & inclination_spring <= inclination_max, # create a column that checks whether inclination from spring lies between max and min of autumn migration
           match_declination = declination_spring >= declination_min & declination_spring <= declination_max,
           match_intensity = intensity_spring >= intensity_min & intensity_spring <= intensity_max, 
           match_all = match_inclination & match_declination & match_intensity) %>% # creates a TRUE if all other columns above are TRUE
    select(-match_inclination, -match_declination, -match_intensity) %>% 
    filter(match_all == TRUE) # remove all grid cells that do not overlap with previously experienced magnetic field values
  
  ## convert data to spatial objects
  magField_sf<-grid_EU %>%
    filter(poly_id %in% unique(experienced_mag_field$poly_id))
  
  ####### EXTRACT SPRING MIGRATION DATA FOR THIS INDIVIDUAL
  data_spring <- data %>%
    filter(id==i) %>% 
    filter(timestamp>=mi$spring_start) %>%
    filter(timestamp<=mi$spring_end) %>% arrange(timestamp) %>%
    mutate(long = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
    st_join(., magField_sf)

  ####### SUMMARISE DATA FOR THIS INDIVIDUAL
  mi$prop_locs_in_mag<-nrow(data_spring %>% filter(is.na(poly_id)))/nrow(data_spring)
  OUTPUT<-bind_rows(OUTPUT,mi)
    
  ####### CREATE MULTILINE STRING FOR BOTH MIGRATIONS
  data_sf<-bind_rows((data_spring %>% mutate(season="spring")),
                     (data_autumn %>% mutate(season="autumn"))) %>%
    group_by(id,season) %>%
    arrange(timestamp)
    ### this leads to very weird plotting patterns
    # summarize() %>%
    # st_cast("MULTILINESTRING")
  
  ####### PLOT INDIVIDUAL
  basemap +
    geom_sf(data = magField_sf, color = 'lightgreen', size=1.5, fill='lightgreen', alpha=0.2) +
    geom_path(data = data_sf,aes(x=long, y=lat, color = season), linewidth=0.7) +
    labs(x = "Longitude",
         y = "Latitude") +
    scale_color_manual(name = 'First migration', values = c('autumn' = 'navyblue', 'spring' = 'firebrick')) +
    geom_text(aes(x=1300000, y=7000000, label=i), colour="gold",size=6, fontface="bold") +
    theme(legend.position="inside",
          legend.position.inside = c(0.17,0.87),  
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12), 
          plot.title = element_text(size = 15), 
          axis.title = element_text(size = 12))
  ggsave(filename = sprintf('output/%s_magnetic_field_migration.png',i), height = 10, width = 12)
  
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. SUMMARIZE DATA ACROSS ALL INDIVIDUALS --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## not a single bird left the magnetic field already experienced in autumn!
mean(OUTPUT$prop_locs_in_mag)
range(OUTPUT$prop_locs_in_mag)
fwrite(OUTPUT,"output/REKI_first_mig_summaries.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## COMPILE ALL GRAPHS INTO ONE FILE ####### --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### resized images with https://imageresizer.com/bulk-resize
library(png)
library(grid)
library(gridExtra)
plots<-list.files("output/.", pattern=glob2rx("REKI*.png"))
plots

pdf("output/REKI_migration_mag_field_all.pdf")
for (i in 1:length(plots)) {
  p1 <- readPNG(paste0("output/",plots[i]), native = FALSE)
  grid.raster(p1, x = unit(0.5, "npc"), width=unit(1, "npc"))
  plot.new()
}
dev.off()