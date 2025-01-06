##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## EXTRACTING EGYPTIAN VULTURE DATA FOR NAVIGATION EXPLORATION #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### facilitate collaboration with Richard Holland on bird navigation
### Steffen Oppel developed a quantitative migration index - in 'EGVU_weekly_migration_index.r'

#### Script for studying magnetic orientation of raptors, here EGVU
#### Script written by Filibert Heim, filibert.heim@posteo.de, following instructions from Steffen Oppel and hints from Will Schneider

### goal is to quantify what proportion of Egyptian Vultures migrate back within magnetic field experienced on first autumn migration

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
library(move2)
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
library(tidyverse)
library(units)
select <- dplyr::select
filter <- dplyr::filter 
rename <- dplyr::rename
library(PWFSLSmoke)  ## necessary for logger.info function


# load movebank filter functions
#source("C:/Users/filib/Documents/Praktika/Sempach/Tracking data/MoveApps_filter_functions.r") ### loads outlier and speed filter functions from MoveApps
source("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/DataPrep/MoveApps_filter_functions.r") ### loads outlier and speed filter functions from MoveApps

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Data loading  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in data - here, we are not working anymore with solely 10 bird but rather with (all/most) data available
locs <- read.csv(file = 'data/EGVU_sample_locations_all.csv', sep = ',', dec = '.', header = T) # exported data from Steffen 
mig <- fread('data/EGVU_mig_predictions.csv') %>% 
  mutate(MIG_MAN = if_else(MIG_MAN == '', NA, MIG_MAN))
birds <- read.csv('data/Neophron percnopterus Bulgaria_Greece-reference-data.csv')
str(locs) # check data structure
str(mig)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Data formatting  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# improve format of BIRDS
birds <- birds[1:118,] # remove rows with nonsense data
birds <- birds %>% rename(bird_name = animal.id, tag_id = tag.id)

# get only juvenile birds and solve issues with duplicate, nearly identical data entries 
juv <- unique(c(birds[grep(pattern = 'captive|headstar', x = birds$animal.group.id),] %>% pull(bird_name),  # get all birds that were tagged as juv. for sure
                birds[grep(pattern = 'juv|juv.|juvenile|chick|nestling', x = birds$animal.comments),] %>% pull(bird_name), 
                birds[grep(pattern = 'juvenile|nestling', x = birds$animal.life.stage),] %>% pull(bird_name)))

# find all duplicates in this data (occurs where tags where removed and replaced)
birds %>% group_by(bird_name) %>% summarise(duplicate = length(bird_name)) %>% filter(duplicate > 1) %>% pull(bird_name) # find all duplictes 
birds_juv <- birds %>% 
  filter(bird_name %in% juv) %>% 
  filter(!tag_id == c('171354', '182250')) %>% # delete duplicates (attention: for our proposes it does not make a difference which rows are deleted, except for Izi, do not delete tag_id 221441), also delete Dalol 182250 as he has been caugt in Africa (no data for autumn migration)
  group_by(bird_name) %>% slice(1) 

birds_juv <- birds_juv %>% filter(!bird_name == 'Adi') # Adi has been caught as ad, but in comments a nstling is mentioned and thus removed manually 

# continue with normal formatting of bird data 
birds_juv <- birds_juv %>%
  mutate(latest_date_born = if_else(animal.exact.date.of.birth == '', animal.latest.date.born, animal.exact.date.of.birth), 
         sex = if_else(animal.sex == '', NA, animal.sex), # remove '' and replace by NA 
         birth_lat = animal.birth.hatch.latitude, 
         birth_long = animal.birth.hatch.longitude, 
         group = animal.group.id, 
         deploy_date = deploy.on.date) %>%
  select(bird_name, tag_id, latest_date_born, sex, birth_lat, birth_long, group, deploy_date)

# improve format for MIG 
mig <- mig %>% 
  mutate(date = as_date(date),
         migration = factor(if_else(MIG_PRED > 0.5, 'migratory', 'stationary')), # create factorial variable migration 
         migration_error = if_else(migration == MIG_MAN,"correct","false"), # create flag indicating whether classification is correct or not depending on manual annotation, if available 
         migration_error = if_else(is.na(MIG_MAN), 'automated_classification', migration_error), # update migration error so that all NAs (non-manullay classified data) dont cause falses
         migration = factor(if_else(migration_error == 'false', MIG_MAN, migration)), # replace all cases where errors occured in the automated classification with manual classification
         bird_name = factor(bird_name), 
         ymo = paste0(year(date), '_', format(date, '%m'))) %>% 
  select(-MIG_PRED, -MIG_MAN)

# improve format for LOCS
locs <- locs %>% mutate(timestamp = as_datetime(timestamp), 
                        bird_name = as.factor(individual_local_identifier), 
                        bird_name2 = as.factor(individual_local_identifier),
                        tag_id = as.factor(tag_local_identifier), 
                        study_id = as.factor(study_id), 
                        ymo = paste0(year(timestamp), '_', format(timestamp, '%m'))) %>%
  rename(long = location_long, 
         lat = location_lat) %>% 
  select(-individual_local_identifier, -tag_local_identifier, -visible)

# make tracking data usable by creating a move2 object
locs_sf <- st_as_sf(locs, coords = c('long', 'lat'), crs = 4326, remove = F)
locs_move2 <- mt_as_move2(locs_sf, time_column = 'timestamp', track_id_column = 'bird_name2')

# # apply filtering functions from moveapps
# # remove outliers with function loaded from local file 
locs_filt <- rmOutlier(locs_move2, maxspeed = 35, MBremove = T)
dim(locs_filt)

# remove duplicate locations
locs_filt <- mt_filter_unique(locs_filt)
dim(locs_filt)# no duplicates seem to exist

# remove missing data and remove probably wrong locations in unlikely parts of the world
mt_track_id(locs_filt) <- NULL  ## converts the move2 object to an sf object
locs_filt <- locs_filt %>% 
  mutate(long = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>% 
  filter(!is.na(timestamp)) %>%
  filter(!is.na(lat)) %>% 
  filter(!is.na(long)) %>% 
  filter(long > -20, long < 60, lat > -35, lat < 56) %>% # this captures basically everything between Denmark and South Africa
  # mutate(age_cy = as.integer(as.integer(as.Date(timestamp)-as.Date(latest_date_born)))/365) %>% # calculate age_cy 
  select(bird_name,tag_id, timestamp, lat, long, ymo) %>% # select needed columns and remove unneeded bird_id's 
  filter(!is.na(bird_name)) %>%
  st_transform(3857)

head(locs_filt)
dim(locs_filt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bring data together from mig, locs and birds  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add automatically annotated migration data and bird data 
data <- locs_filt %>% 
  left_join(mig %>% select(bird_name, ymo, migration), by = join_by(bird_name, ymo)) %>% 
  inner_join(birds_juv, by = join_by(bird_name, tag_id)) %>% # join bird data but keep only these birds where data of birds exists (this removes all ad and immat birds since they are not included in birds_juv)
  filter(!is.na(migration))

dim(data)




# 
# 
# 
# # this code creates consecutive labelling for the migratory periods - 1 is first autumn migration, 2 is second block of migration
# # this does not work and produces only NA
# data <- data %>%
#   group_by(bird_name) %>%
#   mutate(migration_flag = if_else(migration == "migratory" & (dplyr::lag(migration, default = NA_character_) %in% c(NA, "stationary")),1, NA_real_), # finds all the transitions between NA/stationary and migratory and labels them with 1
#          first_migration = if_else(migration == 'migratory', cumsum(coalesce(migration_flag, 0)), NA_real_)) # labels all continuous blocks 
# 
# unique(data$first_migration)
# # this produces an overview table to verify the first_migration column 
# data %>% group_by(bird_name) %>% 
#   summarise(latest_date_born = first(date(latest_date_born)), 
#             deploy_date = first(date(deploy_date)), 
#             first_50_apart = min(date(timestamp)),
#             max_time = max(date(timestamp)), 
#             age_year = as.numeric(difftime(max_time, latest_date_born, units = 'days'))/365, 
#             number_of_migs = max(as.numeric(first_migration), na.rm = T)) %>% print(n = 76)
# 
# # as seen above, I can remove all birds with less than 2 migration periods 
# birds_mig <- data %>% group_by(bird_name) %>% # get all birds with more than 2 migratory periods
#   summarise(number_of_migs = max(as.numeric(first_migration), na.rm = T)) %>% filter(number_of_migs >= 2)  %>% pull(bird_name)
# data_mig <- data %>% filter(bird_name %in% birds_mig)
# 
# # check results - sometimes there are gaps in migratory periods which means that I have to aggregate those
# data_mig %>% drop_na(first_migration) %>% 
#   group_by(bird_name, first_migration) %>%
#   summarise(latest_date_born = first(date(latest_date_born)), 
#             deploy_date = first(date(deploy_date)), 
#             start_mig = min(date(timestamp)),
#             end_mig = max(date(timestamp)))
# 
# # classify the migratory periods depending on time of the year and try to consider individuals which stayed longer in Africa without migrating back
# data_mig <- data_mig %>% 
#   group_by(bird_name) %>% 
#   mutate(first_migration = if_else(year(timestamp) == year(deploy_date) & migration == 'migratory' | ymd(paste(year(deploy_date) + 1, '-01-31')) >= timestamp & migration == 'migratory', 1, first_migration), # aggregate migratory blocks if they were in the same year or the consecutive year until the 31st of Jan, I use deploy_date to consider headstarts
#          first_migration = if_else(ymd(paste(year(deploy_date) + 1, '-01-31')) < timestamp & timestamp <  ymd(paste(year(deploy_date) + 1, '-07-01')) & migration == 'migratory', 2, first_migration), # aggregate all migratory periods from deploy_date +1 and aggregate them to spring migrations // looks like this is enough
#          first_migration = if_else(ymd(paste(year(deploy_date) + 2, '-01-31')) < timestamp & timestamp <  ymd(paste(year(deploy_date) + 2, '-07-01')) & migration == 'migratory' & first_migration != 4, 2, first_migration))
# # attention - the classification of migratory periods (more in detail their consecutive numbering) is inconsistent after the 3rd migration, but since we are only interested in the first migration, this doesn't matter
# 
# # remove all data before the start of the 3rd migration
# data_mig <- data_mig %>%
#   group_by(bird_name) %>%
#   mutate(first_migration_3_time = min(timestamp[first_migration == 3], na.rm = TRUE),  # get min timestamp of 3rd and 4th migratory period 
#          first_migration_4_time = min(timestamp[first_migration == 4], na.rm = TRUE)) %>%  
#   filter(timestamp < first_migration_3_time, timestamp < first_migration_4_time) %>%  # keep only rows before these two timestamps
#   select(-first_migration_3_time, -first_migration_4_time, -migration_flag)
# 
# # relabel first migration and second migration as autumn and spring migrations
# data_mig <- data_mig %>% mutate(first_migration = factor(case_when(first_migration == 1 ~ 'autumn', 
#                                                                    first_migration == 2 ~ 'spring', 
#                                                                    TRUE ~ as.character(first_migration))))
# head(data_mig)
# 
# 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Create a raster for the migration route of EGVU  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create grid across Europe bounded by data extent
flyway <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(3857) %>%
  st_crop(y=st_bbox(data))

grid_flyway <- flyway %>%
    st_make_grid(cellsize = 50000, what = "polygons",
                 square = FALSE) %>% # This statements leads to hexagons
  st_sf() %>% ## turn it into a data frame to add ID column
  mutate(poly_id=seq_along(geometry))

# extract background map tiles - check the provider argument, there are many options
basemap <- basemap_ggplot(ext=st_bbox(data), map_service = "esri", map_type = "world_dark_gray_base")
basemap +
  geom_sf(data = grid_flyway, color = 'lightgreen', size=1.5, fill=NA)


# create a dataframe with lat and long of grid cells
df_rast <- st_centroid(grid_flyway) %>%
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
  group_by(bird_name, mig_year, season) %>%
  summarise(locs=sum(n), start=min(timestamp),end=max(timestamp))

nav_out_spring<-nav_out_summary %>%
  filter(season=="spring") %>%
  select(bird_name,mig_year,locs,start,end) %>%
  rename(spring_locs=locs,spring_start=start,spring_end=end, spring_year=mig_year) %>%
  group_by(bird_name) %>%
  slice_min(spring_year,n=2)
nav_out_autumn<-nav_out_summary %>%
  filter(season=="autumn") %>%
  select(bird_name,mig_year,locs,start,end) %>%
  rename(autumn_locs=locs,autumn_start=start,autumn_end=end,autumn_year=mig_year) %>%
  group_by(bird_name) %>%
  slice_min(autumn_year)

 
### FIND ANIMALS THAT HAVE BOTH AUTUMN AND SPRING MIGRATION AND PROCESS DATA FOR THOSE
nav_out_summary<-nav_out_autumn %>%
  full_join(nav_out_spring, by=c("bird_name")) %>%
  filter(!is.na(spring_year)) %>%
  filter(!is.na(autumn_year)) %>%
  mutate(mean_spring=spring_start+((spring_end-spring_start)/2)) %>%
  mutate(migID=paste(bird_name,spring_year, sep="_"))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. LOOP OVER ALL INDIVIDUALS TO EXTRACT FIRST AUTUMN AND SPRING MIGRATION  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OUTPUT<-data.frame()
for (i in nav_out_summary$migID) {
  
  mi<-nav_out_summary %>% filter(migID==i)
  if(i=="Awash_2020"){mi$autumn_end<-mi$spring_start}
  
  # get magnetic field values for mean date of each birds first spring migration across entire migration flyway raster
  df_rast$inclination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = mi$mean_spring)$inclination
  df_rast$declination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = mi$mean_spring)$declination
  df_rast$intensity <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = mi$mean_spring)$intensity
  
  # select locations in first autumn migration with magnetic field values for the exact timestamp
  data_autumn <- data %>%
    #st_drop_geometry() %>%
    filter(bird_name==mi$bird_name) %>% 
    filter(timestamp>=mi$autumn_start) %>% 
    filter(timestamp<=mi$autumn_end)
  
  # get magnetic field values for actual date of each birds first autumn migration
  data_autumn$inclination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$inclination
  data_autumn$declination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$declination
  data_autumn$intensity <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$intensity
  
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
  experienced_mag_field <- df_rast %>% mutate(bird_name=mi$bird_name) %>%
    full_join(data_autumn, by = 'bird_name', suffix = c('_spring', '_autumn'), relationship = "many-to-many") %>% 
    mutate(match_inclination = inclination_spring >= inclination_min & inclination_spring <= inclination_max, # create a column that checks whether inclination from spring lies between max and min of autumn migration
           match_declination = declination_spring >= declination_min & declination_spring <= declination_max,
           match_intensity = intensity_spring >= intensity_min & intensity_spring <= intensity_max, 
           match_all = match_inclination & match_declination & match_intensity) %>% # creates a TRUE if all other columns above are TRUE
    select(-match_inclination, -match_declination, -match_intensity) %>% 
    filter(match_all == TRUE) # remove all grid cells that do not overlap with previously experienced magnetic field values
  
  ## convert data to spatial objects
  magField_sf<-grid_flyway %>%
    filter(poly_id %in% unique(experienced_mag_field$poly_id))
  
  ####### EXTRACT SPRING MIGRATION DATA FOR THIS INDIVIDUAL
  data_spring <- data %>%
    filter(bird_name==mi$bird_name) %>% 
    filter(timestamp>=mi$spring_start) %>%
    filter(timestamp<=mi$spring_end) %>% arrange(timestamp) %>%
    st_join(., magField_sf)

  ####### SUMMARISE DATA FOR THIS INDIVIDUAL
  mi$prop_locs_in_mag<-nrow(data_spring %>% filter(is.na(poly_id)))/nrow(data_spring)
  OUTPUT<-bind_rows(OUTPUT,mi)
    
  ####### CREATE MULTILINE STRING FOR BOTH MIGRATIONS
  data_sf<-bind_rows((data_spring %>% mutate(season="spring")),
                     (data_autumn %>% mutate(season="autumn"))) %>%
    mutate(long = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
    group_by(bird_name,season) %>%
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
    geom_text(aes(x=4900000, y=6000000, label=i), colour="gold",size=6, fontface="bold") +
    theme(legend.position="inside",
          legend.position.inside = c(0.17,0.87),  
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12), 
          plot.title = element_text(size = 15), 
          axis.title = element_text(size = 12))
  ggsave(filename = sprintf('output/EGVU_%s_magnetic_field_migration.png',i), height = 10, width = 12)
  
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. SUMMARIZE DATA ACROSS ALL INDIVIDUALS --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## not a single bird left the magnetic field already experienced in autumn!
mean(OUTPUT$prop_locs_in_mag)
fwrite(OUTPUT,"output/EGVU_first_mig_summaries.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## COMPILE ALL GRAPHS INTO ONE FILE ####### --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### resized images with https://imageresizer.com/bulk-resize
library(png)
library(grid)
library(gridExtra)
plots<-list.files("output/.", pattern=glob2rx("EGVU*.png"))
plots
pdf("output/EGVU_migration_mag_field_all.pdf")
for (i in 1:length(plots)) {
  p1 <- readPNG(paste0("output/",plots[i]), native = FALSE)
  grid.raster(p1, x = unit(0.5, "npc"), width=unit(1, "npc"))
  plot.new()
}
dev.off()