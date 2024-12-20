#### Script for Creating Magnetic Envelopes of first Autumn Migration and plotting first Spring Migration above - this is the result from a chat with Richard Holland, Will Schneiders and Steffen Oppel
#### Studying magnetic orientation of raptors, here EGVU
#### Script written by Filibert Heim, filibert.heim@posteo.de, in Dez 20124 following instructions from Steffen Oppel and hints from Will Schneider and Richard Holland 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load required packages 
library(oce)
library(sf)
library(terra)
library(readxl)
library(lubridate)
library(move2)
library(data.table)
library(basemaps)
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
library(tidyverse)
select <- dplyr::select
filer <- dplyr::filter 
rename <- dplyr::rename

# load movebank filter functions
source("C:/Users/filib/Documents/Praktika/Sempach/Tracking data/MoveApps_filter_functions.r") ### loads outlier and speed filter functions from MoveApps

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

# find all duplicates in this data (occures where tags where removed and replaced)
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

# apply filtering functions from moveapps
# remove outliers with function loaded from local file 
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
  filter(!is.na(bird_name))
locs_filt <- as.data.frame(locs_filt) %>% select(-geometry) # convert st object to normal data frame
head(locs_filt)
dim(locs_filt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Bring data together from mig, locs and birds  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add automatically annotated migration data and bird data 
data <- locs_filt %>% 
  left_join(mig %>% select(bird_name, ymo, migration, migration_error), by = join_by(bird_name, ymo)) %>% 
  inner_join(birds_juv, by = join_by(bird_name, tag_id)) %>% # join bird data but keep only these birds where data of birds exists (this removes all ad and immat birds since they are not included in birds_juv)
  mutate(migration = factor(migration))

# NAs are introduced, see where they come from 
data %>% filter(is.na(migration)) %>% print() # if possible find a way to deal with them - possibly also in the migration classification script?

# this code creates consecutive labelling for the migratory periods - 1 is first autumn migration, 2 is second block of migration 
data <- data %>%
  group_by(bird_name) %>%
  mutate(migration_flag = if_else(migration == "migratory" & (lag(migration, default = NA_character_) %in% c(NA, "stationary")),1, NA_real_), # finds all the transitions between NA/stationary and migratory and labels them with 1
         first_migration = if_else(migration == 'migratory', cumsum(coalesce(migration_flag, 0)), NA_real_)) # labels all continuous blocks 

# this produces an overview table to verify the first_migration column 
data %>% group_by(bird_name) %>% 
  summarise(latest_date_born = first(date(latest_date_born)), 
            deploy_date = first(date(deploy_date)), 
            first_50_apart = min(date(timestamp)),
            max_time = max(date(timestamp)), 
            age_year = as.numeric(difftime(max_time, latest_date_born, units = 'days'))/365, 
            number_of_migs = max(as.numeric(first_migration), na.rm = T)) %>% print(n = 76)

# as seen above, I can remove all birds with less than 2 migration periods 
birds_mig <- data %>% group_by(bird_name) %>% # get all birds with more than 2 migratory periods
  summarise(number_of_migs = max(as.numeric(first_migration), na.rm = T)) %>% filter(number_of_migs >= 2)  %>% pull(bird_name)
data_mig <- data %>% filter(bird_name %in% birds_mig)

# check results - sometimes there are gaps in migratory periods which means that I have to aggregate those
data_mig %>% drop_na(first_migration) %>% 
  group_by(bird_name, first_migration) %>%
  summarise(latest_date_born = first(date(latest_date_born)), 
            deploy_date = first(date(deploy_date)), 
            start_mig = min(date(timestamp)),
            end_mig = max(date(timestamp)))

# classify the migratory periods depending on time of the year and try to consider individuals which stayed longer in Africa without migrating back
data_mig <- data_mig %>% 
  group_by(bird_name) %>% 
  mutate(first_migration = if_else(year(timestamp) == year(deploy_date) & migration == 'migratory' | ymd(paste(year(deploy_date) + 1, '-01-31')) >= timestamp & migration == 'migratory', 1, first_migration), # aggregate migratory blocks if they were in the same year or the consecutive year until the 31st of Jan, I use deploy_date to consider headstarts
         first_migration = if_else(ymd(paste(year(deploy_date) + 1, '-01-31')) < timestamp & timestamp <  ymd(paste(year(deploy_date) + 1, '-07-01')) & migration == 'migratory', 2, first_migration), # aggregate all migratory periods from deploy_date +1 and aggregate them to spring migrations // looks like this is enough
         first_migration = if_else(ymd(paste(year(deploy_date) + 2, '-01-31')) < timestamp & timestamp <  ymd(paste(year(deploy_date) + 2, '-07-01')) & migration == 'migratory' & first_migration != 4, 2, first_migration))
# attention - the classification of migratory periods (more in detail their consecutive numbering) is inconsistent after the 3rd migration, but since we are only interested in the first migration, this doesn't matter

# remove all data before the start of the 3rd migration
data_mig <- data_mig %>%
  group_by(bird_name) %>%
  mutate(first_migration_3_time = min(timestamp[first_migration == 3], na.rm = TRUE),  # get min timestamp of 3rd and 4th migratory period 
         first_migration_4_time = min(timestamp[first_migration == 4], na.rm = TRUE)) %>%  
  filter(timestamp < first_migration_3_time, timestamp < first_migration_4_time) %>%  # keep only rows before these two timestamps
  select(-first_migration_3_time, -first_migration_4_time, -migration_flag)

# relabel first migration and second migration as autumn and spring migrations
data_mig <- data_mig %>% mutate(first_migration = factor(case_when(first_migration == 1 ~ 'autumn', 
                                             first_migration == 2 ~ 'spring', 
                                            TRUE ~ as.character(first_migration))))
head(data_mig)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Create a raster for the migration route of EGVU's  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create binding box 
bbox <- terra::ext(range(data_mig$long), range(data_mig$lat))

# create raster with geographic CRS for Europe where grid cell size varies with location
rast <- rast(ext = bbox, resolution = 0.5, crs = "EPSG:4326") # not sure if maybe the world EPSG would be better, just replace by 'EPSG:4326'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Get magnetic field values for first spring migration in raster  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate mean autumn migration date a
data_mig <- data_mig %>%
  group_by(bird_name) %>%
  mutate(mean_spring = if_else(first_migration == 'spring', 
                               mean(timestamp[first_migration == 'spring'], na.rm = TRUE), 
                               as.POSIXct(NA)))

# create df to extract magnetic field values for autumn migration for every bird 
birds <- unique(data_mig$bird_name) # get ids from all birds 
mean_spring <- data_mig %>% filter(first_migration == 'spring') %>% group_by(bird_name) %>% summarise(mean_spring = mean(timestamp)) # get mean autumn migration dates 

df_rast <- as.data.frame(xyFromCell(rast, 1:ncell(rast))) # get the coordinated from raster
df_rast <- df_rast %>% rename(long = x, lat = y) # rename raster

df_rast <- data.frame(long = rep(df_rast$long, times = length(birds)), 
                      lat = rep(df_rast$lat, times = length(birds)), 
                      bird_name = rep(birds, each = nrow(df_rast)),
                      inclination = NA, 
                      declination = NA, 
                      intensity = NA)

df_rast <- df_rast %>% left_join(mean_spring, by = join_by(bird_name)) # add the mean_spring migration dates 

# check how many cases there are where no spring migration occurred yet, remove them because otherwise error in magneticField() extraction
df_rast <- df_rast %>% drop_na(mean_spring) # get rid of all cases where no spring migration date could be calculates because of missing data 

# get magnetic field values for mean date of each birds first spring migration
df_rast$inclination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$inclination
df_rast$declination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$declination
df_rast$intensity <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$intensity

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Get magnetic field values for first autumn migration for all locations of each bird --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# filter for autumn migration
data_autumn <- data_mig %>% 
  filter(first_migration == 'autumn') 

# get magnetic field values for all locations in first autumn migration for the exact timestamp
data_autumn$declination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$declination
data_autumn$inclination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$inclination
data_autumn$intensity <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$intensity

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Get information of grid cells where autumn migration values lied within spring migration  --------
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

# create one df of data_autumn and df_rast from spring migration for comparison
overlap <- df_rast %>% 
  full_join(data_autumn %>% select(-long, -lat), by = join_by(bird_name), suffix = c('_spring', '_autumn'), relationship = 'many-to-many') %>% 
  select(-migration, -first_migration) %>% 
  mutate(match_inclination = inclination_spring >= inclination_min & inclination_spring <= inclination_max, # create a column that checks whether inclination from spring lies between max and min of autumn migration
         match_declination = declination_spring >= declination_min & declination_spring <= declination_max,
         match_intensity = intensity_spring >= intensity_min & intensity_spring <= intensity_max, 
         match_all = match_inclination & match_declination & match_intensity) %>% # creates a TRUE if all other columns above are TRUE
  select(-match_inclination, -match_declination, -match_intensity) %>% 
  filter(match_all == TRUE) # romve all grid cells that are not needed

dim(overlap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Plot data --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# order data frames for appropriate plotting with geom_line
overlap <- overlap %>% group_by(bird_name) %>% arrange(bird_name, timestamp)
data_mig <- data_mig %>% group_by(bird_name) %>% arrange(bird_name, timestamp)

# create a plot 
ggplot() + 
  geom_point(data = overlap, 
             mapping = aes(x = long, y = lat, color = 'Experienced Magentic Field in Autumn'), size = 1, shape = 15) +
  geom_point(data = data_mig %>% mutate(first_migration = case_when(first_migration == 'autumn' ~ 'First Autumn Migration', first_migration == 'spring' ~ 'First Spring Migration', is.na(first_migration) ~ 'Stationary Period')), 
            mapping = aes(x = long, y = lat, color = first_migration), size = .4) +
  coord_equal() + 
  facet_wrap(~bird_name) +
  labs(title = "Experienced Magnetic Field During First Autumn Migration and Spring Migration Route in Juvenile Egyptian Vultures",
    x = "Longitude",
    y = "Latitude") +
  scale_color_manual(name = 'Legend', values = c('Experienced Magentic Field in Autumn' = 'grey80', 'First Spring Migration' = 'firebrick', 'First Autumn Migration' = 'blue', 'Stationary Period' = 'grey30')) +
  theme_bw() + 
  theme(legend.position = 'right', 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/EGVU_magnetic_field_autumn_spring.png', height = 9, width = 12)

# check why only a small number of birds appear in plot 
length(unique(data_mig$bird_name)) # its a bit complicated: we have 111 levels, but only 86 of them have data 

# create plots in loop for a few selected individuals 
birds <- c('Tatul', 'Dobromir', 'Sanie', 'Polya')
for(i in 1:length(birds)){
  ggplot() + 
    geom_point(data = overlap %>% filter(bird_name == birds[i]), 
               mapping = aes(x = long, y = lat, color = 'Experienced Magentic Field in Autumn'), size = 3.8, shape = 15) +
    geom_point(data = data_mig %>% filter(bird_name == birds[i]) %>% mutate(first_migration = case_when(first_migration == 'autumn' ~ 'First Autumn Migration', first_migration == 'spring' ~ 'First Spring Migration', is.na(first_migration) ~ 'Stationary Period')), 
               mapping = aes(x = long, y = lat, color = first_migration), size = 2.5) +
    coord_equal() + 
    facet_wrap(~bird_name) +
    labs(title = paste0("Magnetic Field and Migration in Juvenile Egyptian Vulture ", birds[i]),
         x = "Longitude",
         y = "Latitude") +
    scale_color_manual(name = 'Legend', values = c('Experienced Magentic Field in Autumn' = 'grey80', 'First Spring Migration' = 'firebrick', 'First Autumn Migration' = 'blue', 'Stationary Period' = 'grey30')) +
    theme_bw() + 
    theme(legend.position = 'right', 
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12), 
          plot.title = element_text(size = 15), 
          axis.title = element_text(size = 12))
  ggsave(filename = sprintf('output/EGVU_magnetic_field_autumn_spring_%s.png', birds[i]), height = 9, width = 12)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Plot data on a map background --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## convert data to spatial objects - only needed if tmap is used
str(overlap)
magField_sf<-overlap %>%
  st_as_sf(coords = c("long", "lat"), crs=4326) %>%
  st_transform(3857)


springMig_sf<-data_mig %>% filter(first_migration == 'spring') %>%
  st_as_sf( coords = c("long", "lat"), crs=4326) %>%
  group_by(bird_name) %>%
  arrange(timestamp) %>%
  summarize() %>%
  # st_cast("MULTILINESTRING")%>% # this ensures, that no weird connections are made between points that don't belong together - this was the quickest way for troubleshooting and there are probably better solutions! 
  st_transform(3857)

# create sf data 
data_sf<-data_mig %>% filter(first_migration == 'spring') %>% 
  st_as_sf( coords = c("long", "lat"), crs=4326) %>%
  st_transform(3857) %>% arrange(bird_name, timestamp)


# extract background map tiles - check the provider argument, there are many options
basemap <- basemap_ggplot(ext=st_bbox(data_sf), map_service = "esri", map_type = "world_dark_gray_base")



# create a plot 
basemap +
  geom_sf(data = magField_sf, aes(color = 'Experienced Magnetic Field in Autumn'), size=1.5) +
  geom_sf(data = springMig_sf, mapping = aes(color = 'Spring Migration')) +
  coord_sf(default_crs = sf::st_crs(3857)) +
  facet_wrap(~bird_name, nrow = 3) +
  labs(title = "Experienced Magnetic Field During Autumn Migration and Spring Migration route for twenty Egyptian Vultures",
       x = "Longitude",
       y = "Latitude") +
  scale_color_manual(name = 'Legend', values = c('Experienced Magnetic Field in Autumn' = 'forestgreen', 'Spring Migration' = 'firebrick')) +
  scale_x_continuous(breaks = seq(10, 40, 10)) +
  theme(legend.position = 'bottom',  
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/EGVU_magnetic_field_autumn_spring_basemap.png', height = 9, width = 12)


