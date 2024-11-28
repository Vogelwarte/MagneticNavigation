
#### Script for studying magnetic orientation of raptors, here EGVU 
#### Script written by Filibert Heim, filibert.heim@posteo.de, following instructions from Steffen Oppel and hints from Will Schneider 

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
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
library(tidyverse)
select <- dplyr::select
filer <- dplyr::filter 
rename <- dplyr::rename


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Data loading and formatting  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in data 
locs <- read.csv(file = 'data/EGVU_sample_locations10.csv', sep = ',', dec = '.', header = T) # exported data from Steffen 
birds <- readRDS(file = 'data/EGVU_15min_locs.rds') %>% mt_track_data() %>%# other data drom Steffen just to get information on single individuals
  dplyr::rename(bird_id = individual_local_identifier, comments = individual_comments) %>% 
  mutate(latest_date_born = if_else(is.na(latest_date_born), exact_date_of_birth, latest_date_born)) %>% 
  select(individual_id, deployment_id, comments, bird_id, ring_id, sex, latest_date_born, tag_local_identifier) # include certain columns here if more information on the birds is needed 
mig <- read_excel(path = 'data/Mig months.xlsx', sheet = 'Mig months') 
str(locs) # check data structure
str(birds)

# improve format for locs
locs <- locs %>% mutate(timestamp = as_datetime(timestamp), 
                        bird_name = as.factor(individual_local_identifier), 
                        tag_id = as.factor(tag_local_identifier), 
                        study_id = as.factor(study_id), 
                        ymo = paste0(year(timestamp), '_', format(timestamp, '%m'))) %>% 
  select(-individual_local_identifier, -tag_local_identifier, -visible)

#improve format for birds 
birds <- birds %>% rename(bird_name = bird_id, tag_id = tag_local_identifier) %>% 
  select(bird_name, tag_id, sex, latest_date_born)

# improve format of mig
mig[,c(1,3,17:18)] <- NULL
mig <- mig[grep(pattern = paste0(birds$tag_id, collapse = '|'), x = mig$file_list),] # this filters for all bird_ids (names of the birds) but doesnt catch all of them because some were born in 2020 and the data is from 2019
mig <- mig %>%
  separate(file_list, into = c("bird_name", "tag_id", "year"), sep = "_", extra = "drop") %>%
  mutate(year = gsub("\\.png$", "", year)) %>% 
  pivot_longer(cols = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), names_to = 'month', values_to = 'migration') %>% 
  mutate(month = case_when(month == 'Jan' ~ '01', month == 'Feb' ~ '02', month == 'Mar' ~ '03', month == 'Apr' ~ '04', month == 'May' ~ '05', month == 'Jun' ~ '06', month == 'Jul' ~ '07', month == 'Aug' ~ '08', month == 'Sep' ~ '09',
                           month == 'Oct' ~ '10', month == 'Nov' ~ '11', month == 'Dec' ~ '12')) %>% 
  mutate(ymo = paste0(as.character(year), '_', month)) %>% 
  mutate(migration = case_when(migration == 0 ~ NA, migration == 1 ~ 'stationary', migration == 2 ~ 'migratory'), 
         ymo = paste0(year,'_', month)) 

head(mig)
head(birds, 10)
head(locs)

# bring locs and bird data together, calculate age_cy and long format
data <- locs %>% left_join(birds, by = join_by(bird_name, tag_id)) %>% 
  mutate(age_cy = as.integer(as.integer(timestamp - latest_date_born)/365)) %>% 
  rename(long = location_long, lat = location_lat) %>% 
  select(-study_id, -declination, -inclination, -intensity) # remove all unneeded columns 

# add manually annotated migration data 
data <- data %>% left_join(mig %>% select(-tag_id, -Comments), by = join_by(bird_name, ymo)) %>% 
  mutate(migration = factor(migration))
data <- data %>% drop_na(migration) # remove all rows for which no information on migration exists (this is because the classification has been done manually in 2019 but we also use data from afterwards)

# create column which tells whether the current migration was the first autumn migration 
data <- data %>%
  group_by(bird_name) %>%
  mutate(
    first_migration = cumsum(migration == 'migratory' & lag(migration, default = 'stationary') != 'migratory'),
    first_migration = ifelse(migration == 'migratory', first_migration, NA)) %>%
  ungroup() 

data <- data %>% 
  group_by(bird_name) %>% 
  filter(migration == 'migratory') %>% 
  mutate(
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)), 1, first_migration), 
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)) + 1 & timestamp < as.POSIXct(paste0(year(min(timestamp, na.rm = TRUE)) + 1, '-07-15 00:00:00'), format = "%Y-%m-%d %H:%M:%S"), 2, first_migration)) %>%
  filter(first_migration %in% c(1,2)) %>% select(bird_name, timestamp, first_migration) %>% 
  right_join(data %>% select(-first_migration), by = join_by(bird_name, timestamp))

data <- data %>% mutate(first_migration = factor(case_when(first_migration == 1 ~ 'autumn', 
                                             first_migration == 2 ~ 'spring', 
                                            TRUE ~ as.character(first_migration))))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Create a raster for the migration route of EGVU's  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create binding box 
bbox <- terra::ext(range(data$long), range(data$lat))

# create raster with geographic CRS for Europe where grid cell size varies with location
rast <- rast(ext = bbox, resolution = 0.5, crs = "EPSG:4326") # not sure if maybe the world EPSG would be better, just replace by 'EPSG:4326'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Get magnetic field values for first spring migration in raster  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate mean autumn migration date a
data <- data %>%
  group_by(bird_name) %>%
  mutate(mean_spring = if_else(first_migration == 'spring', 
                               mean(timestamp[first_migration == 'spring'], na.rm = TRUE), 
                               as.POSIXct(NA)))

# create df to extract magnetic field values for autumn migration for every bird 
birds <- unique(data$bird_name) # get ids from all birds 
mean_spring <- data %>% filter(first_migration == 'spring') %>% group_by(bird_name) %>% summarise(mean_spring = mean(timestamp)) # get mean autumn migration dates 

df_rast <- as.data.frame(xyFromCell(rast, 1:ncell(rast))) # get the coordinated from raster
df_rast <- df_rast %>% rename(long = x, lat = y) # rename raster

df_rast <- data.frame(long = rep(df_rast$long, times = length(birds)), 
                      lat = rep(df_rast$lat, times = length(birds)), 
                      bird_name = rep(birds, each = nrow(df_rast)),
                      inclination = NA, 
                      declination = NA, 
                      intensity = NA)

df_rast <- df_rast %>% left_join(mean_spring, by = join_by(bird_name)) # add the mean_spring migration dates 
df_rast <- df_rast %>% drop_na(mean_spring) # get rid of all cases where no spring migration date could be calculates because of missing data 

# get magnetic field values for mean date of each birds first spring migration
df_rast$inclination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$inclination
df_rast$declination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$declination
df_rast$intensity <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$intensity

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Get magnetic field values for first autumn migration for all locations of each bird --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# filter for autumn migration
data_autumn <- data %>% 
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
  full_join(data_autumn %>% select(-long, -lat, -mean_spring), by = join_by(bird_name), suffix = c('_spring', '_autumn'), relationship = 'many-to-many') %>% 
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

# create a plot 
ggplot() + 
  geom_point(data = overlap, 
             mapping = aes(x = long, y = lat, color = 'Experienced Magentic Field in Autumn'), size = 1.5) +
  geom_path(data = data %>% filter(first_migration == 'spring'), 
            mapping = aes(x = long, y = lat, color = 'Spring Migration')) +
  coord_equal() + 
  facet_wrap(~bird_name) +
  labs(title = "Experienced Magnetic Field During Autumn Migraion and Spring Migration route for ten Egyptian Vulture",
    x = "Longitude",
    y = "Latitude") +
  scale_color_manual(name = 'Legend', values = c('Experienced Magentic Field in Autumn' = 'grey80', 'Spring Migration' = 'red')) +
  theme_bw() + 
  theme(legend.position = c(0.85,0.2), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/EGVU_magnetic_field_autumn_spring.png', height = 9, width = 12)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Plot data on a map background --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## convert data to spatial objects - only needed if tmap is used
str(overlap)
magField_sf<-overlap %>%
  st_as_sf(coords = c("long", "lat"), crs=4326) %>%
  st_transform(3857)


springMig_sf<-data %>% filter(first_migration == 'spring') %>%
  st_as_sf( coords = c("long", "lat"), crs=4326) %>%
  group_by(bird_name) %>%
  arrange(timestamp) %>%
  summarize() %>%
  # st_cast("MULTILINESTRING")%>% # this ensures, that no weird connections are made between points that don't belong together - this was the quickest way for troubleshooting and there are probably better solutions! 
  st_transform(3857)

# create sf data 
data_sf<-data %>%
  st_as_sf( coords = c("long", "lat"), crs=4326) %>%
  st_transform(3857) %>% arrange(bird_name, timestamp)


# extract background map tiles - check the provider argument, there are many options
basemap <- basemap_ggplot(ext=st_bbox(data_sf), map_service = "esri", map_type = "world_dark_gray_base")



# create a plot 
basemap +
  geom_sf(data = magField_sf, aes(color = 'Experienced Magnetic Field in Autumn'), size=1.5) +
  geom_sf(data = springMig_sf, mapping = aes(color = 'Spring Migration')) +
  coord_sf(default_crs = sf::st_crs(3857)) +
  facet_wrap(~bird_name) +
  labs(title = "Experienced Magnetic Field During Autumn Migration and Spring Migration route for ten Red Kites",
       x = "Longitude",
       y = "Latitude") +
  scale_color_manual(name = 'Legend', values = c('Experienced Magnetic Field in Autumn' = 'forestgreen', 'Spring Migration' = 'firebrick')) +
  theme(legend.position = c(0.85,0.12),  
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 15), 
        axis.title = element_text(size = 12))
ggsave(filename = 'output/EGVU_magnetic_field_autumn_spring.png', height = 9, width = 12)


