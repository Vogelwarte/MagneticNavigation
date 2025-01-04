#~~
##### EGVU classification of migratory periods using a random forerst model--------  ##########################
#~~

## written by Filibert Heim, filibert.heim@posteo.de, in Dzember 2024
## goal is to annotate EGVU tracks automatically - migratory or stationary per month 
## weekly migration index conceived and developed by Steffen Oppel, Patrick Scherler, Florian Orgeret and the Red Kite Team at vogelwarte.ch


### load pre-prepared functions and models ----
setwd('C:/Users/filib/Documents/Praktika/Sempach/MagneticNavigation/')
source("C:/Users/filib/Documents/Praktika/Sempach/Tracking data/MoveApps_filter_functions.r") ### loads outlier and speed filter functions from MoveApps
source("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/DataPrep/MoveApps_filter_functions.r") ### loads outlier and speed filter functions from MoveApps


# load and install packages --- 
library(tidyverse)
library(sf)
library(amt)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(data.table); setDTthreads(percent = 65)
library(foreach)
library(doParallel)
library(randomForest)
library(ggpubr)
library(ranger)
library(move2)
library(units)
library(geosphere)
library(tmap)
library(readxl)
#install.packages('MazamaSpatialUtils')
#devtools::install_github('MazamaScience/PWFSLSmoke', build_vignettes=TRUE)
#install.packages('https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages')
library(PWFSLSmoke)  ## necessary for logger.info function in MoveApps filter functions
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename





#~~
#### Load data and make annotated migration data  -----------------------------------------------
#~~

# load data 
locs <- read.csv(file = 'data/EGVU_sample_locations_all.csv', sep = ',', dec = '.', header = T) # # load tracking data of EGVU which has been exported from Movebank by Steffen 
mig <- read_excel(path = 'data/Mig months.xlsx', sheet = 'Mig months') # this is a subset of annotated data, unfortunately only until 2019 

# make annotated data usable
mig[,c(1,3,17:18)] <- NULL
mig <- mig[grep(pattern = paste(c(paste0(unique(locs$tag_id), collapse = '|'),paste0(unique(locs$bird_name), collapse = '|')), collapse = '|'), x = mig$file_list),] # this filters for all bird_ids (names of the birds) but doesnt catch all of them because some were born in 2020 and the data is from 2019
mig <- mig %>%
  separate(file_list, into = c("bird_name", "tag_id", "year"), sep = "_", extra = "drop") %>%
  mutate(year = gsub("\\.png$", "", year)) %>% 
  pivot_longer(cols = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), names_to = 'month', values_to = 'migration') %>% 
  mutate(month = case_when(month == 'Jan' ~ '01', month == 'Feb' ~ '02', month == 'Mar' ~ '03', month == 'Apr' ~ '04', month == 'May' ~ '05', month == 'Jun' ~ '06', month == 'Jul' ~ '07', month == 'Aug' ~ '08', month == 'Sep' ~ '09',
                           month == 'Oct' ~ '10', month == 'Nov' ~ '11', month == 'Dec' ~ '12')) %>% 
  mutate(ymo = paste0(as.character(year), '_', month)) %>% 
  mutate(migration = case_when(migration == 0 ~ NA, migration == 1 ~ 'stationary', migration == 2 ~ 'migratory'), 
         ymo = paste0(year,'_', month))

# improve data structure of locs 
locs <- locs %>% mutate(timestamp = as_datetime(timestamp), 
                        bird_name = as.factor(individual_local_identifier), 
                        bird_name2 = as.factor(individual_local_identifier), # create dummy column 
                        tag_id = as.factor(tag_local_identifier), 
                        study_id = as.factor(study_id)) %>%
  rename(long = location_long, 
         lat = location_lat) %>% 
  select(-individual_local_identifier, -tag_local_identifier, -visible)

# make tracking data usable by creating a move2 object
locs_sf <- st_as_sf(locs, coords = c('long', 'lat'), crs = 4326, remove = F)
locs_move2 <- mt_as_move2(locs_sf, time_column = 'timestamp', track_id_column = 'bird_name2')


#~~
#### Data cleaning to remove outliers and duplicates----------------------------
#~~

# data time intervals are 2h - no reduction of resolution is needed

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
  select(bird_name,tag_id, timestamp, lat, long) %>% # select needed columns and remove unneeded bird_id's 
  filter(!is.na(bird_name))
head(locs_filt)
dim(locs_filt)

# get rid of all unneeded objects 
rm(locs)
rm(locs_move2)
gc()





#~~
#### Data preparation for migration index --------------------------------------
#~~

## to predict the migration index, we need several input metrics:
## iou = intensity of use = distance / MCP area to provide an area-standardized measure of the directionality of movements
## MOB = mobility metric, calculated from weekly diff(long)^2 + diff(lat)^2 (Scherler et al. 2024)
## dist = total travel distance per week in km
## area = area of the 95% MCP for a week in sq km
## nsd = net squared displacement for a week
## straight = straightness as calculated by direct distance start to end (of week) / total travel distance

# convert to sf for adjusting coordinate system and calculating distances
empty <- st_as_sfc("POINT(EMPTY)", crs = 4087) ## very important to choose the right CRS for equal distance/area projection: 3035 for Europe, 9834,9835, 54034, 4087 for the World - I do not know which is best!

EGVU_sf <- locs_filt %>% 
  mutate(ymo = paste(year(timestamp), format(timestamp, '%m'), sep = '_')) %>% # note that lubridate::month() does not return leading zeroes and will therefore cause problems in ordering
  mutate(ywk = paste(year(timestamp), format(timestamp,"%U"), sep = "_")) %>% # note that lubridate::isoweek() does not return leading zeroes and will therefore cause problems in ordering
  rename(id = bird_name)  %>%
  mutate(ymo_id = paste(id, ymo, sep = '_')) %>% # creates a new variable that consists of id and week of year
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>% 
  st_transform(4087) %>% ## very important to choose the right CRS for equal distance/area projection: 3035 for Europe, 9834,9835, 54034, 4087 for the World - I do not know which is best!
  group_by(id) %>% 
  mutate(
    elapsed_time = lead(timestamp) - timestamp,  ## time difference in seconds
    distance_to_next = sf::st_distance(geometry, lead(geometry, default = empty), by_element = TRUE),  ## distance in metres
    speed=as.vector(distance_to_next)/as.vector(elapsed_time))   ## speed will be in m/s

#### INSPECT THE CURIOUS CASES OF LONG DISTANCES PER WEEK
# outliers<-WHST_sf %>% 
#   st_drop_geometry() %>% # this removes the column geometry 
#   group_by(ywk_id) %>% 
#   summarise(weekdist = sum(distance_to_next)/1000) %>% # convert to kilometers 
#   filter(as.numeric(weekdist)>2000)
# weird<-WHST_sf %>% filter(ywk_id %in% outliers$ywk_id) %>% mutate(size=0.2)
# 
# #### LOOK AT THE DATA TO SPOT OBVIOUS OUTLIERS --------------------------------
# # create bounding box and extract tiles - check the provider argument, there are many options
# bbox <- st_sfc(st_point(c(-20, -35)), st_point(c(55, 56)), crs = 4326) %>% st_bbox()
# basemap <- maptiles::get_tiles(x = bbox, 
#                                zoom = 3,
#                                crop = TRUE, provider = "OpenTopoMap")
# tmap_mode("view")
# tm_shape(basemap)+
#   tm_rgb()+
#   tm_shape(weird %>% st_transform(4326))  +
#   tm_symbols(size="speed",col = "red")


#### CALCULATE DISTANCE PER WEEK AND INDIVIDUAL --------------------------------
dist <- EGVU_sf %>% 
  st_drop_geometry() %>% # this removes the column geometry 
  group_by(ymo_id) %>% 
  summarise(monthdist = as.numeric(sum(distance_to_next)/1000)) %>% # convert to kilometers and remove unit [m]
  filter(!is.na(monthdist)) # filter out all cases (mostly uncomplete months)


#### CALCULATE MOBILE PHASE INDEX PER WEEK AND INDIVIDUAL --------------------------------
# adapted from Scherler et al. 2024
# calculated in WGS84 geographic coordinates

MOB <- EGVU_sf %>% 
  ungroup() %>% 
  st_transform(4326) %>% 
  mutate(lat = sf::st_coordinates(.)[,2], long = sf::st_coordinates(.)[,1]) %>% 
  st_drop_geometry() %>% 
  group_by(ymo_id) %>% 
  summarise(monthlatrange = max(lat) - min(lat),  # calculate weekly latitudinal range traveled
            monthlongrange = max(long) - min(long)) %>% # calculate weekly longitudinal range traveled
  mutate(MOB = sqrt(monthlatrange^2 + monthlongrange^2)) %>% 
  select(ymo_id, MOB)



#### CALCULATE INTENSITY OF USE PER WEEK AND INDIVIDUAL --------------------------------
# based on Mumme et al. 2023 https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16769
# adapted from Almeida et al. (2010) Indices of movement behaviour: conceptual background, effects of scale and location errors. Zoologia (Curitiba), 27, 674-680.

# create an amt track because the amt package has lots of functions to calculate metrics like iou, nsd, area etc.
trackingdata <- EGVU_sf %>% 
  ungroup() %>% 
  mutate(lat = sf::st_coordinates(.)[,2], long = sf::st_coordinates(.)[,1]) %>% 
  st_drop_geometry()

track_amt <- mk_track(tbl=trackingdata, # this function creates a track 
                      .x = long, # provide coordinates
                      .y = lat,
                      .t = timestamp, # provide time information 
                      id = ymo_id, # additional column to group id_weeks
                      crs = 4087) # provide coordinate reference system - same as above

# add n locations and remove all weeks with <3 locations
dist <- track_amt %>% dplyr::count(id) %>%
  rename(ymo_id = id) %>%
  filter(n > 3) %>%   ## select any other threshold that is sensible for your data
  left_join(dist, by = "ymo_id")  %>%
  filter(!is.na(monthdist))

# start parallel loop over all individuals to calculate IOU Metric for every week of tracking data
n.cores <- 8  # this can potentially take a while, so run the function in parallel
registerDoParallel(n.cores)

iououtput <- 
  foreach(s = 1:dim(dist)[1],.combine=rbind, .packages=c('amt','tidyverse'),.inorder=FALSE,.errorhandling="remove",.verbose=FALSE) %dopar% {
    out<-dist[s,]
    mcp1 <- hr_mcp(track_amt %>% filter(id==out$ymo_id), levels = c(0.95,0.99,1))
    hr_area<-hr_area(mcp1)
    out$area<-hr_area$area[hr_area$level==1]/1000000
    out$iou<-as.numeric(intensity_use(track_amt %>% filter(id==out$ymo_id)))
    nsd <- nsd(track_amt %>% filter(id==out$ymo_id))
    out$nsd<-max(nsd)-min(nsd)  ### calculate range of nsd
    out$straight<-straightness(track_amt %>% filter(id==out$ymo_id))
    return(out)
  } ### end loop over all individuals

dim(iououtput)




#~~
#### COMBINE ALL DATA, CREATE TRAINING DATA SET --------------------------------------
#~~

# bring together data from MOB and IOU
EGVU_data <- iououtput %>% 
  rename(dist=monthdist) %>%
  left_join(MOB, by = "ymo_id") %>%
  left_join(EGVU_sf %>%
              st_drop_geometry() %>%
              ungroup() %>%
              group_by(ymo_id) %>%
              summarise(n_locs=length(lat), lat=mean(lat, na.rm=T), long=mean(long, na.rm=T)) %>%
              filter(ymo_id %in% iououtput$ymo_id), by='ymo_id') %>%
  separate_wider_delim(cols=ymo_id,delim="_",names=c("bird_name","year","month")) %>% 
  mutate(date = ymd(paste(year, month, '01', spe = '-')))

# add manually annotated data 
EGVU <- EGVU_data %>% left_join(mig %>% select(bird_name, year, month, migration), by = join_by(bird_name, year, month)) 

# select all annotated data as training data 
EGVU_train <- EGVU %>% filter(!is.na(migration)) # this is roughly 1/3 of the data 


#~~
#### TRAIN MODEL  --------------------------------------
#~~

# check model requirements: balanced categories in data set 
EGVU_train %>% group_by(migration) %>% summarise(n = length(migration)) # teh ratio is 116/190 migration and stationary 


# find optimal number of trees 
tree <- data.frame(num.trees = c(seq(3, 500, 10), seq(500, 5000, 100)), oob = NA)
for(i in 1:nrow(tree)){
  tree$oob[i] <- ranger(formula = as.factor(migration)~dist+area+iou+MOB+nsd+straight, data = EGVU_train, probability = T, num.trees = tree$num.trees[i], mtry = 2, importance = 'permutation')$prediction.error
}
plot(tree$num.trees, tree$oob, type = 'line', col = 'red') # seems like everything is okay at num.trees = 500

# find optimal number of variables to split each node, default and recommendation is sqrt(number.of.predictors) = 2.5, but try values around 2
mtry <- data.frame(mtry = 1:6, oob = NA)
for(i in 1:nrow(tree)){
  mtry$oob[i] <- ranger(formula = as.factor(migration)~dist+area+iou+MOB+nsd+straight, data = EGVU_train, probability = T, num.trees = 500, mtry = mtry$mtry[i], importance = 'permutation')$prediction.error
}
plot(mtry$mtry, mtry$oob, type = 'line', col = 'red') # continue with 2


# fit final model 
rf <- ranger(formula = as.factor(migration)~dist+area+iou+MOB+nsd+straight, data = EGVU_train, probability = T, num.trees = 500, mtry = 2, importance = 'permutation')
rf$prediction.error

# predict rf probability to all training data for accuracy assessment
EGVU_train$migration_pred <- predict(rf, data=EGVU_train, type="response", oob=T)$predictions[,1]
EGVU_train <- EGVU_train %>%
  mutate(dist = as.numeric(dist)) %>%
  mutate(migration_class = ifelse(migration_pred < 0.5,'stationary','migratory')) %>% 
  mutate(class = ifelse(migration_class == migration,"correct","false")) %>% ### this is for RF predictions
  mutate(evaluation = case_when(migration_class == "migratory" & class == 'correct' ~ 'migratory_correct', 
                               migration_class == 'migratory' & class == 'false' ~ 'migratory_false', 
                               migration_class == 'stationary' & class == 'false' ~ 'stationary_false', 
                               migration_class == 'stationary' & class == 'correct' ~ 'stationary_correct'))

# plot predictions against variables
mobpred<-ggplot(EGVU_train,aes(x=MOB, y=migration_pred, colour=evaluation)) + geom_point(size = 2)+theme(legend.position=c(0.85,0.2)) + theme_bw()
distpred<-ggplot(EGVU_train,aes(x=as.numeric(dist), y=migration_pred, colour=evaluation)) + geom_point(size = 2) +guides(colour = FALSE) + theme_bw()
nsdpred<-ggplot(EGVU_train,aes(x=nsd, y=migration_pred, colour=evaluation)) + geom_point(size = 2)+guides(colour = FALSE) + theme_bw()
areapred<-ggplot(EGVU_train,aes(x=area, y=migration_pred, colour=evaluation)) + geom_point(size = 2)+guides(colour = FALSE) + theme_bw()
ggarrange(mobpred, distpred, areapred, nsdpred,ncol = 2, nrow = 2, legend = 'right', common.legend = T)

# filter for mis-classifications to understand erroneous predictions 
EGVU_train %>% filter(class == 'false') # only  misclassifiacations, I think thats okay

#~~
#### APPLY RANDOM FOREST MODEL TO CLASSIFY MIGRATORY PERIODS PER MONTH TO WHOLE EGVU DATA   --------------------------------------
#~~

# predict for complete data set 
EGVU$migration_pred <- predict(rf, data = EGVU, type = 'response', oob = T)$predictions[,1]


#~~
#### FORMAT DATA, CREATE FLAG FOR MIS-CLASSIFICATION AND PLOT RESULTS --------------------------------------
#~~

# create flags to show mis-classifications from manually annotated data (same step as above)
EGVU <- EGVU %>%
  mutate(dist = as.numeric(dist)) %>%
  mutate(migration_class = ifelse(migration_pred < 0.5,'stationary','migratory')) %>% 
  mutate(class = ifelse(migration_class == migration,"correct","false")) %>% ### this is for RF predictions
  mutate(evaluation = case_when(migration_class == "migratory" & class == 'correct' ~ 'migratory_correct', 
                                migration_class == 'migratory' & class == 'false' ~ 'migratory_false', 
                                migration_class == 'stationary' & class == 'false' ~ 'stationary_false', 
                                migration_class == 'stationary' & class == 'correct' ~ 'stationary_correct'))


# plot predictions against variables
mobpred<-ggplot(EGVU,aes(x=MOB, y=migration_pred, colour=evaluation)) + geom_point(size = 2)+theme(legend.position=c(0.85,0.2)) + theme_bw()
distpred<-ggplot(EGVU,aes(x=as.numeric(dist), y=migration_pred, colour=evaluation)) + geom_point(size = 2) +guides(colour = FALSE) + theme_bw()
nsdpred<-ggplot(EGVU,aes(x=nsd, y=migration_pred, colour=evaluation)) + geom_point(size = 2)+guides(colour = FALSE) + theme_bw()
areapred<-ggplot(EGVU,aes(x=area, y=migration_pred, colour=evaluation)) + geom_point(size = 2)+guides(colour = FALSE) + theme_bw()
ggarrange(mobpred, distpred, areapred, nsdpred,ncol = 2, nrow = 2, legend = 'right', common.legend = T)


# plot predicted behaviour against time
EGVU %>%
  ggplot() +
  geom_line(mapping = aes(x = date, y = migration_pred, col = migration_pred)) +
  scale_color_viridis_c(alpha=1,begin=0,end=1,direction=1) +
  facet_wrap(~bird_name, scales="free_x") +
  theme(legend.position="none")
# ggsave("WHST_mig_pred_over_time.jpg", width=20, height=10)


#### MAP PREDICTED BEHAVIOUR --------------------------------
# create bounding box and extract tiles - check the provider argument, there are many options
bbox <- st_sfc(st_point(c(-20, -35)), st_point(c(55, 56)), crs = 4326) %>% st_bbox()
basemap <- maptiles::get_tiles(x = bbox,
                               zoom = 1,
                               crop = TRUE, provider = "OpenTopoMap")
tmap_mode("plot")
outmap<-tm_shape(basemap)+
  tm_rgb()+
  tm_shape(EGVU %>% st_as_sf(coords = c('long', 'lat'), crs = 4326))  +
  tm_symbols(size=0.3,col = "migration_pred", palette="viridis")
outmap
# tmap_save(outmap,"WHST_pred_mig_map.jpg")



#~~
#### EXPORT CLASSIFIED DATA  --------------------------------------
#~~

# bring data in a rather similar structure as the predictions are in other classifications for REKI and WHST
EGVU_export <- EGVU %>% select(-migration_class, -class, -evaluation) %>% rename(MIG_PRED = migration_pred, MIG_MAN = migration)

# export data 
fwrite(EGVU_export, "data/EGVU_mig_predictions.csv")



