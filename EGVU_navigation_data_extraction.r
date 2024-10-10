##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## EXTRACTING RED KITE DATA FOR NAVIGATION EXPLORATION #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### facilitate collaboration with Richard Holland on bird navigation
### Steffen Oppel developed a quantitative migration index - in 'REKI_weekly_migration_index.r'
### this script selects 10 example individuals (5 juv and 5 ad) that have the most migration weeks



####### LIBRARIES REQUIRED
library(tidyverse)
library(sf)
library(tmap)
library(move2)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(data.table); setDTthreads(percent = 65)
library(leaflet)
library(units)
library(oce)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP DOWNLOAD OF TRACKING DATA FROM MOVEBANK --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####### SPECIFY THE MOVEBANK ID OF THE STUDY FOR WHICH DATA ARE SHARED
MYSTUDY<-15869951
MYUSERNAME<-"Steffen"
movebank_store_credentials(username=MYUSERNAME, key_name = getOption("move2_movebank_key_name"), force = TRUE)
movebank_download_study_info(study_id=MYSTUDY[1])$sensor_type_ids

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DOWNLOAD MOVEBANK DATA AND ANIMAL INFO ----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

birds<-movebank_retrieve(study_id=MYSTUDY[1], entity_type="individual") %>%
  dplyr::rename(individual_id=id,bird_id=local_identifier) %>%
  dplyr::select(individual_id,bird_id,sex,latest_date_born) %>%
  filter(bird_id %in% c("Polya","Panteley", "Iliaz","Hedjet","Dobromir","Tatul","Volen","Sanie","Sava","Solomon")) %>%    ## not returned "Levkipos","Lola",
  arrange(bird_id) #%>%
  # mutate(start=c(ymd("2012-09-18"),ymd("2018-07-25"),ymd("2018-09-25"))) %>%
  # mutate(end=c(ymd("2012-10-08"),ymd("2019-01-27"),ymd("2018-12-29"))) 

## download gps data for multiple individuals
locs<-movebank_retrieve(study_id=MYSTUDY, entity_type="event",
                              sensor_type_id = "gps",
                              event_reduction_profile = "EURING_02",  ## this filters it to locations at least 50 km apart
                              individual_id = birds$individual_id
)


# ## filter first autumn migration
# #Iliaz: 19/09/2012	-07/10/2012
# #Panteley 2018-07-27T05:55:49.814610Z	2019-01-26T15:01:04.624955Z
# #Polya	2018-10-03T00:00:00Z	2018-12-28T18:00:00Z
# 
# 
# Iliaz<-locs %>%
#   rename(bird_id=individual_local_identifier) %>%
#   dplyr::select(bird_id,timestamp,location_lat,location_long) %>%
#   dplyr::filter(!is.na(timestamp)) %>%
#   dplyr::filter(bird_id=="Iliaz") %>%
#   #dplyr::filter(timestamp>birds$start[1] & timestamp<birds$end[1])   ### for first migration only
#   dplyr::filter(timestamp>birds$start[1] & timestamp<(birds$end[1]+years(4))) ### for more extended journeys
# Polya<-locs %>%
#   rename(bird_id=individual_local_identifier) %>%
#   dplyr::select(bird_id,timestamp,location_lat,location_long) %>%
#   dplyr::filter(!is.na(timestamp)) %>%
#   dplyr::filter(bird_id=="Polya") %>%
#   #dplyr::filter(timestamp>birds$start[3] & timestamp<birds$end[3])
#   dplyr::filter(timestamp>birds$start[3] & timestamp<(birds$end[3]+years(4)))
# Panteley<-locs %>%
#   rename(bird_id=individual_local_identifier) %>%
#   dplyr::select(bird_id,timestamp,location_lat,location_long) %>%
#   dplyr::filter(!is.na(timestamp)) %>%
#   dplyr::filter(bird_id=="Panteley") %>%
#   #dplyr::filter(timestamp>birds$start[2] & timestamp<birds$end[2])
#   dplyr::filter(timestamp>birds$start[2] & timestamp<(birds$end[2]+years(4)))
# 
# locs<-bind_rows(Iliaz, Polya,Panteley)

head(locs)
dim(locs)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SHOW TRACKS ON MAP 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_tracks<-locs %>%
  st_as_sf(coords=c('location_long','location_lat'), crs=4326)
tmap_mode("view")
c_osm <- tmaptools::read_osm(st_bbox(sample_tracks))
tm_basemap(server="OpenStreetMap") +
  tm_shape(sample_tracks)+
  tm_symbols(col = 'individual_local_identifier', size = 0.1)
tmap_mode("plot")
EGVUmap<-tm_shape(c_osm) +
  tm_rgb() +
  tm_shape(sample_tracks)+
  tm_symbols(col = 'individual_local_identifier', size = 0.1)
tmap_save(EGVUmap,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/EGVU_sample_tracks10.jpg")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT MAGNETIC FIELD FOR ALL TRACKS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


EGVU_mf<-magneticField(longitude=st_coordinates(sample_tracks)[,1],
                       latitude=st_coordinates(sample_tracks)[,2],
                       time=sample_tracks$timestamp,
                       version = 13)


sample_tracks<-sample_tracks %>%
  mutate(declination=EGVU_mf$declination,
         inclination=EGVU_mf$inclination,
         intensity=EGVU_mf$intensity)

locs<-locs %>%
  mutate(declination=EGVU_mf$declination,
         inclination=EGVU_mf$inclination,
         intensity=EGVU_mf$intensity)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE TRACKS AS CSV AND GPKG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


fwrite(locs,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/EGVU_sample_locations10.csv")
st_write(sample_tracks,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/EGVU_sample_locations10.gpkg", append=F)

