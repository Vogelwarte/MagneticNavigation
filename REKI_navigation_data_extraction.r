##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## EXTRACTING RED KITE DATA FOR NAVIGATION EXPLORATION #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### facilitate collaboration with Richard Holland on bird navigation
### Steffen Oppel developed a quantitative migration index - in 'REKI_weekly_migration_index.r'
### this script selects 10 example individuals (5 juv and 5 ad) that have the most migration weeks

### updated script on 10 Oct to take 10 juveniles with outbound and return migration

### updated selection on 23 Dec 2024 to select 100 juveniles for a larger sample size

####### LIBRARIES REQUIRED
library(tidyverse)
library(sf)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(data.table); setDTthreads(percent = 65)
library(tmap)
library(oce)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA ON MIGRATION METRICS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load output of weekly migration index calculation
load("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_weekly_migration_metrics.RData")


### EXTRACT N OF MIGRATION WEEKS AND YEARS FOR EACH INDIVIDUAL
head(export)

IND_DATA<-export %>% separate(ywk, c("year", "week"), sep = "_", remove=FALSE) %>%
  mutate(year=as.numeric(year)) %>%
  mutate(MIG_WEEK=ifelse(MIG_PRED<0.5,0,1)) %>%
  group_by(id) %>%
  summarise(n_mig_weeks=sum(MIG_WEEK),first=min(year), last= max(year)) %>%
  mutate(range=last-first) %>%
  arrange(desc(n_mig_weeks), desc(range))
dim(IND_DATA)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ADD INFORMATION ON AGE AND SEX
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

reki_ind<-fread("C:/Users/sop/Dropbox/MASTER_databases_redkite/ind_life_histories/Individual_life_history_2015-2021.csv") %>%
  select(bird_id, tag_year,sex_compiled,age) %>%
  rename(id=bird_id)


## select individuals ##
SEL_INDS<-IND_DATA %>% left_join(reki_ind, by="id") %>%
  mutate(AD=ifelse(age=="3CY+",1,0)) %>%
  group_by(AD) %>%
  slice_head(n=100) %>%
  filter(!is.na(age)) #%>%
  #filter(AD==0)  ## use only juveniles to compare first autumn and first spring migration

if(!TRUE %in% (c(928,515) %in% SEL_INDS$id)) {
  SEL_INDS<-IND_DATA %>% left_join(reki_ind, by="id") %>%
    mutate(AD=ifelse(age=="3CY+",1,0)) %>%
    filter(id %in% c(928,515)) %>%
    bind_rows(SEL_INDS)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER TRACKING DATA FOR 10 INDIVIDUALS AND SHOW ON MAP 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_tracks<-REKI_sf %>%
  st_transform(4326) %>%
  filter(id %in% SEL_INDS$id) %>%
  mutate(id=as.numeric(id)) %>%
  left_join(export, by=c("id","ywk"))  %>%
  mutate(id=paste("REKI",id,sep="_")) %>%
  mutate(MIGRATION=ifelse(MIG_PRED<0.5,"not migrating","migrating")) %>%
  mutate(MIGRATION=ifelse(is.na(MIGRATION),"not migrating",MIGRATION)) %>%
  select(id, timestamp,geometry,distance_to_next,speed,MIGRATION)
  
### several weird duplicate records, which needs to be removed
duplicates<-export %>% group_by(id, ywk) %>% summarise(N=length(MOB)) %>% filter(N>1)
# sample_tracks[5266642,]
# export %>% filter(id==312) %>% filter(ywk=="2019_10")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SHOW TRACKS ON MAP 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmap_mode("view")
c_osm <- tmaptools::read_osm(st_bbox(sample_tracks))
tm_basemap(server="OpenStreetMap") +
  tm_shape(sample_tracks)+
  tm_symbols(col = 'id', size = 0.1)
tmap_mode("plot")
REKImap<-tm_shape(c_osm) +
  tm_rgb() +
  tm_shape(sample_tracks)+
  tm_symbols(col = 'id', size = 0.1)
tmap_save(REKImap,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_tracks100.jpg")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT MAGNETIC FIELD FOR ALL TRACKS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


REKI_mf<-magneticField(longitude=st_coordinates(sample_tracks)[,1],
                       latitude=st_coordinates(sample_tracks)[,2],
                       time=sample_tracks$timestamp,
                       version = 13)


sample_tracks<-sample_tracks %>%
  mutate(declination=REKI_mf$declination,
         inclination=REKI_mf$inclination,
         intensity=REKI_mf$intensity)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE TRACKS AS CSV AND GPKG or AS GOOGLE EARTH FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fwrite(sample_tracks,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_locations100.csv")
st_write(sample_tracks,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_locations100.gpkg", append=FALSE)
saveRDS(sample_tracks, file = "C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_tracks100.rds", version=3)

sample_tracks %>% st_write("data/REKI_sample_points.kml", append=FALSE)
sample_tracks %>% st_cast("LINESTRING") %>% st_write("data/REKI_sample_lines.kml", append=FALSE)


