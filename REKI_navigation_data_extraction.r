##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## EXTRACTING RED KITE DATA FOR NAVIGATION EXPLORATION #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### facilitate collaboration with Richard Holland on bird navigation
### Steffen Oppel developed a quantitative migration index - in 'REKI_weekly_migration_index.r'
### this script selects 10 example individuals (5 juv and 5 ad) that have the most migration weeks

### updated script on 10 Oct to take 10 juveniles with outbound and return migration


####### LIBRARIES REQUIRED
library(tidyverse)
library(sf)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(data.table); setDTthreads(percent = 65)
library(tmap)


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
  slice_head(n=10) %>%
  filter(!is.na(age)) %>%
  filter(AD==0)  ## use only juveniles to compare first autumn and first spring migration



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
  




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SHOW TRACKS ON MAP 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmap_mode("view")
tm_basemap(server="OpenStreetMap") +
  tm_shape(sample_tracks)+
  tm_symbols(col = 'id', size = 0.1)
tmap_mode("plot")
REKImap<-tm_basemap(server="OpenStreetMap") +
  tm_shape(sample_tracks)+
  tm_symbols(col = 'id', size = 0.1)
tmap_save(REKImap,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_tracks10.jpg")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE TRACKS AS CSV AND GPKG or AS GOOGLE EARTH FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fwrite(sample_tracks,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_locations10.csv")
st_write(sample_tracks,"C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_locations10.gpkg")
saveRDS(sample_tracks, file = "C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/Navigation/data/REKI_sample_tracks.rds", version=3)

sample_tracks %>% st_write("data/REKI_sample_points.kml", append=FALSE)
sample_tracks %>% st_cast("LINESTRING") %>% st_write("data/REKI_sample_lines.kml", append=FALSE)


