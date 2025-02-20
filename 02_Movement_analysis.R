#' ---
#' title : "Movement and Integrated Step-selection analysis for dogs"
#' subtitle: "Mudumalai Tiger Reserve, India"  
#' author: "Nilanjan Chatterjee"
#' date: "`r format(Sys.time(), '%d %B, %Y')`" 
#' output:
#'   html_document:
#'     toc: true
#' ---


#' ## Preamble
#' Load libraries
#+ warning=FALSE, message=FALSE
library(amt)
library(tidyverse)
library(here)
library(terra)
library(readxl)
library(sf)
library(geosphere)

#' Read the data
all_dog <- read.csv(here("Telemetry data", "MB_dog_move.csv"))
glimpse(all_dog)
all_dog$timestamp <-as.POSIXct(all_dog$timestamp, format = "%Y.%m.%d %H:%M:%OS",
                               tz = "Asia/Calcutta")

dat_summary<- read_excel("telemetry summary stats.xlsx")
glimpse(dat_summary)

#' Data wrangling
### Joined the village and some more columns with the movement data which
### will be used further for analysis
names(all_dog)[12] <- "Name" ###change the column name of the movement data
all_dog_joined <- full_join(all_dog, dat_summary[,c(1,7,9:11,20,21)], by  = "Name")
glimpse(all_dog_joined)

#' ########## Movement analysis
###create amt track data from the movement data 
track_dat <- all_dog %>% drop_na(timestamp) %>%
  arrange(Name, timestamp) %>%
  mk_track(location.long, location.lat, timestamp, 
           Name = Name,
           crs = 4326)  # CRS 4326 is the coordinate reference system

### I am looking for the fix interval rates across individual to see if there is 
### any large gap in the data 
mov_datsum <- track_dat %>% nest(data=-"Name") %>% 
  mutate(data = map(data, summarize_sampling_rate)) %>% 
  dplyr::select(Name , data) %>% 
  unnest(cols = data) %>% print(n =Inf)

#### create data for step-length first ----
tmp_amtdat <- track_dat %>% nest(data=- c("Name")) %>%
  mutate(data = map(data, ~ .x %>%
                      steps(lonlat = TRUE))) %>%
  unnest(cols = data) 

###Extract the hours from the data
tmp_amtdat <- tmp_amtdat |> mutate(hour = hour(t1_))

tmp_amtdat <- tmp_amtdat |> 
  mutate(daynight =case_when(hour >6 & hour < 18 ~ "day",
                             TRUE ~ "night"))

###Check the difference between day and night movement 
mov_sum <- tmp_amtdat |> group_by(Name, daynight) |>
  summarise(mean_speed = mean(sl_*60/as.numeric(dt_))) |>
  pivot_wider(names_from = daynight,
              values_from = mean_speed) 

joined_dat <- full_join(mov_sum, dat_summary[,-c(18,19)], by = "Name")
glimpse(joined_dat)

#### Check the results of movement, area and intensity of area use 
joined_dat |> filter( `Data count` != "Poor") |>
  ggplot(aes( x = village, y = intensity_use)) +
  geom_boxplot(fill = "skyblue")+
  labs(y= "Intensity of Use")+
  #ylim(0,35)+
  theme_bw()

joined_dat |> filter( `Data count` != "Poor") |>
  ggplot(aes( x = village, y = est)) +
  geom_boxplot(fill = "skyblue")+
  labs(y= "95% AKDE (sq.km)")+
  ylim(0,35)+
  theme_bw()

#' ### Distance calculation from home and centroid for all dogs

all_dog_joined <- all_dog_joined |> mutate(hour = hour(timestamp)) 

all_dog_joined <- all_dog_joined|> mutate(daynight =case_when(hour >6 & hour < 18 ~ "day",
                                           TRUE ~ "night")) |>
                                   mutate(home_dist =  distHaversine(cbind(location.long, location.lat),
                                                                     cbind(`Home long`,`Home lat`)),
                                          centroid_dist =  distHaversine(cbind(location.long, location.lat),
                                                                         cbind(Centroid_long,Centroid_lat)))


### Check the plot
all_dog_joined |> filter( `Data count` != "Poor") |>
  ggplot(aes(x = village, y = home_dist,fill = daynight)) + 
  geom_boxplot()+
  labs(y = "Distance from home") +
  theme_bw()
  
###################################################################################
#' ##### Step-selection analysis

### Load the habitat raster

hab_rast <- rast("Classified villages.tif")
plot(hab_rast)

sessionInfo()