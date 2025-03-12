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
library(glmmTMB)
library(mixedSSA)
library(msm)

#' ## Read the data
all_dog <- read.csv(here("Telemetry data", "MB_dog_move.csv"))
glimpse(all_dog)
all_dog$timestamp <-as.POSIXct(all_dog$timestamp, format = "%Y.%m.%d %H:%M:%OS",
                               tz = "Asia/Calcutta")

dat_summary<- read_excel("telemetry summary stats.xlsx")
glimpse(dat_summary)

#' #### Data wrangling
### Joined the village and some more columns with the movement data which
### will be used further for analysis
names(all_dog)[12] <- "Name" ###change the column name of the movement data
all_dog_joined <- full_join(all_dog, dat_summary[,c(1,7,9:11,20,21)], by  = "Name")
glimpse(all_dog_joined)

#' #### Movement analysis
### create amt track data from the movement data 
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

all_dog_joined <- all_dog_joined|> 
  mutate(daynight =case_when(hour >6 & hour < 18 ~ "day",
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
  
#' ## 2. Step-selection analysis
### Load the habitat raster

hab_rast <- rast("Classified villages.tif")
plot(hab_rast)
names(hab_rast) <-"habitat"

###create spatial data and change the latlong projection to UTM projection
### also dropped unnecessary columns for further analysis
all_dog_sf <- all_dog_joined |> filter(gps.hdop <10) |>
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) |>
  st_transform(st_crs(hab_rast)) |>
  select(!c(event.id, visible, algorithm.marked.outlier,
            manually.marked.outlier, sensor.type, individual.taxon.canonical.name,
            tag.local.identifier,study.name))

### extract the coordinates to create the tracks and steps dataset
all_dog_sf <- all_dog_sf |>
  mutate(easting = st_coordinates(all_dog_sf)[,1],
         northing = st_coordinates(all_dog_sf)[,2])
glimpse(all_dog_sf)

### Create track data using amt 
dogtrack <- mk_track(all_dog_sf, .x=easting, .y=northing, .t=timestamp, crs=32643, 
            Ind_id = Name,
            village = village)

### Create steps data for ISSF analysis
dog_steps <- dogtrack %>% nest(data=-c("Ind_id", "village")) %>%
  mutate(data = map(data, ~ .x %>%
                      track_resample(rate = minutes(5), tolerance = minutes(1)) %>%
                      filter_min_n_burst(min_n = 3) %>% # drop bursts with less than 3 observations
                      steps_by_burst())) %>%
  unnest(cols=data) %>% # by using unnest before calling random steps, all 
  # random steps will use the same tentative movement kernel
  random_steps(n_control =20) %>%
  extract_covariates(hab_rast, where="both") %>%
  filter(!is.na(ta_)) %>% # remove steps without a turn angle  
  mutate(log_sl_=log(sl_), cos_ta_ = cos(ta_)) 
dog_steps$step_id <-paste0(dog_steps$Ind_id, "_", dog_steps$step_id_)
glimpse(dog_steps)
mosaic::tally(~ habitat_start + case_, data= dog_steps)
saveRDS(dog_steps, "dog_issa_data.rds") ###save the data for further use

#' ### 2.1 Dummy coding for categorical variables
### There are six classes, we want to merge scrubland and open fallow
### to minimize the number of classes
### We also dropped step-lengths less than 500 meter 

dog_steps$habitat_end[dog_steps$habitat_end>10]<-10
dog_steps$habitat_start[dog_steps$habitat_start>10]<-10

dog_steps1 <- fastDummies::dummy_cols(dog_steps, 
                                      select_columns = c("habitat_start", "habitat_end"))
names(dog_steps1)[20:29] <-c("Water_start", "Forest_start", "Agriculture_start", "Settlement_start",
                            "Scrubland_start","Water_end", "Forest_end", "Agriculture_end", "Settlement_end",
                            "Scrubland_end")
dog_steps1 <- dog_steps1 |> filter( sl_ < 500)

#' ### 2.1 Fitting Integrated Step-selection function (ISSF)
#' We are only using the categorical habitat covariate to understand 
#' movement of the dog in the different habitats   
## ---- eval=FALSE----------------------------------------------------------------------------------

### We will use the glmmTMB and individuals and villages as random effects 
### We want to fit models where the habitat use and movement vary across the 
### individuals and village 
### Model 1 

### Final Model  
ssf_hab2 <- glmmTMB(case_ ~  (sl_)+log_sl_ + cos_ta_+ Forest_end + Agriculture_end + Settlement_end +  Scrubland_end + 
                      Forest_start:sl_ + Settlement_start:sl_ + Forest_start:log_sl_ + Settlement_start:log_sl_ +
                      Agriculture_start:sl_ + Scrubland_start:sl_ + Agriculture_start:log_sl_ + Scrubland_start:log_sl_ +
                      Forest_start:cos_ta_ + Settlement_start:cos_ta_ + Agriculture_start:cos_ta_ + Scrubland_start:cos_ta_ +
                     (1|step_id)+ ((0+ Forest_end + Agriculture_end + Settlement_end +  Scrubland_end)| village)+ 
                      (0 + sl_ +  log_sl_ + cos(ta_) |village) -1,
                              family=poisson(), data = dog_steps1, doFit=FALSE)

###Model 3
ssf_hab <- glmmTMB(case_ ~  (sl_)+log_sl_ + cos_ta_+ as.factor(habitat_end) + 
                          (1|step_id)+ (0+ as.factor(habitat_end)| village)+ 
                          as.factor(habitat_start):sl_+ as.factor(habitat_start):log_sl_+
                          as.factor(habitat_start):cos_ta_ +                   
                          (0 + sl_+log_sl_+cos_ta_|village)-1,
                        family=poisson(), data = dog_steps, doFit=FALSE)

ssf_hab2$parameters$theta[1] <- log(1e3)
nvarparm<-length(ssf_hab2$parameters$theta)
ssf_hab2$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
ssf.hab_mod6 <- glmmTMB:::fitTMB(ssf_hab2)
#saveRDS(ssf.hab_mod, "dog_ssf_wvillrf_movparam.rds")

## -------------------------------------------------------------------
#' Model Summary
#' ###  Calculating the updated movement parameters with the mixedSSA package
#' Read the output from the fitted model
ssf.hab_mod6 <- readRDS("Updated_dog_ssf_mod.rds") 
summary(ssf.hab_mod6)
modelsummary::modelplot(ssf.hab_mod6, coef_omit = c("village|SD"))
#ggsave("Coef_plot.jpeg", width = 9, height = 6, units = "in", dpi=300)

### We are calculating the tentative distribution for step-lengths here from the observed data 
sl_dist <- amt::fit_distr(ssf.hab_mod6$frame$sl_[ssf.hab_mod6$frame$case_ == "TRUE"], "gamma",na.rm = T)

### The package mixedSSA can be used for both categorical and continuous data
### We are making sure that the variable of interest is categorical in our case 
ssf.hab_mod6$frame <- ssf.hab_mod6$frame %>% mutate(Agriculture_start = as.factor(Agriculture_start),
                                                    Settlement_start = as.factor(Settlement_start),
                                                    Forest_start = as.factor(Forest_start),
                                                    Scrubland_start = as.factor(Scrubland_start))

#' ### 2.1 Update the movement distribution for different habitat
#' Habitat Categorical covariate 
Agri_mov <- update_dist(model = ssf.hab_mod6,
            dist_name = "gamma",
            beta_sl = "sl_", # the name of the step lengths coefficient in our model
            beta_log_sl = "log_sl_", # the name of the log(sl) coefficient in our model
            interaction_var_name = "Agriculture_start", #or newdist
            random_effects_var_name = "village", # for individual update
            #quantiles = quantiles, # for continuous interaction variables
            tentative_dist = sl_dist)@updated_parameters
Forest_mov <- update_dist(model = ssf.hab_mod6,
                          dist_name = "gamma",
                          beta_sl = "sl_", # the name of the step lengths coefficient in our model
                          beta_log_sl = "log_sl_", # the name of the log(sl) coefficient in our model
                          interaction_var_name = "Forest_start", #or newdist
                          random_effects_var_name = "village", # for individual update
                          #quantiles = quantiles, # for continuous interaction variables
                          tentative_dist = sl_dist)@updated_parameters

Scrubland_mov <- update_dist(model = ssf.hab_mod6,
                        dist_name = "gamma",
                        beta_sl = "sl_", # the name of the step lengths coefficient in our model
                        beta_log_sl = "log_sl_", # the name of the log(sl) coefficient in our model
                        interaction_var_name = "Scrubland_start", #or newdist
                        random_effects_var_name = "village", # for individual update
                        #quantiles = quantiles, # for continuous interaction variables
                        tentative_dist = sl_dist)@updated_parameters
Settlement_mov <- update_dist(model = ssf.hab_mod6,
                          dist_name = "gamma",
                          beta_sl = "sl_", # the name of the step lengths coefficient in our model
                          beta_log_sl = "log_sl_", # the name of the log(sl) coefficient in our model
                          interaction_var_name = "Settlement_start", #or newdist
                          random_effects_var_name = "village", # for individual update
                          #quantiles = quantiles, # for continuous interaction variables
                          tentative_dist = sl_dist)@updated_parameters

all_hab <- rbind(Agri_mov[4:5,], Forest_mov[4:5,], Settlement_mov[4:5,], Scrubland_mov[4:5,])
all_hab$habitat <- rep(c("Agriculture", "Forest", "Settlement", "Scrubland"), each =2)
all_hab <- all_hab |> mutate(mean = scale*shape) ###Gamma distribution 

#' ### 2.2 Delta method to estimate SE for habitat variable 
#'  Categorical Habitat covariate 
estmean <- summary(ssf.hab_mod6)$coefficients$cond[, 1]
estvar <- vcov(ssf.hab_mod6)$cond

sl_df <-data.frame()
sl_df[1,1:2] <- sl_dist$params
s11 <- sl_df[1, 1]
s12 <- sl_df[1, 2]

## delta method
#sl_df$se[1] <- sqrt(deltamethod(~ (s11 ) * (1 / (1 / s12 )), estmean, estvar))
#all_hab$se[2] <- sqrt(deltamethod(~ (s11 +  x2) * (1 / (1 / s12 -  x1)), estmean, estvar))
all_hab$se[1] <- sqrt(deltamethod(~ (s11 +  x14 + x2) * (1 / (1 / s12 - ( x12 + x1))), estmean, estvar))
all_hab$se[2] <- sqrt(deltamethod(~ (s11 +  x14 + x2) * (1 / (1 / s12 - ( x12 + x1))), estmean, estvar))
all_hab$se[3] <- sqrt(deltamethod(~ (s11 +  x10 + x2) * (1 / (1 / s12 - ( x8 + x1))), estmean, estvar))
all_hab$se[4] <- sqrt(deltamethod(~ (s11 +  x10 + x2) * (1 / (1 / s12 - ( x8 + x1))), estmean, estvar))
all_hab$se[5] <- sqrt(deltamethod(~ (s11 +  x11 + x2) * (1 / (1 / s12 - ( x9 + x1))), estmean, estvar))
all_hab$se[6] <- sqrt(deltamethod(~ (s11 +  x11 + x2) * (1 / (1 / s12 - ( x9 + x1))), estmean, estvar))
all_hab$se[7] <- sqrt(deltamethod(~ (s11 +  x15 + x2) * (1 / (1 / s12 - ( x13 + x1))), estmean, estvar))
all_hab$se[8] <- sqrt(deltamethod(~ (s11 +  x15 + x2) * (1 / (1 / s12 - ( x13 + x1))), estmean, estvar))

###check the final output
all_hab

### Plot the results
all_hab |> #filter(habitat != "tentative")|>
  arrange(mean)|>
  ggplot() +
  geom_point(aes(x = habitat, y = mean, col = village), size = 3 , position = position_dodge(0.2)) +
  geom_linerange(aes(x = habitat, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), 
                 linewidth = 1, position = position_dodge2(0.2)) +
  scale_x_discrete(limit = c("Settlement",  "Scrubland" , "Agriculture" ,
                              "Forest"))+
  labs(x = "Habitat class", y = "Mean Step-length (m)", col = "Village") +
  theme_bw() +
  theme(legend.position = c(0.1,0.8),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#ggsave("Village_habitat_movement.jpeg", width=7.5, height = 5, units = "in", dpi = 300)

sessionInfo()