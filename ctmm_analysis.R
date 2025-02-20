library(ctmm)
library(tidyverse)
library(move)
library(amt)
library(move2)

setwd("C:/Users/nchatterjee/Work/Manuscript/Dog_population")
#dir()

indv1 <-read.csv("./Telemetry data/374 Mod.csv")
indv1$time <- as.POSIXct(indv1$Date.Time, format="%d-%m-%Y %H:%M", tz="UTC")
indv1$time <-indv1$time - (indv1$Time.to.fix) ##time to fix is in seconds
glimpse(indv1)
summary(indv1)
indv1 <- indv1 |> filter(Long >70) ##remove outliers with 0 latlong
indv1 <- indv1[order(indv1$time), ]
try_data <- move(x=indv1$Long, y=indv1$Lat, 
                  time=indv1$time, 
                  data=indv1, proj=CRS("+proj=longlat +ellps=WGS84"), 
                  animal="Burger", sensor="GPS")
#plot(try_data, type="b", pch=20)
dat_df <- as.data.frame(try_data)
dat_ctmm <- as.telemetry(dat_df)

###akde analysis
SVF1 <- variogram(dat_ctmm)
# Estimate an initial model
GUESS1 <- ctmm.guess(dat_ctmm, variogram=SVF1, interactive=FALSE)
# Select the best model
FIT1 <- ctmm.select(dat_ctmm, GUESS1, trace=2)
ud <- akde(dat_ctmm, FIT1, weights=TRUE)
summary(ud) ##wil give the HR size

###All individuals together
all_indv <-read.csv("./Telemetry data/MB_Dog_Move.csv")
all_indv$timestamp <- as.POSIXct(all_indv$timestamp, format="%Y.%m.%d %H:%M:%S", tz="UTC")
#indv1$time <-indv1$time - (indv1$Time.to.fix) ##time to fix is in seconds
glimpse(all_indv)
summary(all_indv)
#indv1 <- indv1 |> filter(Long >70) ##remove outliers with 0 latlong
#indv1 <- indv1[order(indv1$time), ]
try_data <- move(x=all_indv$location.long, y=all_indv$location.lat, 
                 time=all_indv$timestamp, 
                 data=all_indv, proj=CRS("+proj=longlat +ellps=WGS84"), 
                 animal=all_indv$individual.local.identifier, 
                 sensor=all_indv$sensor.type)
#plot(try_data, type="b", pch=20)
dat_df <- as.data.frame(try_data)
dat_ctmm <- as.telemetry(dat_df)
plot(dat_ctmm)

indv_1 <-variogram(dat_ctmm$Appu)


