library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(osmdata)
library(readxl)

#### dog UD results ----
dog_UD<-read_excel("C:/Users/nchatterjee/Work/Manuscript/Dog_population/telemetry summary stats.xlsx")
head(dog_UD)
glimpse(dog_UD)

dog_UD |> filter( `Data count` != "Poor") |>
  ggplot(aes( x = village, y = est)) +
  geom_boxplot(fill = "skyblue")+
  labs(y= "95% AKDE (sq.km)")+
  ylim(0,35)+
  theme_bw()

dog_UD |> 
  ggplot(aes( x = village, y = area)) + 
  geom_boxplot(fill = "skyblue")+
  theme_bw()

dog_UD |> filter( `Data count` != "Poor") |>
  ggplot(aes( x = village, y = intensity_use)) +
  geom_boxplot(fill = "skyblue")+
  labs(y= "Intensity of Use")+
  #ylim(0,35)+
  theme_bw()

###################################
### spatial files 
village_shp <-vect("C:/Users/nchatterjee/Work/Manuscript/Dog_population/Moyar and Masinagudi boundaries.kmz")

unzip("C:/Users/nchatterjee/Work/Manuscript/Dog_population/Moyar and Masinagudi boundaries.kmz", 
      exdir = "C:/Users/nchatterjee/Work/Manuscript/Dog_population/vill_kmz_extracted")  # Extract the file

# Identify the extracted KML file
kml_file <- list.files("C:/Users/nchatterjee/Work/Manuscript/Dog_population/vill_kmz_extracted", pattern = "\\.kml$", full.names = TRUE)

# Read the KML file into R
spatial_data <- st_read(kml_file)

# View the first few rows
head(spatial_data)

### download OSM data
bb <- (st_bbox(spatial_data))
temp_osm <- opq(bbox = bb) %>%
  add_osm_feature(key = 'highway') %>%
  add_osm_feature(key = 'name') %>%
  osmdata_sf()

lin_data <- temp_osm$osm_lines[, c("osm_id", "name", "highway", "geometry")]

ggplot()+
  geom_sf(data = lin_data)+
  geom_sf(data = spatial_data, fill = NA) + 
  coord_sf(ylim = c(11.56,11.62), xlim = c(76.62,76.71))+
  theme_bw()

villages <-rast("C:/Users/nchatterjee/Work/Manuscript/Dog_population/Classified villages.tif")
