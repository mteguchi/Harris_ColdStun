
rm(list=ls())
source("cold_stun_functions.R")
library(tidyverse)

save.fig <- FALSE #TRUE

dat0 <- read.csv(file = "data/ColdStun_data_July2019.csv", 
                 header = TRUE)  %>%
  rownames_to_column() %>%
  mutate(Weight_kg = Admit_weight_kg,
         CCL_cm = CCL,
         Sex = as.factor(toupper(Sex)),
         ID = rowname,
         Species = Species_Code, 
         Body_Temp_C = Initial_body_temp_C,
         Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>%
  select(ID, Latitude, Longitude, 
         Sex, Weight_kg, CCL_cm, 
         Body_Temp_C, Date, Species,
         Hypothermic) 

summary(dat0)

water.color <- "lightblue"
land.color <- "darkgray"
border.color <- "gray20"

E.end <- -112
W.end <- -128
N.end <- 51
S.end <- 32
# coast line data are from here: http://openstreetmapdata.com/data/coastlines or
# https://shoreline.noaa.gov/data/datasheets/medres.html and
# http://datapages.com/gis-map-publishing-program/gis-open-files/global-framework/global-heat-flow-database/shapefiles-list
# https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2016-eng.cfm

# https://www.ngdc.noaa.gov/mgg/shorelines/shorelines.html f is the full resolution - too big
# h is high resolution - should be good enough. 
all.coast <- sp::spTransform(rgdal::readOGR(dsn = "~/R/Oceans and Maps/gshhg-shp-2.3.7/GSHHS_shp/h",
                                            layer = "GSHHS_h_L1",
                                            verbose = FALSE),
                             sp::CRS("+proj=longlat +datum=WGS84"))

W.coast <- raster::crop(all.coast, raster::extent(c(W.end, E.end, S.end, N.end)))

W.coast.df <- broom::tidy(W.coast) 
#   filter(long < E.end & long > W.end & lat < N.end & lat > S.end)

# coast.line.df <- do.call(rbind, W.coast.df)
# colnames(coast.line.df) <- c('X', 'Y', 'idx')
# coast.line.Sp <- latlon2sp(coast.line.df, center.latlon)

# the following shape file was obttained from here: #https://gis.stackexchange.com/questions/153514/us-mexican-border-data
# note that the USGS link is broken but direct link "here" is available. 
# http://txpub.usgs.gov/BEHI/Data_download/Boundaries_Layers/International_Boundary_shp.zip
US_MX_border <- sp::spTransform(rgdal::readOGR(dsn = "~/R/Oceans and Maps/International_Boundary/shp",
                                               layer = "International_Boundary_Final",
                                               verbose = FALSE),
                                sp::CRS("+proj=longlat +datum=WGS84"))

US_MX_border <- raster::crop(US_MX_border, raster::extent(c(W.end, E.end, 30, 35)))
US_MX_border.df <- broom::tidy(US_MX_border)

# US-Canada border was obtained from here: 
# https://hifld-geoplatform.opendata.arcgis.com/datasets/canada-and-us-border
US_Canada_border <- sp::spTransform(rgdal::readOGR(dsn = "~/R/Oceans and Maps/Canada_and_US_Border",
                                                  layer = "Canada_and_US_Border",
                                                  verbose = FALSE),
                                   sp::CRS("+proj=longlat +datum=WGS84"))

US_Canada_border <- raster::crop(US_Canada_border, raster::extent(c(W.end, E.end, 47, 51)))
US_Canada_border.df <- broom::tidy(US_Canada_border) 

# state borders from here: https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
state_border <- sp::spTransform(rgdal::readOGR(dsn = "~/R/Oceans and Maps/cb_2017_us_state_500k",
                                                   layer = "cb_2017_us_state_500k",
                                                   verbose = FALSE),
                                    sp::CRS("+proj=longlat +datum=WGS84"))

state_border <- raster::crop(state_border, raster::extent(c(W.end, E.end, 30, 50)))
state_border.df <- broom::tidy(state_border)

p1 <- ggplot() + 
  geom_polygon(data = data.frame(y = c(N.end, S.end, S.end, N.end, N.end),
                                 x = c(W.end, W.end, E.end, E.end, W.end)),
               aes(x = x, y = y),
               fill = water.color,
               alpha = 0.8) + 
  geom_polygon(fill = land.color,
               data = W.coast.df,
               aes(x=long, y=lat, group = id),
               alpha = 0.9) +  
  geom_path(data = US_MX_border.df,
            aes(x = long, y = lat, group = group),
            color = border.color,
            size = 0.5) + 
  geom_path(data = US_Canada_border.df,
            aes(x = long, y = lat, group = group),
            color = border.color,
            size = 0.5) + #coord_map()
  geom_path(data = state_border.df,
            aes(x = long, y = lat, group = group),
            color = border.color,
            size = 0.5) + #coord_map()
  geom_path(data = data.frame(x = c(E.end, E.end), 
                              y = c(S.end, N.end)),
            aes(x = x, y = y),
            color = land.color,
            size = 1.2) +  
  geom_point(data = dat0,
             aes(x = Longitude, 
                 y = Latitude,
                 color = Species),
             size = 2,
             alpha = 0.5)  +
  coord_map() +
  xlim(c(-128, E.end))+
  ylab(expression(paste("Latitude (", degree, "N)"))) +
  xlab(expression(paste("Longitude (", degree, "W)", sep=""))) 
    
if (save.fig)
  ggsave(plot = p1,
         device = "png",
         dpi = 600,
         filename = "figures/stranding_map.png")
p1