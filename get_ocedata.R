# for downloading SST etc files from CoastWatch

rm(list=ls())
source("cold_stun_functions.R")
library(tidyverse)
library(ncdf4)


save.fig <- FALSE #TRUE

dat0 <- read.csv(file = "data/ColdStun_data_June2019.csv", 
                 header = TRUE)  %>%
  rownames_to_column() %>%
  select(-c(Field_ID, Database_ID, NMFS_ID, Other_ID, Location, City, County, Country)) %>%
  mutate(Age = as.factor(toupper(Age)),
         Sex = as.factor(toupper(Sex)),
         Condition = as.factor(toupper(Condition)),
         ID = rowname,
         Species = as.factor(ifelse(Common_Name == "Green sea turtle", "CM",
                                    ifelse(Common_Name == "Loggerhead", "CC", "LO"))),
         Body_Temp_C = Confirmed.body.temp.C,
         Year = Year_Initially_Observed,
         Month = Month_Initially_Observed,
         Day = Day_Initially_Observed,
         Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>%
  select(ID, Age, State, Latitude, Longitude, 
         Sex, Weight_kg, CCL_cm, CCW_cm,
         Body_Temp_C, Date, Species) %>%
  filter(Latitude > 34.45)

ds <- c(3, 5, 10)
d.n <- length(ds)   # this should match the number of dx's below. 

# Get SST for 1, 5, 10 km around each point. 
latlon <- select(dat0, Latitude, Longitude) %>% 
  transmute(X = Longitude, Y = Latitude,
            BeginX_1 = NA, BeginY_1 = NA, EndX_1 = NA, EndY_1 = NA,
            BeginX_2 = NA, BeginY_2 = NA, EndX_2 = NA, EndY_2 = NA,
            BeginX_3 = NA, BeginY_3 = NA, EndX_3 = NA, EndY_3 = NA) 

UTM <- data.frame(X = NA, Y = NA,
                  BeginX_1 = NA, BeginY_1 = NA, EndX_1 = NA, EndY_1 = NA,
                  BeginX_2 = NA, BeginY_2 = NA, EndX_2 = NA, EndY_2 = NA,
                  BeginX_3 = NA, BeginY_3 = NA, EndX_3 = NA, EndY_3 = NA,
                  zone = NA)

# because the range of latitude is quite large, approximate center is moved along each data point.
# convert all lat/lon pairs to UTM coordinates in the appropriate zone. 
centers.UTM <- vector(mode = "list", length = nrow(latlon))
k <- 1
for (k in 1:nrow(dat0)){
  approx.center <- data.frame(X = -120, Y = dat0[k, "Latitude"])
  sp::coordinates(approx.center) <- c("X", "Y")
  sp::proj4string(approx.center) <- CRS("+proj=longlat +datum=WGS84")
  
  # North of ~ Pt Conception is zone 10, and south of it is zone 11
  zone <- ifelse(dat0[k, "Latitude"] > 34.46, 10, 11)
  
  centers.UTM[[k]] <- spTransform(approx.center,
                                  CRS(paste0("+proj=utm +zone=", zone, " ellps=WGS84")))
  
  tmp <- latlon2sp(latlon[k, c("X", "Y")], 
                   center.UTM = centers.UTM[[k]], zone )
  
  UTM[k, ] <- c(tmp@data$newX, tmp@data$newY,
                tmp@data$newX - ds[1], tmp@data$newY - ds[1], 
                tmp@data$newX + ds[1], tmp@data$newY + ds[1],
                tmp@data$newX - ds[2], tmp@data$newY - ds[2], 
                tmp@data$newX + ds[2], tmp@data$newY + ds[2],
                tmp@data$newX - ds[3], tmp@data$newY - ds[3], 
                tmp@data$newX + ds[3], tmp@data$newY + ds[3],
                zone)
  
}

# This loop should go over all dx's (2 each) - Begin and End (X,Y)
k <- k1 <- 1
for (k in 1:(2*d.n)){
  tmp <- data.frame(newX = UTM[, (2*k+1)],
                    newY = UTM[, (2*k+2)])
  for (k1 in 1:nrow(tmp)){
    latlon[k1, (2*k+1):(2*k+2)] <- sp2latlon(tmp[k1,], 
                                             center.UTM = centers.UTM[[k1]], 
                                             zone = UTM[k1, "zone"])@coords[, c("X", "Y")]
    
  }
  
}

latlon$Date <- dat0$Date

k <- k1 <- 1

# for all rows
for (k in 1:nrow(latlon)){
  # wind comes in 0.25 degree resolutions so we just pick one without worrying about other spatial
  # resolutions
  out.file.name.wind <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                               "wind.nc")
  
  out.file.name.wind.30d <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                   "wind_30d.nc")
  
  # for number of distances
  for (k1 in 1:d.n){
    out.file.name <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                            ds[k1], "km_sst.nc")
    out.file.name.anom <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                 ds[k1], "km_sstanom.nc")
    out.file.name.lag30d<- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                  ds[k1], "km_sst_lag30d.nc")
    out.file.name.anom.lag30d <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                        ds[k1], "km_sstanom_lag30d.nc")
    out.file.name.air.temp<- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                    ds[k1], "km_airtemp.nc")
    out.file.name.air.temp.lag30d <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                            ds[k1], "km_airtemp_lag30d.nc")
    out.file.name.sst.0125 <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                     ds[k1], "km_sst0125.nc")
    out.file.name.sst.0125.lag30d <- paste0("data/ncfiles/turtle_", dat0[k, "ID"], "_",
                                            ds[k1], "km_sst0125_lag30d.nc")
    
    
    xlim <- c(latlon[k, (4*k1 - 1)],
              latlon[k, (4*k1 + 1)])
    ylim <- c(latlon[k, (4*k1)],
              latlon[k, (4*k1 + 2)])
    #tlim <- c(latlon[k, "Date"] - 7,
    #          latlon[k, "Date"] + 7)
    
    if(!file.exists(out.file.name.sst.0125)){
      # SST, Aqua MODIS, NPP, 0.0125 degrees, West US, Day time (11 microns), 
      # 2002-present (14 day composite)
      # https://upwell.pfeg.noaa.gov/erddap/griddap/erdMWsstd14day_LonPM180.html
      sst0125URL <- paste0("https://upwell.pfeg.noaa.gov/erddap/griddap/erdMWsstd14day_LonPM180.nc?sst[(",
                           latlon[k, "Date"], "T00:00:00Z):1:(", 
                           latlon[k, "Date"], "T00:00:00Z)][(0.0):1:(0.0)][(", 
                           signif(ylim[1], digits = 6), "):1:(", 
                           signif(ylim[2], digits = 6), ")][(", 
                           signif(xlim[1], digits = 7), "):1:(", 
                           signif(xlim[2], digits = 7), ")]")
      
      test <- download.file(sst0125URL,
                            destfile = out.file.name.sst.0125,
                            mode='wb')
    }
    
    if(!file.exists(out.file.name.sst.0125.lag30d)){
      # SST, Aqua MODIS, NPP, 0.0125 degrees, West US, Day time (11 microns), 
      # 2002-present (14 day composite)
      # https://upwell.pfeg.noaa.gov/erddap/griddap/erdMWsstd14day_LonPM180.html
      sst0125lag30dURL <- paste0("https://upwell.pfeg.noaa.gov/erddap/griddap/erdMWsstd14day_LonPM180.nc?sst[(",
                                 (latlon[k, "Date"]-29), "T00:00:00Z):1:(", 
                                 (latlon[k, "Date"]-29), "T00:00:00Z)][(0.0):1:(0.0)][(", 
                                 signif(ylim[1], digits = 6), "):1:(", 
                                 signif(ylim[2], digits = 6), ")][(", 
                                 signif(xlim[1], digits = 7), "):1:(", 
                                 signif(xlim[2], digits = 7), ")]")
      
      test <- download.file(sst0125lag30dURL,
                            destfile = out.file.name.sst.0125.lag30d,
                            mode='wb')
    }
    
    if(!file.exists(out.file.name)){
      
      # Multi-scale Ultra-high resolution (MUR) SST analysis fv04.1, Monthly
      sstURL <- paste0("https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst[(", 
                       latlon[k, "Date"], "T00:00:00Z):1:(", 
                       latlon[k, "Date"], "T00:00:00Z)][(", 
                       signif(ylim[1], digits = 4), "):1:(", 
                       signif(ylim[2], digits = 4), ")][(", 
                       signif(xlim[1], digits = 5), "):1:(", 
                       signif(xlim[2], digits = 5), ")]")
      
      # Aqua MODIS, NPP, 4km, Daytime 2003-present (8 day composite)
      # sstURL <- paste0("https://upwell.pfeg.noaa.gov/erddap/griddap/erdMH1sstd8day.nc?sst[(", 
      #                  latlon[k, "Date"], "T00:00:00Z):1:(", latlon[k, "Date"], ")][(", 
      #                  ylim[2], "):1:(", ylim[1], "][(", 
      #                  xlim[1], "):1:(", xlim[2], ")]")
      
      test <- download.file(sstURL,
                            destfile = out.file.name,
                            mode='wb')
      
    }
    
    if (!file.exists(out.file.name.anom)){
      # Multi-scale Ultra-high resolution (MUR) SST analysis anomaly fv04.1, Monthly
      sstanomURL <- paste0("https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41anommday.nc?sstAnom[(",
                           latlon[k, "Date"], "T00:00:00Z):1:(", 
                           latlon[k, "Date"], "T00:00:00Z)][(", 
                           signif(ylim[1], 4), "):1:(", 
                           signif(ylim[2], 4), ")][(", 
                           signif(xlim[1], 5), "):1:(", 
                           signif(xlim[2], 5), ")]")
      test <- download.file(sstanomURL,
                            destfile = out.file.name.anom,
                            mode='wb')
      
    }
    
    if (!file.exists(out.file.name.lag30d)){
      # Multi-scale Ultra-high resolution (MUR) SST analysis anomaly fv04.1, Monthly
      sst30dURL <- paste0("https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst[(", 
                          (latlon[k, "Date"] - 29), "T00:00:00Z):1:(", 
                          (latlon[k, "Date"] - 29), "T00:00:00Z)][(", 
                          signif(ylim[1], 4), "):1:(", 
                          signif(ylim[2], 4), ")][(", 
                          signif(xlim[1], 5), "):1:(", 
                          signif(xlim[2], 5), ")]")
      
      test <- download.file(sst30dURL,
                            destfile = out.file.name.lag30d,
                            mode='wb')
      
    }
    
    if (!file.exists(out.file.name.anom.lag30d)){
      # Multi-scale Ultra-high resolution (MUR) SST analysis anomaly fv04.1, Monthly
      sstanomURL <- paste0("https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41anommday.nc?sstAnom[(",
                           (latlon[k, "Date"] - 29), "T00:00:00Z):1:(", 
                           (latlon[k, "Date"] - 29), "T00:00:00Z)][(", 
                           signif(ylim[1], 4), "):1:(", signif(ylim[2], 4), ")][(", 
                           signif(xlim[1], 5), "):1:(", signif(xlim[2], 5), ")]")
      test <- download.file(sstanomURL,
                            destfile = out.file.name.anom.lag30d,
                            mode='wb')
      
    }
    
    # havent' found good airtemp dataset yet.
    # if (!file.exists(out.file.name.air.temp)){
    #   # NOAA NOS SOS, Experimental, 1853-present, Air temperature
    #   # (https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosATemp.html)
    #   airtempURL <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosATemp.nc?longitude%2Clatitude%2Ctime%2Cair_temperature%2Cquality_flags&longitude%3E=", 
    #                        xlim[1], "&longitude%3C=", xlim[2], 
    #                        "&latitude%3E=", ylim[1], "&latitude%3C=", ylim[2], 
    #                        "&time%3E=", latlon[k, "Date"], 
    #                        "T00%3A00%3A00Z&time%3C=", latlon[k, "Date"], "T00%3A00%3A00Z")
    #   
    #   test <- download.file(airtempURL,
    #                         destfile = out.file.name.air.temp,
    #                         mode='wb')
    #   
    #   datafileID <- nc_open(out.file.name.air.temp)
    #   lon <- ncvar_get(datafileID, varid="longitude")
    #   lat <- ncvar_get(datafileID, varid="latitude")
    #   time <- ncvar_get(datafileID, varid="time")
    #   time <- as.POSIXlt(time, origin='1970-01-01',tz= "GMT")  
    #   mat <- ncvar_get(datafileID, varid = '')
    #   nc_close(datafileID)
    #   
    #   air_temp[k, k1] <- mean(mat, na.rm = T)
    #   
    #   
    # }
    
    
  }
  if (!file.exists(out.file.name.wind)){
    # wind speed and direction data end at 2018-10-09: 
    
    # NOAA/NCDC blended daily global 0.25-degree sea surface winds, 1987-2011, Lon +/-180 
    # (https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOwDly_LonPM180.html)
    
    #longitude needs to be increments of 0.25 degrees:
    # because we are always in the negative longitude values...
    cut.off.lon <- c(0.75, 0.5, 0.25, 0)
    dif.lon <- abs(trunc(latlon[k, "X"]) - cut.off.lon - latlon[k, "X"])
    X.cut.off <- trunc(latlon[k, "X"]) - cut.off.lon[dif.lon == min(dif.lon)]
    X.min <- ifelse(X.cut.off >= latlon[k, "X"], X.cut.off - 0.25, X.cut.off)
    X.max <- X.min + 0.25
    
    # and latitude is always positive in our case
    cut.off.lat <- c(0, 0.25, 0.5, 0.75)
    dif.lat <- abs(latlon[k, "Y"] - cut.off.lat - latlon[k, "Y"])
    Y.cut.off <- trunc(latlon[k, "Y"]) - cut.off.lat[dif.lat == min(dif.lat)]
    Y.min <- ifelse(Y.cut.off >= latlon[k, "Y"], Y.cut.off - 0.25, Y.cut.off)
    Y.max <- Y.min + 0.25
    
    if (lubridate::year(latlon[k, "Date"]) < 2011){
      var.name <- "ncdcOwDly"
    } else {
      var.name <-"ncdcOwDlyP" 
      
    }
    
    if (latlon[k, "Date"] < as.Date("2018-10-10")){
      windURL <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/", var.name, "_LonPM180.nc?u[(",
                        latlon[k, "Date"], "T09:00:00Z):1:(", latlon[k, "Date"],
                        "T09:00:00Z)][(10.0):1:(10.0)][(", Y.min, "):1:(", Y.max, ")][(", 
                        X.min, "):1:(", X.max, ")],v[(", latlon[k, "Date"], "T09:00:00Z):1:(", 
                        latlon[k, "Date"], "T09:00:00Z)][(10.0):1:(10.0)][(", Y.min, "):1:(", 
                        Y.max, ")][(", X.min, "):1:(", X.max, ")],w[(", latlon[k, "Date"],
                        "T09:00:00Z):1:(", latlon[k, "Date"], "T09:00:00Z)][(10.0):1:(10.0)][(", 
                        Y.min, "):1:(", Y.max, ")][(", X.min, "):1:(", X.max, ")]")
      
      test <- download.file(windURL,
                            destfile = out.file.name.wind,
                            mode='wb')
      
    }
  }  
  
 
  if (!file.exists(out.file.name.wind.30d)){
    # wind speed and direction data end at 2018-10-09: 
    
    # NOAA/NCDC blended daily global 0.25-degree sea surface winds, 1987-2011, Lon +/-180 
    # (https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOwDly_LonPM180.html)
    
    #longitude needs to be increments of 0.25 degrees:
    # because we are always in the negative longitude values...
    cut.off.lon <- c(0.75, 0.5, 0.25, 0)
    dif.lon <- abs(trunc(latlon[k, "X"]) - cut.off.lon - latlon[k, "X"])
    X.cut.off <- trunc(latlon[k, "X"]) - cut.off.lon[dif.lon == min(dif.lon)]
    X.min <- ifelse(X.cut.off >= latlon[k, "X"], X.cut.off - 0.25, X.cut.off)
    X.max <- X.min + 0.25
    
    # and latitude is always positive in our case
    cut.off.lat <- c(0, 0.25, 0.5, 0.75)
    dif.lat <- abs(latlon[k, "Y"] - cut.off.lat - latlon[k, "Y"])
    Y.cut.off <- trunc(latlon[k, "Y"]) - cut.off.lat[dif.lat == min(dif.lat)]
    Y.min <- ifelse(Y.cut.off >= latlon[k, "Y"], Y.cut.off - 0.25, Y.cut.off)
    Y.max <- Y.min + 0.25
    
    if (latlon[k, "Date"] < as.Date("2011-10-01")){
      var.name <- "ncdcOwDly" # to 2011-09-30
      windURL <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/", 
                        var.name, "_LonPM180.nc?u[(",
                        latlon[k, "Date"]-29, "T09:00:00Z):1:(", 
                        latlon[k, "Date"],
                        "T09:00:00Z)][(10.0):1:(10.0)][(", 
                        Y.min, "):1:(", Y.max, ")][(", 
                        X.min, "):1:(", X.max, ")],v[(", 
                        latlon[k, "Date"]-29, "T09:00:00Z):1:(", 
                        latlon[k, "Date"], "T09:00:00Z)][(10.0):1:(10.0)][(", 
                        Y.min, "):1:(", 
                        Y.max, ")][(", X.min, "):1:(", X.max, ")],w[(", 
                        latlon[k, "Date"]-29,
                        "T09:00:00Z):1:(", latlon[k, "Date"], 
                        "T09:00:00Z)][(10.0):1:(10.0)][(", 
                        Y.min, "):1:(", Y.max, ")][(", X.min, 
                        "):1:(", X.max, ")]")
      
      test <- download.file(windURL,
                            destfile = out.file.name.wind.30d,
                            mode='wb')
    } else {
      var.name <-"ncdcOwDlyP" # from 2011-10-01
      if((latlon[k, "Date"] - 30) < as.Date("2011-10-01")){
        var.name.2 <- "ncdcOwDly"
        windURL <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/", var.name.2,
                          "_LonPM180.nc?u[(",
                          latlon[k, "Date"]-29, 
                          "T09:00:00Z):1:(2011-09-30T09:00:00Z)][(10.0):1:(10.0)][(",
                          Y.min, "):1:(", Y.max, ")][(", 
                          X.min, "):1:(", X.max, ")],v[(", latlon[k, "Date"]-29,
                          "T09:00:00Z):1:(2011-09-30T09:00:00Z)][(10.0):1:(10.0)][(", Y.min, "):1:(", 
                          Y.max, ")][(", X.min, "):1:(", X.max, ")],w[(", latlon[k, "Date"]-29,
                          "T09:00:00Z):1:(2011-09-30T09:00:00Z)][(10.0):1:(10.0)][(", 
                          Y.min, "):1:(", Y.max, ")][(", X.min, "):1:(", X.max, ")]")
        
        test <- download.file(windURL,
                              destfile = paste0(unlist(out.file.name.wind.30d, "\\.nc"), "_dat1.nc"),
                              mode='wb')
        
        windURL <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/", var.name, 
                          "_LonPM180.nc?u[(2011-10-01T09:00:00Z):1:(", latlon[k, "Date"],
                          "T09:00:00Z)][(10.0):1:(10.0)][(", Y.min, "):1:(", Y.max, ")][(", 
                          X.min, "):1:(", X.max, ")],v[(2011-10-01T09:00:00Z):1:(", 
                          latlon[k, "Date"], "T09:00:00Z)][(10.0):1:(10.0)][(", Y.min, "):1:(", 
                          Y.max, ")][(", X.min, "):1:(", X.max, ")],w[(2011-10-01T09:00:00Z):1:(",
                          latlon[k, "Date"], "T09:00:00Z)][(10.0):1:(10.0)][(", 
                          Y.min, "):1:(", Y.max, ")][(", X.min, "):1:(", X.max, ")]")
        
        test <- download.file(windURL,
                              destfile = out.file.name.wind.30d,
                              mode='wb')
      } else {
        if (latlon[k, "Date"] < as.Date("2018-10-10")){
          
          windURL <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/", var.name,
                            "_LonPM180.nc?u[(",
                            latlon[k, "Date"]-29, "T09:00:00Z):1:(", latlon[k, "Date"],
                            "T09:00:00Z)][(10.0):1:(10.0)][(", Y.min, "):1:(", Y.max, ")][(", 
                            X.min, "):1:(", X.max, ")],v[(", latlon[k, "Date"]-29, "T09:00:00Z):1:(", 
                            latlon[k, "Date"], "T09:00:00Z)][(10.0):1:(10.0)][(", Y.min, "):1:(", 
                            Y.max, ")][(", X.min, "):1:(", X.max, ")],w[(", latlon[k, "Date"]-29,
                            "T09:00:00Z):1:(", latlon[k, "Date"], "T09:00:00Z)][(10.0):1:(10.0)][(", 
                            Y.min, "):1:(", Y.max, ")][(", X.min, "):1:(", X.max, ")]")
          
          test <- download.file(windURL,
                                destfile = out.file.name.wind.30d,
                                mode='wb')
        }
      }
      
      
    }
  }
  
}
