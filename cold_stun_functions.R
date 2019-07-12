

latlon2sp <- function(in.df, center.UTM, zone){
  coordinates(in.df) <- c("X", "Y")
  proj4string(in.df) <- CRS("+proj=longlat +datum=WGS84")
  out.df <- spTransform(in.df, CRS(paste0("+proj=utm +zone=", zone, " ellps=WGS84")))
  out.df$newX <- (out.df$X - center.UTM$X)/1000
  out.df$newY <- (out.df$Y - center.UTM$Y)/1000
  return(out.df)
}

# function to convert the new coordinate system back to lat/lon
sp2latlon <- function(in.xy, center.UTM, zone){
  X <- in.xy$newX * 1000 + center.UTM$X
  Y <- in.xy$newY * 1000 + center.UTM$Y
  in.df <- data.frame(X = X, Y = Y)
  coordinates(in.df) <- c('X', 'Y')
  proj4string(in.df) <- CRS(paste0("+proj=utm +zone=", zone, " ellps=WGS84"))
  out.df <- spTransform(in.df, CRS("+proj=longlat +datum=WGS84"))
  return(out.df)
}

get_ocedata_fcn <- function(ds, dat0, dataset = "coldstun"){
  
  if (dataset != "coldstun"){
    id.name <- "turtle_all"
  } else {
    id.name <- "turtle"
  }
  
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
  
  sst01 <- sst01_anom <- sst01_sd <- sst01_lag30d <- sst01_anom_30d <- matrix(nrow = nrow(latlon),
                                                                              ncol = d.n)
  
  sst0125 <- sst0125_min <- sst0125_sd <-  matrix(nrow = nrow(latlon),
                                                  ncol = d.n)
  
  sst01_min <- sst01_min_lag30d <- air_temp <- air_temp_30d <- matrix(nrow = nrow(latlon),
                                                                  ncol = d.n)
  sst0125_lag30d <- sst0125_min_lag30d <- sst0125_sd_lag30d <- matrix(nrow = nrow(latlon),
                                                                      ncol = d.n)
  
  wind <- wind_max <- wind_30d <- wind_max_30d <- matrix(nrow = nrow(latlon),
                                                         ncol = 1)
  k <- k1 <- 1
  
  # for all rows
  for (k in 1:nrow(dat0)){
    # wind comes in 0.25 degree resolutions so we just pick one without worrying about other spatial
    # resolutions
    # out.file.name.wind <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
    #                              "wind.nc")
    # 
    # out.file.name.wind.30d <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
    #                                  "wind_30d.nc")
    
    # for number of distances
    for (k1 in 1:d.n){
      out.file.name <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                              ds[k1], "km_sst.nc")
      out.file.name.anom <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                                   ds[k1], "km_sstanom.nc")
      out.file.name.lag30d<- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                                    ds[k1], "km_sst_lag30d.nc")
      out.file.name.anom.lag30d <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                                          ds[k1], "km_sstanom_lag30d.nc")
      out.file.name.air.temp<- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                                      ds[k1], "km_airtemp.nc")
      out.file.name.air.temp.lag30d <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                                              ds[k1], "km_airtemp_lag30d.nc")
      out.file.name.sst.0125 <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                                       ds[k1], "km_sst0125.nc")
      out.file.name.sst.0125.lag30d <- paste0("data/ncfiles/", id.name, "_", dat0[k, "ID"], "_",
                                              ds[k1], "km_sst0125_lag30d.nc")
      
      
      xlim <- c(latlon[k, (4*k1 - 1)],
                latlon[k, (4*k1 + 1)])
      ylim <- c(latlon[k, (4*k1)],
                latlon[k, (4*k1 + 2)])
      #tlim <- c(latlon[k, "Date"] - 7,
      #          latlon[k, "Date"] + 7)
      
      datafileID <- nc_open(out.file.name.sst.0125)
      lon <- ncvar_get(datafileID, varid="longitude")
      lat <- ncvar_get(datafileID, varid="latitude")
      time <- ncvar_get(datafileID, varid="time")
      time <- as.POSIXlt(time, origin='1970-01-01',tz= "GMT")  
      mat <- ncvar_get(datafileID, varid = 'sst')
      nc_close(datafileID)
      
      sst0125[k, k1] <- mean(mat, na.rm = T)
      sst0125_sd[k, k1] <- sqrt(var(as.vector(mat), na.rm = T))
      sst0125_min[k, k1] <- min(mat, na.rm = T)
      
      datafileID <- nc_open(out.file.name.sst.0125.lag30d)
      lon <- ncvar_get(datafileID, varid="longitude")
      lat <- ncvar_get(datafileID, varid="latitude")
      time <- ncvar_get(datafileID, varid="time")
      time <- as.POSIXlt(time, origin='1970-01-01',tz= "GMT")  
      mat <- ncvar_get(datafileID, varid = 'sst')
      nc_close(datafileID)
      
      sst0125_lag30d[k, k1] <- mean(mat, na.rm = T)
      sst0125_sd_lag30d[k, k1] <- sqrt(var(as.vector(mat), na.rm = T))
      sst0125_min_lag30d[k, k1] <- min(mat, na.rm = T)
      
      datafileID <- nc_open(out.file.name)
      lon <- ncvar_get(datafileID, varid="longitude")
      lat <- ncvar_get(datafileID, varid="latitude")
      time <- ncvar_get(datafileID, varid="time")
      time <- as.POSIXlt(time, origin='1970-01-01',tz= "GMT")  
      mat <- ncvar_get(datafileID, varid = 'sst')
      nc_close(datafileID)
      
      sst01[k, k1] <- mean(mat, na.rm = T)
      sst01_sd[k, k1] <- sqrt(var(as.vector(mat), na.rm = T))
      sst01_min[k, k1] <- min(mat, na.rm = T)

      datafileID <- nc_open(out.file.name.anom)
      lon <- ncvar_get(datafileID, varid="longitude")
      lat <- ncvar_get(datafileID, varid="latitude")
      time <- ncvar_get(datafileID, varid="time")
      time <- as.POSIXlt(time, origin='1970-01-01',tz= "GMT")  # true? check with nc files (5/30/2019)
      mat <- ncvar_get(datafileID, varid = 'sstAnom')
      nc_close(datafileID)
      
      sst01_anom[k, k1] <- mean(mat, na.rm = T)
      
      datafileID <- nc_open(out.file.name.lag30d)
      lon <- ncvar_get(datafileID, varid="longitude")
      lat <- ncvar_get(datafileID, varid="latitude")
      time <- ncvar_get(datafileID, varid="time")
      time <- as.POSIXlt(time, origin='1970-01-01',tz= "GMT")  
      mat <- ncvar_get(datafileID, varid = 'sst')
      nc_close(datafileID)
      
      sst01_lag30d[k, k1] <- mean(mat, na.rm = T)
      sst01_min_lag30d[k, k1] <- min(mat, na.rm = T)
      
      datafileID <- nc_open(out.file.name.anom.lag30d)
      lon <- ncvar_get(datafileID, varid="longitude")
      lat <- ncvar_get(datafileID, varid="latitude")
      time <- ncvar_get(datafileID, varid="time")
      time <- as.POSIXlt(time, origin='1970-01-01',tz= "GMT")  
      mat <- ncvar_get(datafileID, varid = 'sstAnom')
      nc_close(datafileID)
      
      sst01_anom_30d[k, k1] <- mean(mat, na.rm = T)
      
    }
    
  #   if (latlon[k, "Date"] < as.Date("2018-10-10")){
  #     datafileID <- nc_open(out.file.name.wind)
  #     lon <- ncvar_get(datafileID, varid="longitude")
  #     lat <- ncvar_get(datafileID, varid="latitude")
  #     time <- ncvar_get(datafileID, varid="time")
  #     time <- as.POSIXlt(time, origin='1970-01-01', tz= "GMT")  
  #     u <- ncvar_get(datafileID, varid = "u")
  #     v <- ncvar_get(datafileID, varid = "v")
  #     w <- ncvar_get(datafileID, varid = "w")
  #     
  #     nc_close(datafileID)
  #   } else {
  #     U <- v <- w <- NA
  #   }
  #   wind[k, 1] <- mean(w, na.rm = T)
  #   wind_max[k, 1] <- max(w, na.rm = T)
  #   #wind_u[k, 1] <- mean(u, na.rm = T)
  #   #wind_v[k, 1] <- mean(v, na.rm = T)
  #   
  #   if (((latlon[k, "Date"] > as.Date("2011-09-30")) & 
  #        (latlon[k, "Date"] - 30) < as.Date("2011-10-01"))){
  #     
  #     datafileID <- nc_open(out.file.name.wind.30d)
  #     
  #     datafileID2 <- nc_open(paste0(unlist(out.file.name.wind.30d, "\\.nc"), "_dat1.nc"))
  #     
  #     lon <- ncvar_get(datafileID, varid="longitude")
  #     lat <- ncvar_get(datafileID, varid="latitude")
  #     time <- ncvar_get(datafileID, varid="time")
  #     time <- as.POSIXlt(time, origin='1970-01-01', tz= "GMT")  
  #     u <- ncvar_get(datafileID, varid = "u")
  #     v <- ncvar_get(datafileID, varid = "v")
  #     w <- ncvar_get(datafileID, varid = "w")
  #     
  #     lon2 <- ncvar_get(datafileID2, varid="longitude")
  #     lat2 <- ncvar_get(datafileID2, varid="latitude")
  #     time2 <- ncvar_get(datafileID2, varid="time")
  #     time2 <- as.POSIXlt(time2, origin='1970-01-01', tz= "GMT")  
  #     u2 <- ncvar_get(datafileID2, varid = "u")
  #     v2 <- ncvar_get(datafileID2, varid = "v")
  #     w2 <- ncvar_get(datafileID2, varid = "w")
  #     
  #     nc_close(datafileID)
  #     nc_close(datafileID2)
  #     
  #     u <- abind::abind(u, u2)
  #     v <- abind::abind(v, v2)
  #     w <- abind::abind(w, w2)
  #     
  #   } else {
  #     if (latlon[k, "Date"] < as.Date("2018-10-10")){
  #       datafileID <- nc_open(out.file.name.wind.30d)
  #       lon <- ncvar_get(datafileID, varid="longitude")
  #       lat <- ncvar_get(datafileID, varid="latitude")
  #       time <- ncvar_get(datafileID, varid="time")
  #       time <- as.POSIXlt(time, origin='1970-01-01', tz= "GMT")  
  #       u <- ncvar_get(datafileID, varid = "u")
  #       v <- ncvar_get(datafileID, varid = "v")
  #       w <- ncvar_get(datafileID, varid = "w")
  #       
  #       nc_close(datafileID)
  #     } else {
  #       w <- u <- v <- NA
  #     }
  #   }
  #   
  #   wind_30d[k, 1] <- mean(w, na.rm = T)
  #   wind_max_30d[k, 1] <- max(w, na.rm = T)
  #   #wind_u_30d[k, 1] <- mean(u, na.rm = T)
  #   #wind_v_30d[k, 1] <- mean(v, na.rm = T)
  #   
  }
  
  sst01.df <- data.frame(sst01)
  colnames(sst01.df) <- c("SST_1", "SST_2", "SST_3")
  
  sst01.anom.df <- data.frame(sst01_anom)
  colnames(sst01.anom.df) <- c("SST_1_ANOM", "SST_2_ANOM", "SST_3_ANOM")
  
  sst01.sd.df <- data.frame(sst01_sd)
  colnames(sst01.sd.df) <- c("SST_1_SD", "SST_2_SD", "SST_3_SD")
  
  sst01.min.df <- data.frame(sst01_min)
  colnames(sst01.min.df) <- c("SST_1_min", "SST_2_min", "SST_3_min")
  
  sst01.lag30d.df <- data.frame(sst01_lag30d)
  colnames(sst01.lag30d.df) <- c("SST_1_lag30", "SST_2_lag30", "SST_3_lag30")
  
  sst01.min.lag30d.df <- data.frame(sst01_min_lag30d)
  colnames(sst01.min.lag30d.df) <- c("SST_1_min_lag30", "SST_2_min_lag30", "SST_3_min_lag30")
  
  sst01.anom.lag30d.df <- data.frame(sst01_anom_30d)
  colnames(sst01.anom.lag30d.df) <- c("SST_1_ANOM_lag30", "SST_2_ANOM_lag30", "SST_3_ANOM_lag30")
  
  sst0125.df <- data.frame(sst0125)
  colnames(sst0125.df) <- c("SST0125_1", "SST0125_2", "SST0125_3")
  
  sst0125.sd.df <- data.frame(sst0125_sd)
  colnames(sst0125.sd.df) <- c("SST0125_1_SD", "SST0125_2_SD", "SST0125_3_SD")
  
  sst0125.min.df <- data.frame(sst0125_min)
  colnames(sst0125.min.df) <- c("SST0125_1_min", "SST0125_2_min", "SST0125_3_min")
  
  sst0125.lag30d.df <- data.frame(sst0125_lag30d)
  colnames(sst0125.lag30d.df) <- c("SST0125_1_lag30d", "SST0125_2_lag30d", "SST0125_3_lag30d")
  
  sst0125.sd.lag30d.df <- data.frame(sst0125_sd_lag30d)
  colnames(sst0125.sd.lag30d.df) <- c("SST0125_1_SD", "SST0125_2_SD", "SST0125_3_SD")
  
  sst0125.min.lag30d.df <- data.frame(sst0125_min_lag30d)
  colnames(sst0125.min.lag30d.df) <- c("SST0125_1_min_lag30", "SST0125_2_min_lag30", "SST0125_3_min_lag30")
  
  # wind.df <- data.frame(wind)
  # wind.max.df <- data.frame(wind_max)
  # wind.30d.df <- data.frame(wind_30d)
  # wind.max.30d.df <- data.frame(wind_max_30d)
  
  out.list <- list(sst0125.df = sst0125.df, 
                   sst0125.sd.df = sst0125.sd.df, 
                   sst0125.min.df = sst0125.min.df,
                   sst0125.lag30d.df = sst0125.lag30d.df, 
                   sst0125.sd.lag30d.df = sst0125.sd.lag30d.df, 
                   sst0125.min.lag30d.df = sst0125.min.lag30d.df,
                   sst01.df = sst01.df, 
                   sst01.sd.df = sst01.sd.df, 
                   sst01.min.df = sst01.min.df, 
                   sst01.anom.df = sst01.anom.df, 
                   sst01.lag30d.df = sst01.lag30d.df,
                   sst01.min.lag30d.df = sst01.min.lag30d.df, 
                   sst01.anom.lag30d.df = sst01.anom.lag30d.df) 
                   # wind.df = wind.df, 
                   # wind.max.df = wind.max.df, 
                   # wind.30d.df = wind.30d.df, 
                   # wind.max.30d.df = wind.max.30d.df)
  return(out.list)
}
