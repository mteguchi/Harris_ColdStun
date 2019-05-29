

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
