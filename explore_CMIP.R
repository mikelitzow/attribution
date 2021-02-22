## explore CMIP preindustrial downloads

library(tidyverse)
library(ncdf4)
library(ncdf4.helpers)
library(maps)
library(mapdata)
library(fields)
library(chron)
library(zoo)
library(ggpubr)
library(s2dverification)

nc <- nc_open("./CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

nc

d <- nc.get.time.series(nc, v = "tos",
                        time.dim.name = "time")

d <- as.Date(format(d, "%Y-%m-%d"))

lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

SST <- ncvar_get(nc, "tos")

dim(SST) # 329 lons, 384 lats, 6000 months

names(dim(SST)) <- c("lon", "lat", "time")

names(dim(lon)) <- names(dim(lat)) <- c("lon", "lat")

CDORemap(SST, 
         lons=lon,
         lats=lat,
         grid='r180x90',
         method='bilinear',
         crop=c(180, 240, 40, 70))



# lon are constant across rows, so we just need one column!
lon.dim <- x[,1]


# lat are constant across columns, so we just need one row!
lat.dim <- y[1,]

# for the GOA, we need lat 49 - 61 and long 209 - 231
i.min <- which.min(abs(lon.dim-150))
i.max <- which.min(abs(lon.dim-250))
i.count <- 1 + i.max-i.min


j.min <- which.min(abs(lat.dim-20))
j.max <- which.min(abs(lat.dim-68))
j.count <- 1 + j.max-j.min

# and get date vector


# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
x <- lon.dim[i.min:i.max]
y <- lat.dim[j.min:j.max]

lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

lat <- lat+10
y <- y+10

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)


dim(sst) # 320 x 384 x 6000
dim(x) # 320 x 384
dim(y) # 320 x 384
x[,1] # 320 (different longitudes)
y[,1] # 320 (same) latitudes

x[1,] # 384 (slightly decreasing) longitudes
y[1,] # 384 latitudes, pole to pole

# try w/ some help from the web
d <- nc.get.time.series(nc, v = "tos",
                               time.dim.name = "time")

d <- as.Date(format(d, "%Y-%m-%d"))

x.keep <- as.vector(x > 130 & x < 250)
y.keep <- as.vector(y > 20 & y < 64)
# 
# lon_index <- x[x > 130 & x < 250]
# lat_index <- y[y > 20 & y < 64]

nc.get.dim.for.axis(nc, "tos", "X")

sst <- nc.get.var.subset.by.axes(nc, "tos",
                                 axis.indices = list(X = x.keep,
                                                     Y = y.keep))
# frustrating!
# I know this should be easy!

sst <- ncvar_get(nc, "tos")

nc.get.dim.axes(nc)

nc.get.coordinate.axes(nc, "tos")

nc.is.regular.dimension(nc)

nc.get.variable.list(nc)

nc.get.dim.axes.from.names(nc)

nc.get.coordinate.axes(nc, "tos")

nc.get.proj4.string(nc, "tos")

sst <- nc.get.var.subset.by.axes(nc, "tos",
                                 axis.indices = list(lon = lon_index,
                                                     lat = lat_index))
data_frame(time = tas_time, 
           tas = as.vector(tas)) %>%
  mutate(time = as.Date(format(time, "%Y-%m-%d"))) %>%
  ggplot(aes(x = time, y = tas)) + 
  geom_line() + 
  xlab("Date") + ylab("Temperature (K)") + 
  ggtitle("Daily modeled near-surface air temprature, 2071-2075",
          subtitle = "At model grid point nearest Beijing, China") + 
  theme_classic()

################################################################################

library(tidync)


nc <- system.file("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc", package = "stars")
nc <- tidync("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

nc

ncmeta::nc_grids("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

ncmeta::nc_vars("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

ncmeta::nc_dim("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc", 0)

ncmeta::nc_atts("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

ncmeta::nc_axes("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

ncmeta::nc_dims("CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

(nc_data <- nc %>% hyper_array())

names(nc_data)

dim(nc_data[[1]])

image(nc_data[[1]])

str(nc_data[[1]])

str(nc_data[[2]])

lapply(nc_data, dim)
