# using a hybrid ncdf and raster approach for querying CMIP runs

library(raster)
library(tidyverse)
library(ncdf4)
library(ncdf4.helpers)
library(maps)
library(mapdata)
library(fields)
library(chron)
library(zoo)
library(ggpubr)


# get some info using ncdf
nc <- nc_open("./CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc")

nc

# query time
d <- nc.get.time.series(nc, v = "tos",
                        time.dim.name = "time")

d <- as.Date(format(d, "%Y-%m-%d"))
year <- years(d)
month <- months(d)

# now lat/long
# these are on 320 x 384 i,j grid

x <- ncvar_get(nc, "lon")
# lon are constant across rows, so we just need one column!
lon.dim <- x[,1]

y <- ncvar_get(nc, "lat")
# lat are constant across columns, so we just need one row!
lat.dim <- y[1,]

# now get the actual data using raster
fname <- "CMIP/tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc" #You will need to specify your location here
nc.b <- brick(fname, varname="tos")  # ts = sst in degrees K

dim(nc.b) # lon, lat, d
nc.b
plot(nc.b[[1]])



# specify area to extract

# first, north pacific - north limit appears off???
lon.min <- which.min(abs(lon.dim-120))
lon.max <- which.min(abs(lon.dim-250))

lat.min <- which.min(abs(lat.dim-20))
lat.max <- which.min(abs(lat.dim-68))

n.pac <- extent(lon.min, lon.max, lat.min, lat.max)

nc.1 <- nc.b[[1]]

np.crop <- crop(nc.1, n.pac)
plot(np.crop)


# now GOA
# for some reason the latitude seems to be off by ten degrees!!
lon.min <- which.min(abs(lon.dim-209))
lon.max <- which.min(abs(lon.dim-231))

lat.min <- which.min(abs(lat.dim-39))
lat.max <- which.min(abs(lat.dim-51))

goa <- extent(lon.min, lon.max, lat.min, lat.max)

# goa <- extent(205, 240, 315, 330)
nc.1 <- nc.b[[1]]

goa.crop <- crop(nc.1, goa)
plot(goa.crop)

lat.dim

dim(x)
nc.b
