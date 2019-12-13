library(tidyverse)
library(ncdf4)
library(maps)
library(mapdata)
library(fields)
library(chron)
library(zoo)
library(ggpubr)

# 
# This comparison uses output from 5 CMIP5 models that were used in the attribution study published in Walsh et al. 2018 BAMS. 
# Thanks to John Walsh for making these summaries available!
#   
# Here are the preindustrial and historical/RCP8.5 (1987-2005 / 2006-2046) simulations:
  
# load

dat <- read.csv("CMIP5 GOA SST.csv")

dat <- dat %>%
  gather(model, anomaly, -Year, -Era)

dat$Era <- ifelse(dat$Era=="present",
                  "historical / RCP8.5", as.character(dat$Era))

dat$Era <- reorder(dat$Era, desc(dat$Era))

ggplot(dat, aes(Year, anomaly, color=model)) +
  theme_bw() +
  geom_line() + 
  facet_wrap(~Era, scales="free_x")

preindustrial <- dat %>%
  filter(Era=="preindustrial")

ggplot(preindustrial, aes(anomaly, fill=model)) +
  theme_bw() +
  geom_density(alpha=0.3) +
  xlim(-5,5)

# We will compare the preindustrial estimates above with ERSSTv5 observations. The observations are for the same area (50º-60ºN, 150º-130ºW). Here is average SST for that area from ERSSTv5:

# quickly (!) estimate 2019 annual anomaly based on data avaliable so far!
# identify latest year and month needed
year <- 2019
month <- "11"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(50):1:(60)][(210):1:(230)]", sep="")

download.file(URL, "GOA.box.ersst.latest")

# process
nc <- nc_open("GOA.box.ersst.latest")

# extract dates
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,240), ylim=c(40,66))
contour(x, y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# get anomaly wrt 1981:2010 mean (1900:2016 SD)
yr <- as.numeric(as.character(years(d)))
m <- months(d)

f <- function(x) tapply(x, yr, mean)

annual.sst <- apply(SST, 2, f) # annual mean for each cell

# get anomaly for each cell
mu <- apply(annual.sst[rownames(annual.sst) %in% 1981:2010,], 2, mean)	# Compute monthly means for each cell during 1981:2010

# change suggested by Nick Bond - using SD for the 1900:2016 time series when calculating anomalies
# doing this b.c my previous anomalies appear to be too large when compared with Walsh et al. 2018
sd <- apply(annual.sst[rownames(annual.sst) %in% 1900:2016,], 2, sd)

mu <- matrix(mu, nrow=nrow(annual.sst), ncol=length(mu), byrow=T)
sd <- matrix(sd, nrow=nrow(annual.sst), ncol=length(sd), byrow=T) 

anom <- (annual.sst - mu)/sd # now each value is a cell-specific normalized anomaly relative to 1981-2010 

weights <-  sqrt(cos(lat*pi/180))
f <- function(x) weighted.mean(x, w=weights, na.rm=T)
annual.anomaly <- apply(anom, 1, f)

annual.anomaly <- annual.anomaly[names(annual.anomaly) >=1900]

anomaly.plot <- data.frame(year=1900:2018, anomaly=annual.anomaly[names(annual.anomaly) %in% 1900:2018],
                           sign=ifelse(annual.anomaly[names(annual.anomaly) %in% 1900:2018]>0, "positive", "negative"))
anomaly.plot$sign <- reorder(anomaly.plot$sign, desc(anomaly.plot$sign))

##################
# now estimate 2019 value
# limit all years to the range of months available in this year!
m <- as.numeric(m)
keep <- m <= as.numeric(month) # this is the most recent month specified for data query above

short.sst <- SST[keep,]
short.weighted.mean <- apply(short.sst, 1, ff)

short.sst <- tapply(short.weighted.mean, yr[keep], mean)

x <- as.vector(short.sst[names(short.sst) %in% 1900:2018])
y <- annual.anomaly[names(annual.anomaly) %in% 1900:2018]

plot(1900:2018, scale(x), col="blue", type="l")
lines(1900:2018, scale(y), col='red')

md1 <- lm(y ~ x)

# now add 2019 estimate to the annual TS
estimated.2019 <- short.sst[names(short.sst)==2019]*coef(md1)[2] + coef(md1)[1]
annual.anomaly[names(annual.anomaly)==2019] <- estimated.2019

anomaly.plot <- data.frame(year=1900:2019, anomaly=annual.anomaly,
                           sign=as.vector(ifelse(annual.anomaly>0, "positive", "negative")),
                           sm.anom=rollmean(annual.anomaly, 3, fill = NA))

anomaly.plot$sign <- reorder(anomaly.plot$sign, desc(anomaly.plot$sign))

ggplot(anomaly.plot, aes(year, anomaly, fill=sign)) +
  theme_bw() +
  geom_col() +
  geom_line(aes(year, sm.anom, group=1)) +
  theme(legend.position='none')

# Now going ahead with the comparison to preindustrial simulations
preind.obs <- rbind(preindustrial,
                    data.frame(Year=2014:2019, 
                               Era="observation",
                               model=NA, 
                               anomaly=annual.anomaly[names(annual.anomaly) %in% 2014:2019]))

# plot preindustrial distributions with recent
ggplot(preind.obs, aes(anomaly, fill=Era)) +
theme_bw() +
  geom_density(alpha=0.3, color="grey") +
  xlim(-5,5) 

points <- data.frame(labels=c(NA, NA, NA, "2017", "2018", "2014, 2015, 2016, 2019"), 
                     anomaly=annual.anomaly[names(annual.anomaly) %in% 2014:2019],
                     density=0)

# and...run the three-year rolling means on preindustrial simulations and plot that distribution

# now...figure out AR(1) values

f <- function(x) ar(x, aic=F, order.max = 1)$ar

model.ar.1 <- tapply(preindustrial$anomaly, preindustrial$model, f)

observed.ar.1 <- f(annual.anomaly)

model.ar.1; observed.ar.1

# drop CCSM4 as autocorrelation is unreasonably low!
smooth.preindustrial <- preindustrial %>%
  filter(model != "CCSM4")

smooth.preindustrial$sm.anomaly <- NA

# now go through and get smooth anomalies for each model separately
mods <- unique(smooth.preindustrial$model)
  
for(i in 1:length(mods)){
  
  temp <- smooth.preindustrial %>%
    filter(model==mods[i])
  
  smooth.preindustrial$sm.anomaly[smooth.preindustrial$model==mods[i]] <- rollmean(temp$anomaly, 3, fill=NA) 
  
}

# now thin so we don't include overlapping windows!
keepers <- seq(2, 59, 3)

smooth.preindustrial <- smooth.preindustrial %>%
  filter(Year %in% keepers)

sm.obs <- rollmean(annual.anomaly, 3, fill=NA)

points.sm <- data.frame(labels=c("2014-2016", "2017-2019"), 
                     anomaly=sm.obs[names(sm.obs) %in% c(2015, 2018)],
                     density=0)

# and plot individual obs years
# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

a.plot <- ggplot(anomaly.plot, aes(year, anomaly, fill=sign)) +
  theme_bw() +
  geom_col(color="black", size=0.05) +
  scale_fill_manual(values=cb[7:6]) +
  theme(legend.position='none', axis.title.x = element_blank()) +
  ylab("Temperature anomaly") +
  scale_x_continuous(breaks=seq(1900,2020,20)) +
  labs(title="a) SST observations")
  
nudge.b <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.15)

# points <- points %>%
#   arrange(desc(anomaly))

b.plot <- ggplot(preindustrial, aes(anomaly)) +
  theme_bw() +
  geom_density(alpha=0.3, color="grey", fill=cb[3]) +
  xlim(-5,2.5) +
  geom_point(aes(jitter(anomaly), density), points, fill=cb[8], shape=21, size=1.2) +
  geom_text(aes(anomaly, density+nudge.b, label=labels), points, size=3.5) +
  xlab("Temperature anomaly") + ylab("Density") +
  coord_flip() +
  ggtitle("b) Annual values")

nudge.c <- 0.09

c.plot <- ggplot(smooth.preindustrial, aes(anomaly)) +
  theme_bw() +
  geom_density(alpha=0.3, color="grey", fill=cb[3]) +
  xlim(-5,2.5) +
  geom_point(aes(anomaly, density), points.sm, fill=cb[8], shape=21, size=1.2) +
  geom_text(aes(anomaly, density+nudge.c, label=labels), points.sm, size=3.5) +
  xlab("Temperature anomaly") + ylab("Density") +
  coord_flip() +
  ggtitle("c) 3-year running means")

png("joint sst plots.png", 6, 6, units="in", res=300)

ggarrange(a.plot, 
          ggarrange(b.plot, c.plot, ncol=2, nrow=1),
          nrow=2)

dev.off()

# and calculate likelihood/probability

p.annual <- pnorm(points$anomaly, mean(preindustrial$anomaly), sd(preindustrial$anomaly), lower.tail = F)
z.annual <- (points$anomaly - mean(preindustrial$anomaly)) / sd(preindustrial$anomaly)

p.smooth <- pnorm(points.sm$anomaly, mean(preindustrial$sm.anomaly, na.rm=T),
                  sd(preindustrial$sm.anomaly, na.rm=T), lower.tail = F)
z.smooth <- (points.sm$anomaly - mean(preindustrial$sm.anomaly, na.rm = T)) / sd(preindustrial$sm.anomaly, na.rm=T)

# and export
preind.summary <- data.frame(class=c("annual", "3-yr"),
                             mean=c(mean(preindustrial$anomaly), mean(preindustrial$sm.anomaly, na.rm = T)),
                             sd=c(sd(preindustrial$anomaly), sd(preindustrial$sm.anomaly, na.rm = T)),
                             max=c(max(preindustrial$anomaly), max(preindustrial$sm.anomaly, na.rm = T)))
write.csv(preind.summary, "preindustrial summary.csv")

obs.summary <- data.frame(year=c(2014:2019, "2014-2016", "2017-2019"),
                          z = c(z.annual, z.smooth),
                          p = c(p.annual, p.smooth))

obs.summary[,2] <- round(obs.summary[,2], 2)
obs.summary[,3] <- round(obs.summary[,3], 4)

write.csv(obs.summary, "observed summary.csv")

