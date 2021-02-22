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
names(dat)[1] <- "Year"

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

# quickly (!) estimate 2020 annual anomaly based on data avaliable so far!
# identify latest year and month needed
year <- 2020
month <- "10"

# URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(50):1:(60)][(210):1:(230)]", sep="")
# 
# download.file(URL, "GOA.box.ersst.latest")

# process
# nc <- nc_open("GOA.box.ersst.latest")
nc <- nc_open("nceiErsstv5_3a8e_bf19_ba95.nc")

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
# now estimate 2020 value
# limit all years to the range of months available in this year!
m <- as.numeric(m)
keep <- m <= as.numeric(month) # this is the most recent month specified for data query above

short.sst <- SST[keep,]
short.weighted.mean <- apply(short.sst, 1, f)

short.sst <- tapply(short.weighted.mean, yr[keep], mean)

x <- as.vector(short.sst[names(short.sst) %in% 1900:2019])
y <- annual.anomaly[names(annual.anomaly) %in% 1900:2019]

plot(1900:2019, scale(x), col="blue", type="l")
lines(1900:2019, scale(y), col='red')

md1 <- lm(y ~ x)

# now add 2020 estimate to the annual TS
estimated.2020 <- short.sst[names(short.sst)==2020]*coef(md1)[2] + coef(md1)[1]
annual.anomaly[names(annual.anomaly)==2020] <- estimated.2020

anomaly.plot <- data.frame(year=1900:2020, anomaly=annual.anomaly,
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
                    data.frame(Year=1900:2020, 
                               Era="observation",
                               model=NA, 
                               anomaly=annual.anomaly[names(annual.anomaly) %in% 1900:2020]))

# plot preindustrial distributions with recent
ggplot(preind.obs, aes(anomaly, fill=Era)) +
theme_bw() +
  geom_density(alpha=0.3, color="grey") +
  xlim(-5,5) 


# and calculate likelihood/probability compared to each model!
mods <- unique(preindustrial$model)
preind.comp.out <- FAR.out <- data.frame()

for(m in 1:length(mods)){
# m <- 1
temp <- preindustrial %>%
  filter(model==mods[m])

# p.annual <- pnorm(points$anomaly, mean(temp$anomaly), sd(temp$anomaly), lower.tail = F)
#z.annual <- (points$anomaly - mean(preindustrial$anomaly)) / sd(preindustrial$anomaly)

# and FAR

FAR <- obs.prob <- preind.prob <- NA 
# make a new 'points' object for 1900-2020
points <- data.frame(anomaly=annual.anomaly[names(annual.anomaly) %in% 1960:2020])



for(i in 1960:2020){ # loop through each year to calculate FAR
  # i <- 2014
  obs.prob[(i-1959)] <- sum(annual.anomaly[61:120] >= points$anomaly[(i-1959)])/60
  preind.prob[(i-1959)] <- sum(temp$anomaly >= points$anomaly[(i-1959)])/nrow(temp)
  FAR[(i-1959)] <- 1-(preind.prob[(i-1959)]/obs.prob[(i-1959)])
}

  FAR.out <- rbind(FAR.out, 
                   data.frame(model=mods[m], 
                              year=1960:2020,
                              anom=points$anomaly,
                              obs.prob=obs.prob,
                              preind.prob=preind.prob,
                              FAR=FAR))

}

FAR.out <- FAR.out %>%
  arrange(model)
FAR.out

# and collapse into a ts to plot
plot.dat <- FAR.out %>%
  group_by(year) %>%
  summarise(anom=mean(anom),
            observational_prob=mean(obs.prob),
            preindustrial_prob=mean(preind.prob),
            Fraction_Attributable_Risk=mean(FAR)) %>%
  pivot_longer(cols=-year)


ggplot(plot.dat, aes(year, value)) +
  geom_line() +
  facet_wrap(~name, scales="free_y")

# and export
write.csv(FAR.out, "FAR values for cod paper.csv")

## Now calculate FAR values for historical / RCP8.5 runs!

dat <- read.csv("CMIP5 GOA SST.csv")
names(dat)[1] <- "Year"

dat <- dat %>%
  filter(Era=="present") %>% # present is historical / RCP8.5!
  select(-Era) %>%
  pivot_longer(cols = -Year)
  
ggplot(dat, aes(Year, value, color=name)) +
  theme_bw() +
  geom_line() 

## now calculate FAR for these values!
# and calculate likelihood/probability compared to each model!
mods <- unique(preindustrial$model)
preind.comp.out <- FAR_8.5_out <- data.frame()

for(m in 1:length(mods)){
  # m <- 1
  temp.pre <- preindustrial %>%
    filter(model==mods[m])

  temp_8.5 <- dat %>%
    filter(name==mods[m])
  
  # and FAR
  
  FAR_8.5 <- prob_8.5 <- preind.prob <- NA 
 
  for(i in 1987:2046){ # loop through each year to calculate FAR
    # i <- 2014
    prob_8.5[(i-1986)] <- sum(temp_8.5$value >= temp_8.5$value[(i-1986)])/60
    preind.prob[(i-1986)] <- sum(temp.pre$anomaly >= temp_8.5$value[(i-1986)])/60
    FAR_8.5[(i-1986)] <- 1-(preind.prob[(i-1986)]/prob_8.5[(i-1986)])
  }
  
  FAR_8.5_out <- rbind(FAR_8.5_out, 
                   data.frame(model=mods[m], 
                              year=1987:2046,
                              FAR=FAR_8.5))
  
}

## now combine 
FAR.obs <- FAR.out %>%
  group_by(year) %>%
  summarise(FAR=mean(FAR))

write.csv(FAR.out, "ERSST FAR.csv")

FAR.obs$source <- "observed"

names(FAR_8.5_out)[1] <- "source"

FAR_8.5_out <- FAR_8.5_out %>%
  select(year, FAR, source)


write.csv(FAR_8.5_out, "CMIP FAR.csv")

FAR.plot <- rbind(FAR.obs, FAR_8.5_out)

ggplot(FAR.plot, aes(year, FAR, color=source)) +
  geom_line()


# get model / obs FAR spread!

FAR_8.5_mean <- FAR_8.5_out %>%
  group_by(year) %>%
  summarise(FARmu = mean(FAR), FARsd = sd(FAR))

FAR_obs_mean <- FAR.out %>%
  group_by(year) %>%
  summarise(FARmu = mean(FAR), FARsd = sd(FAR))

FAR_obs_mean$source <- "Observed (ERSST)"

FAR_8.5_mean$source <- "Modeled (historical / RCP8.5)"

FAR.plot <- rbind(FAR_obs_mean, FAR_8.5_mean)

FAR.plot$source <- reorder(FAR.plot$source, desc(FAR.plot$source))

FAR.plot$ymin <- FAR.plot$FARmu-1.96*FAR.plot$FARsd/sqrt(5)
FAR.plot$ymax <- FAR.plot$FARmu+1.96*FAR.plot$FARsd/sqrt(5)

FAR.plot$ymin <- ifelse(FAR.plot$ymin < 0, 0, FAR.plot$ymin)
FAR.plot$ymax <- ifelse(FAR.plot$ymax > 1, 1, FAR.plot$ymax)

ggplot(FAR.plot, aes(year, FARmu)) +
  geom_line(color="red") +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.3, linetype=0) +
  facet_wrap(~source, ncol=1) +
  theme_bw() +
  xlim(1960, 2048) +
  geom_hline(yintercept = 0.9, lty=2) +
  theme(axis.title.x = element_blank()) +
  ylab("Fraction Attributable Risk") +
  scale_x_continuous(breaks=seq(1960, 2040, 10))

ggsave("observed and modeled FAR.png", width=5, height=4, units='in')





## old stuff below!


rownames(preind.comp.out) <- preind.comp.out[,1]
preind.comp.out <- preind.comp.out[,-1]

preind.comp.out <- as.data.frame(t(preind.comp.out))

rownames(preind.comp.out) <- 1:6

# and put it together for a table in the paper
FARs <- FAR.out %>%
  select(model, year, anom, FAR) %>%
  spread(model, FAR, -year, -anom)

FARs

# xport <- cbind(FARs, preind.comp.out)

xport <- FARs
# write.csv(xport, "annual anomaly preindustrial P values and FAR.csv", row.names = F)
write.csv(xport, "annual anomaly and FAR.csv", row.names = F)

######################################
# now the statistics for smoothed obs!

mods <- unique(smooth.preindustrial$model)
preind.comp.out <- FAR.out <- data.frame()

for(m in 1:length(mods)){
 # m <- 1
  temp <- smooth.preindustrial %>%
    filter(model==mods[m])
  
  p.annual <- pnorm(points.sm$anomaly, mean(temp$anomaly), sd(temp$anomaly), lower.tail = F)
  #z.annual <- (points$anomaly - mean(preindustrial$anomaly)) / sd(preindustrial$anomaly)
  
  # and FAR
  
  FAR <- obs.prob <- preind.prob <- NA 
  ref <- sm.annual[21:40,3]
  
  for(i in 1:2){ # loop through each year to calculate FAR
    i <- 1
    obs.prob[i] <- sum(ref >= points.sm$anomaly[i])/20
    preind.prob[i] <- sum(temp$sm.anomaly >= points.sm$anomaly[i])/nrow(temp)
    FAR[i] <- 1-(preind.prob[i]/obs.prob[i])
  }
  
  FAR.out <- rbind(FAR.out, 
                   data.frame(model=mods[m], 
                              year=c("2014-2016", "2017-2019"),
                              anom=points.sm$anomaly,
                              obs.prob=obs.prob,
                              preind.prob=preind.prob,
                              FAR=FAR))
  
  preind.comp.out <- rbind(preind.comp.out,
                           data.frame(
                             model=mods[m],
                             p.2014.2016=p.annual[1],
                             p.2017.2019=p.annual[2]
                           ))
  
}

FAR.out <- FAR.out %>%
  arrange(model)
FAR.out

rownames(preind.comp.out) <- preind.comp.out[,1]
preind.comp.out <- preind.comp.out[,-1]

preind.comp.out <- as.data.frame(t(preind.comp.out))

rownames(preind.comp.out) <- 1:2

# and put it together for a table in the paper
FARs <- FAR.out %>%
  select(model, year, anom, FAR) %>%
  spread(model, FAR, -year, -anom)

FARs

xport <- cbind(FARs, preind.comp.out)
write.csv(xport, "smoothed anomaly preindustrial P values and FAR.csv", row.names = F)

##############################################
# old stuff below!
p.smooth <- pnorm(points.sm$anomaly, mean(smooth.preindustrial$sm.anomaly, na.rm=T),
                  sd(smooth.preindustrial$sm.anomaly, na.rm=T), lower.tail = F)
z.smooth <- (points.sm$anomaly - mean(smooth.preindustrial$sm.anomaly, na.rm = T)) / sd(smooth.preindustrial$sm.anomaly, na.rm=T)

# and export
preind.summary <- data.frame(class=c("annual", "3-yr"),
                             mean=c(mean(smooth.preindustrial$anomaly), mean(smooth.preindustrial$sm.anomaly, na.rm = T)),
                             sd=c(sd(smooth.preindustrial$anomaly), sd(smooth.preindustrial$sm.anomaly, na.rm = T)),
                             max=c(max(smooth.preindustrial$anomaly), max(smooth.preindustrial$sm.anomaly, na.rm = T)))
write.csv(preind.summary, "preindustrial summary.csv")

obs.summary <- data.frame(year=c(2014:2019, "2014-2016", "2017-2019"),
                          z = c(z.annual, z.smooth),
                          p = c(p.annual, p.smooth))

obs.summary[,2] <- round(obs.summary[,2], 2)
obs.summary[,3] <- round(obs.summary[,3], 4)

write.csv(obs.summary, "observed summary.csv")

# and calculate fraction of attributable risk

# start with annual FAR
obs.probability <- preind.probability <- NA
ref <- annual.anomaly[91:120]

for(i in 2014:2019){

  # i <- 2014
  obs.probability[(i-2013)] <- sum(ref >= annual.anomaly[names(annual.anomaly)==i])/30
  preind.probability[(i-2013)] <- sum(preindustrial$anomaly >= annual.anomaly[names(annual.anomaly)==i])/nrow(preindustrial)
}

annual.FAR <- 1-(preind.probability/obs.probability)

# now 3-yr running mean FAR
keeper.years <- seq(2018, 1990, -3)
obs.prob.sm <- preind.prob.sm <- NA
ref <- sm.obs[names(sm.obs) %in% keeper.years]

counter <- c(2015, 2018)

for(i in 1:2){
  
  # i <- 2014
  obs.prob.sm[i] <- sum(ref >= sm.obs[names(sm.obs)==counter[i]])/length(ref)
  preind.prob.sm[i] <- sum(smooth.preindustrial$anomaly >= annual.anomaly[names(annual.anomaly)==counter[i]])/nrow(smooth.preindustrial)
}

smooth.FAR <- 1-(preind.prob.sm/obs.prob.sm)

# now a SI plot
# and make a 2-panel version!
short.list.colors <- cb[c(6,4,2,7,8)]

annual <- data.frame(year=as.numeric(as.character(names(annual.anomaly))), anomaly=annual.anomaly)

sm.keep <- seq(2018, 1900, -3)

sm.keep <- sm.obs[names(sm.obs) %in% sm.keep]

annual$anomaly3 <- sm.keep[match(annual$year, names(sm.keep))]

model.means <- preindustrial %>%
  group_by(model) %>%
  summarise(max.anomaly=max(anomaly))

sm.model <- smooth.preindustrial %>%
  group_by(model) %>%
  summarize(max.anomaly3=max(sm.anomaly)) %>%
  filter(model != "MRI.CGCM3") # drop the model with unrealistic AR(1) values

plot.a <- ggplot(annual, aes(year, anomaly)) +
  theme_bw() +
  geom_hline(yintercept=c(model.means$max.anomaly), color=short.list.colors) +
  geom_line() +
  geom_point() +
  theme(axis.title.x = element_blank()) +
  ylab("SST anomaly") +
  scale_x_continuous(breaks=seq(1900, 2020, 20), lim=c(1900,2020)) +
  ggtitle("a) Annual SST") 


for(i in 1:nrow(model.means)){
  plot.a <- plot.a + annotate("text", x=2019, y=-1.5-0.3*i, label=model.means$model[i], color=short.list.colors[i], adj=1, size=3)
}

sm.annual <- na.omit(annual)
  
plot.b <- ggplot(sm.annual, aes(year, anomaly3)) +
  theme_bw() +
  geom_hline(yintercept=c(sm.model$max.anomaly3), color=short.list.colors[2:5]) +
  geom_line() +
  geom_point() +
  theme(axis.title.x = element_blank()) +
  ylab("SST anomaly") +
  scale_x_continuous(breaks=seq(1900, 2020, 20), lim=c(1900,2020)) +
  ggtitle("b) Three-year running mean SST") 

for(i in 1:nrow(sm.model)){
  plot.b <- plot.b + annotate("text", x=2019, y=-1.5-0.3*i, label=sm.model$model[i], color=short.list.colors[(i+1)], adj=1, size=3)
}

png("sst anomalies two-panel.png", 10, 3, units="in", res=300)
ggpubr::ggarrange(plot.a, plot.b)
dev.off()
