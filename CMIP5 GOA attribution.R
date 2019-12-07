library(tidyverse)
library(ncdf4)
library(maps)
library(mapdata)
library(fields)
library(chron)
library(zoo)

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

# and group them all
ggplot(preindustrial, aes(anomaly)) +
  theme_bw() +
  geom_density(alpha=0.3) +
  xlim(-5,5)

# and...run the three-year rolling means on preindustrial simulations and plot that distribution
preindustrial$sm.anomaly <- NA




# We will compare the preindustrial estimates above with ERSSTv5 observations. The observations are for the same area (50º-60ºN, 150º-130ºW). Here is average SST for that area from ERSSTv5:

# quickly (!) estimate 2019 annual anomaly based on data avaliable so far!
# identify latest year and month needed
year <- 2019
month <- "10"

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
ggplot(anomaly.plot, aes(year, anomaly, fill=sign)) +
  theme_bw() +
  geom_col() +
  theme(legend.position='none')

##################
# now estimate 2019 value
# limit all years to the range of months available in this year!
m <- as.numeric(m)
keep <- m <= as.numeric(month) # this is the most recent month specified for data query above

short.sst <- SST[keep,]
short.weighted.mean <- apply(short.sst, 1, ff)

short.sst <- tapply(short.weighted.mean, yr[keep], mean)

# interesting, looks a little curved...
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
                           sm.anom=rollmean(anomaly.plot$anomaly, 3, fill = NA))

anomaly.plot$sign <- reorder(anomaly.plot$sign, desc(anomaly.plot$sign))


anomaly.plot$sm.anom <- rollmean(anomaly.plot$anomaly, 3, fill = NA)

plot1 <- ggplot(anomaly.plot, aes(year, anomaly, fill=sign)) +
  theme_bw() +
  geom_col() +
  theme(legend.position='none')

anomaly.plot <- data.frame(year=1900:2019, anomaly=annual.anomaly,
                           sm.anom = rollmean(anomaly.plot$anomaly, 3, fill = NA))
  
plot2 <- plot1 +
  geom_path(aes(year, sm.anom))

# Now going ahead with the comparison to preindustrial simulations.

# plot preindustrial distributions with recent
ggplot(preindustrial, aes(anomaly, fill=model)) +
theme_bw() +
  geom_density(alpha=0.3) +
  xlim(-5,5) +
  geom_vline(xintercept = annual.anomaly[names(annual.anomaly) %in% 2016:2019], lty=2)


Conclusion: the debiasing is more conservative, presumably because it removes internal variability in the recent observations that is not captured by the preindustrial simulations. The raw envelope produces a conclusion more similar to Walsh et al., with a single year outside the envelope of preindustrial conditions.

Just to confirm that my treatment of three-year rolling means is appropriate, I'll calculate AR(1) values of each preindustrial simulation. Here are the model AR(1) values for 1987-2005:
  
  ```{r}

f <- function(x) ar(x, aic=F, order.max = 1)$ar

model.ar.1 <- tapply(compare.dat$anomaly, compare.dat$model, f)

observed.ar.1 <- f(annual.anomaly[names(annual.anomaly) %in% 1987:2005])

model.ar.1
```

And the observed AR(1):
  
  ```{r}
observed.ar.1
```

So it appears that the simulations were extracted as 60-year blocks, with the autocorrelation intact, and so the 3-yr approach that I use is appropriate.

Now, out of curiosity, I'll calculate the envelope values as the multi-model mean max/min.

Here are the debiased and raw plots using this approach:

```{r,fig.height=3, warnings=F}
model.means <- envelope %>%
group_by(model) %>%
summarize(min.anomaly=fmin(anomaly),
max.anomaly=fmax(anomaly),
min.anomaly3=fmin(anomaly.3.yr),
max.anomaly3=fmax(anomaly.3.yr),
min.debiased=fmin(debiased),
max.debiased=fmax(debiased),
min.debiased3=fmin(debiased.3.yr),
max.debiased3=fmax(debiased.3.yr))

multi.model.mean <- as.data.frame(t(colMeans(model.means[,-1])))

ggplot(annual, aes(year, anomaly)) +
theme_bw() +
geom_line() +
geom_point() +
geom_line(aes(year, anomaly3), color=cb[3], size=1) +
geom_ribbon(aes(ymin=multi.model.mean$min.debiased, ymax=multi.model.mean$max.debiased), alpha=0.2) +
geom_ribbon(aes(ymin=multi.model.mean$min.debiased3, ymax=multi.model.mean$max.debiased3), alpha=0.2, fill=cb[3]) +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly") +
scale_x_continuous(breaks=seq(1900, 2020, 20), lim=c(1900,2020)) +
ggtitle("Observations v. debiased envelope (multi-model mean)")

ggplot(annual, aes(year, anomaly)) +
theme_bw() +
geom_line() +
geom_point() +
geom_line(aes(year, anomaly3), color=cb[3], size=1) +
geom_ribbon(aes(ymin=multi.model.mean$min.anomaly, ymax=multi.model.mean$max.anomaly), alpha=0.2) +
geom_ribbon(aes(ymin=multi.model.mean$min.anomaly3, ymax=multi.model.mean$max.anomaly3), alpha=0.2, fill=cb[3]) +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly") +
scale_x_continuous(breaks=seq(1900, 2020, 20), lim=c(1900,2020)) +
ggtitle("Observations v. raw envelope (multi-model mean)")
```

Try plotting the range of min/max sst values for each model! 

```{r, warning=F, fig.height=3}
model.means <- envelope %>%
group_by(model) %>%
summarize(min.anomaly=fmin(anomaly),
max.anomaly=fmax(anomaly),
min.anomaly3=fmin(anomaly.3.yr),
max.anomaly3=fmax(anomaly.3.yr),
min.debiased=fmin(debiased),
max.debiased=fmax(debiased),
min.debiased3=fmin(debiased.3.yr),
max.debiased3=fmax(debiased.3.yr))

plot.model.names <- model.means$model[order(-model.means$max.anomaly)]

model.means <- model.means %>%
arrange(desc(max.anomaly))

ggplot(annual, aes(year, anomaly)) +
theme_bw() +
geom_line() +
geom_point() +
geom_line(aes(year, anomaly3), color=cb[3], size=1) +
geom_hline(yintercept=c(model.means$max.anomaly), alpha=1, lty=2) +
geom_hline(yintercept=c(model.means$max.anomaly3), alpha=1, lty=2, color=cb[3]) +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly") +
scale_x_continuous(breaks=seq(1900, 2020, 20), lim=c(1900,2020)) +
ggtitle("Observations v. raw envelope (individual model envelopes)")

# and make a 2-panel version!
short.list.colors <- cb[c(6,4,2,7,8)]

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
plot.a <- plot.a + annotate("text", x=2019, y=-1.5-0.3*i, label=plot.model.names[i], color=short.list.colors[i], adj=1, size=3)
}

# now drop MRI.CGCM3 as it has unrealistically low AR(1)
new.means <- model.means[c(1,2,4,5),]
new.names <- plot.model.names[c(1,2,4,5)]
new.colors <- short.list.colors[c(1,2,4,5)]

plot.b <- ggplot(annual, aes(year, anomaly3)) +
theme_bw() +
geom_hline(yintercept=c(new.means$max.anomaly3), color=new.colors) +
geom_line() +
# geom_point() +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly") +
scale_x_continuous(breaks=seq(1900, 2020, 20), lim=c(1900,2020)) +
ggtitle("b) Three-year running mean SST") 

for(i in 1:nrow(new.means)){
plot.b <- plot.b + annotate("text", x=2019, y=-1.5-0.3*i, label=new.names[i], color=new.colors[i], adj=1, size=3)
}

png("sst anomalies two-panel.png", 10, 3, units="in", res=300)
ggpubr::ggarrange(plot.a, plot.b)
dev.off()
```

Now...look at the % of observations that are outside the preindustrial observations for each model.

```{r}

outside.envelope <- matrix(nrow=2, ncol=5)
dimnames(outside.envelope) <- list(c("annual", "3-yr"), model.means$model)

for(j in 1:nrow(model.means)){

outside.envelope[1,j] <- sum(na.omit(annual$anomaly) > model.means$max.anomaly[j])
outside.envelope[2,j] <- sum(na.omit(annual$anomaly3) > model.means$max.anomaly3[j])

}

outside.envelope
```

