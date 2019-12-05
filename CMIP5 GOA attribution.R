library(tidyverse)
library(ncdf4)
library(maps)
library(mapdata)
library(fields)
library(chron)

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
```

We will compare the preindustrial estimates above with ERSSTv5 observations. The observations are for the same area (50º-60ºN, 150º-130ºW). Here is average SST for that area from ERSSTv5:
  
  ```{r, message=F, fig.width=5, fig.height=4}
# identify latest year and month needed
year <- 2018
month <- "12"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv4.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(50):1:(60)][(210):1:(230)]", sep="")

download.file(URL, "GOA.box.ersst")

# process
nc <- nc_open("GOA.box.ersst")

# extract dates
raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
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

```

\pagebreak

And here are the annual anonamlies wrt 1981-2010. Just to make sure that I'm not making any errors, I'll outline the steps I'm using:

1. Calculate annual mean SST for each cell.
2. Change these into anomalies: subtract 1981-2010 mean, divide by 1900-2018 sd.
2. Calculate annual mean anomaly for the entire box, weighting each cell by its area.

The resulting time series looks like this:

```{r, fig.height=3}
# get anomaly for 1981:2010
yr <- as.numeric(as.character(years(d)))
m <- months(d)

f <- function(x) tapply(x, yr, mean)

annual.sst <- apply(SST, 2, f) # annual mean for each cell

# get anomaly for each cell
mu <- apply(annual.sst[rownames(annual.sst) %in% 1981:2010,], 2, mean)	# Compute monthly means for each cell during 1981:2010

# change suggested by Nick Bond - using SD for the entire time series when calculating anomalies
# doing this b.c my previous anomalies appear to be too large when compared with Walsh et al. 2018
sd <- apply(annual.sst[rownames(annual.sst) >= 1900,], 2, sd)

mu <- matrix(mu, nrow=nrow(annual.sst), ncol=length(mu), byrow=T)
sd <- matrix(sd, nrow=nrow(annual.sst), ncol=length(sd), byrow=T) 

anom <- (annual.sst - mu)/sd # now each value is a cell-specific normalized anomaly relative to 1981-2010 

weights <-  sqrt(cos(lat*pi/180))
f <- function(x) weighted.mean(x, w=weights, na.rm=T)
annual.anomaly <- apply(anom, 1, f)

annual.anomaly <- annual.anomaly[names(annual.anomaly) >=1900]

anomaly.plot <- data.frame(year=1900:2018, anomaly=annual.anomaly, sign=ifelse(annual.anomaly>0, "positive", "negative"))
anomaly.plot$sign <- reorder(anomaly.plot$sign, desc(anomaly.plot$sign))
ggplot(anomaly.plot, aes(year, anomaly, fill=sign)) +
theme_bw() +
geom_col() +
theme(legend.position='none')

```

Now going ahead with the comparison to preindustrial simulations.

I'm doing this in *two ways*.

First way ("debiasing"):
  
  1. Average annual anomalies for the historical period (1987-2005) for each model.
2. Subract the mean anomaly from the ERSSTv5 version above to calculate the "bias" for each model, i.e., bias = model.mean - observed.mean
3. Preindustrial simulations are then adjusted with this bias value, i.e.: "debiased" value = simulation anomaly - bias.

```{r}
# now compare thes observations with each model!
compare.dat <- dat %>%
  filter(Era=="historical / RCP8.5", Year <=2005) %>%
  mutate(observed=rep(annual.anomaly[names(annual.anomaly) %in% 1987:2005],5))

model.means <- tapply(compare.dat$anomaly, compare.dat$model, mean)
observed.mean <- mean(annual.anomaly[names(annual.anomaly) %in% 1987:2005])

bias <- model.means - observed.mean 

# so...we subtract this bias from each model's preindustrial simulations
# to calculate unbiased estimates of preindustrial distributions for comparison with observations

envelope <- dat %>%
  filter(Era=="preindustrial") %>%
  mutate(debiased=anomaly-bias[match(model, names(bias))])

envelope$debiased.3.yr <- NA

mods <- unique(envelope$model)

for(i in 1:length(mods)){
  # i <- 1
  temp <- envelope %>%
    filter(model==mods[i])
  
  envelope$debiased.3.yr[envelope$model==mods[i]] <- zoo::rollmean(temp$debiased, 3, fill=NA)
  
  envelope$anomaly.3.yr[envelope$model==mods[i]] <- zoo::rollmean(temp$anomaly, 3, fill=NA)
}

fmin <- function(x) min(x, na.rm=T)
fmax <- function(x) max(x, na.rm=T)

ranges <- data.frame(min.anomaly=fmin(envelope$anomaly),
                     max.anomaly=fmax(envelope$anomaly),
                     min.anomaly3=fmin(envelope$anomaly.3.yr),
                     max.anomaly3=fmax(envelope$anomaly.3.yr),
                     min.debiased=fmin(envelope$debiased),
                     max.debiased=fmax(envelope$debiased),
                     min.debiased3=fmin(envelope$debiased.3.yr),
                     max.debiased3=fmax(envelope$debiased.3.yr))
```

The mean historical values for each model:
  
  ```{r}
model.means
```

And the mean ERSSTv5 anomaly for the same period:
  
  ```{r}
observed.mean
```

Which gives us the following bias values:
  
  ```{r}
bias
```

\pagebreak

I then calculate the envelope of estimated preindustrial values as the range of debiased anomaly values for all five models. 

I also calculate three-year running mean anomaly values for each model and plot the range of these. The idea here is to evaluate the likelihood of a multi-year anomaly under preindustrial conditions. *This step rests on the assumption that the 60 year periods from the model simulation were selected as a block, and therefore include the autocorrelation values from the original.*
  
  Here's the resulting plot. Annual observations in black, 3-year running means in blue, envelope of annual preindustrial values in grey shade, envelope of 3-year preindustrial means in blue shade. (2019 estimated from data through August)

```{r, fig.height=3, warning=F}
# quickly (!) estimate 2019 annual anomaly based on data avaliable so far!
# identify latest year and month needed
year <- 2019
month <- "8"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(50):1:(60)][(210):1:(230)]", sep="")

# download.file(URL, "GOA.box.ersst.latest")

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

yr <- as.numeric(as.character(years(d)))
m <- months(d)
weights <-  sqrt(cos(lat*pi/180))
ff <- function(x) weighted.mean(x, w=weights, na.rm=T)
weighted.mean <- apply(SST, 1, ff)

annual.sst <- tapply(weighted.mean, yr, mean)

# limit all years to the range of months available in this year!
m <- as.numeric(m)
keep <- m <= as.numeric(month) # this is the most recent month specified for data query above

short.sst <- SST[keep,]
short.weighted.mean <- apply(short.sst, 1, ff)

short.sst <- tapply(short.weighted.mean, yr[keep], mean)

# interesting, looks a little curved...
x <- as.vector(short.sst[names(short.sst) %in% 1900:2018])
y <- annual.anomaly
md1 <- lm(y ~ x)

# now add 2019 estimate to the annual anoamly TS
estimated.2019 <- short.sst[names(short.sst)==2019]*coef(md1)[2] + coef(md1)[1]
annual.anomaly <- c(annual.anomaly, estimated.2019)
names(annual.anomaly) <- 1900:2019

# and plot...
annual <- data.frame(year=names(annual.anomaly), anomaly=annual.anomaly, 
anomaly3=zoo::rollmean(annual.anomaly, 3, fill=NA))
annual$year <- as.numeric(as.character(annual$year))

# set pallette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(filter(annual, year >= 1950), aes(year, anomaly)) +
theme_bw() +
geom_line() +
geom_point() +
geom_line(aes(year, anomaly3), color=cb[3]) +
geom_ribbon(aes(ymin=ranges$min.debiased, ymax=ranges$max.debiased), alpha=0.2) +
geom_ribbon(aes(ymin=ranges$min.debiased3, ymax=ranges$max.debiased3), alpha=0.2, fill=cb[3]) +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly") +
ggtitle("Observations v. debiased envelope")

```

And here is the same plot, using the actual preindustrial envelope, without the debiasing described above. 
```{r, fig.height=3, warning=F}
# compare with non-debiased estimates...
ggplot(annual, aes(year, anomaly)) +
theme_bw() +
geom_line() +
geom_point() +
geom_line(aes(year, anomaly3), color=cb[3]) +
geom_ribbon(aes(ymin=ranges$min.anomaly, ymax=ranges$max.anomaly), alpha=0.2) +
geom_ribbon(aes(ymin=ranges$min.anomaly3, ymax=ranges$max.anomaly3), alpha=0.2, fill=cb[3]) +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly") +
ggtitle("Observations v. raw envelope")

<<<<<<< HEAD
# and save a png version
plot <- ggplot(filter(annual, year >= 1900), aes(year, anomaly)) +
=======
# export without title
xtra <- data.frame(year=1880:2040, 
min.anomaly=ranges$min.anomaly,
max.anomaly=ranges$max.anomaly,
min.anomaly3=ranges$min.anomaly3,
max.anomaly3=ranges$max.anomaly3)

annual <- left_join(xtra, annual)

save.plot <- ggplot(annual, aes(year, anomaly)) +
>>>>>>> 51c4ed52122777e59d1ad48d9485bd49b8e2976e
theme_bw() +
geom_line() +
geom_point() +
geom_line(aes(year, anomaly3), color=cb[3]) +
<<<<<<< HEAD
geom_ribbon(aes(ymin=ranges$min.anomaly, ymax=ranges$max.anomaly), alpha=0.2) +
geom_ribbon(aes(ymin=ranges$min.anomaly3, ymax=ranges$max.anomaly3), alpha=0.2, fill=cb[3]) +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly")
ggsave(plot=plot, "goa sst observations and raw envelope.png", width=6, height=4, units="in")
=======
geom_ribbon(aes(ymin=annual$min.anomaly, ymax=annual$max.anomaly), alpha=0.2) +
geom_ribbon(aes(ymin=annual$min.anomaly3, ymax=annual$max.anomaly3), alpha=0.2, fill=cb[3]) +
scale_x_continuous(breaks=seq(1900, 2020, 20), limits = c(1900,2020)) +
theme(axis.title.x = element_blank()) +
ylab("SST anomaly") 

ggsave(plot=save.plot, device="png", "GOA sst vs. raw observations.png", width=6, height=3, units="in")
>>>>>>> 51c4ed52122777e59d1ad48d9485bd49b8e2976e
```

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

