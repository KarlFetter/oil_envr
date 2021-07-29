setwd("/media/max/HDD/My_Documents/McAssey/")

library(ggplot2)
library(ggmap)
library(gridExtra)
library(raster)

# Read in GC data and create latitude, longitude, and sample columns. Also adds row for sample not in phenotype data

GCdata <- read.csv("~/Dropbox/oil_envr/methods/max/oilMap/GC-data.csv")
latitude <- as.numeric(rep("", nrow(GCdata)))
longitude <- as.numeric(rep("", nrow(GCdata)))
sample <- rep(FALSE, nrow(GCdata))
site <- as.character(rep("", nrow(GCdata)))
GCdata <- cbind(GCdata, latitude, longitude, sample, site)
levels(GCdata$site) <- c("", "TX1", "TX2", "TX3", "CAN1", "CAN2", "CAN3")
# Add lat and long data from GRIN and select populations that were sampled. ## Ind 29273 is the same as Ind 664692 in transcriptomics data

GCdata[GCdata$Population == 1, 10] <- 31.18916667
GCdata[GCdata$Population == 1, 11] <- -103.57805556
GCdata[GCdata$Population == 1, 12] <- TRUE
GCdata[GCdata$Population == 1, 13] <- "TX2"

GCdata[GCdata$Population == 2, 10] <- 31.03972222
GCdata[GCdata$Population == 2, 11] <- -104.83027778
GCdata[GCdata$Population == 2, 12] <- TRUE
GCdata[GCdata$Population == 2, 13] <- "TX1"

GCdata[GCdata$Population == 3, 10] <- 31.27277778
GCdata[GCdata$Population == 3, 11] <- -102.69222222
GCdata[GCdata$Population == 3, 12] <- TRUE
GCdata[GCdata$Population == 3, 13] <- "TX3"

 GCdata[GCdata$Population == 4, 10] <- 35.19222222
GCdata[GCdata$Population == 4, 11] <- -102.10722222

GCdata[GCdata$Population == 5, 10] <- 35.19888889
GCdata[GCdata$Population == 5, 11] <- -100.79888889

GCdata[GCdata$Population == 6, 10] <- 35.41194444
GCdata[GCdata$Population == 6, 11] <- -99.40388889

GCdata[GCdata$Population == 7, 10] <- 41.10000000
GCdata[GCdata$Population == 7, 11] <- -97.66666667

GCdata[GCdata$Population == 8, 10] <- 41.21050000
GCdata[GCdata$Population == 8, 11] <- -101.64900000

GCdata[GCdata$Population == 9, 10] <- 41.41777778
GCdata[GCdata$Population == 9, 11] <- -104.09777778

GCdata[GCdata$Population == 10, 10] <- 46.60000000
GCdata[GCdata$Population == 10, 11] <- -108.53333333

GCdata[GCdata$Population == 11, 10] <- 46.66666667
GCdata[GCdata$Population == 11, 11] <- -105.26666667

GCdata[GCdata$Population == 12, 10] <- 46.87916667
GCdata[GCdata$Population == 12, 11] <- -102.78916667

GCdata[GCdata$Population == 13, 10] <- 50.04750000
GCdata[GCdata$Population == 13, 11] <- -104.70722222
GCdata[GCdata$Population == 13, 12] <- TRUE
GCdata[GCdata$Population == 13, 13] <- "CAN3"

GCdata[GCdata$Population == 14, 10] <- 50.39361111
GCdata[GCdata$Population == 14, 11] <- -108.48027778
GCdata[GCdata$Population == 14, 12] <- TRUE
GCdata[GCdata$Population == 14, 13] <- "CAN1"

GCdata[GCdata$Population == 15, 10] <- 50.66000000
GCdata[GCdata$Population == 15, 11] <- -105.66472222
GCdata[GCdata$Population == 15, 12] <- TRUE
GCdata[GCdata$Population == 15, 13] <- "CAN2"

Region2 <- GCdata$Region
Region2[Region2 == 1] <- "Texas"
Region2[Region2 == 5] <- "Canada"
Region2[Region2 == 2] <- "Other"
Region2[Region2 == 3] <- "Other"
Region2[Region2 == 4] <- "Other"

GCdata <- cbind(GCdata, Region2)
GCdata$sat <- GCdata$sat/100
GCdata$unsat <- GCdata$unsat/100

# Use ggmap to plot sampling locations

#register_google(key = "AIzaSyD5Si-SDQ4hjGlRFDlHT2fnLXgSmpfZJ18")

myLocation <- c(-115,27.5,-92,55.2)
myMap <- get_stamenmap(bbox = myLocation, source = "stamen", maptype = "toner-lite", zoom = 4, color = "bw", force = TRUE)


map <- ggmap(myMap) + geom_point(aes(x = longitude, y = latitude, color = Region2), data = GCdata, size = 3) + 
  scale_fill_manual(values = c("Texas" = "#fc8d59", "Canada" = "#91bfdb", "Other" = "#888888")) + 
  scale_color_manual(values = c("Texas" = "#fc8d59", "Canada" = "#91bfdb", "Other" = "#888888")) +
  #scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 19)) +
  theme_light(base_size = 20) +
  theme(legend.position = "none") +
   xlab("Longitude") +
  ylab("Latitude")
  # + geom_text(aes(x=longitude, y=latitude, label=site), data = GCdata, vjust = -1, size = 2.5)


map
# Phenotype fig

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lm_eqn <- function(y, x){
  m <- lm(y ~ x);
  pv <- lmp(m);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pv,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        pv = pv, digits = 3))
  as.character(as.expression(eq));
}

lm_eqn(GCdata$unsat, GCdata$latitude)

# Plot regression

GCdata$unsat <- GCdata$unsat * 100

sat <- ggplot(aes(x = unsat, y = latitude, color = Region2), data = GCdata) +
  geom_point() +
  theme_light(base_size = 20) +
  theme(legend.position = "none") +
  xlab("Unsaturated Fatty Acid %") +
  ylab("Latitude") +
  scale_color_manual(values = c("#91bfdb", "#888888", "#fc8d59")) +
  scale_y_continuous(breaks = seq(25, 55, by = 5), limits = c(28, 55), expand = c(0,0)) +
  geom_smooth(method = "lm", aes(group=1), color = "#555555")

# This code is in Oil_Envr_Analyses_April2021_v1.R
p <- plot_grid(r0, r2, r1, r3, r4, r5, ncol = 2, nrow = 3)

plot_grid(map, sat, p, nrow = 1, labels="auto")

#### Get Bioclim data for Karl

worldclim <- getData("worldclim",var="bio",res=2.5)

coords <- GCdata[, 11:10]
coords <- SpatialPoints(coords, proj4string = worldclim@crs)

oilClim <- extract(worldclim, coords)

BIO <- c("Annual Mean Temperature",
"Mean Diurnal Range (Mean of monthly (max temp - min temp))",
"Isothermality (BIO2/BIO7) (* 100)",
"Temperature Seasonality (standard deviation *100)",
"Max Temperature of Warmest Month",
"Min Temperature of Coldest Month",
"Temperature Annual Range (BIO5-BIO6)",
"Mean Temperature of Wettest Quarter",
"Mean Temperature of Driest Quarter",
"Mean Temperature of Warmest Quarter",
"Mean Temperature of Coldest Quarter",
"Annual Precipitation",
"Precipitation of Wettest Month",
"Precipitation of Driest Month",
"Precipitation Seasonality (Coefficient of Variation)",
"Precipitation of Wettest Quarter",
"Precipitation of Driest Quarter",
"Precipitation of Warmest Quarter",
"Precipitation of Coldest Quarter")

colnames(oilClim) <- BIO

GCdata <- cbind(GCdata, oilClim)

write.csv(GCdata, "oilBioClim.csv")
