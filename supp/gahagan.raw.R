# is gahagan biface morphology correlated with geochemistry/raw material?

# load----
## load packages
library(here)
library(tidyverse)
library(geomorph)
library(wesanderson)

## load pXRF + qual.data----
raw <- read.csv('pxrf1.csv',
                header = TRUE,
                as.is = TRUE)

qual.data <- read.csv("qdata.raw.csv", 
                  header = TRUE,
                  as.is = TRUE)

# join spreadsheets with qual.data on left
qdata <- left_join(qual.data, raw,
                   by = 'specimen')

# output composite csv for GM analysis
write.csv(qdata, 
          file = 'qdata.csv', 
          row.names = FALSE)

## load GM data----
source('readmulti.csv.R')

# read .csv files
setwd("./data")
filelist <- list.files(pattern = ".csv")
coords <- readmulti.csv(filelist)
setwd("../")

# read qualitative data
qdata <- read.csv("qdata.csv", 
                  header = TRUE, 
                  row.names = 1)
qdata <- qdata[match(dimnames(coords)[[3]], rownames(qdata)),]

# GPA----
Y.gpa <- gpagen(coords,
                PrinAxes = TRUE,
                print.progress = FALSE)
## plot gpa
plot(Y.gpa)

# geomorph data frame
gdf <- geomorph.data.frame(shape = Y.gpa$coords,
                           size = Y.gpa$Csize,
                           site = qdata$sitename) 

# PCA----
pca <- gm.prcomp(Y.gpa$coords)
summary(pca)

# set plot parameters to plot by context
pal <- wes_palette("Moonrise2")
pch.gps <- c(15,17)[as.factor(qdata$sitename)]
col.gps <- pal[as.factor(qdata$sitename)]
col.hull <- c("#C27D38","#798E87")

## plot pca----
pc.plot <- plot(pca, asp = 1,
                pch = pch.gps,
                col = col.gps)
shapeHulls(pc.plot,
           groups = qdata$sitename,
           group.cols = col.hull)

# plot x/y maxima/minima
# x - minima
mean.shape <- mshape(Y.gpa$coords)
plotRefToTarget(pca$shapes$shapes.comp1$min, 
                mean.shape)
# x - maxima
plotRefToTarget(pca$shapes$shapes.comp1$max, 
                mean.shape)

# y - minima
plotRefToTarget(pca$shapes$shapes.comp2$min, 
                mean.shape)

# y - maxima
plotRefToTarget(pca$shapes$shapes.comp2$max, 
                mean.shape)

# 2B-PLS----
# geochem data w/o Ca/Fe 
geochem <- qdata %>%
  select(Al.K12:Ba.L1, Ni.K12:Zn.K12)

## shape/log.geochem----
# 2B-PLS (geochem.2) [shape/geochem]
cor <- two.b.pls(A1 = Y.gpa$coords,
                 A2 = geochem,
                 print.progress = FALSE,
                 iter = 9999)
summary(cor)

## plot PLS
plot(cor)
