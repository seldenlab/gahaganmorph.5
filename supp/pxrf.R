## load packages
library(here)
library(tidyverse)
library(EnvStats)
library(ggpubr)
library(ggfortify)
library(cluster)
library(ggExtra)
library(RRPP)
library(wesanderson)

## load pXRF + qual.data----
raw <- read.csv('pxrf1.csv',
                header = TRUE,
                as.is = TRUE)

#palette
pal <- wes_palette("Moonrise2")

#pca
df <- raw[c(4:11)]
df <- log(df)
pch.gps <- c(1:2)[as.factor(raw$Site)]
col.gps <- pal[as.factor(raw$Site)]

## pca plot
pca <- autoplot(prcomp(df),
                data = raw,
                asp = 1,
                shape = pch.gps,
                colour = "Site",
                variance_percentage = TRUE,
                loadings = TRUE, 
                loadings.colour = 'blue',
                loadings.label = TRUE,
                loadings.label.size = 3,
                frame = TRUE,
                frame.type = 't') +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal)

ggMarginal(pca, groupColour = TRUE)


