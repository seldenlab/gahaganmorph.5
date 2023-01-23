# compare colors
library(here)
library(colordistance)
library(recolorize)

# get image paths
# convert from rgb to CIELab
images <- colordistance::getImagePaths('la.gahagan')

# generate color distance matrix
cdm <- imageClusterPipeline(images)

library(ape)
tree <- nj(as.dist(cdm))
plot(tree, direction = "upwards", cex = 0.5)

# define add_image function:
add_image <- function(obj, # an object interpretable by rasterImage
                      x = NULL, # x & y coordinates for the center of the image
                      y = NULL,
                      width = NULL, # width of the image
                      interpolate = TRUE, # method for resizing
                      angle = 0) {
  
  # get current plotting window parameters:
  usr <- graphics::par()$usr # extremes of user coordinates in the plotting region
  pin <- graphics::par()$pin # plot dimensions (in inches)
  
  # image dimensions and scaling factor:
  imdim <- dim(obj)
  sf <- imdim[1] / imdim[2]
  
  # set the width of the image (relative to x-axis)
  w <- width / (usr[2] - usr[1]) * pin[1]
  h <- w * sf # height is proportional to width
  hu <- h / pin[2] * (usr[4] - usr[3]) # scale height to y-axis range
  
  # plot the image
  graphics::rasterImage(image = obj,
                        xleft = x - (width / 2), xright = x + (width / 2),
                        ybottom = y - (hu / 2), ytop = y + (hu/2),
                        interpolate = interpolate,
                        angle = angle)
}

# explaining NMDS is beyond the scope of this post (and is probably best left to the ecologists)
# see this link for more: 
# https://cougrstats.wordpress.com/2019/12/11/non-metric-multidimensional-scaling-nmds-in-r/

# for now, we'll just do it in two lines!
library(vegan)
nmds_scores <- scores(metaMDS(comm = as.dist(cdm)))

plot(nmds_scores)
for (i in 1:length(images)) {
  
  # read the image into R:
  img <- png::readPNG(images[i]) 
  
  # add the image:
  add_image(img, x = nmds_scores[i, 1], y = nmds_scores[i, 2],
            width = 0.025)
}

# plot the tree
plot(tree, show.tip.label = FALSE, direction = "upwards")

# get the parameters from the plotting environment
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# get the xy coordinates of the tips
ntip <- lastPP$Ntip

# first n values are the tips, remaining values are the coordinates of the
# nodes:
xy <- data.frame(x = lastPP$xx[1:ntip],
                 y = lastPP$yy[1:ntip]) 

# we can add points to the tree pretty easily using generic functions:
points(xy[ , 1], xy[ , 2], 
       col = viridisLite::viridis(40), 
       pch = 19, cex = 2)

# get image names
imnames <- tools::file_path_sans_ext(basename(images))

# get tip labels
tipnames <- tree$tip.label

# in my case, the tip labels are identical to the image names, so I can
# use these to check that my images are in the right order:
image_order <- match(tipnames, imnames)
images <- images[image_order]

# and plot!
par(mar = rep(0, 4))
plot(tree, show.tip.label = FALSE, direction = "upward")

# get the parameters from the plotting environment
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# get the xy coordinates of the tips
ntip <- lastPP$Ntip
xy <- data.frame(x = lastPP$xx[1:ntip],
                 y = lastPP$yy[1:ntip]) 

for (i in 1:length(images)) {
  add_image(png::readPNG(images[i]),
            x = xy[i, 1],
            y = xy[i, 2], 
            width = 0.5)
}
