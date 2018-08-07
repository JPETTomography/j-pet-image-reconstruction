############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

############## Example for images ##############

# load 'kde'-object containing reconstructed image
load("ReconstructedImage.RData")
# watch structure
str(image_object)
# show H-matrix
image_object$H    # or image_object[["H"]]

# assign self-describing variables
x <- image_object$eval.points[[1]]
y <- image_object$eval.points[[2]]
z <- image_object$eval.points[[3]]
intensity <- image_object$estimate # 3D array

# choose colour palette (100 grades of topographical one)
colour_palette <- topo.colors(100)
# better choice for the palette, but 'RColorBrewer' package is required
library(RColorBrewer)
r_pal <- colorRampPalette(rev(brewer.pal(11,'BrBG')))
colour_palette <- r_pal(1024)

# plot XY cross-section for pixel No 51 on Z-axis (z[51] = 18.75 cm)
image(x, y, intensity[,,51],
      col = colour_palette,
      xlab = "X [cm]",
      ylab = "Y [cm]",
      main = paste0("XY cross-section at Z = ", z[51], " cm")
      )

# plot XZ cross-section for pixel No 51 on Y-axis (y[51] = 0 cm)
image(x, z, intensity[,51,],
      col = colour_palette,
      xlab = "X [cm]",
      ylab = "Z [cm]",
      main = paste0("XZ cross-section at Y = ", y[51], " cm")
)
