##########################################################################
# Code to generate Figure 1 of the manuscript:
# "Scale- and genotype-dependent spatial concordance of Mycobacterium bovis
#  in cattle and badgers"
#
# Authors: Assel Akhmetova, Liliana C. M. Salvador, Katarina Oravcova,
# Eleanor Presho, Carl McCormick, Suzan Thompson, Lorraine Wright,
# Fraser Menzies, Nigel Trimble, Roland Harwood, Rowland R. Kao,
# Adrian Allen, Theo Pepler, Robin Skuce
#
# Code developers: Assel Akhmetova, Liliana C. M. Salvador.
#
# Description: Code to generate map of Northern Ireland (in grey) with the TVR 
# intervention area shown as a black point and an inset of the TVR area with 
# zoomed locations of Mycobacterium bovis isolates from cattle (black dots) and 
# badgers (light blue dots) within the TVR study area (shown in grey) and the 
# buffer zone (shown in white). Some spatial points fall outside the buffer zone; 
# these represent relevant historical M. bovis isolates that were added to the 
# dataset.
#
# Note: working directory setup was deleted from original version
#
# Input:  Supplementary_Table_2.xls / see Data Availability Statement]
#         Shape files TVR and buffer areas/TVR.shp
#         Shape files TVR and buffer areas/2km_Buffer.shp
#   
# Output: Figure_1.pdf
# Requirements: R 4.4.2; packages: sp, rgdal, readxl, sf, maptools, terra
# Contact: lilianasalvador@arizona.edu
# Last updated: October 1, 2026
##########################################################################

install.packages("sp")
install.packages("rgdal") #outdated 
install.packages("readxl")
library(sp)
library(maptools)
require(sf)
require(terra)
library(readxl)
install.packages("~/Downloads/maptools_1.1-8.tar.gz", repos = NULL, type = "source")

#data load
TVRdata <- read_excel("Supplementary_Table_2.xlsx", sheet = 1,
                      col_names = TRUE, na = c("", "NA"))

# Match read.csv output: plain data.frame with read.csv-style column names
TVRdata <- as.data.frame(TVRdata)
names(TVRdata) <- make.names(names(TVRdata), unique = TRUE)

#Rename 'Zone' and '2km' to 'Cattle' variable; Badger to Badger
TVRdata$DVO <- ifelse(TVRdata$DVO == "Badger", 'Badger', 'Cattle')

# Read the shapefile
shapefile_path <- "Shape files TVR and buffer areas/TVR.shp"
#spatial_data <- st_read(shapefile_path)
#head(spatial_data)

buffer_path <- "Shape files TVR and buffer areas/2km_Buffer.shp"
#buffer_data <-st_read(buffer_path)
#head(buffer_data)

# Load the shapefiles
# Read the shapefile using terra
shapefile_data <- vect(shapefile_path)
buffer_data <- vect(buffer_path)
par(mar = c(4, 4, 2, 6))
# Plot the shapefile data 1-2
plot(buffer_data, col = "white", border = "black", axes = FALSE, box = FALSE)
plot(shapefile_data, col = "gray97", border = "black", add=TRUE)

# Basic plot of the shapefile
point_colors <- c("lightsteelblue3", "black")[as.factor(TVRdata$DVO)]
points(TVRdata$X.coord, TVRdata$y.coord, col = point_colors, pch = 20)

# Add a legend for the points
legend("topright", bty = "n",
       legend = c("Cattle", "Badger"), 
       col = c("black", "lightsteelblue3"), 
       pch = 20, inset = c(1, 0))

