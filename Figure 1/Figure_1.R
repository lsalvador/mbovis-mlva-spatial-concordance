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
# Output: Figure_1.jpexg
# Requirements: R 4.4.2; packages: geodata, readxl, terra
# Contact: lilianasalvador@arizona.edu
# Last updated: October 18, 2026
##########################################################################
# ---------------------------------------------------------------
# Figure 1: TVR zone, 2 km buffer, cattle and badger isolates,
# with inset locator map of Northern Ireland, legend and scale bar
# Packages needed: terra, readxl, geodata
# ---------------------------------------------------------------
install.packages(c("terra", "readxl", "geodata"))  

library(terra)
library(readxl)
library(geodata)


# ---- 1. Files ---------------------------------------------------
data_file   <- "Supplementary_Table_2.xlsx"
tvr_file    <- "TVR.shp"
buffer_file <- "2km_Buffer.shp"
out_file    <- "Figure_1.jpeg"

setGDALconfig("SHAPE_RESTORE_SHX=YES")   # rebuilds a missing .shx index

# ---- 2. Colours -------------------------------------------------
col_cattle <- "black"
col_badger <- "lightsteelblue3"

# ---- 3. Load data -----------------------------------------------
TVRdata <- as.data.frame(read_excel(data_file, sheet = 1,
                                    col_names = TRUE, na = c("", "NA")))
names(TVRdata) <- make.names(names(TVRdata), unique = TRUE)
TVRdata$DVO <- ifelse(TVRdata$DVO == "Badger", "Badger", "Cattle")

tvr <- vect(tvr_file)
buf <- vect(buffer_file)

# Northern Ireland districts from GADM (needs internet the first time)
gbr <- gadm(country = "GBR", level = 2, path = tempdir())
ni  <- gbr[gbr$NAME_1 == "Northern Ireland", ]
ni  <- project(ni, crs(tvr))

# Points: assumed to be in the same CRS as the TVR shapefile.
# If not, use crs = "EPSG:29902" (Irish Grid) or "EPSG:2157" (ITM) and project().
#pts <- vect(TVRdata, geom = c("X.coord", "y.coord"), crs = crs(tvr))
 pts <- vect(TVRdata, geom = c("X.coord", "y.coord"), crs = "EPSG:29902")
# pts <- project(pts, crs(tvr))

cattle <- pts[pts$DVO == "Cattle", ]
badger <- pts[pts$DVO == "Badger", ]

# Remove cattle points that fall outside the TVR zone + 2 km buffer
inside <- is.related(cattle, buf, "intersects") | is.related(cattle, tvr, "intersects")
cat("Cattle points removed (outside buffer):", sum(!inside), "| kept:", sum(inside), "\n")
cattle <- cattle[inside, ]

# ---- 4. Extent with empty space around the buffer ---------------
# (room on the left/bottom for the inset and scale bar, top/right for the legend)
e  <- ext(buf)
bw <- xmax(e) - xmin(e)
bh <- ymax(e) - ymin(e)
pad_left   <- 0.40   # empty space left of the buffer (room for the inset)
pad_bottom <- 0.22   # empty space below the buffer (room for the inset and scale bar)
plot_ext <- ext(xmin(e) - pad_left * bw, xmax(e) + 0.08 * bw,
                ymin(e) - pad_bottom * bh, ymax(e) + 0.04 * bh)

# ---- 5. Open JPEG device (high resolution) ----------------------
jpeg(out_file, width = 9, height = 7, units = "in", res = 600, quality = 100)

# ---- 6. Main map ------------------------------------------------
par(mar = c(0.5, 0.5, 0.5, 0.5), xpd = FALSE)
plot(buf, ext = plot_ext, col = "white", border = "black", lwd = 1.2,
     axes = FALSE, box = FALSE, legend = FALSE)
plot(tvr, col = "gray97", border = "black", lwd = 1, add = TRUE)

points(cattle, col = col_cattle, pch = 20, cex = 1.3)   # cattle first
points(badger, col = col_badger, pch = 20, cex = 1.3)   # badgers on top

u <- par("usr")
w <- u[2] - u[1]
h <- u[4] - u[3]

# Legend (top right, drawn outside the clipping region so nothing is cut)
par(xpd = NA)
legend(x = u[2] - 0.20 * w, y = u[4] - 0.01 * h, bty = "n", cex = 1.1,
       legend = c("Cattle", "Badger"),
       col    = c(col_cattle, col_badger),
       pch    = 20, pt.cex = 1.5)

# Scale bar: 4 km (map units must be metres), bottom right
sx <- u[2] - 0.27 * w     # left end of the bar
sy <- u[3] + 0.06 * h     # bottom of the bar
sh <- 0.012 * h           # bar height
rect(sx,        sy, sx + 2000, sy + sh, col = "black", border = "black")
rect(sx + 2000, sy, sx + 4000, sy + sh, col = "white", border = "black")
text(sx + 2000, sy - 0.03 * h, "4 km", cex = 1)

# ---- 7. Inset (bottom left), drawn directly in the main map's coordinates ----

inset_w   <- 0.29     # box width as a fraction of the map width (bigger = bigger box AND map)
inset_x   <- 0.02     # distance from the left edge (bigger = moves right)
inset_y   <- 0.02     # distance from the bottom edge (bigger = moves up)
inset_pad <- 0.14     # white space between coastline and box edge (smaller = bigger map in the box, 0 = touching)

ne  <- ext(ni)
ew  <- xmax(ne) - xmin(ne)
eh  <- ymax(ne) - ymin(ne)

box_w <- inset_w * w                 # box size in map units
box_h <- box_w * eh / ew             # box height follows NI's shape
xl <- u[1] + inset_x * w
yb <- u[3] + inset_y * h
padx <- inset_pad * box_w
pady <- inset_pad * box_h

# Scale NI to fit the box, then move it into place
sc   <- (box_w - 2 * padx) / ew
ni_s <- rescale(ni, fx = sc, x0 = xmin(ne), y0 = ymin(ne))
ni_s <- shift(ni_s, dx = xl + padx - xmin(ne), dy = yb + pady - ymin(ne))

rect(xl, yb, xl + box_w, yb + box_h, col = "white", border = "gray40")
polys(ni_s, col = "gray95", border = "gray55", lwd = 0.4)

# Dot marking the TVR zone
cen <- crds(centroids(tvr))

dot_override <- c(fx = 0.69, fy = 0.25)

if (is.null(dot_override)) {
  dot_x <- xl + padx + (cen[1, 1] - xmin(ne)) * sc
  dot_y <- yb + pady + (cen[1, 2] - ymin(ne)) * sc
} else {
  dot_x <- xl + padx + dot_override["fx"] * ew * sc
  dot_y <- yb + pady + dot_override["fy"] * eh * sc
}
points(dot_x, dot_y, pch = 20, col = "black", cex = 0.8)

## North arrow (top right inside the box)
##north_size <- 0.07
#north(xy = c(xl + 0.88 * box_w, yb + 0.62 * box_h), type = 1, d = 0.12 * box_h)

north_size <- 0.12     # arrow length (was 0.12)
north_head <- 0.05     # arrow head only (default 0.1) 
north_cex  <- 0.7      # size of the letter N
north_x    <- 0.88     # across the box
north_y    <- 0.80     # up the box (was 0.62)
op <- par(cex = north_cex)
north(xy = c(xl + north_x * box_w, yb + north_y * box_h), type = 1, d = north_size * box_h, head = north_head)
par(op)


dev.off()

