##########################################################################
# Code to generate Figure 6 of the manuscript:
# "Scale- and genotype-dependent spatial concordance of Mycobacterium bovis
#  in cattle and badgers"
#
# Authors: Assel Akhmetova, Liliana C. M. Salvador, Katarina Oravcova,
# Eleanor Presho, Carl McCormick, Suzan Thompson, Lorraine Wright,
# Fraser Menzies, Nigel Trimble, Roland Harwood, Rowland R. Kao,
# Adrian Allen, Theo Pepler, Robin Skuce
#
# Code developers: Assel Akhmetova, Theo Pepler, Liliana C. M. Salvador.
#
# Description: Code to determine the spatial distribution of shared MLVA types 
# at multiple kernel density levels. Kernel density estimate contours at 25%, 
# 50%, 75% and 95% probability for the three shared M. bovis (MLVA004, MLVA006, 
# and MLVA297) the three shared M. bovis MLVA types (orange, teal and purple 
# lines). The geographical positions of bTB-affected farms (black dots) and 
# culture-confirmed bTB TVR/RTA badgers (grey dots) are displayed.
#
# Note: working directory setup was deleted from original version
#
# Input:  Supplementary_Table_2.xls / see Data Availability Statement]
# Output: Figure_6.pdf
# Requirements: R 4.4.2; packages: MASS
# Contact: lilianasalvador@arizona.edu
# Last updated: October 1, 2026
##########################################################################

install.packages("MASS")
install.packages("readxl")

library(MASS)
library(readxl)

TVRdata <- read_excel("Supplementary_Table_2.xlsx", sheet = 1,
                      col_names = TRUE, na = c("", "NA"))
  
# Match read.csv output: plain data.frame with read.csv-style column names
TVRdata <- as.data.frame(TVRdata)
names(TVRdata) <- make.names(names(TVRdata), unique = TRUE)


#change DVO to B and C instead of Badger and 2km and Zone
TVRdata$DVO <- ifelse(TVRdata$DVO == 'Badger', 'B', 'C')

#create data set with badgers location (one sample per badger + RTA (one sample per animal))
UniqueBadgers <- subset(TVRdata, DVO == 'B' & !duplicated(Tno),
                        select = c("Genotype","X.coord", "y.coord"))

#create just cattle data table
UniqueCattle <- subset(TVRdata, DVO == "C",  
                       select=c("Genotype","X.coord", "y.coord"))

#merge unique cattle and badger data together
TVRNew <- rbind(UniqueBadgers, UniqueCattle)

#create data frames with only shared genotypes data (4,6,122,297), select only coordinates
genotypes <- subset(TVRNew, Genotype %in% c(4, 6, 297),
                    select = c(Genotype, X.coord, y.coord))

colnames(genotypes) <- c("Genotype", "x", "y")


# function to get the density threshold (z-level) for a given probability contour
# this takes the already-built kde2d surface (kk) as input, which guarantees the threshold 
# is computed from the exact same surface that gets drawn.

getLevel <- function(kk, prob = 0.95) {
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

# x/y limits for plotting
x <- c(305000, 330000)
y <- c(328000, 350000)

# probability levels requested
probs <- c(0.25, 0.50, 0.75, 0.95)

# Pre-compute the kde2d surfaces once per genotype (the density surface itself
# does not depend on prob, only the contour threshold level does)
kk1 <- kde2d(genotypes[genotypes$Genotype == 4, "x"],
             genotypes[genotypes$Genotype == 4, "y"],
             lims = c(x, y))

kk2 <- kde2d(genotypes[genotypes$Genotype == 6, "x"],
             genotypes[genotypes$Genotype == 6, "y"],
             lims = c(x, y))

kk3 <- kde2d(genotypes[genotypes$Genotype == 297, "x"],
             genotypes[genotypes$Genotype == 297, "y"],
           #  h = 4000,
             lims = c(x, y))

# unique sampling locations by host
UniqueBadgers <- unique(TVRdata[TVRdata$DVO == "B", c("X.coord", "y.coord")])
UniqueCattle <- unique(TVRdata[!(TVRdata$DVO %in% "B"), c("X.coord", "y.coord")])

# Colour scheme for the three genotypes 
genotype_cols <- c("MLVA4" = "#FFA500", "MLVA6" = "#0a7e8c", "MLVA297" = "#B85FD9")

# produce the 4-panel figure
pdf("Figure_6.pdf", width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (p in probs) {
  
  level1 <- getLevel(kk1, prob = p)
  level2 <- getLevel(kk2, prob = p)
  level3 <- getLevel(kk3, prob = p)
  
  MASS::eqscplot(x = x, y = y, type = 'n',
                 xlab = "Eastings", ylab = "Northings",
                 main = paste0(p * 100, "% KDE contour"),
                 cex.lab = 1.1, cex.main = 1.2, las = 1, cex.axis = 0.9)
  
  points(UniqueCattle$X.coord, UniqueCattle$y.coord, col = "black", pch = 19, cex = 0.8)
  points(UniqueBadgers$X.coord, UniqueBadgers$y.coord, col = "darkgrey", pch = 19, cex = 0.8)
  
  contour(kk1, levels = level1, drawlabels = FALSE, add = TRUE, col = genotype_cols["MLVA4"], lwd = 2)
  contour(kk2, levels = level2, drawlabels = FALSE, add = TRUE, col = genotype_cols["MLVA6"], lwd = 2)
  contour(kk3, levels = level3, drawlabels = FALSE, add = TRUE, col = genotype_cols["MLVA297"], lwd = 2)
  
  legend("topright", names(genotype_cols), lty = 1, col = genotype_cols, cex = 0.8, bty = "n")
  legend("topleft", c("Cattle", "Badgers"), pch = 19, col = c("black", "darkgrey"), cex = 0.8, bty = "n")
}

dev.off()

