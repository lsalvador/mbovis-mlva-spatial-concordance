##########################################################################
# Code to generate Figures 2,3 and 4 of the manuscript:
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
# Description: 
# Figure 2: Distribution of M. bovis isolates by year and host-species. 
# Distribution of M. bovis isolates (n = 1136) from cattle (black bars) and 
# badgers (grey bars) by sampling years. 
#
# Figure 3. Distribution of M. bovis MLVA types by host-species a. M. bovis 
# MLVA types found only in cattle (black) and those found only in badgers 
# (grey); b. M. bovis MLVA types shared between cattle (black) and badgers (grey).
#
#Figure 4. Distribution of M. bovis MLVA types in home-bred and purchased cattle. 
# M. bovis MLVA types found only in cattle that were home-bred (black bars) and 
# purchased (grey bars) in the area.
#
#
# Note: working directory setup was deleted from original version
#
# Input:  Supplementary_Table_2.xls / see Data Availability Statement]
# Output: Figure_2.pdf
#         Figure_3.pdf
#         Figure_4.pdf
# Requirements: R 4.4.2; packages: gridExtra, ggplot2, ggbreak, patchwork
# Contact: lilianasalvador@arizona.edu
# Last updated: October 1, 2026
##########################################################################


install.packages("ggplot2")
install.packages("ggbreak")
install.packages("gridExtra")  
install.packages("patchwork")
library(gridExtra)
library(ggplot2)
library(ggbreak)
library(patchwork)


#set working directory
TVRdata <- read_excel("Supplementary_Table_2.xlsx", sheet = 1,
                      col_names = TRUE, na = c("", "NA"))

# Match read.csv output: plain data.frame with read.csv-style column names
TVRdata <- as.data.frame(TVRdata)
names(TVRdata) <- make.names(names(TVRdata), unique = TRUE)

#Rename 'Zone' and '2km' to 'Cattle' variable; Badger to Badger
TVRdata$DVO <- ifelse(TVRdata$DVO == "Badger", 'Badger', 'Cattle')

#Create data set with badgers location (one sample per badger + RTA (one sample per animal))
UniqueBadgers <- subset(TVRdata, DVO == 'Badger' & !duplicated(Tno))
#Create just cattle data table 
UniqueCattle <- subset(TVRdata, DVO == "Cattle")

#Create data frame with unique badger cattle data
obsdata <- data.frame(rbind(UniqueBadgers, UniqueCattle))
#Data table with non-shared genotypes 
obsdata1 <- subset(obsdata, Genotype %in% c(1,2,3,5,7,8,9,10,11,19,20,25,27,49,68,73,114,117,145,146,149,158,169,255,266,293,421,423,464,471,543,997))
obsdata1$col <- paste(ifelse(obsdata1$DVO =="Badger", "grey", "black"))
head(obsdata1)

#Data table with only shared genotypes
shared <- subset(obsdata, Genotype %in% c(4, 6, 297, 122))

######################################################################
####~~~~~~~~~~~~~~~~~~~~~~~Figure 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######################################################################
#********************Figure2a,b,c. a.genotypes by years, b. genotypes by Hosts, c.shared genotypes***************
#Plot samples by years
years <- data.frame(table(obsdata$DVO, obsdata$Year)/length(obsdata$No))
colnames(years) <- c("DVO", "Year", "Freq")

plot1 <- ggplot(years, aes(x = Year, y = Freq, fill = DVO)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.9) +
  scale_fill_manual(values = c("Cattle" = "black", "Badger" = "grey")) +
  labs(x = "Years",    # X-axis title
       y = "Proportion from all isolates"  # Y-axis title
  )  +     
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text=element_text(size=9),
        axis.title = element_text(size = 12),
        # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.title = element_blank()
  )

#Create genotypes table by hosts and plot
# Example data
genotypes<- c("1", "2", "3", "5", "7", "8", "9", "10", "11", "19", "20", "25", "27", "49", "68", "73", "114", "117", "145", "146", "149", "158", "169", "255", "266", "293", "421", "423", "464", "471", "543", "997")
# Genotype labels
gens <- data.frame(table(obsdata1$DVO, obsdata1$Genotype)/length(obsdata$No))
colnames(gens) <- c("DVO", "Genotype", "Freq")

# Create a data frame for plotting
data <- data.frame(
  Genotype = gens$Genotype,
  Frequency = gens$Freq,
  DVO = gens$DVO
)

# Save the combined plot to a PDF file -> Figure2
ggsave("Figure_2.pdf", plot1, width = 8, height = 6)


######################################################################
####~~~~~~~~~~~~~~~~~~~~~~~Figure 3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######################################################################

# Create the barplot with axis break and labels beside bars
plot2 <- ggplot(data, aes(x = Genotype, y = Frequency, fill = DVO)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_y_break(c(0.03, 0.13)) + # Add axis break
  scale_fill_manual(values =c("Cattle" = "black", "Badger" = "grey")) +
  labs(x = "MLVA type",    # X-axis title
       y = "Proportion from all isolates"  # Y-axis title
  )  +     
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text=element_text(size=10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.title = element_blank()
  )

#Data table with only shared genotypes
shared <- subset(obsdata, Genotype %in% c(4, 6, 297, 122))
gens1 <- data.frame(table(shared$DVO, shared$Genotype)/length(obsdata$No))
colnames(gens1) <- c("DVO", "Genotype", "Freq")

# Create the barplot with axis break and labels beside bars
plot3 <- ggplot(gens1, aes(x = Genotype, y = Freq, fill = DVO)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_y_break(c(0.2, 0.4)) +
  scale_fill_manual(values = c("Cattle" = "black", "Badger" = "grey")) +
  labs(x = "MLVA type",    # X-axis title
       y = "Proportion from all isolates"  # Y-axis title
  )  +     # Remove legend title
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text=element_text(size=10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.title = element_blank()
  )

# Combine plots using patchwork

#combine plots 2 and 3 -> Figure 3ab
combined_plot <- plot2 + plot3 +
  plot_layout( ncol=2, nrow = 1,
  widths =  c(2, 1))

# Save the combined plot to a PDF file -> Figure 3ab
ggsave("Figure_3.pdf", combined_plot, width = 10, height = 6)

######################################################################
####~~~~~~~~~~~~~~~~~~~~~~~Figure 4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######################################################################

#home-bredvsBrought-in

#create table with brougt-in and home-bred cattle genotypes
table1 <- data.frame(table(TVRdata$H.B, TVRdata$Genotype, dnn=c("Type", "Gen"))/length(TVRdata$No)) # Creates a contingency table
table1$Type <- ifelse(table1$Type == 'H', "Home-bred", "Purchased")


plot4 <- ggplot(table1, aes(x = Gen, y = Freq, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_y_break(c(0.08, 0.15)) +
  scale_fill_manual(values = c("Home-bred" = "black", "Purchased" = "darkgrey")) +
  labs(x = "MLVA type",    # X-axis title
       y = "Proportion from all isolates"  # Y-axis title
  )  +     # Remove legend title
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text=element_text(size=10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.title = element_blank()
  )

pdf("Figure_4.pdf",
    width = 8, height = 7, onefile = FALSE)
print(plot4)
dev.off()