library(tidyverse)

# create a random ASV table with relative abundance values in percentages
asvTable <- matrix(c(70.93,78.58,78.72,69.24,25.53,43.85,7.49,
                     10.00,78.30,18.11,71.16,63.82,47.37,89.87),ncol=2)
colnames(asvTable) <- c("Before","After")
rownames(asvTable) <- 1:7
rownames(asvTable) <- c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7")
asvTable

# melt the asvTable 
library("reshape2")
asvTable_melt <- melt(asvTable)
names(asvTable_melt) <- c("ASV", "Treatment", "RA_perc")
asvTable_melt

# create a bar plot with default colors
library(ggplot2)
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge", colour="black") 
  
# create a bar plot with colors from the RColorBrewer package
# install.packages("RColorBrewer")
library(RColorBrewer)
library(tmaptools)
library(dichromat)

# display all palettes in the RColorBrewer package
display.brewer.all()
# display colors of the Dark2 color palette in the RColorBrewer package
display.brewer.pal(n = 8, name = 'Dark2')
# display HEX code for the colors of the Dark2 color palette in the RColorBrewer package
brewer.pal(n = 8, name = "Dark2")

# let us explore the color palettes in the RColorBrewer package in a Shiny app
palette_explorer()

# assign your chosen color palette to a variable instead of burying it in the ggplot code
colors4ASVbars <- dichromat(get_brewer_pal("Dark2", n = 7))

# create the bar plot with the RColorBrewer color palette
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  ylab("Relative Abundance (%)") +  # add the title on y axis
  xlab("Incubation Times") +
  scale_fill_manual(values = colors4ASVbars)

########################## some fun color palettes ############################# 
#https://rforpoliticalscience.com/2020/07/26/make-wes-anderson-themed-graphs-with-wesanderson-package-in-r/
install.packages("wesanderson")
library(wesanderson)
myWAcolors <- wes_palette("BottleRocket1", n = 7)


# create the bar plot with the Wes Anderson color palette
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  ylab("Relative Abundance (%)") +  # add the title on y axis
  xlab("Incubation Times") +
  scale_fill_manual(values = myWAcolors)

# install.packages("paletteer")
library(paletteer) 
paletteer_d("beyonce::X4")
paletteer_d("colorBlindness::paletteMartin")
