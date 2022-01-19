# install.packages("tidyverse")
library(tidyverse)

# create a random ASV table with relative abundance values in percentages
asvTable <- matrix(c(70.93,78.58,78.72,69.24,25.53,43.85,7.49,
                     10.00,78.30,18.11,71.16,63.82,47.37,89.87),ncol=2)

# reassign names to columns and rows
colnames(asvTable) <- c("Before","After")
rownames(asvTable) <- c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7")

# check what the matrix looks like
asvTable

# melt the asvTable so that the "y" types of values are in the same column
# install.packages("reshape")
library("reshape2")
asvTable_melt <- melt(asvTable)

# assign new names to columns. Renaming works it easier to call the cols when plotting.
names(asvTable_melt) <- c("ASV", "Treatment", "RA_perc")
asvTable_melt

# create a bar plot with default colors
# install.packages("ggplot")
library(ggplot2)
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge") #, colour="black") 
  
# create a bar plot with colors from the RColorBrewer package
# install.packages("RColorBrewer")
library(RColorBrewer)

# display all palettes in the RColorBrewer package
display.brewer.all()
# display colors of the Dark2 color palette in the RColorBrewer package
display.brewer.pal(n = 8, name = 'Dark2')
# display HEX code for the colors of the Dark2 color palette in the RColorBrewer package
brewer.pal(n = 8, name = "Dark2")

# let us explore the color palettes in the RColorBrewer package in a Shiny app
# install.packages("tmaptools")
library(tmaptools)
palette_explorer()

# assign your chosen color palette to a variable instead of burying it in the ggplot code
colors4ASVbars <- get_brewer_pal("Dark2", n = 7)

# create the bar plot with the RColorBrewer color palette
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  ylab("Relative Abundance (%)") +  # add the title on y axis
  xlab("Incubation Times") +
  scale_fill_manual(values = colors4ASVbars)

# lets see what the above colors appear as they do for folks with color blindness 
# install.packages("dichromat")
library(dichromat)

# assign the RColorBrewer color palette to a variable and run the dichromat function as shown below
colors4ASVbars_dichromat <- dichromat(colors4ASVbars, type = c("deutan", "protan", "tritan"))

# plot with the new color palette
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  ylab("Relative Abundance (%)") +  # add the title on y axis
  xlab("Incubation Times") +
  scale_fill_manual(values = colors4ASVbars_dichromat)

########################## some fun color palettes ############################# 
# install.packages("wesanderson")
library(wesanderson)

# assign a variable for the Wes Anderson Color palette 
myWAcolors <- wes_palette("BottleRocket1", n = 7)

# create the bar plot with the Wes Anderson color palette
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  ylab("Relative Abundance (%)") +  # add the title on y axis
  xlab("Incubation Times") +
  scale_fill_manual(values = myWAcolors)

# copy and paste the url below in a web browser
# https://emilhvitfeldt.github.io/r-color-palettes/discrete.html
# install.packages("paletteer")
library(paletteer) 

# assign a variable for your favorite color palette from Paletteer
paletterColors <- paletteer_d("beyonce::X90")

# create the bar plot with the Wes Anderson color palette
ggplot(asvTable_melt, aes(x=Treatment, y=RA_perc, fill=factor(ASV))) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  ylab("Relative Abundance (%)") +  # add the title on y axis
  xlab("Incubation Times") +
  scale_fill_manual(values = paletterColors)

