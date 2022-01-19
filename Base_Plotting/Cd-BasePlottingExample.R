#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Figures cribbed from:
##  Anderegg et al. 2020 
## "Aridity drives coordinated trait shifts but not decreased trait variance across the geographic range of eight Australian trees"
## New Phytologist
## For UCSB EEMB R Seminar Win 22
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This document reproduces the analyses and figures presented in Anderegg et al. 2020 New Phytolotist
# using the data currently available in the 'trait_data' directory in this github repository
# For questions or to report bugs, please contact Leander Anderegg: leanderegg@gmail.com
# last updated: 03 June 2020


# set working directory if needed:
setwd("/Users/leeanderegg/Desktop/UCSB Classes/W2022 R Seminar/Base Plotting")
# create a directory for generated results
results_dirname <- "Results_YYYYMMDD" # "EucTraits_Github_Repo/Results_20200603
dir.create(results_dirname)



################# Load Required Packages and Create Functions ###################3
require(lme4)
require(lmerTest)
require(nlme)
require(RColorBrewer)
require(dplyr)
require(reshape2)
require(tidyr)
require(lmodel2)
#require(car)
#require(MuMIn)
#source("ggplot_helpers.R") # ggplot helper functions to get rid of pesky background

# set a nice palette
pal <- brewer.pal(9, name = "Set1")
pallight <- paste0(pal, "77")
palette(pallight)



# create a function for plotting Standard Major Axis regression lines (Type II regressions) to data, fitted using lmodel2
plot.MAR <- function(xvar, yvar, data, method="SMA", linecol, lwd=1, lty=1) {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    #return(rep(NA, times=7))
    break()
  }
  else{
    if(var(data[,yvar], na.rm=T)==0){
      break()
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      yhat <- tmp.mod$x * slope + intercept
      lines(yhat[order(tmp.mod$x)]~tmp.mod$x[order(tmp.mod$x)], col=linecol, lwd=lwd, lty=lty)
      
    }
  }
}


# to creat dataframes needed for plotting, run everything between
###### RUN EVERYTHING BETWEEN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: LOAD DATA #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# import soil data from the Soil and Landscape Grid of Australia (Grundy et al. 2015)
soilsums <- read.csv("data/Soils_summaries0-60cm_20200701.csv") [,-1]
# this was created with the 'Cd-SoilsExtraction.R' code, also in this repo


# import raw branch-level traits
traits.b0 <- read.csv("data/TraitsAll_branch_20200701.csv", header=T, row.names = 1 )
# NOTE: these are raw data with problematic points identified by the 'Flag_XXX' columns (outliers, young leaves, lost or damaged samples, etc identified with flag >0)
# kill bad trait values so I don't always have to filter them.
traits.b0$LDMC[which(traits.b0$Flag_LDMC>0)] <- NA
traits.b0$LMA[which(traits.b0$Flag_LDMC>0)] <- NA
traits.b0$Al_As[which(traits.b0$Flag_Area>0)] <- NA
traits.b0$LMA[which(traits.b0$Flag_Area>0)] <- NA

# remove three additional outliers identified during analysis. Results are qualitatively robust to removal of these outliers.
#WD of ACAC-PER-B-5-1 is a weird low outlier that needs removed
traits.b0$WD[which(traits.b0$Branchtag=="ACAC-PER-B-5-1")] <- NA
#LMA  "VIMI-FREY-C-1-c" is also weird low outlier
traits.b0$LMA[which(traits.b0$Branchtag == "VIMI-FREY-C-1-c")] <- NA
# bad ESAL tree for LMA, but LDMC and Al_As seem to be OK
traits.b0$LMA[which(traits.b0$Treetag=="ESAL-BEN-C-1")] <- NA
# calculate Huber value instead of Al_As
traits.b0$hub <- 1/traits.b0$Al_As
# add in scaled climate for use with trait-climate models
traits.b0 <- traits.b0 %>% group_by(Species) %>% mutate(MDc_scaled = scale(MDc), PPTc_scaled = scale(PPTc), PETc_scaled=scale(PETc))

# merge soil data
soils.b <- soilsums[match(traits.b0$Plottag, soilsums$Plottag),-1]
traits.b <- data.frame(traits.b0, soils.b)
# renaming soil for residplots function later on
traits.b$pPC1fert <- traits.b$PC1fert
traits.b$pPC2depth <- traits.b$PC2depth
# ACAC has no DBHs (multistemed in general), so making a dummy value of 0 so my trait-climate model fitting function works for ACAC
traits.b$TreeDBH[which(traits.b$Species=="ACAC")] <- 0 # fill with dummy variable so my function works for fitting trait-env models



# aggregate traits to tree average 
# (need this for most plotting, which has too much overplotting at branch level)
traits.t0 <- data.frame(traits.b  %>% group_by (Species,Site,Plot,Tree,Plottag,Sitetag,Treetag,TreeDBH) %>% summarise(
  MDc = mean(MDc), MATc = mean(MATc),
  PETc=mean(PETc), PPTc=mean(PPTc), Tminc = mean(Tminc), Tmaxc = mean(Tmaxc),
  Pminqc = mean(Pminqc), rhc = mean(rhc),
  Lat = mean(Lat), Lon= mean(Lon), StandBA = mean(StandBA, na.rm=T),
  mVolume = mean(Volume, na.rm=T),  sdWD = sd(WD, na.rm=T), WD = mean(WD, na.rm=T),
  sdLDMC = sd(LDMC, na.rm=T), LDMC = mean(LDMC, na.rm=T),
  sdLMA = sd(LMA, na.rm=T), LMA=mean(LMA, na.rm=T),
  avg_nleaves = mean(nleaves, na.rm=T), avg_tot_area = mean(tot_area, na.rm=T),
  median_area = mean(median_area, na.rm=T),
  max_area = max(max_area, na.rm=T), min_area = min(min_area, na.rm=T),
  sdAl_As = sd(Al_As, na.rm=T), Al_As = mean(Al_As, na.rm=T),
  sdhub = sd(hub, na.rm=T), hub = mean(hub, na.mr=T)
)
)
traits.t0$pchs <- rep(16, times=length(traits.t0$Species))
traits.t0$pchs[which(traits.t0$Species=="ACAC")]<- 3
# add in soils data
soils.t <- soilsums[match(traits.t0$Plottag, soilsums$Plottag),-1]
traits.t <- data.frame(traits.t0, soils.t)



# aggregate to plot level
traits.p0 <- data.frame(traits.t %>% group_by (Species,Site,Plot,Plottag,Sitetag) %>% summarise(
  MDc = mean(MDc), MATc = mean(MATc),
  PETc=mean(PETc), PPTc=mean(PPTc), Tminc = mean(Tminc), Tmaxc = mean(Tmaxc),
  Pminqc = mean(Pminqc), rhc = mean(rhc),
  Lat = mean(Lat), Lon= mean(Lon), StandBA = mean(StandBA, na.rm=T), meanDBH = mean(TreeDBH, na.rm=T),
  avg_nleaves = mean(avg_nleaves, na.rm=T),  max_area = max(max_area, na.rm=T), min_area=min(min_area, na.rm=T),
  sdWD = sd(WD, na.rm=T), medWD = median(WD, na.rm=T), WD = mean(WD, na.rm=T),
  sdLMA = sd(LMA, na.rm=T), medLMA = median(LMA, na.rm=T), LMA = mean(LMA, na.rm=T),
  sdLDMC = sd(LDMC, na.rm=T), medLDMC = median(LDMC, na.rm=T), LDMC = mean(LDMC, na.rm=T),
  sdAl_As = sd(Al_As, na.rm=T), medAl_As = median(Al_As, na.rm=T), Al_As=mean(Al_As, na.rm=T),
  sdhub = sd(hub, na.rm=T), medhub = median(hub, na.rm=T), hub = mean(hub, na.rm=T))
  
)
# add in soils data
soils.p <- soilsums[match(traits.p0$Plottag, soilsums$Plottag),-1]
traits.p <- data.frame(traits.p0, soils.p)


# aggregate to site level
traits.s <- data.frame(traits.p %>% group_by (Species,Site,Sitetag) %>% summarise(
  MDc = mean(MDc), MATc = mean(MATc),
  PETc=mean(PETc), PPTc=mean(PPTc), Tminc = mean(Tminc), Tmaxc = mean(Tmaxc),
  Pminqc = mean(Pminqc), rhc = mean(rhc),
  Lat = mean(Lat), Lon= mean(Lon), StandBA = mean(StandBA, na.rm=T), meanDBH = mean(meanDBH, na.rm=T),
  CLY=mean(CLY),SND=mean(SND), SLT = mean(SLT),
  PTO=mean(PTO), NTO=mean(NTO), AWC=mean(AWC), BDW=mean(BDW),
  ECE=mean(ECE), DES=mean(DES), DER=mean(DER), PC1fert = mean(PC1fert),
  PC2depth=mean(PC2depth), Elev=mean(Elev),
  avg_nleaves = mean(avg_nleaves, na.rm=T), max_area = max(max_area, na.rm=T),
  min_area=min(min_area, na.rm=T),
  sdWD = sd(WD, na.rm=T), medWD = median(WD, na.rm=T), WD = mean(WD, na.rm=T),
  sdLMA = sd(LMA, na.rm=T), medLMA = median(LMA, na.rm=T), LMA = mean(LMA, na.rm=T),
  sdLDMC = sd(LDMC, na.rm=T), medLDMC = median(LDMC, na.rm=T), LDMC = mean(LDMC, na.rm=T),
  sdAl_As = sd(Al_As, na.rm=T), medAl_As = median(Al_As, na.rm=T), Al_As=mean(Al_As, na.rm=T),
  sdhub = sd(hub, na.rm=T), medhub = median(hub, na.rm=T), hub = mean(hub, na.rm=T)
)

)
traits.s$AIc <- traits.s$PPTc/traits.s$PETc

# adding log.Al_As and huber value
traits.b$log.Al_As <- log(traits.b$Al_As, base=10)
traits.t$log.Al_As <- log(traits.t$Al_As, base=10)
traits.p$log.Al_As <- log(traits.p$Al_As, base=10)
traits.s$log.Al_As <- log(traits.s$Al_As, base=10)


traits.b$log.hub <- log(traits.b$hub, base=10)
traits.t$log.hub <- log(traits.t$hub, base=10)
traits.p$log.hub <- log(traits.p$hub, base=10)
traits.s$log.hub <- log(traits.s$hub, base=10)



# climate quantiles for the whole species distributions (for code used to create this, email LDL Anderegg)
quants <- read.csv("data/Climate_Quantiles_allspp_20200701.csv")
quants$Species <- quants$spp
levels(quants$Species) <- list(ACAC="acac",ESAL="esal",EMARG="emarg",COCA="coca",OVAT="ovat",VIMI="vimi",AMYG="amyg",OBLI="obli")


########### Species Means ################################
traits.sp <- traits.b %>% group_by(Species) %>% summarise(WD = mean(WD, na.rm=T), LMA=mean(LMA, na.rm=T), LDMC=mean(LDMC, na.rm=T), hub = mean(hub, na.rm=T))
traits.sp$log.hub <- log(traits.sp$hub, base=10)

traits.sp$MD.50 <- quants$quantCMDcALA.50[match(traits.sp$Species, quants$Species)]
traits.sp$MD.90 <- quants$quantCMDcALA.90[match(traits.sp$Species, quants$Species)]


### load summary of intraspecific residual patterns ########
# respats <- read.csv("Intraspecific_Results/Resid_Patterns_20180811.csv")
#respats <- read.csv("Intraspecific_Results/Resid_Patterns_20200116.csv") # with hub and updated trait-climate relationships with soil and stand
#respats <- read.csv("Intraspecific_Results/Resid_Patterns_20200209.csv") # updated when cleaning code. droped Al_As columns
respats <- read.csv("data/Resid_Patterns_20200526.csv") # updated with NP R2R
# ignore the Warning message. not sure why it throws it.

################# END: LOAD DATA #######################################


# AND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

















##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: VISUALIZE TRAIT-CLIMATE #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## . FIG 1: Trait-Env relationships W.HUBER & MD Chelsa ####################

# set a semi-transparent palette to deal with overplotting problems
palette(pallight)

# make a quartz device of the 1/2 column size (based on journal specs)
quartz(width=3.4, height=6)
#jpeg(file=paste0("./",results_dirname,"/Fig1_TraitClimate_ChelsaMD.jpg"), width=3.4, height=6, units = "in",res = 600)

# set the:
# mfrow= number of panels (mfrow=)
# mar= margin for each panel (mar=)
# oma= outer margin around all panels 
# mgp = location of margin labels and ticks
par(mfrow=c(4,1), mar=c(0,4,0,1), mgp=c(2.2,1,0), oma=c(3.5,0,3.5,0), cex.lab=1.3)

# loop over three traits (WD, LMA, LDMC) that can be handled the same way
for( pp in c(1:3)){
  k <- c("WD","LMA","LDMC")[pp]
  # plot the trait against climate, different color for each species
    # suppress x axis with xaxt="n"
  plot(get(k)~MDc, traits.t, pch=pchs, cex=.9, col=Species, ylab="", xaxt="n")
  
  
  # add in a legend above the very top plot (has to be called after plot created)
  if(k=="WD") {
    legend(x = -1200, y=1, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
           , text.font=3, pch=c(3,rep(16, times=7)), lty=1,lwd=1.5, col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=4, bty="n")
    # add in y axis legend for WD
    mtext(expression(paste("WD ",(g/cm^3))), side=2, line=2.2, font=2)
  }
    # add in y axis legend for LMA
  if(k == "LMA"){
    mtext(expression(paste("LMA ",(g/cm^2))), side=2, line=2.2, font=2)
  }
    # add in y axis legend for LDMC
  if(k == "LDMC"){
    mtext(expression(paste("LDMC ",(g/g))), side=2, line=2.2, font=2)
  }
  
  # calcualte and then overplot each species' linear fit, but only for the extent of the data
  for (j in 1:length(levels(traits.p$Species))){
    i <- levels(traits.p$Species)[j] 
    tmp <- traits.t[which(traits.t$Species==i),]
    tmpmod <- lm(get(k)~MDc, tmp)
    tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
    lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pal[j], lwd=3)
    #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
  }
}

# add in log.huber panel
plot(log.hub~MDc, traits.t, pch=pchs, cex=.9, ylab="", col=Species, xaxt="n", yaxt="n")
k <- "log.hub"
# points(get(k)~MD.50, traits.sp, pch=16, cex=2)
for (j in 1:length(levels(traits.p$Species))){
  i <- levels(traits.p$Species)[j] 
  tmp <- traits.t[which(traits.t$Species==i),]
  tmpmod <- lm(get(k)~MDc, tmp)
  tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
  lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pal[j], lwd=3)
  #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
}
# create a special axis for the logged trait
axis(2, labels = c("0.01","0.05","0.10","0.15"), at=log(c(0.01,0.05,0.10,0.15), base=10))
# add a climate axis to the lowest panel
axis(1)
# add an x-axis label using mtext(side=2)
mtext(text = "MD (PET-PPT, mm)", side=1, line=2.5)
# add a y axis label
mtext(expression(paste("HV ",(mm^2/cm^2))), side=2, line=2.2, font=2)
# if saving as a pdf, run this:
#quartz.save(file=paste0("./",results_dirname,"/Fig1_TraitClimate_ChelsaMD.pdf"),type = "pdf" )
# if saving as a jpeg with jpeg() function, run this to save figure:
#dev.off()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## . FIG S4: Trait-Env Across Species####################

# I want that same figure above, but I want to make the individual level data
# even more transparent and plot species averages on top to show
# within-species and among-species co-gradient patterns.

pallighter <- paste0(pal, "22")
palette(pallighter)
quartz(width=3.4, height=6)
par(mfrow=c(4,1), mar=c(0,4,0,1), mgp=c(2.2,1,0), oma=c(3.5,0,3.5,0), cex.lab=1.3)
for( pp in c(1:3)){
  k <- c("WD","LMA","LDMC")[pp]
  plot(get(k)~MDc, traits.t, pch=pchs, cex=.9, col=Species, ylab="", xaxt="n")
  
  if(k=="WD") {
    legend(x = -1200, y=1, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
           , text.font=3, pch=c(3,rep(16, times=7)), lty=1,lwd=1.5, col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=4, bty="n")
    mtext(expression(paste("WD ",(g/cm^3))), side=2, line=2.2, font=2)
    mod <- lm(WD~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
    abline(mod, lwd=2, lty=2, col="grey")
    mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)
  }
  if(k == "LMA"){
    mtext(expression(paste("LMA ",(g/cm^2))), side=2, line=2.2, font=2)
    mod <- lm(LMA~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
    abline(mod, lwd=2)
    mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)
    
    
  }
  if(k == "LDMC"){
    mtext(expression(paste("LDMC ",(g/g))), side=2, line=2.2, font=2)
    mod <- lm(LDMC~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
    abline(mod, lwd=2)
    mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)
    
    
  }
  
  for (j in 1:length(levels(traits.p$Species))){
    i <- levels(traits.p$Species)[j] 
    tmp <- traits.t[which(traits.t$Species==i),]
    tmpmod <- lm(get(k)~MDc, tmp)
    tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
    lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pallight[j], lwd=3)
    points(get(k)~MD.50, traits.sp[which(traits.sp$Species==i),], pch=ifelse(i=="ACAC",3,16), cex=2.5, col=pal[j])
    #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
  }
}

# add in log.huber panel
plot(log.hub~MDc, traits.t, pch=pchs, cex=.9, ylab="", col=Species, xaxt="n", yaxt="n")
k <- "log.hub"
# points(get(k)~MD.50, traits.sp, pch=16, cex=2)
for (j in 1:length(levels(traits.p$Species))){
  i <- levels(traits.p$Species)[j] 
  tmp <- traits.t[which(traits.t$Species==i),]
  tmpmod <- lm(get(k)~MDc, tmp)
  tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
  lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pallight[j], lwd=3)
  points(get(k)~MD.50, traits.sp[which(traits.sp$Species==i),], pch=ifelse(i=="ACAC",3,16), cex=2.5, col=pal[j])
  #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
}
mod <- lm(log.hub~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
abline(mod, lwd=2)
mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)

axis(2, labels = c("0.01","0.05","0.10","0.15"), at=log(c(0.01,0.05,0.10,0.15), base=10))
axis(1)
mtext(text = "MD (PET-PPT, mm)", side=1, line=2.5)
mtext(expression(paste("HV ",(mm^2/cm^2))), side=2, line=2.2, font=2)
#quartz.save(file=paste0("./",results_dirname,"/FigS4_TraitClimate_ChelsaMD_SpeciesMean.pdf"),type = "pdf" )


### significance of among-species trait relationships
summary(lm(WD~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.05235 .
summary(lm(LMA~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.0309 *
summary(lm(LDMC~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.008513 **
summary(lm(log.hub~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.003756 **
















### [insert 1800 lines of inelegant stats code here about trait-climate relationships]




##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: TRAIT COORDINATION #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### [insert 1300 lines of equally inelegant stats code here about trait coordination]

### read in the only relevant output from the above stats:
  #distribution of correlations coefficients between WD and LMA at different scales
WD.LMA.spp.results <- read.csv("data/WD_LMA_spp.csv")
WD.LMA.branch.results <- read.csv("data/WD_LMA_branch.csv")
WD.LMA.tree.results <- read.csv("data/WD_LMA_tree.csv")

  #distribution of correlations coefficients between LDMC and LMA at different scales
LDMC.LMA.spp.results <- read.csv("data/LDMC_LMA_spp.csv")
LDMC.LMA.branch.results <- read.csv("data/LDMC_LMA_branch.csv")
LDMC.LMA.tree.results <- read.csv("data/LDMC_LMA_tree.csv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######### . FIG 2new: Trait covariation 6 panels (WD, LMA, LDMC, log.hub) ###################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# where to put the panel labels 
letlin <- 0 # move panel labels up and down
letadj <- 0 # move panel labels left and right

palette(pallight)

# Make a big plotting device
quartz(width=7.4, height=4.1)
#jpeg(file=paste0("./",results_dirname,"/Fig2_TraitCorrelations.jpg"), width=7.4, height=4.1, res=600, units="in")

# split the device up into 7 plotting panes
mymat <- matrix(c(1,2,3,7
                  ,4,5,6,7), nrow=2, byrow=T)
layout(mat = mymat,widths = c(.9,.9,.9,1.3))
# visualize how your layout works
layout.show(n = 7)

# set the par for the first six plots
par( mar=c(3.6,3.7,1,.5), mgp=c(2,.7,0), cex.lab=1.4, oma=c(0,.2,2,0))

# plot trait-trait relationships for all teh species
plot(WD~LDMC, traits.t, pch=pchs, cex=.9, col=Species, ylab = expression(paste("WD ",(g/cm^3))), xlab="LDMC (g/g)")
n_sig <- 0 # reset the counter of number of species with significant trait-trait correlations
cors <- rep(NA, times=8) # preallocate a vector of correlations coeffs
  # loop over species and calcualte the trait-trait correlations
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LDMC"], traits.t[which(traits.t$Species==i),"WD"])
  cors[j]<- correlation$estimate # store the correlations coefficient
    # if significant:
  if(correlation$p.value <= 0.05){
    # plot MAR with a solid line
    plot.MAR(xvar = "LDMC", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    # and add one to the significance counter
    n_sig <- n_sig + 1
  }
    # if not significant, plot MAR with dotted line
  else{
    plot.MAR(xvar = "LDMC", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("a)",side=3, line=letlin, adj=letadj)

# add text in corner of plot showing how many of 8 species have significant relationships
mtext(text = paste(n_sig,"8", sep="/"), side = 3, line= -1.3, adj=.1)
# add smaller text showing range of correlation coefficients across species
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=3, line=-2,adj=.05, cex=.7)
# legend(x = 0.1, y=.97, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
#        , text.font=3, pch=c(3, rep(16, times=7)), col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=8, bty="n")

# do it all again for different traits.
plot(WD~LMA, traits.t, pch=pchs, cex=.9, col=Species,ylab=expression(paste("WD ",(g/cm^3))), xlab=expression(paste("LMA ",(g/cm^2))))
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"WD"], traits.t[which(traits.t$Species==i),"LMA"])
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(xvar = "LMA", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(xvar = "LMA", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("b)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 3, line= -1.3, adj=.1)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=3, line=-2,adj=.05, cex=.7)


plot(LMA~LDMC, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste("LMA ",(g/cm^2))), xlab="LDMC (g/g)")
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LDMC"], traits.t[which(traits.t$Species==i),"LMA"])
  #print(paste(i, correlation$p.value))
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(xvar = "LDMC", yvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(xvar = "LDMC", yvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("c)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 3, line= -1.3, adj=.1)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=3, line=-2,adj=.05, cex=.7)




plot(log.hub~LDMC, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste(log[10](HV))), xlab="LDMC (g/g)")
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LDMC"], traits.t[which(traits.t$Species==i),"log.hub"])
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(yvar = "log.hub", xvar = "LDMC",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(yvar = "log.hub", xvar = "LDMC",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("d)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 1, line= -2, adj=.9)# -1.3
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=1, line=-1.1,adj=.95, cex=.7)



plot(log.hub~LMA, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste(log[10](HV))), xlab=expression(paste("LMA ", (g/cm^2))))
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LMA"], traits.t[which(traits.t$Species==i),"log.hub"])
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(yvar = "log.hub", xvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(yvar = "log.hub", xvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("e)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 1, line= -2, adj=.9)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=1, line=-1.1,adj=.95, cex=.7)



plot(log.hub~WD, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste(log[10](HV))), xlab=expression(paste("WD ", (g/cm^3))))
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"log.hub"], traits.t[which(traits.t$Species==i),"WD"])
  print(paste(i, correlation$p.value))
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(yvar = "log.hub", xvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(yvar = "log.hub", xvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("f)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 1, line= -2, adj=.9)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=1, line=-1.1,adj=.95, cex=.7)




## Final panel density plot

# set a different par to leave lots of room on top for legends
par(mar=c(5,3.2,8,1),mgp=c(2,1,0))

# make density plots
plot(density(WD.LMA.spp.results$Rho, na.rm=T, adjust=2),lwd=3, xlim=c(-1,1), main="",xlab='Correlation')
lines(density(WD.LMA.branch.results$Rho, na.rm=T, adjust=2),lwd=1, lty=5)
lines(density(WD.LMA.tree.results$Rho, na.rm=T, adjust=2),lwd=2, lty=3)

lines(density(LDMC.LMA.spp.results$Rho, na.rm=T, adjust=2),lwd=3,col="blue")
lines(density(LDMC.LMA.branch.results$Rho, na.rm=T, adjust=2),lwd=1, lty=5,col="blue")
lines(density(LDMC.LMA.tree.results$Rho, na.rm=T, adjust=2),lwd=2, lty=3,col="blue")

# add a legend for the colors (different traits)
legend('topleft', legend = c("LMA-WD", "LMA-LDMC"), lwd=3, col=c("black","blue"), lty=1, bty="n")
# add a legend for the lines (different scales)
legend('top',inset = -.25,xpd=NA, horiz = F, legend=c("branch","individ","site"), lty=c(5,3,1), lwd=c(1,2,3), bty="n")
mtext("g)",side=3, line=letlin, adj=letadj)
# add a legend for the other plots that didn't have enough room
legend(x=-0.25, y=3.5,xjust=.5, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
       , text.font=3, pch=c(3, rep(16, times=7)), col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=3)
#dev.off()
#quartz.save(file=paste0("./",results_dirname,"/Fig2_TraitCorrelations.pdf"),type = "pdf" )



# END FIG 2new ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


################# END: TRAIT COORDINATION ####################################







##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: VARIANCE DECOMPOSITION #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++




#  Approach: use mixed effects models for each species and trait with only random effects (or only branch diameter as fixed effect for WD). Then put the variance components all in one dataframe for additional analysis and visualization

#Note 1: theoretically it shouldn't matter whether these are nested random effects or not, as long as all the lower-level effects are globally unique IDs (according to Zuur et al. 2009).
#In practice, this seems to be qualitatively true but not precisely quantitatively true. It doesn't change any conclusions to define them as nested, but convergence becomes a MAJOR issue for some traits/species. So here I have opted to stick with independant random effects for each organizational level for computation's sake.

#Note 2: after an update of the {lme4} function, some of these models started throwing convergence errors. These models still appear to be fitting appropriately, and have not been tweaked to remove the error. Also 'boundary (singular) fit' warnings indicate that a variance parameter was estimated to be 0.



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########### . Variance Decomposition for all species, all traits #######
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++ first for WD ++++++++++
testACACvd <- lmer(WD~ Diameter + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(WD~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(WD~ Diameter + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(WD~1 + Diameter + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(WD ~ Diameter + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
# testAMYGvd <- lmer(WD ~ Diameter + (1|Plottag/Sitetag/Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(WD ~ Diameter + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(WD ~ Diameter + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(WD ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesWD <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesWD) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesWD <- variancesWD[c(3,2,1,4),]
climsens <- c("PET", "MD", "ppt", "pet","PPT+PET", "PPT", "PPT", "pet")







#++++ Next for LMA ++++++++++
# note: raw LMAs and logged SLAs look pretty similar in a model with all spp. So I think I'll go raw LMA for the time being...
# The general inferences from raw LMA and logged SLA are pretty similar,
# with the rank orders being unchanged except for AMYG and OBLI.
testACACvd <- lmer(LMA~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(LMA~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(LMA~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(LMA~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(LMA ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(LMA ~ 1 +  (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(LMA ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(LMA ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesLMA <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesLMA) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesLMA <- variancesLMA[c(3,2,1,4),]







#++++ Then LDMC ++++++++++
testACACvd <- lmer(LDMC~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(LDMC~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(LDMC~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(LDMC~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesLDMC <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesLDMC) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesLDMC <- variancesLDMC[c(3,2,1,4),]




#++++ Next for Al_As ++++++++++
testACACvd <- lmer(log.Al_As~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(log.Al_As~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(log.Al_As~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(log.Al_As~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesAl_As <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesAl_As) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesAl_As <- variancesAl_As[c(3,2,1,4),]





#++++ Next for hub ++++++++++
# note: actually mathmatically identical to variance decomp of log.Al_As because log.Al_As = -1* log.hub
testACACvd <- lmer(log.hub~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(log.hub~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(log.hub~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(log.hub~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
varianceshub <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(varianceshub) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
varianceshub <- varianceshub[c(3,2,1,4),]




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######### . Full variance decomp of all species ##############

fullWD <- lmer(WD ~ Diameter + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
fullLDMC <- lmer(LDMC ~ 1 + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
fullLMA <- lmer(LMA ~ 1 + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
fullHV <- lmer(log.hub ~ 1 + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
WDvar <- data.frame(VarCorr(fullWD))
LDMCvar <- data.frame(VarCorr(fullLDMC))
LMAvar <- data.frame(VarCorr(fullLMA))
HVvar <- data.frame(VarCorr(fullHV))

variancesall <- data.frame(WDvar[,4],LMAvar[,4], LDMCvar[,4],HVvar[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","bSpecies","wiTree"))
colnames(variancesall) <- c("WD","LMA","LDMC","HV")
# and reorder variance componentsto make them nice plots
variancesall <- variancesall[c(4,3,2,1,5),]
variancesall.scaled <- apply(variancesall, MARGIN=2, FUN = function(x){x/sum(x)})

##### . scaling all variances #############
scaledvariancesWD <- apply(variancesWD, MARGIN=2, FUN= function(x){x/sum(x)})
scaledvariancesLMA <- apply(variancesLMA, MARGIN=2, FUN= function(x){x/sum(x)})
scaledvariancesLDMC <- apply(variancesLDMC, MARGIN=2, FUN= function(x){x/sum(x)})
scaledvarianceshub <- apply(varianceshub, MARGIN=2, FUN= function(x){x/sum(x)})
totalvariancesWD <- colSums(variancesWD)
totalvariancesLMA <- colSums(variancesLMA)
totalvariancesLDMC <- colSums(variancesLDMC)
totalvarianceshub <- colSums(varianceshub)


### .  combining variance estimates into one dataframe:  #################
# # first make them a wide-form dataframe

varswide <- data.frame(rbind(rbind(variancesWD, colSums(variancesWD)), rbind(variancesLMA, colSums(variancesLMA)),
                             rbind(variancesLDMC, colSums(variancesLDMC)), rbind(variancesAl_As, colSums(variancesAl_As)), rbind(varianceshub, colSums(varianceshub))),"trait"=rep(c("WD","LMA","LDMC","Al_As","hub"), each=5))
varswide$level <- rep(c("bSite","wiSite","wiPlot","wiTree","Total"), times=5)

# then melt them into long form
varslong <-melt(data = varswide, id.vars = c("trait","level") , variable.name="Species")
colnames(varslong)[which(colnames(varslong)=="value")] <- "variance"
# now add the species aridity niche
quants$Species <- toupper(quants$spp) # make species codes upper case to match with varslong

# add in species climate niche info from distribution records
# driest range edge
varslong$CMD.90 <- quants$quantCMDcALA.90[match(varslong$Species,quants$Species)] 
varslong$PPT.10 <- quants$quantPPTcALA.10[match(varslong$Species,quants$Species)]
varslong$PET.90 <- quants$quantCMDcALA.90[match(varslong$Species,quants$Species)]
# median
varslong$CMD.50 <- quants$quantCMDcALA.50[match(varslong$Species,quants$Species)]
varslong$PET.50 <- quants$quantPETcALA.50[match(varslong$Species,quants$Species)]
varslong$PPT.50 <- quants$quantPPTcALA.50[match(varslong$Species,quants$Species)]
varslong$level <- factor(varslong$level)
varslong$Species.trait <- with(varslong, paste(Species, trait, sep="."))





## Some might argue that coefficients of variation are more appropriate measures than raw variances for analyzing the change in variance components with aridity. So this code calculates CVs for plotting SI figures.



trait.means <- traits.b %>% group_by(Species) %>% summarise(WD = mean(WD, na.rm=T), LDMC=mean(LDMC, na.rm=T), LMA=mean(LMA, na.rm=T), Al_As=log(mean(Al_As, na.rm=T), base=10), hub=-1*log(mean(hub, na.rm=T), base=10)) #%>% mutate(log.Al_As = log(Al_As,base=10)) %>% select(-Al_As)
# NOTE: Al_As is actually log10(Al_As) but had to keem name for matching with varslong
trait.means.long <- trait.means %>% gather(trait, mean, -Species)
trait.means.long$Species.trait <- paste(trait.means.long$Species, trait.means.long$trait, sep=".")





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## . FIG4: Barplot of variance Decomposition ###########
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## set color choices
#pal.vardecomp <- brewer.pal(n=9, "Set1")
pal.vardecomp <- c(paste0(brewer.pal(n=3, "Set1")[1],"66"), brewer.pal(n=9, "Blues")[c(8,5,2)])
palette(pal.vardecomp)
# set the colors that wil denote within-species, within-genus, within family and across CWMs
colchoices <- c(1,2,4,3,6)


#quartz(width=4.3,height=6.4) 
jpeg(file=paste0("./",results_dirname,"/Fig4_VarainceDecomp_scaled.jpg"), width=4.3, height=6.4, units="in", res=600)
par(mfrow=c(4,1), mar=c(1,3.6,1,3.6), mgp=c(2.3,1,0), oma=c(3.6,0,3,0), cex.lab=1.4, cex.axis=1.1)
p<-barplot(height=as.matrix(scaledvariancesWD)
           , beside=F, names.arg=rep(NA, times=8)
           , col = pal.vardecomp
           , legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree")
           , args.legend=list(bty="n", x=mean(c(1.1,16)), y=1.5, xpd=NA, cex=1.3,xjust=0.5, ncol=4)
           , ylab="% WD Var", las=3
           , xlim=c(1.1,16.2), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
axis(4,at=c(0,0.002,.004,0.006,.008)/max(totalvariancesWD), labels = c("0","0.002","","0.006",""), xpd=NA,col="#333333", col.axis="#333333")
mtext("Tot WD Var", side = 4, line =2.3)
barplot(height=as.matrix(t(totalvariancesWD/max(totalvariancesWD))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
#text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
mtext(text = "a)", side=3, adj=-0.12, line=.2)

p<-barplot(height=as.matrix(scaledvariancesLMA)
           , beside=F, names.arg=rep(NA, times=8)
           , col = pal.vardecomp
           #, legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree"), args.legend=list(bty="n")
           , ylab="% LMA Var", las=3
           , xlim=c(1.1,16.2), width=1, space=1)#, yaxt="n")# , ylim=c(0,2.5e-5)
barplot(height=as.matrix(t(totalvariancesLMA/max(totalvariancesLMA))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
mtext(text = "b)", side=3, adj=-0.12, line=.2)
mtext("Tot LMA Var", side = 4, line =2.3)
axis(4, at=c(0,.5e-5,1e-5,1.5e-5,2e-5)/max(totalvariancesLMA), labels=c("0","","1e-5","","2e-5"),col="#333333", col.axis="#333333")
p<-barplot(height=as.matrix(scaledvariancesLDMC)
           , beside=F, names.arg=rep(NA, times=8)
           , col = pal.vardecomp
           #, legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree"), args.legend=list(bty="n")
           , ylab=" % LDMC Var", las=3
           , xlim=c(1.1,16.2), width=1, space=1)#, yaxt="n")
barplot(height=as.matrix(t(totalvariancesLDMC/max(totalvariancesLDMC))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
mtext(text = "c)", side=3, adj=-0.12, line=.2)
mtext("Tot LDMC Var", side = 4, line =2.3)
axis(4,at=c(0,0.001,.002,0.003,.004)/max(totalvariancesLDMC), labels = c("0","","2e-3","","4e-3"),col="#333333", col.axis="#333333" )
# p<-barplot(height=as.matrix(scaledvariancesAl_As)
#            , beside=F, names.arg=rep(NA, times=8)
#            , col = pal.vardecomp
#            #, legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree"), args.legend=list(bty="n")
#            , ylab=expression(paste("Var in ",log[10](A[L]:A[S]))), las=3)
p <- barplot(height=as.matrix(scaledvarianceshub)
             , beside=F, names.arg=rep(NA, times=8)
             , col=pal.vardecomp
             , ylab=expression(paste("% ",log[10](HV)," Var"))
             , xlim=c(1.1,16.2), width=1, space=1)
barplot(height=as.matrix(t(totalvarianceshub/max(totalvarianceshub))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
mtext(expression(paste("Tot ",log[10](HV)," Var")), 4, line=2.3)
# note: the var decomp of log.Al_As and log.hub are identical, because they are perfectly negatively correlated
mtext(text = "d)", side=3, adj=-0.12, line=.2)
axis(4,at=seq(from=0,to=0.06, by=0.02)/max(totalvarianceshub), labels=c("0","2e-3",NA,"4e-3"),col="#333333", col.axis="#333333")

axis(1,at=p,labels= c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl"),font=3, cex.axis = 1.4, tick=F,las=2)
dev.off()
#quartz.save(file=paste0("./",results_dirname,"/Fig4_VarainceDecomp_scaled.pdf"),type = "pdf" )



############### FIGURE S6: Scaled barplots for full dataset #################3
quartz(width=5.3,height=3.5) 
par( mar=c(2,3.6,1,6), mgp=c(2.3,1,0), oma=c(3.6,0,0,0), cex.lab=1.4, cex.axis=1.1)
p<-barplot(height=as.matrix(variancesall.scaled)
           , beside=F, names.arg=c("WD","LMA","LDMC","log10(HV)")
           , col = c(pal[1],pal.vardecomp), 
           , legend.text= c("Btw Species","Btw Sites", "Within Site", "Within Plot", "Within Tree")
           , args.legend=list(bty="n", x=5.8, y=1, xpd=NA, cex=1,xjust=0.5, ncol=1)
           , ylab="% Total Variance", las=3
)#,  yaxt="n" # ylim=c(0,0.008),
quartz.save(file=paste0("./",results_dirname,"/FigS7_VarainceDecomp_scaled_fulldataset.pdf"),type = "pdf" )




