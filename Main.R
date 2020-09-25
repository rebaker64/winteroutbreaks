
setwd("~/")

# load functions
source('Functions.R', encoding = 'UTF-8')
loadpackages()

#load R_e estimates for NY
load("NYestimR.RData")

#load case data for NY
load("NYcases.RData")

#load Specific Humidity data for NY
load("NYq.RData")

#vary range of starting number infected 
IList <- rev(seq(0.005,0.5,0.01))

#vary range of proportional reduction in R0 during NPI
RChangeVar <- seq(0.65,1,0.01)

#######set up outcomes vars

# proportion susceptible in July, an outcome of varying start I/N
SOut <- rep(NA, times = length(IList))

# size of wintertime peak in climate scenario
SizeClimWinterI <- matrix(NA, nrow = length(IList), ncol = length(RChangeVar))

# size of wintertime peak in constant scenario
SizeWinterI <- matrix(NA, nrow = length(IList), ncol = length(RChangeVar))

# timing of wintertime peak in climate scenario
TimingClimWinterI <- matrix(NA, nrow = length(IList), ncol = length(RChangeVar))

# timing of wintertime peak in constant scenario
TimingWinterI <- matrix(NA, nrow = length(IList), ncol = length(RChangeVar))

#setting plots = TRUE generates times series I/S and R0 plots i.e. Fig 1a/1b
## fix i and j to generate a specific time series
plotting = FALSE
for(j in 1:length(RChangeVar)){
  SOut <- rep(NA, times = length(IList))
for(i in 1:length(IList)){
    lead_time <- yday(NYestimR$date[1])
    end_data_time <- yday(max(NYestimR$date))
    RchangeUse <- RChangeVar[j]
    popuse <- 8000000
    ISetStart <- IList[i]
    
    #run the model with climate driven R0
    predNPI <- runModelWithR(R0min = 1.5, R0max = 2.5, Immunity = 66.25*7, LatCity = 40.7128, 
                             LonCity = -74.6, Var = -227.5, 
                             SHDAT = qout, pop = popuse, Lead = lead_time, R0list = NYestimR$`Median(R)`,
                             timeLengthSim = 2000 , SSet ="Orig", ISet = ISetStart, Rchange = RchangeUse)
    
    #normalize cases relative to observed I/N
    startcases <- yday(NYcases$time[1])
    casesI <- c(rep(0, times = (startcases - 1)), diff(NYcases$cases))
    ratio <- (mean(predNPI$I[1:length(casesI)])/mean(casesI))
    casesI <- casesI*(mean(predNPI$I[1:length(casesI)])/mean(casesI))
 
    # what is S at the end of the data series?
    minSClim <- min(predNPI$S[1:end_data_time])/popuse
    SOut[i] <- minSClim
    
    #get the winter peak size and timing of peak
    ts <- seq(1,10,1/364)[1:length(predNPI$I)]
    I <- predNPI$I
    WinterI <- I[ts > 1.75 & ts <2.25]
    maxWinterIClim <- max(WinterI)/popuse
    SizeClimWinterI[i,j] <-  maxWinterIClim
    I[1:end_data_time] <- 0
    timeMax <- ts[which.max(I[1:708])]
    TimingClimWinterI[i,j] <-  timeMax
    
    # for constant scenario, take the R0 in July, adjusted by NPIs, and project forward
    R0equivclim <- predNPI$R0[1 + lead_time+nrow(NYestimR)]
    R0listNew <- c(NYestimR$`Median(R)`, rep(R0equivclim,2000))[1:2000]
    predConstant <- runModelWithR(R0min = 1.5, R0max = 2.5, Immunity = 66.25*7, LatCity = 40.7128, 
                                  LonCity = -74.6, Var = -227.5, 
                                  SHDAT = qout, pop = popuse, Lead = lead_time, R0list = R0listNew,
                                  timeLengthSim = 2000 , SSet ="Orig", ISet = ISetStart)
    
    #get the winter peak size and timing of peak
    ts <- seq(1,10,1/364)[1:length(predNPI$I)]
    I <- predConstant$I
    WinterI <- I[ts > 1.75 & ts <2.25]
    maxWinterIConst <- max(WinterI)/popuse
    SizeWinterI[i,j] <-  maxWinterIConst
    I[1:end_data_time] <- 0
    timeMax <- ts[which.max(I[1:708])]
    TimingWinterI[i,j] <-  timeMax
    
    if(plotting == TRUE){
     par(mfrow=c(2,1))
     par(mar=c(3,3,1,3))
     pal <- hp(10,option = "ronweasley2") # highly recommend the harrypotter colour palettes :)
     plot(seq(2020,2024,1/364)[1:728],predNPI$R0[1:728],type="l", lwd = 2, bty = "n",  xlab = "",
          ylab = "",main = "New York", col="dodgerblue3")
     lines(seq(2020,2024,1/364)[1:728],predConstant$R0[1:728],type="l",col="black", lwd = 2, lty = 5)
     abline(h = 1, col="grey64")
     title(xlab = "Year", line = 2)
     title(ylab = "R0", line = 2)
     legend(2021.5, 6, lty=c(1,5), col=c("dodgerblue3","black"), legend=c("Climate","Constant") )
     
    
     plot(seq(2020,2024,1/364)[1:728],
           predNPI$I[1:728]/popuse,type="l",
           lwd = 2, bty = "n",  xlab = "",ylab = "",
           col="dodgerblue3")
     lines(seq(2020,2024,1/364)[1:728],predConstant$I[1:728]/popuse,type="l", lwd = 2, lty = 5, col="black")
     lines(seq(2020,2024,1/364)[1:length(casesI)],casesI/popuse,col="grey64")
     title(xlab = "Year", line = 2)
     title(ylab = "I/N", line = 2)
     par(new = TRUE)
     plot(seq(2020,2024,1/364)[1:728],predNPI$S[1:728]/popuse,type="l", lwd = 2, bty = "n", col=pal[10],xlab = "",ylab = "",
          axes = F)
     axis(side = 4)
     mtext(text = "S/N", line = 2, side = 4, col = pal[10])
     text(2021.5, 1, label = paste0("Cases RR = ",round(1/ratio,2)) )


}
    
  }
}


# Plot figures 1c-h
palette <- hp(25,option = "ronweasley2")
pdf(paste0("PropIClim_Base.pdf"),width=5,height=4)
filled.contour(SOut, seq(1,length(RChangeVar),by = 1),SizeClimWinterI - SizeWinterI,
               col = palette, key.title = "Proportion infected (Winter)", 
               plot.axes ={axis(2, at = rev(seq(0,length(RChangeVar),by = 5)), labels = seq(0,35,5)) ; axis(1)})
title(xlab="Proportion susceptible (July)", line = 2)
title(ylab="NPI reduction in R0 (%)", line = 3)
dev.off()

palette <- hp(18,option = "ronweasley2")
pdf(paste0("PropIClim.pdf"),width=5,height=4)
filled.contour(SOut,seq(1,length(RChangeVar),by = 1), SizeClimWinterI,
               col = palette, key.title = "Proportion infected (Winter)", zlim=c(0,0.0853), 
               plot.axes ={axis(2, at = rev(seq(0,length(RChangeVar),by = 5)), labels = seq(0,35,5)) ; axis(1)})
title(xlab="Proportion susceptible (July)", line = 2)
title(ylab="NPI reduction in R0 (%)", line = 3)
dev.off()

palette <- hp(18,option = "ronweasley2")
pdf(paste0("PropIBase.pdf"),width=5,height=4)
filled.contour(SOut, seq(1,length(RChangeVar),by = 1), SizeWinterI,
               col = palette, nlevel= length(palette), key.title = "Proportion infected (Winter)", zlim=c(0,0.0853), 
               plot.axes ={axis(2, at = rev(seq(0,length(RChangeVar),by = 5)), labels = seq(0,35,5)) ; axis(1)})
title(xlab="Proportion susceptible (July)", line = 2)
title(ylab="NPI reduction in R0 (%)", line = 3)
dev.off()

# set to NA if max is occurring at the limits of the time series - means peak is not occuring in the winter months
TimingClimWinterI[TimingClimWinterI>2.94230] <- NA
TimingWinterI[TimingWinterI>2.94230] <- NA
TimingClimWinterI[TimingClimWinterI<1.56] <- NA
TimingWinterI[TimingWinterI<1.56] <- NA

pdf(paste0("TimingDelta.pdf"),width=5,height=4)
palette <- hp(18,option = "lunalovegood") # another great palette from this package!
filled.contour(SOut, seq(1,length(RChangeVar),by = 1), TimingClimWinterI - TimingWinterI,
               col = palette, key.title = "Proportion infected (SH Summer)", 
               plot.axes ={axis(2, at = rev(seq(0,length(RChangeVar),by = 5)), labels = seq(0,35,5)) ; axis(1)})
title(xlab="Proportion susceptible (July)", line = 2)
title(ylab="NPI reduction in R0 (%)", line = 3)
dev.off()


pdf(paste0("TimingClim.pdf"),width=5,height=4)
palette <- hp(15,option = "lunalovegood")
TimingClimWinterI <- TimingClimWinterI - ts[end_data_time]
filled.contour(SOut,seq(1,length(RChangeVar),by = 1), TimingClimWinterI,
               col = palette, key.title = "Proportion infected (SH Summer)", zlim=c(0,1.5), 
               plot.axes ={axis(2, at = rev(seq(0,length(RChangeVar),by = 5)), labels = seq(0,35,5)) ; axis(1)})
title(xlab="Proportion susceptible (July)", line = 2)
title(ylab="NPI reduction in R0 (%)", line = 3)
dev.off()

pdf(paste0("TimingBase.pdf"),width=5,height=4)
palette <- hp(15,option = "lunalovegood")
TimingWinterI <- TimingWinterI - ts[end_data_time]
filled.contour(SOut,seq(1,length(RChangeVar),by = 1), TimingWinterI,
               col = palette, key.title = "Proportion infected (SH Summer)", zlim=c(0,1.5), 
               plot.axes ={axis(2, at = rev(seq(0,length(RChangeVar),by = 5)), labels = seq(0,35,5)) ; axis(1)})
title(xlab="Proportion susceptible (July)", line = 2)
title(ylab="NPI reduction in R0 (%)", line = 3)
dev.off()




