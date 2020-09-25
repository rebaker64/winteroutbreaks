

## wrapper for running the SIRS model
runModelWithR  <- function(pop = 3000000,
                          R0list = R0list,
                          Immunity = 40, # weeks that correspond to duration of immunity
                          LatCity = 40.7128,
                          LonCity = -74.6, SSet = "Orig", Lead = 0, Var = -180, birthrate=0,
                          timeLengthSim = 364, SHDAT = SHDAT, ISet = "Orig", R0min = 1.5, R0max = 2.5, Rchange = 1){
  
  library(deSolve)
  library(doBy)
  library(viridis)
  library("raster")
  library("lubridate")
  library("ncdf4")
  library("psych")
  library("sp")
  library("maps")
  library("maptools")
  library("rgdal")
 
  qout <- SHDAT
  if(Lead==0){
    qout <- rep(qout, each = 7)}
  if(Lead > 0){
    qout <- rep(qout, each = 7)
    qoutNew <- c(qout[Lead:length(qout)], qout[1:(Lead - 1)])
    qout <- qoutNew
  }
  qout <- rep(qout,length=timeLengthSim)
  R0 = exp(Var*qout + log(R0max - R0min)) + R0min
  R0 = R0*Rchange
  lR0 = length(R0list)
  
  R0[1:lR0] <- R0list
  
  if(SSet == "Orig"){ S = pop - 1}
  if(SSet != "Orig"){
    S  = pop*SSet}
  
  if(ISet == "Orig"){ I = 1}
  if(ISet != "Orig"){
    I = ISet
    S = S - I
  }
  R = pop - S - I
  xstart = c(S = S, I = I, R = R)
  times = seq(1, timeLengthSim, by = 1)
  qList <- rep(qout, length = length(times))
  paras = list(D = 5, L = Immunity, R0list = R0,  var = Var, birthrate = birthrate)
  
  # with all parameters set, run model
  out = as.data.frame(ode(xstart, times, SIRS_InterventionWithR, paras))
  out$R0 <- R0
  if(Lead > 0){ dfadd <- data.frame(time = rep(0,times = Lead), S = rep(pop, times = Lead), I = rep(0, times = Lead), R =rep(0, times = Lead), R0 = rep(0, times = Lead) )
  out <- rbind(dfadd, out)   }
  
  return(out)
}


SIRS_InterventionWithR <- function(time, state ,theta) {
  #browser()
  ## Parameters:
  D <- theta[["D"]]
  L <- theta[["L"]]
  var <- theta[["var"]]
  birthrate <- theta[["birthrate"]]
  timestart <- theta[["timestart"]]
  timeend <- theta[["timeend"]]
  R0list <- theta[["R0list"]]
  
  ## States:
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  N <- S + I + R
  
  ## ODEs:
  R0 <- R0list[time]
  
  beta = R0/D
  dS <- birthrate*N + (R/L) -beta * S * I/N - S*birthrate
  dI <- beta * S * I/N - (I/D) - I*birthrate
  dR <- (I/D) - (R/L) - R*birthrate
  
  return(list(c(dS, dI, dR)))
}

## wrapper for running the SIRS model
runModelControl  <- function(pop = 3000000,
                           Immunity = 40, # weeks that correspond to duration of immunity
                           LatCity = 40.7128,
                           LonCity = -74.6, SSet = "Orig", Lead = 0, Var = -180, birthrate=0,
                           timeLengthSim = 364, SHDAT = SHDAT, ISet = "Orig", R0min = 1.5, R0max = 2.5, Rchange = 1){
  
  
  latlist <- unique(SHDAT$lat)
  lonlist <- unique(SHDAT$lon)
  latid <- latlist[which.min(abs(latlist - LatCity))]
  lonid <- lonlist[which.min(abs(lonlist - LonCity))]
  
  SHDATRow <- as.numeric(SHDAT[SHDAT$lat==latid & SHDAT$lon==lonid,])
  
  qout <- as.numeric(SHDATRow[3:54])
  if(Lead==0){
    qout <- rep(qout, each = 7)}
  if(Lead > 0){
    qout <- rep(qout, each = 7)
    qoutNew <- c(qout[Lead:length(qout)], qout[1:(Lead - 1)])
    qout <- qoutNew
  }
  qout <- rep(qout,length=timeLengthSim)
  R0 = exp(Var*qout + log(R0max - R0min)) + R0min
  R0 = R0*Rchange
 
  
  if(SSet == "Orig"){ S = pop - 1}
  if(SSet != "Orig"){
    S  = pop*SSet}
  
  if(ISet == "Orig"){ I = 1}
  if(ISet != "Orig"){
    I = ISet
    S = S - I
  }
  R = pop - S - I
  xstart = c(S = S, I = I, R = R)
  times = seq(1, timeLengthSim, by = 1)
  qList <- rep(qout, length = length(times))
  paras = list(D = 5, L = Immunity, R0list = R0,  var = Var, birthrate = birthrate)
  
  # with all parameters set, run model
  out = as.data.frame(ode(xstart, times, SIRS_InterventionWithR, paras))
  out$R0 <- R0
  if(Lead > 0){ dfadd <- data.frame(time = rep(0,times = Lead), S = rep(pop, times = Lead), I = rep(0, times = Lead), R =rep(0, times = Lead), R0 = rep(0, times = Lead) )
  out <- rbind(dfadd, out)   }
  
  return(out)
}


loadpackages <- function(){
  library("pals")
  library(geosphere)
  require("plotrix")
  require("raster")
  require("rgdal")
  require("dplyr")
  require("doBy")
  require("ncdf4")
  require("lubridate")
  require(EpiEstim)
  require(dplyr)
  require(ggplot2)
  require(RCurl)
  require(reshape2)
  require(purrr)
  require(lubridate)
  require("harrypotter")
  library(deSolve)
  library(doBy)
  library(viridis)
  library("lubridate")
  library("ncdf4")
  library("psych")
  library("sp")
  library("maps")
  library("maptools")
  
}
