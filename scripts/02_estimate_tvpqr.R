#------------------------------------------------------------------------#
# Belief Shocks Arising from Predictive Distributions                    #
#                                                                        #
# Creating Quantile Nowcasts of RGDP (real-time factors and              #
#                                     rolling-window estimation)         #
#                                                                        #
# Michael Pfarrhofer and Maximilian Boeck                                #
#------------------------------------------------------------------------#
rm(list=ls())

# libraries and auxiliary functions
library(readxl)
library(readr)
library(GIGrvg)
library(lubridate)

# setup
set.p      = c(0.05,0.1,0.16,0.5,0.84,0.9,0.95)
facnum_qd  = 4 # as in companion paper
facnum_md  = 8 # less factors than in companion paper
time       = seq(1972.00,2021.5,by=0.25)
time_date  = as.character(zoo::as.Date.ts(ts(NA_real_, start=c(1972,1), end=c(2021,3), frequency=4)))
time_date1 = as.character(ymd(zoo::as.Date.ts(ts(NA_real_, start=c(1972,2), end=c(2021,4), frequency=4))) %m+% months(2))
Traw       = length(time)

# set option for information set to obtain GDP nowcast quantiles
# --------------------------------------------------------------------------------------
# load data
vintages <- read_excel("./data/routput_first_second_third.xlsx", sheet = "DATA", col_types = c("text", 
                                                                                               "numeric", "numeric", "numeric", 
                                                                                               "numeric"), skip = 3)
vintages <- window(ts(vintages[,-1],start=c(1965,3),frequency=4),end=c(2021,3))

# load realtime factors
load("./data/facs_realtime.rda")
rm(fred_facs_rt)

# extract factors (four factors as in companion FRED-QD paper)
fred_facs_rw = lapply(fred_facs_rt_rw, function(l){
  l           = l[,c(1,2,3,4)]
  l[,c(1,4)]  = l[,c(1,4)] * (-1)
  colnames(l) = paste0("fac",1:facnum_qd)
  return(l)
})

# auxiliary variables not in real-time dataset
NFCI_tmp <- ts(read_excel("./data/NFCI.xls", skip = 10)[,-1], start=c(1971,1),frequency=4)
NFCI     <- ts(0,start=c(1965,1),end=c(2021,3),frequency=4)
window(NFCI,start=c(1971,1)) = window(NFCI_tmp, end=c(2021,3))

# --------------------------------------------------------------------------------------
# TVP-QR
setwd("./functions")
source("!qrdhs.R")
setwd("../")

nburn   = 1000
nsave   = 1000
thinfac = 1

mode    = "qr"
tvp     = "dhs"
sv      = FALSE
fhorz   = 0
P       = 4

spec <- paste0(mode,"_",tvp,"lags",P)
grid.p <- c(0.05,0.1,0.16,0.5,0.84,0.9,0.95)

est_ls     = vector(mode="list", length=Traw)
if(!file.exists(paste0("./results/GDP_quantiles_rt_rw_",spec,"_tmp.rda"))){
  quant_mcmc = array(NA_real_, c(nsave, Traw, length(grid.p)), dimnames=list(NULL,time_date,paste0("p",100*grid.p)))
}else{
  load(paste0("./results/GDP_quantiles_rt_rw_",spec,"_tmp.rda"))
}
for(tt in 1:Traw){
  cat(paste0("Round: ", tt, "/",Traw,".\n"))
  
  Yraw <- ts.intersect(window(vintages[,1], end=time[tt]),
                       window(NFCI, end=time[tt]),
                       fred_facs_rw[[time_date1[tt]]])
  colnames(Yraw) = c("routput", "NFCI", paste0("fac",1:facnum_qd))
  
  date.labs   <- zoo::as.Date.ts(Yraw)
  date.labs.p <- date.labs[-c(1:P)]
  
  Ylags = embed(Yraw,P+1)
  Y     = Ylags[,1,drop=FALSE]
  X     = cbind(Ylags[,(ncol(Yraw)+1):ncol(Ylags)],1)
  rownames(Y) <- rownames(X) <- as.character(date.labs.p)
  
  est_ls_tt     = list()
  quant_mcmc_tt = array(NA_real_,dim=c(nsave,nrow(Y),length(grid.p)))
  
  count <- 0
  for(pp in grid.p){
    count <- count+1
    if(!is.na(quant_mcmc[1,tt,count])) next
    tmp <- est_ls[[paste0("p",pp)]] <- try(tvpqr(Y=Y,X=X,Xout=NULL,mode=mode,p=pp,cons.cons=FALSE,tvp=tvp,sv=sv,fhorz=fhorz,
                                                 nburn=nburn,nsave=nsave,thinfac=thinfac,quiet=TRUE),silent=TRUE)
    if(is(tmp,"try-error")) next
    quant_mcmc_tt[,,count] = tmp$fit[,,1]
    quant_mcmc[,tt,count]  = quant_mcmc_tt[,nrow(Y),count]
  }
}
save(quant_mcmc, file=paste0("./results/GDP_quantiles_rt_rw_",spec,"_tmp.rda"))

# rule out quantile crossing
quant_mcmc_sort = apply(quant_mcmc,c(1,2),sort)
quant_post      = t(apply(quant_mcmc_sort,c(1,3),median))

n_sim       = 10
quant_draws = array(NA_real_,dim=c(nsave,n_sim,Traw))

pb <- txtProgressBar(min = 0, max = nsave, style = 3)
for(irep in 1:nsave){
  for(tt in 1:Traw){
    ytmp                  = quant_mcmc_sort[,irep,tt]
    quant_draws[irep,,tt] = rnum_from_p(ytmp,n=n_sim)
  }
  setTxtProgressBar(pb, irep)
}

quant_mcmc_post <- t(apply(quant_draws,3,quantile,probs=set.p))

Y = as.matrix(window(vintages[,1], start=c(1972,1)))
rownames(Y) = time_date

par(mar=c(4,4,1,1),mfrow=c(1,2))
ts.plot(cbind(quant_post,Y),col=c(rep("black",length(set.p)),"red"),lwd=c(rep(1,length(set.p)),2),ylim=c(-40,40))
ts.plot(cbind(quant_mcmc_post,Y),col=c(rep("black",length(set.p)),"red"),lwd=c(rep(1,length(set.p)),2),ylim=c(-40,40))

# storage
save(est_ls,quant_post,quant_mcmc_post,Y,file=paste0("./results/GDP_quantiles_rt_rw_",spec,".rda"))

