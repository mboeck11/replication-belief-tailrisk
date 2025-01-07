#------------------------------------------------------------------------#
# Belief Shocks Arising from Predictive Distributions                    #
#                                                                        #
# Creating Ensemble Forecast of SPF                                      #
#                                                                        #
# Michael Pfarrhofer and Maximilian Boeck                                #
#------------------------------------------------------------------------#
rm(list=ls())

# load packages
library(scoringRules)
library(readxl)
library(zoo)
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)

# setup
set.p <- c(0.05,0.1,0.16,0.5,0.84,0.9,0.95)#seq(0.05,0.95,by=0.01) # quantiles we consider

# read data
spf <- read_xlsx(path = "./data/meanGrowth.xlsx", sheet="RGDP")
spf$drgdp6[spf$drgdp6 == "#N/A"] = NA_real_
spf$drgdp6 = as.numeric(spf$drgdp6)
spf <- ts(spf[,3:7], start=c(1968,4), frequency=4)
colnames(spf) <- seq(0,4)

spfmicro = read_xlsx(path = "./data/SPFmicrodata.xlsx", sheet="RGDP")
bigT     = nrow(spf)
bigN     = ceiling(nrow(spfmicro)/bigT)

indiv    = list()
years    = unique(spfmicro$YEAR)
quarters = sort(unique(spfmicro$QUARTER))
count    = 1

yearlabs <- NULL
for(tt in 1:length(years)){
  curr_year = years[tt]
  for(qq in 1:length(quarters)){
    if(tt == 1 & qq < 4) next
    if(tt == length(years) & qq > 3) next 
    
    curr_quart = quarters[qq]
    
    idx <- which(spfmicro$YEAR == curr_year & spfmicro$QUARTER == curr_quart)
    temp = spfmicro[idx,paste0("RGDP",seq(1,6))]
    temp[temp=="#N/A"] = NA
    temp = apply(temp,2,as.numeric)
    
    # compute annualized growth rate
    temp <- (((temp[,2:ncol(temp)]/temp[,1])^4)-1)*100
    indiv[[count]] = temp
    
    count <- count + 1
    yearlabs <- c(yearlabs,paste0(curr_year,"Q",qq))
  }
}
names(indiv) = yearlabs
indiv        = indiv[1:212] # cut sample to 1968Q4 to 2021Q3

# delete not necessary objects
rm(spfmicro, temp)

# compute stats
means       = ts(t(sapply(indiv, function(l) apply(l, 2, mean, na.rm=TRUE))),start=as.numeric(unlist(strsplit(names(indiv)[1],"Q"))),frequency=4)
vars        = ts(t(sapply(indiv, function(l) apply(l, 2, sd, na.rm=TRUE))),start=as.numeric(unlist(strsplit(names(indiv)[1],"Q"))),frequency=4)
quants      = t(lapply(indiv, function(l) apply(l, c(2), quantile, probs=set.p, na.rm=TRUE)))
quants.time = as.Date.ts(means)

quants                  = array(as.numeric(unlist(quants)), dim=c(7, 5, length(quants)))
dimnames(quants)        = list(paste0("p",100*set.p),colnames(means),as.character(quants.time))
nc_quant_mean           = ts(cbind(t(quants[,"RGDP2",]),means[,"RGDP2"]), start=c(1968,4), frequency=4)
colnames(nc_quant_mean) = str_extract(colnames(nc_quant_mean),"(?<=\\.).*")
colnames(nc_quant_mean)[ncol(nc_quant_mean)] = "mean"
quant_spf_df            = melt(quants)
quant_spf_df$Var3       = as.character(quant_spf_df$Var3)

save(nc_quant_mean, file ="./results/SPF_mean_quantiles.rda")

# load actual data
vintages = read_excel("./data/routput_first_second_third.xlsx", sheet = "DATA", col_types = c("text", 
                                                                                              "numeric", "numeric", "numeric", 
                                                                                              "numeric"), skip = 3)
vintages = window(ts(vintages[,-1],start=c(1965,3),frequency=4),start=c(1968,4),end=c(2021,3))
bigT     = NROW(vintages)
count    = 1
rw_hh    = 12

# optimize predictive metric for variance
pred_loss <- function(sd_est){
  crps_sum <- 0
  for(i in 1:rw_hh){
    tmp <- as.numeric(na.omit(indiv_sub[[i]][,1]))
    for(j in 1:length(tmp)){
      crps_sum <- crps_sum + crps_norm(vintage_sub[i],mean=tmp[j],sd=sd_est)
    }
  }
  return(crps_sum)
}

theta     = rep(NA,bigT - rw_hh)
theta_lab = NULL
for(hh in (rw_hh+1):length(indiv)){
  id_train  = count:(hh-1)
  indiv_sub = indiv[id_train]
  
  vintage_sub = vintages[id_train]
  theta_lab   = c(theta_lab,names(indiv)[max(id_train)+1])
  
  theta[count] = optim(1,pred_loss,method="Brent",lower=0,upper=100)$par
  count = count+1
}
names(theta) = theta_lab
theta        = ts(theta,start=as.numeric(unlist(strsplit(names(theta)[1],"Q"))),frequency=4)
Dt           = window(vars[,1],start=start(theta))

# ts.plot(ts.intersect(Dt,theta),lwd=2,col=c("black","blue"))
ytrue = window(vintages[,1],start=start(theta))
mu    = window(means[,1],start=start(theta))

true_df <- data.frame("date"=as.character(as.Date.ts(ytrue)),"GDP"=as.numeric(ytrue),
                      "ybar"=as.numeric(mu),"s2"=as.numeric(Dt^2),"sigma2"=theta^2)

# sample from the mixture of individual distributions
indiv_sub <- indiv[names(theta)]
nsave_i   <- 10000

mc_sample <- array(NA,dim=c(length(indiv_sub),nsave_i))
for(tt in 1:length(indiv_sub)){
  indiv_tt <- as.numeric(na.omit(indiv_sub[[tt]][,1]))
  theta_tt <- theta[tt]
  
  fc_id <- sample(1:length(indiv_tt),size=nsave_i,replace=TRUE)
  mc_sample[tt,] <- rnorm(nsave_i,indiv_tt[fc_id],theta_tt)
}

nc_quant <- ts(t(apply(mc_sample,c(1),quantile,probs=set.p)),start=start(theta),frequency=4)
nc_quant <- ts.intersect(nc_quant,ts(apply(mc_sample,c(1),mean),start=start(theta),frequency=4))
colnames(nc_quant) <- c(paste0("p",100*set.p),"mean")

save(nc_quant, file ="./results/SPF_quantiles.rda")
