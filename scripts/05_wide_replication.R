#------------------------------------------------------------------------#
# Wide Replication: Belief Shocks VAR                                    #
#                                                                        #
# Maximilian Boeck, Michael Pfarrhofer                                   #
#------------------------------------------------------------------------#
rm(list=ls())

# libraries
library(stringr)
library(readxl)
library(readr)
library(forecast)

# load VAR functions
source("./scripts/functions/bvar_aux.R")

# use outcomes with realtime and/or rolling-window factors
rt = TRUE
rw = TRUE
if(!rt & rw) stop("not possible")

# load ekm IRF
load("./results/irf_ekm.rda")

# load nowcast errors
load(paste0("./results/nc_errors_DF",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".rda"))

# load data
spf_median = read_xlsx("../02_data/medianGrowth.xlsx", sheet="RGDP")
rgdp       = read_xlsx("../02_data/routput_first_second_third.xlsx", sheet="DATA", skip=3, col_names=TRUE)
outbs      = read_xls("../02_data/FRED/OUTBS.xls", skip=10)
rgdp_first = ts(as.numeric(rgdp$First), start=c(1965,3), frequency=4)
rgdp_third = ts(as.numeric(rgdp$Third), start=c(1965,3), frequency=4)
spf_median = ts(spf_median$drgdp2, start=c(1968,4), frequency=4)
outbs      = ts(log(outbs[,2])*100, start=c(1968,4), frequency=4)

data <- ts.intersect(rgdp_first, rgdp_third, spf_median, outbs)
data <- cbind(data, ne_first=(data[,"rgdp_first"]-data[,"spf_median"])/4, ne_third=(data[,"rgdp_third"]-data[,"spf_median"])/4)
colnames(data) <- str_remove(colnames(data), "data\\.")

rm(spf_median, rgdp, outbs, rgdp_first, rgdp_third)

#-------------------------------------------------------------------------

# estimate the VAR
varNames      = c("ne_constructed", "outbs")
varNames.plot = c("Nowcast Error", "Output")
n             = length(varNames)
set.p         = c(0.05,0.1,0.16,0.5,0.84,0.9,0.95,"Actual")
nhor          = 21
plag          = 4
SV            = FALSE

irf_sign_all = array(NA_real_, dim=c(length(set.p),7,n,n,nhor))
for(pp in 1:length(set.p)){
  
  if(pp == length(set.p)){
    tmp_sub = subset(full_df, variable == set.p[pp])
  }else{
    tmp_sub = subset(full_df, variable == paste0("p=",set.p[pp]))
  }
  tmp_nce = ts(tmp_sub$NCE/4, start=c(1972,1), frequency=4)
  tmp_nce = window(tmp_nce, end=c(2019,4))
  nce_sd  = sd(tmp_nce)
  
  # build Yraw
  Yraw     = ts.union(tmp_nce, window(data[,"outbs"], start=c(1972,1), end=c(2019,4)))
  colnames(Yraw) = varNames
  
  args = list(draws = 10000, burnin = 5000, thin=2, cons = TRUE, trend = TRUE, qtrend = TRUE)
  
  # Bayesian VAR
  run <- bvar_wishart_mn(Yraw, plag = plag, args=args)
  
  #------ compute IRFs
  thindraws = run$args$thindraws
  irf_store = array(NA_real_, c(thindraws, n, n, nhor),
                    dimnames=list(NULL, varNames, varNames, seq(0,nhor-1)))
  for(irep in 1:thindraws){
    
    if(irep%%50==0) print(paste0("Round: ",irep))
    
    temp    <- gen_compMat(A=run$A[irep,,], n=n, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    if(SV){
      SIGMA <- apply(run$SIGMA[irep,,,], c(2,3), median)
    }else{
      SIGMA   <- run$SIGMA[irep,,]
    }
    
    if(max(abs(Re(eigen(compMat)$values))) > 1.00)
      next
    
    shock <- t(chol(SIGMA))
    
    # sign restrictions
    condall = 0
    icounter <- 1
    MaxTries <- 10000
    while(condall == 0 && icounter < MaxTries){
      randMat <- matrix(rnorm(n^2),n,n)
      QR <- qr(randMat)
      Q <- qr.Q(QR)
      dimnames(Q) <- list(varNames, varNames)
      Q_bar <- Q%*%diag(((diag(Q)>0)-(diag(Q)<0)))
      dimnames(Q_bar) <- list(varNames, varNames)
      
      irf.check <- shock%*%Q_bar
      
      signCheck = rep(FALSE,2)
      # non-belief shock / technology shock - positive comovement
      signCheck[1] = irf.check["outbs","ne_constructed"] > 0 & irf.check["ne_constructed","ne_constructed"] > 0
      # belief shock - negative comovement
      signCheck[2] = irf.check["ne_constructed","outbs"] < 0 & irf.check["outbs","outbs"] > 0
      
      condall <- prod(signCheck)
      icounter <- icounter + 1
    }
    
    # switch sign of second shock (belief shock)
    Q_bar[,2] <- Q_bar[,2]*(-1)
    
    if(icounter == MaxTries)
      next
    
    shock = shock %*% Q_bar
    
    impresp2 <- array(NA_real_, c(n, n, nhor))
    impresp2[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp2[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irf_store[irep,,,] <- impresp2
  }
  idx2 <- which(!is.na(irf_store[,1,1,1]))
  cat(paste0("Sign-Restrictions fulfilled: ",length(idx2),"/",thindraws,".\n"))
  irf_store[,,2,]      = irf_store[,,2,]
  irf_store[,,1,]      = irf_store[,,1,] 
  irf_sign_all[pp,,,,] = apply(irf_store[idx2,,,], c(2,3,4), quantile, c(.05, .10, .16, .50, .84, .90, .95))
}

shock <- 2
par(mfrow=c(n,length(set.p)-1),mar=c(2,2,1,1))
for(nn in 1:n){
  if(nn==1) ylim1=c(-1.0,2.0)
  if(nn==2) ylim1=c(-2.2,0.2)
  for(ss in 1:(length(set.p)-1)){
    plot.ts(irf_sign_all[ss,4,nn,shock,], col="black", xlab="", ylab="", ylim=ylim1, main=paste0(varNames[nn]," p=",set.p[ss]), xaxt="n")
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign_all[ss,1,nn,shock,],rev(irf_sign_all[ss,7,nn,shock,])),
            col = "grey80", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign_all[ss,2,nn,shock,],rev(irf_sign_all[ss,6,nn,shock,])),
            col = "grey50", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign_all[ss,3,nn,shock,],rev(irf_sign_all[ss,5,nn,shock,])),
            col = "grey30", border=NA)
    lines(irf_sign_all[length(set.p),3,nn,shock,], col="orange", lwd=1, lty=2)
    lines(irf_sign_all[length(set.p),5,nn,shock,], col="orange", lwd=1, lty=2)
    lines(irf_sign_all[length(set.p),4,nn,shock,], col="orange", lwd=2, lty=2)
    lines(irf_sign_all[ss,4,nn,shock,], col="black", lwd=2, lty=1)
    abline(h=0, col="black", lty=1)
    axis(1, at=seq(1,nhor,by=6), labels=seq(0,nhor,by=6))
  }
}

# only p=0.1 & p=0.5 & p=0.9
set.p.plot = c("0.1","0.5","0.9","Actual")
set.p.idx  = which(set.p %in% set.p.plot)

shock <- 2
if(shock == 1){
  pdf(file=paste0("./nonbeliefshocks_tails_p=",plag,"_vs_ekm",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".pdf"), width=9, height=5)
}else if(shock == 2){
  pdf(file=paste0("./beliefshocks_tails_p=",plag,"_vs_ekm",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".pdf"), width=9, height=5)
}
par(mfrow=c(n,length(set.p.plot)-1),mar=c(2.5,3.5,2,1))
for(nn in 1:n){
  ylim1 <- range(irf_sign_all[set.p.idx,,nn,shock,])
  for(ss in set.p.idx[1:3]){
    plot.ts(irf_ekm[4,nn,shock,], col="black", xlab="", ylab="", 
            ylim=ylim1, main=paste0(varNames.plot[nn]," (p=",set.p[ss],")"), axes=FALSE, cex.main=1.4)
    abline(h=pretty(ylim1), col="grey75", lwd=2)
    abline(v=seq(1,nhor,by=6), col="grey75", lwd=2)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_ekm[1,nn,shock,],rev(irf_ekm[7,nn,shock,])), col = "grey90", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_ekm[2,nn,shock,],rev(irf_ekm[6,nn,shock,])), col = "grey70", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_ekm[3,nn,shock,],rev(irf_ekm[5,nn,shock,])), col = "grey45", border=NA)
    lines(irf_ekm[4,nn,shock,], col="black", lwd=3, lty=2)
    abline(h=0, col="black", lty=1, lwd=2)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign_all[ss,1,nn,shock,],rev(irf_sign_all[ss,7,nn,shock,])),
            col = adjustcolor("red", alpha.f=0.15), border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign_all[ss,2,nn,shock,],rev(irf_sign_all[ss,6,nn,shock,])),
            col = adjustcolor("#E22929",alpha=0.35), border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign_all[ss,3,nn,shock,],rev(irf_sign_all[ss,5,nn,shock,])),
            col = adjustcolor("#D01D1D",alpha=0.25), border=NA)
    lines(irf_sign_all[ss,4,nn,shock,], col="darkred", lwd=3, lty=1)
    axis(1, at=c(1,6,11,16,21), labels=c(0,5,10,15,20), font=2, lwd=2, cex.axis=1.1)
    axis(2, at=pretty(ylim1), las=2, font=2, lwd=2, cex.axis=1.1)
    box(lwd=2)
  }
}
dev.off()

