#------------------------------------------------------------------------#
# Narrow Replication: Estimate Belief Shocks                             #
#                                                                        #
# Maximilian Boeck, Michael Pfarrhofer                                   #
#------------------------------------------------------------------------#
rm(list=ls())


library(readxl)
library(stringr)

# load VAR functions
source("./scripts/functions/bvar_aux.R")

#--------------------------------------------------------------------------------------#
# specification
#--------------------------------------------------------------------------------------#
sample      = "original"   # "original", "extended"
plag        = 4
trend       = TRUE
qtrend      = TRUE

vars        = c("ne_first", "outbs")
varNames    = c("Nowcast Error", "Output")
n           = length(vars)

# irf
nhor        = 21
emp_percs   = c(.05, .10, .16, .50, .84, .90, .95)

SampleStart = c(1968,4)
if(sample == "original") SampleEnd = c(2014,4)
if(sample == "extended") SampleEnd = c(2019,4)

#--------------------------------------------------------------------------------------#
# read data
#--------------------------------------------------------------------------------------#
spf_median = read_xlsx("./data/medianGrowth.xlsx", sheet="RGDP")
rgdp       = read_xlsx("./data/routput_first_second_third.xlsx", sheet="DATA", skip=3, col_names=TRUE)
outbs      = read_xls("../02_data/FRED/OUTBS.xls", skip=10)
rgdp_first = ts(as.numeric(rgdp$First), start=c(1965,3), frequency=4)
rgdp_third = ts(as.numeric(rgdp$Third), start=c(1965,3), frequency=4)
spf_median = ts(spf_median$drgdp2, start=c(1968,4), frequency=4)
outbs      = ts(log(outbs[,2])*100, start=c(1968,4), frequency=4)

data = ts.intersect(rgdp_first, rgdp_third, spf_median, outbs)
data = cbind(data, ne_first=(data[,"rgdp_first"]-data[,"spf_median"])/4, ne_third=(data[,"rgdp_third"]-data[,"spf_median"])/4)
colnames(data) = str_remove(colnames(data), "data\\.")

rm(spf_median, rgdp, outbs, rgdp_first, rgdp_third)

#--------------------------------------------------------------------------------------#
# estimate the VAR
#--------------------------------------------------------------------------------------#
Yraw           = window(data[,vars], start=SampleStart, end=SampleEnd)
colnames(Yraw) = varNames

args_uninf = list(draws = 10000, burnin = 5000, thin=2, cons=TRUE, trend=trend, qtrend=qtrend, kappa=c(100^2,1,1,100^2))
args_mn    = list(draws = 10000, burnin = 5000, thin=2, cons=TRUE, trend=trend, qtrend=qtrend, kappa=NULL)

# Bayesian VAR
run_uninf  = bvar_wishart_mn(Yraw, plag = plag, args=args_uninf)
run_uninf2 = bvar_wishart_mn(Yraw, plag = 2, args=args_uninf)
run_uninf3 = bvar_wishart_mn(Yraw, plag = 6, args=args_uninf)
run_mn     = bvar_wishart_mn(Yraw, plag = plag, args=args_mn)

# add no SV
run_uninf$args$SV = run_uninf2$args$SV = run_uninf3$args$SV = run_mn$args$SV = FALSE

mods = c("uninf", "uninf2", "uninf3", "mn")

#--------------------------------------------------------------------------------------#
# compute IRFs
#--------------------------------------------------------------------------------------#
irf_sign_big = array(NA_real_, c(length(emp_percs), n, n, nhor, length(mods)))
for(mm in 1:length(mods)){
  run = get(paste0("run_",mods[mm]))
  
  thindraws = run$args$thindraws
  bigT      = nrow(run$args$Yraw) - run$args$plag
  plag      = run$args$plag
  
  # other stuff
  MaxTries <- 7500
  
  # horizon
  H.restr <- 1
  N.restr <- H.restr*n
  
  Smat <- vector("list", as.integer(n)); Zmat <- vector("list", as.integer(n))
  # first shock: non-belief shock / technology shock - positive comovement
  Smat[[1]] <- matrix(0, 2, n*H.restr) 
  Smat[[1]][1,1] <- 1; Smat[[1]][2,2] <- 1
  
  # second shock: belief shock / demand shock - negative comovement
  Smat[[2]] <- matrix(0, 2, n*H.restr) 
  Smat[[2]][1,1] <- -1; Smat[[2]][2,2] <- 1
  
  irf_store <- array(NA_real_, c(thindraws, n, n, nhor),
                     dimnames=list(NULL, varNames, varNames, seq(0,nhor-1)))
  eps_store <- array(NA_real_, c(thindraws, bigT, n))
  for(irep in 1:thindraws){
    
    if(irep%%50==0) print(paste0("Round: ",irep))
    
    temp    <- gen_compMat(A=run$A[irep,,], n=n, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    if(run$args$SV){
      SIGMA <- apply(run$SIGMA[irep,,,], c(2,3), median)
    }else{
      SIGMA   <- run$SIGMA[irep,,]
    }
    res <- run$res[irep,,]
    
    if(max(abs(Re(eigen(compMat)$values))) > 1)
      next
    
    A0 <- t(chol(SIGMA))
    irf.restr       <- matrix(NA, N.restr, n)
    irf.restr[1:n,] <- A0
    compMati        <- compMat
    if(H.restr > 1){
      for(hh in 2:H.restr){
        irf.restr[((hh-1)*n+1):(hh*n),] <- t(Jm) %*% compMati %*% Jm %*% A0
        compMati <- compMati %*% compMat
      }
    }
    colnames(irf.restr) <- varNames
    rownames(irf.restr) <- paste(rep(varNames,H.restr),".",
                                 rep(seq(0,H.restr-1),each=n),sep="")
    
    # draw rotation matrix here
    icounter <- 0
    condall <- 0
    while(condall == 0 && icounter < MaxTries){
      randMat = matrix(rnorm(n^2),n,n)
      QR      = qr(randMat)
      Q       = qr.Q(QR)
      Q_bar   = Q%*%diag(((diag(Q)>0)-(diag(Q)<0)))
      
      irf.check <- irf.restr %*% Q_bar
      
      check_all <- 0
      for(nn in 1:n){
        S.temp <- Smat[[nn]]
        if(is.null(S.temp)){
          check_all <- check_all + 1
          next
        }
        nidx <- nrow(S.temp)
        check1 <- S.temp %*% sign(irf.check[,nn])
        check2 <- S.temp %*% sign(irf.check[,nn])*(-1)
        if(sum(check1)  == nidx){
          check_all <- check_all + 1
        }else if(sum(check2) == nidx){
          Q[,nn] = -Q[,nn]
          check_all <- check_all + 1
        }else{
          next
        }
      }
      
      if(check_all == n) condall <- 1
    }
    
    #  switch sign of second shock (belief shock)
    Q_bar[,2] <- Q_bar[,2]*(-1)
    
    if(icounter == MaxTries)
      next
    
    shock = A0 %*% Q_bar
    eps   = res %*% t(solve(shock))
    
    impresp2 <- array(NA_real_, c(n, n, nhor))
    impresp2[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp2[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irf_store[irep,,,] <- impresp2
    eps_store[irep,,]  <- eps
  }
  
  idx2 <- which(!is.na(irf_store[,1,1,1]))
  cat(paste0("Sign-Restrictions fulfilled: ",length(idx2),"/",thindraws,".\n"))
  irf_sign  <- apply(irf_store[idx2,,,], c(2,3,4), quantile, emp_percs)
  eps_mean  <- ts(apply(eps_store[idx2,,], c(2,3), mean), start=c(1970,4), frequency=4)
  
  irf_sign_big[,,,,mm] = irf_sign
  rm(irf_sign)
}

if(sample == "original"){
  irf_ekm = irf_sign_big[,,,,1]
  save(irf_ekm, file="./results/irf_ekm.rda")
}

# get first element
irf_sign = irf_sign_big[,,,,1]

par(mfrow=c(n,n),mar=c(2,3,2,2))
for(nn in 1:n){
  for(shock in 1:n){
    ylim1=range(irf_sign[,nn,shock,])
    plot.ts(irf_sign[4,nn,shock,], col="black", xlab="", ylab="", ylim=ylim1, main=varNames[nn], axes=FALSE)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[1,nn,shock,],rev(irf_sign[7,nn,shock,])),
            col = "grey70", border=NA)
    # polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[2,nn,shock,],rev(irf_sign[6,nn,shock,])),
    #         col = "grey50", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[3,nn,shock,],rev(irf_sign[5,nn,shock,])),
            col = "grey40", border=NA)
    lines(irf_sign[4,nn,shock,], col="black", lwd=2, lty=1)
    abline(h=0, col="red", lty=2)
    axis(1, at=seq(1,nhor,by=6), labels=seq(0,nhor,by=6), lwd=2, font=2)
    axis(2, lwd=2, font=2, las=2)
    box(lwd=2)
  }
}

pdf(file=paste0("./beliefshocks_linear_",sample,".pdf"), width=8, height=6)
par(mfrow=c(3,3),mar=c(2.2,2.9,1,1))
par(fig=c(0,0.1,0.90,1))
plot(-10,-10, axes=FALSE,ylim=c(0,1),xlim=c(0,1))
par(fig=c(0,0.1,0.45,0.90),new=TRUE)
plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1))
text(0.5,0.5, "Nowcast Error",srt=90,cex=2, font=2)
par(fig=c(0,0.1,0,0.45),new=TRUE)
plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1))
text(0.5,0.5, "Output", srt=90, cex=2, font=2)
par(fig=c(0.1,0.55,0.90,1),new=TRUE)
plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1))
text(0.5,0.5, "Non-Belief Shock", cex=2, font=2)
shock<-1
for(nn in 1:n){
  par(fig=c(0.1,0.55,c(0.45,0.00)[nn],c(0.90,0.45)[nn]),new=TRUE)
  ylim1=range(irf_sign[,nn,shock,])
  if(nn == 1) ylim1 = c(-0.5, 0.5)
  if(nn == 2) ylim1 = c(-1.5, 1.5)
  plot.ts(irf_sign[4,nn,shock,], col="black", xlab="", ylab="", ylim=ylim1, main="", axes=FALSE)
  abline(h=pretty(ylim1), col="grey75", lwd=2)
  abline(v=seq(1,nhor,by=6), col="grey75", lwd=2)
  polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[1,nn,shock,],rev(irf_sign[7,nn,shock,])),
          col = "grey90", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[2,nn,shock,],rev(irf_sign[6,nn,shock,])),
          col = "grey70", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[3,nn,shock,],rev(irf_sign[5,nn,shock,])),
          col = "grey45", border=NA)
  lines(irf_sign[4,nn,shock,], col="black", lwd=2, lty=1)
  for(mm in 2:length(mods)){
    lines(irf_sign_big[4,nn,shock,,mm], col="grey20", lwd=2, lty=mm)
  }
  abline(h=0, col="black", lwd=2)
  axis(1, at=c(1,6,11,16,21), labels=c(0,5,10,15,20), lwd=2, font=2, cex.axis=1.5)
  if(nn == 1) axis(2, at=c(-0.5,0,0.5), lwd=2, font=2, las=2, cex.axis=1.5)
  if(nn == 2) axis(2, at=c(-1.5,-1.0,-0.5,0,0.5,1,1.5), lwd=2, font=2, las=2, cex.axis=1.5)
  box(lwd=2)
}
par(fig=c(0.55,1,0.90,1),new=TRUE)
plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1))
text(0.5,0.5, "Belief Shock", cex=2, font=2)
shock<-2
for(nn in 1:n){
  par(fig=c(0.55,1,c(0.45,0.00)[nn],c(0.90,0.45)[nn]),new=TRUE)
  ylim1=range(irf_sign[,nn,shock,])
  if(nn == 1) ylim1 = c(-0.5, 0.5)
  if(nn == 2) ylim1 = c(-1.5, 1.5)
  plot.ts(irf_sign[4,nn,shock,], col="black", xlab="", ylab="", ylim=ylim1, main="", axes=FALSE)
  abline(h=pretty(ylim1), col="grey75", lwd=2)
  abline(v=seq(1,nhor,by=6), col="grey75", lwd=2)
  polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[1,nn,shock,],rev(irf_sign[7,nn,shock,])),
          col = "grey90", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[2,nn,shock,],rev(irf_sign[6,nn,shock,])),
          col = "grey70", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irf_sign[3,nn,shock,],rev(irf_sign[5,nn,shock,])),
          col = "grey45", border=NA)
  lines(irf_sign[4,nn,shock,], col="black", lwd=2, lty=1)
  for(mm in 2:length(mods)){
    lines(irf_sign_big[4,nn,shock,,mm], col="grey20", lwd=2, lty=mm)
  }
  abline(h=0, col="black", lwd=2)
  axis(1, at=c(1,6,11,16,21), labels=c(0,5,10,15,20), lwd=2, font=2, cex.axis=1.5)
  if(nn == 1) axis(2, at=c(-0.5,0,0.5), lwd=2, font=2, las=2, cex.axis=1.5)
  if(nn == 2) axis(2, at=c(-1.5,-1.0,-0.5,0,0.5,1,1.5), lwd=2, font=2, las=2, cex.axis=1.5)
  box(lwd=2)
}
dev.off()

