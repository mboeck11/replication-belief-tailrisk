bvar_wishart_mn <- function(Yraw, plag = 1, args = NULL){
  #----------------------------------------INPUTS----------------------------------------------------#
  # prepare arguments
  draws = burnin = 5000; cons=TRUE; trend=FALSE; qtrend=FALSE; thin=1; Ex=NULL; kappa=NULL; save.prior=FALSE; verbose=TRUE
  argnames = c("draws","burnin","cons","trend","qtrend","thin","Ex","kappa","save.prior")
  if(!is.null(args)){
    for(aa in argnames){
      if(aa%in%names(args)) assign(aa, args[[aa]])
    }
    for(aa in names(args)){
      if(!(aa%in%argnames)) stop(paste0("Argument: ",aa, " not in list of possible arguments!"))
    }
  }
  arglist=list(Yraw=Yraw, plag=plag, draws=draws, burnin=burnin, cons=cons, trend=trend, qtrend=qtrend,
               thin=thin, Ex=Ex)
  
  #----------------------------------------PACKAGES--------------------------------------------------#
  require(GIGrvg, quietly=TRUE)
  require(Rcpp, quietly=TRUE)
  require(MASS, quietly=TRUE)
  require(mvtnorm, quietly=TRUE)
  require(LaplacesDemon, quietly=TRUE)
  
  #-------------------------------------------START--------------------------------------------------#
  Traw  <- nrow(Yraw)
  n     <- ncol(Yraw)
  K     <- n*plag
  Ylag  <- mlag(Yraw,plag)
  varNames <- colnames(Yraw)
  if(is.null(varNames)) varNames <- paste0("Y",seq(n))
  varNameslags <- NULL
  for(pp in 1:plag) varNameslags <- c(varNameslags,paste(varNames,".lag",pp,sep=""))
  colnames(Ylag) <- varNameslags
  
  texo <- FALSE; nex <- 0; Exraw <- NULL
  if(!is.null(Ex)){
    Exraw <- Ex; nex <- ncol(Exraw)
    texo <- TRUE
    ExNames <- paste0("Ex.",colnames(Exraw))
    if(is.null(ExNames)) ExNames <- paste0("Ex.",seq(1,nex))
    varNameslags <- c(varNameslags, ExNames)
  }
  
  X <- cbind(Ylag,Exraw)
  X <- X[(plag+1):nrow(X),,drop=FALSE]
  Y <- Yraw[(plag+1):Traw,,drop=FALSE]
  bigT  <- nrow(X)
  
  if(cons){
    X <- cbind(X,1)
    varNameslags <- c(varNameslags,"cons")
    colnames(X)[ncol(X)] <- "cons"
  }
  if(trend){
    X <- cbind(X,seq(1,bigT))
    varNameslags <- c(varNameslags,"trend")
    colnames(X)[ncol(X)] <- "trend"
  }
  if(qtrend){
    X <- cbind(X,seq(1,bigT)^2)
    varNameslags <- c(varNameslags,"qtrend")
    colnames(X)[ncol(X)] <- "qtrend"
  }
  
  k = ncol(X)
  q = k*n
  v = (n*(n-1))/2
  
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv <- try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- MASS::ginv(crossprod(X))
  A_OLS  <- XtXinv%*%(t(X)%*%Y)
  E_OLS  <- Y - X%*%A_OLS
  #a_OLS <- as.vector(A_OLS)
  #SSE  <-  t((Y - X%*%A_OLS))%*%(Y - X%*%A_OLS)
  SIGMA_OLS  <- crossprod(E_OLS)/(bigT-k)
  #IXY  <-   kronecker(diag(M),(t(X)%*%Y))
  
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw     <- A_OLS
  SIGMA_draw <- SIGMA_OLS
  
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  prmean    = 0
  # prior mean hyperparameter
  a_prior <- matrix(0,k,n)
  diag(a_prior) <- prmean
  
  # compute AR residuals of AR(4)
  arsig2 <- get_resid_ar(Yraw, 4)
  
  # kappa set according to BEAR toolbox
  if(is.null(kappa)){
    kappa <- c(0.3, 0.5, 1, 1000)
  }
  
  # prior variance hyperparameter
  A_prior <- matrix(NA_real_, k, n)
  for(nn in 1:n){
    A_prior[,nn] = get_minnesotaV(n, plag, nn, kappa, arsig2, k)
  }
  
  # prior hyperparameter for variance-covariance matrix
  S_prior  <- diag(c(arsig2))
  nu_prior <- n + 2
  
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  ntot  <- draws+burnin
  
  # thinning
  count <- 0
  thindraws    <- draws/thin
  thin.draws   <- seq(burnin+1,ntot,by=thin)
  arglist      <- c(arglist, thindraws=thindraws)
  
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store       <- array(NA_real_, c(thindraws, k, n))
  SIGMA_store   <- array(NA_real_, c(thindraws, n, n))
  res_store     <- array(NA_real_, c(thindraws, bigT, n))
  
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    
    SIGMAinv_draw <- solve(SIGMA_draw)
    A_priorinv    <- 1/A_prior
    V_post <- try(chol2inv(chol(diag(as.vector(A_priorinv)) + SIGMAinv_draw %x% crossprod(X))),silent=TRUE)
    if(is(V_post,"try-error")) V_post = try(solve(diag(as.vector(A_priorinv)) + SIGMAinv_draw %x% crossprod(X)), silent=TRUE)
    if(is(V_post,"try-error")) V_post = ginv(diag(as.vector(A_priorinv)) + SIGMAinv_draw %x% crossprod(X))
    a_post <- V_post %*% (diag(as.vector(A_priorinv)) %*% as.vector(a_prior) + SIGMAinv_draw %x% crossprod(X) %*% as.vector(A_OLS))
    
    A_draw <- matrix(a_post + t(chol(V_post)) %*% rnorm(q),k,n)
    dimnames(A_draw)<-list(varNameslags, varNames)
    
    #----------------------------------------------------------------------------
    # Step 2: Sample variances
    #----------------------------------------------------------------------------
    res_draw <- Y - X %*% A_draw
    
    SIGMA_draw <- rinvwishart(nu = bigT + nu_prior, S = S_prior + crossprod(res_draw))
    #----------------------------------------------------------------------------
    # Step 3: store draws
    #----------------------------------------------------------------------------
    if(irep %in% thin.draws){
      count <- count+1
      
      A_store[count,,]       <- A_draw
      SIGMA_store[count,,]   <- SIGMA_draw
      res_store[count,,]     <- res_draw
    }
    if(irep%%50==0 && verbose) cat("Round: ",irep, "/", ntot,"\n")
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,varNameslags,varNames)
  dimnames(SIGMA_store)=list(NULL,varNames,varNames)
  ret <- list(Y=Y, X=X, A=A_store, SIGMA=SIGMA_store, res=res_store,
              args=arglist)
  return(ret)
}

gen_compMat <- function(A, n, plag){
  Jm          <- matrix(0, n*plag, n)
  Jm[1:n,1:n] <- diag(n)
  
  A   <- A[1:(n*plag),,drop=FALSE]
  Cm  <- matrix(0, n*plag, n*plag)
  if(plag==1) Cm <- t(A) else {
    for(jj in 1:(plag-1)){
      Cm[(jj*n+1):(n*(jj+1)),(n*(jj-1)+1):(jj*n)] <- diag(n)
    }
  }
  bbtemp <- A[1:(n*plag),]
  splace <- 0
  for(pp in 1:plag){
    for(nn in 1:n) {
      Cm[nn,((pp-1)*n+1):(pp*n)] <- t(bbtemp[(splace+1):(splace+n),nn])
    }
    splace <- splace+n
  }
  return(list(Cm=Cm,
              Jm=Jm))
}

mlag <- function(Y,plag){
  Y <- as.matrix(Y)
  Traw <- nrow(Y)
  n <- ncol(Y)
  Ylag <- matrix(0,Traw,plag*n)
  for (ii in 1:plag){
    Ylag[(plag+1):Traw,(n*(ii-1)+1):(n*ii)] <- Y[(plag+1-ii):(Traw-ii),(1:n)]
  }
  colnames(Ylag) <- paste0(colnames(Y),".lag",rep(seq(plag),each=n))
  return(Ylag)
}

get_resid_ar <- function(Yraw,plag){
  # get dimensions
  n    = ncol(Yraw)
  bigT = nrow(Yraw)-plag
  
  # sample variance of AR(plag) process
  sigma_sq  <- matrix(0,n,1)
  for(nn in 1:n){
    Ylag_nn        = mlag(Yraw[,nn],plag)
    Ylag_nn        = Ylag_nn[(plag+1):nrow(Ylag_nn),,drop=FALSE]
    Y_nn           = Yraw[(plag+1):nrow(Yraw),nn,drop=FALSE]
    Ylag_nn        = cbind(Ylag_nn,1)
    alpha_nn       = solve(crossprod(Ylag_nn))%*%crossprod(Ylag_nn,Y_nn)
    sigma_sq[nn,1] = (1/bigT)*t(Y_nn-Ylag_nn%*%alpha_nn)%*%(Y_nn-Ylag_nn%*%alpha_nn)
  }
  
  return(sigma_sq)
}

get_minnesotaV <- function(n, plag, nn, kappa, sig2, k){
  V_nn = matrix(0, k, 1)
  
  # according to BEAR Technical Guide page 12
  
  # construct V_nn
  for(jj in 1:(n*plag)){
    lag <- ceiling(jj/n)
    
    idx <- jj %% n
    idx <- ifelse(idx==0,n,idx)
    if(idx == nn){ # own lag coefficient
      V_nn[jj,1] = (kappa[1]/(lag^kappa[3]))^2
    }else{ # cross-variable lag coeefficient
      V_nn[jj,1] = (sig2[nn]/sig2[idx])*(kappa[1]*kappa[2]/(lag^kappa[3]))^2
    }
  }
  # deterministics
  if(k>(n*plag)){
    for(jj in (n*plag+1):k){
      V_nn[jj,1] = sig2[nn]*(kappa[1]*kappa[4])^2
    }
  }
  
  return(V_nn)
}
