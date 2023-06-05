library(mclust)
library(mvtnorm)
library(zoo)
library(snipEM)

## Perform Model-based clustering to obtain initial values for all the
## parameters of interest, including the sequence of latent states

getInits=function(y, Z, kmax,b) {
  
  initials=list()
  Ti=dim(y)[2]
  r=dim(y)[3]
  n=dim(y)[1]
  n.cov=dim(Z)[3]
  initials$k=rep(NA,Ti)
  initials$xi=initials$sigma=array(NA,c(kmax,kmax,r))
  initials$be=array(b,c(n.cov,kmax,kmax))
  initials$be[,,1]=0
  initials$Be=array(b,c(n.cov,kmax,kmax,kmax,kmax))
  for(r1 in 1:kmax){
  for(r2 in 1:kmax){
  if(r1==r2){
  for(t1 in 1:r1){
    initials$Be[,r1,r2,t1,t1]=0
  }}else if(r1>r2){
    ref=c(1:r2)
    while (length(ref)<r1){
      ref=c(ref,max(ref))
    }
  for(t1 in 1:r1){
    initials$Be[,r1,r2,ref[t1],t1]=0
    initials$Be[,r1,r2,setdiff(c(1:r1),c(1:r2)),1:r2]=0
  }}else if(r1<r2){
    ref=c(1:r1)
    while(length(ref)<r2){
      ref=c(ref,max(ref))
    }
    for(t1 in 1:r2){
      initials$Be[,r1,r2,ref[t1],t1]=0
    }
  } 
  }}
  
  for(ti in 1:Ti) {
    mc=Mclust(y[,ti,])
    initials$k[ti]=min(c(mc$G,kmax))
  }
  
  for(k in 1:kmax) {
    jnk=NULL
    if(any(initials$k==k)) {
      for(w in which(initials$k==k)) {jnk=rbind(jnk,y[,w,])}}
    if(!any(initials$k==k)) {jnk=y[,Ti,]}
    mc=Mclust(jnk,G=k,modelNames="VVI")
    if(is.null(mc)) {mc=Mclust(jnk,G=k,modelNames="VVV")}
    if(is.null(mc)) {mc=Mclust(jnk,G=k,modelNames="VEV")}
    if(is.null(mc)) {mc=Mclust(jnk,G=k,modelNames="EVV")}
    if(is.null(mc)) {mc=Mclust(jnk,G=k,modelNames="EEV")}
    if(is.null(mc)) {mc=Mclust(jnk,G=k,modelNames="VVE")}
    if(is.null(mc)) {mc=Mclust(jnk,G=k,modelNames="EVE")}
    initials$xi[k,1:k,]=t(mc$parameters$mean)
    for(w in 1:k) {
      initials$sigma[k,w,]=sqrt(diag(mc$parameters$variance$sigma[,,w]))}
    
    o=order(initials$xi[k,1:k,1])
    initials$xi[k,1:k,1]=sort(initials$xi[k,1:k,1])
    initials$sigma[k,1:k,1]=initials$sigma[k,o,1]
    for(h in 2:r) {
      initials$xi[k,1:k,h]=initials$xi[k,o,h]
      initials$sigma[k,1:k,h]=initials$sigma[k,o,h]    
    }}
  
  return(initials)}


## The two following functions are used to map the coefficients and
## the covariates into initial and trasition probabilities


# INITIAL PROBABILITIES #

inprob <- function(K, be, z, thB){
  
  if(K==1){
    ip=1
  }else{
    ip=exp(z%*%be)
  }
  ip = ip/sum(ip)
  ip[ip<thB]=thB
  return((ip/sum(ip)))
}

# TRANSITION PROBABILITIES #

tmat <- function(K, be, z, thB){
  if(length(z)==1){
    # solo intercetta
    if(K[2]==1){
      return(matrix(1,K[1],K[2]))
    }
    if(K[1]==1){
      tmatr <- exp(z%*%be)
      return(tmatr/sum(tmatr))
    }
    tmatr <- exp(be)
    tmatr = sweep(tmatr, 1,apply(tmatr,1,sum), "/")
    tmatr[tmatr<thB]=thB
    return(sweep(tmatr, 1,apply(tmatr,1,sum), "/"))
  }else{
    # con covariate
    if(K[2]==1){
      return(matrix(1,K[1],K[2]))
    }
    else if(K[1]==1){
      tmatr <- exp(z%*%be)
      return(tmatr/sum(tmatr))
    }
    else if(K[1] != K[2]){
      tmatr <- t(exp(apply(be,2,crossprod,z)))
      tmatr = sweep(tmatr, 1,apply(tmatr,1,sum), "/")
      tmatr[tmatr<thB]=thB
      return(sweep(tmatr, 1,apply(tmatr,1,sum), "/"))
    }
    else { 
      tmatr <- exp(apply(be,2,crossprod,z))
      tmatr = sweep(tmatr, 1,apply(tmatr,1,sum), "/")
      tmatr[tmatr<thB]=thB
      return(sweep(tmatr, 1,apply(tmatr,1,sum), "/"))
    }
  }
}

# This function performs the f/b recursion to calculate the likelihood
likco=function(y,w,xi,sigma,pi,PI,k,kmax,n,Ti){
  
  qu=matrix(NA,Ti,kmax)
  
  for(j in 1:k[1])
  {qu[1,j]=w*dmvnorm(y[1,],xi[k[1],j,],diag(sigma[k[1],j,]^2),log=TRUE)+log(pi[k[1],j])}
  
  for(t in 2:Ti) {
    for(j in 1:k[t]) {
      if(k[t-1]>1) {
        qu[t,j]=(qu[t-1,1]+log(PI[(t-1),k[t-1],k[t],1,j]))
        for(d in 2:k[t-1]) {
          cb=c(qu[t,j],(qu[t-1,d]+log(PI[(t-1),k[t-1],k[t],d,j])))
          cb[ cb < -700]=-700
          qu[t,j]=log(sum(exp(cb)))}
        qu[t,j]=qu[t,j]+w*dmvnorm(y[t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)}
      
      
      if(k[t-1]==1) {
        qu[t,j]=w*dmvnorm(y[t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)+qu[t-1,1]+log(PI[(t-1),k[t-1],k[t],1,j])
      }
      
    }}
  
  if(k[Ti]>1) {
    cb=qu[Ti,1:k[Ti]]
    cb[cb < -700]=-700
    liks=log(sum(exp(cb)))}
  if(k[Ti]==1) {liks=qu[Ti,1]}
  res=sum(liks)
}

## This is the main function performing reversible jumps mcmc
## for a latent markov model with unknown number of regime profiles

rlmrjmcmc <- function(y,w=1,Z,n.iter,ki,kmax,inits,initials,verbose,msi=2,vsi=2,mxi=0,vxi=1,Tauxi=rep(0.01,r),mbe=0,vbe=1,mBe=0,vBe=1.5,mu=0,vu=1,mw=2,vw=2,mepsb=0,vepsb=1,mepsB=0,vepsB=1,thB=10^-5,cutOffRj=Inf){
  
  
  n = nrow(y)
  Ti = ncol(y)
  r = dim(y)[3]
  pseq = 2*Ti + 4
  n.cov = dim(Z)[3]
  
  liks <- matrix(NA, n.iter, pseq)
  
  k = matrix(NA, n.iter, Ti)
  Be <- array(NA, dim = c(n.iter,n.cov,kmax,kmax,kmax,kmax))
  BeU <- BeC <- array(NA, c(n.cov,kmax,kmax,kmax,kmax))
  PIs <- array(NA, c(n,Ti-1, kmax, kmax, kmax, kmax))
  PIsU <- PIs
  be <- array(NA, c(n.iter,n.cov,kmax,kmax))
  beU <- beC <- array(NA,c(dim(be)[-1]))
  pis <- pisU <- array(NA, c(n,kmax,kmax))
  xi <- array(NA, dim = c(n.iter,kmax,kmax,r))
  xiU <- xiC <- sigmaU <- sigmaC <- array(NA,c(kmax,kmax,r))
  sigma <- array(NA, dim = c(n.iter,kmax,kmax,r))
  o = matrix(NA, kmax, kmax)
  
  if(!isTRUE(inits)){
    k[1,] = rep(median(1:kmax),Ti)
    centers<-apply(y[,1,],2,mean)
    sds<-apply(y[,1,],2,mad)/max(k[1,])
    for(tt in unique(k[1,])){
      xi[1,tt,1:tt,1] <- centers[1]+sds[1]*(seq(-1.5,1.5,length=tt))
      o[tt,1:tt]=order(xi[1,tt,1:tt,1])
      xi[1,tt, 1:tt,1] = sort(xi[1,tt,1:tt,1])
      for(h in (1:dim(y)[3])[-1]) {
        xi[1,tt, 1:tt,h] <- centers[h]+sds[h]*(seq(-1.5,1.5,length=tt))[o[tt,1:tt]]
      }
      for(h in 1:3){
        sigma[1,tt,1:tt,h] <- sds[h]
      }
    }
    for(tt in unique(k[1,])){
      be[1,,tt,1:tt] <- 0
    }
    for(t1 in unique(k[1,])){
      for(t2 in unique(k[1,])){
        Be[1,,t1,t2,1:t1,1:t2] <- 0
      }
    }
    for(i in 1:n){
      for(tt in 2:length(k[1,])){
        PIs[i,(tt-1),k[1,tt-1],k[1,tt],1:k[1,tt-1],1:k[1,tt]]= tmat(K = c(k[1,tt-1],k[1,tt]), Be[1,,k[1,tt-1],k[1,tt],1:k[1,tt-1],1:k[1,tt]],Z[i,tt,],10^-5)
        pis[i,k[1,1],1:k[1,1]] <- inprob(k[1,1], be=be[1,,k[1,1],1:k[1,1]], Z[i,1,], thB = 10^-5)
      }
    }
  }
  if(isTRUE(inits)){
    k[1,] = ki
    for(t1 in unique(k[1,])){
      xi[1,t1,1:t1,] <- initials$xi[t1,1:t1,]
      sigma[1,t1,1:t1,] <- initials$sigma[t1,1:t1,]
    }
    be[1,,k[1,1],1:k[1,1]] <- initials$be[,k[1,1],1:k[1,1]]
    for(tt in 1:(Ti-1)){
      Be[1,,k[1,tt],k[1,tt+1],1:k[1,tt],1:k[1,tt+1]] = initials$Be[,k[1,tt],k[1,tt+1],1:k[1,tt],1:k[1,tt+1]]
    }
    for(i in 1:n){
      pis[i,k[1,1],1:k[1,1]] <- inprob(k[1,1], be=be[1,,k[1,1],1:k[1,1]], Z[i,1,], thB = 10^-5)
      for(tt in 2:length(k[1,])){
        PIs[i,(tt-1),k[1,tt-1],k[1,tt],1:k[1,tt-1],1:k[1,tt]]= tmat(K = c(k[1,tt-1],k[1,tt]), Be[1,,k[1,tt-1],k[1,tt],1:k[1,tt-1],1:k[1,tt]],Z[i,tt,],10^-5)
      }
    }
  }
  ilik = sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xi[1,,,], sigma = sigma[1,,,], pi = pis[i,,], PI = PIs[i,,,,,], k=k[1,], kmax=kmax, Ti=Ti)) 
  liks[1,] = sum(ilik)    
  
  for(iter in 2:n.iter){
    
    k[iter,] <- k[iter-1,]
    
    if(iter>cutOffRj) { 
      
      xi[iter,,,] <- xi[iter-1,,,]
      sigma[iter,,,] <- sigma[iter-1,,,]
      Be[iter,,,,,] <- Be[iter-1,,,,,]
      be[iter,,,] <- be[iter-1,,,]
      liks[iter,1:(2*Ti)] <- liks[iter-1,2*Ti+4]
      
    }
    if(iter<=cutOffRj) {
      ###################  
      ### BIRTH/DEATH ###
      ###################
      if(verbose){print("## birth/death")}
      
      ku = kc = k[iter,]
      xiC <- xi[iter-1,,,]
      sigmaC <- sigma[iter-1,,,]
      BeC <- array(Be[iter-1,,,,,],c(dim(Be)[-1]))
      xiU = sigmaU <- array(NA,c(kmax,kmax,r))
      BeU <- array(NA, c(n.cov,kmax,kmax,kmax,kmax))
      beU <- array(NA,c(dim(be)[-1]))
      # t = 1 #
      if(k[iter,1]==1){ 
        bd1 = 0
      }else if(k[iter,1]==kmax){
        bd1=1
      }else{bd1 = runif(1)}
      
      if(bd1 <=0.5){
        # BIRTH #
        xij0 = sij0 = bej0 = Bej0 = NA
        ku[1] = kc[1]+1 
        j0 = sample(c(1:ku[1]),1) # sample the location
        beC = array(be[iter-1,,,],c(dim(be)[-1]))
        beU[,ku[1],(1:ku[1])[-j0]] = beC[,kc[1],1:kc[1]] 
        bej0 = beU[,ku[1],j0] = rnorm(n.cov,mbe,vbe) # sample parameters from priors
        ksub = ku[-1] 
        L = !any(ksub==ku[1]) # check need for update
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        if(L==1){ # the update step for centroids and sds
          xiU[ku[1],(1:ku[1])[-j0],] = xiC[kc[1],1:kc[1],]
          xij0 = xiU[ku[1],j0,] = rnorm(r,mxi,vxi)
          sigmaU[ku[1],(1:ku[1])[-j0],] = sigmaC[kc[1],1:kc[1],]
          sij0 = sigmaU[ku[1],j0,] = rgamma(r,msi,vsi)
          o[ku[1],1:ku[1]] = order(xiU[ku[1],1:ku[1],1])
          xiU[ku[1],1:ku[1],1]=sort(xiU[ku[1],1:ku[1],1])
          for(h in 2:dim(y)[3]) {
            xiU[ku[1],1:ku[1],h] <- xiU[ku[1],o[ku[1],1:ku[1]],h]
          }
        }
        W = !any(rollapply(kc, 2, identical, c(ku[1],ku[2]))) # check if transitions need updates
        for(tt in 1:(Ti-2)){
          BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
        }
        if(W==1){ # update step for coefficients modulating transition k_1 to k_2
          BeU[,ku[1],ku[2],(1:ku[1])[-j0],1:ku[2]] = BeC[,kc[1],kc[2],1:kc[1],1:kc[2]]
          Bej0 = BeU[,ku[1],ku[2],j0,1:ku[2]] = rnorm(n.cov*kc[2],mBe,vBe)
          
          for(j in 1:ku[1]) {
            BeU[,ku[1],ku[2],j,1:ku[2]] = BeU[,ku[1],ku[2],j,1:ku[2]] - BeU[,ku[1],ku[2],j,min(j,ku[2])]}}
      }else{
        # DEATH #
        xij0 = sij0 = bej0 = Bej0 = NA
        ku[1] = kc[1]-1
        j0 = sample(c(1:kc[1]),1) # sample the location
        beC = array(be[iter-1,,,],c(dim(be)[-1]))
        bej0 = beC[,kc[1],j0]
        beU[,ku[1],(1:ku[1])] = beC[,kc[1],(1:kc[1])[-j0]] # delete the sampled coefficients
        ksub = ku[-1] 
        L = !any(ksub==ku[1]) # check need for update
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        if(L==1){ # update step for centroids and sds
          xij0 = xiC[kc[1],j0,]
          sij0 = sigmaC[kc[1],j0,]
          xiU[ku[1],1:ku[1],] <- xiC[kc[1],(1:kc[1])[-j0],]
          sigmaU[ku[1],1:ku[1],] <- sigmaC[kc[1],(1:kc[1])[-j0],]
          o[ku[1],1:ku[1]] = order(xiU[ku[1],1:ku[1],1])
          xiU[ku[1],1:ku[1],1]=sort(xiU[ku[1],1:ku[1],1])
          for(h in 2:dim(y)[3]) {
            xiU[ku[1],1:ku[1],h] <- xiU[ku[1],o[ku[1],1:ku[1]],h]
          }
        }
        W =!any(rollapply(kc, 2, identical, c(ku[1],ku[2]))) # check if transitions need updates
        for(tt in 1:(Ti-2)){
          BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
        }
        if(W==1){ # update step for coefficients modulating transition k_1 to k_2
          Bej0 = BeC[,kc[1],kc[2],j0,1:kc[2]]
          BeU[,ku[1],ku[2],1:ku[1],1:ku[2]] = BeC[,kc[1],kc[2],(1:kc[1])[-j0],1:kc[2]]
          for(j in 1:ku[1]) {
            
            BeU[,ku[1],ku[2],j,1:ku[2]] = BeU[,ku[1],ku[2],j,1:ku[2]] -
              BeU[,ku[1],ku[2],j,min(j,ku[2])]}
        }
      }
      Be[iter,,,,,] <- BeU
      xi[iter,,,] <- xiU
      be[iter,,,] <- beU
      sigma[iter,,,] <- sigmaU
      k[iter,] = ku
      for(i in 1:n){ # update initial and transition probs
        pisU[i,k[iter,1],1:k[iter,1]] <- inprob(k[iter,1], be=beU[,k[iter,1],1:k[iter,1]], Z[i,1,], thB = 10^-5)
        for(tt in 2:length(k[1,])){
          PIsU[i,(tt-1),k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]]= tmat(K = c(k[iter,tt-1],k[iter,tt]), Be[iter,,k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]],Z[i,tt,],10^-5)        
        }
      }
      ilik = sapply(1:n,function(i) likco(y = y[i,,], w=1, xi = xiU[,,], sigma = sigmaU[,,], pi = pisU[i,,], PI = PIsU[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
      liks[iter,1] = sum(ilik) # calculate likelihood
      if(bd1 <= 0.5){
        # accept a birth
        if(!is.finite(liks[iter,1])){
          Ab1 = 0
        }else{
          Ab1 = (liks[iter,1]-liks[iter-1,pseq])+((sum(dnorm(beU,mbe,vbe,log = T),na.rm = T) + sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T)) - (sum(dnorm(beC,mbe,vbe,log = T),na.rm = T) + sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
            log((1/(kc[1]-1))/(1/ku[1])) +     
            (log(ku[1]) - (sum(dnorm(xij0,mxi,vxi,log = T),na.rm = T)*L +(sum(dgamma(sij0,msi,vsi,log = T),na.rm = T))*L+(sum(dnorm(bej0,mbe,vbe,log = T),na.rm = T))+(sum(dnorm(Bej0,mBe,vBe,log = T),na.rm = T))*W))
          Ab1 = exp(Ab1)}
        ab = runif(1)
        if(ab > min(1,Ab1)){
          liks[iter,1] = liks[iter-1,pseq]
          ku = kc
          k[iter,] = kc
          xiU = xiC
          xi[iter,,,] = xiC
          sigmaU = sigmaC
          Be[iter,,,,,] <- BeC
          sigma[iter,,,] = sigmaC
          be[iter,,,] = beC
          BeU = BeC
        }else{
          be[iter,,,] = beU
          pis = pisU
          PIs = PIsU
        }
      }else{
        # accept a death
        if(!is.finite(liks[iter,1])){
          Ab1 = 0
        }else{
          Ab1 = (liks[iter,1]-liks[iter-1,pseq])+((sum(dnorm(beU,mbe,vbe,log = T),na.rm = T) + sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T)) - (sum(dnorm(beC,mbe,vbe,log = T),na.rm = T) + sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
            log((1/ku[1])/(1/(kc[1]-1))) +     
            (-(sum(dnorm(xij0,mxi,vxi,log = T),na.rm = T)*L +(sum(dgamma(sij0,msi,vsi,log = T),na.rm = T))*L+(sum(dnorm(bej0,mbe,vbe,log = T),na.rm = T))+(sum(dnorm(Bej0,mBe,vBe,log = T),na.rm = T))*W)+log(ku[1]))
          Ab1 = exp(Ab1)}
        ab = runif(1)
        if(ab>min(1,Ab1)){
          liks[iter,1] = liks[iter-1,pseq]
          ku = kc
          xiU = xiC
          xi[iter,,,] = xiC
          k[iter,] = kc
          sigmaU = sigmaC
          sigma[iter,,,] = sigmaC
          be[iter,,,] = beC
          BeU = BeC
          Be[iter,,,,,] <- BeC
        }else{
          be[iter,,,] = beU
          pis = pisU
          PIs = PIsU
        }
      }
      # 1 < t < Ti
      for(t in 2:(Ti-1)){
        BeC = array(Be[iter,,,,,],c(dim(Be)[-1]))
        xiC = xi[iter,,,]
        sigmaC = sigma[iter,,,]
        xiU = sigmaU <- array(NA,c(kmax,kmax,r))
        BeU <- array(NA, c(n.cov,kmax,kmax,kmax,kmax))
        ku = kc <- k[iter,]
        if(kc[t]==1){
          bdt = 0
        }else if(kc[t]==kmax){
          bdt=1
        }else{bdt = runif(1)}
        if(bdt <= 0.5){
          # BIRTH
          xij0 = sij0 = Bej0F = Bej0B = NA
          ku[t] = kc[t]+1
          j0 = sample(c(1:ku[t]),1) # sample the location
          ksub = ku[-t] 
          L = !any(ksub==ku[t])
          for(t1 in unique(ku)){
            xiU[t1,1:t1,] = xiC[t1,1:t1,]
            sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
          }
          if(L==1){ # update step for centroids and sds
            xiU[ku[t],(1:ku[t])[-j0],] = xiC[kc[t],1:kc[t],]
            xij0 = xiU[ku[t],j0,] = rnorm(r,mxi,vxi)
            sigmaU[ku[t],(1:ku[t])[-j0],] = sigmaC[kc[t],1:kc[t],]
            sij0 = sigmaU[ku[t],j0,] = rgamma(r,msi,vsi)
            o[ku[t],1:ku[t]] = order(xiU[ku[t],1:ku[t],1])
            xiU[ku[t],1:ku[t],1]=sort(xiU[ku[t],1:ku[t],1])
            for(h in 2:dim(y)[3]) {
              xiU[ku[t],1:ku[t],h] <- xiU[ku[t],o[ku[t],1:ku[t]],h]
            }
          }
          # check if the transition k_{t} to k_{t+1} need update
          W =!any(rollapply(kc, 2, identical, c(ku[t],ku[t+1])))
          for(tt in 1:(Ti-2)){
            BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
          }
          if(W==1){ # update step for coefficients modulating transition k_{t} to k_{t+1}
            BeU[,ku[t],ku[t+1],(1:ku[t])[-j0],1:ku[t+1]] = BeC[,kc[t],kc[t+1],1:kc[t],1:kc[t+1]]
            Bej0F = BeU[,ku[t],ku[t+1],j0,1:ku[t+1]] = rnorm(n.cov*ku[t+1],mBe,vBe)
            for(j in 1:ku[t]) {
              BeU[,ku[t],ku[t+1],j,1:ku[t+1]] = BeU[,ku[t],ku[t+1],j,1:ku[t+1]] - BeU[,ku[t],ku[t+1],j,min(j,ku[t+1])]}
          }
          # check if the transition k_{t-1} to k_{t} needs update
          H = !any(rollapply(kc, 2, identical, c(ku[t-1],ku[t])))
          if(H==1){ # update step for coefficients modulating transition k_{t-1} to k_{t}
            BeU[,ku[t-1],ku[t],1:ku[t-1],(1:ku[t])[-j0]] = BeC[,kc[t-1],kc[t],1:kc[t-1],1:kc[t]]
            Bej0B = BeU[,ku[t-1],ku[t],1:ku[t-1],j0] = rnorm(n.cov*ku[t-1],mBe,vBe)
            for(j in 1:ku[t-1]) {
              BeU[,ku[t-1],ku[t],j,1:ku[t]] = BeU[,ku[t-1],ku[t],j,1:ku[t]] - BeU[,ku[t-1],ku[t],j,min(j,ku[t])]}
          }
        }else{
          # DEATH 
          xij0 = sij0 = Bej0F = Bej0B = NA
          ku[t] = kc[t]-1
          ksub = ku[-t] 
          j0 = sample(c(1:kc[t]),1) # sample the location
          L = !any(ksub==ku[t])
          for(t1 in unique(ku)){
            xiU[t1,1:t1,] = xiC[t1,1:t1,]
            sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
          }
          if(L==1){ # update step for centroids and sds
            xij0 = xiC[kc[t],j0,]
            sij0 = sigmaC[kc[t],j0,]
            xiU[ku[t],1:ku[t],] <- xiC[kc[t],(1:kc[t])[-j0],]
            sigmaU[ku[t],1:ku[t],] <- sigmaC[kc[t],(1:kc[t])[-j0],]
            o[ku[t],1:ku[t]] = order(xiU[ku[t],1:ku[t],1])
            xiU[ku[t],1:ku[t],1]=sort(xiU[ku[t],1:ku[t],1])
            for(h in 2:dim(y)[3]) {
              xiU[ku[t],1:ku[t],h] <- xiU[ku[t],o[ku[t],1:ku[t]],h]
            }
          }
          # check if transition k_{t-1} to k_{t} needs update
          W =!any(rollapply(kc, 2, identical, c(ku[t],ku[t+1])))
          for(tt in 1:(Ti-2)){
            BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
          }
          if(W==1){ # update step for coefficients modulating transition k_{t} to k_{t+1}
            Bej0F = BeC[,kc[t],kc[t+1],j0,1:kc[t+1]]
            BeU[,ku[t],ku[t+1],1:ku[t],1:ku[t+1]] = BeC[,kc[t],kc[t+1],(1:kc[t])[-j0],1:kc[t+1]]
            for(j in 1:ku[t]) {
              BeU[,ku[t],ku[t+1],j,1:ku[t+1]] = BeU[,ku[t],ku[t+1],j,1:ku[t+1]] -BeU[,ku[t],ku[t+1],j,min(j,ku[t+1])]}
          }
          H = !any(rollapply(kc, 2, identical, c(ku[t-1],ku[t])))
          if(H==1){ # update step for coefficients modulating transition k_{t-1} to k_{t}
            BeU[,ku[t-1],ku[t],1:ku[t-1],1:ku[t]] = BeC[,kc[t-1],kc[t],1:kc[t-1],(1:kc[t])[-j0]]
            Bej0B = BeC[,kc[t-1],kc[t],1:kc[t-1],j0]
            for(j in 1:ku[t-1]){
              BeU[,ku[t-1],ku[t],j,1:ku[t]] = BeU[,ku[t-1],ku[t],j,1:ku[t]] - BeU[,ku[t-1],ku[t],j,min(j,ku[t])]}
          }
        }
        Be[iter,,,,,] <- BeU
        xi[iter,,,] <- xiU
        sigma[iter,,,] <- sigmaU
        k[iter,] = ku
        for(i in 1:n){ # update transition probs
          for(tt in 2:length(k[1,])){
            PIsU[i,(tt-1),k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]]= tmat(K = c(k[iter,tt-1],k[iter,tt]), Be[iter,,k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]],Z[i,tt,],10^-5)
          }
        }
        ilik = sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xiU[,,], sigma = sigmaU[,,], pi = pis[i,,], PI = PIsU[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
        liks[iter,t] = sum(ilik) # calculate likelihood
        if(bdt <= 0.5){
          if(!is.finite(liks[iter,t])){
            Abt = 0
          }else{
            # accept a birth
            Abt = (liks[iter,t] - liks[iter,t-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T)) - (sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
              log((1/(kc[t]-1))/(1/ku[t])) +     
              (log(ku[t]) - (sum(dnorm(xij0,mxi,vxi,log = T),na.rm = T)*L +(sum(dgamma(sij0,msi,vsi,log = T),na.rm = T)*L+(sum(dnorm(Bej0F,mBe,vBe,log = T),na.rm = T))*W + sum(dnorm(Bej0B,mBe,vBe,log = T),na.rm = T))*H))
            Abt = exp(Abt)}
          ab = runif(1)
          if(ab > min(1,Abt)){
            liks[iter,t] = liks[iter,t-1]
            ku = kc
            k[iter,] = kc
            xiU = xiC
            sigmaU = sigmaC
            BeU = BeC
            Be[iter,,,,,] = BeC
            xi[iter,,,]=xiC
            sigma[iter,,,]=sigmaC}else{PIs=PIsU}
        }else{
          if(!is.finite(liks[iter,t])){
            Abt = 0
          }else{
            # accept a death
            Abt = (liks[iter,t] - liks[iter,t-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T)) - (sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
              log((1/ku[t])/(1/(kc[t]-1))) +     
              (-(sum(dnorm(xij0,mxi,vxi,log = T),na.rm = T)*L +(sum(dgamma(sij0,msi,vsi,log = T),na.rm = T)*L+(sum(dnorm(Bej0F,mBe,vBe,log = T),na.rm = T))*W + sum(dnorm(Bej0B,mBe,vBe,log = T),na.rm = T))*H)+log(ku[t]))
            Abt = exp(Abt)}
          ab = runif(1)
          if(ab >min(1,Abt)){
            liks[iter,t] = liks[iter,t-1]
            ku = kc
            k[iter,] = kc
            xiU = xiC
            sigmaU = sigmaC
            BeU = BeC
            Be[iter,,,,,] = BeC
            xi[iter,,,]=xiC
            sigma[iter,,,]=sigmaC}else{PIs=PIsU}
        }
      }
      # t = T
      BeC = array(Be[iter,,,,,],c(dim(Be)[-1]))
      xiC = xi[iter,,,]
      sigmaC = sigma[iter,,,]
      xiU = sigmaU <- array(NA,c(kmax,kmax,r))
      BeU <- array(NA, c(n.cov,kmax,kmax,kmax,kmax))
      ku = kc <- k[iter,]
      xij0 = sij0 = Bej0B = NA
      if(kc[Ti]==1){
        bdT = 0
      }else if(kc[Ti]==kmax){
        bdT=1
      }else{bdT = runif(1)}
      if(bdT<=0.5){
        # BIRTH
        ku[Ti] = kc[Ti]+1
        ksub = ku[-Ti] 
        j0 = sample(c(1:ku[Ti]),1) # sample the location
        L = !any(ksub==ku[Ti])
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        if(L==1){ # update step for centroids and sds
          xiU[ku[Ti],(1:ku[Ti])[-j0],] = xiC[kc[Ti],1:kc[Ti],]
          xij0 = xiU[ku[Ti],j0,] = rnorm(r,mxi,vxi)
          sigmaU[ku[Ti],(1:ku[Ti])[-j0],] = sigmaC[kc[Ti],1:kc[Ti],]
          sij0 = sigmaU[ku[Ti],j0,] = rgamma(r,msi,vsi)
          o[ku[Ti],1:ku[Ti]] = order(xiU[ku[Ti],1:ku[Ti],1])
          xiU[ku[Ti],1:ku[Ti],1]=sort(xiU[ku[Ti],1:ku[Ti],1])
          for(h in 2:dim(y)[3]) {
            xiU[ku[Ti],1:ku[Ti],h] <- xiU[ku[Ti],o[ku[Ti],1:ku[Ti]],h]
          }
        }
        # check if the transition k_{t-1} to k_{t} needs update
        H = !any(rollapply(kc, 2, identical, c(ku[Ti-1],ku[Ti])))
        for(tt in 1:(Ti-2)){
          BeU[,ku[tt],ku[tt+1],1:ku[tt],1:ku[tt+1]] = BeC[,ku[tt],ku[tt+1],1:ku[tt],1:ku[tt+1]]
        }
        if(H==1){ # update step for transition k_{t-1} to k_{t}
          BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],(1:ku[Ti])[-j0]] = BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],1:kc[Ti]]
          Bej0B = BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],j0] = rnorm(n.cov*ku[Ti-1],mBe,vBe)
          for(j in 1:ku[Ti-1]) {
            BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] = BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] - BeU[,ku[Ti-1],ku[Ti],j,min(j,ku[Ti])]}
        }
      }else{
        # DEATH
        xij0 = sij0 = Bej0B = NA
        ku[Ti] = kc[Ti]-1
        ksub = ku[-Ti]
        j0 = sample(c(1:kc[Ti]),1) # sample the location
        L = !any(ksub==ku[Ti])
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        if(L==1){ # update step for centroids and sds
          xij0 = xiC[kc[Ti],j0,]
          sij0 = sigmaC[kc[Ti],j0,]
          xiU[ku[Ti],1:ku[Ti],] <- xiC[kc[Ti],(1:kc[Ti])[-j0],]
          sigmaU[ku[Ti],1:ku[Ti],] <- sigmaC[kc[Ti],(1:kc[Ti])[-j0],]
          o[ku[Ti],1:ku[Ti]] = order(xiU[ku[Ti],1:ku[Ti],1])
          xiU[ku[Ti],1:ku[Ti],1]=sort(xiU[ku[Ti],1:ku[Ti],1])
          for(h in 2:dim(y)[3]) {
            xiU[ku[Ti],1:ku[Ti],h] <- xiU[ku[Ti],o[ku[Ti],1:ku[Ti]],h]
          }
        }
        # check if transition k_{t-1} to k_{t} needs update
        H = !any(rollapply(kc, 2, identical, c(ku[Ti-1],ku[Ti])))
        for(tt in 1:(Ti-2)){
          BeU[,ku[tt],ku[tt+1],1:ku[tt],1:ku[tt+1]] = BeC[,ku[tt],ku[tt+1],1:ku[tt],1:ku[tt+1]]
        }
        if(H==1){ # update step for transition k_{t-1} to k_{t}
          BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],1:ku[Ti]] = BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],(1:kc[Ti])[-j0]]
          Bej0B = BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],j0]
          for(j in 1:ku[Ti-1]){
            BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] = BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] - BeU[,ku[Ti-1],ku[Ti],j,min(j,ku[Ti])]}
        }
      }
      Be[iter,,,,,] <- BeU
      xi[iter,,,] <- xiU
      sigma[iter,,,] <- sigmaU
      k[iter,] = ku
      for(i in 1:n){ # update transition probs
        for(tt in 2:length(k[1,])){
          PIsU[i,(tt-1),k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]]= tmat(K = c(k[iter,tt-1],k[iter,tt]), Be[iter,,k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]],Z[i,tt,],10^-5)
        }
      }
      ilik = sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xiU[,,], sigma = sigmaU[,,], pi = pis[i,,], PI = PIsU[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
      liks[iter,Ti] = sum(ilik) # calculate likelihood
      if(bdT <= 0.5){
        if(!is.finite(liks[iter,Ti])){
          AbT = 0
        }else{
          # accept a birth
          AbT = (liks[iter,Ti] - liks[iter,Ti-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T)) 
                                                   - (sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T))) +log((1/(kc[Ti]-1))/(1/(ku[Ti]))) +
            (log(ku[Ti]) - (sum(dnorm(xij0,mxi,vxi,log = T),na.rm = T)*L + sum(dgamma(sij0,msi,vsi,log = T),na.rm = T)*L + sum(dnorm(Bej0B,mBe,vBe,log = T),na.rm = T)*H))
          AbT = exp(AbT)}
        ab = runif(1)
        if(ab > min(1,AbT)){
          liks[iter,Ti] = liks[iter,Ti-1]
          ku = kc
          k[iter,] = kc
          xiU = xiC
          sigmaU = sigmaC
          BeU = BeC
          Be[iter,,,,,] = BeC
          xi[iter,,,]=xiC
          sigma[iter,,,]=sigmaC}else{PIs=PIsU}
      }else{
        if(!is.finite(liks[iter,Ti])){
          AbT = 0
        }else{
          # accept a death
          AbT = (liks[iter,Ti] - liks[iter,Ti-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T)) 
                                                   - (sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T))) +log((1/(ku[Ti]))/(1/(kc[Ti]-1))) +
            (-(sum(dnorm(xij0,mxi,vxi,log = T),na.rm = T)*L + sum(dgamma(sij0,msi,vsi,log = T),na.rm = T)*L + sum(dnorm(Bej0B,mBe,vBe,log = T),na.rm = T)*H)+log(ku[Ti]) )
          AbT = exp(AbT)}
        ab = runif(1)
        if(ab >min(1,AbT)){
          liks[iter,Ti] = liks[iter,Ti-1]
          ku = kc
          k[iter,] = kc
          xiU = xiC
          sigmaU = sigmaC
          BeU = BeC
          Be[iter,,,,,] = BeC
          xi[iter,,,]=xiC
          sigma[iter,,,]=sigmaC}else{PIs=PIsU}
      }
      
      #####################
      ### SPLIT/COMBINE ###
      #####################
      
      if(verbose){print("## split/combine")}
      ku = kc = k[iter,]
      xiC <- xi[iter,,,]
      sigmaC <- sigma[iter,,,]
      BeC <- array(Be[iter,,,,,],c(dim(Be)[-1]))
      beC = array(be[iter,,,],c(dim(be)[-1]))
      xiU = sigmaU <- array(NA,c(kmax,kmax,r))
      BeU <- array(NA, c(n.cov,kmax,kmax,kmax,kmax))
      # t = 1 # 
      u = w = e1 = e12 = NA
      
      if(kc[1]==1){ 
        sc1 = 0
      }else if(kc[1]==kmax){
        sc1=1
      }else{sc1 = runif(1)}
      
      if(sc1 <= 0.5){
        # SPLIT
        ku[1] = kc[1]+1
        j0=sample(c(1:kc[1]),1) # sample j_0
        j1=j0+1
        beU[,ku[1],(1:ku[1])[-c(j0,j1)]] = beC[,kc[1],(1:kc[1])[-j0]]
        e1 = rnorm(n.cov,mepsb,vepsb)
        beU[,ku[1],j0] = beC[,kc[1],j0] - e1 # split with auxiliary e1
        beU[,ku[1],j1] = beC[,kc[1],j0] + e1
        ksub=kc[-1]
        L = !any(ksub==ku[1]) # check if centroids and sds need split
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        if(L==1){ # split step for centroids and sds
          xiU[ku[1],(1:ku[1])[-c(j0,j1)],] = xiC[kc[1],(1:kc[1])[-j0],]
          sigmaU[ku[1],(1:ku[1])[-c(j0,j1)],] = sigmaC[kc[1],(1:kc[1])[-j0],]
          u = rnorm(r,mu,vu)
          w = rgamma(r,mw,vw)
          xiU[ku[1],j0,] = xiC[kc[1],j0,] - u*sigmaC[kc[1],j0,] # split centroids
          xiU[ku[1],j1,] = xiC[kc[1],j0,] + u*sigmaC[kc[1],j0,] # with auxiliary u
          sigmaU[ku[1],j0,] = sigmaC[kc[1],j0,]*w # split sds
          sigmaU[ku[1],j1,] = sigmaC[kc[1],j0,]/w # with auxiliary w
          o[ku[1],1:ku[1]] = order(xiU[ku[1],1:ku[1],1])
          xiU[ku[1],1:ku[1],1]=sort(xiU[ku[1],1:ku[1],1])
          for(h in 2:dim(y)[3]) {
            xiU[ku[1],1:ku[1],h] <- xiU[ku[1],o[ku[1],1:ku[1]],h]
          }
        }
        # check if transition k_{1} to k_{2} needs update
        W = !any(rollapply(kc, 2, identical, c(ku[1],ku[2])))
        for(tt in 1:(Ti-2)){
          BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
        }
        if(W==1){ # split step for coefficients modulating transition k_{1} to k_{2}
          BeU[,ku[1],ku[2],(1:ku[1])[-c(j0,j1)],1:ku[2]] = BeC[,kc[1],kc[2],(1:kc[1])[-j0],1:kc[2]]
          e12 = rnorm(n.cov*kc[2],mepsB,vepsB)
          BeU[,ku[1],ku[2],j0,1:ku[2]] = BeC[,kc[1],kc[2],j0,1:kc[2]] - e12 # split with
          BeU[,ku[1],ku[2],j1,1:ku[2]] = BeC[,kc[1],kc[2],j0,1:kc[2]] + e12 # auxiliary e12
          for(j in 1:ku[1]){
            BeU[,ku[1],ku[2],j,1:ku[2]] = BeU[,ku[1],ku[2],j,1:ku[2]] - BeU[,ku[1],ku[2],j,min(j,ku[2])]}
        }
      }else{
        # COMBINE
        ku[1] = kc[1]-1
        j0=sample(c((1:(kc[1]))[-kc[1]]),1) 
        j1 = j0+1 # sample a pair of adjacent states to merge
        beU[,ku[1],(1:ku[1])[-j0]] = beC[,kc[1],(1:kc[1])[-c(j0,j1)]]
        beU[,ku[1],j0] = (beC[,kc[1],j0] + beC[,kc[1],j1])/2 # combine step for coefficients modulating initial probs
        e1 = (beC[,kc[1],j0] - beC[,kc[1],j1])/2            # auxialiry variable
        ksub=kc[-1]
        L = !any(ksub==ku[1]) # check if centroids and sds need combine
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        if(L==1){
          xiU[ku[1],(1:ku[1])[-j0],] = xiC[kc[1],(1:kc[1])[-c(j0,j1)],]
          sigmaU[ku[1],(1:ku[1])[-j0],] = sigmaC[kc[1],(1:kc[1])[-c(j0,j1)],]
          xiU[ku[1],j0,] = (xiC[kc[1],j0,] + xiC[kc[1],j1,])/2 # combine centroids
          u = (xiC[kc[1],j0,] - xiC[kc[1],j1,])/2
          sigmaU[ku[1],j0,] = sqrt(sigmaC[kc[1],j0,] * sigmaC[kc[1],j1,]) # combine sds
          w = sqrt(sigmaC[kc[1],j0,]/sigmaC[kc[1],j1,])
          o[ku[1],1:ku[1]] = order(xiU[ku[1],1:ku[1],1])
          xiU[ku[1],1:ku[1],1]=sort(xiU[ku[1],1:ku[1],1])
          for(h in 2:dim(y)[3]) {
            xiU[ku[1],1:ku[1],h] <- xiU[ku[1],o[ku[1],1:ku[1]],h]
          }
        }
        # check if transition k_{1} to k_{2} needs merging
        W = !any(rollapply(kc, 2, identical, c(ku[1],ku[2])))
        for(tt in 1:(Ti-2)){
          BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
        }
        if(W==1){ # combine step for coefficients modulating transition k_{1} to k_{2}
          BeU[,ku[1],ku[2],(1:ku[1])[-j0],1:ku[2]] = BeC[,kc[1],kc[2],(1:kc[1])[-c(j0,j1)],1:kc[2]]
          BeU[,ku[1],ku[2],j0,1:ku[2]] = (BeC[,kc[1],kc[2],j0,1:kc[2]] + BeC[,kc[1],kc[2],j1,1:kc[2]])/2
          e12 = (BeC[,kc[1],kc[2],j0,1:kc[2]] - BeC[,kc[1],kc[2],j1,1:kc[2]])/2
          for(j in 1:ku[1]){
            BeU[,ku[1],ku[2],j,1:ku[2]] = BeU[,ku[1],ku[2],j,1:ku[2]] - BeU[,ku[1],ku[2],j,min(j,ku[2])]}
        }
      }
      Be[iter,,,,,] <- BeU
      xi[iter,,,] <- xiU
      sigma[iter,,,] <- sigmaU
      k[iter,] = ku
      for(i in 1:n){ # update initial and transition probs
        pisU[i,k[iter,1],1:k[iter,1]] <- inprob(k[iter,1], be=beU[,k[iter,1],1:k[iter,1]], Z[i,1,], thB = 10^-5)
        for(tt in 2:length(k[1,])){
          PIsU[i,(tt-1),k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]]= tmat(K = c(k[iter,tt-1],k[iter,tt]), Be[iter,,k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]],Z[i,tt,],10^-5)
        }
      }
      ilik = sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xiU[,,], sigma = sigmaU[,,], pi = pisU[i,,], PI = PIsU[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
      liks[iter,(Ti+1)] = sum(ilik) # calculate likelihood
      if(sc1 <= 0.5){
        if(!is.finite(liks[iter,(Ti+1)])){
          As1 = 0
        }else{
          # accept a split
          J = sum((log(4^r * prod((sigmaC[kc[1],j0,]^2)/w))^L)*L,log(2*n.cov),log(2*kc[2]*n.cov)*W)
          As1 = (liks[iter,(Ti+1)]-liks[iter,Ti])+((sum(dnorm(beU,mbe,vbe,log = T),na.rm = T) + sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T))-(sum(dnorm(beC,mbe,vbe,log = T),na.rm = T)+sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
            log((1/(kc[1]-1))/(1/ku[1])) + log(ku[1]) + (J - (sum(dnorm(e1,mepsb,vepsb,log = T)) + sum(dnorm(e12,mepsB,vepsB,log = T),na.rm = T)*W + sum(dnorm(u,mu,vu,log = T),na.rm = T)*L + sum(dgamma(w,mw,vw,log = T),na.rm = T)*L))
          As1=exp(As1)}
        as = runif(1)
        if(as > min(1,As1)){
          liks[iter,(Ti+1)] = liks[iter,Ti]
          ku = kc
          k[iter,] = kc
          xiU = xiC
          xi[iter,,,] = xiC
          sigmaU = sigmaC
          Be[iter,,,,,] <- BeC
          sigma[iter,,,] = sigmaC
          be[iter,,,] = beC
          BeU = BeC}else{
            pis = pisU
            be[iter,,,] = beU
            PIs = PIsU
          }
      }else{
        if(!is.finite(liks[iter,(Ti+1)])){
          As1 = 0
        }else{
          # accept a combine
          J = sum((log(4^r * prod((sigmaC[kc[1],j0,]^2)/w))^L)*L,log(2*n.cov),log(2*kc[2]*n.cov)*W)
          As1 = (liks[iter,(Ti+1)]-liks[iter,Ti])+((sum(dnorm(beU,mbe,vbe,log = T),na.rm = T) + sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T))-(sum(dnorm(beC,mbe,vbe,log = T),na.rm = T)+sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
            log((1/ku[1])/(1/(kc[1]-1))) - log(ku[1]) + (-J + (sum(dnorm(e1,mepsb,vepsb,log = T)) + sum(dnorm(e12,mepsB,vepsB,log = T),na.rm = T)*W + sum(dnorm(u,mu,vu,log = T),na.rm = T)*L + sum(dgamma(w,mw,vw,log = T),na.rm = T)*L))
          As1 = exp(As1)}
        as = runif(1)
        if(as > min(1,As1)){
          liks[iter,(Ti+1)] = liks[iter,Ti]
          ku = kc
          k[iter,] = kc
          xiU = xiC
          xi[iter,,,] = xiC
          sigmaU = sigmaC
          Be[iter,,,,,] <- BeC
          sigma[iter,,,] = sigmaC
          be[iter,,,] = beC
          BeU = BeC}else{
            pis = pisU
            be[iter,,,] = beU
            PIs = PIsU
          }   
      }
      # 1 < t < T
      for(t in 2:(Ti-1)){
        ev = eb = u = w = NA
        BeC = array(Be[iter,,,,,],c(dim(Be)[-1]))
        xiC = xi[iter,,,]
        sigmaC = sigma[iter,,,]
        xiU = sigmaU <- array(NA,c(kmax,kmax,r))
        BeU <- array(NA, c(n.cov,kmax,kmax,kmax,kmax))
        ku = kc <- k[iter,]
        
        if(kc[t]==1){
          sct = 0
        }else if(kc[t]==kmax){
          sct=1
        }else{sct = runif(1)}
        if(sct <= 0.5){
          # SPLIT
          ksub=ku[-t]  
          j0 = sample(c(1:kc[t]),1) # sample the state j_0
          j1 = j0+1
          ku[t]=kc[t]+1
          for(t1 in unique(ku)){
            xiU[t1,1:t1,] = xiC[t1,1:t1,]
            sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
          }
          L = !any(ksub==ku[t]) # check if centroids and sds need splitting
          if(L==1){
            xiU[ku[t],(1:ku[t])[-c(j0,j1)],] = xiC[kc[t],(1:kc[t])[-j0],]
            sigmaU[ku[t],(1:ku[t])[-c(j0,j1)],] = sigmaC[kc[t],(1:kc[t])[-j0],]
            u = rnorm(r,mu,vu)
            w = dgamma(r,mw,vw)
            xiU[ku[t],j0,] = xiC[kc[t],j0,] - u*sigmaC[kc[t],j0,] # split centroids
            xiU[ku[t],j1,] = xiC[kc[t],j0,] + u*sigmaC[kc[t],j0,] # with auxiliary u
            sigmaU[ku[t],j0,] = sigmaC[kc[t],j0,]*w # split sds
            sigmaU[ku[t],j1,] = sigmaC[kc[t],j0,]/w # with auxiliary w
            o[ku[t],1:ku[t]] = order(xiU[ku[t],1:ku[t],1])
            xiU[ku[t],1:ku[t],1]=sort(xiU[ku[t],1:ku[t],1])
            for(h in 2:dim(y)[3]) {
              xiU[ku[t],1:ku[t],h] <- xiU[ku[t],o[ku[t],1:ku[t]],h]
            }
          }
          # check if transition k_{t} to k_{t+1} needs splitting
          W = !any(rollapply(kc, 2, identical, c(ku[t],ku[t+1])))
          for(tt in 1:(Ti-2)){
            BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
          }
          if(W==1){ # split step for coefficients modulating transition k_{t} to k_{t+1}
            BeU[,ku[t],ku[t+1],(1:ku[t])[-c(j0,j1)],1:ku[t+1]] = BeC[,kc[t],kc[t+1],(1:kc[t])[-j0],1:kc[t+1]]
            eb = rnorm(n.cov*kc[t+1])
            BeU[,ku[t],ku[t+1],j0,1:ku[t+1]] = BeC[,kc[t],kc[t+1],j0,1:kc[t+1]] - eb
            BeU[,ku[t],ku[t+1],j1,1:ku[t+1]] = BeC[,kc[t],kc[t+1],j0,1:kc[t+1]] + eb
            for(j in 1:ku[t]){
              BeU[,ku[t],ku[t+1],j,1:ku[t+1]] = BeU[,ku[t],ku[t+1],j,1:ku[t+1]] - BeU[,ku[t],ku[t+1],j,min(j,ku[t+1])]}
          }
          # check if transition k_{t-1} to k_{t} needs updating
          H = !any(rollapply(kc, 2, identical, c(ku[t-1],ku[t])))
          if(H==1){
            BeU[,ku[t-1],ku[t],1:ku[t-1],(1:ku[t])[-c(j0,j1)]] = BeC[,kc[t-1],kc[t],1:kc[t-1],(1:kc[t])[-j0]]
            ev = rnorm(n.cov*kc[t-1],mepsB,vepsB)
            BeU[,ku[t-1],ku[t],1:ku[t-1],j0] = BeC[,kc[t-1],kc[t],1:kc[t-1],j0] - ev # split step for coefficients
            BeU[,ku[t-1],ku[t],1:ku[t-1],j1] = BeC[,kc[t-1],kc[t],1:kc[t-1],j0] + ev # with auxiliary ev
            for(j in 1:ku[t-1]){
              BeU[,ku[t-1],ku[t],j,1:ku[t]] = BeU[,ku[t-1],ku[t],j,1:ku[t]] - BeU[,ku[t-1],ku[t],j,min(j,ku[t])]}
          }
        }else{
          # COMBINE
          ku[t] = kc[t]-1
          j0 = sample(c((1:kc[t])[-kc[t]]),1)
          j1 = j0+1 # sample a pair of adjacent states to merge
          ksub = ku[-t] 
          L = !any(ksub==ku[t]) # check if centroids and sds need combining
          for(t1 in unique(ku)){
            xiU[t1,1:t1,] = xiC[t1,1:t1,]
            sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
          }
          if(L==1){
            xiU[ku[t],(1:ku[t])[-j0],] = xiC[kc[t],(1:kc[t])[-c(j0,j1)],]
            sigmaU[ku[t],(1:ku[t])[-j0],] = sigmaC[kc[t],(1:kc[t])[-c(j0,j1)],]
            xiU[ku[t],j0,] = (xiC[kc[t],j0,] + xiC[kc[t],j1,])/2 # combine centroids
            sigmaU[ku[t],j0,] = sqrt(sigmaC[kc[t],j0,]*sigmaC[kc[t],j1,])
            u = (xiC[kc[t],j0,] - xiC[kc[t],j1,])/2
            w = sqrt(sigmaC[kc[t],j0,]/sigmaC[kc[t],j1,]) # combine sds
            o[ku[t],1:ku[t]] = order(xiU[ku[t],1:ku[t],1])
            xiU[ku[t],1:ku[t],1]=sort(xiU[ku[t],1:ku[t],1])
            for(h in 2:dim(y)[3]) {
              xiU[ku[t],1:ku[t],h] <- xiU[ku[t],o[ku[t],1:ku[t]],h]
            }
          }
          # check if transition k_{t} to k_{t+1} needs updating
          W = !any(rollapply(kc, 2, identical, c(ku[t],ku[t+1])))
          for(tt in 1:(Ti-2)){
            BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
          }
          if(W==1){ # combine coefficients modulating transition k_{t} to k_{t+1}
            BeU[,ku[t],ku[t+1],(1:ku[t])[-j0],1:ku[t+1]] = BeC[,kc[t],kc[t+1],(1:kc[t])[-c(j0,j1)],1:kc[t+1]]
            BeU[,ku[t],ku[t+1],j0,1:ku[t+1]] = (BeC[,kc[t],kc[t+1],j0,1:kc[t+1]] + BeC[,kc[t],kc[t+1],j1,1:kc[t+1]])/2
            eb = (BeC[,kc[t],kc[t+1],j0,1:kc[t+1]] - BeC[,kc[t],kc[t+1],j1,1:kc[t+1]])/2
            for(j in 1:ku[t]){
              BeU[,ku[t],ku[t+1],j,1:ku[t+1]] = BeU[,ku[t],ku[t+1],j,1:ku[t+1]] - BeU[,ku[t],ku[t+1],j,min(j,ku[t+1])]}
          }
          # check if transition k_{t-1} to k_{t} needs updating
          H = !any(rollapply(kc, 2, identical, c(ku[t-1],ku[t])))
          if(H==1){ # combine coefficients modulating transition k_{t-1} to k_{t}
            BeU[,ku[t-1],ku[t],1:ku[t-1],(1:ku[t])[-j0]] = BeC[,kc[t-1],kc[t],1:kc[t-1],(1:kc[t])[-c(j0,j1)]]
            BeU[,ku[t-1],ku[t],1:ku[t-1],j0] = (BeC[,kc[t-1],kc[t],1:kc[t-1],j0] + BeC[,kc[t-1],kc[t],1:kc[t-1],j1])/2
            ev = (BeC[,kc[t-1],kc[t],1:kc[t-1],j0] - BeC[,kc[t-1],kc[t],1:kc[t-1],j1])/2
            for(j in 1:ku[t-1]){
              BeU[,ku[t-1],ku[t],j,1:ku[t]] = BeU[,ku[t-1],ku[t],j,1:ku[t]] - BeU[,ku[t-1],ku[t],j,min(j,ku[t])]}
          }
        }
        Be[iter,,,,,] <- BeU
        xi[iter,,,] <- xiU
        sigma[iter,,,] <- sigmaU
        k[iter,] = ku
        for(i in 1:n){ # update transition probs
          for(tt in 2:length(k[1,])){
            PIsU[i,(tt-1),k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]]= tmat(K = c(k[iter,tt-1],k[iter,tt]), Be[iter,,k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]],Z[i,tt,],10^-5)
          }
        }
        ilik = sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xiU[,,], sigma = sigmaU[,,], pi = pis[i,,], PI = PIsU[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
        liks[iter,(t+Ti)] = sum(ilik) # calculate likelihood
        if(sct <= 0.5){
          if(!is.finite(liks[iter,(t+Ti)])){
            Ast = 0
          }else{
            # accept a split
            J = sum((log(4^r * prod((sigmaC[kc[t],j0,]^2)/w))^L)*L,log(2*kc[t-1]*n.cov)*H,log(2*kc[t+1]*n.cov)*W)
            Ast = (liks[iter,(t+Ti)]-liks[iter,(t+Ti)-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T))-(sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
              log((1/(kc[t]-1))/(1/ku[t])) + log(ku[t]) + (J - (sum(dnorm(eb,mepsB,vepsB,log = T),na.rm = T)*W + sum(dnorm(ev,mepsB,vepsB,log = T),na.rm = T)*H + sum(dnorm(u,mu,vu,log = T),na.rm = T)*L + sum(dgamma(w,mw,vw,log = T),na.rm = T)*L))
            Ast = exp(Ast)}
          as = runif(1)
          if(as > min(1,Ast)){
            liks[iter,(t+Ti)] = liks[iter,(t+Ti)-1]
            ku = kc
            k[iter,] = kc
            xiU = xiC
            xi[iter,,,] = xiC
            sigmaU = sigmaC
            Be[iter,,,,,] <- BeC
            sigma[iter,,,] = sigmaC
            BeU = BeC}else{
              PIs = PIsU}
        }else{
          if(!is.finite(liks[iter,(t+Ti)])){
            Ast = 0
          }else{
            # accept a combine
            J = sum((log(4^r * prod((sigmaC[kc[t],j0,]^2)/w))^L)*L,log(2*kc[t-1]*n.cov)*H,log(2*kc[t+1]*n.cov)*W)
            Ast = (liks[iter,(t+Ti)]-liks[iter,(t+Ti)-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T))-(sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
              log((1/(kc[t]-1))/(1/ku[t])) + log(ku[t]) + (-J + (sum(dnorm(eb,mepsB,vepsB,log = T),na.rm = T)*W + sum(dnorm(ev,mepsB,vepsB,log = T),na.rm = T)*H + sum(dnorm(u,mu,vu,log = T),na.rm = T)*L + sum(dgamma(w,mw,vw,log = T),na.rm = T)*L))
            Ast = exp(Ast)}
          as = runif(1)
          if(as > min(1,Ast)){
            liks[iter,(t+Ti)] = liks[iter,(t+Ti)-1]
            ku = kc
            k[iter,] = kc
            xiU = xiC
            xi[iter,,,] = xiC
            sigmaU = sigmaC
            Be[iter,,,,,] <- BeC
            sigma[iter,,,] = sigmaC
            BeU = BeC}else{PIs=PIsU}
        }
      }
      # t = T
      u = w = eT = NA
      BeC = array(Be[iter,,,,,],c(dim(Be)[-1]))
      xiC = xi[iter,,,]
      sigmaC = sigma[iter,,,]
      xiU = sigmaU <- array(NA,c(kmax,kmax,r))
      BeU <- array(NA, c(n.cov,kmax,kmax,kmax,kmax))
      ku = kc <- k[iter,]
      labt = c(1:kc[Ti])
      j0 = sample(labt[-kc[Ti]], 1) # sample a state j_0
      j1 = j0+1
      if(kc[Ti]==1){
        scT = 0
      }else if(kc[Ti]==kmax){
        scT=1
      }else{scT = runif(1)}
      if(scT <= 0.5){
        # SPLIT  
        ku[Ti] = kc[Ti] + 1  
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        ksub = ku[-Ti]
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        L = !any(ksub==ku[Ti]) # check if centroids and sds need update
        if(L==1){
          xiU[ku[Ti],(1:ku[Ti])[-c(j0,j1)],] = xiC[kc[Ti],(1:kc[Ti])[-j0],]  
          sigmaU[ku[Ti],(1:ku[Ti])[-c(j0,j1)],] = sigmaC[kc[Ti],(1:kc[Ti])[-j0],]
          u = rnorm(r,mu,vu)
          w = rgamma(r,mw,vw)
          xiU[ku[Ti],j0,] = xiC[kc[Ti],j0,] - u*sigmaC[kc[Ti],j0,] # split centroids
          xiU[ku[Ti],j1,] = xiC[kc[Ti],j0,] + u*sigmaC[kc[Ti],j0,] # with auxiliary u
          sigmaU[ku[Ti],j0,] = sigmaC[kc[Ti],j0,]*w # split sds
          sigmaU[ku[Ti],j1,] = sigmaC[kc[Ti],j0,]/w # with auxiliary w
          o[ku[Ti],1:ku[Ti]] = order(xiU[ku[Ti],1:ku[Ti],1])
          xiU[ku[Ti],1:ku[Ti],1]=sort(xiU[ku[Ti],1:ku[Ti],1])
          for(h in 2:dim(y)[3]) {
            xiU[ku[Ti],1:ku[Ti],h] <- xiU[ku[Ti],o[ku[Ti],1:ku[Ti]],h]
          }
        }
        for(tt in 1:(Ti-2)){
          BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
        }
        # check if transition k_{Ti-1} to k_{Ti} needs update
        H = !any(rollapply(kc, 2, identical, c(ku[Ti-1],ku[Ti])))
        if(H==1){ # split coefficients modulating transition k_{Ti-1} to k_{Ti}
          BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],(1:ku[Ti])[-c(j0,j1)]] = BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],(1:kc[Ti])[-j0]]
          eT = rnorm(n.cov*kc[Ti-1],mepsB,vepsB)
          BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],j0] = BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],j0] - eT # auxiliary eT
          BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],j1] = BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],j0] + eT
          for(j in 1:ku[Ti-1]){
            BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] = BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] - BeU[,ku[Ti-1],ku[Ti],j,min(j,ku[Ti])]
          }
        }
      }else{
        # COMBINE
        ku[Ti] = kc[Ti]-1
        ksub = ku[-Ti]
        for(t1 in unique(ku)){
          xiU[t1,1:t1,] = xiC[t1,1:t1,]
          sigmaU[t1,1:t1,] = sigmaC[t1,1:t1,]
        }
        L = !any(ksub==ku[Ti]) # check if centroids and sds need update
        if(L==1){
          xiU[ku[Ti],(1:ku[Ti])[-j0],] = xiC[kc[Ti],(1:kc[Ti])[-c(j0,j1)],]
          sigmaU[ku[Ti],(1:ku[Ti])[-j0],] = sigmaC[kc[Ti],(1:kc[Ti])[-c(j0,j1)],]
          xiU[ku[Ti],j0,] = (xiC[kc[Ti],j0,] + xiC[kc[Ti],j1,])/2 # combine centroids
          u = (xiC[kc[Ti],j0,] - xiC[kc[Ti],j1,])/2
          sigmaU[ku[Ti],j0,] = sqrt(sigmaC[kc[Ti],j0,]*sigmaC[kc[Ti],j1,]) # combine sds
          w = sqrt(sigmaC[kc[Ti],j0,]/sigmaC[kc[Ti],j1,])
        }
        o[ku[Ti],1:ku[Ti]] = order(xiU[ku[Ti],1:ku[Ti],1])
        xiU[ku[Ti],1:ku[Ti],1]=sort(xiU[ku[Ti],1:ku[Ti],1])
        for(h in 2:dim(y)[3]) {
          xiU[ku[Ti],1:ku[Ti],h] <- xiU[ku[Ti],o[ku[Ti],1:ku[Ti]],h]
        }
        for(tt in 1:(Ti-2)){
          BeU[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]] = BeC[,ksub[tt],ksub[tt+1],1:ksub[tt],1:ksub[tt+1]]
        }
        # check if transition k_{Ti-1} to k_{Ti} needs updating
        H = !any(rollapply(kc, 2, identical, c(ku[Ti-1],ku[Ti])))
        if(H==1){ # combine step for coefficients modulating transition k_{Ti-1} to k_{Ti}
          BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],(1:ku[Ti])[-j0]] = BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],(1:kc[Ti])[-c(j0,j1)]]
          BeU[,ku[Ti-1],ku[Ti],1:ku[Ti-1],j0] = (BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],j0] + BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],j1])/2
          eT = (BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],j0] - BeC[,kc[Ti-1],kc[Ti],1:kc[Ti-1],j1])/2
          for(j in 1:ku[Ti-1]){
            BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] = BeU[,ku[Ti-1],ku[Ti],j,1:ku[Ti]] - BeU[,ku[Ti-1],ku[Ti],j,min(j,ku[Ti])]
          }
        }
      }
      Be[iter,,,,,] <- BeU
      xi[iter,,,] <- xiU
      sigma[iter,,,] <- sigmaU
      k[iter,] = ku
      for(i in 1:n){ # update transition probs
        for(tt in 2:length(k[1,])){
          PIsU[i,(tt-1),k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]]= tmat(K = c(k[iter,tt-1],k[iter,tt]), Be[iter,,k[iter,tt-1],k[iter,tt],1:k[iter,tt-1],1:k[iter,tt]],Z[i,tt,],10^-5)
        }
      }
      ilik=sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xiU[,,], sigma = sigmaU[,,], pi = pis[i,,], PI = PIsU[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
      liks[iter,2*Ti] = sum(ilik) # calculate likelihood
      if(scT <= 0.5){
        if(!is.finite(liks[iter,(2*Ti)])){
          AsT = 0
        }else{
          # accept a split
          J = sum((log(4^r * prod((sigmaC[kc[Ti],j0,]^2)/w))^L)*L,log(2*kc[Ti-1]*n.cov)*H)
          AsT = (liks[iter,(2*Ti)]-liks[iter,(2*Ti)-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T))-(sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
            log((1/(kc[Ti]-1))/(1/ku[Ti])) + log(ku[Ti]) + (J - (sum(dnorm(ev,mepsB,vepsB,log = T),na.rm = T)*H + sum(dnorm(u,mu,vu,log = T),na.rm = T)*L + sum(dgamma(w,mw,vw,log = T),na.rm = T)*L))
          AsT = exp(AsT)}
        as = runif(1)
        if(as > min(1,AsT)){
          liks[iter,2*Ti] = liks[iter,(2*Ti)-1]
          ku <- kc
          k[iter,] <- kc
          xiU = xiC
          xi[iter,,,] = xiC
          sigmaU = sigmaC
          Be[iter,,,,,] <- BeC
          sigma[iter,,,] = sigmaC
          BeU = BeC  
        }else{PIs=PIsU}
      }else{
        if(!is.finite(liks[iter,(2*Ti)])){
          AsT = 0
        }else{
          # accept a combine
          J = sum((log(4^r * prod((sigmaC[kc[Ti],j0,]^2)/w))^L)*L,log(2*kc[Ti-1]*n.cov)*H)
          AsT = (liks[iter,(2*Ti)]-liks[iter,(2*Ti)-1])+((sum(dnorm(xiU,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeU,mBe,vBe,log = T),na.rm = T))-(sum(dnorm(xiC,mxi,vxi,log = T),na.rm = T) + sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T) + sum(dnorm(BeC,mBe,vBe,log = T),na.rm = T)))+
            log((1/(kc[Ti]-1))/(1/ku[Ti])) + log(ku[Ti]) + (-J + (sum(dnorm(ev,mepsB,vepsB,log = T),na.rm = T)*H + sum(dnorm(u,mu,vu,log = T),na.rm = T)*L + sum(dgamma(w,mw,vw,log = T),na.rm = T)*L))
          AsT = exp(AsT)}
        as = runif(1)
        if(as > min(1,AsT)){
          liks[iter,(2*Ti)] = liks[iter,(2*Ti-1)]
          ku = kc
          k[iter,] = kc
          xiU = xiC
          xi[iter,,,] = xiC
          sigmaU = sigmaC
          Be[iter,,,,,] <- BeC
          sigma[iter,,,] = sigmaC
          BeU = BeC  
        }else{PIs=PIsU}
      }
      
    }
    ################
    ### MH STEPS ###  
    ################
    
    # update xi
    if(verbose){print("## update xi")}
    xiU <- xiC <- xi[iter,,,]
    for(tt in unique(k[iter,])){
      xiU[tt,1:tt,1] <- xiC[tt,1:tt,1] + rnorm(length(c(xiC[tt,1:tt,1])),0,Tauxi[1])
      for(h in 2:(dim(y)[3])){
        xiU[tt,1:tt,h] <- xiC[tt,1:tt,h] + rnorm(length(c(xiC[tt,1:tt,h])),0,Tauxi[h])
      }
    }
    for(tt in unique(k[iter,])){
      o[tt,1:tt]<-order(xiU[tt,1:tt,1])
      xiU[tt,1:tt,1] = sort(xiU[tt,1:tt,1])
      for(h in 2:dim(y)[3]) {
        xiU[tt,1:tt,h] <- xiU[tt,o[tt,1:tt],h]
      }
    }
    xi[iter,,,] <- xiU
    
    ilik <- sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xiU, sigma = sigma[iter,,,], pi = pis[i,,], PI = PIs[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
    liks[iter,(2*Ti)+1] = sum(ilik)
    if(is.finite(liks[iter,(2*Ti)+1])){
      D = sum(dnorm(xiU,mxi,vxi,log = T), na.rm = T) - sum(dnorm(xiC,mxi,vxi,log = T), na.rm = T)
      Ax = (liks[iter,(2*Ti)+1] - liks[iter,(2*Ti)])+D
      Ax = exp(Ax)}
    if(!is.finite(liks[iter,2*Ti])){Ax = 0}
    axi = runif(1)
    if(axi > Ax){
      xi[iter,,,] <- xiC
      liks[iter,(2*Ti)+1] <- liks[iter,(2*Ti)]
    }
    
    # update sigma
    if(verbose){print("## update sigma")}
    sigmaC <- sigmaU <- sigma[iter,,,]
    for(tt in unique(k[iter,])){
      sigmaU[tt,1:tt,1] <- exp(log(sigmaC[tt,1:tt,1]) + rnorm(length(c(sigmaC[tt,1:tt,1])),0,0.05))
      for(h in 2:(dim(y)[3])){
        sigmaU[tt,1:tt,h] <- exp(log(sigmaC[tt,1:tt,h]) + rnorm(length(c(sigmaC[tt,1:tt,h])),0,0.05))
      }
    }
    sigma[iter,,,] <- sigmaU
    ilik <- sapply(1:n, function(i) likco(y = y[i,,], w=1, xi = xi[iter,,,], sigma = sigma[iter,,,], pi = pis[i,,], PI = PIs[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
    liks[iter,(2*Ti)+2] = sum(ilik)
    if(is.finite(liks[iter,(2*Ti)+2])){
      As = (liks[iter,(2*Ti)+2] - liks[iter,(2*Ti)+1])+(sum(dgamma(sigmaU,msi,vsi,log = T),na.rm = T)-sum(dgamma(sigmaC,msi,vsi,log = T),na.rm = T))+(sum(log(sigmaC), na.rm = T)-sum(log(sigmaU),na.rm=TRUE))
      As = exp(As)
    }
    if(!is.finite(liks[iter,(2*Ti)+2])){As=0}
    asi = runif(1)
    if(asi > min(1, As)){
      sigma[iter,,,] <- sigmaC
      liks[iter,(2*Ti)+2] = liks[iter,(2*Ti)+1]
    }
    # update initial probabilities
    if(verbose){print("## update be")}
    beU <- beC <- array(be[iter,,,],c(dim(be)[-1]))
    beU[,k[iter,1],1:k[iter,1]] <- beC[,k[iter,1],1:k[iter,1]] + rnorm(length(beC[,k[iter,1],1:k[iter,1]]), 0, 0.1)
    beU[,k[iter,1],1]=0 
    for(i in 1:n){
      pisU[i,k[iter,1],1:k[iter,1]] <- inprob(k[iter,1], be=beU[,k[iter,1],1:k[iter,1]], Z[i,1,], thB = 10^-5)
    }
    be[iter,,,] <- beU
    
    ilik <- sapply(1:n,function(i) likco(y = y[i,,], w=1, xi = xi[iter,,,], sigma = sigma[iter,,,], pi = pisU[i,,], PI = PIs[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti))
    liks[iter,(2*Ti)+3] = sum(ilik)
    if(is.finite(liks[iter,(2*Ti)+3])) {
      Ap = (liks[iter,(2*Ti)+3]-liks[iter,(2*Ti)+2])+(sum(dnorm(beU, mbe, vbe, log = T), na.rm = T)-sum(dnorm(beC,mbe,vbe,log = T), na.rm = T))
      Ap = exp(Ap)}
    if(!is.finite(liks[iter,(2*Ti)+3])) {Ap=0}
    abe <- runif(1)
    if(abe > Ap){
      be[iter,,,] = beC
      liks[iter,(2*Ti)+3] = liks[iter,(2*Ti)+2]
    }else{
      pis = pisU
    }
    
    # update transition probabilities
    if(verbose){print("## update Be")}
    BeC <- array(Be[iter,,,,,],c(dim(Be)[-1]))
    for(j in 1:kmax) {
      for(h in 1:kmax) {
        trans=k[iter,-Ti]==j & k[iter,-1]==h
        
        if(any(trans)) {
          PIsU=PIs
          BeU[,j,h,1:j,1:h] <- BeC[,j,h,1:j,1:h]
          for(s in 1:h){
            if(n.cov>1) {
              BeU[,j,h,1:j,1:h][,-s,s] <- BeC[,j,h,1:j,1:h][,-s,s] + rnorm(length(BeC[,j,h,1:j,1:h][,-s,s]),0,0.1)}
            if(n.cov==1) {BeU[,j,h,1:j,1:h][-s,s] <- BeC[,j,h,1:j,1:h][-s,s] + rnorm(length(BeC[,j,h,1:j,1:h][-s,s]),0,0.1)
            }
          }
          Be[iter,,,,,] <- BeU
          for(i in 1:n){
            for(t in which(trans)){
              PIsU[i,t,j,h,1:j,1:h]= tmat(K = c(j,h), Be[iter,,j,h,1:j,1:h],Z[i,t+1,],10^-5)
            }
            ilik[i] <- likco(y = y[i,,], w=1, xi = xi[iter,,,], sigma = sigma[iter,,,], pi = pis[i,,], PI = PIsU[i,,,,,], k=k[iter,], kmax=kmax, Ti=Ti) 
          }
          liks[iter,pseq] = sum(ilik)
          if(is.finite(liks[iter,(2*Ti)+4])) {
            AP = (liks[iter,pseq]-liks[iter,pseq-1])+(sum(dnorm(c(BeU), mBe, vBe, log = T), na.rm = T)-sum(dnorm(c(BeC),mBe,vBe,log = T), na.rm = T))
            AP=exp(AP)}
          if(!is.finite(liks[iter,pseq])) {AP=0}
          aBe = runif(1)
          if(aBe <= AP){
            PIs <- PIsU
          }else{
            Be[iter,,,,,] <- BeC
            liks[iter,pseq] <- liks[iter,pseq-1]
          }
        }}}
    
    if(verbose) {print(c(iter/n.iter, liks[iter,],k[iter,]))}
  }
  
  return(list(liks=liks,k=k,be=be,Be=Be,xi=xi,sigma=sigma))
}

# This function is used to generate data
genData=function(n=200,k=c(4,4,4,3,3,3),r=3,sepp=4,be,Be) {
  
  Ti=length(k)
  kmax=max(k)
  n.cov=dim(be)[1]-1
  
  pi=array(NA,c(n,kmax,kmax))
  PI=array(NA,c(n,kmax,kmax,kmax,kmax))
  pi[,1,1]=PI[,,1,,1]=1
  
  Z = array(NA,c(n,Ti,n.cov+1))
  Z[,,1] = 1
  Z[,,-1] = rnorm(length(Z[,,-1]))
  
  for(i in 1:n){
    pi[i,k[1],1:k[1]] <- inprob(k[1], be=be[,k[1],1:k[1]], Z[i,1,], thB = 10^-5)
    for(tt in 2:Ti){
      PI[i,k[tt-1],k[tt],1:k[tt-1],1:k[tt]]= tmat(K = c(k[tt-1],k[tt]), Be[,k[tt-1],k[tt],1:k[tt-1],1:k[tt]],Z[i,tt,],10^-5)
    }
  }
  
  xi=array(1,c(kmax,kmax,3)) 
  sigma=xi
  for(j in 1:kmax) {
    for(u in 1:r) {
      xi[j,1:j,u]=seq(0,sepp*j,length=j)}}
  ## generate data ###
  
  y=array(NA,c(n,Ti,r))
  u=matrix(NA,n,Ti) 
  for(i in 1:n) {
    prob=pi[i,k[1],1:k[1]]
    u[i,1]=sample(k[1],1,replace=TRUE,prob=prob)
    for(j in 1:r) {y[i,1,j]=rnorm(1,xi[k[1],u[i,1],j])}
    for(ti in 2:Ti) {
      prob=PI[i,k[ti-1],k[ti],u[i,ti-1],1:k[ti]]
      u[i,ti]=sample(k[ti],1,replace=TRUE,prob=prob)
      for(j in 1:r) {y[i,ti,j]=rnorm(1,xi[k[ti],u[i,ti],j])}}
  }
  
  
  truePars=list(u=u,pi=pi,PI=PI,xi=xi,sigma=sigma,be=be,Be=Be)
  
  return(list(y=y,Z=Z,truePars=truePars))}

