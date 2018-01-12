# Halkidi, Vazirgiannis, Pat Rec Let 2008

# Find representative objects and within-cluster variance

findrep <- function(x,xcen,clustering,cluster,r,p=ncol(x),n=nrow(x),
                    nc=sum(clustering==cluster)){
  repxx <- matrix(0,nrow=r,ncol=p)
  xc <- x[clustering==cluster,,drop=FALSE]
  drx <- matrix(0,nrow=r,ncol=nc)
  drx[1,] <- mahalanobis(xc,xcen,diag(p))
  wvar <- sum(drx[1,])/(nc-1)
  if (is.na(wvar)) wvar <- 0
#   browser()
  repxi <- which.max(drx[1,])
  repxx[1,] <- xc[repxi,]
  maxr <- r
#   browser()
  if (r>1){
    for (ri in 2:r){
      drx[ri,] <- mahalanobis(xc,repxx[ri-1,],diag(p))
      di <- numeric(0)
      for (i in 1:nc)
        di[i] <- min(drx[2:ri,i])
      if (max(di)>0){
        repxi[ri] <- which.max(di)
        repxx[ri,] <- xc[repxi[ri],]
      }
      else{
        maxr <- ri-1
        break
      }
    }
  }
  list(repc=(1:n)[clustering==cluster][repxi],repx=repxi,
           maxr=maxr,wvar=wvar)
}
# repc: indexes relative to x, repx: indexes relative to cluster

cdbw <- function(x,clustering,r=10,s=seq(0.1,0.8,by=0.1),
                 clusterstdev=TRUE,trace=FALSE){
  p <- ncol(x)
  n <- nrow(x)
  x <- as.matrix(x)
  cn <- max(clustering)
  repc <- list()
  repr <- rep(0,cn)
  mrepr <- vrepc <- numeric(0)
  minrepr <- 1
  xcc <- matrix(0,nrow=cn,ncol=p)
  nc <- numeric(0)
  wvar <- numeric(0)
  repx <- list()
  if (trace)
    print("Find representatives")
  for (i in 1:cn){
    nc[i] <- sum(clustering==i)
    xcc[i,] <- colMeans(x[clustering==i,,drop=FALSE])
    rrx <- findrep(x,xcc[i,],clustering,i,r,p,n,nc[i])
    repc[[i]] <- rrx$repc
    repx[[i]] <- rrx$repx
    repr[i] <- rrx$maxr
    wvar[i] <- rrx$wvar
    mrepr[i] <- sum(repr)
    if(i>1)
      minrepr[i] <- mrepr[i-1]+1
    vrepc <- c(vrepc,repc[[i]]) 
  }
  if (trace)
    print(repc)
  stdev <- mean(sqrt(wvar))
# Find the closest representatives for all representatives and the sets rcr
# for all cluster pairs
  dv <- as.matrix(dist(x[vrepc,]))
  dijmin <- list()
  rcr <- list()
  for (i in 1:cn){
    dijmin[[i]] <- list()
    rcr[[i]] <- list()
  }
#  browser()
  for (i in 1:(cn-1)){
    for (j in (i+1):cn){
      dij <- dv[minrepr[i]:mrepr[i],minrepr[j]:mrepr[j],drop=FALSE]
      ii <- ij <- numeric(0)
      dijmin[[i]][[j]] <- dijmin[[j]][[i]] <- numeric(0)
      for (k in 1:repr[i]){
        ii[k] <- which.min(dij[k,])
        dijmin[[i]][[j]][k] <- repc[[j]][ii[k]]
      }
      for (k in 1:repr[j]){
        ij[k] <- which.min(dij[,k])
        dijmin[[j]][[i]][k] <- repc[[i]][ij[k]]
      }
      rcr[[i]][[j]] <- numeric(0)
      for (k in 1:repr[i]){
        if (k==ij[ii[k]])
          rcr[[i]][[j]] <- rbind(rcr[[i]][[j]],
                                 c(dijmin[[i]][[j]][k],dijmin[[j]][[i]][ii[k]]))
      }
    }
  }
  if(trace){
    print("rcr")
    print(rcr)
    print("wvar")
    print(wvar)
    print("stdev")
    print(stdev)
  }
# Find dens(C_i,C_j), dist
  dens <- dk <- matrix(0,ncol=cn,nrow=cn)
  for (i in 1:(cn-1)){
    for (j in (i+1):cn){      
      nrcr <- nrow(rcr[[i]][[j]])
      wsdij <- mean(sqrt(wvar[i]),sqrt(wvar[j])) 
      for (k in 1:nrcr){
        u <- (x[rcr[[i]][[j]][k,1],]+x[rcr[[i]][[j]][k,2],])/2
#        browser()
        ud <- sqrt(mahalanobis(x[clustering==i | clustering==j,,drop=FALSE],u,diag(p)))
        dkd <-  sqrt(sum((x[rcr[[i]][[j]][k,1],]-x[rcr[[i]][[j]][k,2],])^2))
        dk[i,j] <- dk[i,j]+dkd
        dens[i,j] <- dens[i,j]+dkd*sum(ud<wsdij)/(2*wsdij*(nc[i]+nc[j]))
      }
      dens[j,i] <- dens[i,j] <- dens[i,j]/nrcr
      dk[j,i] <- dk[i,j] <- dk[i,j]/nrcr
    }
  }
  if(trace){
    print("dens")
    print(dens)
  }
# Interdens and Sep
  maxd <- mind <- numeric(0)
  for (i in 1:cn){
    maxd[i] <- max(dens[i,])
    mind[i] <- min(dk[i,-i])
  }
  interdens <- mean(maxd)
  sep <- mean(mind/(1+interdens))
  if (trace)
    cat("sep= ",sep," interdens=",interdens," mind=",mind,"\n")
# Intradens and compactness
  ns <- length(s)
  intradens <- numeric(0)
  denscl <- matrix(0,nrow=cn,ncol=ns) 
  for (i in 1:ns){
    for (j in 1:cn){
#     browser()
      xcj <- x[clustering==j,,drop=FALSE]
      if (clusterstdev)
        stdevj <- sqrt(wvar[j])
      else
        stdevj <- stdev
#      dxj <- as.matrix(dist(x))[clustering==j,clustering==j]
      for (k in 1:repr[j]){
        srep <- (1-s[i])*xcj[repx[[j]][k],]+s[i]*xcc[j,]
        dsjk <- mahalanobis(xcj,srep,diag(p))
        denscl[j,i] <- denscl[j,i]+sum(dsjk<stdevj)/nc[j]
      }
      denscl[j,i] <- denscl[j,i]/repr[j]
      if (trace)
        cat("denscl cluster ",j," s ",i,": ",denscl[j,i],"\n")
    }
    intradens[i] <- sum(denscl[,i])/(cn*stdev)
  }
  compactness <- mean(intradens)
  if (trace){
    print(intradens)
    cat("compactness= ",compactness,"\n")
  }
# Intrachange and Cohesion
  ic <- intradens[2:ns]-intradens[1:(ns-1)]
  intrachange <- sum(ic)/(ns-1)
  cohesion <- compactness/(1+intrachange)
  sc <- sep*compactness
  cdbw <- cohesion*sc
  if (trace){
    print(intrachange)
    cat("cohesion= ",cohesion," sc=",sc," cdbw=",cdbw,"\n")
  }
  list(cdbw=cdbw,cohesion=cohesion,compactness=compactness,sep=sep)
}
   
      
# Negentropy increment,  Lago-Fernandez and Corbacho (PatRec 2010)

neginc <- function(x,clustering){  
  p <- ncol(x)
  n <- nrow(x)
  x <- as.matrix(x)
  cn <- max(clustering)
  dcovx <- det(cov(x))
  pri <- numeric(0)
  dcovxi <- numeric(0)
  for (i in 1:cn){
    pri[i] <- sum(clustering==i)/n
    dcovxi[i] <- det(cov(x[clustering==i,,drop=FALSE]))
    if (is.na(dcovxi[i])) dcovxi[i] <- 0
  }
  out <- 0.5*sum(pri*log(dcovxi))-0.5*log(dcovx)-sum(pri*log(pri))
  out
}

# CVNN index by Liu et al. 2013, in Aggarwals book, but "corrected", see handbook

# clusterings should be a list of clustering vectors where the max is the nc.
cvnn <- function(d=NULL,clusterings,k=5){
  d <- as.matrix(d)
  n <- nrow(d)
  neighb <- matrix(0,nrow=n,ncol=k)
  for (i in 1:n)
    neighb[i,] <- order(d[i,])[2:(k+1)]
  cnums <- sep <- comp <- sepf <- compf <- numeric(0)
  lc <- length(clusterings)
  for (i in 1:lc){
    cnums[i] <- max(clusterings[[i]])
    comp[i] <- nji <- 0
    sepj <- numeric(0)
    for (j in 1:cnums[i]){
      nj <- sum(clusterings[[i]]==j)
      compj <- sum(d[clusterings[[i]]==j,clusterings[[i]]==j])
      nji <- nji+nj*(nj-1)
      comp[i] <- comp[i]+compj
#      cat("i=",i," j=",j,"compj=",compj,"\n")
      sepj[j] <-  sum(clusterings[[i]][as.vector(neighb[(1:n)[clusterings[[i]]==j],])]!=j)/(k*nj)
    }
    sep[i] <- max(sepj)
    comp[i] <- comp[i]/nji
  }
  maxsep <- max(sep)
  maxcomp <- max(comp)
  sepf <- sep/maxsep
  compf <- comp/maxcomp
  cvnnindex <- sepf+compf
  out <- list(cvnnindex=cvnnindex,sep=sep,comp=comp)
}
    
    
  

