tdecomp <- function(m){
  wm <- eigen(m, symmetric=TRUE)
  p <- ncol(m)
  wmd <- wm$values
  for (i in 1:p){
    if (abs(wmd[i])<1e-6)
      wmd[i] <- 1e-6
  }
  out <- t(wm$vectors %*% diag(sqrt(wmd)))
  out
}

discrcoord <- function(xd, clvecd, pool="n", ...) {
  x <- as.matrix(xd)
  clvec <- as.integer(clvecd)
  n <- nrow(x)
  p <- ncol(x)
  clf <- factor(clvec)
  cll <- as.integer(levels(clf)) 
  clnum <- length(cll)
  cln <- rep(0, times=clnum)
  for (i in 1:clnum){
    cln[i] <- sum(clvec==cll[i])
  }
  W <- rep(0, times=p*p)
  dim(W) <- c(p,p)
  for (i in 1:clnum){
    clx <- rep(0, times=p*cln[i])
    dim(clx) <- c(cln[i],p)
    for (j in 1:p){
      clx[,j] <- x[,j][clvec==cll[i]]
    }
    if (pool=="n")
      W <- W + ((cln[i]-1)*cov(clx))
    else
      W <- W + (n-1)*cov(clx)/clnum
  }
  Tm <- tdecomp(W)
  Tinv <- solve(Tm)
  S <- (n-1)*cov(x)
  B <- S-W
  Z <- t(Tinv) %*% B %*% Tinv
  dc <- eigen(Z, symmetric=TRUE)
  units <- Tinv %*% dc$vectors * sqrt(n-clnum)
  proj <- x %*% units    
  list(ev=dc$values, units=units, proj=proj, W=W)
}

batvarcoord <- function(xd, clvecd, clnum=1){
  x <- as.matrix(xd)
  clvec <- as.integer(as.integer(clvecd)==as.integer(clnum))
  n <- nrow(x)
  p <- ncol(x)
  cll <- c(0,1)
  cln <- rep(0, times=2)
  for (i in 1:2)
    cln[i] <- sum(clvec==cll[i])
  clx <- list()
  for (i in 1:2){
    clx[[i]] <- rep(0, times=p*cln[i])
    dim(clx[[i]]) <- c(cln[i],p)
    for (j in 1:p)
      clx[[i]][,j] <- x[,j][clvec==cll[i]]
  }
  S1 <- cov(clx[[1]])
  S2 <- cov(clx[[2]])
  W <- solve(S2) %*% S1
  Weigen <- eigen(W)
  rev <- Weigen$values + 1/Weigen$values + 2
  dw <- diag(t(Weigen$vectors) %*% S2 %*% Weigen$vectors)
  svw <- matrix(ncol=p, nrow=p)
  for (i in 1:p)
    svw[,i] <- Weigen$vectors[,i]/sqrt(dw[i])
  units <- svw[,(order(-rev))]
  proj <- x %*% units    
  list(ev=Weigen$values[order(-rev)], rev=rev[(order(-rev))], 
       units=units, proj=proj, W=W, S1=S1, S2=S2)
}

batcoord <- function(xd, clvecd, clnum=1, dom="mean"){
  x <- as.matrix(xd)
  clvec <- as.integer(clvecd)
  clf <- factor(clvec)
  cll <- as.integer(levels(clf)) 
  clz <- length(cll)
  p <- ncol(x)
  if (clz!=2){
    clvec <- as.integer(clvecd==clnum)
    cll <- c(0,1)
    print("Cluster indicator has more than 2 values")
  }
  if (dom=="mean"){
    dcx <- discrcoord(x, clvec, pool="equal")
    x2 <- dcx$proj[,2:p]
    batx <- batvarcoord(x2, as.integer(clvecd==clnum))
    ev <- c(dcx$ev[1],batx$ev)
    rev <- c(max(batx$rev)+1, batx$rev)
    units <- cbind(dcx$units[,1],dcx$units[,2:p] %*% batx$units)
    proj <- cbind(dcx$proj[,1],batx$proj)
  }
  else{
    batx <- batvarcoord(x, as.integer(clvecd==clnum))
    cln <- rep(0, times=2)
    mx <- matrix(nrow=nrow(x),ncol=p)
    for (i in 1:2)
      cln[i] <- sum(clvec==cll[i])
    for (i in 1:2){
      clx <- rep(0, times=p*cln[i])
      dim(clx) <- c(cln[i],p)
      for (j in 1:p){
        clx[,j] <- x[,j][clvec==cll[i]]
      }
      mx[i,] <- colMeans(clx)
    }
    mdiff <- mx[2,]-mx[1,]
    ev <- batx$ev
    rev <- rep(0, times=p)
    for (i in 1:p)
      rev[i] <- (batx$units[,i] %*% mdiff)^2/(1+ev[i])+log(ev[i]+1/ev[i]+2)
    units <- batx$units[,(order(-rev))]
    proj <- batx$proj[,(order(-rev))]
  }    
  list(ev=ev[order(-rev)], rev=rev[order(-rev)], 
       units=units, proj=proj)
}

# discriminant plot
# bw: black/white
plotcluster <- function(x, clvecd, clnum=1,
                        method=ifelse(identical(range(as.integer(clvecd)),
                          as.integer(c(0,1))),"awc","dc"),bw=FALSE, xlab=NULL,
                        ylab=NULL, pch=NULL, col=NULL, ...){
  asym <- any(method==c("bc","vbc","adc","awc","arc","anc"))
  if (asym)
    clvec <- as.integer(as.integer(clvecd)==as.integer(clnum))
  else
    clvec <- as.integer(clvecd)
  cx <- discrproj(x, clvecd, method, clnum, ...)$proj
  if (is.null(xlab))
    xlab <- paste(method,"1")
  if (is.null(ylab))
    ylab <- paste(method,"2")
  if (is.null(pch))
    pch <- if (bw){
             1+clvec 
           }
           else 1
  if (is.null(col))
    col <- if (bw) 1
           else{
             1+clvec
           }
  plot(cx, xlab=xlab, ylab=ylab, pch=pch, col=col, ...)
}

discrproj <- function(x, clvecd, method="awc", clnum=1, ...){
  result <- switch(method,
                   dc=discrcoord(x, clvecd, ...),
                   bc=batcoord(x, clvecd, clnum),
                   vbc=batcoord(x, clvecd, clnum, dom="var"),
                   mvdc=mvdcoord(x, clvecd, clnum, ...),
                   adc=adcoord(x, clvecd, clnum),
                   awc=awcoord(x, clvecd, clnum, ...),
                   arc=awcoord(x, clvecd, clnum, method="mcd", ...),
                   nc=ncoord(x, clvecd, ...),
                   wnc=ncoord(x, clvecd, weighted=TRUE, ...),
                   anc=ancoord(x, clvecd, clnum, ...))
  result
}
                   
                   







