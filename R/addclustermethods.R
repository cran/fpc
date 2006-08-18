kmeansruns <- function(data,k,iter.max=100,runs=100,
                       scaledata=FALSE,plot=FALSE){
  data <- as.matrix(data)
  if (scaledata) data <- scale(data)
  options(show.error.messages = FALSE)
  minSS <- Inf
  kmopt <- NULL
  for (i in 1:runs){
    repeat{
      km <- try(kmeans(data,k,iter.max=iter.max))
      if (class(km) != "try-error") break 
    }
    options(show.error.messages = TRUE)  
    swss <- sum(km$withinss)
    if (swss<minSS){
      kmopt <- km
      minSS <- swss
    }
    if (plot){
      par(ask=TRUE)
      pairs(data,col=km$cluster,main=swss)
    }
  }
  kmopt
}
    
pamk <- function(data,krange=2:10,scaling=FALSE,diss=inherits(data, "dist"),
                 ...){
  require(cluster)
  if (scaling)
    sdata <- scale(data,scale=scaling)
  else
    sdata <- data
  asw <- numeric(max(krange))
  pams <- list()
  for (k in krange){
    pams[[k]] <- pam(sdata, k,...)
    asw[k] <- pams[[k]]$silinfo$avg.width
  }
  k.best <- which.max(asw)
  out <- list(pamobject=pams[[k.best]],nc=k.best)
  out
}

  
