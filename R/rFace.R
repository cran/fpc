#
# face benchmark dataset by M. Maechler and C. Hennig
#
## MM: -  function(n, p)  {where `n' was a bit tricky}
##     -  return grouping as well

rFace <- function(n, p = 6, nrep.top = 2, smile.coef = 0.6,
                  dMoNo = 1.2, dNoEy = 1)
{
    ## Purpose: Generate random "Face" data set -- to be "Hard for Clustering"
    ## -------------------------------------------------------------------------
    ## Arguments: (n,p)    : dimension of result
    ##            nrep.top : #{repetitions} of the top point / "hair tip"
    ##            dMoNo    : distance{Mouth, Nose}
    ##            dNoEy    : distance{Nose,  Eyes} {vertically only}
    ##
    ## MM- TODO's:   o  provide more arguments (face shape)
    ##			    --> separation of clusters
    ## -------------------------------------------------------------------------
    ## Author: Christian Hennig & Martin Maechler, 26 Jun 2002

    if((p <- as.integer(p)) < 2) stop("number of variables p must be at least 2")
    if((n <- as.integer(n)) < 10) stop("number of points  n  must be at least 10")
    if((nrep.top <- as.integer(nrep.top)) < 1) stop("`nrep.top' must be positive")

    ntips <- nrep.top + 2
    n0 <- n - ntips
    m <- n0 %/% 5
    n.5 <- n0 %/% 2
    n.7 <- n.5+m
    n.9 <- n.7+m
    ## shouldn't happen:
    if(m < 1 || n.9 >= n0) stop("number of points n is too small")

    ## Indices of the different groups :
    Gr <- list(chin = 1:m,
               mouth= (m+1):n.5,
               nose = (n.5+1):n.7,
               rEye = (n.7+1):n.9,
               lEye = (n.9+1):n0,
               tips = (n0+1):n)

    face <- matrix(nrow = n, ncol = p)
    ## chin :
    face[Gr$chin, 1] <- U <- runif(m, -3, 3)
    face[Gr$chin, 2] <- rnorm(m, mean = U^2, sd=0.1)
    ## mouth:
    m0m <- 3 # lower mouth mean
    face[Gr$mouth, 1] <- Z <- rnorm(n.5-m, sd= 0.5)
    face[Gr$mouth, 2] <- rnorm(n.5-m, mean = m0m + smile.coef * Z^2, sd= 0.2)
    ## nose :
    face[Gr$nose, 1] <-  rnorm(m, mean=0, sd=0.2)
    yEye <- 17
    ##face[Gr$nose, 2] <- rnorm(m, mean=9, sd= 2.5)
    face[Gr$nose, 2] <- pmin(yEye - dNoEy,
                             (m0m + dMoNo) + rgamma(m, shape= 0.8, scale= 2.5))
    ## right eye  U[circle] :
    rangle <- runif(m, 0, 2*pi)
    rpos   <- runif(m)
    face[Gr$rEye, 1] <- rpos*cos(rangle) +  2
    face[Gr$rEye, 2] <- rpos*sin(rangle) + yEye
    ## left eye:
    face[Gr$lEye, 1] <- rnorm(n0-n.9, mean= -2,   sd= 0.5)
    face[Gr$lEye, 2] <- rnorm(n0-n.9, mean= yEye, sd= 0.5)
    ## `ntips'  ``hair tips'', the last one `nrep.top' times:
    face[n0+1,1:2] <- c(-4.5, 25)
    face[n0+2,1:2] <- c( 4.5, 25)
    for(k in 1:nrep.top)
        face[n-k+1, 1:2] <- c(0,32)

    ##-- Extra coordinates with noise ---

    if(p >= 3) {
        face[, p] <- rexp(n)
        if(p >= 4) {
            face[, p-1] <- rt(n,df=1)
            if(p >= 5)
                for(k in 3:(p-2))
                    face[,k] <- rnorm(n)
        }
    }
    gr <- character(n)
    ng <- names(Gr)
    for(i in seq(ng)) gr[Gr[[i]]] <- ng[i]
    structure(face, grouping = as.factor(gr), indexlist = Gr)
}

