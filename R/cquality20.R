# cquality 19 implements Serhat's stupid clustering methods
# cquality 20 allows for bootstrap stability

# Re-defined average.within with reweighting!
# Changed print.valstat output to be in line with documentation, but
# can make documentation clearar

# standardisation is for ave distance, separation, 
# largest within gap, 
# pamcrit, max.diameter, min.separation
# and can be "none", "max" (max overall distance), "ave" (ave overall distance),
# "q90" (90% quantile).
# Unless stan="none", densityindex, parsimony and
# entropy will be standardised by the maximum possible for given number
# clusters (entropy), maxk (parsimony); within and
# between cluster ss will be standardised by overall.ss, mnnd will be
# standardised by maximum distance to nnkth nearest neighbour
# standardisation can also be a number. pearsongamma will be standardised
# +1/2. dindex, denscut are automatically standardised, dhgap by maxd.
# cvnn by sqrt(n)

# sepall: Should every cluster be taken into account for sepindex? Or just best
# sepprob overall?

# nndist: compute average distance to nnkth nearest neighbour within cluster (mnnd)
# nnk will then also govern the computation of coefficient of variation of
# kth within-cluster nn distance; clusters with <nnk+1 points are ignored
# averagegap: wgap is average within cluster gaps, otherwise max
# pamcrit: compute pam criterion pamc (looking for optimal centroids pamcentroids)
# densityindex: density estimate is points close to
# dfactor*within cluster median distance, weighted down linearly

# maxk max number of clusters for parsimony
# cvstan: standardisation for coeffeicient of variation 
cqcluster.stats <- function (d = NULL, clustering, alt.clustering = NULL,
                             noisecluster = FALSE, 
    silhouette = TRUE, G2 = FALSE, G3 = FALSE, wgap = TRUE, sepindex = TRUE, 
    sepprob = 0.1, sepwithnoise = TRUE, compareonly = FALSE, 
    aggregateonly = FALSE,
# new:                             
    averagegap=FALSE, pamcrit=TRUE,
#                             densityindex=TRUE,
    dquantile=0.1,
    nndist=TRUE, nnk=2, standardisation="max", sepall=TRUE, maxk=10,
    cvstan=sqrt(length(clustering)))
# dk=5, doweight=0.25,    
{
  lweight <- function(x,md)
    (x<md)*(-x/md+1)

## preparations    
    if (!is.null(d)) 
        d <- as.dist(d)
    cn <- max(clustering)
    clusteringf <- as.factor(clustering)
    clusteringl <- levels(clusteringf)
    cnn <- length(clusteringl)
    if (cn != cnn) {
        warning("clustering renumbered because maximum != number of clusters")
        for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
        cn <- cnn
    }
    n <- length(clustering)
    noisen <- 0
    cwn <- cn
    if (noisecluster) {
        noisen <- sum(clustering == cn)
        cwn <- cn - 1
    }
    parsimony <- cn/maxk
# cn: number of clusters including noise; cwn: number of clusters w/o noise
    diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <- within.dist <- between.dist <- numeric(0)
    for (i in 1:cn) cluster.size[i] <- sum(clustering == i)
## standardisation    
    if (is.numeric(standardisation)) stan <- standardisation
    else
      stan <- switch(standardisation,
                   max=max(d),
                   ave=mean(d),
                   q90=quantile(d,0.9),
                   1)
## entropy, corrected rand, vi    
    pk1 <- cluster.size/n
    pk10 <- pk1[pk1 > 0]
    h1 <- -sum(pk10 * log(pk10))
    if (!(standardisation=="none")){
      pkmax <- rep(1,cn)/cn
      h1 <- -h1/sum(pkmax * log(pkmax))
    }
    corrected.rand <- vi <- NULL
    if (!is.null(alt.clustering)) {
        choose2 <- function(v) {
            out <- numeric(0)
            for (i in 1:length(v)) out[i] <- ifelse(v[i] >= 2, 
                choose(v[i], 2), 0)
            out
        }
        cn2 <- max(alt.clustering)
        clusteringf <- as.factor(alt.clustering)
        clusteringl <- levels(clusteringf)
        cnn2 <- length(clusteringl)
        if (cn2 != cnn2) {
            warning("alt.clustering renumbered because maximum != number of clusters")
            for (i in 1:cnn2) alt.clustering[clusteringf == clusteringl[i]] <- i
            cn2 <- cnn2
        }
        nij <- table(clustering, alt.clustering)
        dsum <- sum(choose2(nij))
        cs2 <- numeric(0)
        for (i in 1:cn2) cs2[i] <- sum(alt.clustering == i)
        sum1 <- sum(choose2(cluster.size))
        sum2 <- sum(choose2(cs2))
        pk2 <- cs2/n
        pk12 <- nij/n
        corrected.rand <- (dsum - sum1 * sum2/choose2(n))/((sum1 + 
            sum2)/2 - sum1 * sum2/choose2(n))
        pk20 <- pk2[pk2 > 0]
        h2 <- -sum(pk20 * log(pk20))
        icc <- 0
        for (i in 1:cn) for (j in 1:cn2) if (pk12[i, j] > 0) 
            icc <- icc + pk12[i, j] * log(pk12[i, j]/(pk1[i] * 
                pk2[j]))
        vi <- h1 + h2 - 2 * icc
    }
    if (compareonly) {
        out <- list(corrected.rand = corrected.rand, vi = vi)
    }
    else { # begin !if compareonly
## functions of within/between cluster distances        
#        if (silhouette) 
#            require(cluster)
        dmat <- as.matrix(d)
        within.cluster.ss <- 0
        overall.ss <- nonnoise.ss <- sum(d^2)/n
        if (noisecluster) 
            nonnoise.ss <- sum(as.dist(dmat[clustering <= cwn, 
                clustering <= cwn])^2)/sum(clustering <= cwn)
        ave.between.matrix <- separation.matrix <- matrix(0, 
            ncol = cn, nrow = cn)
        nnd <- numeric(0)
        cvnndc <- rep(NA,cn)
        mnnd <- cvnnd <- NULL
        di <- list()
        for (i in 1:cn) { 
            cluster.size[i] <- sum(clustering == i)
            di <- as.dist(dmat[clustering == i, clustering == 
                i])
            if (i <= cwn) {
                within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
                within.dist <- c(within.dist, di)
            }
            if (length(di) > 0){ 
                diameter[i] <- max(di)                
                average.distance[i] <- mean(di)
                median.distance[i] <- median(di)
            }
            else diameter[i] <- average.distance[i] <- median.distance[i] <- NA
            bv <- numeric(0)
            for (j in 1:cn) {
                if (j != i) {
                  sij <- dmat[clustering == i, clustering == 
                    j]
                  bv <- c(bv, sij)
                  if (i < j) {
                    separation.matrix[i, j] <- separation.matrix[j, 
                      i] <- min(sij)
                    ave.between.matrix[i, j] <- ave.between.matrix[j, 
                      i] <- mean(sij)
                    if (i <= cwn & j <= cwn) 
                      between.dist <- c(between.dist, sij)
                  } # if i<j
                } # if j!=i
            } # for j
            separation[i] <- min(bv)
            average.toother[i] <- mean(bv)
        } # for i
## nndist  and coef of var
        if (nndist){
          kenough <- cluster.size>nnk
#          if (min(cluster.size)<nnk+1)
#            warning("min cluster size smaller than nnk. No distance to neighbours computed.") 
          for (i in (1:cn)[kenough]){
              nndi <- apply(dmat[clustering==i,clustering==i],
                     1,sort,partial=nnk+1)[nnk+1,]
#              print(nndi)
              nnd <- c(nnd,nndi)
              cvnndc[i] <- sd(nndi)/mean(nndi)
          }
          cvnnd <- weighted.mean(cvnndc,pk1,na.rm=TRUE)
          mnnd <- mean(nnd,na.rm=TRUE)
#          print("cluster.size")
#          print(cluster.size)
#          print(nnk)
#          print(kenough)
#          print(nnd)
#          print(mnnd)
          if (!standardisation=="none"){
                maxnnd <- max(apply(dmat,1,sort,partial=nnk+1)[nnk+1,])
#                print("stan")
#                print(maxnnd)
                mnnd <- mnnd/maxnnd
                cvnnd <- cvnnd/cvstan
          }
        }
## end nndist         
        average.between <- mean(between.dist)/stan
#        if (stanbound) average.between <- min(1,average.between)
#        average.within <- mean(within.dist)/stan
        average.within <- weighted.mean(average.distance,cluster.size,na.rm=TRUE)/stan
        nwithin <- length(within.dist)
        nbetween <- length(between.dist)
        between.cluster.ss <- nonnoise.ss - within.cluster.ss
        ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss * 
            (cwn - 1))
        if (!(standardisation=="none")){
          within.cluster.ss <- within.cluster.ss/overall.ss
          between.cluster.ss <- between.cluster.ss/overall.ss
        }
#        if (stanbound) between.cluster.ss <- min(1,between.cluster.ss)
## silhouette widths        
        clus.avg.widths <- avg.width <- NULL
        if (silhouette) {
            sii <- silhouette(clustering, dmatrix = dmat)
            sc <- summary(sii)
            clus.avg.widths <- sc$clus.avg.widths
            if (noisecluster) 
                avg.width <- mean(sii[clustering <= cwn, 3])
            else avg.width <- sc$avg.width
        } # if silhouette
## g2        
        g2 <- g3 <- cn2 <- cwidegap <- widestgap <- sindex <- NULL
        if (G2) {
            splus <- sminus <- 0
            for (i in 1:nwithin) {
                splus <- splus + sum(within.dist[i] < between.dist)
                sminus <- sminus + sum(within.dist[i] > between.dist)
            }
            g2 <- (splus - sminus)/(splus + sminus)
        }
## g3        
        if (G3) {
            sdist <- sort(c(within.dist, between.dist))
            sr <- nwithin + nbetween
            dmin <- sum(sdist[1:nwithin])
            dmax <- sum(sdist[(sr - nwithin + 1):sr])
            g3 <- (sum(within.dist) - dmin)/(dmax - dmin)
        }
## pearsongamma, dunn, dunn2        
        pearsongamma <- cor(c(within.dist, between.dist), c(rep(0, 
            nwithin), rep(1, nbetween)))
        dunn <- min(separation[1:cwn])/max(diameter[1:cwn], na.rm = TRUE)
        acwn <- ave.between.matrix[1:cwn, 1:cwn]
        dunn2 <- min(acwn[upper.tri(acwn)])/max(average.distance[1:cwn])
        if (!standardisation=="none")
          pearsongamma <- (pearsongamma+1)/2
## widest within gap         
        if (wgap) {
            cwidegap <- rep(0, cwn)
            for (i in 1:cwn) if (sum(clustering == i) > 1) 
                cwidegap[i] <- max(hclust(as.dist(dmat[clustering == 
                  i, clustering == i]), method = "single")$height)
            if (averagegap)
              widestgap <- mean(cwidegap)/stan
            else    
              widestgap <- max(cwidegap)/stan
        }
## separation index        
        if (sepindex) {
          if (sepall){
            psep <- rep(NA,n)
            svec <- numeric(0)
            for (i in 1:n) psep[i] <- min(dmat[i, clustering != 
                  clustering[i]])
            lcn <- cwn+(sepwithnoise & noisecluster)
            for (j in 1:lcn){
              scj <- psep[clustering==j]
              qcj <- quantile(scj,sepprob)
              svec <- c(svec,scj[scj<=qcj])
            } # for j
            sindex <- mean(svec)/stan
          } # if sepall
          else{
            psep <- rep(NA, n)
            if (sepwithnoise | !noisecluster) {
                for (i in 1:n) psep[i] <- min(dmat[i, clustering != 
                  clustering[i]])
            } # if (sepwithnoise | !noisecluster)
            else {
                dmatnn <- dmat[clustering <= cwn, clustering <= 
                  cwn]
                clusteringnn <- clustering[clustering <= cwn]
                for (i in 1:(n - noisen)) psep[i] <- min(dmatnn[i, 
                  clusteringnn != clusteringnn[i]])
            } # else (!if (sepwithnoise | !noisecluster))
            sindex <- mean(psep[psep<=quantile(psep, sepprob)])/stan
          } # else (!if sepall)
        } # if sepindex
## pam criterion
        pamc <- pamcentroids <- NULL
        if (pamcrit){
          pamc <- 0
          pamcentroids <- rep(0,cwn)
          for (i in 1:cwn){
            mindsum <- Inf
            for (j in (1:n)[clustering==i]){
              dsum <- sum(dmat[j,clustering==i])
              if (dsum<=mindsum){
                mindsum <-  dsum
                pamcentroids[i] <- j
              } # if dsum best up to now
            } # for j
            pamc <- pamc+mindsum 
          } # for i
          pamc <- pamc/(n*stan)
        } # if pamcrit
## densityindex
        withindensp <- densoc <- percdensoc <- percwdens <- rep(0,n)
        thgap <- numeric(0)
  # withindensp: density
  # densoc: contribution to density that comes from other clusters
  # percdensoc: densoc/max(withindensp)
  # percwdens: withindensp/max(withindensp)      
        npenalty <- dpenalty <- rep(0,cwn)
        cutdist <- quantile(d,dquantile)
        for (i in 1:cwn){
          for (j in (1:n)[clustering==i]){
#              if (cluster.size[i]>1){
               withindensp[j] <- sum(lweight(dmat[j,],cutdist))
               densoc[j] <- sum(lweight(dmat[j,clustering!=i],cutdist))
#              }
#              else{
#                mdd <- min(d[d>0])/2
#                withindensp[j] <- sum(lweight(dmat[j,],mdd))
#                densoc[j] <- sum(lweight(dmat[j,clustering!=i],mdd))
#              }
          } # for j
        } # for i
        percwdens <- withindensp/max(withindensp)
        percdensoc <- densoc/max(withindensp)
        for (i in 1:cwn){
          npenalty[i] <- sum(percdensoc[clustering==i]*percwdens[clustering==i])
        } # for i
        pdistto <- distto <- pclosetomode <- list()
  # In the following:
  # modei: Mode of cluster i
  # closetomode: sequence of points with lowest distance to any point already
  # in closetomode (first is modei)
  # distto: difference in relative densities (percwdens) between new point
  # in closetomode and the one to which it is closest.
  # If this is positive, it is added to dpenalty.
  # pdistto is the number of the old point in closetomode that is closest.
        for (i in 1:cwn){
          di <- dmat[clustering==i,clustering==i,drop=FALSE]
          oindexes <- (1:n)[clustering==i]
          modei <- which.max(withindensp[clustering==i])
          di[modei,] <- Inf
          closetomode <- modei
          pdistto[[i]] <- distto[[i]] <- NA
          if (cluster.size[i]>1){
            for (j in 2:cluster.size[i]){
              minval <- apply(di[,closetomode,drop=FALSE],1,min)
              remainindexes <- oindexes[-closetomode]
              closetomode[j] <- which.min(minval)
              minctm <- which.min(di[closetomode[j],closetomode[1:(j-1)]])
              minctm <- oindexes[closetomode[1:(j-1)][minctm]]
              newp <- oindexes[closetomode[j]]
              thgap <- c(thgap,min(minval)*max(percwdens[remainindexes]))
              distto[[i]][j] <- percwdens[newp]-percwdens[minctm]
              if (percwdens[newp]>percwdens[minctm])
                dpenalty[i] <- dpenalty[i]+(distto[[i]][j])^2
              di[closetomode[j],] <- Inf
              pdistto[[i]][j] <- minctm
            } # for j
          } # if
          pclosetomode[[i]] <- oindexes[closetomode]
        } # for i
        dindex <- sqrt(sum(dpenalty)/n)
        denscut <- sum(npenalty)/n
        highdgap <- max(thgap)
        if (!standardisation=="none"){
          highdgap <- highdgap/stan
        }
        
# ## new density index, no dismissed
#         withindensp <- densoc <- percdensoc <- percwdens <- rep(0,n)
# # withindensp: density
# # densoc: contribution to density that comes from other clusters
# # percdensoc: densoc/max(withindensp within cluster)
# # percwdens: withindensp/max(withindensp within cluster)      
#         denscutc <- rep(0,cwn)
# # denscutc: clusterwise penalty for density contributions from other clusters
# # denscut: sum over all clusters      
#         cutdist <- quantile(d,dquantile)
#         for (i in 1:cwn){
#           for (j in (1:n)[clustering==i]){
# #            if (cluster.size[i]>1){
#               withindensp[j] <- sum(lweight(dmat[j,],cutdist))/n
#               densoc[j] <- sum(lweight(dmat[j,clustering!=i],cutdist))/n
# #            }
# #            else{
# #              mdd <- min(d[d>0])/2
# #              withindensp[j] <- sum(lweight(dmat[j,],mdd))
# #              densoc[j] <- sum(lweight(dmat[j,clustering!=i],mdd))
# #            }
#           } # for j
#           percwdens[clustering==i] <-
#             withindensp[clustering==i]/max(withindensp[clustering==i])
#           percdensoc[clustering==i] <-
#             densoc[clustering==i]/max(withindensp[clustering==i])
#           denscutc[i] <- sum(percdensoc[clustering==i]*percwdens[clustering==i])
#         } # for i
#         denscut <- sum(denscutc)/n
#         dindexc <- rep(0,cwn)
# # Clusterwise sum of squared differences between ordered density and
# # density ordered according to distance from mode.      
#         for (i in 1:cwn){
#           di <- dmat[clustering==i,clustering==i,drop=FALSE]
#           oindexes <- (1:n)[clustering==i]
#           modei <- which.max(withindensp[clustering==i])
#           orderfrommode <- order(di[modei,])
#           distordereddens <- withindensp[clustering==i][orderfrommode]
#           ordereddens <- sort(withindensp[clustering==i],decreasing=TRUE)
#           dindexc[i] <- sum((distordereddens-ordereddens)^2)
#         } # for i
#         dindex <- sum(dindexc)/n
## output        
        if (!aggregateonly) 
            out <- list(n = n, cluster.number = cn, cluster.size = cluster.size, 
                min.cluster.size = min(cluster.size[1:cwn]), 
                noisen = noisen, diameter = diameter, average.distance = average.distance, 
                median.distance = median.distance, separation = separation, 
                average.toother = average.toother, separation.matrix = separation.matrix, 
                ave.between.matrix = ave.between.matrix, avebetween = average.between, 
                avewithin = average.within, n.between = nbetween, 
                n.within = nwithin, maxdiameter = max(diameter[1:cwn], 
                  na.rm = TRUE)/stan, minsep = (sepwithnoise * 
                  min(separation) +
                      (!sepwithnoise) * min(separation[1:cwn]))/stan, 
                withinss = within.cluster.ss,
#                        ave.within.cluster.ss = within.cluster.ss/(n - noisen),
                        clus.avg.silwidths = clus.avg.widths, 
                asw = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma, 
                dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within/average.between, 
                ch = ch, cwidegap = cwidegap, widestgap = widestgap, 
                corrected.rand = corrected.rand, 
                vi = vi, sindex = sindex,
# new
            svec=svec, psep=psep, stan=stan, nnk=nnk,mnnd=mnnd,pamc=pamc,pamcentroids=pamcentroids,
#            dindex=dindex,dsepi=dsepi, dratio=dratio,
#                        pdensity=pdensity,borderp=borderp,
            dindex=dindex,denscut=denscut,highdgap=highdgap,
#            dindexc=dindexc,denscutc=denscutc,
            npenalty=npenalty,dpenalty=dpenalty,
            withindensp=withindensp, densoc=densoc,
            pdistto=pdistto, pclosetomode=pclosetomode, distto=distto,
            percwdens=percwdens, percdensoc=percdensoc, 
                        parsimony=parsimony,cvnnd=cvnnd,cvnndc=cvnndc)
# end if !aggregateonly        
        else out <- list(n = n, cluster.number = cn, min.cluster.size = min(cluster.size[1:cwn]), 
            noisen = noisen, avebetween = average.between, 
            avewithin = average.within, maxdiameter = max(diameter[1:cwn], 
                na.rm = TRUE)/stan, minsep = (sepwithnoise * 
                min(separation) + (!sepwithnoise) * min(separation[1:cwn]))/stan, 
            withinss = within.cluster.ss, # was ave 
            asw = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma, 
            dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within/average.between, 
            ch = ch, widestgap = widestgap, 
            corrected.rand = corrected.rand, vi = vi, sindex = sindex, 
# new
            svec=svec, psep=psep, stan=stan, nnk=nnk,mnnd=mnnd,pamc=pamc,pamcebtroids=pamcentroids,
            dindex=dindex,denscut=denscut,highdgap=highdgap,
#            cpenalty=cpenalty,dpenalty=dpenalty,
#            npenalty=npenalty,
#            dsepi=dsepi, dratio=dratio,
#                        pdensity=pdensity,borderp=borderp,
                         parsimony=parsimony,cvnnd=cvnnd,cvnndc=cvnndc)
# end else !if aggregateonly        
    } # else (!if compareonly)   
    class(out) <- "cquality"
    out
}

# old densityindex
#         dsep <- despi <- dindex <- dratio <- pdensity <- borderp <- NULL
#         if (densityindex){
#           borderp <- matrix(TRUE,nrow=n,ncol=cwn)
#           pdensity <- rep(0,n)
#           for (i in 1:n){
#             sdi <- sort(dmat[i,],partial=dk+1)
#             pdensity[i] <- dk/(2*mean(sdi[2:(dk+1)]))
#             for (j in (1:cwn)[-clustering[i]])
#               borderp[i,j] <- any(dmat[i,clustering==j]<=sdi[dk+1])
#           } # for i
#           pdensity[pdensity==Inf] <- 2*max(pdensity[pdensity<Inf])
#           odensity <- mean(pdensity)
#           dsepi <- nni <- rep(1,cwn)
#           aidensity <- numeric(0)
# #          dsep <- matrix(1,nrow=cwn,ncol=cwn)          
# #           for (i in 1:(cwn-1)){
# #             for (j in (i+1):cwn)
# #               if (any((clustering==i & borderp[,j])|
# #                       ( clustering==j & borderp[,i])))
# #                 dsep[i,j] <- dsep[j,i] <-
# #                   min(c(1,
# #                     max(c(pdensity[clustering==i & borderp[,j]],
# #                         pdensity[clustering==j & borderp[,i]]))/
# #                     min(c(mean(pdensity[clustering==i]),
# #                         mean(pdensity[clustering==j])))))
# #           }
#           for (i in 1:cwn){
#             bdensity <- c(pdensity[clustering!=i & borderp[,i]],
#                           pdensity[clustering==i & apply(borderp,1,sum)>1])
#             idensity <- pdensity[clustering==i & apply(borderp,1,sum)==1]
#             nni[i] <- length(idensity)
#             aidensity <- c(aidensity,idensity)
#             if (nni[i]>0){
#               if (length(bdensity)>0)  
#                 d1 <- mean(bdensity)/mean(idensity)
#               else
#                 d1 <- 0
#               dsepi[i] <- min(c(1,d1))
#             if (length(bdensity)==0) dsepi[i] <- 0
#             } # nni>0
#           } # for i 
#           dratio <- odensity/mean(aidensity)
#           dindex1 <- 1-weighted.mean(dsepi,nni/sum(nni))
#           dindex2 <- 1-min(c(1,dratio))
#           nonoise <- n-noisen
#           doweight <- 0.5-abs((nonoise-sum(nni))/nonoise-0.5)
#           dindex <- (1-doweight)*dindex1+doweight*dindex2
# #          fdsep <- cwn-sum(dsep)
# #          dindex <- fdsep/sqrt(cwn)+(sqrt(maxk)-sqrt(cwn))*(1-dborder)
# #          if (standardisation!="none")
# #            dindex <- dindex/sqrt(maxk)




# largeisgood will transform ave within distance, within.cluster.ss, dindex,
# largest within gap, nndist, pam criterion by 1-x, so that large is good
# all will be bounded by 1 if stanbound, this also bounds
# pearsongamma by 0 from below.
summary.cquality <- function(object,stanbound=TRUE,largeisgood=TRUE,...){
  x <- object
  if (stanbound){
      x$mnnd <- min(1,x$mnnd)
      x$avewithin <- min(1,x$avewithin)
      x$pearsongamma <- max(0,x$pearsongamma)
      x$dindex <- min(1,x$dindex)
      x$denscut <- min(1,x$denscut)
      x$highdgap <- min(1,x$highdgap)
      x$asw <- max(0,x$asw)
      x$widestgap <- min(1,x$widestgap)
      if(!is.null(x$pamc))
        x$pamc <- min(1,x$pamc)
      x$maxdiameter <- min(1,x$maxdiameter)
      x$minsep <- min(1,x$minsep)
      x$sindex <- min(1,x$sindex)
      x$cvnnd <- min(1,x$cvnnd)
  }
  if (largeisgood){
      x$avewithin <- 1-x$avewithin
      x$mnnd <- 1-x$mnnd
      x$widestgap <- 1-x$widestgap
      x$denscut <- 1-x$denscut
      x$withinss <- 1-x$withinss
      if(!is.null(x$pamc))
        x$pamc <- 1-x$pamc
      x$maxdiameter <- 1-x$maxdiameter
      x$dindex <- 1-x$dindex
      x$highdgap <- 1-x$highdgap
      x$cvnnd <- 1-x$cvnnd
  }
  out <- list(avewithin=x$avewithin,nnk=x$nnk,mnnd=x$mnnd,
              asw=x$asw,
              widestgap=x$widestgap,sindex=x$sindex,
              pearsongamma=x$pearsongamma,entropy=x$entropy,pamc=x$pamc,
              withinss=x$withinss,
              dindex=x$dindex,denscut=x$denscut,highdgap=x$highdgap,
              parsimony=x$parsimony,maxdiameter=x$maxdiameter,
              minsep=x$minsep,cvnnd=x$cvnnd)
  class(out) <- "summary.cquality"
  out
}

print.summary.cquality <- function(x,...){
  cat("\n")  
  cat("average.within= ",x$avewithin,"\n")
  cat(x$nnk," nearest neighbour distance= ", x$mnnd,"\n")
  cat("coefficient of variation of ",x$nnk," nearest neighbour distance=",
      x$cvnnd,"\n")
  cat("max diameter= ",x$maxdiameter,"\n")
  cat("within-cluster widest gap= ",x$widestgap,"\n\n")
  cat("separation index= ",x$sindex,"\n")
  cat("min separation= ",x$minsep,"\n\n")
  cat("density index= ",x$dindex,"\n")
  cat("density cut= ",x$denscut,"\n")
  cat("density gap= ",x$highdgap,"\n")
  cat("average silhouette width= ",x$asw,"\n\n")
  cat("Pearson Gamma= ",x$pearsongamma,"\n")
  cat("pam criterion= ",x$pamc,"\n")
  cat("within ss= ",x$withinss,"\n\n")
  cat("entropy= ",x$entropy,"\n")
  cat("parsimony= ",x$parsimony,"\n")
}
  



# Similarity to "normal" or "uniform" distribution 
# x: data, clustering: clustering (has to be 1:k), noisecluster: last cluster
# is to be ignored, nnk: nearest neighbour k for uniform
# The uniform thing is misleading; this will indicate uniformity
# even on disconncted sets.
distrsimilarity <- function(x,clustering,noisecluster = FALSE,
                            distribution=c("normal","uniform"),nnk=2,
                           largeisgood=FALSE,messages=FALSE){
  cn <- max(clustering)
  clusteringf <- as.factor(clustering)
  clusteringl <- levels(clusteringf)
  cnn <- length(clusteringl)
  if (cn != cnn) {
        warning("clustering renumbered because maximum != number of clusters")
        for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
  }
  x <- as.matrix(x)
  nn <- n <- nrow(x)
  p <- ncol(x)
  k <- cnn
  if (noisecluster){
    nn <- n-sum(clustering==k)
    k <- k-1
#    x <- x[clustering<=k,]
  }
  xmahal <- xdknn <- rep(NA,n)
  kdnorm <- kdunif <- NA
  kdnormc <- kdunifc <- wkdnormc <- wkdunifc <- rep(NA,k)
  if ("normal" %in% distribution){
    for (i in 1:k){
#      print(k)
#      print(str(x))      
      cm <- colMeans(x[clustering==i,,drop=FALSE])
      ccov <- cov(x[clustering==i,,drop=FALSE])
#      cat("Cluster",i,"\n")
#      print(sum(clustering==i))
#      print(nn)
#      print(cm)
#      print(ccov)
      cmah <- try(mahalanobis(x[clustering==i,,drop=FALSE],cm,ccov),silent=TRUE)
      if (!is.null(attr(cmah,"class"))){
        cmah <- rep(0,sum(clustering==i))
        if (messages)
          print("cov-matrix not invertible")
      }
      if (any(is.na(cmah))){
        cmah <- rep(0,sum(clustering==i))
        if (messages)
          print("cov-matrix not invertible")
      }
#      print("after try")
      emahal <- ecdf(cmah)
      emahalm <- emahal(cmah)
      memahalm <- min(emahalm)
#      print(emahalm)
#      print(emahalm-pchisq(cmah,df=p))
#      print(emahalm-memahalm-pchisq(cmah,df=p))      
#      print(memahalm)
      kdnormc[i] <- max(abs(emahalm-pchisq(cmah,df=p)))
      kdnormc[i] <- max(c(abs(emahalm-memahalm-pchisq(cmah,df=p)),kdnormc[i]))
      wkdnormc[i] <- kdnormc[i]*sum(clustering==i)/nn
      xmahal[clustering==i] <- cmah
    } # for i
    kdnorm <- sum(wkdnormc)
    if (largeisgood) kdnorm <- 1-kdnorm
#    print(kdnormc)
#    print(wkdnormc)
  } # normal distribution: Mahalanobis distances vs. chisq
#  print("normal end")
  if ("uniform" %in% distribution){
    dmat <- as.matrix(dist(x))
    for (i in 1:k){
#      cat("Cluster",i,"\n")
      dmati <- dmat[clustering==i,clustering==i]
      nnki <- nnk
      if(sum(clustering==i)<nnk+1) nnki <- 1
      if(sum(clustering==i)>1){
        knnd <- numeric(0)
        for (j in 1:sum(clustering==i))
          knnd[j] <- sort(dmati[j,],partial=nnki+1)[nnki+1]
      }
      else
        knnd <- 0
      alphap <- (2 * pi^(p/2))/(p * gamma(p/2))
      if (mean(knnd^p)>0)
        lambdap <- nnki/(alphap * mean(knnd^p))
      else
        lambdap <- 1/alphap
      ekn <- ecdf(knnd^p)
      eknm <- ekn(knnd^p)
      meknm <- min(eknm)
#      print(alphap)
#      print(lambdap)
#      print(knnd^p)
#      print(eknm)
#      print(eknm-pgamma(knnd^p,shape=nnki,rate=lambdap*alphap))
      kdunifc[i] <- max(abs(eknm-pgamma(knnd^p,shape=nnki,rate=lambdap*alphap)))
      kdunifc[i] <- max(c(abs(eknm-
            meknm-pgamma(knnd^p,shape=nnki,rate=lambdap*alphap)),kdunifc[i]))
      wkdunifc[i] <- kdunifc[i]*sum(clustering==i)/nn
      xdknn[clustering==i] <- knnd
    }
    kdunif <- sum(wkdunifc)
    if (largeisgood) kdunif <- 1-kdunif
  }
  out <- list(kdnorm=kdnorm,kdunif=kdunif,kdnormc=kdnormc,kdunifc=kdunifc,
              xmahal=xmahal,xdknn=xdknn)
# kdnorm,kdunif: Kolmogorov distance (norm: Mahal to chisq,
# unif: knndist to Gamma, clusterwise, weighted by cluster size),
# kdnormc, kdunifc: clusterwise values, xmahal: Mahalanobis distances
# xdknn: distance to kth nearest neighbour within cluster
  out
}

# This computes k-centroids clustering to k random centroids (Serhat's version)
stupidkcentroids <- function(xdata, k, distances = inherits(xdata, "dist")){
  if (distances){
    cdist <- as.matrix(xdata)
  }
  else{
    cdist <- as.matrix(dist(xdata))
  }
  n <- ncol(cdist)
  kcent <- sample(n,k)
  clustering <- rep(0,n)
  clustering[kcent] <- 1:k
  topredict <- (1:n)[clustering==0]
  clustering[topredict] <- apply(cdist[topredict,kcent], 1, which.min)
  
  if (distances){
    centroids <- kcent
  }
  else{
    centroids <- xdata[kcent,]
  }
# print(str(centroids))
  out <- list(partition = clustering, 
              centroids = centroids,
              distances=distances)
  out
}


# stupidkcentroids <- function(d,k){
#   cdist <- as.matrix(d)
#   n <- ncol(cdist)
#   kcent <- sample(n,k)
#   clustering <- rep(0,n)
#   clustering[kcent] <- 1:k
#   topredict <- (1:n)[clustering==0]
#   clustering[topredict] <- apply(cdist[topredict,kcent], 1, which.min)
#   clustering
# }

stupidkcentroidsCBI <- function(dmatrix,k,distances=TRUE){
  c1 <- stupidkcentroids(dmatrix,k,distances=distances)
  partition <- c1$partition
  cl <- list()
  for (i in 1:k) cl[[i]] <- partition == i
  out <- list(result = c1, nc = max(c1$partition), clusterlist = cl, partition = partition, 
        clustermethod = "randomkcentroids")
#  print(str(out))
#  print(out$result["centroids"])
  out
}

# This computes a clustering starting from k centroids in which in each step
# the nearest neighbour of an existing cluster is added to that cluster.
stupidknn <- function(d,k){
  cdist <- as.matrix(d)
  n <- ncol(cdist)
  kcent <- sample(n,k)
  clustering <- rep(0,n)
  clustering[kcent] <- 1:k
  classified <- clustering>0
  repeat{
    topredict <- clustering==0
    cdistx <- cdist[topredict,classified,drop=FALSE]
#    print(str(cdistx))
    opt.ind <- which(cdistx==min(cdistx),arr.ind=TRUE)[1,,drop=FALSE]
    clustering[topredict][opt.ind[1,1]] <- clustering[classified][opt.ind[1,2]]
    classified <- clustering>0
    if(sum(classified)==n) break
  }
  clustering
}

stupidknnCBI <- function(dmatrix,k){
  c1 <- stupidknn(dmatrix,k)
  partition <- c1
  cl <- list()
  for (i in 1:k) cl[[i]] <- partition == i
  out <- list(result = c1, nc = max(c1), clusterlist = cl, partition = partition, 
        clustermethod = "randomkcentroids")
  out
}


# This computes a clustering starting from k centroids in which in each step
# the point that is nearest to the furthest point of an existing cluster is added to that cluster.
stupidkfn <- function(d,k){
  cdist <- as.matrix(d)
  n <- ncol(cdist)
  kcent <- sample(n,k)
  clustering <- rep(0,n)
  clustering[kcent] <- 1:k
  classified <- clustering>0
  repeat{
    cclustering  <-  clustering[clustering>0]
    topredict <- clustering==0
    cdistx <- cdist[topredict,classified,drop=FALSE]
#    print(str(cdistx))
    opt.ind <- minmaxd <- numeric(0)
    for (i in 1:k){
      cmk <- cdistx[,cclustering==i,drop=FALSE]
      maxd <- apply(cmk,1,max)
      minmaxd[i] <- min(maxd)
      opt.ind[i] <- which.min(maxd)
#      opt.ind2[i] <- which(cmk==minmaxd[i],arr.ind=TRUE)[1,,drop=FALSE]
    }
    kopt <- which.min(minmaxd)
    clustering[topredict][opt.ind[kopt]] <- kopt
    classified <- clustering>0
    if(sum(classified)==n) break
  }
  clustering
}

stupidkfnCBI <- function(dmatrix,k){
  c1 <- stupidkfn(dmatrix,k)
  partition <- c1
  cl <- list()
  for (i in 1:k) cl[[i]] <- partition == i
  out <- list(result = c1, nc = max(c1), clusterlist = cl, partition = partition, 
        clustermethod = "randomkcentroids")
  out
}


# This computes a clustering starting from k centroids in which in each step
# the point that has minimum average distiance to an existing cluster is added to that cluster.
stupidkaven <- function(d,k){
  cdist <- as.matrix(d)
  n <- ncol(cdist)
  kcent <- sample(n,k)
  clustering <- rep(0,n)
  clustering[kcent] <- 1:k
  classified <- clustering>0
  repeat{
    cclustering  <-  clustering[clustering>0]
    topredict <- clustering==0
    cdistx <- cdist[topredict,classified,drop=FALSE]
#    print(str(cdistx))
    opt.ind <- minmaxd <- numeric(0)
    for (i in 1:k){
      cmk <- cdistx[,cclustering==i,drop=FALSE]
      maxd <- apply(cmk,1,mean)
      minmaxd[i] <- min(maxd)
      opt.ind[i] <- which.min(maxd)
#      opt.ind2[i] <- which(cmk==minmaxd[i],arr.ind=TRUE)[1,,drop=FALSE]
    }
    kopt <- which.min(minmaxd)
    clustering[topredict][opt.ind[kopt]] <- kopt
    classified <- clustering>0
    if(sum(classified)==n) break
  }
  clustering
}

stupidkavenCBI <- function(dmatrix,k){
  c1 <- stupidkaven(dmatrix,k)
  partition <- c1
  cl <- list()
  for (i in 1:k) cl[[i]] <- partition == i
  out <- list(result = c1, nc = max(c1), clusterlist = cl, partition = partition, 
        clustermethod = "randomkcentroids")
  out
}


# This prepares the output of cqcluster.stats for some other
# ways of plotting and printing, plot.clustatsum, print.clustatsum
# cbmethod: clustering method to be used for stability if useboot
clustatsum <- function(datadist=NULL,clustering,noisecluster=FALSE,
                       datanp=NULL,npstats=FALSE, useboot=FALSE,
                              bootclassif=NULL,
                              bootmethod="nselectboot",
                              bootruns=25, cbmethod=NULL,methodpars=NULL,
                       distmethod=NULL, dnnk=2,
                       pamcrit=TRUE,...){
  out <- list()
  if (is.null(datadist))
    datadist <- dist(datanp)
  outcstat <- cqcluster.stats(datadist,clustering,pamcrit=pamcrit,...)
#  out$sum <- summary(out$cstat)
  outsum <- summary(outcstat)
  out <- list()
  out$avewithin <- outsum$avewithin
  out$mnnd <- outsum$mnnd
  out$cvnnd <- outsum$cvnnd
  out$maxdiameter <- outsum$maxdiameter
  out$widestgap <- outsum$widestgap
  out$sindex <- outsum$sindex
  out$minsep <- outsum$minsep
  out$asw <-  outsum$asw
  out$dindex <- outsum$dindex
  out$denscut <- outsum$denscut
  out$highdgap <- outsum$highdgap
  out$pearsongamma <- outsum$pearsongamma
  out$withinss <- outsum$withinss
  out$entropy <- outsum$entropy
  if (pamcrit)
    out$pamc <- outsum$pamc
  if (npstats){
    outd <- distrsimilarity(datanp,clustering,nnk=dnnk,
                            noisecluster=noisecluster,
                            largeisgood=TRUE)
    out$kdnorm <- outd$kdnorm
    out$kdunif <- outd$kdunif
  }
  if (useboot){
    if (distmethod){
      if (bootmethod=="nselectboot"){
#        print(cbmethod)
        out$boot <- do.call(nselectboot,c(list(data=datadist,B=bootruns,
                                          distances=TRUE,
                              clustermethod=get(cbmethod),
                              classification=bootclassif,
                              krange=max(clustering),largeisgood=TRUE),
                              methodpars))$stabk[max(clustering)]
#        cat(cbmethod,"nselectboot end\n")
      }
      else
        out$boot <- do.call(prediction.strength,c(list(xdata=datadist,
                                                  M=bootruns,distances=TRUE,
                              clustermethod=get(cbmethod),
                              classification=bootclassif,
                              Gmin=max(clustering),
                              Gmax=max(clustering)),methodpars))$mean.pred[max(clustering)]
    }
    else{
      if (bootmethod=="nselectboot"){
#        print(cbmethod)
        out$boot <- do.call(nselectboot,c(list(data=datanp,B=bootruns,
                                          distances=FALSE,
                              clustermethod=get(cbmethod),
                              classification=bootclassif,
                              krange=max(clustering),largeisgood=TRUE),
                              methodpars))$stabk[max(clustering)]
#        cat(cbmethod,"nselectboot end\n")
      }
      else
        out$boot <- do.call(prediction.strength,c(list(xdata=datanp,M=bootruns,
                                                  distances=FALSE,
                              clustermethod=get(cbmethod),
                              classification=bootclassif,
                              Gmin=max(clustering),
                              Gmax=max(clustering)),methodpars))$mean.pred[max(clustering)]
    }
  } # if useboot
  out
}

    
# This applies stupidkcentroids, stupidknn, stupidkfn and stupidkaven
# lots of times to data.
randomclustersim <- function(datadist,datanp=NULL,npstats=FALSE,useboot=FALSE,
                              bootmethod="nselectboot",
                              bootruns=25, 
                      G,nnruns=100,kmruns=100,fnruns=100,avenruns=100,
                      nnk=4,dnnk=2,
                      pamcrit=TRUE, 
                      multicore=FALSE,cores=detectCores()-1,monitor=TRUE){
  nnsim <- function(k,g){
    if(monitor) cat(g," clusters; nn run ",k,"\n")
    solution <- data.frame(NA)
    names(solution) <- "avewithin"
    gcl <- stupidknn(datadist,g)
#    print(gcl)
    grs <- cqcluster.stats(datadist,gcl,nnk=nnk,pamcrit=pamcrit)
    grss <- summary(grs)
    solution$avewithin <- grss$avewithin
    solution$mnnd <- grss$mnnd
    solution$cvnnd <- grss$cvnnd
    solution$maxdiameter <- grss$maxdiameter
    solution$widestgap <- grss$widestgap
    solution$sindex <- grss$sindex
    solution$minsep <- grss$minsep
    solution$asw <-  grss$asw
    solution$dindex <- grss$dindex
    solution$denscut <- grss$denscut
    solution$highdgap <- grss$highdgap
    solution$pearsongamma <- grss$pearsongamma
    solution$withinss <- grss$withinss
    solution$entropy <- grss$entropy
    if (pamcrit)
      solution$pamc <- grss$pamc
    if (npstats){
      geysdist <- distrsimilarity(datanp,gcl,nnk=dnnk,largeisgood=TRUE)
      solution$kdnorm <- geysdist$kdnorm
      solution$kdunif <- geysdist$kdunif
    }
    if (useboot){
      if (bootmethod=="nselectboot")
        solution$boot <- nselectboot(datadist,B=bootruns,distances=TRUE,
                               clustermethod=stupidknnCBI,
                               classification="knn",
                               krange=g,largeisgood=TRUE)$stabk[g]
      else
        out$boot <- prediction.strength(datadist,M=bootruns,distances=TRUE,
                              clustermethod=stupidknnCBI,
                              classification="knn",
                              Gmin=g,
                              Gmax=g)$mean.pred[g]
    }
    solution
  }
  
  fnsim <- function(k,g){
    if(monitor) cat(g," clusters; fn run ",k,"\n")
#    print("start fnsim")
    solution <- data.frame(NA)
    names(solution) <- "avewithin"
    gcl <- stupidkfn(datadist,g)
#    print(gcl)
    grs <- cqcluster.stats(datadist,gcl,nnk=nnk,pamcrit=pamcrit)
    grss <- summary(grs)
    solution$avewithin <- grss$avewithin
    solution$mnnd <- grss$mnnd
    solution$cvnnd <- grss$cvnnd
    solution$maxdiameter <- grss$maxdiameter
    solution$widestgap <- grss$widestgap
    solution$sindex <- grss$sindex
    solution$minsep <- grss$minsep
    solution$asw <-  grss$asw
    solution$dindex <- grss$dindex
    solution$denscut <- grss$denscut
    solution$highdgap <- grss$highdgap
    solution$pearsongamma <- grss$pearsongamma
    solution$withinss <- grss$withinss
    solution$entropy <- grss$entropy
    if (pamcrit)
      solution$pamc <- grss$pamc
    if (npstats){
      geysdist <- distrsimilarity(datanp,gcl,nnk=dnnk,largeisgood=TRUE)
      solution$kdnorm <- geysdist$kdnorm
      solution$kdunif <- geysdist$kdunif
    }
    if (useboot){
#      print("boot")
      if (bootmethod=="nselectboot")
        solution$boot <- nselectboot(datadist,B=bootruns,distances=TRUE,
                               clustermethod=stupidkfnCBI,
                               classification="fn",
                               krange=g,largeisgood=TRUE)$stabk[g]
      else
        out$boot <- prediction.strength(datadist,M=bootruns,distances=TRUE,
                              clustermethod=stupidkfnCBI,
                              classification="fn",
                              Gmin=g,
                              Gmax=g)$mean.pred[g]
    }
    print("end fnsim")
    solution
  }
  
  avensim <- function(k,g){
    if(monitor) cat(g," clusters; aven run ",k,"\n")
    solution <- data.frame(NA)
    names(solution) <- "avewithin"
    gcl <- stupidkaven(datadist,g)
#    print(gcl)
    grs <- cqcluster.stats(datadist,gcl,nnk=nnk,pamcrit=pamcrit)
    grss <- summary(grs)
    solution$avewithin <- grss$avewithin
    solution$mnnd <- grss$mnnd
    solution$cvnnd <- grss$cvnnd
    solution$maxdiameter <- grss$maxdiameter
    solution$widestgap <- grss$widestgap
    solution$sindex <- grss$sindex
    solution$minsep <- grss$minsep
    solution$asw <-  grss$asw
    solution$dindex <- grss$dindex
    solution$denscut <- grss$denscut
    solution$highdgap <- grss$highdgap
    solution$pearsongamma <- grss$pearsongamma
    solution$withinss <- grss$withinss
    solution$entropy <- grss$entropy
    if (pamcrit)
      solution$pamc <- grss$pamc
    if (npstats){
      geysdist <- distrsimilarity(datanp,gcl,nnk=dnnk,largeisgood=TRUE)
      solution$kdnorm <- geysdist$kdnorm
      solution$kdunif <- geysdist$kdunif
    }
    if (useboot){
      if (bootmethod=="nselectboot")
        solution$boot <- nselectboot(datadist,B=bootruns,distances=TRUE,
                               clustermethod=stupidkavenCBI,
                               classification="averagedist",
                               krange=g,largeisgood=TRUE)$stabk[g]
      else
        out$boot <- prediction.strength(datadist,M=bootruns,distances=TRUE,
                              clustermethod=stupidkavenCBI,
                              classification="averagedist",
                              Gmin=g,
                              Gmax=g)$mean.pred[g]
    }
     solution
  }
      
  kmsim <- function(k,g){
    if(monitor) cat(g," clusters; km run ",k,"\n")
    solution <- data.frame(NA)
    names(solution) <- "avewithin"
    gcl <- stupidkcentroids(datadist,g)$partition
#    print(gcl)
    grs <- cqcluster.stats(datadist,gcl,nnk=nnk,pamcrit=pamcrit)
    grss <- summary(grs)
    solution$avewithin <- grss$avewithin
    solution$mnnd <- grss$mnnd
    solution$cvnnd <- grss$cvnnd
    solution$maxdiameter <- grss$maxdiameter
    solution$widestgap <- grss$widestgap
    solution$sindex <- grss$sindex
    solution$minsep <- grss$minsep
    solution$asw <-  grss$asw
    solution$dindex <- grss$dindex
    solution$denscut <- grss$denscut
    solution$highdgap <- grss$highdgap
    solution$pearsongamma <- grss$pearsongamma
    solution$withinss <- grss$withinss
    solution$entropy <- grss$entropy
    if (pamcrit)
      solution$pamc <- grss$pamc
    if (npstats){
      geysdist <- distrsimilarity(datanp,gcl,nnk=dnnk,largeisgood=TRUE)
      solution$kdnorm <- geysdist$kdnorm
      solution$kdunif <- geysdist$kdunif
    }
    if (useboot){
      if (bootmethod=="nselectboot")
        solution$boot <- nselectboot(datadist,B=bootruns,distances=TRUE,
                               clustermethod=stupidkcentroidsCBI,
                               classification="centroid",
                               centroidname="centroids",
                               krange=g,largeisgood=TRUE)$stabk[g]
      else
        out$boot <- prediction.strength(datadist,M=bootruns,distances=TRUE,
                              clustermethod=stupidkcentroidsCBI,
                              classification="centroid",
                              centroidname="centroids",
                              Gmin=g,
                              Gmax=g)$mean.pred[g]
    }
     solution
  }
      
  out <- list()
  out$nn <- out$fn <- out$aven <- out$km <- list()
  out$nnruns <- nnruns
  out$fnruns <- fnruns
  out$avenruns <- avenruns
  out$kmruns <- kmruns
  for (g in G){
#    out$nn[[g]] <- out$km[[g]] <- data.frame()
    if (multicore){
      if (nnruns>0)
        out$nn[[g]] <- do.call(rbind,
                               mclapply(1:nnruns,nnsim,g=g,mc.cores=cores))
      if (fnruns>0)
        out$fn[[g]] <- do.call(rbind,
                               mclapply(1:fnruns,fnsim,g=g,mc.cores=cores))
      if (avenruns>0)
        out$aven[[g]] <- do.call(rbind,
                               mclapply(1:avenruns,avensim,g=g,mc.cores=cores))
      if(kmruns>0)
        out$km[[g]] <- do.call(rbind,
                               mclapply(1:kmruns,kmsim,g=g,mc.cores=cores))
    }
    else{
      if (nnruns>0)
        out$nn[[g]] <- do.call(rbind,lapply(1:nnruns,nnsim,g=g))
      if (fnruns>0)
        out$fn[[g]] <- do.call(rbind,lapply(1:fnruns,fnsim,g=g))
      if (avenruns>0)
        out$aven[[g]] <- do.call(rbind,lapply(1:avenruns,avensim,g=g))
      if (kmruns>0)
        out$km[[g]] <- do.call(rbind,lapply(1:kmruns,kmsim,g=g))
    }
  }
  out
}



######################################################
                                        # CBI-version
######################################################


# Standardise output of clustatsum by results of randomclustersim,
# to be used within clusterbenchstats
# percentage: quantile calibration method rather than sd
# useallmethods: use not only stupid clusterings
# useallg: Use all studpid clusterings with all numbers of clusters
# for calibration
# othernc: list of pairs of clustering method number and
# number of clusters larger than max(G).
# Only if useallg=TRUE, these are standardised as well.
# will add density mode index, 0.75*dindex+0.25*highdgap
cgrestandard <- function(clusum,clusim,G,percentage=FALSE,
                               useallmethods=FALSE,
                             useallg=FALSE, othernc=list()){

  qstan <- function(x,svec,nmethods){
#    print(x)
    ls <- length(svec[!is.na(svec)])
    out <- rep(NA,length(x))
    for (i in (1:nmethods)[!is.na(x)])
      out[i] <- (sum(x[i]>svec,na.rm=TRUE)+1)/(ls+1)
    out
  }
               
  sdstan <- function(x,svec){
    sm <- mean(svec,na.rm=TRUE)
    ssd <- sd(svec,na.rm=TRUE)
    (x-sm)/ssd
  }
    
  nmethods <- length(clusum$method)               
  statistics <- c("avewithin","mnnd","cvnnd","maxdiameter",
                  "widestgap","sindex","minsep","asw","dindex","denscut","highdgap",
                  "pearsongamma","withinss","entropy")
  if(!is.null(clusum[[1]][[min(G)]]$pamc))
     statistics <- c(statistics,"pamc")
  if(!is.null(clusum[[1]][[min(G)]]$kdnorm))
    statistics <- c(statistics,"kdnorm","kdunif")
  if(!is.null(clusum[[1]][[min(G)]]$boot))
    statistics <- c(statistics,"boot")

  out <- clusum
  lon <- length(othernc)
  for (st in statistics){
    if (useallg){
      svec <- numeric(0)
      sumvec <- list()
      for (g in G){
        sumvec[[g]] <- rep(NA,nmethods)
        for (i in 1:nmethods){
#          print(st)
#          print(g)
#          print(i)
          if (!identical(clusum[[i]][[g]],NA))
            sumvec[[g]][i] <- clusum[[i]][[g]][[st,exact=FALSE]]
          }
        if (clusim$kmruns>0)
          svec <- c(svec, clusim$km[[g]][[st,exact=FALSE]])
        if (clusim$nnruns>0)
          svec <- c(svec, clusim$nn[[g]][[st,exact=FALSE]])
        if (clusim$fnruns>0)
          svec <- c(svec, clusim$fn[[g]][[st,exact=FALSE]])
        if (clusim$avenruns>0)
          svec <- c(svec, clusim$aven[[g]][[st,exact=FALSE]])
        if (useallmethods)
          svec <- c(svec,sumvec[[g]])
      }
      for (g in G){
#        if(percentage){
#          print("cgrestandard")
#          print(str(sumvec))
#          print(str(svec))
#          print(nmethods)
#          rvec <- qstan(sumvec[[g]],svec,nmethods)
#        }
        if(percentage)
          rvec <- qstan(sumvec[[g]],svec,nmethods)
        else
          rvec <- sdstan(sumvec[[g]],svec)
        for (i in 1:nmethods)
          if (!identical(clusum[[i]][[g]],NA))
            out[[i]][[g]][st] <- rvec[i]
      }
      if (lon>0)
        for(j in 1:lon){
          i <- othernc[[j]][1]
          g <- othernc[[j]][2]
          if (g>1){
            if (percentage){
              ls <- length(svec[!is.na(svec)])
              out[[i]][[g]][st] <-
                (sum(clusum[[i]][[g]][[st,exact=FALSE]]>
                     svec,na.rm=TRUE)+1)/(ls+1)
            }
            else
              out[[i]][[g]][st] <-
                sdstan(clusum[[i]][[g]][[st,exact=FALSE]],svec)
          }
        }  
    } # if useallg
    else{
      for (g in G){
        svec <- numeric(0)
        sumvec <- rep(NA,nmethods)
#        browser()
        for (i in 1:nmethods){
#          print(st)
#          print(g)
#          print(i)
          if (!identical(clusum[[i]][[g]],NA))
            sumvec[i] <- clusum[[i]][[g]][[st,exact=FALSE]]
        }
        if (clusim$kmruns>0)
          svec <- c(svec, clusim$km[[g]][[st,exact=FALSE]])
        if (clusim$nnrun>0)
          svec <- c(svec, clusim$nn[[g]][[st,exact=FALSE]])
        if (clusim$fnruns>0)
          svec <- c(svec, clusim$fn[[g]][[st,exact=FALSE]])
        if (clusim$avenrun>0)
          svec <- c(svec, clusim$aven[[g]][[st,exact=FALSE]])
        if (useallmethods)
          svec <- c(svec,sumvec)
#        if(percentage){
#          print("cgres. not useallg")
#          print(str(sumvec))
#          print(str(svec))
#          print(nmethods)
#          rvec <- qstan(sumvec,svec,nmethods)
#        }
        if(percentage)
          rvec <- qstan(sumvec,svec,nmethods)
        else
          rvec <- sdstan(sumvec,svec)
        for (i in 1:nmethods)
          if (!identical(clusum[[i]][[g]],NA))
            out[[i]][[g]][st] <- rvec[i]
      } # for g
      if (lon>0)
        for(j in 1:lon){
          i <- othernc[[j]][1]
          g <- othernc[[j]][2]
          if (g>1)
              out[[i]][[g]][st] <- NA
        } 
    } # else (!useallg)
  } # for st
# Computation of dmode
  statistics <- c(statistics,"dmode")
#  print(nmethods)
  for (i in 1:nmethods){
    for (g in G){
      if (!identical(out[[i]][[g]],NA)){
#        print(i)
#        print(g)
#        print(out[[i]][[g]])
        out[[i]][[g]]$dmode <- 0.75*out[[i]][[g]]$dindex+0.25*out[[i]][[g]]$highdgap
      }
    }
  } # for i
  if (lon>0)
    for(j in 1:lon){
      i <- othernc[[j]][1]
      g <- othernc[[j]][2]
      if (g>1)
        out[[i]][[g]]$dmode <- 0.75*out[[i]][[g]]$dindex+0.25*out[[i]][[g]]$highdgap
    }
  class(out) <- "valstat"
  out
}
      




# clustermethod: string vector with CBI functions
# clustermethodpars: list of length length(clustermethod)
# that specifies parameters for all clustermethods
# distmethod: Will the method interpret input data as distances?
# ncinput: Does the method run with number of clusters as input?

# Note that different parameter choices of, e.g., dbscan, can be
# compared providing dbscanCBI several times as clustermethod
# with different clustermethodpars.
cluster.magazine <- function(data,G,diss = inherits(data, "dist"),
                             scaling=TRUE, clustermethod,
                             distmethod=rep(TRUE,length(clustermethod)),
                             ncinput=rep(TRUE,length(clustermethod)),
                             clustermethodpars,
#                             nstart=10,iter.max=100,
#                             nnk=0, emModelNames=NULL, mdsmethod="classical",
#                             mdsdim=4,
#                             summary.out=FALSE,points.out=FALSE,
#                             usepam=TRUE, samples=100,
#                             usepdf=TRUE, graphtype="pairs",
#                             lambda=c(0.1,0.01,0.001), 
                             trace=TRUE){
# scaling kmeans, linkage methods, pam/clara
# nstart kmeans    
# iter.max kmeans
# nnk mclust
# emModelNames mclust
# mdsmethod distnoisemclust
# mdsdim distnoisemclust 
# summary.out mclust
# points.out distnoisemclust
# samples clara
# graphtype pdfCluster  
# lambda pdfCluster, only effect if graphtype="pairs"    
  out <- list()
  out$output <- list()
  out$clustering <- list()
  out$noise <- list()
  out$othernc <- list()
  otherncc <- 1
#  out$method <- character(0)
#  else sdata <- data
  if(diss){
    if(all(distmethod))
      datadist <- as.dist(data)
    else
      stop("If the input is distance data, all methods must be distmethods.")
  }
  else{
    if (!identical(scaling, FALSE)) 
      sdata <- scale(data, center = TRUE, scale = scaling)
    else
      sdata <- data
    if (any(distmethod))
      datadist <- dist(sdata)
  }
  lmo <- length(clustermethod)
  for (i in 1:lmo){
    if (trace) print(clustermethod[i])
    out$output[[i]] <- out$clustering[[i]] <- out$noise[[i]] <- as.list(rep(NA,max(G)))
    clusterfunction <- get(clustermethod[i])
    if (ncinput[[i]]){
      for (g in G){
        cpars <- clustermethodpars[[i]]
        cpars$k <- g
        if (distmethod[i]){
          if (clustermethod[[i]] %in% c("disthclustCBI","distnoisemclustCBI",
                                        "disttrimkmeansCBI"))
            cpars$dmatrix <- datadist
          else
            cpars$data <- datadist
        }
        else
          cpars$data <- sdata
        out$output[[i]][[g]] <- do.call(clusterfunction, cpars)
        out$output[[i]][[g]]$clusterlist <- NULL
        out$clustering[[i]][[g]] <- out$output[[i]][[g]]$partition
        if (is.null(out$output[[i]][[g]]$nccl))
          out$noise[[i]][[g]] <- FALSE
        else{
          if (out$output[[i]][[g]]$nc>out$output[[i]][[g]]$nccl)
            out$noise[[i]][[g]] <- TRUE
          else
            out$noise[[i]][[g]] <- FALSE
        }
      }
    } # if ncinput
    else{
      cpars <- clustermethodpars[[i]]
      if (distmethod[i]){
        if (clustermethod[[i]] %in% c("disthclustCBI","distnoisemclustCBI",
                                      "disttrimkmeansCBI"))
          cpars$dmatrix <- datadist
        else
          cpars$data <- datadist
      }
      else
        cpars$data <- sdata
      coutput <- do.call(clusterfunction, cpars)
      coutput$clusterlist <- NULL
      if (is.null(coutput$nccl)){
        cnum <- coutput$nc
        out$noise[[i]][[cnum]] <- FALSE
      }
      else{
        cnum <- max(coutput$nccl,1)      
        if (coutput$nc>cnum)
          out$noise[[i]][[cnum]] <- TRUE
        else
          out$noise[[i]][[cnum]] <- FALSE
      }
      if (cnum>max(G) | cnum<min(G)){
        out$othernc[[otherncc]] <- c(i,cnum)
        otherncc <- otherncc+1
      }
      out$output[[i]][[cnum]] <- do.call(clusterfunction, cpars)
      out$clustering[[i]][[cnum]] <- out$output[[i]][[cnum]]$partition
    }
  }
  out  
}
# othernc will be a list of pairs of method number and number of clusters
# whenever the number of clusters is larger max(G)



# useallmethods, useallg: These are for restandardisation;
# useallmethods means that all method results are used, not only stupid
# clusterings; useallg: standardisation uses all values of G everywhere
# (not only for the same value)
# ncinput and others: see cluster.magazine.
# if ncinput=FALSE for some methods and number of clusters is
# estimated > max(G), see cgrestandard for the handling.
# This requires useallg=TRUE, otherwise these results are ignored
# for standardisation.
#
# useboot will involve a stability resampling method.
# bootclassif is a vector of classification methods to be used with
# the different clustering methods.
# bootmethod is either "nselectboot" or "prediction.strength"
# bootruns is the number of bootstrap replicates, all only
# active if useboot=TRUE.

# Help page: clustermethodpars last entry must be specified
# Mentions output list member statistics that doesn't exist
# (it's stat$statistics).
# Can't use methods like dbscan that don't take fixed nc with useboot!
clusterbenchstats <- function(data,G,diss = inherits(data, "dist"),
                                  scaling=TRUE, clustermethod,
                                  methodnames=clustermethod,
                              distmethod=rep(TRUE,length(clustermethod)),
                              ncinput=rep(TRUE,length(clustermethod)),
                              clustermethodpars,
#                              nstart=10,iter.max=100,
#                              nnk=0, emModelNames=NULL, mdsmethod="classical",
#                              mdsdim=4,summary.out=FALSE,points.out=FALSE,
#                              usepam=TRUE, samples=100,
                              npstats=FALSE,
                              useboot=FALSE,
                              bootclassif=NULL,
                              bootmethod="nselectboot",
                              bootruns=25,
#                              usepdf=TRUE, graphtype="pairs",
#                              lambda=c(0.1,0.01,0.001), 
                              trace=TRUE,
                              pamcrit=TRUE,snnk=2,
                              dnnk=2,
#                              noisecluster=FALSE,
                              nnruns=100,kmruns=100,fnruns=100,avenruns=100,
                              multicore=FALSE,cores=detectCores()-1,
                              useallmethods=TRUE,
                              useallg=FALSE,...){

  comsum <- function(i){
    if (trace)
      cat("comsum ",i,"\n")
    statl <- as.list(rep(NA,max(G)))
    for (g in G)
      if(!identical(out$cm$clustering[[i]][[g]],NA))
        statl[[g]] <- clustatsum(datadist,out$cm$clustering[[i]][[g]],
                                         noisecluster=out$cm$noise[[i]][[g]],
                                         datanp=data,npstats=npstats,
                                 useboot=useboot,
                                 bootclassif=bootclassif[i],
                                 bootmethod=bootmethod,
                                 bootruns=bootruns, cbmethod=clustermethod[i],
                                 methodpars=clustermethodpars[[i]],
                                 distmethod=distmethod[i],
                                 nnk=snnk,pamcrit=pamcrit,...)
    lon <- length(out$cm$othernc)
    if (lon>0)
      for (j in 1:lon)
        if(out$cm$othernc[[j]][1]==i){
          g <- out$cm$othernc[[j]][2]
          if (g>1)
            statl[[g]] <-
              clustatsum(datadist,out$cm$clustering[[i]][[g]],
                                         noisecluster=out$cm$noise[[i]][[g]],
                                         datanp=data,npstats=npstats,
                                 useboot=useboot,
                                 bootclassif=bootclassif[i],
                                 bootmethod=bootmethod,
                                 cbmethod=clustermethod[i],
                                 methodpars=clustermethodpars[[i]],
                                 bootruns=bootruns, 
                       nnk=snnk,pamcrit=pamcrit,...)
        }        
    statl
  }
   
  
  if(diss){
    datadist <- as.dist(data)
    data <- datadist
    datanp <- NULL
  }
  else{
    datanp <- data
    datadist <- dist(data)
  }
  out <- list()
  out$cm <- cluster.magazine(data,G=G,diss = diss,
                             scaling=scaling,
                             clustermethod=clustermethod,
                             distmethod=distmethod,
                             ncinput=ncinput,
                             clustermethodpars=clustermethodpars,
#                             nstart=nstart,iter.max=iter.max,
#                              nnk=nnk, emModelNames=emModelNames,
#                              mdsmethod=mdsmethod,
#                              mdsdim=mdsdim,summary.out=summary.out,
#                              points.out=points.out,
#                              usepam=usepam, samples=samples,
#                              usepdf=usepdf, graphtype=graphtype,
#                              lambda=lambda, 
                             trace=trace)
  nmethods <- length(clustermethod)
#  browser()
  if (trace)
    print("Computation of validity statistics")
  if (multicore)
    out$stat <- mclapply(1:nmethods,comsum,mc.cores=cores)
  else
    out$stat <- lapply(1:nmethods,comsum)
#  print(str(out$stat))
  out$stat$maxG <- max(G)
  out$stat$minG <- min(G)
  out$stat$method <-  clustermethod
  out$stat$name <-  methodnames
  class(out$stat) <- "valstat"
#  browser()
  if (trace)
    print("Simulation")
  out$sim <- randomclustersim(datadist,datanp,npstats=npstats,
                                 useboot=useboot,
                                 bootmethod=bootmethod,
                                 bootruns=bootruns, G=G,
                       nnruns=nnruns,kmruns=kmruns,fnruns=fnruns,
                       avenruns=avenruns,
                       nnk=snnk,dnnk=dnnk,pamcrit=pamcrit,
                       multicore=multicore,
                       cores=cores,monitor=trace)
  if (trace)
    print("Simulation quantile re-standardisation")
#  browser()
  out$qstat <- cgrestandard(out$stat,out$sim,G,percentage=TRUE,
                                useallmethods=useallmethods,useallg=useallg,
                                out$cm$othernc)
  if (trace)
    print("Simulation sd re-standardisation")
  out$sstat <- cgrestandard(out$stat,out$sim,G,percentage=FALSE,
                            useallmethods=useallmethods,useallg=useallg,
                                out$cm$othernc)
  ostatistics <- c("avewithin","mnnd","cvnnd","maxdiameter",
                  "widestgap","sindex","minsep","asw","dindex","denscut","highdgap",
                  "pearsongamma","withinss","entropy")
  if(pamcrit)
     ostatistics <- c(ostatistics,"pamc")
  if(npstats)
    ostatistics <- c(ostatistics,"kdnorm","kdunif")
  if(useboot)
    ostatistics <- c(ostatistics,"boot")
  ostatistics <- c(ostatistics,"dmode")
  out$stat$statistics <- out$qstat$statistics <- out$sstat$statistics <- ostatistics
  class(out) <- "clusterbenchstats"
  out
}
#> str(ds,max.level=1)
#List of 6
# $ cm        :List of 4
# $ stat      :List of 11
#  ..- attr(*, "class")= chr "clustatsum"
# $ sim       :List of 4
# $ qstat     :List of 11
#  ..- attr(*, "class")= chr "clustatsum"
# $ sstat     :List of 11
#  ..- attr(*, "class")= chr "clustatsum"
# $ statistics: chr [1:17] "avewithin" "mnnd" "cvnnd" "maxdiameter" ...

# out$cm has results of cluster.magazine
# > str(ds$cm,max.level=1)
#List of 4
# $ output    :List of 8
# $ clustering:List of 8
# $ noise     :List of 8
# $ othernc   :List of 1
# out$cm$output[[i]][[j]]: method i number of clusters j, CBI method output
# out$cm$clustering[[i]][[j]]: method i number of clusters j, clustering
# out$cm$noise[[i]][[j]]: method i number of clusters j is there noise?
# out$cm$othernc: numbers of clusters for methods that estimate this
# (list with methodsnumber, number of clusters) outside the specified range.

# out$stat has validity statistics
# out$stat[[i]][[j]]: method i number of clusters j, all statistics 
# out$stat$method has methods, out$stat$name has method names to be used for
# listing and plotting
# also out$stat$maxG and minG

# out$sim$km and out$sim$nn[[j]] have the stupid clustering results for j clusters
# There's also out$sim$kmruns and nnruns, number of runs

# out$qstat and out$sstat are organised like $stat, with quantile and
# stdev-standardised statistics


print.clusterbenchstats <- function(x,...){
  cat("Output object of clusterbenchstats.\n")
  cat("Clustering methods: ",x$stat$name," \n")
  cat("Cluster validation statistics: ",x$stat$statistics,"\n") 
  cat("Numbers of clusters minimum: ",x$stat$minG," maximum: ",
      x$stat$maxG,"\n")
  cat("Output components are cm, stat, sim, qstat, sstat.")
  cat("stat, qstat, and sstat are valstat-objects.")
  cat("Use plot.valstat and print.valstat on these to get more information.\n")
}


# include.othernc will overwrite xlim default.
# It should be clusterbenchoutput$cm$othernc otherwise. 
plot.valstat <- function(x,simobject=NULL,statistic="sindex",
                            xlim=NULL,ylim=c(0,1),
                            nmethods=length(x)-5,
                            col=1:nmethods,cex=1,pch=c("c","f","a","n"),
                            simcol=rep(grey(0.7),4),
                         shift=c(-0.1,-1/3,1/3,0.1),include.othernc=NULL,...){
  ion <- !is.null(include.othernc)
  if (ion) ion <- length(include.othernc)>0
  othernc <- numeric(0)
  q <- x$minG:x$maxG
  if (ion){
    for (i in 1:length(include.othernc))
      othernc[i] <- include.othernc[[i]][2]
    q <- sort(c(othernc,q))
  }
  if (is.null(xlim))
    xlim <- c(min(q)-0.5,max(q)+0.5)  
  plot(1,type="n",xlim=xlim,xlab="Number of clusters",
       ylim=ylim,ylab=statistic,xaxt="n")
  axis(1,at=q)
  if(!is.null(simobject)){
    for(g in x$minG:x$maxG){
      if(simobject$kmruns>0)
        points(rep(g+shift[1],simobject$kmruns),
               unlist(simobject$km[[g]][statistic]),
               pch=pch[1],col=simcol[1],cex=cex)
      if(simobject$nnruns>0)
        points(rep(g+shift[4],simobject$nnruns),
               unlist(simobject$nn[[g]][statistic]),
               pch=pch[4],col=simcol[4], cex=cex)
      if(simobject$fnruns>0)
        points(rep(g+shift[2],simobject$fnruns),
               unlist(simobject$fn[[g]][statistic]),
               pch=pch[2],col=simcol[2], cex=cex)
      if(simobject$avenruns>0)
        points(rep(g+shift[3],simobject$avenruns),
               unlist(simobject$aven[[g]][statistic]),
               pch=pch[3],col=simcol[3], cex=cex)
    }
  }    
  for(i in 1:nmethods)
    for(g in x$minG:x$maxG){
#      print(i)
#      print(g)
#      print(x[[i]][[g]])
      if(!identical(x[[i]][[g]],NA)){
#        points(g,x[[i]][[g]]$stats[statistic],pch=pch,col=col[i])
        text(g,x[[i]][[g]][statistic],labels=x$name[i],
             col=col[i],cex=cex)
      }
    }
    if (ion)
      for(h in 1:length(include.othernc)){
        mn <- include.othernc[[h]][1]
        ch <- include.othernc[[h]][2]
        if (ch>1)
          text(ch,x[[mn]][[ch]][statistic],labels=x$name[mn],
             col=col[mn],cex=cex)
      }        
}

print.valstat <- function(x,statistics=x$statistics,
                          nmethods=length(x)-5,aggregate=FALSE,
                          weights=NULL,digits=2,
                          include.othernc=NULL,...){

  q <- x$minG:x$maxG
  ion <- !is.null(include.othernc)
  if (ion) ion <- length(include.othernc)>0
  othernc <- numeric(0)
  if (ion){
    for (i in 1:length(include.othernc))
      othernc[i] <- include.othernc[[i]][2]
    q <- sort(c(othernc,q))
    q <- q[q>1]
  }
  lq <- length(q)
  if (aggregate){
    aggmatrix <- matrix(NA,nrow=nmethods,ncol=lq)
    for(j in 1:nmethods){
      for(i in 1:lq)
        if(length(x[[j]])>=q[i]){
          if(!identical(x[[j]][[q[i]]],NA)){
#            print(as.vector(as.matrix(x[[j]][[q[i]]])))
            x[[j]][[q[i]]]$aggregate <- weighted.mean(unlist(x[[j]][[q[i]]]),weights)
          }
        }
    }
    statistics <- c(statistics,"aggregate")
  } # if aggregate
  printobject <- list()
  l <- 1
  for(statistic in statistics){
    cat(statistic,"\n")
    printobject[[l]] <- data.frame(x$name)
    q <- x$minG:x$maxG
    ion <- !is.null(include.othernc)
    if (ion) ion <- length(include.othernc)>0
    othernc <- numeric(0)
    if (ion){
      for (i in 1:length(include.othernc))
        othernc[i] <- include.othernc[[i]][2]
      q <- sort(c(othernc,q))
      q <- q[q>1]
    }
    lq <- length(q)
    for(i in 1:lq){
      printobject[[l]][,i+1] <- rep(NA,nmethods)
      for (j in 1:nmethods){
          if(length(x[[j]])>=q[i])
            if(!identical(x[[j]][[q[i]]],NA)){
              if(is.null(unlist(x[[j]][[q[i]]][statistic])))
                printobject[[l]][j,i+1] <- NA
              else
                printobject[[l]][j,i+1] <- x[[j]][[q[i]]][statistic]
            } 
      } # for j
    } # for i
    names(printobject[[l]]) <- c("method",sapply(q,toString))
    print(printobject[[l]],digits=digits)
    l <- l+1
  } # for statistic
  invisible(printobject)
}



