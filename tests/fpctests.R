# This tests a few things that are not run in the examples.

library(fpc)

set.seed(4634)
face <- rFace(300,dMoNo=2,dNoEy=0,p=3)
grface <- as.integer(attr(face,"grouping"))
discrproj(face,grface, clnum=1, method="bc")$units
discrproj(face,grface, clnum=1, method="anc")$units
discrproj(face,grface, clnum=1, method="awc")$units


pamk(face,krange=1:5,criterion="ch",usepam=FALSE,critout=TRUE)

set.seed(20000)
face50 <- rFace(50,dMoNo=2,dNoEy=0,p=2)
pamk(dist(face50),krange=1:5,criterion="asw",critout=TRUE)

x <- c(1,2,3,6,6,7,8,120)
ff8 <- fixmahal(x)
summary(ff8)
  # ...dataset a bit too small for the defaults...
ff9 <- fixmahal(x, mnc=3, startn=3)
summary(ff9)

set.seed(776655)
v1 <- rnorm(100)
v2 <- rnorm(100)
d1 <- sample(1:5,100,replace=TRUE)
d2 <- sample(1:4,100,replace=TRUE)
ldata <- cbind(v1,v2,d1,d2)
fr <- flexmixedruns(ldata,
    continuous=2,discrete=2,simruns=1,initial.cluster=c(rep(1,5),rep(2,45),
                                        rep(3,50)),
                    control=list(minprior=0.1),
                    n.cluster=3,allout=FALSE)
print(fr$optsummary)

dface <- dist(face50)


hclusttreeCBI(face50,minlevel=2,method="complete",scaling=TRUE)

noisemclustCBI(face50,G=1:5,emModelNames="VVV",nnk=2)

distnoisemclustCBI(dface,G=5,emModelNames="EEE",nnk=2,
                        mdsmethod="classical",
                        mdsdim=2)

mahalCBI(face50,clustercut=0.5)

library(trimcluster)
disttrimkmeansCBI(dface,k=4,scaling=TRUE,trim=0.1,
                         mdsmethod="classical",
                         mdsdim=2)


complete3 <- cutree(hclust(dface),3)

cluster.stats(dface,complete3,G2=TRUE)

library(diptest)

data(crabs)
dc <- crabs[,4:8]
cmo <- mclustBIC(crabs[,4:8],G=9,modelNames="EEE")
set.seed(12345)
cm <- mclustBIC(crabs[,4:8],G=9,modelNames="EEE",
                initialization=list(noise=(1:200)[sample(200,50)]))


scm <- summary(cm,crabs[,4:8])
scmo <- summary(cmo,crabs[,4:8])

summary(mergenormals(crabs[,4:8],scm,method="ridge.ratio",by=0.05))
summary(mergenormals(crabs[,4:8],scmo,method="ridge.uni",by=0.05))
summary(mergenormals(crabs[,4:8],scm,method="diptantrum",by=0.05))
summary(mergenormals(crabs[,4:8],scmo,method="dipuni",by=0.05))
summary(mergenormals(crabs[,4:8],scm,method="predictive",M=2))

