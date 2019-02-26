################################################################################
#December 09, 2018, R script for data analysis

################################################################################
#load the Virginia Lyme disease data from the SpatialVS 1.1 package, and 
#compute the distance matrix for each eco-region
library(SpatialVS)
#load ecoregion 0 data
data("lyme.svs.eco0")
lyme.svs.eco0.dat$dist=distmat.compute(location=lyme.svs.eco0.dat$location, dist.min=0.4712249)
#load ecoregion 1 data
data("lyme.svs.eco1")
lyme.svs.eco1.dat$dist=distmat.compute(location=lyme.svs.eco1.dat$location, dist.min=0.2821849)
detach(package:SpatialVS)


################################################################################
source("Functions.r")
#Note: not all R functions in this paper are implemented in the R package SpatialVS 1.1.
#Thus it is recommended to use the functions in "Functions.r" to reproduce the 
#results in this paper, in case that the SpatialVS package may update in the future.
################################################################################
#plot of pairwise correlations
cov.corr.plot()
################################################################################
#Variable selection 
################################################################################
#Conduct variable selection, run on Linux machine
#Note: for the final results, we ran a fine grid of 20x40 for alpha and lambda, which
#could take multiple days on a Linux server.
##Ecoregion 0
svs.res.eco0.obj<-SpatialVS(dat.obj=lyme.svs.eco0.dat, alpha.vec=seq(0.05, 1,,20), lambda.vec=seq(0.05, 7, len=40), method="PQL")
lyme.svs.eco0.fit.obj=svs.res.eco0.obj$PQL.best.obj
##Ecoregion 1
svs.res.eco1.obj<-SpatialVS(dat.obj=lyme.svs.eco1.dat, alpha.vec=seq(.05, 1,, 20), lambda.vec=seq(3, 15,, 40), method="PQL")
lyme.svs.eco1.fit.obj=svs.res.eco1.obj$PQL.best.obj
#End of run on Linux machine
################################################################################

################################################################################
#create data objects for bootstrap
#one can use the saved intermediate objects for examples
load("lyme.svs.eco0.fit.obj.RData") 
load("lyme.svs.eco1.fit.obj.RData") 

#Ecoregion 0
lyme.eco0.boot.dat=list(dat.obj=lyme.svs.eco0.dat, fit.obj=lyme.svs.eco0.fit.obj)
#Ecoregion 1
lyme.eco1.boot.dat=list(dat.obj=lyme.svs.eco1.dat, fit.obj=lyme.svs.eco1.fit.obj)

save(lyme.eco0.boot.dat, file="lyme.eco0.boot.dat")
save(lyme.eco1.boot.dat, file="lyme.eco1.boot.dat")

################################################################################
#check overall fitting
overall.exp.vs.obs(obj0=lyme.eco0.boot.dat, obj1=lyme.eco1.boot.dat)
################################################################################
 
################################################################################
#parameteric bootstrap, run on Linux machine, using snowfall
library(snowfall)
no.nodes=1
no.cores=32*no.nodes-1
per.core=35
sfInit(parallel=TRUE, cpus=no.cores,type="SOCK")
sfSource("Functions.r")
load(file="lyme.eco0.boot.dat")
load(file="lyme.eco1.boot.dat")
control.default=list(maxIter=20,iwls=10^(-4),tol1=10^(-3),tol2=10^(-3)) 
sfExportAll()
sfClusterSetupRNG(seed=round(2^32*runif(1)))

##Ecoregion 0
tmp=sfClusterApplyLB(1:(per.core*no.cores), lyme.fit.boot, boot.dat=lyme.eco0.boot.dat, alpha.vec=seq(0.7, .9,, 5), lambda.vec=seq(4.5, 6.5,,5)) 
save(tmp,file="tmp.lyme.eco0.fit.boot.1")

##Ecoregion 1
tmp=sfClusterApplyLB(1:(per.core*no.cores), lyme.fit.boot, boot.dat=lyme.eco1.boot.dat, alpha.vec=seq(0.7, .8,, 5), lambda.vec=seq(3.5, 5,,5)) #
save(tmp,file="tmp.lyme.eco1.fit.boot.2")
sfStop()

##collect results
boot.est.eco0.mat=lyme.boot.res.collection(filename="tmp.lyme.eco0.fit.boot", seqs=c(1)) 
boot.est.eco1.mat=lyme.boot.res.collection(filename="tmp.lyme.eco1.fit.boot", seqs=c(2)) 

save(boot.est.eco0.mat, file="boot.est.eco0.mat")
save(boot.est.eco1.mat, file="boot.est.eco1.mat")
################################################################################
##End of run on Linux machine
################################################################################
#one can use the saved intermediate objects for examples
load("boot.est.eco0.mat.RData")
load("boot.est.eco1.mat.RData")

eco0.ci.mat=lyme.boot.est.ci(boot.dat=lyme.eco0.boot.dat, est.mat=boot.est.eco0.mat) 
eco1.ci.mat=lyme.boot.est.ci(boot.dat=lyme.eco1.boot.dat, est.mat=boot.est.eco1.mat) 

##Table 9
#general table for parameter estimates and confidence intervals
xx=round(cbind(eco0.ci.mat[,c(1, 4,5)], eco1.ci.mat[,c(1, 4,5)]), 3)
print(xtable(xx, digits=3),sanitize.text.function = function(x) {x})

################################################################################
##End########################################################################### 
