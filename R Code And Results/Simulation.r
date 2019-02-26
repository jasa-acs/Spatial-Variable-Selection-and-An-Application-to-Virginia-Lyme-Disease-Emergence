################################################################################
#May 26, 2018, R script for simulation study

#######################################################################
#run on Linux machine using snowfall, see readme for me details
#######################################################################

library(snowfall)
no.nodes=1
no.cores=32*no.nodes-1
per.core=10

sfInit(parallel=TRUE, cpus=no.cores,type="SOCK")
sfSource("Functions.r")

sfExportAll()
sfClusterSetupRNG(seed=round(2^32*runif(1)))

#######################################################################
#Table 2, PQL
#######################################################################
result<-sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=0, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05))
save(result, file="result_case1_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05))
save(result, file="result_case2_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.8), adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05))
save(result, file="result_case3_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.5), adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05))
save(result, file="result_case4_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, noeff_corr=T, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05))
save(result, file="result_case5_adaptive")


#######################################################################
#Table 3, PQL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=0, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.5, 5))
save(result, file="result_case1_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.5, 5))
save(result, file="result_case2_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.8), adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.5, 5))
save(result, file="result_case3_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.5), adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.5, 5))
save(result, file="result_case4_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, noeff_corr=T, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.5, 5))
save(result, file="result_case5_adaptive.cov1")

#######################################################################
#Table 4, PQL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=0, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.1, 10))
save(result, file="result_case1_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.1, 10))
save(result, file="result_case2_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.8), adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.1, 10))
save(result, file="result_case3_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.5), adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.1, 10))
save(result, file="result_case4_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, noeff_corr=T, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), cov.pars=c(0.1, 10))
save(result, file="result_case5_adaptive.cov3")

#######################################################################
#Table 5, PQL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=0, adaptive=T, nocor=T, betaeff=c(0.2, 0.3, 0.4, 0.5, 0.7, 0.8, -0.1, -0.6, -0.9, -1), lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=20)
save(result, file="result_case1_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, adaptive=T, nocor=T, betaeff=c(0.2, 0.3, 0.4, 0.5, 0.7, 0.8, -0.1, -0.6, -0.9, -1), lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=20)
save(result, file="result_case2_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.8), adaptive=T, nocor=T, betaeff=c(0.2, 0.3, 0.4, 0.5, 0.7, 0.8, -0.1, -0.6, -0.9, -1), lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=20)
save(result, file="result_case3_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.5), adaptive=T, nocor=T, betaeff=c(0.2, 0.3, 0.4, 0.5, 0.7, 0.8, -0.1, -0.6, -0.9, -1), lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=20)
save(result, file="result_case4_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, noeff_corr=T, adaptive=T, nocor=T, betaeff=c(0.2, 0.3, 0.4, 0.5, 0.7, 0.8, -0.1, -0.6, -0.9, -1), lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=20)
save(result, file="result_case5_adaptive.p1")

#######################################################################
#Table 6, PQL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=0, adaptive=T, nocor=T,  lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=25)
save(result, file="result_case1_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, adaptive=T, nocor=T,  lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=25)
save(result, file="result_case2_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.8), adaptive=T, nocor=T,  lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=25)
save(result, file="result_case3_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=c(3,2), rho=c(0.8, 0.5), adaptive=T, nocor=T,  lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=25)
save(result, file="result_case4_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=5, rho=0.8, noeff_corr=T, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), p.total=25)
save(result, file="result_case5_adaptive.p2")


#######################################################################
#Table 2, APL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case1_adaptive", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case1_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case2_adaptive", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case2_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case3_adaptive", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case3_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case4_adaptive", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case4_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case5_adaptive", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case5_adaptive")

#######################################################################
#Table 3, APL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case1_adaptive.cov1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case1_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case2_adaptive.cov1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case2_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case3_adaptive.cov1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case3_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case4_adaptive.cov1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case4_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case5_adaptive.cov1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case5_adaptive.cov1")

#######################################################################
#Table 4, APL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case1_adaptive.cov3", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case1_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case2_adaptive.cov3", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case2_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case3_adaptive.cov3", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case3_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case4_adaptive.cov3", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case4_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case5_adaptive.cov3", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case5_adaptive.cov3")

#######################################################################
#Table 5, APL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case1_adaptive.p1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case1_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case2_adaptive.p1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case2_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case3_adaptive.p1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case3_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case4_adaptive.p1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case4_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case5_adaptive.p1", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case5_adaptive.p1")

#######################################################################
#Table 6, APL
#######################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case1_adaptive.p2", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case1_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case2_adaptive.p2", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case2_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case3_adaptive.p2", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case3_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case4_adaptive.p2", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case4_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_case5_adaptive.p2", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_case5_adaptive.p2")


################################################################################
#Table S2, comparisons with existing method.
################################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case1_adaptive", lasso.formula=lasso.formula)
save(result, file="resultExM_case1_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case2_adaptive", lasso.formula=lasso.formula)
save(result, file="resultExM_case2_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case3_adaptive", lasso.formula=lasso.formula)
save(result, file="resultExM_case3_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case4_adaptive", lasso.formula=lasso.formula)
save(result, file="resultExM_case4_adaptive")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case5_adaptive", lasso.formula=lasso.formula)
save(result, file="resultExM_case5_adaptive")

################################################################################
#Table S3, comparisons with existing method.
################################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case1_adaptive.cov1", lasso.formula=lasso.formula)
save(result, file="resultExM_case1_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case2_adaptive.cov1", lasso.formula=lasso.formula)
save(result, file="resultExM_case2_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case3_adaptive.cov1", lasso.formula=lasso.formula)
save(result, file="resultExM_case3_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case4_adaptive.cov1", lasso.formula=lasso.formula)
save(result, file="resultExM_case4_adaptive.cov1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case5_adaptive.cov1", lasso.formula=lasso.formula)
save(result, file="resultExM_case5_adaptive.cov1")


################################################################################
#Table S4, comparisons with existing method.
################################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case1_adaptive.cov3", lasso.formula=lasso.formula)
save(result, file="resultExM_case1_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case2_adaptive.cov3", lasso.formula=lasso.formula)
save(result, file="resultExM_case2_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case3_adaptive.cov3", lasso.formula=lasso.formula)
save(result, file="resultExM_case3_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case4_adaptive.cov3", lasso.formula=lasso.formula)
save(result, file="resultExM_case4_adaptive.cov3")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case5_adaptive.cov3", lasso.formula=lasso.formula)
save(result, file="resultExM_case5_adaptive.cov3")

################################################################################
#Table S5, comparisons with existing method.
################################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case1_adaptive.p1", lasso.formula=lasso.formula)
save(result, file="resultExM_case1_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case2_adaptive.p1", lasso.formula=lasso.formula)
save(result, file="resultExM_case2_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case3_adaptive.p1", lasso.formula=lasso.formula)
save(result, file="resultExM_case3_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case4_adaptive.p1", lasso.formula=lasso.formula)
save(result, file="resultExM_case4_adaptive.p1")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case5_adaptive.p1", lasso.formula=lasso.formula)
save(result, file="resultExM_case5_adaptive.p1")


################################################################################
#Table 7, comparisons with existing method.
################################################################################
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case1_adaptive.p2", lasso.formula=lasso.formula)
save(result, file="resultExM_case1_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case2_adaptive.p2", lasso.formula=lasso.formula)
save(result, file="resultExM_case2_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case3_adaptive.p2", lasso.formula=lasso.formula)
save(result, file="resultExM_case3_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case4_adaptive.p2", lasso.formula=lasso.formula)
save(result, file="resultExM_case4_adaptive.p2")

result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_case5_adaptive.p2", lasso.formula=lasso.formula)
save(result, file="resultExM_case5_adaptive.p2")


################################################################################
#Table 8, comparision with Covariate simulated from real data.
################################################################################
#PQL
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu, numofcor=0, adaptive=T, nocor=T, lambda.vec.ext=seq(0.01, 1, by=0.05), realX=T, realX.file="lyme.eco0.boot.dat")
save(result, file="result_realX_adaptive")
#APL
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_LP, name="result_realX_adaptive", adaptive=T, nocor=T, lambda.vec=c(0.1, seq(0.5, 5, by=0.5), 6:10))
save(result, file="resultLP_realX_adaptive")
#Existing methods
result=sfClusterApplyLB(1:(per.core*no.cores), spatialVS_simu_existing_methods, name="result_realX_adaptive")
save(result, file="resultExM_realX_adaptive")

sfStop()

###########################################################################
#summarize all the simulation results into tables
#######################################################################
#Table 2
spatialVS_sim_coll.APL.PQL.print.tab(file.name.ext="", betaeff=c(-0.5, 0.75, 1, -0.75, -1), cov.pars=c(0.1, 5), p.zero=10)
#######################################################################
#Table 3
spatialVS_sim_coll.APL.PQL.print.tab(file.name.ext=".cov1", betaeff=c(-0.5, 0.75, 1, -0.75, -1), cov.pars=c(0.5, 5), p.zero=10)
#######################################################################
#Table 4
spatialVS_sim_coll.APL.PQL.print.tab(file.name.ext=".cov3", betaeff=c(-0.5, 0.75, 1, -0.75, -1), cov.pars=c(0.1, 10), p.zero=10)
#######################################################################
#Table 5
spatialVS_sim_coll.APL.PQL.print.tab(file.name.ext=".p1", betaeff=c(0.2, 0.3, 0.4, 0.5, 0.7, 0.8, -0.1, -0.6, -0.9, -1), cov.pars=c(0.5, 5), p.zero=10)
#######################################################################
#Table 6
spatialVS_sim_coll.APL.PQL.print.tab(file.name.ext=".p2", betaeff=c(-0.5, 0.75, 1, -0.75, -1), cov.pars=c(0.1, 5), p.zero=20) 
################################################################################
#Table S2
spatialVS_sim_coll.ExM.print.tab(file.name.ext="")
################################################################################
#Table S3
spatialVS_sim_coll.ExM.print.tab(file.name.ext=".cov1")
################################################################################
#Table S4
spatialVS_sim_coll.ExM.print.tab(file.name.ext=".cov3")
################################################################################
#Table S5
spatialVS_sim_coll.ExM.print.tab(file.name.ext=".p1", betaeff=c(0.2, 0.3, 0.4, 0.5, 0.7, 0.8, -0.1, -0.6, -0.9, -1), p.zero=10)
################################################################################
#Table 7
spatialVS_sim_coll.ExM.print.tab(file.name.ext=".p2", p.zero=20)
################################################################################
#Table 8,
spatialVS_sim_coll.realX.print.tab()
################################################################################
###End of simulation############################################################










