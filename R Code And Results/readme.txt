--------------------------------------------------------------------------------
Readme:
--------------------------------------------------------------------------------

This zip file contains the R code for simulation and analysis in the paper "Spatial Variable Selection and An Application to Virginia Lyme Disease Emergence," authored by Yimeng Xie, Li Xu, Jie Li, Xinwei Deng, Yili Hong, Korine Kolivras, and David N. Gaines.

There are eight files contained in this zip file, including this readme file. The file names are:

1. readme.txt
2. Functions.r
3. Data-Analysis.r
4. Simulation.r
5. lyme.svs.eco0.fit.obj.RData
6. lyme.svs.eco1.fit.obj.RData
7. boot.est.eco0.mat.RData
8. boot.est.eco1.mat.RData


The code files consist of all the information needed to replicate the key results in the paper.

"Functions.r" contains all the R functions.

"Simulation.r" contains the R codes for generating the results in Tables 2-8.

Note that parallel computing is used in the simulation study, and high-performance computers are required. For this paper, the computing for the simulation was done using the R parallel package 'snowfall'. The Linux system used is SMP x86_64 GNU/Linux, with 2 x E5-2683v4 2.1GHz CPU (32 cores). The running time for each simulation scenario varies but it could take up to 60 hours. 

"Data-Analysis.r" contains the R code for the data analysis results in Section 6. Besides the code files, additional files are available that might be helpful for readers to reproduce the result in Section 6 of the paper. In particular,

"lyme.svs.eco0.fit.obj.RData" contains the final model fitting results for ecoregion 0.
"lyme.svs.eco1.fit.obj.RData" contains the final model fitting results for ecoregion 1.
"boot.est.eco0.mat.RData" contains the bootstrap estimates for ecoregion 0.
"boot.est.eco1.mat.RData" contains the bootstrap estimates for ecoregion 1.


--------------------------------------------------------------------------------
Note:
--------------------------------------------------------------------------------

There is an R package named "SpatialVS" (version 1.1) that contains the Virginia Lyme disease data and the major functions for spatial variable selections.

--------------------------------------------------------------------------------
End of Readme
--------------------------------------------------------------------------------
