# Spatial-Variable-Selection-and-An-Application-to-Virginia-Lyme-Disease-Emergence

# Author Contributions Checklist Form

## Data

### Abstract 

The Lyme disease dataset used for this paper contains the number of disease cases occurred in Virginia from 2006 to 2011, and demographic data and land cover data in Virginia. The demographic data (e.g., population density, median income, average age) are from the 2010 US census. Land cover data were obtained from the Multi-Resolution Land Cover Consortium for 2006.

### Availability 

The data are available via the R package “SpatialVS” (version 1.1), which can be downloaded from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).

### Description

In the R package, the data object named “lyme.svs.eco0.dat” consists of the disease counts for each census tract from 2006 to 2011, covariate information, population counts, and the location of the centroid for each census tract for ecoregion 0 (i.e., Northern Piedmont, Blue Ridge, Ridge and Valley and Central Appalachian areas). The data object named “lyme.svs.eco1.dat” consists of similar information for ecoregion 1 (i.e., in Piedmont, Middle Atlantic Coastal Plain, and Southeastern Plains areas). A detailed description of the data is available in the documentation of the R package.

## Code

### Abstract

The directory file named “R Code And Results” consists of all the information needed to replicate the key results in the paper.

### Description 

The code files are in .R format, which can be edited by common text editors. In particular, the simulation code is available in “Simulation.r”, and the data analysis code is available in “Data- Analysis.r”.

## Instructions for Use

### Reproducibility 

The files in the subdirectory can be used to reproduce results in Sections 5 and 6 in the paper. Parallel computing is used in the simulation study and high-performance computers are required. Code for Tables 2-8 can be found in “Simulation.r.” The data analysis results in Section 6 can be reproduced by using the code in “Data-Analysis.r.” All the functions are contained in “Functions.r.” Besides the code files, additional files are available that might be helpful for readers to reproduce the result in Section 6 of the paper. In particular,

* “lyme.svs.eco0.fit.obj.RData” contains the final model fitting results for ecoregion 0. 
* “lyme.svs.eco1.fit.obj.RData” contains the final model fitting results for ecoregion 1. 
* “boot.est.eco0.mat.RData” contains the bootstrap estimates for ecoregion 0.
* “boot.est.eco1.mat.RData” contains the bootstrap estimates for ecoregion 1.

The readme file in the folder contains the details. The R version used for computing is R 3.3.1. The required packages and their versions are MASS 7.3-45, nlme 3.1-128, glmmLasso 1.5.1, xtable 1.8-2, and SpatialVS 1.1.

### Notes

There is an R package named "SpatialVS" (version 1.1) that contains the Virginia Lyme disease data and the major functions for spatial variable selection.
