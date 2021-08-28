# RVMMAT

RVMMAT is an R program that performs retrospective association testing for longitudinal binary traits in population samples. 


## Installation Instructions:

### Required software and packages
    
1. R (http://www.r-project.org/)
    
2. Package    expm, stats, mvtnorm, MASS, utils
    
3. PLINK 1.0 or 1.9 (https://www.cog-genomics.org/plink2)

Please install the required R package before you install LBRAT package. Please install the **RVMMAT** as following steps.

 
### Install RVMMAT on LUNIX or Mac OSX

```
git clone https://github.com/ZWang-Lab/RVMMAT.git

R CMD INSTALL RVMMAT

```

Alternatively, one can install it in R using the following code.
### Install LBRAT in R
```
library(devtools)
devtools::install_github(repo = 'ZWang-Lab/RVMMAT')

```
## Usage instructions

All functions and examples in the RVMMAT are available in the manual (https://github.com/ZWang-Lab/RVMMAT/blob/master/RVMMAT_0.0.1.0.pdf).


