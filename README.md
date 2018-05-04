# General-MSE
A 2-box, generalized management strategy evaluation framework that allows time-varying processes.
See the wiki!

### Installation
To install the `GeMS` package, the [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) package is required.

```
install.packages("devtools")
devtools::install_github("szuwalski/GeMS")
```

### Code Example
```
##----------------
## Run, GeMS, run!
##----------------
# Load package
library(GeMS)

# Set directory where control files are
dir.MSE <- system.file("Cod_1_Production",package="GeMS") #Located in the Examples/ folder of the GeMS directory

# Names of control files to be tested
OMNames<-c("Cod_LowProd_CTL", "Cod_Base_CTL", "Cod_HighProd_CTL")

# Beginning the MSE
run_GeMS(OMNames,dir.MSE)

# Plots can then be found at file.path(dir.MSE, "plots")

##-----------------
## Parallel example
##-----------------
# run_GeMS(OMNames,dir.MSE,
#	        runparallel = T, cores = 3)
```
