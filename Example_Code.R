##----------------
## Run, GeMS, run!
##----------------
# Load package
library(GeMS)

# Set directory where control files are
dir.MSE <- system.file("extdata/Cod_1_Production",package="GeMS") #Located in the extdata/ folder

# Names of control files to be tested
OMNames<-c("Cod_LowProd_CTL", "Cod_Base_CTL", "Cod_HighProd_CTL")

# Beginning the MSE
run_GeMS(OMNames,dir.MSE)

# Plots can then be found at file.path(dir.MSE, "plots")

##-----------------
## Parallel example
##-----------------
#run_GeMS(OMNames,dir.MSE,
#	        runparallel = T, cores = 3)
