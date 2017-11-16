#----------------
# Run, GeMS, run!
#----------------

# Set working directory to where GeMS is located
setwd("~/Box\ Sync/GeMS/General-MSE")

# Directory where control files are
tempdir <- file.path(getwd(), "Cod_5_AgeStructure")

# Names of control files to be tested
OMNames<-c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL")

# Loading the wrapper function
source("run_GeMS.R")

# Beginning the MSE
# CurDir is the directory where control files are
# CreateFolderNameList is the vector of control file names
# GeMSops is a named list of options used in the GeMS() function
#	silent = T sends console output into a console.log text file
#	Adoptions = character string of ADMB options, e.g. -nohess, -gbs, -cbs etc.
#	runparallel = T runs control files in parallel*
#		* Requires packages foreach and doParallel, and number of cores must be registered
#		  E.g.:
#		  doParallel::registerDoParallel(parallel::detectCores()-2);  # where 2 is the number of cores to keep free
#
run_GeMS(CurDir=tempdir, CreateFolderNameList=OMNames,
		 GeMSops = list(silent = T, ADoptions = "-gbs 2000000000"))

#-----------------
# Parallel example
#-----------------
setwd("~/Box\ Sync/GeMS/General-MSE")

nCores <- parallel::detectCores() - 4
doParallel::registerDoParallel(nCores)

tempdir <- file.path(getwd(), "Cod_5_AgeStructure")
OMNames<-c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL")

source("run_GeMS.R")
run_GeMS(CurDir=tempdir, CreateFolderNameList=OMNames,
	     runparallel = T,
	     GeMSops = list(silent = T, ADoptions = "-gbs 2000000000"))
