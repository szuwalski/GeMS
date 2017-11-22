#----------------
# Run, GeMS, run!
#----------------

# Set working directory to where GeMS is located
setwd("~/Box\ Sync/GeMS/General-MSE")

# Directory where control files are
tempdir <- file.path(getwd(), "TestHCRs")

# Names of control files to be tested
OMNames<-c("Cod_AgeStructure_CTL_1","Cod_AgeStructure_CTL_2","Cod_AgeStructure_CTL_3","Cod_AgeStructure_CTL_4")

# Loading the wrapper function
source("run_GeMS.R")

# Beginning the MSE
# CurDir is the directory where control files are
# CreateFolderNameList is the vector of control file names
# GeMSops is a named list of options used in the GeMS() function
#	silent = T sends console output into a console.log text file
#	Adoptions = character string of ADMB options, e.g. -nohess, -gbs, -cbs etc.
#	runparallel = T runs control files in parallel
#	cores = integer stating number of cores to use in parallel processing

run_GeMS(CurDir=tempdir, CreateFolderNameList=OMNames,
		 GeMSops = list(silent = T, ADoptions = "-gbs 2000000000"))

#-----------------
# Parallel example
#-----------------
setwd("~/Box\ Sync/GeMS/General-MSE")

tempdir <- file.path(getwd(), "TestHCRs")
OMNames<-c("Cod_AgeStructure_CTL_1","Cod_AgeStructure_CTL_2","Cod_AgeStructure_CTL_3","Cod_AgeStructure_CTL_4")

source("run_GeMS.R")
run_GeMS(CurDir=tempdir, CreateFolderNameList=OMNames,
	     runparallel = T, cores = 2,
	     GeMSops = list(silent = T, ADoptions = "-gbs 2000000000"))
