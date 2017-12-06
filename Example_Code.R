##----------------
## Run, GeMS, run!
##----------------

# Set directory where GeMS is located
dir.GeMS <- "~/Box\ Sync/GeMS/General-MSE"

# Set directory where control files are
#dir.MSE <- ""~/Box\ Sync/GeMS/General-MSE"/TestHCRs"
dir.MSE <- file.path(dir.GeMS, "TestHCRs")

# Names of control files to be tested
OMNames<-c("Cod_AgeStructure_CTL_1","Cod_AgeStructure_CTL_2","Cod_AgeStructure_CTL_3","Cod_AgeStructure_CTL_4")

# Loading the wrapper function
source(file.path(dir.GeMS,"run_GeMS.R"))

# Beginning the MSE
# GeMSDir is the directory where GeMS files are
# CurDir is the directory where control files are
# CreateFolderNameList is the vector of control file names
# GeMSops is an optional named list of options used in the GeMS() function
#	silent = T sends ADMB console output into a console.log text file
#	Adoptions = character string of ADMB options, e.g. -nohess, -gbs, -cbs etc.
#	runparallel = T runs control files in parallel
#	cores = integer stating number of cores to use in parallel processing
# shutup is a T/F option to display run_GeMS() messages

run_GeMS(GeMSDir=dir.GeMS, CurDir=dir.MSE,
		 CreateFolderNameList=OMNames)

##----------------------
## Example using GeMSops
##----------------------
#run_GeMS(GeMSDir=dir.GeMS, CurDir=dir.MSE,
#		 CreateFolderNameList=OMNames,
#		 GeMSops = list(silent = T, ADoptions = "-gbs 2000000000"))

##-----------------
## Parallel example
##-----------------
#run_GeMS(GeMSDir=dir.GeMS, CurDir=dir.MSE,
#		 CreateFolderNameList=OMNames,
#		 GeMSops = list(silent = T, ADoptions = "-gbs 2000000000")
#	     runparallel = T, cores = 4)
