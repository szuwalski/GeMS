run_GeMS <- function(CurDir,CreateFolderNameList,runparallel=F,cores = 1,GeMSops=NA,shutup=F) {

	if(!file.exists("General_MSE_main.R")) {
		stop(print("Set working directory to folder where GeMS is located."))
	}

	dir.main <- getwd()
	
	source("General_MSE_helpers.R")
	source("General_MSE_main.R")

	if(runparallel) {
		if(!requireNamespace("foreach")) {
			install.packages(c("foreach","doParallel"))
			message("Installed packages foreach and doParallel.")
		}
		suppressPackageStartupMessages(library(foreach))
		if(cores == 1) {message(paste0("Only one core registered. Runs will not be done in parallel.\n",
									  "Set cores=2 or greater, depending on number of cores to be registered."))
						runparallel <- F
					}
		if(cores>1) {
			if(!shutup) message("Running scenarios in parallel.")
			cl <- parallel::makeCluster(cores)
			doParallel::registerDoParallel(cl)

			foreach(mod=1:length(CreateFolderNameList),.export=ls(envir=globalenv())) %dopar% {
				setwd(dir.main)
				Inout<-ReadCTLfile(file.path(CurDir,paste0(CreateFolderNameList[mod],".csv")))
				temp <- list(MSEdir = CurDir, out=Inout,
							 CreateFolderName = CreateFolderNameList[mod])
				if(!is.na(GeMSops)) {GeMSvars <- append(temp,GeMSops)}
				if(is.na(GeMSops)) {GeMSvars <- temp}
		
				do.call(GeMS, args = GeMSvars)
			}
		}
	}

	if(!runparallel) {
		if(!shutup) message("Running scenarios in sequence.")
		for(mod in 1:length(CreateFolderNameList)) {
				setwd(dir.main)
				Inout<-ReadCTLfile(file.path(CurDir,paste0(CreateFolderNameList[mod],".csv")))
				temp <- list(MSEdir = CurDir, out=Inout,
							 CreateFolderName = CreateFolderNameList[mod])
				if(!is.na(GeMSops)) {GeMSvars <- append(temp,GeMSops)}
				if(is.na(GeMSops)) {GeMSvars <- temp}
		
				do.call(GeMS, args = GeMSvars)
		 }
	 }

	if(!shutup) message("End of simulations.")
	if(length(CreateFolderNameList)>1) {
		if(!shutup) message("Creating combined figures.")
		OMfile <- ReadCTLfile(file.path(CurDir,paste0(rev(CreateFolderNameList)[1],".csv")))
		if(OMfile$OM$AssessmentType==1) {
			ProductionModelOutput(Inout=OMfile,CTLNames=CreateFolderNameList,MSEdir = CurDir)
		}

		if(OMfile$OM$AssessmentType==2) {
			if(OMfile$OM$SimYear-OMfile$OM$InitYear < 6) {
				retpeels <- OMfile$OM$SimYear-OMfile$OM$InitYear
				message(paste0("Number of years for retrospective peels == ", retpeels, "."))
				AgeStructureComp(Inout=OMfile,CTLNames=CreateFolderNameList,RetroPeels=retpeels,MSEdir=CurDir)
			}
			if(OMfile$OM$SimYear-OMfile$OM$InitYear >= 6) {
				message("Number of years for retrospective peels == 6.")
				AgeStructureComp(Inout=OMfile,CTLNames=CreateFolderNameList,MSEdir=CurDir)
			}
		}
	}

	setwd(dir.main)
	if(!shutup) message("End of GeMS run. Returned to original GeMS directory.")

}
