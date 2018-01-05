run_GeMS <- function(CreateFolderNameList,GeMSDir,CurDir,
					 runparallel=F,cores = 1,GeMSops=NA,shutup=F) {

	if(!file.exists(file.path(GeMSDir,"General_MSE_main.R"))) {
		stop(print("Set GeMSdir to folder where GeMS is located. Full file path required."))
	}

	if(!file.exists(file.path(CurDir,paste0(CreateFolderNameList[1],".csv")))) {
		stop(print("Set CurDir to folder where CTL files are located. Full file path required."))
	}
	
	if(sum(grep("ReadCTLfile",ls()))==0) {
		source(file.path(GeMSDir,"General_MSE_helpers.R"))
		source(file.path(GeMSDir,"General_MSE_main.R"))
	}

	if(runparallel) {
		if(!requireNamespace("foreach")) {
			stop(print(paste0("Packages foreach and doParallel need to be installed, and foreach needs to be loaded. Example:\n",
							  "install.packages(c('foreach','doParallel')\n",
							  "library(foreach)")))
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
				setwd(CurDir)
				Inout<-ReadCTLfile(paste0(CreateFolderNameList[mod],".csv"))
				temp <- list(MSEdir = CurDir, GeMSdir = GeMSDir, out=Inout,
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
				setwd(CurDir)
				Inout<-ReadCTLfile(paste0(CreateFolderNameList[mod],".csv"))
				temp <- list(MSEdir = CurDir, GeMSdir = GeMSDir, out=Inout,
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
			ProductionModelOutput(Inout=OMfile,CTLNames=CreateFolderNameList,MSEdir=CurDir)
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

	setwd(CurDir)
	if(!shutup) message("End of GeMS run. Returned to MSE directory.")

}
