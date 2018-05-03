#' Function to run GeMS
#' 
#' @description Main function to run a full Generalized MSE.
#' 
#' @param CTLNameList Vector of CTL files
#' @param MSEdir Directory containing CTL files
#' @param runparallel Logical; operating models should be run in parallel. If you're doing this, make sure that \pkg{foreach} is loaded, using \code{library(foreach)}
#' @param cores Number of cores to be used for parallel runs
#' @param silent Logical; Show output on console
#' @param ... Anything to be passed to \code{\link{GeMS}}, including \code{ADoptions} and \code{ADsilent} for ADMB options.
#'
#' @return \code{\link{GeMS}} output
#'
#' @seealso \code{\link{GeMS}}
#' @export
#' @examples
#' \dontrun{
#' OMNames <- c("Cod_LowProd_CTL","Cod_Base_CTL","Cod_HighProd_CTL")
#' MSEdir <- "~/GeneralMSE/Examples/Cod_1_Production"
#' run_GeMS(OMNames, MSEdir)
#' }
run_GeMS <- function(CTLNameList,MSEdir=getwd(),
					 runparallel=F,cores = 1,silent=T,...) {

  if(!file.exists(file.path(MSEdir,paste0(CTLNameList[1],".csv")))) stop("Set MSEdir to directory containing CTL files.")
	if(runparallel) {
		if(!requireNamespace("foreach",quietly=T)) {
			stop(print(paste0("Packages foreach and doParallel need to be installed, and foreach needs to be loaded. Example:\n",
							  "install.packages(c('foreach','doParallel')\n",
							  "library(foreach)")))
		}
		if(cores == 1) {message(paste0("Only one core registered. Runs will not be done in parallel.\n",
									  "Set cores=2 or greater, depending on number of cores to be registered."))
						runparallel <- F
		}
		if(cores>1) {
			if(!silent) message("Running scenarios in parallel.")
			cl <- parallel::makeCluster(cores)
			doParallel::registerDoParallel(cl)

			foreach(mod=1:length(CTLNameList),.export=ls(envir=globalenv())) %dopar% {
			  Inout<-ReadCTLfile(file.path(MSEdir,CTLNameList[mod]))
			  do.call(GeMS, args = list(out=Inout,CTLName=CTLNameList[mod],
			                            MSEdir=MSEdir,silent=silent,...))
			}
			doParallel::stopImplicitCluster()
		}
	}

	if(!runparallel) {
		if(!silent) message("Running scenarios in sequence.")
		for(mod in 1:length(CTLNameList)) {
				Inout<-ReadCTLfile(file.path(MSEdir,CTLNameList[mod]))
				do.call(GeMS, args = list(out=Inout,CTLName=CTLNameList[mod],
				                          MSEdir=MSEdir,silent=silent,...))
		 }
	 }

	OMfile <- ReadCTLfile(file.path(MSEdir,paste0(rev(CTLNameList)[1])))
	if(OMfile$OM$AssessmentType==1) {
		ProductionModelOutput(out=OMfile,CTLNameList=CTLNameList,MSEdir=MSEdir)
	}
	
	if(!silent) message("End of simulations.")
#	if(length(CTLNameList)>1) {
		if(!silent) message("Creating combined figures.")
		if(OMfile$OM$AssessmentType==2) {
			if(OMfile$OM$SimYear-OMfile$OM$InitYear < 6) {
				retpeels <- OMfile$OM$SimYear-OMfile$OM$InitYear
				message(paste0("Number of years for retrospective peels == ", retpeels, "."))
				AgeStructureComp(out=OMfile,CTLNameList=CTLNameList,RetroPeels=retpeels,MSEdir=MSEdir)
			}
			if(OMfile$OM$SimYear-OMfile$OM$InitYear >= 6) {
				message("Number of years for retrospective peels == 6.")
				AgeStructureComp(out=OMfile,CTLNameList=CTLNameList,MSEdir=MSEdir)
			}
		}
#	}
	if(!silent) message("End of GeMS run.")

}
