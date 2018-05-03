#' Pull True Projections
#' Not used anywhere... yet.
#' 
#' @param CTLNameList Vector of CTL file names
#' @param out List output from \code{\link{ReadCTLfile}}
#' @param MSEdir Directory containing CTL files
#'
#' @return Plots comparing age-structured estimating models
#'
#' @export
#' @examples
#' \dontrun{
#' MSEdir <- "~/GeneralMSE/Examples/Cod_5_AgeStructure"
#' OMNames <- c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")
#' out <- ReadCTLfile(OMNames[1])
#' AgeStructureComp(out=out,
#'                  CTLNameList=OMNames,
#'                  MSEdir=MSEdir,
#'                  plotNames=c("Base","Fixed M","Estimate M"))
#' }
#' 
PullTrueProj<-function(CTLNameList,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
MaxAge		<-out$OM$MaxAge
Ages			<-seq(1,MaxAge)
t0			<-out$OM$t0s		# check this later when going to spatial

TrueRec	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))
TrueFmort	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))
TrueCatch	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))
TrueSpbio	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))
TrueOFL	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))

for(p in seq_along(CTLNameList))
{
DrawDir		<-CTLNameList[p]
Inout			<-ReadCTLfile(file.path(MSEdir,paste0(CTLNameList[p],"_CTL")))

for(n in 1:Nsim)
{
IndSimFolder	<-file.path(DrawDir,n,SimYear)
TRU			<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))

#==true
temp		<-grep("fishing mortality",TRU)
TrueFmort[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]

#==true Catch
temp		<-grep("Catch",TRU)
TrueCatch[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]

#==true Spawning biomass
temp		<-grep("spawning biomass",TRU)
TrueSpbio[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]


#==recruitment
temp		<-grep("recruitment",TRU)
TrueRec[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]

#==OFL 
temp			<-grep("OFL",TRU)
TrueOFL[n,,p]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]
}
}
list(TrueRec,TrueFmort,TrueOFL,TrueCatch,TrueSpbio)
}