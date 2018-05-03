#' Reference Points
#' Extracts reference points from age-structured model runs for use in \code{\link{AgeStructureComp}}
#'
#' @param CTLName Name of CTL file
#' @param out List from \code{\link{ReadCTLfile}}
#' @param MSEDir Directory containing CTL file
#'
#' @return List of reference points list(B35,F35,OFL,tB35,tF35,tOFL)
#'
#' @export
#' @examples
#' \dontrun{
#' MSEdir <- "~/GeneralMSE/Examples/Cod_5_AgeStructure"
#' OMName <- "Cod_AgeStructure_CTL"
#' out <- ReadCTLfile(OMName)
#' CheckReferencePoints(OMName, out, MSEdir)
#' }

CheckReferencePoints<-function(CTLName,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)

B35	<-matrix(ncol=Nsim,nrow=SimYear)
F35	<-matrix(ncol=Nsim,nrow=SimYear)
OFL	<-matrix(ncol=Nsim,nrow=SimYear) 
tB35	<-matrix(ncol=Nsim,nrow=SimYear)
tOFL	<-matrix(ncol=Nsim,nrow=SimYear)

for(n in 1:Nsim)
for(v in (InitYear+1):SimYear)
{
IndSimFolder<-file.path(MSEdir,CTLName,n,v)
REP<-readLines(file.path(IndSimFolder,"simass.REP"))
TRU<-readLines(file.path(IndSimFolder,"TrueQuantities.DAT"))

temp<-grep("B35",REP)[1]
B35[v,n]<-as.numeric((unlist(strsplit(REP[temp+1],split=" "))))

temp<-grep("F35",REP)[2]
F35[v,n]<-as.numeric((unlist(strsplit(REP[temp+1],split=" "))))

temp<-grep("OFL",REP)[2]
OFL[v,n]<-as.numeric((unlist(strsplit(REP[temp+1],split=" "))))

if(v==SimYear)
{
temp<-grep("B35",TRU)
tB35[,n] <-as.numeric((unlist(strsplit(TRU[temp+n],split=" "))))

temp<-grep("OFL",TRU)
tOFL[,n]<-as.numeric((unlist(strsplit(TRU[temp+n],split=" ")))) 
}
temp<-grep("F35",TRU)
tF35 <-as.numeric((unlist(strsplit(TRU[temp+1],split=" ")))) 


}
list(B35,F35,OFL,tB35,tF35,tOFL)
}
