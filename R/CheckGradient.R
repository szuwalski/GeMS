#' Extracts gradients from age-structured model runs for use in \code{\link{AgeStructureComp}}
#'
#' @param CTLName Name of CTL file
#' @param out List from \code{\link{ReadCTLfile}}
#' @param MSEdir Directory containing CTL file
#'
#' @return Matrix of gradients
#'
#' @export
#' @examples
#' \dontrun{
#' MSEdir <- "~/GeMS/Examples/Cod_1_Production"
#' CTLName <- "Cod_Base_CTL"
#' out <- ReadCTLfile(CTLName)
#' CheckGradient(CTLName,out,MSEdir)
#' }

CheckGradient<-function(CTLName,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)

Gradients	<-matrix(ncol=Nsim,nrow=SimYear)

for(n in 1:Nsim)
for(v in (InitYear+1):SimYear)
{
IndSimFolder<-file.path(MSEdir,CTLName,n,v)
PAR<-readLines(file.path(IndSimFolder,"SimAss.par"))

temp<-grep("gradient",PAR)[1]
GradientLine<-(unlist(strsplit(PAR[temp],split=" ")))
Gradients[v,n]<-as.numeric(GradientLine[length(GradientLine)])
}
return(Gradients)
}