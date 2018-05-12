#' Check Retrospective
#' @description Extracts calculated retrospective bias from age-structured model runs for use in \code{\link{AgeStructureComp}}
#'
#' @param RetroPeels Number of years to conduct retrospective analysis
#' @param CTLName Name of CTL file
#' @param PlotRetro 1/0 value to plot retrospective analyses
#' @param out List from \code{\link{ReadCTLfile}}
#' @param MSEDir Directory containing CTL file
#'
#' @return Predicted spawning biomass and true spawning biomass from retrospective analysis.
#'
#' @export
#' @examples
#' \dontrun{
#' CTLName <- "Cod_Base_CTL"
#' out <- ReadCTLfile(CTLName)
#' MSEdir <- "~/GeMS/Examples/Cod_1_Production"
#' CheckRetro(CTLName,out,MSEdir)
#' }

CheckRetro<-function(CTLName,out,MSEdir,RetroPeels=6,PlotRetro=0)
{
  Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
  SimYear		<-out$OM$SimYear			# total number of years in simulation
  InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
  
  predSpBioRetro	<-array(dim=c(RetroPeels,SimYear,Nsim))
  trueSpBioRetro    <-matrix(ncol=SimYear,nrow=Nsim)
  
  for(n in 1:Nsim)
    for(v in 0:(RetroPeels-1))
    {
    IndSimFolder<-file.path(CTLName,n,SimYear-v)
    REP<-readLines(file.path(MSEdir,IndSimFolder,"simass.rep"))
    DAT<-readLines(file.path(MSEdir,IndSimFolder,"simass.DAT"))
    
    temp<-grep("survey years",DAT)[1]
    yearsDat<-suppressWarnings(as.numeric(unlist(strsplit(DAT[temp+1],split=" "))))
    temp<-grep("spawning biomass",REP)
    predSpBioRetro[(v+1),(out$OM$start_assessment:(SimYear-v-1)),n]<-suppressWarnings(as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)])
  
  }
  
  for(n in 1:Nsim)
  {
    IndSimFolder<-file.path(CTLName,n,SimYear)
    TRU<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))
    temp<-grep("spawning biomass",TRU)
    trueSpBioRetro[n,]<-suppressWarnings(as.numeric(unlist(strsplit(TRU[temp+1],split=" "))))
  }
  
  plotSegment<-predSpBioRetro[,((SimYear-RetroPeels-4):SimYear),1]
  if(PlotRetro==1)
  {
    dev.new()
    par(mar=c(4,5,1,1))
    plot(y=plotSegment[1,],x=((SimYear-RetroPeels-4):SimYear),type="l",ylim=c(min(plotSegment,na.rm=T),max(plotSegment,na.rm=T)),lty=2,las=1,
     ylab="",xlab="Year")
    for(x in 2:RetroPeels)
     lines(y=plotSegment[x,],x=((SimYear-RetroPeels-4):SimYear),lty=x+1)
    
     lines(y=trueSpBioRetro[1,((SimYear-RetroPeels-4):SimYear)],x=((SimYear-RetroPeels-4):SimYear),col=2)
    
    legend("topleft",bty='n',lty=1,col=2,"Truth")
  }
  list(predSpBioRetro,trueSpBioRetro)
}
