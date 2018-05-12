#' Compare model estimates
#' 
#' @description Used within the \code{\link{AgeStructureComp}} function to compare true values with estimated values.
#' 
#' @param input1 Array of true quantities
#' @param input2 Array of estimated quantities 
#' @param title Quantity being plotted, e.g. "Catch", "Survey Index", etc.
#' @param CTLNameList Vector of CTL file names
#' @param plotNames Vector of the same length as CTLNameList with "prettier" names for plotting
#' @param nSim Number of replicates per OM, specified Nsim in CTL file
#'
#' @return Plot comparing age-structured estimating models
#'
#' @export
#' 
plot_stuff_inside<-function(input1,input2,title,CTLNameList,plotNames=NA,nSim=Nsim)
{ 
  if(is.na(plotNames[1])) plotNames<-CTLNameList
  par(mfcol=c(1,length(CTLNameList)),mar=c(.1,.1,.1,.1),oma=c(4,6,1,4))
  if(nSim>1) 
  {
    for(y in seq_along(CTLNameList))
    {
      boxplot(input1[,,y],type="l",ylim=c(0,max(input1,input2,na.rm=T)),
              las=1,xaxt='n',ylab='',yaxt='n')
      for(x in 1:nrow(  input2))
        lines(input2[x,,y],col='#ff000066')
      if(y==1)
      {
        axis(side=2,las=1)
        mtext(side=2,title,line=5,cex=.9)
      }
      axis(side=1)
      mtext(side=3,bty='n',text=plotNames[y])
    }
  }
  
  if(nSim==1) {
    for(y in seq_along(CTLNameList))
    {
      plot(input1[,,y],type="l",ylim=c(0,max(input1,input2,na.rm=T)),
           las=1,xaxt='n',ylab='',yaxt='n')
      for(x in 1:nrow(  input2))
        lines(input2[x,,y],col='#ff000066')
      if(y==1)
      {
        axis(side=2,las=1)
        mtext(side=2,title,line=5,cex=.9)
      }
      mtext(side=3,bty='n',text=plotNames[y])
      axis(side=1)  
    }
  }
  
  legend('topleft',col=c(1,2),pch=c(15,NA),lty=c(NA,1),legend=c("True","Estimated"),bty='n')
}