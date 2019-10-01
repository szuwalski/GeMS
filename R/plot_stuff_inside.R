#' Compare model estimates
#' 
#' @description Used within the \code{\link{AgeStructureComp}} function to compare true values with estimated values.
#' 
#' @param input1 Array of true quantities
#' @param input2 Array of estimated quantities 
#' @param title Quantity being plotted, e.g. "Catch", "Survey Index", etc.
#' @param CTLNameList Vector of CTL file names
#' @param plotNames Vector of the same length as CTLNameList with "prettier" names for plotting
#' @param runs2plot Matrix of runs to plot
#'
#' @return Plot comparing age-structured estimating models
#'
#' @export
#' 
plot_stuff_inside<-function(input1,input2,title,CTLNameList,runs2plot,plotNames=CTLNameList,ProjYr=out$OM$InitYear)
{ 
  if(length(dim(input1))<3) {
    temp<-input1
    input1<-array(temp,dim=c(dim(temp),length(CTLNameList)))
  }
  if(length(dim(input2))<3) {
    temp<-input2
    input2<-array(temp,dim=c(dim(temp),length(CTLNameList)))
  }
  
  par(mar=c(2,.5,.1,.1),oma=c(2,6,1,2))

  plotProj<- (dim(input1)[2]-ProjYr)>2

  if(plotProj) {
  inmat <- matrix(nrow=3,ncol=length(CTLNameList),byrow=T,
                  c((seq_along(CTLNameList)*2-1),
                    (seq_along(CTLNameList)*2-1),
                    (seq_along(CTLNameList)*2)))

  layout(inmat)
  }

  if(!plotProj) par(mfcol=c(1,length(CTLNameList)))

  if(dim(runs2plot)[2]>1) 
  {
    for(y in seq_along(CTLNameList))
    {
      boxplot(input1[,,y],type="l",ylim=c(0,max(input1,input2,na.rm=T)),
              las=1,xaxt='n',ylab='',yaxt='n')
      for(x in runs2plot[,y])
        lines(input2[x,,y],col='#ff000066')
      if(y==1)
      {
        axis(side=2,las=1)
        mtext(side=2,title,line=5,cex=.9)
      }
      axis(side=1, at=c(1,ProjYr,dim(input1[,,y])[2]))
      mtext(side=3,bty='n',text=plotNames[y])
      abline(v=ProjYr,lty=3)

      if(y==length(CTLNameList))   legend('topleft',col=c(1,2,1),pch=c(15,NA,NA),lty=c(NA,1,3),legend=c("Observed","Predicted","Projection begins"),bty='n')

      #### Zoom in to projected years

    if(plotProj) {
     boxplot(input1[,ProjYr:dim(input1)[2],y],type="l",ylim=c(0,max(input1[,ProjYr:dim(input1)[2],],input2[,ProjYr:dim(input2)[2],],na.rm=T)),
              las=1,xaxt='n',ylab='',yaxt='n')
      for(x in runs2plot[,y])
        lines(input2[x,ProjYr:dim(input1)[2],y],col='#ff000066')
      if(y==1)
      {
        axis(side=2,las=1)
        mtext(side=2,paste("Projected",title,sep=" "),line=5,cex=.9)
      }
       axis(side=1,labels=c((ProjYr+1),round((length(ProjYr:dim(input1)[2]))/2)+ProjYr,dim(input1)[2]),at=c(1,round((length(ProjYr:dim(input1)[2])+1)/2),length(ProjYr:dim(input1)[2])))
      }
    }
  }
  
  if(dim(runs2plot)[2]==1) {
    for(y in seq_along(CTLNameList))
    {
      boxplot(input1[,,y],type="l",ylim=c(0,max(input1,input2,na.rm=T)),
           las=1,xaxt='n',ylab='',yaxt='n')
      for(x in runs2plot[,y])
        lines(input2[x,,y],col='#ff000066')
      if(y==1)
      {
        axis(side=2,las=1)
        mtext(side=2,title,line=5,cex=.9)
      }
      mtext(side=3,bty='n',text=plotNames[y])
      axis(side=1)  
      abline(v=ProjYr,lty=3)

    if(plotProj) {
      boxplot(input1[,ProjYr:dim(input1)[2],y],type="l",ylim=c(0,max(input1[,ProjYr:dim(input1)[2],],input2[,ProjYr:dim(input2)[2],],na.rm=T)),
              las=1,xaxt='n',ylab='',yaxt='n')
      for(x in runs2plot[,y])
        lines(input2[x,ProjYr:dim(input1)[2],y],col='#ff000066')
      if(y==1)
      {
        axis(side=2,las=1)
        mtext(side=2,paste("Projected",title,sep=" "),line=5,cex=.9)
      }
       axis(side=1,labels=c((ProjYr+1),round((length(ProjYr:dim(input1)[2]))/2)+ProjYr,dim(input1)[2]),at=c(1,round((length(ProjYr:dim(input1)[2])+1)/2),length(ProjYr:dim(input1)[2])))
     }
     }
    }
    par(mfcol=c(1,1))
    mtext(side=1,line=3,"Time")
  
  }


