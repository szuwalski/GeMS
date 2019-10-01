#' Compare Production Model Outputs
#' 
#' Analyse and compare output from production estimating models
#' 
#' @param out List ouput from \code{\link{ReadCTLfile}}
#' @param CTLNameList List of CTL files
#' @param MSEdir Directory containing CTL files
#' @param ylimIN Vector of optional manual limits for y-axis
#' @param plotNames Vector of the same length as CTLNameList with "prettier" names for plotting
#'
#' @return GeMS Object # To be implemented
#'
#' @export
#' @examples
#' \dontrun{
#' MSEdir <- "~/GeneralMSE/Examples/Cod_1_Production"
#' OMNames <- c("Cod_LowProd_CTL","Cod_Base_CTL","Cod_HighProd_CTL")
#' out <- ReadCTLfile(OMNames[1])
#' ProductionModelOutput(out,OMNames,MSEdir,plotNames=c("Low Productivity", "Base", "High Productivity"))
#' }
ProductionModelOutput<-function(out,CTLNameList,MSEdir,ylimIN=NA,plotNames=NA)
{
if(is.na(plotNames[1])) plotNames<-CTLNameList
Data<-rep(list(list()))
for(x in seq_along(CTLNameList))
 Data[[x]]<-readLines(file.path(MSEdir,CTLNameList[x],"ProdOutputs.csv"))

SimYear<-out$OM$SimYear

#==estimated and true CPUE
estCPUE<-array(dim=c(out$OM$Nsim,out$OM$SimYear,length(CTLNameList)))
trueCPUE<-array(dim=c(out$OM$Nsim,out$OM$SimYear,length(CTLNameList)))
trueBMSY<-rep(0,length(CTLNameList))
estBMSY<-array(dim=c(out$OM$Nsim,out$OM$SimYear,length(CTLNameList)))
trueFMSY<-rep(0,length(CTLNameList))
estFMSY<-array(dim=c(out$OM$Nsim,out$OM$SimYear,length(CTLNameList)))
estTAC<-array(dim=c(out$OM$Nsim,out$OM$SimYear,length(CTLNameList)))
trueTAC<-array(dim=c(out$OM$Nsim,out$OM$SimYear,length(CTLNameList)))

for(y in seq_along(CTLNameList))
{
 tmp<-grep("est cpue",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+out$OM$Nsim)]
for(x in seq_along(Quant))
 estCPUE[x,,y]<-suppressWarnings(as.numeric(unlist(strsplit(Quant[x],split=" "))))

tmp<-grep("true CPUE",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+out$OM$Nsim)]
for(x in seq_along(Quant))
 trueCPUE[x,,y]<-suppressWarnings(as.numeric(unlist(strsplit(Quant[x],split=" "))))

tmp<-grep("true BMSY",Data[[y]])
trueBMSY[y]<-suppressWarnings(as.numeric(Data[[y]][(tmp+1)])[1])

tmp<-grep("est BMSY",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+out$OM$Nsim)]
for(x in seq_along(Quant))
 estBMSY[x,,y]<-suppressWarnings(as.numeric(unlist(strsplit(Quant[x],split=" "))))

tmp<-grep("true FMSY",Data[[y]])
trueFMSY[y]<-suppressWarnings(as.numeric(Data[[y]][(tmp+1)])[1])

tmp<-grep("est FMSY",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+out$OM$Nsim)]
for(x in seq_along(Quant))
 estFMSY[x,,y]<-suppressWarnings(as.numeric(unlist(strsplit(Quant[x],split=" "))))

 tmp<-grep("est total allowable catch",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+out$OM$Nsim)]
for(x in seq_along(Quant))
 estTAC[x,,y]<-suppressWarnings(as.numeric(unlist(strsplit(Quant[x],split=" "))))

tmp<-grep("true total allowable catch",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+out$OM$Nsim)]
for(x in seq_along(Quant)) {
 trueTAC[x,,y]<-suppressWarnings(as.numeric(unlist(strsplit(Quant[x],split=" "))))
 trueTAC[x,,y][trueTAC[x,,y]%in%0]<-NA
}

}

png(file.path(MSEdir,"plots",paste0("ProductionFits_",paste(CTLNameList,sep="_",collapse=""),".png")),res=500,width=6,height=4,units='in')
par(mfcol=c(1,length(CTLNameList)),mar=c(.1,.1,.1,.1),oma=c(4,6,1,1))

for(y in seq_along(CTLNameList))
{
#  boxplot(trueCPUE[,out$OM$start_assessment:out$OM$SimYear-1,y],type="l",ylim=c(0,max(trueCPUE,na.rm=T)),
#  las=1,xaxt='n',ylab='',yaxt='n')
#   input<-trueCPUE[,out$OM$start_assessment:out$OM$SimYear-1,y]
# for(x in 1:nrow(estCPUE))
#   lines(estCPUE[x,out$OM$start_assessment:out$OM$SimYear-1,y],col='#ff000033')
#  use_ylim<-c(0,max(trueCPUE,estCPUE,na.rm=T))
#  axis(1)
#  if(y==1)
#  {
#    #PolygonPlots(Truth=trueCPUE[,,y],Estimated=estCPUE[,,y],SimYear=out$OM$SimYear,Nsim=out$OM$Nsim,ylimIN = use_ylim)
#    mtext(side=3,plotNames[y],cex=.7)
#    axis(side=2,las=1)
#    mtext(side=2,"Biomass",line=4,cex=.9)
#    legend("topright",col=c(1,2,"#0000ff99","#00800077"),pch=c(15,NA,NA,NA),lty=c(NA,1,1,2),
#           legend=c("Observations","Estimates","True BMSY","Estimated BMSY range"),bty='n',cex=.7)
#  }
#  if(y>1)
#    #PolygonPlots(Truth=trueCPUE[,,y],Estimated=estCPUE[,,y],SimYear=out$OM$SimYear,Nsim=out$OM$Nsim,ylimIN = use_ylim)
#  mtext(side=3,plotNames[y],cex=.7)
#  abline(h=trueBMSY[y],col="#0000ff99",lty=1)
#  abline(h=max(estBMSY[,,y],na.rm=T),col="#00800099",lty=2)
#  abline(h=min(estBMSY[,,y],na.rm=T),col="#00800099",lty=2)
#  
# }
# par(mfrow=c(1,1))
# mtext(side=1,"Year",line=3,cex=.9)
# dev.off()
 use_ylim<-c(0,max(trueCPUE,estCPUE,na.rm=T))
  if(y==1)
  {
    PolygonPlots(Truth=trueCPUE[,,y],Estimated=estCPUE[,,y],SimYear=out$OM$SimYear,Nsim=out$OM$Nsim,ylimIN = use_ylim)
    axis(side=2,las=1)
    mtext(side=2,"Biomass",line=5,cex=.9)
    legend("topright",col=c(1,2,"#0000ff99","#00800077"),lty=c(1,1,1,2),
           legend=c("Observations","Estimates","True BMSY","Estimated BMSY range"),bty='n',cex=.7)
  }
  if(y>1) 
   PolygonPlots(Truth=trueCPUE[,,y],Estimated=estCPUE[,,y],SimYear=out$OM$SimYear,Nsim=out$OM$Nsim,ylimIN = use_ylim)
 mtext(side=3,plotNames[y],cex=.7)
 abline(h=trueBMSY[y],col="#0000ff99",lty=1)
 abline(h=max(estBMSY[,,y],na.rm=T),col="#00800099",lty=2)
 abline(h=min(estBMSY[,,y],na.rm=T),col="#00800099",lty=2)
 axis(side=1)

# use_ylim<-c(0,max(trueTAC,estTAC,na.rm=T))
#  
#  if(y==1)
#  {
#    PolygonPlots(Truth=trueTAC[,,y],Estimated=estTAC[,,y],SimYear=out$OM$SimYear,Nsim=out$OM$Nsim,ylimIN = use_ylim,bottom=T)
#    mtext(side=2,"Total allowable catch",line=4.5,cex=.9)
#    legend('topleft',col=c(1,2),pch=c(15,NA),lty=c(NA,1),legend=c("True","Estimated"),bty='n',cex=.7)
#  }
#  if(y>1)
#   {
#    PolygonPlots(Truth=trueTAC[,,y],Estimated=estTAC[,,y],SimYear=out$OM$SimYear,Nsim=out$OM$Nsim,ylimIN = use_ylim)
#    axis(side=1)
#   }
# 
 }
par(mfcol=c(1,1),xpd=NA)
  mtext(side=1,"Year",line=2.5,cex=.9)
dev.off()

png(file.path(MSEdir,"plots",paste0("ProductionRefPoints_",paste(CTLNameList,sep="_",collapse=""),".png")),res=500,width=5,height=4.5,units='in')
RelativeErrorBMSY<-list(list())
RelativeErrorFMSY<-list(list())
RelativeErrorTAC<-list(list())

for(y in seq_along(CTLNameList))
{
  temp<-sweep(estBMSY[,,y],MAR=2,trueBMSY[y],FUN="-")
  RelativeErrorBMSY[[y]]<-sweep(temp,MAR=2,trueBMSY[y],FUN="/")
  temp<-sweep(estFMSY[,,y],MAR=2,trueFMSY[y],FUN="-")
  RelativeErrorFMSY[[y]]<-sweep(temp,MAR=2,trueFMSY[y],FUN="/")
  RelativeErrorTAC[[y]]<-(estTAC[,,y]-trueTAC[,,y])/trueTAC[,,y]
}
temp<-c(unlist(RelativeErrorBMSY),unlist(RelativeErrorFMSY),unlist(RelativeErrorTAC))
temp[!is.finite(temp)]<-NA
yUp   <-quantile(temp,na.rm=T,probs=c(0.975))
if(yUp>10) yUp<-5
ydown <-quantile(temp,na.rm=T,probs=c(0.025))
if(ydown<(-10)) ydown<-(-5)
use_ylim<-c(ydown,yUp)

if((SimYear-out$OM$InitYear)>2) {
  par(mfcol=c(3,length(CTLNameList)),mar=c(.1,.1,.1,.1),oma=c(4,6,1,1))
  for(y in seq_along(CTLNameList))
  {
    inShape<-apply(RelativeErrorBMSY[[y]][,1:(ncol(RelativeErrorBMSY[[y]]))],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)
    plot(-100000000000,ylim=c(ydown,yUp),las=1,ylab="",xlab="Year",xaxt='n',xlim=c(out$OM$InitYear,out$OM$SimYear),yaxt='n')

    polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(inShape[1,1:SimYear-1],rev(inShape[4,1:SimYear-1])),col='darkgrey',border=F)
    polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(inShape[2,1:SimYear-1],rev(inShape[3,1:SimYear-1])),col='lightgrey',border=F)
    
  
    mtext(side=3,plotNames[y],cex=.7)
  
    abline(h=0,lty=2)
  
    if(y==1)
    {
      axis(side=2,las=1)
      mtext(side=2,"Relative error",line=2.5,cex=.7)
      mtext(side=2,"Target biomass",line=3.5,cex=.7)
    }
  
  # Plot FMSY
    inShape<-apply(RelativeErrorFMSY[[y]][,1:(ncol(RelativeErrorFMSY[[y]]))],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)
    
    plot(-100000000000,ylim=c(ydown,yUp),las=1,ylab="",xlab="Year",xlim=c(out$OM$InitYear,out$OM$SimYear),yaxt='n',xaxt='n')
    polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(inShape[1,1:SimYear-1],rev(inShape[4,1:SimYear-1])),col='darkgrey',border=F)
    polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(inShape[2,1:SimYear-1],rev(inShape[3,1:SimYear-1])),col='lightgrey',border=F)
    axis(side=1)
    
    abline(h=0,lty=2)
    if(y==1)
    {
    # axis(side=1)
      axis(side=2,las=1)
      mtext(side=2,"Relative error",line=2.5,cex=.7)
      mtext(side=2,"Target fishing mortality",line=3.5,cex=.7)
    }
    
    # Plot TAC
    inShape<-apply(RelativeErrorTAC[[y]][,1:(ncol(RelativeErrorTAC[[y]]))],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)
    
    plot(-100000000000,ylim=c(ydown,yUp),las=1,ylab="",xlab="Year",xlim=c(out$OM$InitYear,SimYear),yaxt='n')
    polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(inShape[1,1:SimYear-1],rev(inShape[4,1:SimYear-1])),col='darkgrey',border=F)
    polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(inShape[2,1:SimYear-1],rev(inShape[3,1:SimYear-1])),col='lightgrey',border=F)
    
    
    abline(h=0,lty=2)
    if(y==1)
    {
    # axis(side=1)
     axis(side=2,las=1)
     mtext(side=2,"Relative error",line=2.5,cex=.7)
     mtext(side=2,"Total Allowable Catch",line=3.5,cex=.7)
    }
}


if((SimYear-out$OM$InitYear)==1) {
  REBMSY <- REFMSY <- RETAC <- rep(NA,out$OM$Nsim)

  for(y in seq_along(CTLNameList)) {
      REBMSY <- cbind(REBMSY,RelativeErrorBMSY[[y]][,(out$OM$InitYear+1):(out$OM$SimYear)])
      REFMSY <- cbind(REFMSY,RelativeErrorFMSY[[y]][,(out$OM$InitYear+1):(SimYear)])
      RETAC <-cbind(RETAC,RelativeErrorTAC[[y]][,(out$OM$InitYear+1):(SimYear)])
  }

  REBMSY <- REBMSY[,-1]
  REFMSY <- REFMSY[,-1]
  RETAC <-  RETAC[,-1]

  inYlim<-c(min(REBMSY,REFMSY,RETAC,na.rm=T),max(REBMSY[REBMSY<6],
                                                 REFMSY[REFMSY<6],
                                                 RETAC[RETAC<6],na.rm=T))
  
  par(mfrow=c(3,1))
  
  boxplot(REBMSY,xaxt='n',ylim=inYlim,las=1)
  legend("topleft",c("(a) Target Biomass"),bty='n')
  mtext(side=2,"Relative Error",line=3,cex=.7)
  abline(h=0,lty=2)
  
  boxplot(REFMSY,xaxt='n',ylim=inYlim,las=1)
  legend("topleft",c("(b) Target Fishing Mortality"),bty='n')
  mtext(side=2,"Relative Error",line=3,cex=.7)
  abline(h=0,lty=2)
  
  boxplot(RETAC,xaxt='n',ylim=inYlim,las=1)
  legend("topleft",c("(c) Total Allowable Catch"),bty='n')
  mtext(side=2,"Relative Error",line=3,cex=.7)
  abline(h=0,lty=2)
  text(seq_along(plotNames), par("usr")[3]-.05, labels = plotNames, srt = 45, adj = c(1.1,1.1), xpd = NA, cex=1.1)

}

par(mfcol=c(1,1),xpd=NA)
if((SimYear-out$OM$InitYear)>=2) mtext(side=1,"Year",line=2.5,cex=.9)
dev.off()

#RefPointsError <- cbind(seq_along(CTLNameList),RelativeErrorFMSY[,yrs])
#yrs <- (out$OM$InitYear+1):out$OM$SimYear
#names(RefPointsError) <- c("OM","RefPoint","Year")

}

