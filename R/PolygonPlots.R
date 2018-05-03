#' Plots polygons
#' 
#' @return Plots plots plots
#'
#' @export
#' 
PolygonPlots<-function(Truth=NA,Estimated=NA,Observed=NA,AddLegend=F,bottom=F,quantity=NA,SimYear=NA,Nsim=NA,ylimIN=NA)
{
 EstShape<-apply(Estimated[,1:(SimYear)],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)
 TrueShape<-apply(Truth[,1:(SimYear)],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)

 DataColor<-"#6699ff11"
 EstimateColor25<-"darkgrey"
 EstimateColor5<-"lightgrey"
 TruthColor<-"#FF000033"

 if(!is.na(ylimIN[1]))
   use_ylim<-ylimIN
 if(is.na(ylimIN[1]))
   use_ylim<-c(0,max(EstShape,TrueShape,Observed,na.rm=t))
 
if(bottom==F)
 plot(-100000000000,ylim=use_ylim,xlim=c(1,SimYear),las=1,ylab="",xlab="Year",xaxt='n',yaxt='n')
if(bottom==T)
 plot(-100000000000,ylim=use_ylim,xlim=c(1,SimYear),las=1,ylab="",xlab="Year")

if(!is.na(quantity))
 legend("bottomleft",bty='n',legend=quantity)

 if(length(Observed)>1)
 {
 if(AddLegend==T)
  legend("topright",bty='n',pch=c(15,16,NA),lty=c(NA,NA,1),col=c(EstimateColor5,DataColor,TruthColor),legend=c("Estimates","Observed Data","Truth"))
 for(z in 1:Nsim)
  points(Observed[z,]~seq(1,(SimYear-1)),col=DataColor,pch=16)
 }
 polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(TrueShape[1,1:SimYear-1],rev(TrueShape[4,1:SimYear-1])),col=EstimateColor5,border=F)
 polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(TrueShape[2,1:SimYear-1],rev(TrueShape[3,1:SimYear-1])),col=EstimateColor25,border=F)
 if(AddLegend==T&length(Observed)<1)
 legend("topright",bty='n',pch=c(15,NA),lty=c(NA,1),col=c(EstimateColor5,TruthColor),legend=c("Truth","Estimates"))

 for(z in 1:Nsim)
  lines(Estimated[z,]~seq(1,(SimYear)),col=TruthColor,lty=1)
}