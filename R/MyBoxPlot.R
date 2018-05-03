#' Customised Box Plot
#' Used within the \code{\link{GeMSplots}} function, not quite implemented.
#' 
#' @param Inputs Matrix of values to be plotted
#' @param Title Title of the plot
#' @param bottomplot 1/0 if figure is at the bottom of plot, includes the legend
#' @param refVal Reference value to set maximum y limit
#' 
#' @return Boxplots
#'
#' @export
#' @examples
#' \dontrun{
#' MyBoxPlot(FMSY[,InitYear:SimYear],"FMSY")
#' }
#' 
MyBoxPlot<-function(Inputs,Title,bottomplot=0,refVal=0)
{
 maxVal<-max(Inputs,na.rm=T)
 minVal<-min(Inputs,na.rm=T)

 ylimTop<-maxVal
 if(ylimTop<refVal)
  ylimTop<-refVal

 ylimBot<-minVal
 if(ylimBot>refVal)
  ylimBot<-refVal

 coltemp<-apply(Inputs,2,median,na.rm=T)
 colIn<-coltemp
 colIn[coltemp>refVal]<-"green"
 colIn[coltemp<=refVal]<-'red'

 if(bottomplot==0)
 boxplot(Inputs,ylim=c(ylimBot,ylimTop),las=1,xaxt='n',col=colIn)
 if(bottomplot==1)
 boxplot(Inputs,ylim=c(ylimBot,ylimTop),las=1,col=colIn)
 abline(h=refVal,lty=2)
 legend("topleft",bty='n',legend=Title)
}