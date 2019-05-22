#' Age-Structure Comparison
#' @description Compare and plot age-structured models
#' @description Used within the \code{\link{run_GeMS}} function, can also be run on its own.
#' 
#' @param out List output from \code{\link{ReadCTLfile}}
#' @param RetroPeels Number of years to peel for retrospective analysis
#' @param CTLNameList Vector of CTL file names
#' @param MSEdir Directory containing CTL files
#' @param plotNames Vector of the same length as CTLNameList with "prettier" names for plotting
#' @param Nruns Number of runs to plot in Population Processes plots (selectivity and growth)
#' @param plottiff Logical; make plots in TFF instead of PNG?
#' @param GradientTolerance Tolerance level for convergence; runs with gradient levels greater than or equal to this value will be discarded
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
AgeStructureComp<-function(out,RetroPeels=6,CTLNameList,MSEdir,
                            plotNames=CTLNameList,Nruns=out$OM$Nsim,
                            plottiff=F,GradientTolerance=1E-3)
{
  TakeRows<-(out$OM$SimYear-RetroPeels+1):out$OM$SimYear
  GradientSave<-array(dim=c(RetroPeels,out$OM$Nsim,length(CTLNameList)))
  
  for(x in seq_along(CTLNameList))
  {
    out<-ReadCTLfile(file.path(MSEdir,paste0(CTLNameList[x])))
    GradientSave[,,x]<-CheckGradient(CTLName=CTLNameList[x],out=out,MSEdir)[TakeRows,]
  }
  #ScenCols<-c("grey",as.numeric(seq(2,length(CTLNameList),1)))
  ScenCols<-colorspace::rainbow_hcl(length(CTLNameList))

  convRuns<-GradientSave[RetroPeels,,]<GradientTolerance

  if(!plottiff)png(file.path(MSEdir,"plots",paste0("ConvergenceRates_",paste(CTLNameList,sep="_",collapse=""),".png")),height=6,width=9,units='in',res=1200)
  if(plottiff)tiff(file.path(MSEdir,"plots",paste0("ConvergenceRates_",paste(CTLNameList,sep="_",collapse=""),".tiff")),height=2,width=3,units='in',res=1200)  

  par(oma=c(3,2,.5,.5),mar=c(4,1,3.5,.5))
  temp<-barplot(colSums(convRuns),main=paste0("Convergence Rates based on Gradients < ",GradientTolerance),xaxt='n')
  text(temp, par("usr")[3]-.05, labels = plotNames, srt = 45, adj = c(1.1,1.1), xpd = NA, cex=1.1)
  text(temp, colSums(convRuns)+2.5 ,.01*colSums(convRuns),cex=1.1,xpd=NA)
  dev.off()

  #==reference points
  TakeRows	<-(out$OM$InitYear+1):out$OM$SimYear
  B35save	<-array(dim=c(length(TakeRows),out$OM$Nsim,length(CTLNameList)))
  F35save	<-array(dim=c(length(TakeRows),out$OM$Nsim,length(CTLNameList)))
  OFLsave	<-array(dim=c(length(TakeRows),out$OM$Nsim,length(CTLNameList)))
  tOFLsave	<-array(dim=c(length(TakeRows),out$OM$Nsim,length(CTLNameList)))
  tF35save	<-array(dim=c(length(TakeRows),out$OM$Nsim,length(CTLNameList)))
  tB35save	<-array(dim=c(length(TakeRows),out$OM$Nsim,length(CTLNameList)))
  
  for(x in seq_along(CTLNameList))
  {
    temp<-CheckReferencePoints(CTLNameList[x],out,MSEdir)
    B35save[,,x]	<-temp[[1]][TakeRows,]
    F35save[,,x]	<-temp[[2]][TakeRows,]
    OFLsave[,,x]	<-temp[[3]][TakeRows,]
    tB35save[,,x]	<-temp[[4]][TakeRows,]
    tF35save[,,x]	<-temp[[5]][TakeRows]
    tOFLsave[,,x]	<-temp[[6]][TakeRows,]
  } 
  
  tB35save[tB35save<=0]<-1E-10
  tF35save[tF35save<=0]<-1E-10
  tOFLsave[tOFLsave<=0]<-1E-10

  BigB35<-matrix(nrow=out$OM$Nsim,ncol=length(CTLNameList))
  BigF35<-matrix(nrow=out$OM$Nsim,ncol=length(CTLNameList))
  BigOFL<-matrix(nrow=out$OM$Nsim,ncol=length(CTLNameList))
  if(out$OM$Nsim>1 & (out$OM$SimYear-out$OM$InitYear)>1) {
    for(x in seq_along(CTLNameList))
    {
      BigB35[,x]<-apply((B35save[,,x]-tB35save[,,x])/tB35save[,,x],2,median,na.rm=T)
      BigF35[,x]<-apply((F35save[,,x]-tF35save[,,x])/tF35save[,,x],2,median,na.rm=T)
      BigOFL[,x]<-apply((OFLsave[,,x]-tOFLsave[,,x])/tOFLsave[,,x],2,median,na.rm=T)
    }
  }
  
  if(out$OM$Nsim==1 | (out$OM$SimYear-out$OM$InitYear)==1) {
    for(x in seq_along(CTLNameList))
    {
      BigB35[,x]<-median((B35save[,,x]-tB35save[,,x])/tB35save[,,x],na.rm=T)
      BigF35[,x]<-median((F35save[,,x]-tF35save[,,x])/tF35save[,,x],na.rm=T)
      BigOFL[,x]<-median((OFLsave[,,x]-tOFLsave[,,x])/tOFLsave[,,x],na.rm=T)
    }
  }
  
  ReB35<-ifelse(convRuns,BigB35,NA)
  ReF35<-ifelse(convRuns,BigF35,NA)
  ReOFL<-ifelse(convRuns,BigOFL,NA)
  
  #=plotting output
  # GeMSplots(out=out,DrawDir=CTLNameList[x],MSEdir=MSEdir)
  MohnsRho<-array(dim=c(RetroPeels-1,out$OM$Nsim,length(CTLNameList)))
  SSBbias<-array(dim=c(RetroPeels,out$OM$Nsim,length(CTLNameList)))
  TrueOFL<-array(dim=c(RetroPeels,out$OM$Nsim,length(CTLNameList)))
  
  for(x in seq_along(CTLNameList))
  {
    out			<-ReadCTLfile(file.path(MSEdir,CTLNameList[x]))
    DrawDir		<-CTLNameList[x]
    RetroOuts		<-CheckRetro(RetroPeels=RetroPeels,DrawDir,PlotRetro=0,out=out,MSEdir)
    if(RetroPeels>1) {
      MohnsRho[,,x]	<-CalculateMohn(RetroOuts)

    }
    if(RetroPeels==1) {
      MohnsRho<-NA
    }
    SSBbias[,,x]  <-CalculateBias(RetroOuts)
  }
  
  BigMohn<-matrix(nrow=out$OM$Nsim,ncol=length(CTLNameList))
  if(RetroPeels==1) BigMohn<-NA
  if(out$OM$Nsim > 1 & (RetroPeels-1)>1) {
    for(x in seq_along(CTLNameList))
    {
      temp<-MohnsRho[,,x]
      BigMohn[,x]<-apply(temp,2,mean) 
    }
  }
  
  if(out$OM$Nsim == 1 | (RetroPeels-1) == 1) {
    for(x in seq_along(CTLNameList))
    {
      temp<-MohnsRho[,,x]
      BigMohn[,x]<-mean(temp) 
    }
  }

  BigBias<-matrix(nrow=out$OM$Nsim,ncol=length(CTLNameList))
  
  for(x in seq_along(CTLNameList))
    {
      temp<-SSBbias[,,x]

      if(out$OM$Nsim > 1 & RetroPeels > 1) {
        BigBias[,x]<-apply(temp,2,mean) 
      }
      if(out$OM$Nsim == 1 & RetroPeels > 1) {
        BigBias[,x]<-mean(temp) 
      }
      if(RetroPeels == 1) {
        BigBias[,x]<-temp 
      }
  }

  convBigMohn<-ifelse(convRuns,BigMohn,NA)
  convBigBias<-ifelse(convRuns,BigBias,NA)

  if(!plottiff)png(file.path(MSEdir,"plots",paste0("CompareRefPoints",paste(CTLNameList,sep="_",collapse=""),".png")),height=7,width=3.5,units='in',res=1200)
  if(plottiff)tiff(file.path(MSEdir,"plots",paste0("CompareRefPoints",paste(CTLNameList,sep="_",collapse=""),".tiff")),height=5,width=2.5,units='in',res=1200)
  par(mfrow=c(5,1),mar=c(.1,.1,.3,.1),oma=c(10,6,1,1))
  
  inYlim<-c(min(convBigMohn,convBigBias,ReB35,ReF35,ReOFL,na.rm=T),max(convBigMohn[convBigMohn<6],
                                                        convBigBias[convBigBias<6],ReB35[ReB35<6],
                                                        ReF35[ReF35<6],ReOFL[ReOFL<6],na.rm=T))
  boxplot(convBigMohn,col=ScenCols,ylim=inYlim,xaxt='n',las=1)
  legend("topleft",c("(a) Retrospective bias"),bty='n')
  mtext(side=2,"Mohn's Rho",line=3,cex=.7)
  abline(h=0,lty=2)
  
  boxplot(convBigBias,col=ScenCols,xaxt='n',ylim=inYlim,las=1)
  legend("topleft",c("(b) Spawning biomass"),bty='n')
  mtext(side=2,"Relative Error",line=3,cex=.7)
  abline(h=0,lty=2)
  boxplot(ReB35,col=ScenCols,ylim=inYlim,xaxt='n',las=1)
  legend("topleft",expression("(c) B"[35]),bty='n')
  mtext(side=2,"Relative Error",line=3,cex=.7)
  abline(h=0,lty=2)
  boxplot(ReF35,col=ScenCols,ylim=inYlim,xaxt='n',las=1)
  legend("topleft",expression("(d) F"[35]),bty='n')
  mtext(side=2,"Relative Error",line=3,cex=.7)
  abline(h=0,lty=2)
  boxplot(ReOFL,col=ScenCols,ylim=inYlim,xaxt='n',las=2)
  axis(1, at = seq_along(plotNames), labels = FALSE)
  text(seq_along(plotNames), par("usr")[3]-.05, labels = plotNames, srt = 45, adj = c(1.1,1.1), xpd = NA, cex=1.1)
  abline(h=0,lty=2)
  legend("topleft",c("(e) OFL"),bty='n')
  mtext(side=2,"Relative Error",line=3,cex=.7)
  dev.off()
  
  quants<-PullTimevary(out,MSEdir,CTLNameList,convRuns)
  
  if(!plottiff)png(file.path(MSEdir,"plots",paste0("ComparePopulationProcess_",paste(CTLNameList,sep="_",collapse=""),".png")),height=5,width=7.5,units='in',res=1200)
  if(plottiff)tiff(file.path(MSEdir,"plots",paste0("ComparePopulationProcess_",paste(CTLNameList,sep="_",collapse=""),".tiff")),height=4,width=6,units='in',res=1200)
  inmat<-matrix(c(1,1,2,2,3,3,
                  4,4,4,5,5,5),nrow=2,byrow=T)

  convRunNos<-ifelse(convRuns,row(convRuns),NA)
  runs2plot<-apply(convRunNos,2,function(x,runs) {
                                  if(runs<=length(x[!is.na(x)])) return(sample(x[!is.na(x)],runs))
                                  if(runs>length(x[!is.na(x)])) return(rep(x[!is.na(x)],length=runs))},
                                    runs=Nruns)

#  if(Nruns<out$OM$Nsim) cat(paste0("Plotting population processes from run numbers: ",paste(runs2plot,collapse=" "),".\n"))
  layout(inmat)
  par(mar=c(.1,.1,.1,.1),oma=c(4,5,4,4))
  
  #==RECRUITMENT
  input<-quants[[1]]
  tInput<-quants[[2]] 
  
  plot(-100000,xlim=c(1,dim(input)[2]),ylim=c(0,max(input,na.rm=T)),las=1,xaxt='n')
  for(x in seq_along(CTLNameList))
  {
    #color<-seq(1,length(CTLNameList)+1)
    color<-colorspace::rainbow_hcl(length(CTLNameList))
    incol<-adjustcolor(color,alpha.f=.2)
    for(y in runs2plot[,x]) 
      lines(input[y,,x],col=incol[x])
  } 
  medians<- apply(input,c(2,3),median,na.rm=T)
  for(x in seq_along(CTLNameList)) lines(medians[,x],col=color[x],cex=2)
  
  abline(v=out$OM$InitYear,lty=3)
  
  #==true R
  plotIn<-apply(tInput,2,median,na.rm=T)
  plotIn[1]<-NA
  lines(plotIn,col="black",lwd=1.5,lty=2)
  
  legend("bottomleft",bty='n',"(a) Recruitment")
  mtext(side=2,"Numbers",line=4,cex=.7)

  #==FISHING MORTALITY
  input<-quants[[3]]
  tInput<-quants[[4]]
  
  plot(-1000,xlim=c(1,dim(input)[2]),ylim=range(c(0,input,tInput),na.rm=T),las=1,xaxt='n',yaxt='n')
  axis(side=3)
  mtext(side=3,line=2,"Time")
  for(x in seq_along(CTLNameList))
  {
    color<-colorspace::rainbow_hcl(length(CTLNameList))
    incol<-adjustcolor(color,alpha.f=.2)
    temp<-input[,,x]
    for(y in runs2plot[,x]) 
      lines(input[y,,x],col=incol[x])
  } 
  medians<- apply(input,c(2,3),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=color[x],cex=2)
  
  abline(v=out$OM$InitYear,lty=3)
  
  #==true F
  lines(apply(tInput,2,median),col="black",lwd=1.5,lty=2)
  legend("bottomleft",bty='n',"(b) Fishing mortality")
  
  #==NATURAL MORTALITY
  input<-quants[[7]]
  tInput<-quants[[8]] 
  
  plot(-1000,xlim=c(1,dim(input)[2]),ylim=range(c(0,input,tInput),na.rm=T),las=1,xaxt='n',yaxt='n')
  axis(side=4,las=1)
  for(x in seq_along(CTLNameList))
  {
    #color<-seq(1,length(CTLNameList)+1)
    color<-colorspace::rainbow_hcl(length(CTLNameList))
    incol<-adjustcolor(color,alpha.f=.2)
    for(y in runs2plot[,x]) 
      lines(input[y,,x],col=incol[x])
  } 
  medians<- apply(input,c(2,3),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=color[x],cex=2)
  
  #==true M
  #dim(tInput)
  lines(tInput,col="black",lwd=1.5,lty=2)
  
  abline(v=out$OM$InitYear,lty=3)
  
  legend("bottomleft",bty='n',"(c) Natural mortality")
  mtext(side=4,"rate (/year)",line=2,cex=.7)
  
  #==SELECTIVITY
  input<-quants[[9]]
  tInput<-quants[[10]] 

  plot(-1000,xlim=c(1,which(input[1,,1,1]>0.99)[3]),ylim=c(0,max(input,na.rm=T)),las=1)
  for(x in seq_along(CTLNameList))
  {
    #color<-seq(1,length(CTLNameList)+1)
    color<-colorspace::rainbow_hcl(length(CTLNameList))
    incol<-adjustcolor(color,alpha.f=.2)
    tCol<-rgb(0,0,0,0.2)
    for(y in 1:dim(input)[1])
    {
      #lines(tInput[y,,x],col=tCol)
      for(z in runs2plot)
        lines(input[y,,z,x],col=incol[x])
    }   
  } 
  
  medians<- apply(input,c(2,4),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=color[x],cex=2)
  
  lines(tInput[1,,1],col="black",lwd=1.5,lty=2)
  lines(tInput[out$OM$SimYear,,1],col="black",lwd=1.5,lty=2)
  mtext(side=2,"Selectivity",line=2.3,cex=.7)
  legend("topleft",bty='n',"(d) Fishery selectivity")
  
  #==GROWTH
  input<-quants[[5]]
  tInput<-quants[[6]] 

  plot(-1000,xlim=c(1,out$OM$MaxAge),ylim=c(0,max(input,na.rm=T)),las=1,yaxt='n')
  axis(side=4,las=1)
  for(x in seq_along(CTLNameList))
  {
    #color<-seq(1,length(CTLNameList)+1)
    color<-colorspace::rainbow_hcl(length(CTLNameList))
    incol<-adjustcolor(color,alpha.f=.2)
    tCol<-rgb(0,0,0,0.2)
    for(y in 1:dim(input)[1])
    {
      #lines(tInput[y,,x],col=tCol)
      for(z in seq_along(runs2plot))
        lines(input[y,,z,x],col=incol[x])
    }   
  } 
  
  medians<- apply(input,c(2,4),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=color[x],cex=2)
  mtext(side=4,"Length (cm)",line=2.3,cex=.7)
  
  #==TRUE GROWTH
  lines(tInput[1,,1],col="black",lwd=1.5,lty=2)
  lines(tInput[out$OM$SimYear,,1],col="black",lwd=2,lty=2)
  legend("topleft",bty='n',"(e) Length at age")
  
  legend("bottomright",bty='n',col=c(1,color,"black"),lty=c(3,rep(1,length(ScenCols)),2),
         legend=c("Projection begins",plotNames,"Truth"),lwd=2)
  mtext(side=1,outer=T,"Age",line=2.3)
  
  dev.off()
  
  
  #================================
  # Plot the estimates of spawning biomass
  #====================================
  if(!plottiff)png(file.path(MSEdir,"plots",paste0("Spbio_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".png")),height=5,width=7.5,units='in',res=1200)
  if(plottiff)tiff(file.path(MSEdir,"plots",paste0("Spbio_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".tiff")),height=5,width=7.5,units='in',res=1200)
  plot_stuff_inside(input1=quants[[14]][,out$OM$start_assessment:out$OM$SimYear,],input2=quants[[15]][runs2plot,out$OM$start_assessment:out$OM$SimYear,],title="Biomass",CTLNameList=CTLNameList,plotNames=plotNames,nSim=out$OM$Nsim)
  dev.off()
  if(!plottiff)png(file.path(MSEdir,"plots",paste0("SurveyInd_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".png")),height=5,width=7.5,units='in',res=1200)
  if(plottiff)tiff(file.path(MSEdir,"plots",paste0("SurveyInd_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".tiff")),height=5,width=7.5,units='in',res=1200)
  plot_stuff_inside(input1=quants[[17]],input2=quants[[16]][runs2plot,,],title="Survey",CTLNameList=CTLNameList,plotNames=plotNames,nSim=out$OM$Nsim)
  dev.off()
  if(!plottiff)png(file.path(MSEdir,"plots",paste0("CPUE_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".png")),height=5,width=7.5,units='in',res=1200)
  if(plottiff)tiff(file.path(MSEdir,"plots",paste0("CPUE_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".tiff")),height=5,width=7.5,units='in',res=1200)
  plot_stuff_inside(input1=quants[[19]],input2=quants[[18]][runs2plot,,],title="CPUE",CTLNameList=CTLNameList,plotNames=plotNames,nSim=out$OM$Nsim)
  dev.off()
  if(!plottiff)png(file.path(MSEdir,"plots",paste0("Catch_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".png")),height=5,width=7.5,units='in',res=1200)
  if(plottiff)tiff(file.path(MSEdir,"plots",paste0("Catch_true_vs_est",paste(CTLNameList,sep="_",collapse=""),".tiff")),height=5,width=7.5,units='in',res=1200)
  plot_stuff_inside(input1=quants[[21]],input2=quants[[20]][runs2plot,,],title="Catch",CTLNameList=CTLNameList,plotNames=plotNames,nSim=out$OM$Nsim)
  dev.off()
 

  
  
  #==================================
  # relative error in OFL over time
  # NEED TO GET THIS ONCE FIGURED OUT AT SOME POINT
  #====================================
  OFLrel<-0
  if(OFLrel>0)
  {
  InitYear<-out$OM$InitYear
  SimYear<-out$OM$SimYear
  
  REofl<-rep(list(list()))
  for(x in seq_along(SaveTrueOFL))
    REofl[[x]]<-(SaveTrueOFL[[x]][,(InitYear+1):(SimYear-1),]-SaveTrueOFL[[x]][,(InitYear+1):(SimYear-1),])/SaveTrueOFL[[x]][,(InitYear+1):(SimYear-1),]
  
  emNames<-c("Base","M vary","Growth vary","Sel vary")
  omNames<-c("Base","M vary","Sel vary","Growth vary")
  hexDec<-c("#ff000055","#00ff0055","#0000ff55")
  par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
  for(x in seq_along(REofl))
  {
    #boxplot(REofl[[x]][,,1],col="darkgrey",ylim=c(min(REofl[[x]]),max(REofl[[x]])),xaxt='n',las=1)
    
    plot(NA,ylim=c(min(unlist(REofl)),max(unlist(REofl))),xlim=c(1,9),xaxt='n',las=1)
    plotQuant<-REofl[[x]][,,1]
    sortQuant<-apply(plotQuant,2,sort)
    #polygon(x=c(seq(1,9),seq(9,1)),y=c(sortQuant[3,],rev(sortQuant[18,])),col="#cfcfcf99",border=F)
    #lines(sortQuant[3,],col='darkgrey',lty=2)
    #lines(sortQuant[18,],col='darkgrey',lty=2)
    lines(apply(sortQuant,2,median),col='darkgrey',lwd=2)
    for(y in 2:4)
    {
      plotQuant<-REofl[[x]][,,y]
      sortQuant<-apply(plotQuant,2,sort)
      # lines(sortQuant[3,],col=y,lty=2)
      #lines(sortQuant[18,],col=y,lty=2)
      lines(apply(sortQuant,2,median,na.rm=T),col=y,lwd=2)
      
      #polygon(x=c(seq(1,9),seq(9,1)),y=c(sortQuant[3,],rev(sortQuant[18,])),col=hexDec[y-1],border=F)
      #lines(apply(sortQuant,2,median),col=hexDec[y-1])
    } 
    abline(h=0,lty=2)
    legend('topright',bty='n',omNames[x])
    if(x==1)
      legend('topleft',emNames,col=c('darkgrey',2,3,4),pch=15,bty='n')
  }
  axis(side=1)
  mtext(outer=T,side=2,"Relative error in TAC",line=2.5)
  mtext(outer=T,side=1,"Projection year",line=2)
  
  dev.new()
  boxplot(SaveTrueOFL[[1]][,(InitYear+1):(SimYear-1),],col="#00ff0033")
  boxplot(SaveEstOFL[[1]][,(InitYear+1):(SimYear-1),1],add=T,col="#ff000033")
  
  dev.new()
  par(mfrow=c(4,2),mar=c(.1,.1,.1,.1),oma=c(4,5,1,4))
  for(x in seq_along(REofl))
  {
    plot(NA,ylim=c(min(unlist(SaveEstOFL),na.rm=T),max(unlist(SaveEstOFL),na.rm=T)),xlim=c(1,9),xaxt='n',las=1)
    plotQuant<-SaveEstOFL[[x]][,51:60,1]
    sortQuant<-apply(plotQuant,2,sort)
    lines(apply(sortQuant,2,median),col='darkgrey',lwd=2)
    for(y in 2:4)
    {
      plotQuant<-SaveEstOFL[[x]][,51:60,y]
      sortQuant<-apply(plotQuant,2,sort)
      lines(apply(sortQuant,2,median),col=y,lwd=2)
    } 
    abline(h=0,lty=2)
    legend('topright',bty='n',paste("(",letters[2*x-1],") ",omNames[x],sep=""))
    
    if(x==1)
      legend('topleft',emNames,col=c('darkgrey',2,3,4),pch=15,bty='n')
    if(x==4)
      axis(side=1)
    #==plot the status true
    temp<-lapply(SaveSpBio,'[',,51:60,)
    tempSB<-t(as.matrix(temp[[x]][,,1]))
    tempB35<-tB35save[1,1,1,x]
    bigStat<-tempSB/tempB35
    plot(NA,ylim=c(.75,1.75),xlim=c(1,9),xaxt='n',las=1,yaxt='n',las=1)
    axis(side=4,las=1)
    plotQuant<-bigStat
    lines(apply(plotQuant,1,median),col='darkgrey',lwd=2)
    for(y in 2:4)
    {
      tempSB<-t(as.matrix(temp[[x]][,,y]))
      tempB35<-tB35save[1,1,y,x]
      bigStat<-tempSB/tempB35
      plotQuant<-bigStat
      lines(apply(plotQuant,1,median),col=y,lwd=2)
    } 
    legend('topright',bty='n',paste("(",letters[2*x],") ",sep=""))
    abline(h=1,lty=2)
    if(x==4)
      axis(side=1)
  }
  mtext(side=4,outer=T,text=expression("B/B"[35]),line=2.5)
  mtext(side=2,outer=T,"Catch",line=3.25)
  mtext(side=1,outer=T,"Projection year",line=2.5)
  }
}
