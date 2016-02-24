#' StockRecCheck checks stock recruitment
#'
#' @param rec recruitment
#' @param spbio spawning biomass
#' @param plotSR plot spawners per recruit
#'
#' @return list of checks
#' @export
StockRecCheck<-function(rec,spbio,plotSR=F)
{
  Stx		<-srStarts(rec~spbio,type="BevertonHolt",param=2)
  x		<-c(Stx$a,Stx$Rp,1)
  outs		<-optim(x,BevHolt)
  counter	<-1
  while(outs$convergence!=0 & counter<20)
  {
    x		<-outs$par
    outs		<-optim(x,BevHolt)
    counter	<-counter+1
  }

  Stx1<-srStarts(rec~spbio,type="Ricker",param=2)
  x1		<-c(Stx1$a,Stx1$b,1)
  outs1		<-optim(x1,Ricker)
  x1		<-outs1$par
  outs1		<-optim(x1,Ricker)
  x1		<-outs1$par
  outs1		<-optim(x1,Ricker)

  x2		<-c(mean(rec),sd(rec))
  outs2		<-optim(x2,StrLine)

  AvgAIC	<-outs2$value
  RickerAIC	<-outs1$value
  BHaic		<-outs$value

  AvgPar	<-outs2$par
  RickerPars	<-outs1$par
  BHpars	<-outs$par

  #===check parameters
  AICR		<-6+2*RickerAIC
  AICB		<-6+2*BHaic
  AICA		<-4+2*AvgAIC

  #==report decision about stock recruit relationship
  SRchoice	<-"Average"
  returnPar	<-AvgPar

  if(AICR<AICA-2 & AICR<AICB)
  {
    SRchoice	<-"Rick"
    returnPar	<-RickerPars
  }

  if(AICB<AICA-2 & AICB<AICR)
  {
    SRchoice	<-"BevH"
    returnPar	<-BHpars
  }

  if(plotSR==T)
  {
    predSpbio<-seq(1,max(spbio),max(spbio)/40)
    alpha<-abs(BHpars[1])
    beta<-abs(BHpars[2])
    plotPred<-alpha*predSpbio/(1+(predSpbio/beta))
    plot(rec~spbio,ylim=c(0,max(rec)),xlim=c(0,max(spbio)))
    lines(plotPred~predSpbio)

    alpha<-RickerPars[1]
    beta<-RickerPars[2]
    plotPred<-predSpbio*exp(alpha-beta*predSpbio)
    lines(plotPred~predSpbio,lty=2,col=2)


    abline(h=AvgPar[1],lty=4,col=4)

    legend("topright",bty='n',col=c(1,2,4),lty=c(1,2,4),legend=c(paste("BH ",round(AICB,2)),
                                                                 paste("Rick ", round(AICR,2)),paste("Avg ",round(AICA,2))))
  }
  list(SRchoice,returnPar)
}