#' Plot Life History
#' 
#' Plots life history based on inputs from CTL file, used in \code{\link{GeMS}}
#' 
#' @param LenAtAgeN Matrix of lengths at age for northern population
#' @param LenAtAgeS Matrix of lengths at age for southern population
#' @param matureN Proportions mature for northern population
#' @param matureS Proportions mature for southern population
#' @param vulnN Vulnerability matrix for northern population
#' @param vulnS Vulnerability matrix for southern population
#' @param survSelN Survey selectivity matrix for northern population
#' @param survSelS Survey selectivity matrix for southern population
#' @param WeightAtAgeN Weight-at-age matrix for northern population
#' @param WeightAtAgeS Weight-at-age matrix for southern population
#' @param MovementN Proportion of emigration to the southern population
#' @param MovementS Proportion of immigration into the northern population
#' @param NatMs Natural Mortality for southern population
#' @param NatMn Natural Mortality for northern population
#' @param VirBioN Initial virgin biomass for northern population
#' @param VirBioS Initial virgin biomass for southern population
#' @param RzeroN Virgin recruitment parameter for northern population
#' @param steepnessN Steepness parameter for northern population
#' @param steepnessS Steepness parameter for southern population
#' @param RzeroS Virgin recruitment parameter for southern population
#' @param RecErrS Recruitment deviations
#' @param sigmaRn Recruitment variability for northern population
#' @param sigmaRs Recruitment variability for southern population
#' @param HistoricalFn Vector of historical fishing mortalities for northern population
#' @param HistoricalFs Vector of historical fishing mortalities for southern population
#' @param SimYear Year of simulation
#' @param MaxAge Maximum age
#' 
#' @return Plots plots plots
#'
#' @export
#'
PlotLifeHistory<-function(LenAtAgeN,LenAtAgeS,matureN,matureS,vulnN,vulnS,survSelN,survSelS,WeightAtAgeN,
	WeightAtAgeS,MovementN,MovementS,NatMs,NatMn,VirBioN,VirBioS,RzeroN,RecErrN,steepnessN,steepnessS,
	RzeroS,RecErrS,sigmaRn,sigmaRs,HistoricalFn,HistoricalFs,SimYear,MaxAge)
  {
 plot(LenAtAgeN[1,],type="l",las=1,ylim=c(0,max(LenAtAgeN)))
 for(x in 2:SimYear)
  lines(LenAtAgeN[x,])
 mtext(side=3,"Length at Age (N)",cex=.7)

 plot(LenAtAgeS[1,],type="l",las=1,ylim=c(0,max(LenAtAgeS)))
 for(x in 2:SimYear)
  lines(LenAtAgeS[x,])
 mtext(side=3,"Length at Age (S)",cex=.7)

 plot(matureN[1,],las=1,type="l",ylim=c(0,max(matureN)))
 for(x in 2:SimYear)
  lines(matureN[x,])
 mtext(side=3,"Maturity at age (N)",cex=.7)

 plot(matureS[1,],las=1,type="l",ylim=c(0,max(matureS)))
 for(x in 2:SimYear)
  lines(matureS[x,])
 mtext(side=3,"Maturity at age (S)",cex=.7)

 plot(vulnN[1,],las=1,type="l",ylim=c(0,max(vulnN)))
 for(x in 2:SimYear)
  lines(vulnN[x,])
 mtext(side=3,"Fishery selectivity at age(N)",cex=.7)

 plot(vulnS[1,],las=1,type="l",ylim=c(0,max(vulnS)))
 for(x in 2:SimYear)
  lines(vulnS[x,])
 mtext(side=3,"Fishery selectivity at age (S)",cex=.7)

 plot(survSelN[1,],las=1,type="l",ylim=c(0,max(survSelN)))
 for(x in 2:SimYear)
  lines(survSelN[x,])
 mtext(side=3,"Index selectivity at age(N)",cex=.7)

 plot(survSelS[1,],las=1,type="l",ylim=c(0,max(survSelS)))
 for(x in 2:SimYear)
  lines(survSelS[x,])
 mtext(side=3,"Index selectivity at age (S)",cex=.7)

 plot(WeightAtAgeN[1,],las=1,type="l",ylim=c(0,max(WeightAtAgeN)))
 for(x in 2:SimYear)
  lines(WeightAtAgeN[x,])
 mtext(side=3,"Weight at age (N)",cex=.7)

 plot(WeightAtAgeS[1,],las=1,type="l",ylim=c(0,max(WeightAtAgeS)))
 for(x in 2:SimYear)
  lines(WeightAtAgeS[x,])
 mtext(side=3,"Weight at age (S)",cex=.7)

 plot(MovementN[1,],las=1,type="l",ylim=c(0,max(MovementN)))
 for(x in 2:SimYear)
  lines(MovementN[x,])
  mtext(side=3,"Movement at age (N)",cex=.7)

 plot(MovementS[1,],las=1,type="l",ylim=c(0,max(MovementS)))
 for(x in 2:SimYear)
  lines(MovementS[x,])
  mtext(side=3,"Movement at age (S)",cex=.7)

 plot(NatMs,las=1,type="l",ylim=c(0,max(NatMs)))
 lines(NatMn,lty=2)
 legend("bottomleft",bty='n',lty=c(1,2),legend=c("South","North"))
  mtext(side=3,"Natural mortality by year",cex=.7)

ssb<-seq(1,VirBioN,VirBioN/100)
record<-ssb
  for(x in seq_along(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,],MaxAge=MaxAge,sigmaRin=sigmaRn[1])
   }
plot(record/10000~ssb,type='l',las=1)
mtext(side=3,"SSB vs Rec (N)",cex=.7)
ssb<-seq(1,VirBioS,VirBioS/100)
record<-ssb
  for(x in seq_along(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,],MaxAge=MaxAge,sigmaRin=sigmaRs[1])
   }
plot(record/10000~ssb,type='l',las=1)
mtext(side=3,"SSB vs Rec (S)",cex=.7)


plot(HistoricalFn,type='l',las=1)
lines(HistoricalFs,lty=2,col=2)
mtext(side=3,"Historical F",cex=.7)
 }
 