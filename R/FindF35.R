#' Calculates F35 for \code{\link{GeMS}}
#' 
#' @param MaxAge Maximum age
#' @param vulnN Vulnerability matrix for northern population
#' @param vulnS Vulnerability matrix for southern population
#' @param NatMn Natural Mortality for northern population
#' @param NatMs Natural Mortality for southern population
#' @param matureN Proportions mature for northern population
#' @param matureS Proportions mature for southern population
#' @param WeightAtAgeN Weight-at-age matrix for northern population
#' @param WeightAtAgeS Weight-at-age matrix for southern population
#' @param steepnessN Steepness parameter for northern population
#' @param steepnessS Steepness parameter for southern population
#' @param RzeroN Virgin recruitment parameter for northern population
#' @param RzeroS Virgin recruitment parameter for southern population
#' @param LenAtAgeN Matrix of lengths at age for northern population
#' @param LenAtAgeS Matrix of lengths at age for southern population
#' @param inRec Historical recruitment values
#' @param RefYear Year of forecast
#' @param fmortSin Historical fishing mortality
#' @param MovementN Proportion of emigration to the southern population
#' @param MovementS Proportion of immigration into the northern population
#'
#' @return Vector of historical and calculated future Fs.
#' 
#' @export

FindF35<-function(MaxAge,vulnN,vulnS,NatMn,NatMs,matureN,matureS,WeightAtAgeN,
 WeightAtAgeS,steepnessN,steepnessS,RzeroN,RzeroS,LenAtAgeN,LenAtAgeS,inRec,RefYear,
 fmortSin=0,MovementN,MovementS)
{
tempNn		<-matrix(ncol=MaxAge,nrow=100)
tempNn[1,]		<-initialN(Rzero=RzeroN[RefYear],NatM=NatMn[RefYear],inAge=MaxAge)
tempNs		<-matrix(ncol=MaxAge,nrow=100)
tempNs[1,]		<-initialN(Rzero=RzeroS[RefYear],NatM=NatMs[RefYear],inAge=MaxAge)
tempCatchN		<-rep(0,100)
tempCatchS		<-rep(0,100)
tempRecN		<-rep(0,100)
tempRecS		<-rep(0,100)
tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=100)
tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=100)

Target<-0.35

Bzero     <-ProjPopDym(fmortN=0,MaxAge=MaxAge,vulnN=vulnN,
                      vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
                      WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
                      steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,
                      LenAtAgeS=LenAtAgeS,ProxyRec=inRec,RefYear=RefYear,
                      MovementN=MovementN,MovementS=MovementS)[[2]]

maxF	<-3
minF	<-0.01
for(k in 1:20)
{
inFn	<-(maxF+minF)/2
inFs	<-fmortSin
for (j in 2:100)
{
 for (i in 2:(MaxAge-1))
  {
   tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-inFn*vulnN[RefYear,i-1])*exp(-NatMn[RefYear])
   tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-inFs*vulnS[RefYear,i-1])*exp(-NatMs[RefYear])
  }
   tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-inFn*vulnN[RefYear,MaxAge])*exp(-NatMn[RefYear])+ tempNn[j-1,MaxAge]*exp(-inFn*vulnN[RefYear,MaxAge])*exp(-NatMn[RefYear])
   tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-inFs*vulnS[RefYear,MaxAge])*exp(-NatMs[RefYear])+ tempNs[j-1,MaxAge]*exp(-inFs*vulnS[RefYear,MaxAge])*exp(-NatMs[RefYear])

   EggsN			<-sum(tempNn[j-1,]*matureN[RefYear,]*WeightAtAgeN[RefYear,])
   EggsS			<-sum(tempNs[j-1,]*matureS[RefYear,]*WeightAtAgeS[RefYear,])

   moveFromN		<-tempNn[j,]*MovementN[RefYear,]
   moveFromS		<-tempNs[j,]*MovementS[RefYear,]

   tempNn[j,]		<-tempNn[j,]-moveFromN+moveFromS
   tempNs[j,]		<-tempNs[j,]-moveFromS+moveFromN

   tempNn[j,1]		<-inRec
   tempNs[j,1]		<-inRec

   tempRecN[j]		<-tempNn[j,1]
   tempRecS[j]		<-tempNs[j,1]

   tempCatchAtAgeN[j,]	<-((vulnN[RefYear,]*inFn)/(vulnN[RefYear,]*inFn+NatMn[RefYear])) * (1-exp(-(vulnN[RefYear,]*inFn+NatMn[RefYear]))) * tempNn[j-1,]
   tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[RefYear,])

   tempCatchAtAgeS[j,]	<-((vulnS[RefYear,]*inFs)/(vulnS[RefYear,]*inFs+NatMs[RefYear])) * (1-exp(-(vulnS[RefYear,]*inFs+NatMs[RefYear]))) * tempNs[j-1,]
   tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[RefYear,])

 }
if(EggsN>Target*Bzero)
 minF	<-inFn
if(EggsN<Target*Bzero)
 maxF	<-inFn
}

list(inFn,(EggsN)/(inRec))
}