#' Calculates yield at F, to be passed to optimizing function
#'
#' @param x Historical fishing mortality for nothern population
#' @param VirInitN Initial biomass of northern population
#' @param VirInitS Initial biomass of southern population
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
#' @param sigmaRn Recruitment variability for northern population
#' @param sigmaRs Recruitment variability for southern population
#' @param RefYear Year of forecast
#' @param inputFs Historical fishing mortality for southern population
#' @param MovementN Proportion of emigration to the southern population
#' @param MovementS Proportion of immigration into the northern population
#'
#' @return Value of yield at F35%
#'
#' @export
#' @examples
#' \dontrun{
#' 
#' }

FindFMSY<-function(x,MaxAge,VirInitN,VirInitS,vulnN,vulnS,NatMn,NatMs,matureN,matureS,WeightAtAgeN,
 WeightAtAgeS,steepnessN,steepnessS,RzeroN,RzeroS,LenAtAgeN,LenAtAgeS,sigmaRn,sigmaRs,RefYear,
 inputFs=0,MovementN,MovementS)
{
tempNn		<-matrix(ncol=MaxAge,nrow=100)
tempNn[1,]		<-VirInitN
tempNs		<-matrix(ncol=MaxAge,nrow=100)
tempNs[1,]		<-VirInitS
tempCatchN		<-rep(0,100)
tempCatchS		<-rep(0,100)
tempRecN		<-rep(0,100)
tempRecS		<-rep(0,100)
tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=100)
tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=100)

inFs<-inputFs
inFn<-x

for (j in 2:100)
{
 for (i in 2:(MaxAge-1))
  {
   tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-inFn*vulnN[RefYear,i-1])*exp(-NatMn[RefYear])
   tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-inFs*vulnS[RefYear,i-1])*exp(-NatMs[RefYear])
  }
   tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-inFn*vulnN[RefYear,MaxAge])*exp(-NatMn[RefYear])+ tempNn[j-1,MaxAge]*exp(-inFn*vulnN[RefYear,MaxAge])*exp(-NatMn[RefYear])
   tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-inFs*vulnS[RefYear,MaxAge])*exp(-NatMs[RefYear])+ tempNs[j-1,MaxAge]*exp(-inFs*vulnS[RefYear,MaxAge])*exp(-NatMs[RefYear])

   moveFromN		<-tempNn[j,]*MovementN[RefYear,]
   moveFromS		<-tempNs[j,]*MovementS[RefYear,]

   tempNn[j,]		<-tempNn[j,]-moveFromN+moveFromS
   tempNs[j,]		<-tempNs[j,]-moveFromS+moveFromN

   EggsN			<-sum(tempNn[j-1,]*matureN[RefYear,]*WeightAtAgeN[RefYear,])
   EggsS			<-sum(tempNs[j-1,]*matureS[RefYear,]*WeightAtAgeS[RefYear,])

   tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[RefYear],RzeroIN=RzeroN[RefYear],RecErrIN=0,recType="BH",NatMin=NatMn[RefYear],
							vulnIN=vulnN[RefYear,],matureIN=matureN[RefYear,],weightIN=WeightAtAgeN[RefYear,],LenAtAgeIN=LenAtAgeN[RefYear,],MaxAge=MaxAge,sigmaRin=sigmaRn[RefYear])
   tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[RefYear],RzeroIN=RzeroS[RefYear],RecErrIN=0,recType="BH",NatMin=NatMs[RefYear],
							vulnIN=vulnS[RefYear,],matureIN=matureS[RefYear,],weightIN=WeightAtAgeS[RefYear,],LenAtAgeIN=LenAtAgeS[RefYear,],MaxAge=MaxAge,sigmaRin=sigmaRs[RefYear])
   tempRecN[j]		<-tempNn[j,1]
   tempRecS[j]		<-tempNs[j,1]

   tempCatchAtAgeN[j,]	<-((vulnN[RefYear,]*inFn)/(vulnN[RefYear,]*inFn+NatMn[RefYear])) * (1-exp(-(vulnN[RefYear,]*inFn+NatMn[RefYear]))) * tempNn[j-1,]
   tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[RefYear,])

   tempCatchAtAgeS[j,]	<-((vulnS[RefYear,]*inFs)/(vulnS[RefYear,]*inFs+NatMs[RefYear])) * (1-exp(-(vulnS[RefYear,]*inFs+NatMs[RefYear]))) * tempNs[j-1,]
   tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[RefYear,])

 }
obj <-  -1*(tempCatchS[j]+tempCatchN[j])
return(obj)
}