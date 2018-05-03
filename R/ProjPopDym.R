#' Project population dynamics
#' 
#' @return List of projected yield and biomass
#'
#' @export
#' 
ProjPopDym<-function(fmortN,fmortS=0,MaxAge,vulnN,vulnS,NatMn,NatMs,matureN,matureS,WeightAtAgeN,
 WeightAtAgeS,steepnessN,steepnessS,RzeroN,RzeroS,LenAtAgeN,LenAtAgeS,ProxyRec=0,sigmaRn,sigmaRs,
 RefYear,MovementN,MovementS)
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

inFs<-fmortS
inFn<-fmortN

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
if(ProxyRec==0)
{
   tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[RefYear],RzeroIN=RzeroN[RefYear],RecErrIN=0,recType="BH",NatMin=NatMn[RefYear],
							vulnIN=vulnN[RefYear,],matureIN=matureN[RefYear,],weightIN=WeightAtAgeN[RefYear,],LenAtAgeIN=LenAtAgeN[RefYear,],MaxAge,sigmaRin=sigmaRn[RefYear])
   tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[RefYear],RzeroIN=RzeroS[RefYear],RecErrIN=0,recType="BH",NatMin=NatMs[RefYear],
							vulnIN=vulnS[RefYear,],matureIN=matureS[RefYear,],weightIN=WeightAtAgeS[RefYear,],LenAtAgeIN=LenAtAgeS[RefYear,],MaxAge,sigmaRin=sigmaRs[RefYear])
}
if(ProxyRec>0)
{
   tempNn[j,1]		<-ProxyRec
   tempNs[j,1]		<-ProxyRec
}
   tempRecN[j]		<-tempNn[j,1]
   tempRecS[j]		<-tempNs[j,1]

   tempCatchAtAgeN[j,]	<-((vulnN[RefYear,]*inFn)/(vulnN[RefYear,]*inFn+NatMn[RefYear])) * (1-exp(-(vulnN[RefYear,]*inFn+NatMn[RefYear]))) * tempNn[j-1,]
   tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[RefYear,])

   tempCatchAtAgeS[j,]	<-((vulnS[RefYear,]*inFs)/(vulnS[RefYear,]*inFs+NatMs[RefYear])) * (1-exp(-(vulnS[RefYear,]*inFs+NatMs[RefYear]))) * tempNs[j-1,]
   tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[RefYear,])

 }
 outYield	<-tempCatchN[j]
 outBiomass	<-EggsN
list(outYield,outBiomass)
}
