###################################################
# Generalized spatial management strategy evaluation
# Two box, populations processes in each box can
# vary over time, movement between boxes can vary over 
# time, effort in boxes can be allocated asymmetrically
# production and age-structured estimation models can be tested
# written by Cody Szuwalski, 11/2015

#==source a file with helper functions
source("C:/Users/Cody/Desktop/Publications/BlueShark/Helpers.R")

#==this is where the output and files will be stored
CurDir<-"C:/Users/Cody/Desktop/Publications/BlueShark/"
setwd(CurDir)


#=======================
#==simulation controls==
#=======================
Nsim			<-2	# number of simulations to do in the MSE
SimYear		<-60	# total number of years in simulation
InitYear		<-50	# year in which MSE starts (i.e. the number of years of data available)
AssessmentType	<-"Non-spatial age"	# this can be "Production" or "Non-spatial age" ("Spatial age" in the works)
FisheryIndepentDat<-0
LifeHistoryPlots	<-1
 if(LifeHistoryPlots==1)
 {
  dev.new()
  par(mfrow=c(4,4),mar=c(2,3,2,2))
 }

EstimationPlots	<-1

depletion	<-0.3		# model is conditioned to a specific depletion given a trajectory of effort
CatchShareN	<-0.5		# the amount of effort allocated to a given area
CatchShareS	<-1-CatchShareN

#==sampling uncertainty
CatchCVn	<-0.1
CatchCVs	<-0.1

IndexCVn	<-0.2
IndexCVs	<-0.2

LenSampleN	<-200
LenSampleS	<-200

GrowthSDn<-10
GrowthSDs<-10
#==========================================================
#=================population dynamics======================
#=========================================================
#==When making dynamics time-varying, only make them vary after
#=='InitYear'.  The idea is that this simulates a relatively steady
#==environment until climate change, then things change.
#===========================================================
#===========================================================
MaxAge<-30
Ages<-seq(1,MaxAge)

#==natural mortality================
NatMn		<-rep(.2,SimYear)
NatMs		<-rep(.2,SimYear)

#==Length at age====================
#==north==
VonKn		<-rep(.2,SimYear)
LinfN		<-rep(300,SimYear)
t0n		<-rep(.9,SimYear)
LenAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 LenAtAgeN[i,]<-LinfN[i]*(1-exp(-VonKn[i]*(Ages-t0n[i])))

if(LifeHistoryPlots==1)
 {
 plot(LenAtAgeN[1,],type="l",las=1)
 for(x in 2:SimYear)
  lines(LenAtAgeN[x,])
 mtext(side=3,"Length at Age (N)",cex=.7)
 }

#==south==
VonKs		<-rep(.2,SimYear)
LinfS		<-rep(300,SimYear)
t0s		<-rep(.9,SimYear)
LenAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 LenAtAgeS[i,]<-LinfS[i]*(1-exp(-VonKs[i]*(Ages-t0s[i])))

if(LifeHistoryPlots==1)
 {
 plot(LenAtAgeS[1,],type="l",las=1)
 for(x in 2:SimYear)
  lines(LenAtAgeS[x,])
 mtext(side=3,"Length at Age (S)",cex=.7)
 }

#==specify the number of length bins
 BinWidth<-(max(LenAtAgeS,LenAtAgeN)*1.05)/MaxAge
 LengthBins<-seq(0,max(LenAtAgeS,LenAtAgeN)*1.05,BinWidth)
 LengthBinsMid<-LengthBins[1:(length(LengthBins)-1)] + mean(LengthBins[1:2])

#==maturity==========================
#==by age now, can be changed to by length==
#==north==
mat50n	<-rep(6,SimYear)
mat95n	<-rep(8,SimYear)

matureN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 matureN[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50n[i])/(mat95n[i]-mat50n[i])))
 #matureN[i,]	<-1/(1+exp(-1*log(19)*(Ages-mat50n)/(mat95n-mat50n)))

if(LifeHistoryPlots==1)
 {
 plot(matureN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(matureN[x,])
 mtext(side=3,"Maturity at age (N)",cex=.7)
 }

#==south==
mat50s	<-rep(6,SimYear)
mat95s	<-rep(8,SimYear)
matureS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 matureS[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50s[i])/(mat95s[i]-mat50s[i])))
 #matureS[i,]	<-1/(1+exp(-1*log(19)*(Ages-mat50s[i])/(mat95s[i]-mat50s[i])))

if(LifeHistoryPlots==1)
 {
 plot(matureS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(matureS[x,])
 mtext(side=3,"Maturity at age (S)",cex=.7)
 }

#=====================================
#==fishery selectivitiy============== 
#==north==
sel50n	<-rep(190,SimYear)
sel95n	<-rep(225,SimYear)
vulnN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50n[i])/(sel95n[i]-sel50n[i])))

vulnN[vulnN<0.01]<-0

if(LifeHistoryPlots==1)
 {
 plot(vulnN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(vulnN[x,])
 mtext(side=3,"Fishery selectivity at len(N)",cex=.7)

 }

#==south==
sel50s	<-rep(190,SimYear)
sel95s	<-rep(225,SimYear)
vulnS		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50s[i])/(sel95s[i]-sel50s[i])))

vulnS[vulnS<0.01]<-0

if(LifeHistoryPlots==1)
 {
 plot(vulnS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(vulnS[x,])
 mtext(side=3,"Fishery selectivity at len (S)",cex=.7)
 }

#==index selectivity=====================
#==north==
surv50n	<-rep(90,SimYear)
surv95n	<-rep(150,SimYear)
survSelN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50n[i])/(surv95n[i]-surv50n[i])))

if(LifeHistoryPlots==1)
 {
 plot(survSelN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(survSelN[x,])
 mtext(side=3,"Index selectivity at len(N)",cex=.7)
 }

#==south==
surv50s	<-rep(90,SimYear)
surv95s	<-rep(150,SimYear)
survSelS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50s[i])/(surv95s[i]-surv50s[i])))

if(LifeHistoryPlots==1)
 {
 plot(survSelS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(survSelS[x,])
 mtext(side=3,"Index selectivity at len (S)",cex=.7)
 }

#==weight at length==========================
alphaN		<-0.0000017
betaN			<-3
WeightAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtAgeN[i,]	<-alphaN*LenAtAgeN[i,]^betaN

if(LifeHistoryPlots==1)
 {
 plot(WeightAtAgeN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(WeightAtAgeN[x,])
 mtext(side=3,"Weight at length (N)",cex=.7)
 }

alphaS		<-0.0000017
betaS			<-3
WeightAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtAgeS[i,]	<-alphaS*LenAtAgeS[i,]^betaS

if(LifeHistoryPlots==1)
 {
 plot(WeightAtAgeS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(WeightAtAgeS[x,])
 mtext(side=3,"Weight at length (S)",cex=.7)
 }

#==movement from box to box==
#==north
MaxMovingN	<-0
Move50n	<-rep(6,SimYear)
Move95n	<-rep(10,SimYear)
MovementN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 MovementN[i,]	<-MaxMovingN/(1+exp(-1*log(19)*(Ages-Move50n[i])/(Move95n[i]-Move50n[i])))

if(LifeHistoryPlots==1)
 {
 plot(MovementN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(MovementN[x,])
  mtext(side=3,"Movement at age (N)",cex=.7)
 }


#==south
MaxMovingS	<-0
Move50s	<-rep(6,SimYear)
Move95s	<-rep(10,SimYear)
MovementS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 MovementS[i,]	<-MaxMovingS/(1+exp(-1*log(19)*(Ages-Move50s[i])/(Move95s[i]-Move50s[i])))

if(LifeHistoryPlots==1)
 {
 plot(MovementS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(MovementS[x,])
  mtext(side=3,"Movement at age (S)",cex=.7)
 }

#==recruitment (these should have low steepness because of low fecundity)
#==these parameters are hard to specify or estimate with the given data
#==they largely determine the dynamics of the system though...
steepnessN	<-rep(.7,SimYear)
sigmaRn	<-rep(.001,SimYear)
RzeroN	<-rep(100000,SimYear)

steepnessS	<-rep(.7,SimYear)
sigmaRs	<-rep(.001,SimYear)
RzeroS	<-rep(100000,SimYear)

RecErrN	<-matrix(rnorm(SimYear*Nsim,0,sigmaRn[1]),ncol=SimYear)
RecErrS	<-matrix(rnorm(SimYear*Nsim,0,sigmaRs[1]),ncol=SimYear)

#===================================================
#==Virgin numbers at age, biomass, initial depletion, recruitment
#===================================================
VirInitN<-initialN(Rzero=RzeroN[1],NatM=NatMn[1],inAge=MaxAge)
VirInitS<-initialN(Rzero=RzeroS[1],NatM=NatMs[1],inAge=MaxAge)

VirBioN<-sum(VirInitN*matureN[1,]*WeightAtAgeN[1,])
VirBioS<-sum(VirInitS*matureS[1,]*WeightAtAgeS[1,])

ExploitBioN<-sum(VirInitN*vulnN[1,]*WeightAtAgeN[1,])
ExploitBioS<-sum(VirInitS*vulnS[1,]*WeightAtAgeS[1,])

if(LifeHistoryPlots==1)
 {
ssb<-seq(1,VirBioN,VirBioN/100)
record<-ssb
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,])
   }
plot(record~ssb)
mtext(side=3,"SSB vs Rec (N)",cex=.7)
ssb<-seq(1,VirBioS,VirBioS/100)
record<-ssb
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,])
   }
plot(record~ssb)
mtext(side=3,"SSB vs Rec (S)",cex=.7)
}  

#================================
# Set the initial conditions by specifying an F
# that will put the population at the designated
# depletion after the specified time period
# (may have to make this variable--fish up, back off--to get contrast)
#===============================
targetDep	<-(VirBioN+VirBioS)*depletion
maxExp	<-0.999
minExp	<-0.001

#==defines the exploitation pattern
maxMult	<-1.5
levelMult	<-1
peak		<-round(InitYear/3)
levelOff	<-2*peak
ascending	<-seq(0.01,maxMult,(maxMult-0.01)/(peak))
descending	<-seq(maxMult,levelMult,-(maxMult-levelMult)/(peak))
steady	<-rep(levelMult,peak)
ExpSeries	<-c(ascending,descending,steady)
ExpSeries	<-ExpSeries[1:InitYear]

#==finds an exploitation rate that returns a given depletion at the end of the 
#==INCORPORATE THE CATCHSHARE BY AREA INTO THIS CHUNK OF CODE============
for(x in 1:30)
{
tempNn	<-matrix(ncol=MaxAge,nrow=InitYear)
tempNn[1,]	<-VirInitN

tempNs	<-matrix(ncol=MaxAge,nrow=InitYear)
tempNs[1,]	<-VirInitS

tempCatchN	<-rep(0,InitYear)
tempCatchS	<-rep(0,InitYear)

tempRecN	<-tempCatchN
tempRecS	<-tempCatchS

tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=InitYear)
tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=InitYear)

tempExp	<-(maxExp+minExp)/2
tempF		<- -log(1-tempExp)
tempFsrsN	<-ExpSeries*tempF
tempFsrsS	<-ExpSeries*tempF
#tmpExpSrs	<-1-exp(-tempFsrs)

for (j in 2:InitYear)
{
 for (i in 2:(MaxAge-1))
  {
   tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-tempFsrsN[j]*vulnN[1,i-1])*exp(-NatMn[1])
   tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-tempFsrsS[j]*vulnS[1,i-1])*exp(-NatMs[1])

  }
   tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-tempFsrsN[j]*vulnN[1,MaxAge])*exp(-NatMn[1])+ tempNn[j-1,MaxAge]*exp(-tempFsrsN[j]*vulnN[1,MaxAge])*exp(-NatMn[1])
   tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-tempFsrsS[j]*vulnS[1,MaxAge])*exp(-NatMs[1])+ tempNs[j-1,MaxAge]*exp(-tempFsrsS[j]*vulnS[1,MaxAge])*exp(-NatMs[1])

   EggsN			<-sum(tempNn[j-1,]*matureN[1,]*WeightAtAgeN[1,])
   EggsS			<-sum(tempNs[j-1,]*matureS[1,]*WeightAtAgeS[1,])

   tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1,j],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,])
   tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1,j],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,])
   tempRecN[j]		<-tempNn[j,1]
   tempRecS[j]		<-tempNs[j,1]

   tempCatchAtAgeN[j,]	<-((vulnN[1,]*tempFsrsN[j])/(vulnN[1,]*tempFsrsN[j]+NatMn[1])) * (1-exp(-(vulnN[1,]*tempFsrsN[j]+NatMn[1]))) * tempNn[j-1,]
   tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[1,])

   tempCatchAtAgeS[j,]	<-((vulnS[1,]*tempFsrsS[j])/(vulnS[1,]*tempFsrsS[j]+NatMs[1])) * (1-exp(-(vulnS[1,]*tempFsrsS[j]+NatMs[1]))) * tempNs[j-1,]
   tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[1,])

 }

tempDepN	<-sum(tempNn[InitYear,]*matureN[1,]*WeightAtAgeN[1,])
tempDepS	<-sum(tempNs[InitYear,]*matureS[1,]*WeightAtAgeS[1,])
tempDep	<-tempDepN+tempDepS

if(tempDep>targetDep)
 minExp	<-tempExp
if(tempDep<targetDep)
 maxExp	<-tempExp
initExp	<-tempExp
}

HistoricFmort<-ExpSeries*tempExp
mtext(side=3,"Exploitation history",cex=.7)
round((-log(1-HistoricFmort)),3)
#===============================================================
# BEGIN SIMULATION OF ASSESSMENT AND HARVEST
#===============================================================

#==tempNn and tempNs from above are the starting points
projNn	<-array(dim=c(SimYear,MaxAge,Nsim))
projNs	<-array(dim=c(SimYear,MaxAge,Nsim))

projCatchAtAgeN	<-array(dim=c(SimYear,MaxAge,Nsim))
projCatchAtAgeS	<-array(dim=c(SimYear,MaxAge,Nsim))

projCatchN	<-matrix(nrow=Nsim,ncol=SimYear)
projCatchS	<-matrix(nrow=Nsim,ncol=SimYear)
projSSBn	<-matrix(nrow=Nsim,ncol=SimYear)
projSSBs	<-matrix(nrow=Nsim,ncol=SimYear)
projExpBn	<-matrix(nrow=Nsim,ncol=SimYear)
projExpBs	<-matrix(nrow=Nsim,ncol=SimYear)
projSurvN	<-matrix(nrow=Nsim,ncol=SimYear)
projSurvS	<-matrix(nrow=Nsim,ncol=SimYear)
projRecN	<-matrix(nrow=Nsim,ncol=SimYear)
projRecS	<-matrix(nrow=Nsim,ncol=SimYear)
projFmortN	<-matrix(nrow=Nsim,ncol=SimYear)
projFmortS	<-matrix(nrow=Nsim,ncol=SimYear)

projCatLenFreqN	<-array(dim=c(SimYear,MaxAge,Nsim))
projCatLenFreqS	<-array(dim=c(SimYear,MaxAge,Nsim))
projSurvLenFreqN	<-array(dim=c(SimYear,MaxAge,Nsim))
projSurvLenFreqS	<-array(dim=c(SimYear,MaxAge,Nsim))

#==============================================================
# calculate the catch proportion at (length) for assessment
#==============================================================
library(splines)

 tempCatAtLenN<-matrix(nrow=MaxAge,ncol=MaxAge)
 tempCatAtLenS<-matrix(nrow=MaxAge,ncol=MaxAge)
 GrowthSDn<-10
 GrowthSDs<-10

for(x in 1:Nsim)
 for(y in 2:InitYear)
 {

 #==make length frequencies for catch==
 for(w in 1:nrow(tempCatAtLenN))
 {
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
  ProbN<-probtemp/sum(probtemp)
  tempCatAtLenN[w,]<-tempCatchAtAgeN[y,w]*ProbN
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
  ProbS<-probtemp/sum(probtemp)
  tempCatAtLenS[w,]<-tempCatchAtAgeS[y,w]*ProbS
 }

 projCatLenFreqN[y,,x]<-apply(tempCatAtLenN,2,sum)
 projCatLenFreqS[y,,x]<-apply(tempCatAtLenS,2,sum)
}

if(LifeHistoryPlots==1)
 {
 dev.new()
 par(mfrow=c(ceiling(sqrt(InitYear)),ceiling(sqrt(InitYear))),mar=c(.1,.1,.1,.1))
 for(i in 2:InitYear)
  barplot(rbind(projCatLenFreqN[i,,1],projCatLenFreqS[i,,1]),beside=T,xaxt='n',yaxt='n')
 }

#==============================================================
# calculate the survey proportion at age (length) for assessment
#==============================================================
 tempSurvAtLenN<-matrix(nrow=MaxAge,ncol=MaxAge)
 tempSurvAtLenS<-matrix(nrow=MaxAge,ncol=MaxAge)
for(x in 1:Nsim)
 for(y in 2:InitYear)
 {
 #==make length frequencies for s==
 for(w in 1:nrow(tempCatAtLenN))
 {
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
  ProbN<-probtemp/sum(probtemp)
  tempSurvAtLenN[w,]<-tempNn[y,w]*survSelN[y,w]*ProbN
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
  ProbS<-probtemp/sum(probtemp)
  tempSurvAtLenS[w,]<-tempNs[y,w]*survSelS[y,w]*ProbS
 }

 projSurvLenFreqN[y,,x]<-apply(tempSurvAtLenN,2,sum)
 projSurvLenFreqS[y,,x]<-apply(tempSurvAtLenS,2,sum)
}

if(LifeHistoryPlots==1)
 {
 dev.new()
 par(mfrow=c(ceiling(sqrt(InitYear)),ceiling(sqrt(InitYear))),mar=c(.1,.1,.1,.1))
 for(i in 2:InitYear)
  barplot(rbind(projSurvLenFreqN[i,,1],projSurvLenFreqS[i,,1]),beside=T,xaxt='n',yaxt='n')
 }

#==transfer historic time series from above
for(x in 1:Nsim)
{
projNn[1:InitYear,,x]			<-tempNn
projNs[1:InitYear,,x]			<-tempNs
projCatchN[x,1:InitYear]		<-tempCatchN
projCatchS[x,1:InitYear]		<-tempCatchS
projRecN[x,1:InitYear]			<-tempRecN
projRecS[x,1:InitYear]			<-tempRecS
projFmortN[x,1:InitYear]		<-HistoricFmort
projFmortS[x,1:InitYear]		<-HistoricFmort
projSurvN[x,1:InitYear]			<-apply(tempNn*survSelN[1:InitYear,]*WeightAtAgeN[1:InitYear,],1,sum)
projSurvS[x,1:InitYear]			<-apply(tempNs*survSelS[1:InitYear,]*WeightAtAgeS[1:InitYear,],1,sum)
}

for(x in 1:InitYear)
{
 projSSBn[,x]	<-sum(projNn[x,,1]*matureN[1,]*WeightAtAgeN[1,])
 projSSBs[,x]	<-sum(projNs[x,,1]*matureS[1,]*WeightAtAgeS[1,])
 projExpBn[,x]	<-sum(projNn[x,,1]*vulnN[1,]*WeightAtAgeN[1,])
 projExpBs[,x]	<-sum(projNs[x,,1]*vulnS[1,]*WeightAtAgeS[1,])
}

CatchErrorN		<-matrix(rnorm(Nsim*SimYear,0,CatchCVn),nrow=Nsim,ncol=SimYear)
CatchErrorS		<-matrix(rnorm(Nsim*SimYear,0,CatchCVs),nrow=Nsim,ncol=SimYear)
CatchAssessN	<-projCatchN*exp(CatchErrorN)
CatchAssessS	<-projCatchS*exp(CatchErrorS)

CPUEErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
CPUEErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
CPUEAssessN		<-projExpBn*exp(CPUEErrorN)
CPUEAssessS		<-projExpBs*exp(CPUEErrorS)

SurvErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
SurvErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
SurvAssessN		<-projSurvN*exp(SurvErrorN)
SurvAssessS		<-projSurvS*exp(SurvErrorS)

#==save management and assessment quantities
FMSY<-matrix(nrow=Nsim,ncol=SimYear)
BMSY<-matrix(nrow=Nsim,ncol=SimYear)
TAC<-matrix(nrow=Nsim,ncol=SimYear)
CurBio<-matrix(nrow=Nsim,ncol=SimYear)

#=================================================================================
#==BEGIN SIMULATION===========================================================

for(z in 1:Nsim)
{
 for(y in InitYear:(SimYear-1))
 {
 #==pull data for assessment
 CatchData<-CatchAssessN[z,!is.na(CatchAssessN[z,])] + CatchAssessS[z,!is.na(CatchAssessS[z,])]
 CPUEData<-CPUEAssessN[z,!is.na(CPUEAssessN[z,])] + CPUEAssessS[z,!is.na(CPUEAssessS[z,])] 
 SurvData<-SurvAssessN[z,!is.na(SurvAssessN[z,])] + SurvAssessS[z,!is.na(SurvAssessS[z,])] 

 #=================================================================================
 #===ASSESS THE STOCK, CHOOSE A METHOD AT THE TOP=================================
 #=================================================================================

 #==Production model==
 if(AssessmentType == "Production")
 {
 x		<-c(1.5*(VirBioN+VirBioS),0.2)	#initial values
 outs		<-nlminb(start=x,objective=ProdMod,CatchData=CatchData,IndexData=CPUEData )
 PredBio	<-ProdModPlot(outs$par,CatchData,CPUEData,plots=EstimationPlots)
 FMSY[z,y] 	<-abs(outs$par[2])/2
 BMSY[z,y]	<-abs(outs$par[1])/2
 MSY		<-BMSY[z,y]*FMSY[z,y]
 CurBio[z,y]<-PredBio[length(PredBio)]
 TAC[z,y]	<-CurBio[z,y]*FMSY[z,y]
 }

 #==Age-structured model
 if(AssessmentType == "Non-spatial age")
 {
 #==make folder if it isn't there
 setwd(CurDir)
 dir.create("AgeAssOutputs",showWarnings=F)
 IndSimFolder<-paste("AgeAssOutputs/",y,sep="")
 dir.create(IndSimFolder,showWarnings=F)

 #==copy the .exe into it (where does this come from? github?)
 file.copy(from="GenAss/simass.exe",to=IndSimFolder)
 
 #==write the true values
 setwd(IndSimFolder)
 file.create("TrueQuantities.DAT")
 cat("#True quantities","\n",file="TrueQuantities.DAT")
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#recruitment","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(0,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#Catch","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(projCatchN[z,]+projCatchS[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#survey selectivity","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(c(surv50n,surv95n),"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#fishery selectivity","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(c(sel50n,sel95n),"\n",file="TrueQuantities.DAT",append=TRUE)
# cat("#fishing mortality","\n",file="TrueQuantities.DAT",append=TRUE)
# cat(TrueFmort,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#spawning biomass","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(projSSBn[z,]+projSSBs[z,],"\n",file="TrueQuantities.DAT",append=TRUE)

 selAtAgeFunc<-function(inputLen,K,Linf,t0)
 {
  selAtAge<-((log(1-inputLen/Linf)/-K))+t0
 return(selAtAge)
 }
 selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
 selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])
 selAtAgeFunc(surv50n[1],VonKn[1],LinfN[1],t0n[1])
 selAtAgeFunc(surv95n[1],VonKn[1],LinfN[1],t0n[1])

 #==.CTL file  !!!!!!!!!!!!!!!!!how should we deal with space?!!!!!!!!!!!!!!!!!!!!!!
 file.create("SimAss.CTL")
 cat("#Simulated assessment control file","\n",file="SimAss.CTL")
 cat("#","\n",file="SimAss.CTL",append=TRUE)
 cat("#weight parameters","\n",file="SimAss.CTL",append=TRUE)
 cat(c(alphaN,betaN) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#length parameters","\n",file="SimAss.CTL",append=TRUE)
 cat(c(VonKn[1],LinfN[1],t0n[1]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#natural mortality","\n",file="SimAss.CTL",append=TRUE)
 cat(NatMn[length(NatMn)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#maturity","\n",file="SimAss.CTL",append=TRUE)
 cat(c(mat50n[length(mat50n)],mat95n[length(mat95n)]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#steepness","\n",file="SimAss.CTL",append=TRUE)
 cat(steepnessN[length(steepnessN)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#small number","\n",file="SimAss.CTL",append=TRUE)
 cat(0.0001,"\n",file="SimAss.CTL",append=TRUE)
 cat("#initial smoothness weight","\n",file="SimAss.CTL",append=TRUE)
 cat(0.1,"\n",file="SimAss.CTL",append=TRUE) 
 cat("#Fishery Independent data available?","\n",file="SimAss.CTL",append=TRUE)
 cat(FisheryIndepentDat,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Length at age standard deviation","\n",file="SimAss.CTL",append=TRUE)
 cat(GrowthSDn,"\n",file="SimAss.CTL",append=TRUE)

 #==.DAT file !!!!!!!!!!!!!!!!!FIX THIS!!!!!!!!!!!!!!!!!!!!!!
 file.create("SimAss.DAT")
 cat("#Simulated generalized assessment","\n",file="SimAss.DAT")
 cat("#","\n",file="SimAss.DAT",append=TRUE)
 cat("#start/end year","\n",file="SimAss.DAT",append=TRUE)
 cat(c(1,y),"\n",file="SimAss.DAT",append=TRUE)

 #==aggregate the length frequencies from the two areas
 cat("#number of ages","\n",file="SimAss.DAT",append=TRUE)
 cat(MaxAge,"\n",file="SimAss.DAT",append=TRUE)
 cat("#ages","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,MaxAge),"\n",file="SimAss.DAT",append=TRUE)

 #==aggregate the index of abundance from two areas
 cat("#CPUE years","\n",file="SimAss.DAT",append=TRUE)
 cat(y,"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for CPUE","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,y),"\n",file="SimAss.DAT",append=TRUE)

 cat("#survey years","\n",file="SimAss.DAT",append=TRUE)
 cat(y,"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for survey","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,y),"\n",file="SimAss.DAT",append=TRUE)

 cat("#length freq from catch years","\n",file="SimAss.DAT",append=TRUE)
 cat(y,"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for length freq catch","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,y),"\n",file="SimAss.DAT",append=TRUE)
 cat("#catch length sample sizes","\n",file="SimAss.DAT",append=TRUE)
 cat(LenSampleN/10,"\n",file="SimAss.DAT",append=TRUE) 

 cat("#length freq from survey years","\n",file="SimAss.DAT",append=TRUE)
 cat(y,"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for length freq survey","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,y),"\n",file="SimAss.DAT",append=TRUE)
 cat("#survey sample sizes","\n",file="SimAss.DAT",append=TRUE)
 cat(LenSampleN/10,"\n",file="SimAss.DAT",append=TRUE) 

 cat("#total survey biomass","\n",file="SimAss.DAT",append=TRUE) 
 cat(SurvData,"\n",file="SimAss.DAT",append=TRUE)
 cat("#survey CV","\n",file="SimAss.DAT",append=TRUE)
 cat(IndexCVn,"\n",file="SimAss.DAT",append=TRUE)

 cat("#total CPUE","\n",file="SimAss.DAT",append=TRUE) 
 cat(CPUEData,"\n",file="SimAss.DAT",append=TRUE)
 cat("#CPUE CV","\n",file="SimAss.DAT",append=TRUE)
 cat(IndexCVn,"\n",file="SimAss.DAT",append=TRUE)

 cat("#total catch biomass","\n",file="SimAss.DAT",append=TRUE)
 cat(CatchData,"\n",file="SimAss.DAT",append=TRUE)
 cat("#catch CV","\n",file="SimAss.DAT",append=TRUE)
 cat(CatchCVn,"\n",file="SimAss.DAT",append=TRUE)
 cat("#Length Bins","\n",file="SimAss.DAT",append=TRUE)
 cat(LengthBinsMid,"\n",file="SimAss.DAT",append=TRUE)

projCatLenFreqN[is.na(projCatLenFreqN)]<-0
projCatLenFreqS[is.na(projCatLenFreqS)]<-0
projSurvLenFreqN[is.na(projSurvLenFreqN)]<-0
projSurvLenFreqS[is.na(projSurvLenFreqS)]<-0

 cat("#catch length counts","\n",file="SimAss.DAT",append=TRUE)
 for(i in 1:y)
  cat(projCatLenFreqN[i,,z]+projCatLenFreqS[i,,z],"\n",file="SimAss.DAT",append=TRUE)

 cat("#survey length counts","\n",file="SimAss.DAT",append=TRUE)
 for(i in 1:y)
  cat(projSurvLenFreqN[i,,z]+projSurvLenFreqS[i,,z],"\n",file="SimAss.DAT",append=TRUE)

 #==run the code
 shell("simass -nohess")

#==================================
# PLOT THE OUTPUT (if instructed)
#==================================
REP<-readLines(paste(CurDir,IndSimFolder,"/simass.rep",sep=""))
CTL<-readLines(paste(CurDir,IndSimFolder,"/simass.CTL",sep=""))
DAT<-readLines(paste(CurDir,IndSimFolder,"/simass.DAT",sep=""))
TRU<-readLines(paste(CurDir,IndSimFolder,"/TrueQuantities.DAT",sep=""))
data2<-scan(paste(CurDir,IndSimFolder,"/simass.PAR",sep=""),what="character")

 if(EstimationPlots==1)
  AgeAssPlot(REP,CTL,DAT,TRU,data2)

temp<-grep("OFL",REP)[2]
OFL<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

 #==calculate the TAC based on a SRR or a proxy, depending on data
 TAC[z,y]	<-OFL
 }

 if(EstimationPlots==1 & AssessmentType == "Production")
 {
  lines(BMSY[z,],col=3,lty=3)
  legend("topright",bty='n',col=c(NA,1,2,3),pch=c(15,NA,NA,NA),lty=c(NA,1,2,3),
          legend=c("Index obs","Index est","Catch","BMSY"))
 }
 #==================================================================================
 #==find the F that removes the TAC given the fishery selectivity===================
 #==this TAC is not necessarily the 'true' TAC--there is error in the biomass and FMSY
 #==================================================================================

 minExp<-0.00001
 maxExp<-0.99999

 for(p in 1:30)
 {
 tempExp	<-(maxExp+minExp)/2
 tempF	<- -log(1-tempExp)
 tempFn	<-tempF*CatchShareN
 tempFs	<-tempF*CatchShareS

   finCatLenFn	<-((tempFn*vulnN[y,])/(vulnN[y,]*tempFn+NatMn[y])) * (1-exp(-((tempFn*vulnN[y,])+NatMn[y]))) * projNn[y,,z]
   finCatFn		<-sum(finCatLenFn*WeightAtAgeN[y,])

   finCatLenFs	<-((tempFs*vulnS[y,])/(vulnS[y,]*tempFs+NatMs[y])) * (1-exp(-((tempFs*vulnS[y,])+NatMs[y]))) * projNn[y,,z]
   finCatFs		<-sum(finCatLenFs*WeightAtAgeS[y,])

   TempTotCatch	<-finCatFs+finCatFn

 if(TempTotCatch<TAC[z,y])
  minExp	<-tempExp
 if(TempTotCatch>TAC[z,y])
  maxExp	<-tempExp
 initExp	<-tempExp
 }

 ApplyF	<- -log(1-tempExp)
 ApplyFn	<-ApplyF*CatchShareN
 ApplyFs	<-ApplyF*CatchShareS
 
 #==apply F that would take a calculated TAC====================================================
 #==update population dynamics
  for (i in 2:(MaxAge-1))
  {
   projNn[y+1,i,z]		<-projNn[y,i-1,z]*exp(-ApplyFn*vulnN[y,i-1])*exp(-NatMn[y])
   projNs[y+1,i,z]		<-projNs[y,i-1,z]*exp(-ApplyFs*vulnS[y,i-1])*exp(-NatMs[y])
  }
   projNn[y+1,MaxAge,z]		<-(projNn[y,(MaxAge-1),z])*exp(-ApplyFn*vulnN[y,MaxAge])*exp(-NatMn[y])+ projNn[y,MaxAge,z]*exp(-ApplyFn*vulnN[y,MaxAge])*exp(-NatMn[y])
   projNs[y+1,MaxAge,z]		<-(projNs[y,(MaxAge-1),z])*exp(-ApplyFs*vulnS[y,MaxAge])*exp(-NatMs[y])+ projNs[y,MaxAge,z]*exp(-ApplyFn*vulnS[y,MaxAge])*exp(-NatMs[y])

   EggsN				<-sum(projNn[y-1,,z]*matureN[y,]*WeightAtAgeN[y,])
   EggsS				<-sum(projNs[y-1,,z]*matureS[y,]*WeightAtAgeS[y,])

   projNn[y+1,1,z]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[y],RzeroIN=RzeroN[y],RecErrIN=RecErrN[z,y],recType="BH",NatMin=NatMn[y],
							vulnIN=vulnN[y,],matureIN=matureN[y,],weightIN=WeightAtAgeN[y,],LenAtAgeIN=LenAtAgeN[y,])
   projNs[y+1,1,z]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[y],RzeroIN=RzeroS[y],RecErrIN=RecErrS[z,y],recType="BH",NatMin=NatMs[y],
							vulnIN=vulnS[y,],matureIN=matureS[y,],weightIN=WeightAtAgeS[y,],LenAtAgeIN=LenAtAgeS[y,])

   projCatchAtAgeN[y+1,,z]	<-((ApplyFn*vulnN[y,])/(vulnN[y,]*ApplyFn+NatMn[y])) * (1-exp(-((ApplyFn*vulnN[y,])+NatMn[y]))) * projNn[y+1,,z]
   projCatchN[z,y+1]		<-sum(projCatchAtAgeN[y+1,,z]*WeightAtAgeN[y,])

   projCatchAtAgeS[y+1,,z]	<-((ApplyFs*vulnS[y,])/(vulnS[y,]*ApplyFs+NatMs[y])) * (1-exp(-((ApplyFs*vulnS[y,])+NatMs[y]))) * projNs[y+1,,z]
   projCatchS[z,y+1]		<-sum(projCatchAtAgeS[y+1,,z]*WeightAtAgeS[y,])

   #==UPDATE CATCH DATA AND INDEX (survey and CPUE) DATA
   projSSBn[z,y+1]	<-sum(projNn[y+1,,z]*matureN[y,]*WeightAtAgeN[y,])
   projSSBs[z,y+1]	<-sum(projNs[y+1,,z]*matureS[y,]*WeightAtAgeS[y,])
   projExpBn[z,y+1]	<-sum(projNn[y+1,,z]*vulnN[y,]*WeightAtAgeN[y,])
   projExpBs[z,y+1]	<-sum(projNs[y+1,,z]*vulnS[y,]*WeightAtAgeS[y,])
   projSurvN[z,y+1]	<-sum(projNn[y+1,,z]*survSelN[y,]*WeightAtAgeN[y,])
   projSurvS[z,y+1]	<-sum(projNs[y+1,,z]*survSelS[y,]*WeightAtAgeS[y,])

   CatchAssessN[z,y+1]	<-projCatchN[z,y+1]*exp(CatchErrorN[z,y+1])
   CatchAssessS[z,y+1]	<-projCatchS[z,y+1]*exp(CatchErrorS[z,y+1])

   CPUEAssessN[z,y+1]	<-projExpBn[z,y+1]*exp(CPUEErrorN[z,y+1])
   CPUEAssessS[z,y+1]	<-projExpBs[z,y+1]*exp(CPUEErrorS[z,y+1])

   SurvAssessN[z,y+1]	<-projSurvN[z,y+1]*exp(CPUEErrorN[z,y+1])
   SurvAssessS[z,y+1]	<-projSurvS[z,y+1]*exp(CPUEErrorS[z,y+1])
   
  #==UPDATE CATCH LENGTH FREQS====
    for(w in 1:nrow(tempCatAtLenN))
	 {
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
	  ProbN<-probtemp/sum(probtemp)
	  tempCatAtLenN[w,]<-projCatchAtAgeN[y+1,w,z]*ProbN
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
	  ProbS<-probtemp/sum(probtemp)
	  tempCatAtLenS[w,]<-projCatchAtAgeS[y+1,w,z]*ProbS
	 }
	 projCatLenFreqN[y+1,,z]<-apply(tempCatAtLenN,2,sum)
	 projCatLenFreqS[y+1,,z]<-apply(tempCatAtLenS,2,sum)
	

  #==============================================================
  # calculate the survey proportion at age (length) for assessment
  #==============================================================
	 #==make length frequencies for s==
   for(w in 1:nrow(tempCatAtLenN))
    {
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
	ProbN<-probtemp/sum(probtemp)
	tempSurvAtLenN[w,]<-projNn[y,w,z]*survSelN[y,w]*ProbN
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
	ProbS<-probtemp/sum(probtemp)
	tempSurvAtLenS[w,]<-projNs[y,w,z]*survSelS[y,w]*ProbS
     }

    projSurvLenFreqN[y+1,,z]<-apply(tempSurvAtLenN,2,sum)
    projSurvLenFreqS[y+1,,z]<-apply(tempSurvAtLenS,2,sum)

 } # end y
} # end z


