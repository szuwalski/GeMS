###################################################
# Generalized spatial management strategy evaluation
# Two box, populations processes in each box can
# vary over time, movement between boxes can vary over 
# time, effort in boxes can be allocated asymmetrically
# production and age-structured estimation models can be tested
# written by Cody Szuwalski, 11/2015

#==source a file with helper functions
source("C:/General-MSE/General MSE helpers.R")

#==this is where the output and files will be stored
CurDir<-"C:/General-MSE/"
setwd(CurDir)
input<-"C:/General-MSE/CTLfile_smallyellowcroaker.csv"
input<-"C:/General-MSE/CTLfile.csv"
out<-ReadCTLfile(input)

CalcFMSY<-0
#=======================
#==simulation controls==
#=======================
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
AssessmentType	<-out$OM$AssessmentType		# this can be "Production" or "Non-spatial age" ("Spatial age" and "Time-varying production" in the works)
FisheryIndepenDat <-out$OM$FisheryIndepenDat
LifeHistoryPlots	<-out$OM$LifeHistoryPlots
EstimationPlots	<-out$OM$EstimationPlots	# this plots the diagnostics for each run of the estimation model (will be a lot of plots!)
AssessmentData	<-out$OM$AssessmentData		# this plots the diagnostics for each run of the estimation model (will be a lot of plots!)

depletion		<-out$OM$depletion		# model is conditioned to a specific depletion given a trajectory of effort
CatchShareN		<-out$OM$CatchShareN		# the amount of effort allocated to a given area
CatchShareS		<-1-CatchShareN

#==sampling uncertainty
CatchCVn	<-out$OM$CatchCVn
CatchCVs	<-out$OM$CatchCVs

IndexCVn	<-out$OM$IndexCVn
IndexCVs	<-out$OM$IndexCVs

LenSampleN	<-out$OM$LenSampleN
LenSampleS	<-out$OM$LenSampleS

GrowthSDn	<-out$OM$GrowthSDn
GrowthSDs	<-out$OM$GrowthSDs

#==========================================================
#=================population dynamics======================
#=========================================================
#==When making dynamics time-varying, only make them vary after
#=='InitYear'.  The idea is that this simulates a relatively steady
#==environment until climate change, then things change. MEH.
#===========================================================
#===========================================================
MaxAge<-out$OM$MaxAge
Ages<-seq(1,MaxAge)

#==natural mortality================
 NatMn		<-CleanInput(out$OM$NatMn,SimYear)
 NatMs		<-CleanInput(out$OM$NatMn,SimYear)

#==Length at age====================
VonKn		<-CleanInput(out$OM$VonKn,SimYear)
LinfN		<-CleanInput(out$OM$LinfN,SimYear)
t0n		<-CleanInput(out$OM$t0n,SimYear)
LenAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 LenAtAgeN[i,]<-LinfN[i]*(1-exp(-VonKn[i]*(Ages-t0n[i])))

VonKs		<-CleanInput(out$OM$VonKs,SimYear)
LinfS		<-CleanInput(out$OM$LinfS,SimYear)
t0s		<-CleanInput(out$OM$t0s,SimYear)
LenAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 LenAtAgeS[i,]<-LinfS[i]*(1-exp(-VonKs[i]*(Ages-t0s[i])))

#==specify the number of length bins
 BinWidth<-(max(LenAtAgeS,LenAtAgeN)*1.05)/MaxAge
 LengthBins<-seq(0,max(LenAtAgeS,LenAtAgeN)*1.05,BinWidth)
 LengthBinsMid<-LengthBins[1:(length(LengthBins)-1)] + mean(LengthBins[1:2])

#==maturity at age==========================
mat50n	<-CleanInput(out$OM$mat50n,SimYear)
mat95n	<-CleanInput(out$OM$mat95n,SimYear)
matureN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 matureN[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50n[i])/(mat95n[i]-mat50n[i])))

mat50s	<-CleanInput(out$OM$mat50s,SimYear)
mat95s	<-CleanInput(out$OM$mat95s,SimYear)
matureS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 matureS[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50s[i])/(mat95s[i]-mat50s[i])))

#=====================================
#==fishery selectivitiy============== 
sel50n	<-CleanInput(out$OM$sel50n,SimYear)
sel95n	<-CleanInput(out$OM$sel95n,SimYear)
vulnN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50n[i])/(sel95n[i]-sel50n[i])))

vulnN[vulnN<0.01]<-0

sel50s	<-CleanInput(out$OM$sel50s,SimYear)
sel95s	<-CleanInput(out$OM$sel95s,SimYear)
vulnS		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50s[i])/(sel95s[i]-sel50s[i])))

vulnS[vulnS<0.01]<-0

#==index selectivity=====================
surv50n	<-CleanInput(out$OM$surv50n,SimYear)
surv95n	<-CleanInput(out$OM$surv95n,SimYear)
survSelN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50n[i])/(surv95n[i]-surv50n[i])))

surv50s	<-CleanInput(out$OM$surv50s,SimYear)
surv95s	<-CleanInput(out$OM$surv95s,SimYear)
survSelS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50s[i])/(surv95s[i]-surv50s[i])))

#==weight at age==========================
alphaN		<-CleanInput(out$OM$alphaN,SimYear)
betaN			<-CleanInput(out$OM$betaN,SimYear)
WeightAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtAgeN[i,]	<-alphaN[i]*LenAtAgeN[i,]^betaN[i]

alphaS		<-CleanInput(out$OM$alphaS,SimYear)
betaS			<-CleanInput(out$OM$betaS,SimYear)
WeightAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtAgeS[i,]	<-alphaS[i]*LenAtAgeS[i,]^betaS[i]

#==movement from box to box==
MaxMovingN	<-CleanInput(out$OM$MaxMovingN,SimYear)
Move50n	<-CleanInput(out$OM$Move50n,SimYear)
Move95n	<-CleanInput(out$OM$Move95n,SimYear)
MovementN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 MovementN[i,]	<-MaxMovingN[i]/(1+exp(-1*log(19)*(Ages-Move50n[i])/(Move95n[i]-Move50n[i])))

MaxMovingS	<-CleanInput(out$OM$MaxMovingS,SimYear)
Move50s	<-CleanInput(out$OM$Move50s,SimYear)
Move95s	<-CleanInput(out$OM$Move95s,SimYear)
MovementS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 MovementS[i,]	<-MaxMovingS[i]/(1+exp(-1*log(19)*(Ages-Move50s[i])/(Move95s[i]-Move50s[i])))

#==recruitment parameters==
steepnessN	<-CleanInput(out$OM$steepnessN,SimYear)
sigmaRn	<-CleanInput(out$OM$sigmaRn,SimYear)
RzeroN	<-CleanInput(out$OM$RzeroN,SimYear)

steepnessS	<-CleanInput(out$OM$steepnessS,SimYear)
sigmaRs	<-CleanInput(out$OM$sigmaRs,SimYear)
RzeroS	<-CleanInput(out$OM$RzeroS,SimYear)

RecErrN	<-matrix(rnorm(SimYear*Nsim,0,sigmaRn[1]),ncol=SimYear)
RecErrS	<-matrix(rnorm(SimYear*Nsim,0,sigmaRs[1]),ncol=SimYear)

#==historic fishing mortality
HistoricalF		<-CleanInput(out$OM$HistoricalF,SimYear)[1:InitYear]
HarvestControl	<-out$OM$HarvestControl
ConstantCatch	<-out$OM$ConstantCatch
ConstantF		<-out$OM$ConstantF
HCalpha		<-out$OM$HCalpha
HCbeta		<-out$OM$HCbeta

SmallNum		<-out$OM$SmalNum
InitSmooth		<-out$OM$InitSmooth	
FmortPen		<-out$OM$FmortPen
RecruitPen		<-out$OM$RecruitPen		

#===================================================
#==Virgin numbers at age, biomass, initial depletion, recruitment
#===================================================
VirInitN<-initialN(Rzero=RzeroN[1],NatM=NatMn[1],inAge=MaxAge)
VirInitS<-initialN(Rzero=RzeroS[1],NatM=NatMs[1],inAge=MaxAge)

VirBioN<-sum(VirInitN*matureN[1,]*WeightAtAgeN[1,])
VirBioS<-sum(VirInitS*matureS[1,]*WeightAtAgeS[1,])

ExploitBioN<-sum(VirInitN*vulnN[1,]*WeightAtAgeN[1,])
ExploitBioS<-sum(VirInitS*vulnS[1,]*WeightAtAgeS[1,])

#========================================================================
# FIND MSY FOR THE POPULATION (should just make a popdym function, dummy)=
#========================================================================
if(CalcFMSY==1)
{
SearchFmort		<-seq(0.01,3*NatMn[1],(NatMn[1]-0.01)/100)
SearchYield		<-rep(0,length(SearchFmort))
SearchBiomass	<-rep(0,length(SearchFmort))
for(p in 1:length(SearchFmort))
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

inFs<-SearchFmort[p]*CatchShareS
inFn<-SearchFmort[p]*CatchShareN

for (j in 2:100)
{
 for (i in 2:(MaxAge-1))
  {
   tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-inFn*vulnN[1,i-1])*exp(-NatMn[1])
   tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-inFs*vulnS[1,i-1])*exp(-NatMs[1])

  }
   tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-inFn*vulnN[1,MaxAge])*exp(-NatMn[1])+ tempNn[j-1,MaxAge]*exp(-inFn*vulnN[1,MaxAge])*exp(-NatMn[1])
   tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-inFs*vulnS[1,MaxAge])*exp(-NatMs[1])+ tempNs[j-1,MaxAge]*exp(-inFs*vulnS[1,MaxAge])*exp(-NatMs[1])

   EggsN			<-sum(tempNn[j-1,]*matureN[1,]*WeightAtAgeN[1,])
   EggsS			<-sum(tempNs[j-1,]*matureS[1,]*WeightAtAgeS[1,])

   tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1,j],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,])
   tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1,j],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,])
   tempRecN[j]		<-tempNn[j,1]
   tempRecS[j]		<-tempNs[j,1]

   tempCatchAtAgeN[j,]	<-((vulnN[1,]*inFn)/(vulnN[1,]*inFn+NatMn[1])) * (1-exp(-(vulnN[1,]*inFn+NatMn[1]))) * tempNn[j-1,]
   tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[1,])

   tempCatchAtAgeS[j,]	<-((vulnS[1,]*inFs)/(vulnS[1,]*inFs+NatMs[1])) * (1-exp(-(vulnS[1,]*inFs+NatMs[1]))) * tempNs[j-1,]
   tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[1,])

 }
 SearchYield[p]	<-tempCatchS[j]+tempCatchN[j]
 SearchBiomass[p]	<-EggsN+EggsS
}

trueFMSY	<-SearchFmort[which.max(SearchYield)]
trueUMSY	<- 1-(exp(-trueFMSY))
trueBMSY	<-SearchBiomass[which.max(SearchYield)]
trueMSY	<-max(SearchYield)

dev.new()
par(mfrow=c(1,2))
plot(SearchYield~SearchFmort)
plot(SearchYield~SearchBiomass)
}
#=========================================================================
# INITALIZE THE POPULATION
# Option 1: Set the initial conditions by scaling the specified F series
# that will put the population at the designated depletion 
#
# Option 2: Set depletion to 0 and use the input time series of F to 
# initialize the population
#==========================================================================
if(depletion>0)
{
targetDep	<-(VirBioN+VirBioS)*depletion
maxExp	<-0.999
minExp	<-0.001

#==defines the exploitation pattern
ExpSeries	<-HistoricalF

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

tempRecN	<-rep(0,InitYear)
tempRecS	<-rep(0,InitYear)

tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=InitYear)
tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=InitYear)

tempExp	<-(maxExp+minExp)/2
tempF		<- -log(1-tempExp)
tempFsrsN	<-ExpSeries*tempF
tempFsrsS	<-ExpSeries*tempF

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
HistoricalF<-ExpSeries*tempExp
print(c("Final Depletion",tempDep))
}

if(depletion==0)
{
tempNn	<-matrix(ncol=MaxAge,nrow=InitYear)
tempNn[1,]	<-VirInitN
tempNs	<-matrix(ncol=MaxAge,nrow=InitYear)
tempNs[1,]	<-VirInitS
tempCatchN	<-rep(0,InitYear)
tempCatchS	<-rep(0,InitYear)
tempRecN	<-rep(0,InitYear)
tempRecS	<-rep(0,InitYear)
tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=InitYear)
tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=InitYear)

HistoricalFs<-HistoricalF
HistoricalFn<-HistoricalF

for (j in 2:InitYear)
{
 for (i in 2:(MaxAge-1))
  {
   tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-HistoricalFn[j]*vulnN[1,i-1])*exp(-NatMn[j])
   tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-HistoricalFs[j]*vulnS[1,i-1])*exp(-NatMs[j])

  }
   tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-HistoricalFn[j]*vulnN[1,MaxAge])*exp(-NatMn[j])+ tempNn[j-1,MaxAge]*exp(-HistoricalFn[j]*vulnN[1,MaxAge])*exp(-NatMn[j])
   tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-HistoricalFs[j]*vulnS[1,MaxAge])*exp(-NatMs[j])+ tempNs[j-1,MaxAge]*exp(-HistoricalFs[j]*vulnS[1,MaxAge])*exp(-NatMs[j])

   EggsN			<-sum(tempNn[j-1,]*matureN[1,]*WeightAtAgeN[1,])
   EggsS			<-sum(tempNs[j-1,]*matureS[1,]*WeightAtAgeS[1,])

   tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1,j],recType="BH",NatMin=NatMn[j],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,])
   tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1,j],recType="BH",NatMin=NatMs[j],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,])
   tempRecN[j]		<-tempNn[j,1]
   tempRecS[j]		<-tempNs[j,1]

   tempCatchAtAgeN[j,]	<-((vulnN[1,]*HistoricalFn[j])/(vulnN[1,]*HistoricalFn[j]+NatMn[j])) * (1-exp(-(vulnN[1,]*HistoricalFn[j]+NatMn[j]))) * tempNn[j-1,]
   tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[1,])

   tempCatchAtAgeS[j,]	<-((vulnS[1,]*HistoricalFs[j])/(vulnS[1,]*HistoricalFs[j]+NatMs[j])) * (1-exp(-(vulnS[1,]*HistoricalFs[j]+NatMs[j]))) * tempNs[j-1,]
   tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[1,])

 }
}

if(LifeHistoryPlots==1)
 PlotLifeHistory()

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
 tempCatAtLenN<-matrix(nrow=MaxAge,ncol=MaxAge)
 tempCatAtLenS<-matrix(nrow=MaxAge,ncol=MaxAge)

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



if(AssessmentData==1)
 {
 dev.new()
 par(mfrow=c(ceiling(sqrt(InitYear)),ceiling(sqrt(InitYear))),mar=c(.1,.1,.1,.1))
 for(i in 2:InitYear)
  barplot(rbind(projCatLenFreqN[i,,1],projCatLenFreqS[i,,1]),beside=T,xaxt='n',yaxt='n')
 }

#==============================================================
# calculate the survey proportion at length for assessment
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

if(AssessmentData==1)
 {
 dev.new()
 par(mfrow=c(ceiling(sqrt(InitYear)),ceiling(sqrt(InitYear))),mar=c(.1,.1,.1,.1))
 for(i in 2:InitYear)
  barplot(rbind(projSurvLenFreqN[i,,1],projSurvLenFreqS[i,,1]),beside=T,xaxt='n',yaxt='n')
 }

#==transfer historical time series from above
for(x in 1:Nsim)
{
projNn[1:InitYear,,x]			<-tempNn
projNs[1:InitYear,,x]			<-tempNs
projCatchN[x,1:InitYear]		<-tempCatchN
projCatchS[x,1:InitYear]		<-tempCatchS
projRecN[x,1:InitYear]			<-tempRecN
projRecS[x,1:InitYear]			<-tempRecS
projFmortN[x,1:InitYear]		<-HistoricalF
projFmortS[x,1:InitYear]		<-HistoricalF
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

#==storage for management and assessment quantities
FMSY<-matrix(nrow=Nsim,ncol=SimYear)
BMSY<-matrix(nrow=Nsim,ncol=SimYear)
TAC<-matrix(nrow=Nsim,ncol=SimYear)
CurBio<-matrix(nrow=Nsim,ncol=SimYear)

#=================================================================================
#==BEGIN PROJECTIONS===========================================================

for(z in 1:Nsim)
{
 for(y in (InitYear+1):SimYear)
 {
 #==pull data for assessment
 CatchData<-CatchAssessN[z,!is.na(CatchAssessN[z,])] + CatchAssessS[z,!is.na(CatchAssessS[z,])]
 CPUEData<-CPUEAssessN[z,!is.na(CPUEAssessN[z,])] + CPUEAssessS[z,!is.na(CPUEAssessS[z,])] 
 SurvData<-SurvAssessN[z,!is.na(SurvAssessN[z,])] + SurvAssessS[z,!is.na(SurvAssessS[z,])] 

 #=================================================================================
 #===ASSESS THE STOCK, CHOOSE A METHOD IN CTL FILE=================================
 #=================================================================================

 #==Production model==
 if(AssessmentType == 1)
 {
 if(y==InitYear)
  x		<-c(1.2*(VirBioN+VirBioS),0.3)	#initial values
 if(y>InitYear)
  x		<-outs$par
 outs		<-nlminb(start=x,objective=ProdMod,CatchData=CatchData[2:(y-1)],IndexData=CPUEData[2:(y-1)])
 PredBio	<-ProdModPlot(outs$par,CatchData[2:(y-1)],CPUEData[2:(y-1)],plots=EstimationPlots)
 FMSY[z,y] 	<-abs(outs$par[2])/2
 BMSY[z,y]	<-abs(outs$par[1])/2
 MSY		<-BMSY[z,y]*FMSY[z,y]
 CurBio[z,y]<-PredBio[length(PredBio)-1]
 TAC[z,y]	<-HarvestControlRule(FMSY=FMSY[z,y],BMSY=BMSY[z,y],ExploitBio=CurBio[z,y],SpawnBio=CurBio[z,y],
			alpha=HCalpha,beta=HCbeta,HarvestControl=HarvestControl,ConstantCatch=ConstantCatch,
  			ConstantF=ConstantF)
 }

 #==Age-structured model
 if(AssessmentType == 2)
 {
 #==make folder if it isn't there
 setwd(CurDir)
 dir.create("AgeAssOutputs")
 IndSimFolder<-paste("AgeAssOutputs/",y,sep="")
 dir.create(IndSimFolder)

 #==copy the .exe into it (where does this come from? github?)
 file.copy(from="GenAss/simass.exe",to=IndSimFolder)
 
 #==write the true values
 #==Probably don't need this if it is stored in the MSE object
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
 cat(c(alphaN[1],betaN[1]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#length parameters","\n",file="SimAss.CTL",append=TRUE)
 cat(c(VonKn[1],LinfN[1],t0n[1]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#natural mortality","\n",file="SimAss.CTL",append=TRUE)
 cat(NatMn[length(NatMn)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#maturity","\n",file="SimAss.CTL",append=TRUE)
 cat(c(mat50n[length(mat50n)],mat95n[length(mat95n)]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#steepness","\n",file="SimAss.CTL",append=TRUE)
 cat(steepnessN[length(steepnessN)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#small number","\n",file="SimAss.CTL",append=TRUE)
 cat(SmallNum,"\n",file="SimAss.CTL",append=TRUE)
 cat("#initial smoothness weight","\n",file="SimAss.CTL",append=TRUE)
 cat(InitSmooth,"\n",file="SimAss.CTL",append=TRUE) 
 cat("#Fishery Independent data available?","\n",file="SimAss.CTL",append=TRUE)
 cat(FisheryIndepenDat,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Length at age standard deviation","\n",file="SimAss.CTL",append=TRUE)
 cat(GrowthSDn,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on fishing mortality penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(FmortPen,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on recruit penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(RecruitPen,"\n",file="SimAss.CTL",append=TRUE)

 cat("#harvest control rule selected","\n",file="SimAss.CTL",append=TRUE)
 cat(HarvestControl,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Constantcatch if harvest control ==2, set catch level","\n",file="SimAss.CTL",append=TRUE)
 cat(ConstantCatch,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Constantcatch if harvest control ==3, set fmort level","\n",file="SimAss.CTL",append=TRUE)
 cat(ConstantF,"\n",file="SimAss.CTL",append=TRUE)
 cat("#40 10 parameters for harvest control 4","\n",file="SimAss.CTL",append=TRUE)
 cat(HCalpha,"\n",file="SimAss.CTL",append=TRUE)
 cat("#40 10 parameters for harvest control 4","\n",file="SimAss.CTL",append=TRUE)
 cat(HCbeta,"\n",file="SimAss.CTL",append=TRUE)
 



 #==.DAT file !!!!!!!!!!!!!!!!!FIX THIS!!!!!!!!!!!!!!!!!!!!!!
 file.create("SimAss.DAT")
 cat("#Simulated generalized assessment","\n",file="SimAss.DAT")
 cat("#","\n",file="SimAss.DAT",append=TRUE)
 cat("#start/end year","\n",file="SimAss.DAT",append=TRUE)
 cat(c(1,(y-1)),"\n",file="SimAss.DAT",append=TRUE)

 #==aggregate the length frequencies from the two areas
 cat("#number of ages","\n",file="SimAss.DAT",append=TRUE)
 cat(MaxAge,"\n",file="SimAss.DAT",append=TRUE)
 cat("#ages","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,MaxAge),"\n",file="SimAss.DAT",append=TRUE)

 #==aggregate the index of abundance from two areas
 cat("#CPUE years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-1),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for CPUE","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,(y-1)),"\n",file="SimAss.DAT",append=TRUE)

 cat("#survey years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-1),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for survey","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,(y-1)),"\n",file="SimAss.DAT",append=TRUE)

 cat("#length freq from catch years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-1),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for length freq catch","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,(y-1)),"\n",file="SimAss.DAT",append=TRUE)
 cat("#catch length sample sizes","\n",file="SimAss.DAT",append=TRUE)
 cat(LenSampleN,"\n",file="SimAss.DAT",append=TRUE) 

 cat("#length freq from survey years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-1),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for length freq survey","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,(y-1)),"\n",file="SimAss.DAT",append=TRUE)
 cat("#survey sample sizes","\n",file="SimAss.DAT",append=TRUE)
 cat(LenSampleN,"\n",file="SimAss.DAT",append=TRUE) 

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
 for(i in 1:(y-1))
  cat(projCatLenFreqN[i,,z]+projCatLenFreqS[i,,z],"\n",file="SimAss.DAT",append=TRUE)

 cat("#survey length counts","\n",file="SimAss.DAT",append=TRUE)
 for(i in 1:(y-1))
  cat(projSurvLenFreqN[i,,z]+projSurvLenFreqS[i,,z],"\n",file="SimAss.DAT",append=TRUE)

 #==run the code
 shell("simass -nohess")

#==================================
# PLOT THE ASSESSMENT OUTPUT (if instructed)
#==================================

#IndSimFolder<-"AgeAssOutputs/51"

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

 if(EstimationPlots==1 & AssessmentType == 1)
 {
  lines(BMSY[z,],col=3,lty=3)
  legend("topright",bty='n',col=c(NA,1,2,3),pch=c(15,NA,NA,NA),lty=c(NA,1,2,3),
          legend=c("Index obs","Index est","Catch","BMSY"))
 }
 #==================================================================================
 #==find the F that removes the TAC given the fishery selectivity===================
 #==this TAC is not necessarily the 'true' TAC--there is error in the biomass and FMSY
 #==================================================================================

 if(CatchShareN==0.5)
 {
 minExp<-0.00001
 maxExp<-0.99999

 for(p in 1:30)
 {
 tempExp	<-(maxExp+minExp)/2
 tempF	<- -log(1-tempExp)
 #tempFn	<-tempF*CatchShareN
 #tempFs	<-tempF*CatchShareS
 tempFn	<-tempF
 tempFs	<-tempF

   finCatLenFn	<-((tempFn*vulnN[y,])/(vulnN[y,]*tempFn+NatMn[y])) * (1-exp(-((tempFn*vulnN[y,])+NatMn[y]))) * projNn[y-1,,z]
   finCatFn		<-sum(finCatLenFn*WeightAtAgeN[y,])

   finCatLenFs	<-((tempFs*vulnS[y,])/(vulnS[y,]*tempFs+NatMs[y])) * (1-exp(-((tempFs*vulnS[y,])+NatMs[y]))) * projNs[y-1,,z]
   finCatFs		<-sum(finCatLenFs*WeightAtAgeS[y,])

   TempTotCatch	<-finCatFs+finCatFn

 if(TempTotCatch<TAC[z,y])
  minExp	<-tempExp
 if(TempTotCatch>TAC[z,y])
  maxExp	<-tempExp
 initExp	<-tempExp
 }
 ApplyF	<- -log(1-tempExp)

 #ApplyFn	<-ApplyF*CatchShareN
 #ApplyFs	<-ApplyF*CatchShareS
 ApplyFn	<-ApplyF
 ApplyFs	<-ApplyF
 }

 if(CatchShareN!=0.5)
 {
 }
 #===UPDATE DATA FOR ASSESSMENT====================================
   projCatchAtAgeN[y,,z]	<-((ApplyFn*vulnN[y,])/(vulnN[y,]*ApplyFn+NatMn[y])) * (1-exp(-((ApplyFn*vulnN[y,])+NatMn[y]))) * projNn[y-1,,z]
   projCatchN[z,y]		<-sum(projCatchAtAgeN[y,,z]*WeightAtAgeN[y,])

   projCatchAtAgeS[y,,z]	<-((ApplyFs*vulnS[y,])/(vulnS[y,]*ApplyFs+NatMs[y])) * (1-exp(-((ApplyFs*vulnS[y,])+NatMs[y]))) * projNs[y-1,,z]
   projCatchS[z,y]		<-sum(projCatchAtAgeS[y,,z]*WeightAtAgeS[y,])

   #==UPDATE CATCH DATA AND INDEX (survey and CPUE) DATA
   projSSBn[z,y]	<-sum(projNn[y-1,,z]*matureN[y,]*WeightAtAgeN[y,])
   projSSBs[z,y]	<-sum(projNs[y-1,,z]*matureS[y,]*WeightAtAgeS[y,])
   projExpBn[z,y]	<-sum(projNn[y-1,,z]*vulnN[y,]*WeightAtAgeN[y,])
   projExpBs[z,y]	<-sum(projNs[y-1,,z]*vulnS[y,]*WeightAtAgeS[y,])
   projSurvN[z,y]	<-sum(projNn[y-1,,z]*survSelN[y,]*WeightAtAgeN[y,])
   projSurvS[z,y]	<-sum(projNs[y-1,,z]*survSelS[y,]*WeightAtAgeS[y,])

   CatchAssessN[z,y]	<-projCatchN[z,y]*exp(CatchErrorN[z,y])
   CatchAssessS[z,y]	<-projCatchS[z,y]*exp(CatchErrorS[z,y])

   CPUEAssessN[z,y]	<-projExpBn[z,y]*exp(CPUEErrorN[z,y])
   CPUEAssessS[z,y]	<-projExpBs[z,y]*exp(CPUEErrorS[z,y])

   SurvAssessN[z,y]	<-projSurvN[z,y]*exp(CPUEErrorN[z,y])
   SurvAssessS[z,y]	<-projSurvS[z,y]*exp(CPUEErrorS[z,y])
   
  #==UPDATE CATCH LENGTH FREQS====
    for(w in 1:nrow(tempCatAtLenN))
	 {
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
	  ProbN<-probtemp/sum(probtemp)
	  tempCatAtLenN[w,]<-projCatchAtAgeN[y,w,z]*ProbN
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
	  ProbS<-probtemp/sum(probtemp)
	  tempCatAtLenS[w,]<-projCatchAtAgeS[y,w,z]*ProbS
	 }
	 projCatLenFreqN[y,,z]<-apply(tempCatAtLenN,2,sum)
	 projCatLenFreqS[y,,z]<-apply(tempCatAtLenS,2,sum)

  #==============================================================
  # calculate the survey proportion at age (length) for assessment
  #==============================================================
	 #==make length frequencies for s==
   for(w in 1:nrow(tempCatAtLenN))
    {
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
	ProbN<-probtemp/sum(probtemp)
	tempSurvAtLenN[w,]<-projNn[y-1,w,z]*survSelN[y,w]*ProbN
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
	ProbS<-probtemp/sum(probtemp)
	tempSurvAtLenS[w,]<-projNs[y-1,w,z]*survSelS[y,w]*ProbS
     }

    projSurvLenFreqN[y,,z]<-apply(tempSurvAtLenN,2,sum)
    projSurvLenFreqS[y,,z]<-apply(tempSurvAtLenS,2,sum)
 #==apply F that would take a calculated TAC====================================================
 #==update population dynamics
  for (i in 2:(MaxAge-1))
  {
   projNn[y,i,z]		<-projNn[y-1,i-1,z]*exp(-ApplyFn*vulnN[y,i-1])*exp(-NatMn[y])
   projNs[y,i,z]		<-projNs[y-1,i-1,z]*exp(-ApplyFs*vulnS[y,i-1])*exp(-NatMs[y])
  }
   projNn[y,MaxAge,z]		<-(projNn[y-1,(MaxAge-1),z])*exp(-ApplyFn*vulnN[y,MaxAge])*exp(-NatMn[y])+ projNn[y-1,MaxAge,z]*exp(-ApplyFn*vulnN[y,MaxAge])*exp(-NatMn[y])
   projNs[y,MaxAge,z]		<-(projNs[y-1,(MaxAge-1),z])*exp(-ApplyFs*vulnS[y,MaxAge])*exp(-NatMs[y])+ projNs[y-1,MaxAge,z]*exp(-ApplyFn*vulnS[y,MaxAge])*exp(-NatMs[y])

   EggsN				<-sum(projNn[y-1,,z]*matureN[y,]*WeightAtAgeN[y,])
   EggsS				<-sum(projNs[y-1,,z]*matureS[y,]*WeightAtAgeS[y,])

   projNn[y,1,z]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[y],RzeroIN=RzeroN[y],RecErrIN=RecErrN[z,y],recType="BH",NatMin=NatMn[y],
							vulnIN=vulnN[y,],matureIN=matureN[y,],weightIN=WeightAtAgeN[y,],LenAtAgeIN=LenAtAgeN[y,])
   projNs[y,1,z]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[y],RzeroIN=RzeroS[y],RecErrIN=RecErrS[z,y],recType="BH",NatMin=NatMs[y],
							vulnIN=vulnS[y,],matureIN=matureS[y,],weightIN=WeightAtAgeS[y,],LenAtAgeIN=LenAtAgeS[y,])



 } # end y
} # end z


