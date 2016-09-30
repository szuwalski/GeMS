###################################################################
# Generalized spatial management strategy evaluation
# Two box, populations processes in each box can
# vary over time, movement between boxes can vary over 
# time, effort in boxes can be allocated asymmetrically
# production and age-structured estimation models can be tested
# written by Cody Szuwalski, 11/2015
###################################################################

GeMS<-function(out,CreateFolderName)
{
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
PlotYieldCurve	<-out$OM$PlotYieldCurve

depletion		<-out$OM$depletion		# model is conditioned to a specific depletion given a trajectory of effort
CatchShareN		<-out$OM$CatchShareN		# the amount of effort allocated to a given area
CatchShareS		<-1-CatchShareN

#==sampling uncertainty (change this to be time-varying as well)
CatchCVn	<-CleanInput(out$OM$CatchCVn,SimYear)
CatchCVs	<-CleanInput(out$OM$CatchCVs,SimYear)

IndexCVn	<-CleanInput(out$OM$IndexCVn,SimYear)
IndexCVs	<-CleanInput(out$OM$IndexCVs,SimYear)

LenSampleN	<-CleanInput(out$OM$LenSampleN,SimYear)
LenSampleS	<-CleanInput(out$OM$LenSampleS,SimYear)

GrowthSDn	<-CleanInput(out$OM$GrowthSDn,SimYear)
GrowthSDs	<-CleanInput(out$OM$GrowthSDs,SimYear)

#==========================================================
#=================population dynamics processes============
#==============================================+===========
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
 BinWidth		<-(max(LenAtAgeS,LenAtAgeN)*1.05)/(out$OM$LengthBinN)
 LengthBins		<-seq(0,max(LenAtAgeS,LenAtAgeN)*1.05,BinWidth)
 LengthBinsMid	<-LengthBins[1:(length(LengthBins)-1)] + mean(LengthBins[1:2])
 LengthBinN		<-length(LengthBinsMid)

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

#==needs to be changed to allow for all 
RecErrN	<-matrix(rnorm(SimYear*Nsim,0,sigmaRn),ncol=SimYear,byrow=T)
RecErrS	<-matrix(rnorm(SimYear*Nsim,0,sigmaRs),ncol=SimYear,byrow=T)

#==historic fishing mortality
HistoricalF		<-CleanInput(out$OM$HistoricalF,SimYear)[1:InitYear]
PastFsd		<-out$OM$PastFsd
HarvestControl	<-out$OM$HarvestControl
ConstantCatch	<-out$OM$ConstantCatch
ConstantF		<-out$OM$ConstantF
HCalpha		<-out$OM$HCalpha
HCbeta		<-out$OM$HCbeta

SmallNum		<-out$OM$SmalNum
InitSmooth		<-out$OM$InitSmooth	
FmortPen		<-out$OM$FmortPen
RecruitPen		<-out$OM$RecruitPen	
Mpenalty		<-out$OM$Mpenalty	
Growthpenalty	<-out$OM$Growthpenalty	
SelPenalty		<-out$OM$SelPenalty

EstM			<-out$OM$EstM
TimeVaryM		<-out$OM$TimeVaryM	
EstGrowthK		<-out$OM$EstGrowthK
TimeVaryGrowthK	<-out$OM$TimeVaryGrowthK
EstLinf		<-out$OM$EstLinf
TimeVaryLinf	<-out$OM$TimeVaryLinf
TimeVarySel50	<-out$OM$TimeVarySel50
TimeVarySel95	<-out$OM$TimeVarySel95
ProjectTimeVary	<-out$OM$ProjectTimeVary
InitValSig		<-out$OM$InitValSig

InitBzeroMod	<-out$OM$InitBzeroMod
InitGrowthRate	<-out$OM$InitGrowthRate

#==expand historical fishing mortality==

#===================================================
#==Virgin numbers at age, biomass, initial depletion, recruitment
#===================================================
VirInitN		<-initialN(Rzero=RzeroN[1],NatM=NatMn[1],inAge=MaxAge)
VirInitS		<-initialN(Rzero=RzeroS[1],NatM=NatMs[1],inAge=MaxAge)

VirBioN		<-sum(VirInitN*matureN[1,]*WeightAtAgeN[1,])
VirBioS		<-sum(VirInitS*matureS[1,]*WeightAtAgeS[1,])

ExploitBioN		<-sum(VirInitN*vulnN[1,]*WeightAtAgeN[1,])
ExploitBioS		<-sum(VirInitS*vulnS[1,]*WeightAtAgeS[1,])

#=========================================================================
# INITALIZE THE POPULATION
#==========================================================================
tempNn		<-array(dim=c(InitYear,MaxAge,Nsim))
tempNn[1,,]		<-VirInitN
tempNs		<-array(dim=c(InitYear,MaxAge,Nsim))
tempNs[1,,]		<-VirInitS
tempCatchN		<-matrix(ncol=InitYear,nrow=Nsim)
tempCatchS		<-matrix(ncol=InitYear,nrow=Nsim)
tempRecN		<-matrix(ncol=InitYear,nrow=Nsim)
tempRecS		<-matrix(ncol=InitYear,nrow=Nsim)
tempCatchAtAgeN	<-array(dim=c(InitYear,MaxAge,Nsim))
tempCatchAtAgeS	<-array(dim=c(InitYear,MaxAge,Nsim))

#==Make a matrix of past Fs based on a sd and an input time series==
HistoricalFs<-matrix(HistoricalF,ncol=InitYear,nrow=Nsim,byrow=T)
HistoricalFn<-matrix(HistoricalF,ncol=InitYear,nrow=Nsim,byrow=T)
Ferr<-matrix(rnorm(InitYear*Nsim,1,PastFsd),ncol=InitYear)
HistoricalFs<-HistoricalFs*Ferr
HistoricalFn<-HistoricalFn*Ferr

for(k in 1:Nsim)
{
 for (j in 2:InitYear)
 {
  for (i in 2:(MaxAge-1))
  {
   tempNn[j,i,k]		<-tempNn[j-1,i-1,k]*exp(-HistoricalFn[k,j]*vulnN[1,i-1])*exp(-NatMn[j])
   tempNs[j,i,k]		<-tempNs[j-1,i-1,k]*exp(-HistoricalFs[k,j]*vulnS[1,i-1])*exp(-NatMs[j])
  }
   #tempNn[,,k]
   tempNn[j,MaxAge,k]	<-(tempNn[j-1,(MaxAge-1),k])*exp(-HistoricalFn[k,j]*vulnN[j,MaxAge])*exp(-NatMn[j])+ tempNn[j-1,MaxAge,k]*exp(-HistoricalFn[k,j]*vulnN[j,MaxAge])*exp(-NatMn[j])
   tempNs[j,MaxAge,k]	<-(tempNs[j-1,(MaxAge-1),k])*exp(-HistoricalFs[k,j]*vulnS[j,MaxAge])*exp(-NatMs[j])+ tempNs[j-1,MaxAge,k]*exp(-HistoricalFs[k,j]*vulnS[j,MaxAge])*exp(-NatMs[j])

   EggsN			<-sum(tempNn[j-1,,k]*matureN[j,]*WeightAtAgeN[j,])
   EggsS			<-sum(tempNs[j-1,,k]*matureS[j,]*WeightAtAgeS[j,])

   tempNn[j,1,k]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[j],RzeroIN=RzeroN[j],RecErrIN=RecErrN[k,j],recType="BH",NatMin=NatMn[j],
							vulnIN=vulnN[j,],matureIN=matureN[j,],weightIN=WeightAtAgeN[j,],LenAtAgeIN=LenAtAgeN[j,],MaxAge=MaxAge,sigmaRin=sigmaRn[j])
   tempNs[j,1,k]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[j],RzeroIN=RzeroS[j],RecErrIN=RecErrS[k,j],recType="BH",NatMin=NatMs[j],
							vulnIN=vulnS[j,],matureIN=matureS[j,],weightIN=WeightAtAgeS[j,],LenAtAgeIN=LenAtAgeS[j,],MaxAge=MaxAge,sigmaRin=sigmaRs[j])
   tempRecN[k,j]		<-tempNn[j,1,k]
   tempRecS[k,j]		<-tempNs[j,1,k]

   tempCatchAtAgeN[j,,k]<-((vulnN[j,]*HistoricalFn[k,j])/(vulnN[j,]*HistoricalFn[k,j]+NatMn[j])) * (1-exp(-(vulnN[j,]*HistoricalFn[k,j]+NatMn[j]))) * tempNn[j-1,,k]
   tempCatchN[k,j]	<-sum(tempCatchAtAgeN[j,,k]*WeightAtAgeN[j,])

   tempCatchAtAgeS[j,,k]<-((vulnS[j,]*HistoricalFs[k,j])/(vulnS[j,]*HistoricalFs[k,j]+NatMs[j])) * (1-exp(-(vulnS[j,]*HistoricalFs[k,j]+NatMs[j]))) * tempNs[j-1,,k]
   tempCatchS[k,j]	<-sum(tempCatchAtAgeS[j,,k]*WeightAtAgeS[j,])
 }
}
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

projCatLenFreqN	<-array(dim=c(SimYear,LengthBinN,Nsim))
projCatLenFreqS	<-array(dim=c(SimYear,LengthBinN,Nsim))
projSurvLenFreqN	<-array(dim=c(SimYear,LengthBinN,Nsim))
projSurvLenFreqS	<-array(dim=c(SimYear,LengthBinN,Nsim))

#=storage for the true quantities
trueRec		<-matrix(ncol=SimYear,nrow=Nsim)
trueCatch		<-matrix(ncol=SimYear,nrow=Nsim)
trueFmort		<-matrix(ncol=SimYear,nrow=Nsim)
for(x in 1:Nsim)
 trueFmort[x,1:InitYear]<-HistoricalF

trueSpbio		<-matrix(ncol=SimYear,nrow=Nsim)
trueSurvInd		<-matrix(ncol=SimYear,nrow=Nsim)
trueCPUEind		<-matrix(ncol=SimYear,nrow=Nsim)

#========================================================================
# FIND REFERENCE POINTS FOR THE POPULATION 
#========================================================================
if(PlotYieldCurve==1)
{
SearchFmort		<-seq(0.01,3*NatMn[1],(NatMn[1]-0.01)/100)
SearchYield		<-rep(0,length(SearchFmort))
SearchBiomass	<-rep(0,length(SearchFmort))
for(p in 1:length(SearchFmort))
{
tempOut<-ProjPopDym(SearchFmort[p])
SearchYield[p]<-tempOut[[1]]
SearchBiomass[p]<-tempOut[[2]]
}
dev.new()
par(mfrow=c(1,2))
plot(SearchYield~SearchFmort)
plot(SearchYield~SearchBiomass)
}

#============================================
# calculate 'true' reference points
# figure out a way to determine when to do this for all sims and when for just 1
# as long as population parametesr are from the same control file, this should be the same 
# for all populations regardless of fishing history or recruitment
#===================================
FindFMSYin<-1
if(FindFMSYin==1)
{
x<-0.3
FmsyOut	<-nlminb(x,FindFMSY,MaxAge=MaxAge,VirInitN=VirInitN,VirInitS=VirInitS,vulnN=vulnN,
                     	vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
				WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
				steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,LenAtAgeS=LenAtAgeS,sigmaRn=sigmaRn,sigmaRs=sigmaRs)
trueFMSY	<-FmsyOut$par
trueUMSY	<- 1-(exp(-trueFMSY))
trueBMSY	<-ProjPopDym(fmort=trueFMSY,MaxAge=MaxAge,vulnN=vulnN,
                     	vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
				WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
				steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,LenAtAgeS=LenAtAgeS,sigmaRn=sigmaRn,sigmaRs=sigmaRs)[[2]]
trueMSY	<- -FmsyOut$objective
}

#======================================
# calculate reference point proxies
#======================================
ConstRec	<-100000
outsF35	<-FindF35(MaxAge=MaxAge,vulnN=vulnN,
                     	vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
				WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
				steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,
				LenAtAgeS=LenAtAgeS,inRec=ConstRec)

 
trueF35	<-outsF35[[1]]
trueB35	<-ProjPopDym(fmort=trueFMSY,MaxAge=MaxAge,vulnN=vulnN,
                     	vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
				WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
				steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,
				LenAtAgeS=LenAtAgeS,sigmaRn=sigmaRn,sigmaRs=sigmaRs)[[2]]


#==============================================================
# calculate the catch proportion at (length) for assessment
#==============================================================
 tempCatAtLenN<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))
 tempCatAtLenS<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))

for(x in 1:Nsim)
 for(y in 2:InitYear)
 {
 #==make length frequencies for catch==
 for(w in 1:ncol(LenAtAgeN))
 {
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
  ProbN<-probtemp/sum(probtemp)
  tempCatAtLenN[w,,x]<-tempCatchAtAgeN[y,w,x]*ProbN
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
  ProbS<-probtemp/sum(probtemp)
  tempCatAtLenS[w,,x]<-tempCatchAtAgeS[y,w,x]*ProbS
 }
 projCatLenFreqN[y,,x]<-apply(tempCatAtLenN[,,x],2,sum)
 projCatLenFreqS[y,,x]<-apply(tempCatAtLenS[,,x],2,sum)
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
 tempSurvAtLenN<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))
 tempSurvAtLenS<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))

for(x in 1:Nsim)
 for(y in 2:InitYear)
 {
 #==make length frequencies for s==
 for(w in 1:ncol(LenAtAgeS))
 {
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
  ProbN<-probtemp/sum(probtemp)
  tempSurvAtLenN[w,,x]<-tempNn[y,w,x]*survSelN[y,w]*ProbN
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
  ProbS<-probtemp/sum(probtemp)
  tempSurvAtLenS[w,,x]<-tempNs[y,w,x]*survSelS[y,w]*ProbS
 }
 projSurvLenFreqN[y,,x]<-apply(tempSurvAtLenN[,,x],2,sum)
 projSurvLenFreqS[y,,x]<-apply(tempSurvAtLenS[,,x],2,sum)
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
projNn[1:InitYear,,x]			<-tempNn[,,x]
projNs[1:InitYear,,x]			<-tempNs[,,x]
projCatchN[x,1:InitYear]		<-tempCatchN[x,]
projCatchS[x,1:InitYear]		<-tempCatchS[x,]
projRecN[x,1:InitYear]			<-tempRecN[x,]
projRecS[x,1:InitYear]			<-tempRecS[x,]
projFmortN[x,1:InitYear]		<-HistoricalFn[x,]
projFmortS[x,1:InitYear]		<-HistoricalFs[x,]
projSurvN[x,1:InitYear]			<-apply(tempNn[,,x]*survSelN[1:InitYear,]*WeightAtAgeN[1:InitYear,],1,sum)
projSurvS[x,1:InitYear]			<-apply(tempNs[,,x]*survSelS[1:InitYear,]*WeightAtAgeS[1:InitYear,],1,sum)

for(y in 1:InitYear)
{
 projSSBn[x,y]	<-sum(projNn[y,,x]*matureN[y,]*WeightAtAgeN[y,])
 projSSBs[x,y]	<-sum(projNs[y,,x]*matureS[y,]*WeightAtAgeS[y,])
 projExpBn[x,y]	<-sum(projNn[y,,x]*vulnN[y,]*WeightAtAgeN[y,])
 projExpBs[x,y]	<-sum(projNs[y,,x]*vulnS[y,]*WeightAtAgeS[y,])
}
}


#==fill in true storage arrays
for(x in 1:Nsim)
{
trueRec[x,1:InitYear]		<-tempRecN[x,]+tempRecS	[x,]	
trueCatch[x,1:InitYear]		<-tempCatchN[x,]+tempCatchS[x,]		
trueSpbio[x,1:InitYear]		<-projSSBn[x,1:InitYear]+projSSBs[x,1:InitYear]	
trueSurvInd[x,1:InitYear]	<-projSurvN[x,1:InitYear]+projSurvS[x,1:InitYear]
trueCPUEind[x,1:InitYear]	<-projExpBn[x,1:InitYear]+projExpBs[x,1:InitYear]	
}

CatchErrorN		<-matrix(rnorm(Nsim*SimYear,0,CatchCVn),nrow=Nsim,ncol=SimYear)
CatchErrorS		<-matrix(rnorm(Nsim*SimYear,0,CatchCVs),nrow=Nsim,ncol=SimYear)
CatchAssessN	<-projCatchN*exp(CatchErrorN-(CatchCVn^2/2))
CatchAssessS	<-projCatchS*exp(CatchErrorS-(CatchCVs^2/2))

CPUEErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
CPUEErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
CPUEAssessN		<-projExpBn*exp(CPUEErrorN-(IndexCVn^2/2))
CPUEAssessS		<-projExpBs*exp(CPUEErrorS-(IndexCVs^2/2))

SurvErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
SurvErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
SurvAssessN		<-projSurvN*exp(SurvErrorN-(IndexCVn^2/2))
SurvAssessS		<-projSurvS*exp(SurvErrorS-(IndexCVs^2/2))

#==storage for management and assessment quantities
FMSY		<-matrix(nrow=Nsim,ncol=SimYear)
BMSY		<-matrix(nrow=Nsim,ncol=SimYear)
TAC		<-matrix(nrow=Nsim,ncol=SimYear)
trueTAC	<-matrix(nrow=Nsim,ncol=SimYear)
trueOFL	<-matrix(nrow=Nsim,ncol=SimYear)
CurBio	<-matrix(nrow=Nsim,ncol=SimYear)
EstBio	<-matrix(nrow=Nsim,ncol=SimYear)
Converge	<-matrix(nrow=Nsim,ncol=SimYear)

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
 setwd(CurDir)
 dir.create(CreateFolderName)
 setwd(CreateFolderName)
# InitBzeroMod<-1
# InitGrowthRate<-.2
 if(y==(InitYear+1))
  x		<-c(InitBzeroMod*(VirBioN+VirBioS),InitGrowthRate)	#initial values
 if(y>(InitYear+1))
  x		<-outs$par
 outs		<-nlminb(start=x,objective=ProdMod,CatchData=CatchData[2:(y-1)],IndexData=CPUEData[2:(y-1)])
 Converge[z,y]<-outs$convergence
 PredBio	<-ProdModPlot(outs$par,CatchData[2:(y-1)],CPUEData[2:(y-1)],plots=EstimationPlots)
 FMSY[z,y] 	<-abs(outs$par[2])/2
 BMSY[z,y]	<-abs(outs$par[1])/2
 MSY		<-BMSY[z,y]*FMSY[z,y]
 CurBio[z,y]<-PredBio[length(PredBio)-1]
 inBio	<-projExpBn[z,y-1]+projExpBs[z,y-1]
 TAC[z,y]	<-HarvestControlRule(FMSY=FMSY[z,y],BMSY=BMSY[z,y],ExploitBio=CurBio[z,y],SpawnBio=CurBio[z,y],
			alpha=HCalpha,beta=HCbeta,HarvestControl=HarvestControl,ConstantCatch=ConstantCatch,
  			ConstantF=ConstantF)
 trueTAC[z,y]<-HarvestControlRule(FMSY=trueFMSY,BMSY=trueBMSY,ExploitBio=inBio,SpawnBio=inBio,
			alpha=HCalpha,beta=HCbeta,HarvestControl=HarvestControl,ConstantCatch=ConstantCatch,
  			ConstantF=ConstantF)
 if(y==SimYear)
 {
  EstBio[z,1:length(PredBio)]<-PredBio
 }
 if(y==SimYear & z==Nsim)
 {
 inFile<-paste("ProdOutputs.csv",sep="")
 file.create(inFile)

 cat("#True quantities","\n",file=inFile)
 cat("#","\n",file=inFile,append=TRUE)
 cat("#true Catch","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueCatch[m,],"\n",file=inFile,append=TRUE)
 cat("#true fishing mortality","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueFmort[m,],"\n",file=inFile,append=TRUE)
 cat("#true spawning biomass","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueSpbio[m,],"\n",file=inFile,append=TRUE)
 cat("#true CPUE","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueCPUEind[m,],"\n",file=inFile,append=TRUE)
 cat("#true total allowable catch","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueTAC[m,],"\n",file=inFile,append=TRUE)
 cat("#true BMSY","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueBMSY,"\n",file=inFile,append=TRUE)
 cat("#true FMSY","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueFMSY,"\n",file=inFile,append=TRUE)
 cat("#est cpue ind","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(EstBio[m,],"\n",file=inFile,append=TRUE)
 cat("#est spawning biomass by year (terminal year)","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(CurBio[m,],"\n",file=inFile,append=TRUE)
 cat("#est total allowable catch","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(TAC[m,],"\n",file=inFile,append=TRUE)
 cat("#est BMSY","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(BMSY[m,],"\n",file=inFile,append=TRUE)
 cat("#est FMSY","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(FMSY[m,],"\n",file=inFile,append=TRUE)

 cat("#data cpue","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(CPUEData,"\n",file=inFile,append=TRUE)

 }
}

#==Projection based model--no assessment
 if(AssessmentType == 0)
 {
 #==make folder if it isn't there
 setwd(CurDir)
 dir.create(CreateFolderName)
 dir.create(paste(CreateFolderName,"/",z,sep=""))
 IndSimFolder<-paste(CreateFolderName,"/",z,"/",y,sep="")
 dir.create(IndSimFolder)

 #==write the true values
 setwd(IndSimFolder)
 file.create("TrueQuantities.DAT")
 cat("#True quantities","\n",file="TrueQuantities.DAT")
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#recruitment","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueRec[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#Catch","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueCatch[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#fishing mortality","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueFmort[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#spawning biomass","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueSpbio[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#survey index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueSurvInd[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#cpue index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueCPUEind[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 #cat("#total allowable catch","\n",file="TrueQuantities.DAT",append=TRUE)
 #cat(trueTAC[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#BMSY","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueBMSY,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#FMSY","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueFMSY,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#B35","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueB35,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#F35","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueF35,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#OFL","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueOFL[z,],"\n",file="TrueQuantities.DAT",append=TRUE)

 #==calculate the TAC based on a SRR or a proxy, depending on data
 #==allow this to be a user defined function of any of the things that can be input
 TAC[z,y]	<-trueSpbio[z,y-1]*(1-exp(-ConstantF))
 }

 #==Age-structured assessment model for setting the TAC
 if(AssessmentType == 2)
 {
 #==make folder if it isn't there
 setwd(CurDir)
 dir.create(CreateFolderName)
 dir.create(paste(CreateFolderName,"/",z,sep=""))
 IndSimFolder<-paste(CreateFolderName,"/",z,"/",y,sep="")
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
 for(z in 1:Nsim)
 cat(trueRec[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#Catch","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueCatch[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#fishing mortality","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueFmort[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#spawning biomass","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueSpbio[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#survey index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueSurvInd[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#cpue index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(z in 1:Nsim)
 cat(trueCPUEind[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 #cat("#total allowable catch","\n",file="TrueQuantities.DAT",append=TRUE)
 #cat(trueTAC[z,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#BMSY","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueBMSY,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#FMSY","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueFMSY,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#B35","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueB35,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#F35","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueF35,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#OFL","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueOFL[z,],"\n",file="TrueQuantities.DAT",append=TRUE)

 #selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
 #selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])
 #selAtAgeFunc(surv50n[1],VonKn[1],LinfN[1],t0n[1])
 #selAtAgeFunc(surv95n[1],VonKn[1],LinfN[1],t0n[1])

 #==.CTL file  !!!!!!!!!!!!!!!!!how should we deal with space?!!!!!!!!!!!!!!!!!!!!!!
 file.create("SimAss.CTL")
 cat("#Simulated assessment control file","\n",file="SimAss.CTL")
 cat("#","\n",file="SimAss.CTL",append=TRUE)
 cat("#weight parameters","\n",file="SimAss.CTL",append=TRUE)
 cat(c(alphaN[1],betaN[1]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#how many years to average over timevarying process for projectsion","\n",file="SimAss.CTL",append=TRUE)
 cat(ProjectTimeVary,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Estimate growth? if positive, indicates phase, no estimate if negative","\n",file="SimAss.CTL",append=TRUE)
 cat(EstGrowthK ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#growth time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVaryGrowthK ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Estimate growth? if positive, indicates phase, no estimate if negative","\n",file="SimAss.CTL",append=TRUE)
 cat(EstLinf ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#growth time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVaryLinf ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#length parameters","\n",file="SimAss.CTL",append=TRUE)
 cat(c(VonKn[1],LinfN[1],t0n[1]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#estimate natural mortality? if positive, number indicates phase, if negative, not estimated","\n",file="SimAss.CTL",append=TRUE)
 cat(EstM,"\n",file="SimAss.CTL",append=TRUE)
 cat("#natural mortality (used as starting value if estimated)","\n",file="SimAss.CTL",append=TRUE)
 cat(NatMn[length(NatMn)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#Allow M to vary over time? 1 = yes, 0 = no","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVaryM,"\n",file="SimAss.CTL",append=TRUE) 
 cat("#selectivity 50% time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVarySel50 ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#selectivity 95% time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVarySel95 ,"\n",file="SimAss.CTL",append=TRUE)
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
 cat(GrowthSDn[length(steepnessN)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on fishing mortality penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(FmortPen,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on recruit penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(RecruitPen,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on time-varying natural mortality penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(Mpenalty,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on time-varying growth penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(Growthpenalty,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on time-varying selelctivity penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(SelPenalty,"\n",file="SimAss.CTL",append=TRUE)

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

#==write .PIN file to get initial values close for estimation
 cat("# Simulated assessment pin file","\n",file="SimAss.PIN")
 cat("#","\n",file="SimAss.PIN",append=TRUE)

 InputSurvSel50S<-selAtAgeFunc(surv50s,VonKs,LinfS,t0s)
 InputSurvSel50N<-selAtAgeFunc(surv50n,VonKn,LinfN,t0n)
 InputSurvSel95S<-selAtAgeFunc(surv95s,VonKs,LinfS,t0s)
 InputSurvSel95N<-selAtAgeFunc(surv95n,VonKn,LinfN,t0n)

 cat("# srv_sel50","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(c(InputSurvSel50S,InputSurvSel50N))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# srv_sel95","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(c(InputSurvSel95S,InputSurvSel95N))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# stNatLen","\n",file="SimAss.PIN",append=TRUE)
 cat(log(VirInitN+VirInitS)*rnorm(length(log(VirInitN+VirInitS)),1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# log_avg_fmort_dir","\n",file="SimAss.PIN",append=TRUE)
 cat((log(mean(trueFmort,na.rm=T)))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# fmort_dir_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-1),"\n",file="SimAss.PIN",append=TRUE)
 cat("# mean_log_rec","\n",file="SimAss.PIN",append=TRUE)
 cat((log(mean(trueRec,na.rm=T)))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# rec_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-1),"\n",file="SimAss.PIN",append=TRUE)
 cat("# log_avg_NatM","\n",file="SimAss.PIN",append=TRUE)
 cat((log(mean(c(NatMs,NatMn),na.rm=T)))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# rec_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-1),"\n",file="SimAss.PIN",append=TRUE)

 cat("# log_avg_GrowthK","\n",file="SimAss.PIN",append=TRUE)
 cat((log(mean(c(VonKn,VonKs),na.rm=T)))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# log_avg_Linf","\n",file="SimAss.PIN",append=TRUE)
 cat((log(mean(c(LinfN,LinfS),na.rm=T)))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)

 cat("# GrowthK_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-1),"\n",file="SimAss.PIN",append=TRUE)
 cat("# Linf_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-1),"\n",file="SimAss.PIN",append=TRUE)

 InputSel50S<-selAtAgeFunc(sel50s,VonKs,LinfS,t0s)
 InputSel50N<-selAtAgeFunc(sel50n,VonKn,LinfN,t0n)
 InputSel95S<-selAtAgeFunc(sel95s,VonKs,LinfS,t0s)
 InputSel95N<-selAtAgeFunc(sel95n,VonKn,LinfN,t0n)

 cat("# SelPars50","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(c(InputSel50S,InputSel50N),na.rm=T)*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# SelPars95","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(c(InputSel95S,InputSel95N),na.rm=T)*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)

 cat("# SelPars_dev50","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-1),"\n",file="SimAss.PIN",append=TRUE)
 cat("# SelPars_dev95","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-1),"\n",file="SimAss.PIN",append=TRUE)
 	
 #==.DAT file 
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
 cat("#Number of length bins","\n",file="SimAss.DAT",append=TRUE)
 cat(LengthBinN,"\n",file="SimAss.DAT",append=TRUE)

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
 cat(LenSampleN[seq(1,(y-1))],"\n",file="SimAss.DAT",append=TRUE) 

 cat("#length freq from survey years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-1),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for length freq survey","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,(y-1)),"\n",file="SimAss.DAT",append=TRUE)
 cat("#survey sample sizes","\n",file="SimAss.DAT",append=TRUE)
 cat(LenSampleN[seq(1,(y-1))],"\n",file="SimAss.DAT",append=TRUE) 

 cat("#total survey biomass","\n",file="SimAss.DAT",append=TRUE) 
 cat(SurvData,"\n",file="SimAss.DAT",append=TRUE)
 cat("#survey CV","\n",file="SimAss.DAT",append=TRUE)
 cat(IndexCVn[seq(1,(y-1))],"\n",file="SimAss.DAT",append=TRUE)

 cat("#total CPUE","\n",file="SimAss.DAT",append=TRUE) 
 cat(CPUEData,"\n",file="SimAss.DAT",append=TRUE)
 cat("#CPUE CV","\n",file="SimAss.DAT",append=TRUE)
 cat(IndexCVn[seq(1,(y-1))],"\n",file="SimAss.DAT",append=TRUE)

 cat("#total catch biomass","\n",file="SimAss.DAT",append=TRUE)
 cat(CatchData,"\n",file="SimAss.DAT",append=TRUE)
 cat("#catch CV","\n",file="SimAss.DAT",append=TRUE)
 cat(CatchCVn[seq(1,(y-1))],"\n",file="SimAss.DAT",append=TRUE)
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

 ApplyFn	<-ApplyF
 ApplyFs	<-ApplyF
 }
 trueFmort[z,y]	<-ApplyF

 if(CatchShareN!=0.5)
 {
 #==FIGURE THIS OUT SOON==
 #==i THINK JUST MAKE there be an area that has an F designated relative
 #==to the other area, so the user inputs a fraction of the fishing pressure
 #==that is goign on in one area
 #==0.0000001, would be an MPA, for example
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

   trueCatch[z,y]		<-projCatchN[z,y]+projCatchS[z,y]
   CatchAssessN[z,y]	<-projCatchN[z,y]*exp(CatchErrorN[z,y]-(CatchCVn[y]^2/2))
   CatchAssessS[z,y]	<-projCatchS[z,y]*exp(CatchErrorS[z,y]-(CatchCVs[y]^2/2))

   trueCPUEind[z,y]	<-projExpBn[z,y]+projExpBs[z,y]
   CPUEAssessN[z,y]	<-projExpBn[z,y]*exp(CPUEErrorN[z,y]-(IndexCVn[y]^2/2))
   CPUEAssessS[z,y]	<-projExpBs[z,y]*exp(CPUEErrorS[z,y]-(IndexCVs[y]^2/2))

   trueSurvInd[z,y]	<-projSurvN[z,y]+projSurvS[z,y]
   SurvAssessN[z,y]	<-projSurvN[z,y]*exp(SurvErrorN[z,y]-(IndexCVn[y]^2/2))
   SurvAssessS[z,y]	<-projSurvS[z,y]*exp(SurvErrorS[z,y]-(IndexCVs[y]^2/2))
   
   trueSpbio[z,y]		<-projSSBn[z,y]+projSSBs[z,y]

  #==UPDATE CATCH LENGTH FREQS====
    for(w in 1:ncol(LenAtAgeN))
	 {
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
	  ProbN<-probtemp/sum(probtemp)
	  tempCatAtLenN[w,,z]<-projCatchAtAgeN[y,w,z]*ProbN
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
	  ProbS<-probtemp/sum(probtemp)
	  tempCatAtLenS[w,,z]<-projCatchAtAgeS[y,w,z]*ProbS
	 }
	 projCatLenFreqN[y,,z]<-apply(tempCatAtLenN[,,z],2,sum)
	 projCatLenFreqS[y,,z]<-apply(tempCatAtLenS[,,z],2,sum)

  #==============================================================
  # calculate the survey proportion at age (length) for assessment
  #==============================================================
	 #==make length frequencies for s==
   for(w in 1:ncol(LenAtAgeS))
    {
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
	ProbN<-probtemp/sum(probtemp)
	tempSurvAtLenN[w,,z]<-projNn[y-1,w,z]*survSelN[y,w]*ProbN
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
	ProbS<-probtemp/sum(probtemp)
	tempSurvAtLenS[w,,z]<-projNs[y-1,w,z]*survSelS[y,w]*ProbS
     }

    projSurvLenFreqN[y,,z]<-apply(tempSurvAtLenN[,,z],2,sum)
    projSurvLenFreqS[y,,z]<-apply(tempSurvAtLenS[,,z],2,sum)

   # sum(projSurvLenFreqN[y,,z])
   # sum(projNn[y,,z])

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
							vulnIN=vulnN[y,],matureIN=matureN[y,],weightIN=WeightAtAgeN[y,],LenAtAgeIN=LenAtAgeN[y,],MaxAge=MaxAge,sigmaRin=sigmaRn[y])
   projNs[y,1,z]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[y],RzeroIN=RzeroS[y],RecErrIN=RecErrS[z,y],recType="BH",NatMin=NatMs[y],
							vulnIN=vulnS[y,],matureIN=matureS[y,],weightIN=WeightAtAgeS[y,],LenAtAgeIN=LenAtAgeS[y,],MaxAge=MaxAge,sigmaRin=sigmaRs[y])
   trueRec[z,y]		<-projNn[y,1,z]+projNs[y,1,z]

 #==calculate true OFL
 #==find FOFL given spawning biomass
  tempSpBio			<-EggsN+EggsS
  FutMort			<-trueF35
  if(tempSpBio<trueB35)
  {
   FutMort<-0
   if(tempSpBio>HCbeta*trueB35)
    FutMort = trueF35*(tempSpBio/trueB35-HCalpha)/(1-HCalpha);
   }
  #==find the catch for the FOFL
   trueCatchAtAgeN	<-((vulnN[y,]*FutMort)/(vulnN[y,]*FutMort+NatMn[y])) * (1-exp(-(vulnN[y,]*FutMort+NatMn[y]))) * projNn[y-1,,z]
   trueCatchN		<-sum(trueCatchAtAgeN*WeightAtAgeN[y,])

   trueCatchAtAgeS	<-((vulnS[y,]*FutMort)/(vulnS[y,]*FutMort+NatMs[y])) * (1-exp(-(vulnS[y,]*FutMort+NatMs[y]))) * projNs[y-1,,z]
   trueCatchS		<-sum(trueCatchAtAgeS*WeightAtAgeS[y,])
   trueOFL[z,y]		<-trueCatchS+trueCatchN

  print(paste("Year ",y," of ",SimYear," in simulation ",z," of ",Nsim,sep=""))

 } # end y
} # end z

} # end of GeMS


