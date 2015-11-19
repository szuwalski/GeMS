###################################################
# Generalized spatial management strategy evaluation
# Two box, populations processes in each box can
# vary over time, movement between boxes can vary over 
# time, effort in boxes can be allocated asymmetrically
# written by Cody Szuwalski, 11/2015

#==source a file with helper functions
library(devtools)
inURL<-"https://raw.githubusercontent.com/szuwalski/General-MSE/master/General%20MSE%20helpers.R"
source_url(inURL)

#=======================
#==simulation controls==
#=======================
Nsim			<-1	# number of simulations to do in the MSE
SimYear		<-100	# total number of years in simulation
InitYear		<-50	# year in which MSE starts (i.e. the number of years of data available)

LifeHistoryPlots	<-1
 if(LifeHistoryPlots==1)
 {
  dev.new()
  par(mfrow=c(4,4),mar=c(2,3,2,2))
 }

EstimationPlots	<-1

depletion	<-0.3
CatchShareN	<-0.5
CatchShareS	<-1-CatchShareN

#==sampling uncertainty
CatchCVn	<-0.1
CatchCVs	<-0.1

IndexCVn	<-0.2
IndexCVs	<-0.2

LenSampleN	<-200
LenSampleS	<-200

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
 mtext(side=3,"Maturity at length (N)",cex=.7)
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
 mtext(side=3,"Maturity at length (S)",cex=.7)
 }

#=====================================
#==fishery selectivitiy============== 
#==north==
sel50n	<-rep(50,SimYear)
sel95n	<-rep(150,SimYear)
vulnN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50n[i])/(sel95n[i]-sel50n[i])))

if(LifeHistoryPlots==1)
 {
 plot(vulnN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(vulnN[x,])
 mtext(side=3,"Fishery selectivity (N)",cex=.7)

 }

#==south==
sel50s	<-rep(50,SimYear)
sel95s	<-rep(150,SimYear)
vulnS		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50s[i])/(sel95s[i]-sel50s[i])))

if(LifeHistoryPlots==1)
 {
 plot(vulnS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(vulnS[x,])
 mtext(side=3,"Fishery selectivity (S)",cex=.7)
 }

#==index selectivity=====================
#==north==
surv50n	<-rep(50,SimYear)
surv95n	<-rep(150,SimYear)
survSelN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50n[i])/(surv95n[i]-surv50n[i])))

if(LifeHistoryPlots==1)
 {
 plot(survSelN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(survSelN[x,])
 mtext(side=3,"Index selectivity (N)",cex=.7)
 }

#==south==
surv50s	<-rep(50,SimYear)
surv95s	<-rep(150,SimYear)
survSelS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50s[i])/(sel95s[i]-sel50s[i])))

if(LifeHistoryPlots==1)
 {
 plot(survSelS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(survSelS[x,])
 mtext(side=3,"Index selectivity (S)",cex=.7)
 }

#==weight at length==========================
alphaN		<-0.0000017
betaN			<-3
WeightAtLenN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtLenN[i,]	<-alphaN*LenAtAgeN[i,]^betaN

if(LifeHistoryPlots==1)
 {
 plot(WeightAtLenN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(WeightAtLenN[x,])
 mtext(side=3,"Weight at length (N)",cex=.7)
 }

alphaS		<-0.0000017
betaS			<-3
WeightAtLenS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtLenS[i,]	<-alphaS*LenAtAgeS[i,]^betaS

if(LifeHistoryPlots==1)
 {
 plot(WeightAtLenS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(WeightAtLenS[x,])
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
steepnessN	<-rep(.5,SimYear)
sigmaRn	<-rep(.001,SimYear)
RzeroN	<-rep(1,SimYear)

steepnessS	<-rep(.5,SimYear)
sigmaRs	<-rep(.001,SimYear)
RzeroS	<-rep(1,SimYear)

RecErrN  <-matrix(rnorm(SimYear*Nsim,0,sigmaRn[1]),ncol=SimYear)
RecErrS	<-matrix(rnorm(SimYear*Nsim,0,sigmaRs[1]),ncol=SimYear)


#===================================================
#==Virgin numbers at age, biomass, initial depletion, recruitment
#===================================================
source("C:/Users/Cody/Desktop/BlueShark/Helpers.R")
VirInitN<-initialN(Rzero=RzeroN[1],NatM=NatMn[1],inAge=MaxAge)
VirInitS<-initialN(Rzero=RzeroS[1],NatM=NatMs[1],inAge=MaxAge)

VirBioN<-sum(VirInitN*matureN[1,]*LenAtAgeN[1,]*WeightAtLenN[1,])
VirBioS<-sum(VirInitS*matureS[1,]*LenAtAgeS[1,]*WeightAtLenS[1,])

ExploitBioN<-sum(VirInitN*vulnN[1,]*LenAtAgeN[1,]*WeightAtLenN[1,])
ExploitBioS<-sum(VirInitS*vulnS[1,]*LenAtAgeS[1,]*WeightAtLenS[1,])

if(LifeHistoryPlots==1)
 {
ssb<-seq(1,VirBioN,VirBioN/100)
record<-ssb
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtLenN[1,],LenAtAgeIN=LenAtAgeN[1,])
   }
plot(record~ssb)
mtext(side=3,"SSB vs Rec (N)",cex=.7)
ssb<-seq(1,VirBioS,VirBioS/100)
record<-ssb
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtLenS[1,],LenAtAgeIN=LenAtAgeS[1,])
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

tempCatchAtLenN	<-matrix(ncol=MaxAge,nrow=InitYear)
tempCatchAtLenS	<-matrix(ncol=MaxAge,nrow=InitYear)

tempExp	<-(maxExp+minExp)/2
tempF		<- -log(1-tempExp)
tempFsrs	<-ExpSeries*tempF
tmpExpSrs	<-1-exp(-tempFsrs)

for (j in 2:InitYear)
{
 for (i in 2:(MaxAge-1))
  {
  # tempNn[j,i]		<-tempNn[j-1,i-1]*(1-tmpExpSrs[j]*vulnN[1,i-1])*exp(-NatMn[1])
  # tempNs[j,i]		<-tempNs[j-1,i-1]*(1-tmpExpSrs[j]*vulnS[1,i-1])*exp(-NatMs[1])

   tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-tempFsrs[j]*vulnN[1,i-1])*exp(-NatMn[1])
   tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-tempFsrs[j]*vulnS[1,i-1])*exp(-NatMs[1])
  }
  # tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*(1-tmpExpSrs[j]*vulnN[1,MaxAge])*exp(-NatMn[1])+ tempNn[j-1,MaxAge]*(1-tmpExpSrs[j]*vulnN[1,MaxAge])*exp(-NatMn[1])
  # tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*(1-tmpExpSrs[j]*vulnS[1,MaxAge])*exp(-NatMs[1])+ tempNs[j-1,MaxAge]*(1-tmpExpSrs[j]*vulnS[1,MaxAge])*exp(-NatMs[1])
   tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-tempFsrs[j]*vulnN[1,MaxAge])*exp(-NatMn[1])+ tempNn[j-1,MaxAge]*exp(-tempFsrs[j]*vulnN[1,MaxAge])*exp(-NatMn[1])
   tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-tempFsrs[j]*vulnS[1,MaxAge])*exp(-NatMs[1])+ tempNs[j-1,MaxAge]*exp(-tempFsrs[j]*vulnS[1,MaxAge])*exp(-NatMs[1])

   EggsN			<-sum(tempNn[j-1,]*matureN[1,]*LenAtAgeN[1,]*WeightAtLenN[1,])
   EggsS			<-sum(tempNs[j-1,]*matureS[1,]*LenAtAgeS[1,]*WeightAtLenS[1,])

   tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1,j],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtLenN[1,],LenAtAgeIN=LenAtAgeN[1,])
   tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1,j],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtLenS[1,],LenAtAgeIN=LenAtAgeS[1,])

   tempCatchAtLenN[j,]	<-(tempFsrs[j]/(tempFsrs[j]+NatMn[1])) * (1-exp(-((tempFsrs[j]*vulnN[1,])+NatMn[1]))) * tempNn[j,]
   tempCatchN[j]		<-sum(tempCatchAtLenN[j,]*LenAtAgeN[1,]*WeightAtLenN[1,])

   tempCatchAtLenS[j,]	<-(tempFsrs[j]/(tempFsrs[j]+NatMs[1])) * (1-exp(-((tempFsrs[j]*vulnS[1,])+NatMs[1]))) * tempNs[j,]
   tempCatchS[j]		<-sum(tempCatchAtLenS[j,]*LenAtAgeS[1,]*WeightAtLenS[1,])

 }

tempDepN	<-sum(tempNn[InitYear,]*matureN[1,]*LenAtAgeN[1,]*WeightAtLenN[1,])
tempDepS	<-sum(tempNs[InitYear,]*matureS[1,]*LenAtAgeS[1,]*WeightAtLenS[1,])
tempDep	<-tempDepN+tempDepS

if(tempDep>targetDep)
 minExp	<-tempExp
if(tempDep<targetDep)
 maxExp	<-tempExp
initExp	<-tempExp
}
plot(ExpSeries*tempExp)
mtext(side=3,"Exploitation history",cex=.7)
#===============================================================
# BEGIN SIMULATION OF ASSESSMENT AND HARVEST
#===============================================================
 if(EstimationPlots==1)
 {
  dev.new()
  par(mfrow=c(1,1))
  }
#==tempNn and tempNs from above are the starting points
projNn	<-array(dim=c(SimYear,MaxAge,Nsim))
projNs	<-array(dim=c(SimYear,MaxAge,Nsim))

projCatchAtLenN	<-array(dim=c(SimYear,MaxAge,Nsim))
projCatchAtLenS	<-array(dim=c(SimYear,MaxAge,Nsim))

projCatchN	<-matrix(nrow=Nsim,ncol=SimYear)
projCatchS	<-matrix(nrow=Nsim,ncol=SimYear)
projSSBn	<-matrix(nrow=Nsim,ncol=SimYear)
projSSBs	<-matrix(nrow=Nsim,ncol=SimYear)
projExpBn	<-matrix(nrow=Nsim,ncol=SimYear)
projExpBs	<-matrix(nrow=Nsim,ncol=SimYear)


#==INSERT EXPLOITABLE BIOMASS HERE AND USE THAT FOR PRODUCTION MODEL??

#==transfer historic time series from above
for(x in 1:Nsim)
{
projNn[1:InitYear,,x]			<-tempNn
projNs[1:InitYear,,x]			<-tempNs
projCatchN[x,1:InitYear]		<-tempCatchN
projCatchS[x,1:InitYear]		<-tempCatchS

projCatchAtLenN[1:InitYear,,x]	<-tempCatchAtLenS	#numbers at age
projCatchAtLenS[1:InitYear,,x]	<-tempCatchAtLenN	#numbers at age
}

for(x in 1:InitYear)
{
 projSSBn[,x]	<-sum(projNn[x,,1]*matureN[1,]*LenAtAgeN[1,]*WeightAtLenN[1,])
 projSSBs[,x]	<-sum(projNs[x,,1]*matureS[1,]*LenAtAgeS[1,]*WeightAtLenS[1,])
 projExpBn[,x]	<-sum(projNn[x,,1]*vulnN[1,]*LenAtAgeN[1,]*WeightAtLenN[1,])
 projExpBs[,x]	<-sum(projNs[x,,1]*vulnS[1,]*LenAtAgeS[1,]*WeightAtLenS[1,])
}

CatchErrorN		<-matrix(rnorm(Nsim*SimYear,0,CatchCVn),nrow=Nsim,ncol=SimYear)
CatchErrorS		<-matrix(rnorm(Nsim*SimYear,0,CatchCVs),nrow=Nsim,ncol=SimYear)
CatchAssessN	<-projCatchN*exp(CatchErrorN)
CatchAssessS	<-projCatchS*exp(CatchErrorS)

IndexErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
IndexErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
IndexAssessN	<-projExpBn*exp(IndexErrorN)
IndexAssessS	<-projExpBs*exp(IndexErrorS)

#===DO LENGTH FREQUENCIES EVENTUALLY==
#CatchLenAsses
LenSampleN	<-200
LenSampleS	<-200

#==save management and assessment quantities
FMSY<-matrix(nrow=Nsim,ncol=SimYear)
BMSY<-matrix(nrow=Nsim,ncol=SimYear)
TAC<-matrix(nrow=Nsim,ncol=SimYear)
CurBio<-matrix(nrow=Nsim,ncol=SimYear)

#==into the simulation===========================================================
for(z in 1:Nsim)
{
 for(y in InitYear:(SimYear-1))
 {
 #==collect data
 CatchData<-CatchAssessN[z,!is.na(CatchAssessN[z,])] + CatchAssessS[z,!is.na(CatchAssessS[z,])]
 IndexData<-IndexAssessN[z,!is.na(IndexAssessN[z,])] + IndexAssessS[z,!is.na(IndexAssessS[z,])] 

 #==assess with production model==================================================
 #==add other assessment methods later
 #==fox, PT, age-structured
 x		<-c(.75*(VirBioN+VirBioS),0.2)	#initial values
 outs		<-nlminb(start=x,objective=ProdMod,CatchData=CatchData,IndexData=IndexData )
 PredBio	<-ProdModPlot(outs$par,CatchData,IndexData,plots=EstimationPlots)
 FMSY[z,y] 	<-abs(outs$par[2])/2
 BMSY[z,y]	<-outs$par[1]/2
 MSY		<-BMSY[z,y]*FMSY[z,y]
 CurBio[z,y]<-PredBio[length(PredBio)]
 TAC[z,y]	<-CurBio[z,y]*FMSY[z,y]

 if(EstimationPlots==1)
 {
  lines(BMSY[z,],col=3,lty=3)
  legend("topright",bty='n',col=c(NA,1,2,3),pch=c(15,NA,NA,NA),lty=c(NA,1,2,3),
          legend=c("Index obs","Index est","Catch","BMSY"))
 }
 #==find the F that removes the TAC given the fishery selectivity===================
 #==this TAC is not necessarily the 'true' TAC--there is error in the biomass and FMSY
 minExp<-0.00001
 maxExp<-0.99999

 for(p in 1:30)
 {
 tempExp	<-(maxExp+minExp)/2
 tempF	<- -log(1-tempExp)
 tempFn	<-tempF*CatchShareN
 tempFs	<-tempF*CatchShareS

   finCatLenFn	<-(tempFn/(tempFn+NatMn[y])) * (1-exp(-((tempFn*vulnN[y,])+NatMn[y]))) * projNn[y,,z]
   finCatFn		<-sum(finCatLenFn*LenAtAgeN[y,]*WeightAtLenN[y,])

   finCatLenFs	<-(tempFs/(tempFs+NatMs[y])) * (1-exp(-((tempFs*vulnS[y,])+NatMs[y]))) * projNn[y,,z]
   finCatFs		<-sum(finCatLenFs*LenAtAgeS[y,]*WeightAtLenS[y,])

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

   EggsN				<-sum(projNn[y-1,,z]*matureN[y,]*LenAtAgeN[y,]*WeightAtLenN[y,])
   EggsS				<-sum(projNs[y-1,,z]*matureS[y,]*LenAtAgeS[y,]*WeightAtLenS[y,])

   projNn[y+1,1,z]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[y],RzeroIN=RzeroN[y],RecErrIN=RecErrN[z,y],recType="BH",NatMin=NatMn[y],
							vulnIN=vulnN[y,],matureIN=matureN[y,],weightIN=WeightAtLenN[y,],LenAtAgeIN=LenAtAgeN[y,])
   projNs[y+1,1,z]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[y],RzeroIN=RzeroS[y],RecErrIN=RecErrS[z,y],recType="BH",NatMin=NatMs[y],
							vulnIN=vulnS[y,],matureIN=matureS[y,],weightIN=WeightAtLenS[y,],LenAtAgeIN=LenAtAgeS[y,])

   projCatchAtLenN[y+1,,z]	<-(ApplyFn/(ApplyFn+NatMn[y])) * (1-exp(-((ApplyFn*vulnN[y,])+NatMn[y]))) * projNn[y+1,,z]
   projCatchN[z,y+1]		<-sum(projCatchAtLenN[y+1,,z]*LenAtAgeN[y,]*WeightAtLenN[y,])

   projCatchAtLenS[y+1,,z]	<-(ApplyFs/(ApplyFs+NatMs[y])) * (1-exp(-((ApplyFs*vulnS[y,])+NatMs[y]))) * projNs[y+1,,z]
   projCatchS[z,y+1]		<-sum(projCatchAtLenS[y+1,,z]*LenAtAgeS[y,]*WeightAtLenS[y,])

   #==UPDATE CATCH DATA AND INDEX DATA
   projSSBn[z,y+1]	<-sum(projNn[y+1,,z]*matureN[y,]*LenAtAgeN[y,]*WeightAtLenN[y,])
   projSSBs[z,y+1]	<-sum(projNs[y+1,,z]*matureS[y,]*LenAtAgeS[y,]*WeightAtLenS[y,])
   projExpBn[z,y+1]	<-sum(projNn[y+1,,z]*vulnN[y,]*LenAtAgeN[y,]*WeightAtLenN[y,])
   projExpBs[z,y+1]	<-sum(projNs[y+1,,z]*vulnS[y,]*LenAtAgeS[y,]*WeightAtLenS[y,])

   CatchAssessN[z,y+1]	<-projCatchN[z,y+1]*exp(CatchErrorN[z,y+1])
   CatchAssessS[z,y+1]	<-projCatchS[z,y+1]*exp(CatchErrorS[z,y+1])

   IndexAssessN[z,y+1]	<-projExpBn[z,y+1]*exp(IndexErrorN[z,y+1])
   IndexAssessS[z,y+1]	<-projExpBs[z,y+1]*exp(IndexErrorS[z,y+1])
 } # end y
} # end z















