###############################################################
SBPRfunc<-function(F,Rzero,NatM,vuln,mature,weight,LenAtAge,MaxAge)
{
  #initial conditions
  N<-matrix(ncol=MaxAge,nrow=50)
  N[1,]<-initialN(Rzero,NatM,MaxAge)
  
  for (j in 2:50)
  {
 for (i in 2:(MaxAge-1))
  N[j,i]<-N[j-1,i-1]*exp(-F*vuln[i-1])*exp(-NatM)

 N[j,MaxAge]<-(N[j-1,MaxAge-1]+N[j-1,MaxAge])*exp(-F*vuln[MaxAge])*exp(-NatM)
 N[j,1]<-Rzero
  }
  SBP<-sum(N[50,]*mature*weight)
  SBPR<-sum(N[50,]*mature*weight)/Rzero
  YPR<-sum(N[50,]*(1-exp(-F*vuln)) *weight)
  list(SBP,SBPR,YPR)
  #return(SBPR)
}

###############################################################

Recruitment<-function(EggsIN,steepnessIN,RzeroIN,RecErrIN,recType,NatMin,vulnIN,matureIN,weightIN,LenAtAgeIN,MaxAge,sigmaRin)
{
  if(recType=="BH")
  {
 Szero<-unlist(SBPRfunc(0,RzeroIN,NatMin,vulnIN,matureIN,weightIN,LenAtAgeIN,MaxAge)[1])	#since there are only a couple possible values, code could be sped up by changing this
 alpha<-(Szero/RzeroIN)*(1-steepnessIN)/(4*steepnessIN)
 beta<-(5*steepnessIN-1)/(4*steepnessIN*RzeroIN)
 Recruit<-(EggsIN/(alpha+beta*EggsIN))*exp(RecErrIN-(sigmaRin^2)/2)
  }
  if(recType=="Rick")
  {
 Szero<-unlist(SBPRfunc(0,RzeroIN,NatMin,vulnIN,matureIN,weightIN,LenAtAgeIN,MaxAge)[1])
 beta<-log(5*steepnessIN)/(0.8*Szero)
 alpha<-log(RzeroIN/Szero)+(beta*Szero)
 Recruit<-(EggsIN*exp(alpha-beta*EggsIN))*exp(RecErrIN-(sigmaRin^2)/2)
  }
  return(Recruit)
}

###############################################################

initialN<-function(Rzero,NatM,inAge)
{
  test		<-rep(0,100)
  test[1]	<-Rzero[1]
  for (i in 2:100)
   {test[i]	<-test[i-1]*exp(-NatM)}
  virInit2	<-c(test[1:(inAge-1)],sum(test[inAge:100]))
  return(virInit2)
}

###############################################################

BevHolt<-function(x)
{
alpha		<-abs(x[1])
beta		<-abs(x[2])
sigma		<-abs(x[3])
Pred		<-alpha*spbio/(1+(spbio/beta))
loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
return(sum(-1*loglike))
}

###############################################################

Ricker<-function(x)
{
alpha		<-x[1]
beta		<-x[2]
sigma		<-abs(x[3])
Pred		<-spbio*exp(alpha-beta*spbio)
loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
return(sum(-1*loglike))
}

###############################################################

StrLine<-function(x)
{
 Pred		<-x[1]
 sigma	<-x[2]
 loglike	<- log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
 return(sum(-1*loglike))
}

###############################################################

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

###############################################################

calcBHrec<-function(SRouts,spbioIn)
{
 alpha	<-abs(SRouts[[2]][1])
 beta		<-abs(SRouts[[2]][2])
 Pred		<-alpha*spbioIn/(1+(alpha*spbioIn/beta)) 
 return(Pred)
}

###############################################################

calcRICKrec<-function(SRouts,spbioIn)
{
 alpha	<-SRouts[[2]][1]
 beta		<-SRouts[[2]][2]
 Pred		<-spbioIn*exp(alpha-beta*spbioIn)
 return(Pred)
}

###############################################################

calcREF<-function(SRouts)
{
if(SRouts[[1]][1]!="Average")
{
 testExp	<-seq(0.001,1,0.01)
 recYield	<-rep(0,length(testExp))
 projYrs	<-100
 tempN	<-matrix(nrow=projYrs,ncol=maxAge)
 tempN[1,]	<-virInit
 tempYield	<-rep(0,length(testExp))
for(f in 1:length(testF))
{
 tempExp<-testExp[f]
 for(j in 2:projYrs)
 {
  for (i in 2:(maxAge-1))
   tempN[j,i]		<-tempN[j-1,i-1]*(1-tempExp*vuln[i-1])*survival
   tempN[j,maxAge]	<-(tempN[j-1,(maxAge-1)])*(1-tempExp*vuln[maxAge])*survival + tempN[j-1,maxAge]*(1-tempExp*vuln[maxAge])*survival
   spbioIn			<-sum(tempN[j-1,]*mature*weight)
   if(SRouts[[1]][1]=="BevH")
    tempN[j,1]		<-calcBHrec(SRouts,spbioIn) 
   if(SRouts[[1]][1]=="Rick")
    tempN[j,1]		<-calcRICKrec(SRouts,spbioIn) 
  }
   tempYield[f]		<-sum(tempN[j,]*tempExp*vuln*weight)
 }
}
calcRICKrec(SRouts,virBio)
if(SRouts[[1]][1]=="Average")
{
#STUFF
}

}

###############################################################
#==produces population and yield from a given fishing mortality
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

##############################################################
#==calculates Yield at F, to be passed to optimizing function
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
return(-1*(tempCatchS[j]+tempCatchN[j]))
}
#================================================================================
#==calculates F35%
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
Bzero			<-ProjPopDym(fmortN=0,MaxAge=MaxAge,vulnN=vulnN,
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

##############################################################
ProdMod<-function(x,CatchData,IndexData)
{
K<-abs(x[1])
r<-abs(x[2])
predBio<-rep(0,length(IndexData))
predBio[1]<-K
for(i in 2:length(CatchData))
{
 predBio[i]<-predBio[i-1]+r*predBio[i-1]*(1-predBio[i-1]/K)-CatchData[i]
}
SSQ<-sum((predBio-IndexData)^2,na.rm=T)
return(SSQ)
}

#################################################################
ProdModPlot<-function(x,CatchData,IndexData,plots=0)
{
K<-abs(x[1])
r<-abs(x[2])
predBio<-rep(0,length(IndexData)+1)
predBio[1]<-K
CatchData<-c(CatchData,0)
for(i in 2:length(CatchData))
{
 predBio[i]<-predBio[i-1]+r*predBio[i-1]*(1-predBio[i-1]/K)-CatchData[i]
}
if(plots==1)
{
plot(IndexData[1:(length(IndexData)-1)],ylim=c(0,max(IndexData,na.rm=T)),las=1,ylab="")
lines(predBio[1:(length(predBio)-1)])
lines(CatchData[1:(length(CatchData)-1)],lty=2,col=2)
}
return(predBio)
}

#################################################################
AgeAssPlot<-function(REP,CTL,DAT,TRU,data2,MSEdir,CTLName)
{
temp<-grep("survey years",DAT)[1]
yearsDat<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
temp<-grep("number of ages",DAT)[1]
maxAge<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
AgeBin<-seq(1,maxAge)
temp<-grep("Length Bin",DAT)[1]
LengthBin<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))

temp<-grep("recruitment",TRU)
trueRec<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("Catch",TRU)
trueCatch<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("survey selectivity",TRU)
trueSurvSelPar<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("fishery selectivity",TRU)
trueFishSelPar<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("fishing mortality",TRU)
trueFmort<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("spawning biomass",TRU)
trueSpbio<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))

temp<-grep("spawning biomass",REP)
predSpBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("B35",REP)
B35<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))


temp<-grep("pred survey bio",REP)
predSurvBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs survey bio",REP)
obsSurvBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred cpue bio",REP)
predCpueBio <-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs cpue bio",REP)
cpueIndex <-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred catch bio",REP)
predCatchBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs catch bio",REP)
obsCatchBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("obs surv len freq",REP)
obsSurlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]
temp<-grep("pred surv len freq",REP)
predSurlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]

temp<-grep("obs catch len freq",REP)
obsCatchlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]
temp<-grep("pred catch len freq",REP)
predCatchlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]

#==================================================
#   Par file plots
#================================================
cnt<-grep("mean_log_rec",data2)
meanREC<-as.numeric(data2[cnt+1])
cnt<-grep("rec_dev",data2)
recDevs<-as.numeric(data2[(cnt+1):(cnt+yearsDat)])
TotRec<-exp(meanREC+recDevs)

fmortYrs<-seq(1,yearsDat)
cnt<-grep("log_avg_fmort_dir",data2)
logFdir<-as.numeric(data2[cnt+1])
cnt<-grep("fmort_dir_dev",data2)
fmort_dir_dev<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])
TotFdir<-exp(logFdir+fmort_dir_dev)

cnt<-grep("stNatLen",data2)
stNatLen<-as.numeric(data2[(cnt+1):(cnt+maxAge)])

png(file.path(MSEdir,"plots",paste0("PopulationProcessesEst_",CTLName,".png")),res=1200,width=6,height=6,units="in")
par(mfrow=c(3,1),mar=c(1,4,1,1),oma=c(3,1,1,1))
plot(TotRec[1:(yearsDat-1)],type="l",las=1,ylab="Recruitment")
lines(trueRec[1:(yearsDat-1)],lty=2,col=2)
plot(TotFdir,type="l",las=1,ylab="Directed F")
lines(trueFmort,lty=2,col=2)
barplot(exp(stNatLen))
#==plot estimated selectivities

cnt<-grep("srv_sel50",data2)
sel50surv<-as.numeric(data2[cnt+1])
cnt<-grep("srv_sel95",data2)
sel95surv<-as.numeric(data2[cnt+1])

cnt<-grep("log_avg_SelPars50",data2)
sel50fish<-exp(as.numeric(data2[cnt+1]))
cnt<-grep("log_avg_SelPars95",data2)
sel95fish<-exp(as.numeric(data2[cnt+1]))

cnt<-grep("SelPars_dev50",data2)
sel50_dir_dev<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])
cnt<-grep("SelPars_dev95",data2)
sel95_dir_dev<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])

SurvSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50surv)/(sel95surv-sel50surv)))
FishSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50fish)/(sel95fish-sel50fish)))

trueSurvSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueSurvSelPar[1])/(trueSurvSelPar[2]-trueSurvSelPar[1])))
#plot(trueSurvSel~LengthBin,type="l",xlab="Age",ylab="Selectivity")
trueFishSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueFishSelPar[1])/(trueFishSelPar[2]-trueFishSelPar[1])))
#plot(FishSel~AgeBin,type="l",xlab="Age",ylab="Selectivity")

#lines(FishSel~AgeBin,lty=2,col=2)
#lines(SurvSel~AgeBin,lty=2,col=2)
#lines(trueFishSel~LengthBin,lty=2,col=1)
#legend("bottomright",bty='n',col=c(1,1,2,2),lty=c(1,2,1,2),legend=c("True Survey","Est Survey","True Fishery","Est Fishery"))
dev.off()

png(file.path(MSEdir,"plots", paste0("FitsToDataAgeStruc_",CTLName,".png")),res=1200,width=6,height=6,units="in")
par(mfrow=c(2,2),mar=c(.1,.1,.1,.1),oma=c(5,5,1,5))
#plot(obsSurvBio,yaxt='n',ylim=c(0,max(obsSurvBio,predSurvBio,trueSpbio,na.rm=T)))
#axis(side=2,las=1)
#lines(predSurvBio/2)
#lines(trueSpbio,col=2,lty=2)
#lines(predSpBio,col=1,lty=1)
#legend("topright",bty='n',"Survey")

#plot(obsSurvBio,yaxt='n',ylim=c(0,max(obsSurvBio,predSurvBio,trueSpbio,na.rm=T)))
plot(obsSurvBio,yaxt='n',xaxt='n')
axis(side=2,las=1)
lines(predSurvBio)
legend("topright",bty='n',"Survey")

#plot(cpueIndex,yaxt='n',ylim=c(0,max(predCpueBio ,cpueIndex)))
plot(cpueIndex,yaxt='n',xaxt='n')
axis(side=4,las=1)
lines(predCpueBio  )
#lines(trueCatch,col=2,lty=2)
legend("topright",bty='n',"CPUE")

plot(obsCatchBio ,yaxt='n')
lines(predCatchBio,yaxt='n')
axis(side=2,las=1)
legend("topright",bty='n',"Catch")

plot(predSpBio,yaxt='n',type="l",ylim=c(0,max(predSpBio,na.rmm=T)))
lines(trueSpbio,col=2,lty=2)
abline(h=B35,lty=3,col=3)
axis(side=4,las=1)
legend("topright",bty='n',c("Spawning biomass","True","Estimated"),lty=c(NA,1,2))
dev.off()


png(file.path(MSEdir,"plots",paste0("FitsToSurvLengths_",CTLName,".png")),res=1200,width=6,height=6,units="in")
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsSurlen[i,],axes=F)
lines(predSurlen[i,])
}
plot.new()
legend("center","Survey")
dev.off()

png(file.path(MSEdir,"plots",paste0("FitsToCatchLengths_",CTLName,".png")),res=1200,width=6,height=6,units="in")
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsCatchlen[i,],axes=F)
lines(predCatchlen[i,])
}
plot.new()
legend("center","Catch")
dev.off()
}
##################################################

TakeOut<-function(SearchTerm,SearchPool)
{
 TakeIndex<-grep(SearchTerm,SearchPool[,3])[1]
 temp<-suppressWarnings(as.numeric(SearchPool[TakeIndex,1]))
 if(is.na(temp))
  temp<-eval(parse(text=SearchPool[TakeIndex,1]))
 return(temp)
}
#input<-"C:/Users/Cody/Desktop/Publications/BlueShark/CTLfile.csv"

#########################################################
ReadCTLfile<-function(input)
{
 OM				<-NULL
 SearchPool			<-as.matrix(read.csv(input))
 OM$Nsim			<-TakeOut("Nsim",SearchPool)
 OM$SimYear			<-TakeOut("SimYear",SearchPool)
 OM$InitYear		<-TakeOut("InitYear",SearchPool)
 OM$AssessmentType	<-TakeOut("AssessmentType",SearchPool)
 OM$FisheryIndepenDat	<-TakeOut("FisheryIndepenDat",SearchPool)
 OM$LifeHistoryPlots	<-TakeOut("LifeHistoryPlots",SearchPool)
 OM$EstimationPlots	<-TakeOut("EstimationPlots",SearchPool)
 OM$AssessmentData	<-TakeOut("AssessmentData",SearchPool)
 OM$PlotYieldCurve	<-TakeOut("PlotYieldCurve",SearchPool)

 OM$CatchCVn		<-TakeOut("CatchCVn",SearchPool)
 OM$CatchCVs		<-TakeOut("CatchCVs",SearchPool)
 OM$CPUEcvN			<-TakeOut("CPUEcvN",SearchPool)
 OM$CPUEcvS			<-TakeOut("CPUEcvS",SearchPool)
 OM$IndexCVn		<-TakeOut("IndexCVn",SearchPool)
 OM$IndexCVs		<-TakeOut("IndexCVs",SearchPool)
 OM$LenSampleN		<-TakeOut("LenSampleN",SearchPool)
 OM$LenSampleS		<-TakeOut("LenSampleS",SearchPool)
 OM$GrowthSDn		<-TakeOut("GrowthSDn",SearchPool)
 OM$GrowthSDs		<-TakeOut("GrowthSDs",SearchPool)
 OM$LengthBinN		<-TakeOut("LengthBinN",SearchPool)

 OM$MaxAge			<-TakeOut("MaxAge",SearchPool)
 OM$NatMn			<-TakeOut("NatMn",SearchPool)
 OM$VonKn			<-TakeOut("VonKn",SearchPool)
 OM$LinfN			<-TakeOut("LinfN",SearchPool)
 OM$t0n			<-TakeOut("t0n",SearchPool)
 OM$mat50n			<-TakeOut("mat50n",SearchPool)
 OM$mat95n			<-TakeOut("mat95n",SearchPool)
 OM$alphaN			<-TakeOut("alphaN",SearchPool)
 OM$betaN			<-TakeOut("betaN",SearchPool)
 OM$MaxMovingN		<-TakeOut("MaxMovingN",SearchPool)
 OM$Move50n			<-TakeOut("Move50n",SearchPool)
 OM$Move95n			<-TakeOut("Move95n",SearchPool)
 OM$steepnessN		<-TakeOut("steepnessN",SearchPool)
 OM$RzeroN			<-TakeOut("RzeroN",SearchPool)
 OM$sigmaRn			<-TakeOut("sigmaRn",SearchPool)

 OM$sel50n			<-TakeOut("sel50n",SearchPool)
 OM$sel95n			<-TakeOut("sel95n",SearchPool)
 OM$surv50n			<-TakeOut("surv50n",SearchPool)
 OM$surv95n			<-TakeOut("surv95n",SearchPool)

 OM$NatMs			<-TakeOut("NatMs",SearchPool)
 OM$VonKs			<-TakeOut("VonKs",SearchPool)
 OM$LinfS			<-TakeOut("LinfS",SearchPool)
 OM$t0s			<-TakeOut("t0s",SearchPool)
 OM$mat50s			<-TakeOut("mat50s",SearchPool)
 OM$mat95s			<-TakeOut("mat95s",SearchPool)
 OM$alphaS			<-TakeOut("alphaS",SearchPool)
 OM$betaS			<-TakeOut("betaS",SearchPool)
 OM$MaxMovingS		<-TakeOut("MaxMovingS",SearchPool)
 OM$Move50s			<-TakeOut("Move50s",SearchPool)
 OM$Move95s			<-TakeOut("Move95s",SearchPool)
 OM$steepnessS		<-TakeOut("steepnessS",SearchPool)
 OM$RzeroS			<-TakeOut("RzeroS",SearchPool)
 OM$sigmaRs			<-TakeOut("sigmaRs",SearchPool)

 OM$sel50s			<-TakeOut("sel50s",SearchPool)
 OM$sel95s			<-TakeOut("sel95s",SearchPool)
 OM$surv50s			<-TakeOut("surv50s",SearchPool)
 OM$surv95s			<-TakeOut("surv95s",SearchPool)

 OM$HistoricalFn		<-TakeOut(SearchTerm="HistoricalFn",SearchPool=SearchPool)
 OM$HistoricalFs		<-TakeOut(SearchTerm="HistoricalFs",SearchPool=SearchPool)
 OM$PastFsdN		<-TakeOut("PastFsdN",SearchPool)
 OM$PastFsdS		<-TakeOut("PastFsdS",SearchPool)

 OM$HarvestControlN	<-TakeOut("HarvestControlN",SearchPool)
 OM$ConstantCatchN	<-TakeOut("ConstantCatchN",SearchPool)
 OM$ConstantFn		<-TakeOut("ConstantFn",SearchPool)
 OM$HCalphaN		<-TakeOut("HCalphaN",SearchPool)
 OM$HCbetaN			<-TakeOut("HCbetaN",SearchPool)

 OM$HarvestControlS	<-TakeOut("HarvestControlS",SearchPool)
 OM$ConstantCatchS	<-TakeOut("ConstantCatchS",SearchPool)
 OM$ConstantFs		<-TakeOut("ConstantFs",SearchPool)
 OM$HCalphaS		<-TakeOut("HCalphaS",SearchPool)
 OM$HCbetaS			<-TakeOut("HCbetaS",SearchPool)

 OM$SmalNum			<-TakeOut("SmallNum",SearchPool)
 OM$InitSmooth		<-TakeOut("InitSmooth",SearchPool)
 OM$FmortPen		<-TakeOut("FmortPen",SearchPool)
 OM$RecruitPen		<-TakeOut("RecruitPen",SearchPool)
 OM$Mpenalty		<-TakeOut("Mpenalty",SearchPool)
 OM$EstM			<-TakeOut("EstM",SearchPool)
 OM$TimeVaryM		<-TakeOut("TimeVaryM",SearchPool)

 OM$EstGrowthK		<-TakeOut("EstGrowthK",SearchPool)
 OM$TimeVaryGrowthK	<-TakeOut("TimeVaryGrowthK",SearchPool)
 OM$EstLinf			<-TakeOut("EstLinf",SearchPool)
 OM$TimeVaryLinf		<-TakeOut("TimeVaryLinf",SearchPool)
 OM$Growthpenalty		<-TakeOut("Growthpenalty",SearchPool)
 OM$TimeVarySel50		<-TakeOut("TimeVarySel50",SearchPool)
 OM$TimeVarySel95		<-TakeOut("TimeVarySel95",SearchPool)
 OM$SelPenalty		<-TakeOut("SelPenalty",SearchPool)
 OM$ProjectTimeVary	<-TakeOut("ProjectTimeVary",SearchPool)
 OM$InitValSig		<-TakeOut("InitValSig",SearchPool)

 OM$InitBzeroMod		<-TakeOut("InitBzeroMod",SearchPool)
 OM$InitGrowthRate	<-TakeOut("InitGrowthRate",SearchPool)

 OM$TwoPop			<-TakeOut("TwoPop",SearchPool)

 list(OM=OM)
}
#====================================================
CleanInput<-function(input,SimYear)
{
if(length(input)==1)
 output		<-rep(input,SimYear)
if(length(input)>1)
 output		<-input
 return(output)
}
#=============================================
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
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,],MaxAge=MaxAge,sigmaRin=sigmaRn[1])
   }
plot(record~ssb,type='l')
mtext(side=3,"SSB vs Rec (N)",cex=.7)
ssb<-seq(1,VirBioS,VirBioS/100)
record<-ssb
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,],MaxAge=MaxAge,sigmaRin=sigmaRs[1])
   }
plot(record~ssb,type='l')
mtext(side=3,"SSB vs Rec (S)",cex=.7)


plot(HistoricalFn,type='l')
lines(HistoricalFs,lty=2,col=2)
mtext(side=3,"Historical F",cex=.7)
 }

#=======================================================

HarvestControlRule<-function(FMSY,BMSY,ExploitBio,SpawnBio,alpha,beta,HarvestControl,ConstantCatch=0,ConstantF=0)
{
 if(HarvestControl==1)
  outTAC<-ExploitBio*(FMSY)
 if(HarvestControl==2)
  outTAC<-ConstantCatch
 if(HarvestControl==3)
  outTAC<-ExploitBio*(ConstantF)

 if(HarvestControl==4)
  {
   useF	<-FMSY
   if(SpawnBio<BMSY)
    {
    useF	<-0
    if(SpawnBio>(BMSY*beta)) 
      useF	<-FMSY-((FMSY/(BMSY-beta*BMSY))*(BMSY-SpawnBio))
     }
   outTAC<-ExploitBio*useF
  }
 return(outTAC)
}
#####################################################################
   
PolygonPlots<-function(Truth=NA,Estimated=NA,Observed=NA,AddLegend=F,bottom=F,quantity=NA,SimYear=NA,Nsim=NA)
{
 EstShape<-apply(Estimated[,1:(ncol(Estimated))],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)
 TrueShape<-apply(Truth[,1:(ncol(Truth))],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)

 DataColor<-"#6699ff11"
 EstimateColor25<-"darkgrey"
 EstimateColor5<-"lightgrey"
 TruthColor<-"#FF000033"
  
  if(bottom==F)
   plot(-100000000000,ylim=c(0,max(EstShape,TrueShape,Observed,na.rm=t)),xlim=c(1,SimYear),las=1,ylab="",xlab="Year",xaxt='n')
  if(bottom==T)
   plot(-100000000000,ylim=c(0,max(EstShape,TrueShape,Observed,na.rm=t)),xlim=c(1,SimYear),las=1,ylab="",xlab="Year")
  
  if(!is.na(quantity))
 legend("bottomleft",bty='n',legend=quantity)

 if(length(Observed)>1)
 {
 if(AddLegend==T)
  legend("topright",bty='n',pch=c(15,16,NA),lty=c(NA,NA,1),col=c(EstimateColor5,DataColor,TruthColor),legend=c("Estimates","Observed Data","Truth"))
 for(z in 1:Nsim)
  points(Observed[z,]~seq(1,(SimYear-1)),col=DataColor,pch=16)
 }
 polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(TrueShape[1,1:SimYear-1],rev(TrueShape[4,1:SimYear-1])),col=EstimateColor5,border=F)
 polygon(x=c(seq(1,(SimYear-1)),rev(seq(1,(SimYear-1)))),y=c(TrueShape[2,1:SimYear-1],rev(TrueShape[3,1:SimYear-1])),col=EstimateColor25,border=F)
 if(AddLegend==T&length(Observed)<1)
 legend("topright",bty='n',pch=c(15,NA),lty=c(NA,1),col=c(EstimateColor5,TruthColor),legend=c("Truth","Estimates"))

 for(z in 1:Nsim)
  lines(Estimated[z,]~seq(1,(SimYear)),col=TruthColor,lty=1)
}

#####################################################################
selAtAgeFunc<-function(inputLen,K,Linf,t0)
 {
  selAtAge<-((log(1-inputLen/Linf)/-K))+t0
 return(selAtAge)
 }
#####################################################################
CheckRetro<-function(RetroPeels=6,DrawDir,PlotRetro=0,out,MSEdir)
{
  Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
  SimYear		<-out$OM$SimYear			# total number of years in simulation
  InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
  
  predSpBioRetro	<-array(dim=c(RetroPeels,SimYear,Nsim))
  trueSpBioRetro    <-matrix(ncol=SimYear,nrow=Nsim)
  
  for(n in 1:Nsim)
    for(v in 0:(RetroPeels-1))
    {
    IndSimFolder<-file.path(DrawDir,n,SimYear-v)
    REP<-readLines(file.path(MSEdir,IndSimFolder,"simass.rep"))
    DAT<-readLines(file.path(MSEdir,IndSimFolder,"simass.DAT"))
    
    temp<-grep("survey years",DAT)[1]
    yearsDat<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
    temp<-grep("spawning biomass",REP)
    predSpBioRetro[(v+1),(1:(SimYear-v-1)),n]<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
  
  }
  
  for(n in 1:Nsim)
  {
    IndSimFolder<-file.path(DrawDir,n,SimYear)
    TRU<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))
    temp<-grep("spawning biomass",TRU)
    trueSpBioRetro[n,]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
  }
  
  plotSegment<-predSpBioRetro[,((SimYear-RetroPeels-4):SimYear),1]
  if(PlotRetro==1)
  {
    dev.new()
    par(mar=c(4,5,1,1))
    plot(y=plotSegment[1,],x=((SimYear-RetroPeels-4):SimYear),type="l",ylim=c(min(plotSegment,na.rm=T),max(plotSegment,na.rm=T)),lty=2,las=1,
     ylab="",xlab="Year")
    for(x in 2:RetroPeels)
     lines(y=plotSegment[x,],x=((SimYear-RetroPeels-4):SimYear),lty=x+1)
    
     lines(y=trueSpBioRetro[1,((SimYear-RetroPeels-4):SimYear)],x=((SimYear-RetroPeels-4):SimYear),col=2)
    
    legend("topleft",bty='n',lty=1,col=2,"Truth")
  }
  list(predSpBioRetro,trueSpBioRetro)
}

#################################################
CheckGradient<-function(DrawDir,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)

Gradients	<-matrix(ncol=Nsim,nrow=SimYear)

for(n in 1:Nsim)
for(v in (InitYear+1):SimYear)
{
IndSimFolder<-file.path(DrawDir,n,v)
PAR<-readLines(file.path(MSEdir,IndSimFolder,"simass.par"))

temp<-grep("gradient",PAR)[1]
GradientLine<-(unlist(strsplit(PAR[temp],split=" ")))
Gradients[v,n]<-as.numeric(GradientLine[length(GradientLine)])
}
return(Gradients)
}
#################################################
CheckReferencePoints<-function(DrawDir,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)

B35	<-matrix(ncol=Nsim,nrow=SimYear)
F35	<-matrix(ncol=Nsim,nrow=SimYear)
OFL	<-matrix(ncol=Nsim,nrow=SimYear) 
tB35	<-matrix(ncol=Nsim,nrow=SimYear)
tOFL	<-matrix(ncol=Nsim,nrow=SimYear)

for(n in 1:Nsim)
for(v in (InitYear+1):SimYear)
{
IndSimFolder<-file.path(DrawDir,n,v)
REP<-readLines(file.path(MSEdir,IndSimFolder,"simass.REP"))
TRU<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))

temp<-grep("B35",REP)[1]
B35[v,n]<-as.numeric((unlist(strsplit(REP[temp+1],split=" "))))

temp<-grep("F35",REP)[2]
F35[v,n]<-as.numeric((unlist(strsplit(REP[temp+1],split=" "))))

temp<-grep("OFL",REP)[2]
OFL[v,n]<-as.numeric((unlist(strsplit(REP[temp+1],split=" "))))

if(v==SimYear)
{
temp<-grep("B35",TRU)
tB35[,n] <-as.numeric((unlist(strsplit(TRU[temp+n],split=" "))))

temp<-grep("OFL",TRU)
tOFL[,n]<-as.numeric((unlist(strsplit(TRU[temp+n],split=" ")))) 
}
temp<-grep("F35",TRU)
tF35 <-as.numeric((unlist(strsplit(TRU[temp+1],split=" ")))) 


}
list(B35,F35,OFL,tB35,tF35,tOFL)
}
###############################################################
CalculateBias<-function(RetroOuts)
{
 preds<-RetroOuts[[1]]
 BiasSSB<-matrix(ncol=dim(preds)[3],nrow=dim(preds)[1])
 for(x in 1:dim(preds)[3])
  {
   tempPreds<-preds[,,x]
   tempTrue<-RetroOuts[[2]][x,]
  for(y in 1:dim(preds)[1])
   BiasSSB[y,x]<-(tempPreds[y,(ncol(tempPreds)-y)]-tempTrue[(ncol(tempPreds)-y)])/tempTrue[(ncol(tempPreds)-y)]
   }
 return(BiasSSB)
}

###############################################################
CalculateMohn<-function(RetroOuts)
{
 preds<-RetroOuts[[1]]
 MohnsRho<-matrix(ncol=dim(preds)[3],nrow=(dim(preds)[1]-1))
 for(x in 1:dim(preds)[3])
  {
   tempPreds<-preds[,,x]
   tempTrue<-preds[1,,x]
  for(y in 2:dim(preds)[1])
   MohnsRho[y-1,x]<-(tempPreds[y,(ncol(tempPreds)-y)]-tempTrue[(ncol(tempPreds)-y)])/tempTrue[(ncol(tempPreds)-y)]
   }
 return(MohnsRho)
}


########################################################
GeMSplots<-function(SourceFolder,out)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
AssessmentType	<-out$OM$AssessmentType		# this can be "Production" or "Non-spatial age" ("Spatial age" and "Time-varying production" in the works)
MaxAge		<-out$OM$MaxAge

if(AssessmentType==1)
{
 PROD<-readLines(file.path(SourceFolder,"ProdOuts.csv"))

temp<-grep("true Catch",PROD)[1]
TrueCatchOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true fishing mortality",PROD)[1]
TrueFmortOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true spawning biomass",PROD)[1]
TrueSpBioOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true CPUE",PROD)[1]
TrueCPUEOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true total allowable catch",PROD)[1]
TrueTACOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

temp<-grep("true BMSY",PROD)[1]
TrueBMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true FMSY",PROD)[1]
TrueFMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

temp<-grep("est cpue ind",PROD)[1]
EstCPUEOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est spawning biomass by year",PROD)[1]
EstCurbioOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est total allowable catch",PROD)[1]
EstTACOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est BMSY",PROD)[1]
EstBMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est FMSY",PROD)[1]
EstFMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

temp<-grep("data cpue",PROD)[1]
DateCPUEout<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

pdf(file.path(SourceFolder,"RelativeErrors.pdf"),height=3,width=4)
PolygonPlots(Truth=TrueCPUEOut,Estimated=EstCPUEOut,Observed=DateCPUEout)
dev.off()

#==plots
REtac<-(EstTACOut-TrueTACOut)/TrueTACOut
REbmsy<-(EstBMSYOut-TrueBMSYOut[1])/TrueBMSYOut[1]
REfmsy<-(EstFMSYOut-TrueFMSYOut[1])/TrueFMSYOut[1]
REbio<-(EstCurbioOut-TrueSpBioOut)/TrueSpBioOut
LostYield<-TrueTACOut-EstTACOut
PercLostYield<-LostYield/TrueTACOut
PercOverfishedTrue<-TrueSpBioOut/TrueBMSYOut[1]
PercOverfishedEst<-EstCurbioOut/TrueBMSYOut[1]

MyBoxPlot<-function(Inputs,Title,bottomplot=0,refVal=0)
{
 maxVal<-max(Inputs,na.rm=T)
 minVal<-min(Inputs,na.rm=T)

 ylimTop<-maxVal
 if(ylimTop<refVal)
  ylimTop<-refVal

 ylimBot<-minVal
 if(ylimBot>refVal)
  ylimBot<-refVal

 coltemp<-apply(Inputs,2,median,na.rm=T)
 colIn<-coltemp
 colIn[coltemp>refVal]<-"green"
 colIn[coltemp<=refVal]<-'red'

 if(bottomplot==0)
 boxplot(Inputs,ylim=c(ylimBot,ylimTop),las=1,xaxt='n',col=colIn)
 if(bottomplot==1)
 boxplot(Inputs,ylim=c(ylimBot,ylimTop),las=1,col=colIn)
 abline(h=refVal,lty=2)
 legend("topleft",bty='n',legend=Title)
}

pdf(file.path(SourceFolder,"RelativeErrors.pdf"),height=8,width=4)
par(mfrow=c(4,1),mar=c(0.1,4,.1,1),oma=c(3,0,3,0))
MyBoxPlot(REtac[,InitYear:SimYear],"Total allowable catch")
MyBoxPlot(REbmsy[,InitYear:SimYear],"BMSY")
MyBoxPlot(REfmsy[,InitYear:SimYear],"FMSY")
MyBoxPlot(REbio[,InitYear:SimYear],"Estimated biomass",bottomplot=1)
mtext(outer=T,side=3,"Relative error")

pdf(file.path(SourceFolder,"YieldStats.pdf"),height=8,width=4)
par(mfrow=c(4,1),mar=c(0.1,4,.1,1),oma=c(3,0,0,0))
MyBoxPlot(LostYield[,InitYear:SimYear],"Lost yield")
MyBoxPlot(PercLostYield[,InitYear:SimYear],"Lost yield/yield")
MyBoxPlot(PercOverfishedTrue[,InitYear:SimYear],"True B/BMSY",refVal=1)
MyBoxPlot(PercOverfishedEst[,InitYear:SimYear],"Estimated B/BMSY",bottomplot=1,refVal=1)

}

if(AssessmentType==2)
{


#==predicted derived quantities storage
predSpBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
predSurvBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
obsSurvBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
predCpueBio 	<-matrix(ncol=SimYear-1,nrow=Nsim)
obsCpueIndex 	<-matrix(ncol=SimYear-1,nrow=Nsim)
predCatchBio	<-matrix(ncol=SimYear-1,nrow=Nsim)
obsCatchBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
obsSurlen		<-array(dim=c(SimYear-1,MaxAge,Nsim))
predSurlen		<-array(dim=c(SimYear-1,MaxAge,Nsim))
obsCatchlen		<-array(dim=c(SimYear-1,MaxAge,Nsim))
predCatchlen	<-array(dim=c(SimYear-1,MaxAge,Nsim))
trueRecPlot		<-matrix(ncol=SimYear-1,nrow=Nsim)
trueCatchPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueSpbioPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueCPUEindPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueSurvIndPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueFmortPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)

trueB35Plot		<-rep(0,Nsim)
trueBMSYPlot	<-rep(0,Nsim)
trueF35Plot		<-rep(0,Nsim)
trueFMSYPlot	<-rep(0,Nsim)

predB35Plot		<-matrix(ncol=SimYear-1,nrow=Nsim)
predBMSYPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
predF35Plot		<-matrix(ncol=SimYear-1,nrow=Nsim)
predFMSYPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)

predB35PlotEnd	<-rep(0,Nsim)
predBMSYPlotEnd	<-rep(0,Nsim)
predF35PlotEnd	<-rep(0,Nsim)
predFMSYPlotEnd	<-rep(0,Nsim)


#===parameters===
meanREC		<-rep(0,Nsim)
recDevs		<-matrix(ncol=SimYear-1,nrow=Nsim)
TotRec		<-matrix(ncol=SimYear-1,nrow=Nsim)
logFdir		<-rep(0,Nsim)
fmort_dir_dev	<-matrix(ncol=SimYear-1,nrow=Nsim)
TotFdir		<-matrix(ncol=SimYear-1,nrow=Nsim)
sel50surv		<-rep(0,Nsim)
sel95surv		<-rep(0,Nsim)
sel50fish		<-rep(0,Nsim)
sel95fish		<-rep(0,Nsim)

OutFolder<-file.path(MSEdir,SourceFolder)
pdf(file.path(OutFolder,"AgeStrucurePlots.pdf"))
#==cycle through final .rep files, pull out data
for(n in 1:Nsim)
{
IndSimFolder<-file.path(SourceFolder,n,SimYear)

REP<-readLines(file.path(MSEdir,IndSimFolder,"simass.rep"))
CTL<-readLines(file.path(MSEdir,IndSimFolder,"simass.CTL"))
DAT<-readLines(file.path(MSEdir,IndSimFolder,"simass.DAT"))
TRU<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))
data2<-scan(file.path(MSEdir,IndSimFolder,"simass.PAR"),what="character")

if(n==1)
{
temp<-grep("survey years",DAT)[1]
yearsDat<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
temp<-grep("number of ages",DAT)[1]
maxAge<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
AgeBin<-seq(1,maxAge)
temp<-grep("Length Bin",DAT)[1]
LengthBin<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
}

temp<-grep("recruitment",TRU)
trueRecPlot[n,]		<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("Catch",TRU)
trueCatchPlot[n,]		<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("fishing mortality",TRU)
trueFmortPlot[n,]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("spawning biomass",TRU)
trueSpbioPlot[n,]		<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("cpue index",TRU)
trueCPUEindPlot[n,]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("survey index",TRU)
trueSurvIndPlot[n,]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]

#temp<-grep("B35",TRU)
#trueB35Plot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
#temp<-grep("F35",TRU)
#trueF35Plot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("BMSY",TRU)
trueBMSYPlot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("FMSY",TRU)
trueFMSYPlot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))

temp<-grep("spawning biomass",REP)
predSpBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("B35",REP)
predB35PlotEnd[n]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("F35",REP)[2]
predF35PlotEnd[n]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

temp<-grep("pred survey bio",REP)
predSurvBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs survey bio",REP)
obsSurvBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred cpue bio",REP)
predCpueBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs cpue bio",REP)
obsCpueIndex[n,] 	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred catch bio",REP)
predCatchBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs catch bio",REP)
obsCatchBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("obs surv len freq",REP)
obsSurlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]
temp<-grep("pred surv len freq",REP)
predSurlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]

temp<-grep("obs catch len freq",REP)
obsCatchlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]
temp<-grep("pred catch len freq",REP)
predCatchlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]

#==================================================
#   Par file plots
#================================================
cnt<-grep("mean_log_rec",data2)
meanREC[n]<-as.numeric(data2[cnt+1])
cnt<-grep("rec_dev",data2)
recDevs[n,]<-as.numeric(data2[(cnt+1):(cnt+yearsDat)])
TotRec[n,]<-exp(meanREC[n]+recDevs[n,])

fmortYrs<-seq(1,yearsDat)
cnt<-grep("log_avg_fmort_dir",data2)
logFdir[n]<-as.numeric(data2[cnt+1])
cnt<-grep("fmort_dir_dev",data2)
fmort_dir_dev[n,]<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])
TotFdir[n,]<-exp(logFdir[n]+fmort_dir_dev[n,])
}

#==plot the fitted data
par(mfrow=c(3,1),mar=c(.1,3,.1,.1),oma=c(2,3,0,0))
PolygonPlots(Truth=trueCatchPlot,Estimated=predCatchBio,Observed=obsCatchBio,AddLegend=T,quantity="Catch")
PolygonPlots(Truth=trueCPUEindPlot,Estimated=predCpueBio,Observed=obsCpueIndex,quantity="CPUE")
PolygonPlots(Truth=trueSurvIndPlot,Estimated=predSurvBio,Observed=obsSurvBio,quantity="Survey")

#==plot the estimated processes
par(mfrow=c(2,1),mar=c(.1,3,.1,.1),oma=c(2,3,0,0))
PolygonPlots(Truth=trueRecPlot,Estimated=TotRec,AddLegend=T,quantity="Recruitment")
PolygonPlots(Truth=trueFmortPlot,Estimated=TotFdir,quantity="Fishing mortality")

#==SELECTIVITY==
# selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(surv50n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(surv95n[1],VonKn[1],LinfN[1],t0n[1])
# pull pars from OM
# translate to age with function
# calculate
# plot

#==plot the reference points
PolygonPlots(Truth=trueSpbioPlot,Estimated=predSpBio)
  plot(trueRecPlot~trueSpbioPlot,ylim=c(0,max(trueRecPlot,na.rm=T)))

#==errors in reference points
REbmsy<-(predB35PlotEnd-trueBMSYPlot)/trueBMSYPlot
REfmsy<-(predF35PlotEnd-trueFMSYPlot)/trueFMSYPlot
par(mfrow=c(2,1),mar=c(.1,3,.1,.1),oma=c(2,3,0,0))
boxplot(REbmsy)
boxplot(REfmsy)
## relative errors in everything
## summary table to compare between other scenarios
PlotSel<-0
if(PlotSel==1)
{
#dev.new()
par(mfrow=c(4,1),mar=c(1,4,1,1),oma=c(3,1,1,1))
#==plot estimated selectivities

cnt<-grep("srv_sel50",data2)
sel50surv<-as.numeric(data2[cnt+1])
cnt<-grep("srv_sel95",data2)
sel95surv<-as.numeric(data2[cnt+1])

cnt<-grep("SelPars50",data2)
sel50fish<-as.numeric(data2[cnt+1])
cnt<-grep("SelPars95",data2)
sel95fish<-as.numeric(data2[cnt+1])

trueSel50<-out$OM$sel50s
trueSel95<-out$OM$sel95s
trueSurv50<-out$OM$surv50s
trueSurv95<-out$OM$surv95s

SurvSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50surv)/(sel95surv-sel50surv)))
FishSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50fish)/(sel95fish-sel50fish)))

trueSurvSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueSurv50)/(trueSurv95-trueSurv50)))
trueFishSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueSel50)/(trueSel95-trueSel50)))

plot(trueSurvSel~LengthBin,type="l",xlab="Length",ylab="Selectivity")

lines(FishSel~AgeBin,lty=2,col=2)
lines(SurvSel~AgeBin,lty=2,col=1)
lines(trueFishSel~LengthBin,lty=1,col=2)
legend("bottomright",bty='n',col=c(1,1,2,2),lty=c(1,2,1,2),legend=c("True Survey","Est Survey","True Fishery","Est Fishery"))
}

#dev.new()
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsSurlen[i,,1],axes=F)
lines(predSurlen[i,,1])
}
plot.new()
legend("center","Survey")
#dev.new()
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsCatchlen[i,,1],axes=F)
lines(predCatchlen[i,,1])
}
plot.new()
legend("center","Catch")
dev.off()
}

}
########################################################
PullTimevary<-function(inFolders,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
MaxAge		<-out$OM$MaxAge
Ages			<-seq(1,MaxAge)
t0			<-out$OM$t0s		# check this later when going to spatial

Recruitment	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))
FishMort	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))
TrueRec	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))
TrueFmort	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))
TrueCatch	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))
TrueSpbio	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))
EstSpbio	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))
NatMvary	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))
SelVary	<-array(dim=c(SimYear,MaxAge,Nsim,length(inFolders))) 
GrowthVary	<-array(dim=c(SimYear,MaxAge,Nsim,length(inFolders)))
TrueGrow	<-array(dim=c(nrow=SimYear,ncol=MaxAge,length(inFolders)))
TrueSel	<-array(dim=c(nrow=SimYear,ncol=MaxAge,length(inFolders)))
TrueM   <-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))
TrueOFL	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))
EstOFL	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(inFolders)))

for(p in 1:length(inFolders))
{
DrawDir		<-inFolders[p]
Inout			<-ReadCTLfile(file.path(MSEdir,paste0(inFolders[p],".csv")))

IndSimFolder	<-file.path(DrawDir,Nsim,SimYear)
TRU			<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))

#==true
temp		<-grep("fishing mortality",TRU)
TrueFmort[,,p]<-matrix(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" "))),ncol=SimYear,byrow=T)

#==true Catch
temp		<-grep("Catch",TRU)
TrueCatch[,,p]<-matrix(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" "))),ncol=SimYear,byrow=T)

#==true Spawning biomass
temp		<-grep("spawning biomass",TRU)
TrueSpbio[,,p]<-matrix(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" "))),ncol=SimYear,byrow=T)

temp		<-grep("recruitment",TRU)
TrueRec[,,p]<-matrix(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" "))),ncol=SimYear,byrow=T)

#==OFL 
temp			<-grep("OFL",TRU)
TrueOFL[,,p]	<-matrix(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" "))),ncol=SimYear,byrow=T)

for(n in 1:Nsim)
{
IndSimFolder	<-file.path(DrawDir,n,SimYear)
PAR			<-readLines(file.path(MSEdir,IndSimFolder,"simass.par"))

#==fishing mortality
temp		<-grep("log_avg_fmort_dir",PAR)[1]
AvgF		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp		<-grep("fmort_dir_dev",PAR)[1]
fDevs		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
FishMort[n,,p]<-exp(AvgF+fDevs)

#==est spawning biomass
IndSimFolder2			<-file.path(DrawDir,n,SimYear)
pullREP				<-readLines(file.path(MSEdir,IndSimFolder2,"simass.rep"))
temp					<-grep("spawning biomass",pullREP)
EstSpbio[n,,p]			<-as.numeric((unlist(strsplit(pullREP[temp+1],split=" "))))[2:(SimYear+1)]

#==recruitment
temp		<-grep("mean_log_rec",PAR)[1]
AvgRec	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp		<-grep("rec_dev",PAR)[1]
RecDevs	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Recruitment[n,,p]<-exp(AvgRec+RecDevs)

for(w in (InitYear+1):SimYear)
{
IndSimFolder2			<-file.path(DrawDir,n,w)
pullREP				<-readLines(file.path(MSEdir,IndSimFolder2,"simass.rep"))
temp					<-grep("OFL",pullREP)[2]
EstOFL[n,w,p]			<-as.numeric((unlist(strsplit(pullREP[temp+1],split=" "))))
}

#==timevary
#==growth combos
#==write a function already that pulls the variable...

if(Inout$OM$TimeVaryLinf>0 & Inout$OM$TimeVaryGrowthK <=0)
{
temp			<-grep("log_avg_Linf",PAR)[1]
log_avg_Linf	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("Linf_dev",PAR)[1]
Linf_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Linf			<-exp(log_avg_Linf+Linf_dev)

GrowthK		<-rep(Inout$OM$VonKn,SimYear)
}

if(Inout$OM$TimeVaryLinf<=0 & Inout$OM$TimeVaryGrowthK >0)
{
temp			<-grep("log_avg_GrowthK",PAR)[1]
log_avg_GrowthK	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("GrowthK_dev",PAR)[1]
GrowthK_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
GrowthK		<-exp(log_avg_GrowthK+GrowthK_dev)

Linf			<-rep(Inout$OM$LinfN,SimYear)
}

if(Inout$OM$TimeVaryLinf>0 & Inout$OM$TimeVaryGrowthK>0)
{
temp			<-grep("log_avg_Linf",PAR)[1]
log_avg_Linf	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("Linf_dev",PAR)[1]
Linf_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Linf			<-exp(log_avg_Linf+Linf_dev)

temp			<-grep("log_avg_GrowthK",PAR)[1]
log_avg_GrowthK	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("GrowthK_dev",PAR)[1]
GrowthK_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
GrowthK		<-exp(log_avg_GrowthK+GrowthK_dev)
}

if(Inout$OM$TimeVaryLinf<=0 & Inout$OM$TimeVaryGrowthK<=0)
{
Linf			<-rep(Inout$OM$LinfN,SimYear)
GrowthK		<-rep(Inout$OM$VonKn,SimYear)
}

for(x in 1:SimYear)
 GrowthVary[x,,n,p]<-Linf[x]*(1-exp(-GrowthK[x]*(Ages-t0)))

#==natural mortality
if(Inout$OM$TimeVaryM >0)
{
temp			<-grep("log_avg_NatM",PAR)[1]
avgM			<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("m_dev",PAR)[1]
mDevs			<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
NatMvary[n,,p]	<-exp(avgM+mDevs)
}

if(Inout$OM$TimeVaryM <=0 )
{
  temp			<-grep("log_avg_NatM",PAR)[1]
  avgM			<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
  NatMvary[n,,p]	<-exp(avgM)
}

#==fishing selectivity combos
# selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])


if(Inout$OM$TimeVarySel50 >0 & Inout$OM$TimeVarySel95 <=0)
{
temp			<-grep("SelPars50",PAR)[1]
SelPars50		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev50",PAR)[1]
SelPars_dev50	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel50			<-SelPars50+SelPars_dev50

temp			<-grep("SelPars95",PAR)[1]
SelPars95		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
}

if(Inout$OM$TimeVarySel50 <=0 & Inout$OM$TimeVarySel95 >0)
{
temp			<-grep("SelPars95",PAR)[1]
SelPars95		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev95",PAR)[1]
SelPars_dev95	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel95			<-SelPars95+SelPars_dev95

temp			<-grep("SelPars50",PAR)[1]
SelPars50		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
}

if(Inout$OM$TimeVarySel50>0 & Inout$OM$TimeVarySel95 >0)
{
temp			<-grep("SelPars95",PAR)[1]
SelPars95		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev95",PAR)[1]
SelPars_dev95	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel95			<-SelPars95+SelPars_dev95

temp			<-grep("SelPars50",PAR)[1]
SelPars50		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev50",PAR)[1]
SelPars_dev50	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel50			<-SelPars50+SelPars_dev50
}

if(Inout$OM$TimeVarySel50<=0 & Inout$OM$TimeVarySel95 <=0)
{
temp			<-grep("SelPars50",PAR)[1]
Sel50		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
temp			<-grep("SelPars95",PAR)[1]
Sel95		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
}

for(x in 1:SimYear)
 SelVary[x,,n,p]	<-1/( 1 + exp( -1*log(19)*(Ages-Sel50[x])/(Sel95[x]-Sel50[x])))

}
#==============================================
# Pull the 'true' values for each process
#==============================================
tGrowthK		<-Inout$OM$VonKn
if(length(tGrowthK)==1)
 tGrowthK		<-rep(tGrowthK,SimYear)
tLinf			<-Inout$OM$LinfN
if(length(tLinf)==1)
 tLinf		<-rep(tLinf,SimYear)

for(x in 1:SimYear)
 TrueGrow[x,,p]	<-tLinf[x]*(1-exp(-tGrowthK[x]*(Ages-t0)))

 tSel50<-selAtAgeFunc(Inout$OM$sel50n,Inout$OM$VonKn,Inout$OM$LinfN,Inout$OM$t0n)
 tSel95<-selAtAgeFunc(Inout$OM$sel95n,Inout$OM$VonKn,Inout$OM$LinfN,Inout$OM$t0n)

if(length(tSel50)==1)
 tSel50		<-rep(tSel50,Inout$OM$SimYear)
if(length(tSel95)==1)
 tSel95		<-rep(tSel95,Inout$OM$SimYear)

for(x in 1:SimYear)
 TrueSel[x,,p]	<-1/( 1 + exp( -1*log(19)*(Ages-tSel50[x])/(tSel95[x]-tSel50[x])))

TrueM			<-Inout$OM$NatMn
if(length(TrueM)==1)
 TrueM		<-rep(TrueM,SimYear)
}

list(Recruitment,TrueRec,FishMort,TrueFmort,GrowthVary,TrueGrow,NatMvary,TrueM,SelVary,TrueSel,EstOFL,TrueOFL,TrueCatch,TrueSpbio,EstSpbio)
}

########################################################
PullTrueProj<-function(inFolders,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
MaxAge		<-out$OM$MaxAge
Ages			<-seq(1,MaxAge)
t0			<-out$OM$t0s		# check this later when going to spatial

TrueRec	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))
TrueFmort	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))
TrueCatch	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))
TrueSpbio	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))
TrueOFL	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(inFolders)))

for(p in 1:length(inFolders))
{
DrawDir		<-inFolders[p]
Inout			<-ReadCTLfile(file.path(MSEdir,paste0(inFolders[p],"_CTL.csv")))

for(n in 1:Nsim)
{
IndSimFolder	<-file.path(DrawDir,n,SimYear)
TRU			<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))

#==true
temp		<-grep("fishing mortality",TRU)
TrueFmort[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]

#==true Catch
temp		<-grep("Catch",TRU)
TrueCatch[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]

#==true Spawning biomass
temp		<-grep("spawning biomass",TRU)
TrueSpbio[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]


#==recruitment
temp		<-grep("recruitment",TRU)
TrueRec[n,,p]<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]

#==OFL 
temp			<-grep("OFL",TRU)
TrueOFL[n,,p]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(SimYear-1)]
}
}
list(TrueRec,TrueFmort,TrueOFL,TrueCatch,TrueSpbio)
}

#=============================================================
# Production model plots
#============================================================
ProductionModelOutput<-function(Inout,CTLNames,MSEdir)
{
Data<-rep(list(list()))
for(x in 1:length(CTLNames))
 Data[[x]]<-readLines(file.path(MSEdir,CTLNames[x],"ProdOutputs.csv"))

#==estimated and true CPUE
estCPUE<-array(dim=c(Inout$OM$Nsim,Inout$OM$SimYear,length(CTLNames)))
trueCPUE<-array(dim=c(Inout$OM$Nsim,Inout$OM$SimYear,length(CTLNames)))
trueBMSY<-rep(0,length(CTLNames))
estBMSY<-array(dim=c(Inout$OM$Nsim,Inout$OM$SimYear,length(CTLNames)))
trueFMSY<-rep(0,length(CTLNames))
estFMSY<-array(dim=c(Inout$OM$Nsim,Inout$OM$SimYear,length(CTLNames)))
estTAC<-array(dim=c(Inout$OM$Nsim,Inout$OM$SimYear,length(CTLNames)))
trueTAC<-array(dim=c(Inout$OM$Nsim,Inout$OM$SimYear,length(CTLNames)))

for(y in 1:length(CTLNames))
{
 tmp<-grep("est cpue",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estCPUE[x,,y]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

tmp<-grep("true CPUE",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 trueCPUE[x,,y]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

tmp<-grep("true BMSY",Data[[y]])
trueBMSY[y]<-as.numeric(Data[[y]][(tmp+1)])

tmp<-grep("est BMSY",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estBMSY[x,,y]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

tmp<-grep("true FMSY",Data[[y]])
trueFMSY[y]<-as.numeric(Data[[y]][(tmp+1)])

tmp<-grep("est FMSY",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estFMSY[x,,y]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

 tmp<-grep("est total allowable catch",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estTAC[x,,y]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

tmp<-grep("true total allowable catch",Data[[y]])
Quant<-Data[[y]][(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 trueTAC[x,,y]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

}

png(file.path(MSEdir,"plots",paste0("ProductionFits_",paste(CTLNames,collapse="_"),".png")),res=1200,width=5,height=4.5,units='in')
par(mfcol=c(2,length(CTLNames)),mar=c(.1,.1,.1,.1),oma=c(4,6,1,4))

for(y in 1:length(CTLNames))
{
# boxplot(trueCPUE[,,y],type="l",ylim=c(0,max(trueCPUE,na.rm=T)),
#  las=1,xaxt='n',ylab='',yaxt='n')
#   input<-trueCPUE[,,y]
  # for(x in 1:nrow(estCPUE))
  #   lines(estCPUE[x,,y],col='#ff000033')
  
PolygonPlots(Truth=trueCPUE[,,y],Estimated=estCPUE[,,y],SimYear=Inout$OM$SimYear,Nsim=Inout$OM$Nsim)

if(y==1)
{
 axis(side=2,las=1)
 mtext(side=2,"Biomass",line=5,cex=.9)
}

abline(h=trueBMSY,col="#0000ff99",lty=1)
abline(h=max(estBMSY,na.rm=T),col="#00800099",lty=2)
abline(h=min(estBMSY,na.rm=T),col="#00800099",lty=2)

legend("topright",col=c(1,2,"#0000ff99","#00800077"),pch=c(15,NA,NA,NA),lty=c(NA,1,1,2),
       legend=c("Observations","Estimates","True BMSY","Estimated BMSY range"),bty='n',cex=.7)


PolygonPlots(Truth=trueTAC[,,y],Estimated=estTAC[,,y],SimYear=Inout$OM$SimYear,Nsim=Inout$OM$Nsim)

# boxplot(trueTAC[,,y],type="l",ylim=c(0,max(trueTAC,na.rm=T)),
#         las=1,xaxt='n',ylab='',yaxt='n')
# for(x in 1:nrow(estCPUE))
#   lines(estTAC[x,,y],col='#ff000033')
if(y==1)
{
  # axis(side=2,las=1)
  mtext(side=2,"Total allowable catch",line=4.5,cex=.9)
}
}
axis(side=1)

legend('topleft',col=c(1,2),pch=c(15,NA),lty=c(NA,1),legend=c("True","Estimated"),bty='n',cex=.7)
dev.off()

png(file.path(MSEdir,"plots",paste0("ProductionRefPoints_",paste(CTLNames,collapse="_"),".png")),res=1200,width=5,height=4.5,units='in')
par(mfcol=c(2,length(CTLNames)),mar=c(.1,.1,.1,.1),oma=c(4,6,1,4))
for(y in 1:length(CTLNames))
{
temp<-sweep(estBMSY[,,y],MAR=2,trueBMSY[y],FUN="-")
RelativeErrorBMSY<-sweep(temp,MAR=2,trueBMSY[y],FUN="/")
temp<-sweep(estFMSY[,,y],MAR=2,trueFMSY[y],FUN="-")
RelativeErrorFMSY<-sweep(temp,MAR=2,trueFMSY[y],FUN="/")

yUp   <-max(RelativeErrorBMSY,RelativeErrorFMSY,na.rm=T)
ydown <-min(RelativeErrorBMSY,RelativeErrorFMSY,na.rm=T)

inShape<-apply(RelativeErrorBMSY[,1:(ncol(RelativeErrorBMSY))],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)
plot(-100000000000,ylim=c(ydown,yUp),las=1,ylab="",xlab="Year",xaxt='n',xlim=c(Inout$OM$InitYear,Inout$OM$SimYear))
polygon(x=c(seq(1,(Inout$OM$SimYear-1)),rev(seq(1,(Inout$OM$SimYear-1)))),y=c(inShape[1,1:Inout$OM$SimYear-1],rev(inShape[4,1:Inout$OM$SimYear-1])),col='darkgrey',border=F)
polygon(x=c(seq(1,(Inout$OM$SimYear-1)),rev(seq(1,(Inout$OM$SimYear-1)))),y=c(inShape[2,1:Inout$OM$SimYear-1],rev(inShape[3,1:Inout$OM$SimYear-1])),col='lightgrey',border=F)
abline(h=0,lty=2)

if(y==1)
{
 axis(side=2,las=1)
 mtext(side=2,"Relative error",line=2.5,cex=.9)
 mtext(side=2,"Target biomass",line=3.5,cex=.9)
}

inShape<-apply(RelativeErrorFMSY[,1:(ncol(RelativeErrorFMSY))],2,quantile,probs=c(.05,.25,.75,.95),na.rm=T)
plot(-100000000000,ylim=c(ydown,yUp),las=1,ylab="",xlab="Year",xaxt='n',xlim=c(Inout$OM$InitYear,Inout$OM$SimYear))
polygon(x=c(seq(1,(Inout$OM$SimYear-1)),rev(seq(1,(Inout$OM$SimYear-1)))),y=c(inShape[1,1:Inout$OM$SimYear-1],rev(inShape[4,1:Inout$OM$SimYear-1])),col='darkgrey',border=F)
polygon(x=c(seq(1,(Inout$OM$SimYear-1)),rev(seq(1,(Inout$OM$SimYear-1)))),y=c(inShape[2,1:Inout$OM$SimYear-1],rev(inShape[3,1:Inout$OM$SimYear-1])),col='lightgrey',border=F)
abline(h=0,lty=2)

if(y==1)
{
 axis(side=2,las=1)
 mtext(side=2,"Relative error",line=2.5,cex=.9)
 mtext(side=2,"Target fishing mortality",line=3.5,cex=.9)
}
}

dev.off()
}

#############################################################
#==cOMPARISON OF AGE STRUCTURED MODELS
############################################################
AgeStructureComp<-function(Inout,RetroPeels=6,CTLNames,MSEdir)
{
  TakeRows<-(Inout$OM$SimYear-RetroPeels+1):Inout$OM$SimYear
  GradientSave<-array(dim=c(RetroPeels,Inout$OM$Nsim,length(CTLNames)))
  
  for(x in 1:length(CTLNames))
  {
    Inout<-ReadCTLfile(file.path(MSEdir,paste0(CTLNames[x],".csv")))
    GradientSave[,,x]<-CheckGradient(DrawDir=CTLNames[x],out=Inout,MSEdir)[TakeRows,]
  }
  ScenCols<-c("grey",as.numeric(seq(2,length(CTLNames),1)))
  
  #==reference points
  TakeRows	<-(Inout$OM$InitYear+1):Inout$OM$SimYear
  B35save	<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CTLNames)))
  F35save	<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CTLNames)))
  OFLsave	<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CTLNames)))
  tOFLsave	<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CTLNames)))
  tF35save	<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CTLNames)))
  tB35save	<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CTLNames)))
  
  for(x in 1:length(CTLNames))
  {
    temp<-CheckReferencePoints(DrawDir=CTLNames[x],out=Inout,MSEdir)
    B35save[,,x]	<-temp[[1]][TakeRows,]
    F35save[,,x]	<-temp[[2]][TakeRows,]
    OFLsave[,,x]	<-temp[[3]][TakeRows,]
    tB35save[,,x]	<-temp[[4]][TakeRows,]
    tF35save[,,x]	<-temp[[5]][TakeRows]
    tOFLsave[,,x]	<-temp[[6]][TakeRows,]
  } 
  
  BigB35<-matrix(nrow=Inout$OM$Nsim,ncol=length(CTLNames))
  BigF35<-matrix(nrow=Inout$OM$Nsim,ncol=length(CTLNames))
  BigOFL<-matrix(nrow=Inout$OM$Nsim,ncol=length(CTLNames))
  if(Inout$OM$Nsim>1) {
    for(x in 1:length(CTLNames))
    {
      BigB35[,x]<-apply((B35save[,,x]-tB35save[,,x])/tB35save[,,x],2,median,na.rm=T)
      BigF35[,x]<-apply((F35save[,,x]-tF35save[,,x])/tF35save[,,x],2,median,na.rm=T)
      BigOFL[,x]<-apply((OFLsave[,,x]-tOFLsave[,,x])/tOFLsave[,,x],2,median,na.rm=T)
    }
  }
  
  if(Inout$OM$Nsim==1) {
    for(x in 1:length(CTLNames))
    {
      BigB35[,x]<-median((B35save[,,x]-tB35save[,,x])/tB35save[,,x],na.rm=T)
      BigF35[,x]<-median((F35save[,,x]-tF35save[,,x])/tF35save[,,x],na.rm=T)
      BigOFL[,x]<-median((OFLsave[,,x]-tOFLsave[,,x])/tOFLsave[,,x],na.rm=T)
    }
  }
  
  ReB35<-BigB35
  ReF35<-BigF35
  ReOFL<-BigOFL
  
  #=plotting output
  # GeMSplots(out=Inout,SourceFolder=CTLNames[x])
  MohnsRho<-array(dim=c(RetroPeels-1,Inout$OM$Nsim,length(CTLNames)))
  SSBbias<-array(dim=c(RetroPeels,Inout$OM$Nsim,length(CTLNames)))
  TrueOFL<-array(dim=c(RetroPeels,Inout$OM$Nsim,length(CTLNames)))
  
  for(x in 1:length(CTLNames))
  {
    Inout			<-ReadCTLfile(file.path(MSEdir,paste0(CTLNames[x],".csv")))
    DrawDir		<-CTLNames[x]
    RetroOuts		<-CheckRetro(RetroPeels=RetroPeels,DrawDir,PlotRetro=0,out=Inout,MSEdir)
    MohnsRho[,,x]	<-CalculateMohn(RetroOuts)
    SSBbias[,,x]	<-CalculateBias(RetroOuts)
  }
  
  BigMohn<-matrix(nrow=Inout$OM$Nsim,ncol=length(CTLNames))
  if(Inout$OM$Nsim > 1) {
    for(x in 1:length(CTLNames))
    {
      temp<-MohnsRho[,,x]
      BigMohn[,x]<-apply(temp,2,mean) 
    }
  }
  
  if(Inout$OM$Nsim == 1) {
    for(x in 1:length(CTLNames))
    {
      temp<-MohnsRho[,,x]
      BigMohn[,x]<-mean(temp) 
    }
  }
  
  BigBias<-matrix(nrow=Inout$OM$Nsim,ncol=length(CTLNames))
  if(Inout$OM$Nsim > 1) {
    for(x in 1:length(CTLNames))
    {
      temp<-SSBbias[,,x]
      BigBias[,x]<-apply(temp,2,mean) 
    }
  }
  
  if(Inout$OM$Nsim == 1) {
    for(x in 1:length(CTLNames))
    {
      temp<-SSBbias[,,x]
      BigBias[,x]<-mean(temp) 
    }
  }
  
  png(file.path(MSEdir,"plots","CompareRefPoints.png"),height=7,width=3.5,units='in',res=1200)
  par(mfrow=c(5,1),mar=c(.1,.1,.3,.1),oma=c(4,6,1,1))
  
  inYlim<-c(min(BigMohn,BigBias,ReB35,ReF35,BigOFL),max(BigMohn,BigBias,ReB35,ReF35,BigOFL))
  boxplot(BigMohn,col=ScenCols,ylim=inYlim,xaxt='n',las=1)
  legend("topleft",c("(a) Retrospective bias"),bty='n')
  mtext(side=2,"Mohn's rho",line=3,cex=.7)
  abline(h=0,lty=2)
  
  boxplot(BigBias,col=ScenCols,xaxt='n',ylim=inYlim,las=1)
  legend("topleft",c("(b) Spawning biomass"),bty='n')
  mtext(side=2,"relative error",line=3,cex=.7)
  abline(h=0,lty=2)
  boxplot(ReB35,col=ScenCols,,ylim=inYlim,xaxt='n',las=1)
  legend("topleft",expression("(c) B"[35]),bty='n')
  mtext(side=2,"relative error",line=3,cex=.7)
  abline(h=0,lty=2)
  boxplot(ReF35,col=ScenCols,ylim=inYlim,xaxt='n',las=1)
  legend("topleft",expression("(c) F"[35]),bty='n')
  mtext(side=2,"relative error",line=3,cex=.7)
  abline(h=0,lty=2)
  boxplot(BigOFL,col=ScenCols,ylim=inYlim,las=1,names=CTLNames)
  abline(h=0,lty=2)
  legend("topleft",c("(e) OFL"),bty='n')
  mtext(side=2,"relative error",line=3,cex=.7)
  dev.off()
  
  quants<-PullTimevary(inFolders=CTLNames,out=Inout,MSEdir)
  
  png(file.path(MSEdir,"plots","ComparePopulationProcess.png"),height=5,width=7.5,units='in',res=1200)
  inmat<-matrix(c(1,1,2,2,3,3,
                  4,4,4,5,5,5),nrow=2,byrow=T)
  layout(inmat)
  par(mar=c(.1,.1,.1,.1),oma=c(4,5,4,4))
  
  #==RECRUITMENT
  input<-quants[[1]]
  tInput<-quants[[2]] 
  
  plot(-100000,xlim=c(1,dim(input)[2]),ylim=c(0,max(input,na.rm=T)),las=1,xaxt='n')
  for(x in 1:length(CTLNames))
  {
    color<-seq(1,length(CTLNames)+1)
    incol<-adjustcolor(color,alpha.f=.08)
    tCol<-rgb(0,0,0,0.2)
    temp<-input[,,x]
    for(y in 1:dim(input)[1])
      lines(input[y,,x],col=incol[x])
  } 
  medians<- apply(input,c(2,3),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=x,cex=2)
  
  abline(v=Inout$OM$InitYear,lty=3)
  
  #==true R
  plotIn<-apply(tInput,2,median)
  plotIn[1]<-NA
  lines(plotIn,col="yellow",lwd=2,lty=2)
  
  legend("bottomleft",bty='n',"(a) Recruitment")
  mtext(side=2,"Numbers",line=4,cex=.7)
  #==FISHING MORTALITY
  input<-quants[[3]]
  NatM<-quants[[7]]
  tInput<-quants[[4]] 
  
  plot(-1000,xlim=c(1,dim(input)[2]),ylim=c(0,max(input,NatM,na.rm=T)),las=1,xaxt='n',yaxt='n')
  axis(side=3)
  mtext(side=3,line=2,"Time")
  for(x in 1:length(CTLNames))
  {
    color<-seq(1,length(CTLNames)+1)
    incol<-adjustcolor(color,alpha.f=.08)
    tCol<-rgb(0,0,0,0.2)
    temp<-input[,,x]
    for(y in 1:dim(input)[1])
      lines(input[y,,x],col=incol[x])
  } 
  medians<- apply(input,c(2,3),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=x,cex=2)
  
  abline(v=Inout$OM$InitYear,lty=3)
  
  #==true F
  lines(apply(tInput,2,median),col="yellow",lwd=2,lty=2)
  legend("bottomleft",bty='n',"(b) Fishing mortality")
  
  #==NATURAL MORTALITY
  input<-quants[[7]]
  tInput<-quants[[8]] 
  
  plot(-1000,xlim=c(1,dim(input)[2]),ylim=c(0,max(input,NatM,na.rm=T)),las=1,xaxt='n',yaxt='n')
  axis(side=4,las=1)
  for(x in 1:length(CTLNames))
  {
    color<-seq(1,length(CTLNames)+1)
    incol<-adjustcolor(color,alpha.f=.08)
    tCol<-rgb(0,0,0,0.2)
    temp<-input[,,x]
    for(y in 1:dim(input)[1])
      lines(input[y,,x],col=incol[x])
  } 
  medians<- apply(input,c(2,3),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=x,cex=2)
  
  #==true M
  dim(tInput)
  lines(tInput,col="yellow",lwd=2,lty=2)
  
  abline(v=Inout$OM$InitYear,lty=3)
  
  legend("bottomleft",bty='n',"(c) Natural mortality")
  mtext(side=4,"rate (/year)",line=2,cex=.7)
  
  #==SELECTIVITY
  input<-quants[[9]]
  tInput<-quants[[10]] 
  plot(-1000,xlim=c(1,which(input[1,,1,1]>0.99)[3]),ylim=c(0,max(input,na.rm=T)),las=1)
  for(x in 1:length(CTLNames))
  {
    color<-seq(1,length(CTLNames)+1)
    incol<-adjustcolor(color,alpha.f=.008)
    tCol<-rgb(0,0,0,0.2)
    temp<-input[,,,x]
    for(y in 1:dim(input)[1])
    {
      #lines(tInput[y,,x],col=tCol)
      for(z in 1:dim(input)[3])
        lines(input[y,,z,x],col=incol[x])
    }   
  } 
  
  medians<- apply(input,c(2,4),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=x,cex=2)
  
  lines(tInput[1,,1],col="yellow",lwd=2,lty=2)
  lines(tInput[Inout$OM$SimYear,,1],col="yellow",lwd=2,lty=2)
  mtext(side=2,"Selectivity",line=2.3,cex=.7)
  legend("topleft",bty='n',"(d) Fishery selectivity")
  
  #==GROWTH
  input<-quants[[5]]
  tInput<-quants[[6]] 
  plot(-1000,xlim=c(1,Inout$OM$MaxAge),ylim=c(0,max(input,na.rm=T)),las=1,yaxt='n')
  axis(side=4,las=1)
  for(x in 1:length(CTLNames))
  {
    color<-seq(1,length(CTLNames)+1)
    incol<-adjustcolor(color,alpha.f=.008)
    tCol<-rgb(0,0,0,0.2)
    temp<-input[,,,x]
    for(y in 1:dim(input)[1])
    {
      #lines(tInput[y,,x],col=tCol)
      for(z in 1:dim(input)[3])
        lines(input[y,,z,x],col=incol[x])
    }   
  } 
  
  medians<- apply(input,c(2,4),median,na.rm=T)
  for(x in 1:ncol(medians))
    lines(medians[,x],col=x,cex=2)
  mtext(side=4,"Length (cm)",line=2.3,cex=.7)
  
  #==TRUE GROWTH
  lines(tInput[1,,1],col="yellow",lwd=2,lty=2)
  lines(tInput[Inout$OM$SimYear,,1],col="yellow",lwd=2,lty=2)
  legend("topleft",bty='n',"(e) Length at age")
  
  legend("bottomright",bty='n',col=c(1,ScenCols,"yellow"),lty=c(3,rep(1,length(ScenCols)),2),
         legend=c("Projection begins",CTLNames,"Truth"),lwd=2)
  mtext(side=1,outer=T,"Age",line=2.3)
  
  dev.off()
  
  
  #================================
  # Plot the comparisons of fits like the production models
  #====================================
  
  
  png(file.path(MSEdir,"plots",paste0("AgeStructuredFits_",CTLNames,".png")),height=600,width=900)
  par(mfcol=c(2,length(CTLNames)),mar=c(.1,.1,.1,.1),oma=c(4,6,1,4))
  
  if(Inout$OM$Nsim>1) {
    for(y in 1:length(CTLNames))
    {
    boxplot(quants[[14]][,,y],type="l",ylim=c(0,max(quants[[14]],na.rm=T)),
            las=1,xaxt='n',ylab='',yaxt='n')
    for(x in 1:nrow(  quants[[14]]))
      lines(quants[[15]][x,,y],col='#ff000033')
    if(y==1)
    {
      axis(side=2,las=1)
      mtext(side=2,"Biomass",line=5,cex=.9)
    }
    legend("topright",bty='n',,legend=CTLNames[y])
   # abline(h=trueBMSY,col="#0000ff99",lty=1)
  #  abline(h=max(estBMSY,na.rm=T),col="#00800099",lty=2)
   # abline(h=min(estBMSY,na.rm=T),col="#00800099",lty=2)
    
    #legend("topright",col=c(1,2,"#0000ff99","#00800077"),pch=c(15,NA,NA,NA),lty=c(NA,1,1,2),
    #       legend=c("Observations","Estimates","True BMSY","Estimated BMSY range"),bty='n')
    
    boxplot(quants[[12]][,,y],type="l",ylim=c(0,max(quants[[12]],na.rm=T)),
            las=1,xaxt='n',ylab='',yaxt='n')
    for(x in 1:nrow(quants[[11]]))
      lines(quants[[11]][x,,y],col='#ff000033')
    if(y==1)
    {
      axis(side=2,las=1)
      mtext(side=2,"Total allowable catch",line=4.5,cex=.9)
    }
    axis(side=1)  
  }
  }
  
    if(Inout$OM$Nsim==1) {
    for(y in 1:length(CTLNames))
    {
    plot(quants[[14]][,,y],type="l",ylim=c(0,max(quants[[14]],na.rm=T)),
            las=1,xaxt='n',ylab='',yaxt='n')
    for(x in 1:nrow(  quants[[14]]))
      lines(quants[[15]][x,,y],col='#ff000033')
    if(y==1)
    {
      axis(side=2,las=1)
      mtext(side=2,"Biomass",line=5,cex=.9)
    }
    legend("topright",bty='n',,legend=CTLNames[y])
   # abline(h=trueBMSY,col="#0000ff99",lty=1)
  #  abline(h=max(estBMSY,na.rm=T),col="#00800099",lty=2)
   # abline(h=min(estBMSY,na.rm=T),col="#00800099",lty=2)
    
    #legend("topright",col=c(1,2,"#0000ff99","#00800077"),pch=c(15,NA,NA,NA),lty=c(NA,1,1,2),
    #       legend=c("Observations","Estimates","True BMSY","Estimated BMSY range"),bty='n')
    
    plot(quants[[12]][,,y],type="l",ylim=c(0,max(quants[[12]],na.rm=T)),
            las=1,xaxt='n',ylab='',yaxt='n')
    for(x in 1:nrow(quants[[11]]))
      lines(quants[[11]][x,,y],col='#ff000033')
    if(y==1)
    {
      axis(side=2,las=1)
      mtext(side=2,"Total allowable catch",line=4.5,cex=.9)
    }
    axis(side=1)  
  }
  }
  
  legend('topleft',col=c(1,2),pch=c(15,NA),lty=c(NA,1),legend=c("True","Estimated"),bty='n')
  dev.off()
  
  
  #==================================
  # relative error in OFL over time
  # NEED TO GET THIS ONCE FIGURED OUT AT SOME POINT
  #====================================
  OFLrel<-0
  if(OFLrel>0)
  {
  InitYear<-Inout$OM$InitYear
  SimYear<-Inout$OM$SimYear
  
  REofl<-rep(list(list()))
  for(x in 1:length(SaveTrueOFL))
    REofl[[x]]<-(SaveTrueOFL[[x]][,(InitYear+1):(SimYear-1),]-SaveTrueOFL[[x]][,(InitYear+1):(SimYear-1),])/SaveTrueOFL[[x]][,(InitYear+1):(SimYear-1),]
  
  emNames<-c("Base","M vary","Growth vary","Sel vary")
  omNames<-c("Base","M vary","Sel vary","Growth vary")
  hexDec<-c("#ff000055","#00ff0055","#0000ff55")
  par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
  for(x in 1:length(REofl))
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
      lines(apply(sortQuant,2,median),col=y,lwd=2)
      
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
  for(x in 1:length(REofl))
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
