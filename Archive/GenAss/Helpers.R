###############################################################
SBPRfunc<-function(F,Rzero)
{
#initial conditions
N<-matrix(ncol=maxAge,nrow=50)
N[1,]<-initialN(Rzero)

for (j in 2:50)
{
 for (i in 2:(maxAge-1))
  N[j,i]<-N[j-1,i-1]*exp(-F*vuln[i-1])*survival

 N[j,maxAge]<-(N[j-1,maxAge-1]+N[j-1,maxAge])*exp(-F*vuln[maxAge])*survival
 N[j,1]<-Rzero
}
SBP<-sum(N[50,]*mature*weight)
SBPR<-sum(N[50,]*mature*weight)/Rzero
YPR<-sum(N[50,]*(1-exp(-F*vuln)) *weight)
list(SBP,SBPR,YPR)
#return(SBPR)
}

###############################################################

Recruitment<-function()
{
if(recType=="BH")
{
 Szero<-unlist(SBPRfunc(0,Rzero[j])[1])	#since there are only a couple possible values, code could be sped up by changing this
 alpha<-(Szero/Rzero[j])*(1-steepness)/(4*steepness)
 beta<-(5*steepness-1)/(4*steepness*Rzero[j])
 Recruit<-(Eggs/(alpha+beta*Eggs))*exp(RecErr[j])
}
if(recType=="Rick")
{
 Szero<-unlist(SBPRfunc(0,Rzero[j])[1])
 beta<-log(5*steepness)/(0.8*Szero)
 alpha<-log(Rzero[j]/Szero)+(beta*Szero)
 Recruit<-(Eggs*exp(alpha-beta*Eggs))*exp(RecErr[j])
}
return(Recruit)
}

###############################################################

initialN<-function(Rzero)
{
test		<-rep(0,100)
test[1]	<-Rzero[1]
for (i in 2:100)
 {test[i]	<-test[i-1]*survival}
virInit2	<-c(test[1:(maxAge-1)],sum(test[maxAge:100]))
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
 
calcRemoval<-function( calcREFouts ) 
{
#==this should have options for how the control rule will change
 }
