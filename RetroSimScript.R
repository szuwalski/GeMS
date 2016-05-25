#################################################
#=='Fixing' retrospective biases in stock assessment
#==Szuwalski et al.
##################################################
#set.seed(001)


#==MAKE COMPARISON PLOTS OF ALL THE PARAMETERS

#==this is where the output and files will be stored
CurDir<-"C:/Users/cody_admin/Desktop/Retro/"
setwd(CurDir)
source("General_MSE_helpers.R")
source("General_MSE_main.R")

#==Names of the CTL files (which will also be the names of the folders for output)
#CreateFolderNameList<-c("NoEstBase01","EstMbase01","EstGrowthBase01","EstSelBase01")
#CreateFolderNameList<-c("NoEstSelvary","EstMselvary","EstGrowthSelvary","EstSelSelvary")
#CreateFolderNameList<-c("NoEstGrowvary","EstMgrowvary","EstGrowthGrowvary","EstSelGrowvary")
CreateFolderNameList<-c("NoEstMvary","EstMmvary","EstGrowthMvary","EstSelMvary")
#CreateFolderNameList<-c("NoEstDet","EstMDet","EstGrowthDet","EstSelDet")
#CreateFolderNameList<-c("EstGrowthDet")

#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
 Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],"_CTL.csv",sep=""))
 GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
}

#Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[1],"_CTL.csv",sep=""))
#==Save gradients
RetroPeels<-6
TakeRows<-(Inout$OM$SimYear-RetroPeels+1):Inout$OM$SimYear
GradientSave<-array(dim=c(RetroPeels,Inout$OM$Nsim,length(CreateFolderNameList)))

for(x in 1:length(CreateFolderNameList))
{
 Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[1],"_CTL.csv",sep=""))
 GradientSave[,,x]<-CheckGradient(DrawDir=CreateFolderNameList[x],out=Inout)[TakeRows,]
}
ScenCols<-c("grey",2,3,4)
#==reference points
TakeRows<-(Inout$OM$InitYear+1):Inout$OM$SimYear
B35save<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CreateFolderNameList)))
F35save<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CreateFolderNameList)))
OFLsave<-array(dim=c(length(TakeRows),Inout$OM$Nsim,length(CreateFolderNameList)))

for(x in 1:length(CreateFolderNameList))
{
 temp<-CheckReferencePoints(DrawDir=CreateFolderNameList[x],out=Inout)
 B35save[,,x]<-temp[[1]][TakeRows,]
 F35save[,,x]<-temp[[2]][TakeRows,]
 OFLsave[,,x]<-temp[[3]][TakeRows,]
} 


BigB35<-matrix(nrow=Inout$OM$Nsim,ncol=length(CreateFolderNameList))
BigF35<-matrix(nrow=Inout$OM$Nsim,ncol=length(CreateFolderNameList))
BigOFL<-matrix(nrow=Inout$OM$Nsim,ncol=length(CreateFolderNameList))
for(x in 1:length(CreateFolderNameList))
{
 BigB35[,x]<-apply(B35save[,,x],2,median) 
 BigF35[,x]<-apply(F35save[,,x],2,median) 
 BigOFL[,x]<-apply(OFLsave[,,x],2,median) 
}
dev.new()
par(mfrow=c(length(CreateFolderNameList),3),mar=c(.1,2,.1,.1),oma=c(3,6,3,4))
for(x in 1:length(CreateFolderNameList))
{
 #boxplot(c(MohnsRho[,,x]),col=x)
par(mar=c(.1,2,.1,.1))
 boxplot(B35save[,,x],col=ScenCols[x],xaxt='n',las=1,ylim=c(quantile(B35save,0.1),quantile(B35save,0.9)))
 abline(h=median(B35save),lty=2)
if(x==1)
mtext(side=3,"Target biomass")

par(mar=c(.1,4,.1,.1))
 boxplot(F35save[,,x],col=ScenCols[x],xaxt='n',las=1,,ylim=c(quantile(F35save,0.1),quantile(F35save,0.9)))
 #,ylim=c(0.09,0.2))
 abline(h=median(F35save),lty=2)
if(x==1)
mtext(side=3,"Target fishing mortality")

par(mar=c(.1,.1,.1,.1))
 boxplot(OFLsave[,,x],col=ScenCols[x],xaxt='n',yaxt='n',ylim=c(quantile(OFLsave,0.1),quantile(OFLsave,0.9)))
 axis(side=4,las=1)
 abline(h=median(OFLsave),lty=2)
if(x==1)
mtext(side=3,"OFL")
}
mtext(outer=T,side=2,"No time varying",adj=.9,line=3,cex=.7)
mtext(outer=T,side=2,"M time varying",adj=.65,line=3,cex=.7)
mtext(outer=T,side=2,"Growth time varying",adj=.35,line=3,cex=.7)
mtext(outer=T,side=2,"Sel time varying",adj=.1,line=3,cex=.7)


#=plotting output
# GeMSplots(out=Inout,SourceFolder=CreateFolderNameList[x])

#==
MohnsRho<-array(dim=c(RetroPeels-1,Inout$OM$Nsim,length(CreateFolderNameList)))
SSBbias<-array(dim=c(RetroPeels,Inout$OM$Nsim,length(CreateFolderNameList)))
for(x in 1:length(CreateFolderNameList))
{
Inout			<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],"_CTL.csv",sep=""))
DrawDir		<-paste(CreateFolderNameList[x],sep="")
RetroOuts		<-CheckRetro(RetroPeels=RetroPeels,DrawDir,PlotRetro=0,out=Inout)
MohnsRho[,,x]	<-CalculateMohn(RetroOuts)
SSBbias[,,x]	<-CalculateBias(RetroOuts)
}
dev.new()
par(mfrow=c(length(CreateFolderNameList),1),mar=c(.1,.1,.1,.1),oma=c(3,3,3,3))
for(x in 1:length(CreateFolderNameList))
{
 #boxplot(c(MohnsRho[,,x]),col=x)
 boxplot(MohnsRho[,,x],col=x)
 abline(h=0,lty=2)
}

dev.new()
par(mfrow=c(length(CreateFolderNameList),1),mar=c(.1,.1,.1,.1),oma=c(3,3,3,3))
for(x in 1:length(CreateFolderNameList))
{
 #boxplot(c(MohnsRho[,,x]),col=x)
 boxplot(SSBbias[,,x],col=ScenCols[x])
 abline(h=0,lty=2)
}

BigMohn<-matrix(nrow=Inout$OM$Nsim,ncol=length(CreateFolderNameList))
for(x in 1:length(CreateFolderNameList))
{
 temp<-MohnsRho[,,x]
 BigMohn[,x]<-apply(temp,2,mean) 
}

BigBias<-matrix(nrow=Inout$OM$Nsim,ncol=length(CreateFolderNameList))
for(x in 1:length(CreateFolderNameList))
{
 temp<-SSBbias[,,x]
 BigBias[,x]<-apply(temp,2,mean) 
}

dev.new()
par(mfrow=c(5,1),mar=c(.1,.1,.1,.1),oma=c(4,6,1,1))
boxplot(BigMohn,col=ScenCols,xaxt='n',las=1)
mtext(side=2,"Retrospective bias",line=4,cex=.7)
abline(h=0,lty=2)

boxplot(BigBias,col=ScenCols,xaxt='n',las=1)
mtext(side=2,"Bias",line=4,cex=.7)
abline(h=0,lty=2)

boxplot(BigB35,col=ScenCols,ylim=c(quantile(BigB35,0.1),quantile(BigB35,0.9)),xaxt='n',las=1)
mtext(side=2,"Estimated B35%",line=4.2,cex=.7)
boxplot(BigF35,col=ScenCols,ylim=c(quantile(BigF35,0.1),quantile(BigF35,0.9)),xaxt='n',las=1)
mtext(side=2,"Estimated F35%",line=4,cex=.7)
boxplot(BigOFL,col=ScenCols,ylim=c(quantile(BigOFL,0.1),quantile(BigOFL,0.9)),las=1,names=c("No vary","M vary","Growth vary","Sel vary"))
mtext(side=2,"Estimated OFL",line=4,cex=.7)

#=====compare estimated vs. time varying processes====
quants<-PullTimevary(inFolders=CreateFolderNameList,out=Inout)
dev.new(width=7,height=4)
inmat<-matrix(c(1,1,2,2,3,3,
		    4,4,4,5,5,5),nrow=2,byrow=T)
layout(inmat)
par(mar=c(.1,.1,.1,.1),oma=c(4,4,4,4))

#==RECRUITMENT
input<-quants[[1]]
tInput<-quants[[2]] 

plot(-100000,xlim=c(1,dim(input)[2]),ylim=c(0,max(input,na.rm=T)),las=1,xaxt='n')
for(x in 1:length(CreateFolderNameList))
{
 color<-seq(1,length(CreateFolderNameList)+1)
 incol<-adjustcolor(color,alpha.f=.08)
 tCol<-rgb(0,0,0,0.2)
 temp<-input[,,x]
 for(y in 1:dim(input)[1])
   lines(input[y,,x],col=incol[x])
} 
medians<- apply(input,c(2,3),median,na.rm=T)
for(x in 1:ncol(medians))
 lines(medians[,x],col=x,cex=2)
#==true R
plotIn<-apply(tInput,2,median)
plotIn[1]<-NA
lines(plotIn,col="yellow",lwd=2)

legend("bottomleft",bty='n',"Recruitment")

#==FISHING MORTALITY
input<-quants[[3]]
NatM<-quants[[7]]
tInput<-quants[[4]] 

plot(-1000,xlim=c(1,dim(input)[2]),ylim=c(0,max(input,NatM,na.rm=T)),las=1,xaxt='n',yaxt='n')
axis(side=3)
mtext(side=3,line=2,"Time")
for(x in 1:length(CreateFolderNameList))
{
 color<-seq(1,length(CreateFolderNameList)+1)
 incol<-adjustcolor(color,alpha.f=.08)
 tCol<-rgb(0,0,0,0.2)
 temp<-input[,,x]
 for(y in 1:dim(input)[1])
   lines(input[y,,x],col=incol[x])
} 
medians<- apply(input,c(2,3),median,na.rm=T)
for(x in 1:ncol(medians))
 lines(medians[,x],col=x,cex=2)
#==true F
dim(tInput)
lines(apply(tInput,2,median),col="yellow",lwd=2)
legend("bottomleft",bty='n',"Fishing mortality")

#==NATURAL MORTALITY
input<-quants[[7]]
tInput<-quants[[8]] 
dim(input)
plot(-1000,xlim=c(1,dim(input)[2]),ylim=c(0,max(input,NatM,na.rm=T)),las=1,xaxt='n',yaxt='n')
axis(side=4,las=1)
for(x in 1:length(CreateFolderNameList))
{
 color<-seq(1,length(CreateFolderNameList)+1)
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
lines(tInput,col="yellow",lwd=2)

legend("bottomleft",bty='n',"Natural mortality")

#==SELECTIVITY
input<-quants[[9]]
tInput<-quants[[10]] 
plot(-1000,xlim=c(1,which(input[1,,1,1]>0.99)[3]),ylim=c(0,max(input,na.rm=T)),las=1)
for(x in 1:length(CreateFolderNameList))
{
 color<-seq(1,length(CreateFolderNameList)+1)
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
legend("bottomright",bty='n',"Fishery selectivity")

#==GROWTH
input<-quants[[5]]
tInput<-quants[[6]] 
plot(-1000,xlim=c(1,MaxAge),ylim=c(0,max(input,na.rm=T)),las=1,yaxt='n')
axis(side=4,las=1)
for(x in 1:length(CreateFolderNameList))
{
 color<-seq(1,length(CreateFolderNameList)+1)
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

legend("bottomright",bty='n',"Length at age")
mtext(side=1,outer=T,"Age",line=2.3)
