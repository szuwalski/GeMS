rm(list=ls())

#==name where the GeMS folder lives on your computer
TopDir<-"C:/GeMS/"

#==source important scripts
source(paste(TopDir,"General_MSE_helpers.R",sep=""))
source(paste(TopDir,"General_MSE_main.R",sep=""))

#==identify the directory where current analysis will be performed and CTL files
CurDir<-"C:/GeMS/Cod_5_AgeStructure/"
setwd(CurDir)
CreateFolderNameList<-c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")

#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
  Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],".csv",sep=""))
  GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
}

#==compare estimating M vs. not estimating M when M is varying
CreateFolderNameList<-c("Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")
AgeStructureComp(Inout=Inout,CreateFolderNameList=CreateFolderNameList)
