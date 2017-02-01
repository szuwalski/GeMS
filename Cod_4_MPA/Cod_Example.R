rm(list=ls())
source("C:/GeMS/General_MSE_helpers.R")
source("C:/GeMS/General_MSE_main.R")

CurDir<-"C:/GeMS/Cod_4_MPA/"
setwd(CurDir)
CreateFolderNameList<-c("Cod_MPA_CTL")

#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
  Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],".csv",sep=""))
  GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
  ProductionModelOutput(CreateFolderNameList=CreateFolderNameList[x])
}

