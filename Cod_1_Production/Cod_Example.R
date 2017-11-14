rm(list=ls())
source("C:/GeMS/General_MSE_helpers.R")
source("C:/GeMS/General_MSE_main.R")

CurDir<-"C:/GeMS/Cod_1_Production/"
setwd(CurDir)
CreateFolderNameList<-c("Cod_Base_CTL","Cod_LowProd_CTL","Cod_veryLowProd_CTL")
CreateFolderNameList<-c("Cod_veryLowProd_CTL")

#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
 Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],".csv",sep=""))
 GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
 ProductionModelOutput(CreateFolderNameList[x])
}


