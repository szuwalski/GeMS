rm(list=ls())

#==name where the GeMS folder lives on your computer
GeMS.dir <- "C:/GeMS"
#==identify the directory where current analysis will be performed and CTL files
dir.MSE<-"C:/GeMS/Cod_5_AgeStructure"

#==source important scripts
source(file.path(GeMS.dir,"run_GeMS.R"))

OMNames<-c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")

#==Loop that reads in CTL files stored in dir.MSE and executes code
run_GeMS(GeMSDir=GeMS.dir, CurDir=dir.MSE,
		 CreateFolderNameList=OMNames)

##==compare just two scenarios: estimating M vs. not estimating M when M is varying
##==WARNING: Overwrites plots made from run_GeMS
#CreateFolderNameList<-c("Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")
#AgeStructureComp(Inout=Inout,CreateFolderNameList=CreateFolderNameList,MSEdir=dir.MSE)
