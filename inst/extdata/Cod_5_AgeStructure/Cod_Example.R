library(GeMS)
#==identify the directory where current analysis will be performed and CTL files
MSEdir<-system.file("extdata/Cod_5_AgeStructure",package="GeMS")

CTLNameList<-c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")

#==Loop that reads in CTL files stored in dir.MSE and executes code in sequence
run_GeMS(CTLNameList=CTLNameList,
         MSEdir=MSEdir)

#==Loop that reads in CTL files stored in dir.MSE and executes code in parallel
# run_GeMS(CTLNameList=CTLNameList,
#          MSEdir=MSEdir,
#          runparallel=T, cores=3)

##==compare just two scenarios: estimating M vs. not estimating M when M is varying
# CTLNameList<-c("Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")
# Inout<-ReadCTLfile(file.path(MSEdir,CTLNameList[1]))
# AgeStructureComp(out=Inout,CTLNameList=CTLNameList,MSEdir=MSEdir)
