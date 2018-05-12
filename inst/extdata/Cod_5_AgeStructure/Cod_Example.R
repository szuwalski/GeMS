library(GeMS)
#==identify the directory where current analysis will be performed and CTL files
Cur.dir<-system.file("extdata/Cod_5_AgeStructure",package="GeMS")

OMNames<-c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")

#==Loop that reads in CTL files stored in dir.MSE and executes code in sequence
run_GeMS(CTLNameList=OMNames,
         MSEdir=Cur.dir)

#==Loop that reads in CTL files stored in dir.MSE and executes code in parallel
# run_GeMS(CTLNameList=OMNames,
#         MSEdir=Cur.dir,
#         runparallel=T, cores=3)

##==compare just two scenarios: estimating M vs. not estimating M when M is varying
#OMNames<-c("Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")
# Inout<-ReadCTLfile(file.path(Cur.dir,OMNames[1]))
# AgeStructureComp(out=Inout,CTLNameList=OMNames,MSEdir=Cur.dir)
