library(GeMS)

#==Assumes your working directory is set to where the control file is located.

CTLNameList<-c("Cod_Base","Cod_F35","Cod_F35_Slope","Cod_LenComp","Cod_EM","Cod_EM_LenComp","Cod_EM_smooth")

#==Loop that reads in CTL files stored in MSEdir and executes code in sequence
run_GeMS(CTLNameList=CTLNameList[7])

#==Loop that reads in CTL files stored in MSEdir and executes code in parallel
 run_GeMS(CTLNameList=CTLNameList,
          runparallel=T, cores=7)

##==compare just two scenarios: estimating M vs. not estimating M when M is varying
 Inout<-ReadCTLfile(CTLNameList[1])
 AgeStructureComp(out=Inout,CTLNameList=CTLNameList[1:3],MSEdir=getwd(),Nruns=20,
 					plotNames=c("Constant Catch","F35","Sloped F35"))


##==compare just two scenarios: estimating M vs. not estimating M when M is varying
 Inout<-out<-ReadCTLfile(CTLNameList[2])
 AgeStructureComp(out=Inout,CTLNameList=CTLNameList[c(2,5,6,7,8)],MSEdir=getwd(),Nruns=20,
 					plotNames=c("Fixed L50","Vary L50","Better Fut. Data","Better Hist. Data","Better Hist. Data and Smoother"))
 
 CTLName<-"Cod_EM_SuperLenComp_smooth"
 AgeStructureComp(out=ReadCTLfile(CTLName),CTLNameList=CTLName,MSEdir=getwd())

