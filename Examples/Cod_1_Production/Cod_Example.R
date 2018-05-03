rm(list=ls())
library(GeMS)
Cur.dir<-"C:/GeMS/Examples/Cod_1_Production/"

OMNames<-c("Cod_LowProd_CTL","Cod_Base_CTL","Cod_HighProd_CTL")

run_GeMS(MSEdir=Cur.dir,
         CTLNameList=OMNames)

##-----------------
## Parallel example
##-----------------
#run_GeMS(MSEdir=Cur.dir,
#         CTLNameList=OMNames,
#         runparallel = T, cores = 3)
