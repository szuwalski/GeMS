library(GeMS)
Cur.dir<-system.file("extdata/Cod_1_Production",package="GeMS")

OMNames<-c("Cod_LowProd_CTL","Cod_Base_CTL","Cod_HighProd_CTL")

run_GeMS(MSEdir=Cur.dir,
         CTLNameList=OMNames)

##-----------------
## Parallel example
##-----------------
#run_GeMS(MSEdir=Cur.dir,
#         CTLNameList=OMNames,
#         runparallel = T, cores = 3)
