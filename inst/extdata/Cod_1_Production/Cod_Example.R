library(GeMS)
MSEdir<-system.file("extdata/Cod_1_Production",package="GeMS")

OMNames<-c("Cod_LowProd_CTL","Cod_Base_CTL","Cod_HighProd_CTL")

run_GeMS(MSEdir=MSEdir,
         CTLNameList=OMNames)

##-----------------
## Parallel example
##-----------------
#run_GeMS(MSEdir=MSEdir,
#         CTLNameList=OMNames,
#         runparallel = T, cores = 3)
