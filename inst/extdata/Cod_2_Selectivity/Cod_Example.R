library(GeMS)
Cur.dir<-system.file("Cod_2_Selectivity",package="GeMS")

OMNames<-c("Cod_Selectivity_CTL")

run_GeMS(CTLNameList=OMNames,
		 MSEdir=Cur.dir)


