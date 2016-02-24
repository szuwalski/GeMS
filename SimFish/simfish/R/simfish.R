#' \code{simfish} simulates codys fishery
#'
#' @param Nsim number of simulations
#' @param SimYear number of years to run simulation
#' @param depletion target depletion
#' @param CatchShareN share of the catch in the north parth
#' @param CalcFMSY 1 or 0 to calculate Fmsy doesn't work now
#' @param SmallNum constant to add to zeros
#' @param InitSmooth no idea
#' @param FmortPen no idea
#' @param RecruitPen no idea
#' @param CatchCVn CV of the catch in the north
#' @param CatchCVs CV of the catch in the south
#' @param IndexCVn CV of the survey index in the north
#' @param IndexCVs CV of the survey index in the south
#' @param LenSampleN number of lengths samples in the north
#' @param LenSampleS number of length samples in the south
#' @param GrowthSDn growth standard deviation in the north
#' @param GrowthSDs growth standard deviation in the south
#' @param surv50n lower survey selectivity in the north
#' @param surv95n upper survey selectivity in the north
#' @param surv50s lower survey selectivity in the south
#' @param surv95s upper survey selectivity in the south
#' @param MaxAge maximum age of the species
#' @param NatMn natural mortality in the north
#' @param NatMs natural mortality in the south
#' @param VonKn von bert k in the north
#' @param VonKs von bert k in the south
#' @param LinfN linfinity in the north
#' @param LinfS linfinity in the south
#' @param t0n t0 in the north
#' @param t0s t0 in the south
#' @param mat50n lower maturity in the north
#' @param mat50s lower maturity in the south
#' @param mat95n upper maturity in the north
#' @param mat95s upper maturity in the north
#' @param alphaN alpha of weight function in the north
#' @param betaN beta of weight function in the north
#' @param alphaS alpha of weight function in the south
#' @param betaS beta of weight function in the south
#' @param MaxMovingN maximum moving numbers I think in the north
#' @param MaxMovingS maximum moving numbers I think in the south
#' @param Move50n lower movement ogive north
#' @param Move50s lower movement ogive south
#' @param Move95n upper movement ogive north
#' @param Move95s upper movement ogive south
#' @param steepnessN recruitment steepness in the north
#' @param steepnessS recruitment steepness in the south
#' @param sigmaRn recruitment deviations in the north
#' @param sigmaRs recruitment deviations in the south
#' @param RzeroN eq recruitment in the north
#' @param RzeroS eq recruitment in the south
#' @param sel50n lower fishing selectivity ogive in the north
#' @param sel50s lower fishing selectivity ogive in the south
#' @param sel95n upper fishing selectivity ogive in the north
#' @param sel95s upper fishing selectivity ogive in the south
#' @param HarvestControl no idea
#' @param HistoricalF historical pattern of fishing mortality
#' @param ConstantCatch no idea
#' @param ConstantF no idea
#' @param HCalpha no idea
#' @param HCbeta mo idea
#'
#' @return list of real and observed data and plots
#' @export
simfish <- function(Nsim = 1, SimYear = 50, depletion = 0, CatchShareN = 1, CalcFMSY = 0,
                                    SmallNum = 1e-4, InitSmooth = 3, FmortPen = 3, RecruitPen = 3,
                    CatchCVn = 0.01, CatchCVs = 0.01, IndexCVn = 0.05, IndexCVs = 0.05,
                                    LenSampleN= 200, LenSampleS = 200, GrowthSDn = 20, GrowthSDs = 20,
                                    surv50n = 90, surv95n =  150, surv50s = 90, surv95s = 150,
                    MaxAge = 30, NatMn = 0.2, NatMs = 0.2, VonKn = 0.2, VonKs = 0.2, LinfN = 300,
                                 LinfS = 300, t0n = 0.9, t0s = 0.9, mat50n = 6, mat50s = 6, mat95n = 8, mat95s = 8,
                                 alphaN = 1.7e-06, betaN = 3, alphaS = 1.7e-06, betaS = 3,
                                 MaxMovingN = 0, MaxMovingS = 0, Move50n = 6, Move50s = 6, Move95n = 10, Move95s = 10,
                                 steepnessN = 0.7, steepnessS = 0.7, sigmaRn = 0.001, sigmaRs =  0.001, RzeroN = 1e5,
                                 RzeroS = 1e5,
                    sel50n = 190, sel50s = 190, sel95n = 225, sel95s = 225, HarvestControl = 3,
                                 HistoricalF = 0.3, ConstantCatch = 0, ConstantF = 0.1, HCalpha = 0.05, HCbeta = 0.5 )
  {

  #=======================
  #==simulation controls==
  #=======================

  plot_theme <- theme_economist() + theme(text = element_text(size = 12))



  Nsim <- Nsim   #out$OM$Nsim			# number of simulations to do in the MSE
  SimYear <- SimYear #  out$OM$SimYear			# total number of years in simulation
  InitYear <- SimYear #out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
  depletion	 <- depletion		# model is conditioned to a specific depletion given a trajectory of effort
  CatchShareN	 <- CatchShareN		# the amount of effort allocated to a given area
  CatchShareS	 <- 1-CatchShareN
  SmallNum <- SmallNum
  InitSmooth <- InitSmooth
  FmortPen <- FmortPen
  RecruitPen <- RecruitPen
  CalcFMSY <- CalcFMSY
  #==========================================================
  #=================population dynamics======================
  #=========================================================
  #==When making dynamics time-varying, only make them vary after
  #=='InitYear'.  The idea is that this simulates a relatively steady
  #==environment until climate change, then things change. MEH.
  #===========================================================
  #===========================================================
  MaxAge <- MaxAge

  Ages <- seq(1,MaxAge)

  #==natural mortality================
  NatMn	<-  CleanInput(NatMn, SimYear)

  NatMs		<-  CleanInput(NatMs, SimYear)

  #==Length at age====================
  VonKn		<-  CleanInput(VonKn, SimYear)
  LinfN		<-  CleanInput(LinfN,SimYear)
  t0n		<-  CleanInput(t0n, SimYear)
  LenAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    LenAtAgeN[i,]<-LinfN[i]*(1-exp(-VonKn[i]*(Ages-t0n[i])))

  VonKs		<-  CleanInput(VonKs, SimYear)
  LinfS		<-  CleanInput(LinfS, SimYear)
  t0s		<-  CleanInput(t0s, SimYear)
  LenAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear){
    LenAtAgeS[i,]<- LinfS[i]*(1-exp(-VonKs[i]*(Ages-t0s[i])))
  }

  #==specify the number of length bins
  BinWidth<-(max(LenAtAgeS,LenAtAgeN)*1.05)/MaxAge
  LengthBins<-seq(0,max(LenAtAgeS,LenAtAgeN)*1.05,BinWidth)
  LengthBinsMid<-LengthBins[1:(length(LengthBins)-1)] + mean(LengthBins[1:2])

  #==maturity at age==========================
  mat50n	<- CleanInput(mat50n, SimYear)
  mat95n	<- CleanInput(mat95n, SimYear)
  matureN	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    matureN[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50n[i])/(mat95n[i]-mat50n[i])))

  mat50s	<- CleanInput(mat50s, SimYear)
  mat95s	<-  CleanInput(mat95s,SimYear)
  matureS	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    matureS[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50s[i])/(mat95s[i]-mat50s[i])))

  #==weight at age==========================
  alphaN		<-   CleanInput(alphaN, SimYear)
  betaN			<-  CleanInput(betaN, SimYear)
  WeightAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    WeightAtAgeN[i,] <- alphaN[i]*LenAtAgeN[i,]^betaN[i]

  alphaS <- CleanInput(alphaS, SimYear)
  betaS <-  CleanInput(betaS, SimYear)
  WeightAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    WeightAtAgeS[i,]	<-alphaS[i]*LenAtAgeS[i,]^betaS[i]

  #==movement from box to box==
  MaxMovingN	<- CleanInput(MaxMovingN, SimYear)
  Move50n	<- CleanInput(Move50n, SimYear)
  Move95n	<- CleanInput(Move50n, SimYear)
  MovementN	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    MovementN[i,]	<-MaxMovingN[i]/(1+exp(-1*log(19)*(Ages-Move50n[i])/(Move95n[i]-Move50n[i])))

  MaxMovingS <- CleanInput(MaxMovingS,SimYear)
  Move50s	<- CleanInput(Move50s, SimYear)
  Move95s	<- CleanInput(Move95s, SimYear)
  MovementS	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    MovementS[i,]	<-MaxMovingS[i]/(1+exp(-1*log(19)*(Ages-Move50s[i])/(Move95s[i]-Move50s[i])))

  #==recruitment parameters==
  steepnessN	<-  CleanInput(steepnessN, SimYear)
  sigmaRn	<- CleanInput(sigmaRn, SimYear)
  RzeroN	<- CleanInput(RzeroN, SimYear)

  steepnessS	<-  CleanInput(steepnessS, SimYear)
  sigmaRs	<-  CleanInput(sigmaRs, SimYear)
  RzeroS	<- CleanInput(RzeroS ,SimYear)

  RecErrN	<-matrix(rnorm(SimYear*Nsim,0,sigmaRn[1]),ncol=SimYear)
  RecErrS	<-matrix(rnorm(SimYear*Nsim,0,sigmaRs[1]),ncol=SimYear)

 # Fishing Fleet ----

  sel50n	<-  CleanInput(sel50n, SimYear)
  sel95n	<-  CleanInput(sel95n, SimYear)
  vulnN		<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    vulnN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50n[i])/(sel95n[i]-sel50n[i])))

  vulnN[vulnN<0.01]<-0

  sel50s	<-  CleanInput(sel50s, SimYear)
  sel95s	<- CleanInput(sel95s, SimYear)
  vulnS		<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    vulnS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50s[i])/(sel95s[i]-sel50s[i])))

  vulnS[vulnS<0.01]<-0

  #==historic fishing mortality
  HistoricalF		<-   CleanInput(HistoricalF, SimYear)[1:InitYear]
  HistoricalF[is.na(HistoricalF)] <- mean(HistoricalF, na.rm = T)
  HarvestControl	<- HarvestControl #out$OM$HarvestControl
  ConstantCatch	<- ConstantCatch #out$OM$ConstantCatch
  ConstantF		<- ConstantF #out$OM$ConstantF
  HCalpha		<- HCalpha #out$OM$HCalpha
  HCbeta		<- HCbeta #out$OM$HCbeta

  # Surveys -----------------------------------------------------------------

  #==sampling uncertainty
  CatchCVn	<- CatchCVn
  CatchCVs	<- CatchCVs

  IndexCVn	<- IndexCVn
  IndexCVs	<- IndexCVs

  LenSampleN	<- LenSampleN
  LenSampleS	<- LenSampleS

  GrowthSDn	<- GrowthSDn
  GrowthSDs	<- GrowthSDs


  #==index selectivity=====================
  surv50n	<- CleanInput(surv50n, SimYear)
  surv95n	<-  CleanInput(surv95n, SimYear)
  survSelN		<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    survSelN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50n[i])/(surv95n[i]-surv50n[i])))


  surv50s	<-  CleanInput(surv50s,SimYear)
  surv95s	<-  CleanInput(surv95s,SimYear)
  survSelS	<-matrix(nrow=SimYear,ncol=MaxAge)
  for(i in 1:SimYear)
    survSelS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50s[i])/(surv95s[i]-surv50s[i])))


  #===================================================
  #==Virgin numbers at age, biomass, initial depletion, recruitment
  #===================================================
  VirInitN<-initialN(Rzero=RzeroN[1],NatM=NatMn[1],inAge=MaxAge)
  VirInitS<-initialN(Rzero=RzeroS[1],NatM=NatMs[1],inAge=MaxAge)

  VirBioN<- sum(VirInitN*matureN[1,]*WeightAtAgeN[1,])
  VirBioS<-sum(VirInitS*matureS[1,]*WeightAtAgeS[1,])

  ExploitBioN<-sum(VirInitN*vulnN[1,]*WeightAtAgeN[1,])
  ExploitBioS<-sum(VirInitS*vulnS[1,]*WeightAtAgeS[1,])

  #========================================================================
  # FIND MSY FOR THE POPULATION (should just make a popdym function, dummy)=
  #========================================================================

  if(CalcFMSY==1)
  {
    SearchFmort		<-seq(0.01,3*NatMn[1],(NatMn[1]-0.01)/100)
    SearchYield		<-rep(0,length(SearchFmort))
    SearchBiomass	<-rep(0,length(SearchFmort))
    for(p in 1:length(SearchFmort))
    {
      tempNn		<-matrix(ncol=MaxAge,nrow=100)
      tempNn[1,]		<-VirInitN
      tempNs		<-matrix(ncol=MaxAge,nrow=100)
      tempNs[1,]		<-VirInitS
      tempCatchN		<-rep(0,100)
      tempCatchS		<-rep(0,100)
      tempRecN		<-rep(0,100)
      tempRecS		<-rep(0,100)
      tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=100)
      tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=100)

      inFs<-SearchFmort[p]*CatchShareS
      inFn<-SearchFmort[p]*CatchShareN

      for (j in 2:100)
      {
        for (i in 2:(MaxAge-1))
        {
          tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-inFn*vulnN[1,i-1])*exp(-NatMn[1])
          tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-inFs*vulnS[1,i-1])*exp(-NatMs[1])

        }
        tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-inFn*vulnN[1,MaxAge])*exp(-NatMn[1])+ tempNn[j-1,MaxAge]*exp(-inFn*vulnN[1,MaxAge])*exp(-NatMn[1])
        tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-inFs*vulnS[1,MaxAge])*exp(-NatMs[1])+ tempNs[j-1,MaxAge]*exp(-inFs*vulnS[1,MaxAge])*exp(-NatMs[1])

        EggsN			<-sum(tempNn[j-1,]*matureN[1,]*WeightAtAgeN[1,])
        EggsS			<-sum(tempNs[j-1,]*matureS[1,]*WeightAtAgeS[1,])
        tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1,j],recType="BH",NatMin=NatMn[1],
                                   vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,], MaxAge = MaxAge)
        tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1,j],recType="BH",NatMin=NatMs[1],
                                   vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,], MaxAge = MaxAge)
        tempRecN[j]		<-tempNn[j,1]
        tempRecS[j]		<-tempNs[j,1]

        tempCatchAtAgeN[j,]	<-((vulnN[1,]*inFn)/(vulnN[1,]*inFn+NatMn[1])) * (1-exp(-(vulnN[1,]*inFn+NatMn[1]))) * tempNn[j-1,]
        tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[1,])

        tempCatchAtAgeS[j,]	<-((vulnS[1,]*inFs)/(vulnS[1,]*inFs+NatMs[1])) * (1-exp(-(vulnS[1,]*inFs+NatMs[1]))) * tempNs[j-1,]
        tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[1,]) #catch in weight

      }
      SearchYield[p]	<-tempCatchS[j]+tempCatchN[j]
      SearchBiomass[p]	<-EggsN+EggsS
    }

    trueFMSY	<-SearchFmort[which.max(SearchYield)]
    trueUMSY	<- 1-(exp(-trueFMSY))
    trueBMSY	<-SearchBiomass[which.max(SearchYield)]
    trueMSY	<-max(SearchYield)

    dev.new()
    par(mfrow=c(1,2))
    plot(SearchYield~SearchFmort)
    plot(SearchYield~SearchBiomass)
  }
  #=========================================================================
  # INITALIZE THE POPULATION
  # Option 1: Set the initial conditions by scaling the specified F series
  # that will put the population at the designated depletion
  #
  # Option 2: Set depletion to 0 and use the input time series of F to
  # initialize the population
  #==========================================================================
  if(depletion>0)
  {
    targetDep	<-(VirBioN+VirBioS)*depletion
    maxExp	<-0.999
    minExp	<-0.001

    #==defines the exploitation pattern
    ExpSeries	<-HistoricalF

    #==finds an exploitation rate that returns a given depletion at the end of the
    #==INCORPORATE THE CATCHSHARE BY AREA INTO THIS CHUNK OF CODE============
    for(x in 1:30)
    {
      tempNn	<-matrix(ncol=MaxAge,nrow=InitYear)
      tempNn[1,]	<-VirInitN

      tempNs	<-matrix(ncol=MaxAge,nrow=InitYear)
      tempNs[1,]	<-VirInitS

      tempCatchN	<-rep(0,InitYear)
      tempCatchS	<-rep(0,InitYear)

      tempRecN	<-rep(0,InitYear)
      tempRecS	<-rep(0,InitYear)

      tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=InitYear)
      tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=InitYear)

      tempExp	<-(maxExp+minExp)/2
      tempF		<- -log(1-tempExp)
      tempFsrsN	<-ExpSeries*tempF
      tempFsrsS	<-ExpSeries*tempF

      for (j in 2:InitYear)
      {
        for (i in 2:(MaxAge-1))
        {
          tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-tempFsrsN[j]*vulnN[1,i-1])*exp(-NatMn[1])
          tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-tempFsrsS[j]*vulnS[1,i-1])*exp(-NatMs[1])

        }
        tempNn[j,MaxAge]	<-(tempNn[j-1,(MaxAge-1)])*exp(-tempFsrsN[j]*vulnN[1,MaxAge])*exp(-NatMn[1])+ tempNn[j-1,MaxAge]*exp(-tempFsrsN[j]*vulnN[1,MaxAge])*exp(-NatMn[1])
        tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-tempFsrsS[j]*vulnS[1,MaxAge])*exp(-NatMs[1])+ tempNs[j-1,MaxAge]*exp(-tempFsrsS[j]*vulnS[1,MaxAge])*exp(-NatMs[1])

        EggsN			<-sum(tempNn[j-1,]*matureN[1,]*WeightAtAgeN[1,])
        EggsS			<-sum(tempNs[j-1,]*matureS[1,]*WeightAtAgeS[1,])

        tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1,j],recType="BH",NatMin=NatMn[1],
                                   vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,])
        tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1,j],recType="BH",NatMin=NatMs[1],
                                   vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,])
        tempRecN[j]		<-tempNn[j,1]
        tempRecS[j]		<-tempNs[j,1]

        tempCatchAtAgeN[j,]	<-((vulnN[1,]*tempFsrsN[j])/(vulnN[1,]*tempFsrsN[j]+NatMn[1])) * (1-exp(-(vulnN[1,]*tempFsrsN[j]+NatMn[1]))) * tempNn[j-1,]
        tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[1,])

        tempCatchAtAgeS[j,]	<-((vulnS[1,]*tempFsrsS[j])/(vulnS[1,]*tempFsrsS[j]+NatMs[1])) * (1-exp(-(vulnS[1,]*tempFsrsS[j]+NatMs[1]))) * tempNs[j-1,]
        tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[1,])

      }

      tempDepN	<-sum(tempNn[InitYear,]*matureN[1,]*WeightAtAgeN[1,])
      tempDepS	<-sum(tempNs[InitYear,]*matureS[1,]*WeightAtAgeS[1,])
      tempDep	<-tempDepN+tempDepS

      if(tempDep>targetDep)
        minExp	<-tempExp
      if(tempDep<targetDep)
        maxExp	<-tempExp
      initExp	<-tempExp
    }
    HistoricalF<-ExpSeries*tempExp
    print(c("Final Depletion",tempDep))
  }

  if(depletion==0)
  {
    tempNn	<-matrix(ncol=MaxAge,nrow=InitYear)
    tempNn[1,]	<-VirInitN
    tempNs	<-matrix(ncol=MaxAge,nrow=InitYear)
    tempNs[1,]	<-VirInitS
    tempCatchN	<-rep(0,InitYear)
    tempCatchS	<-rep(0,InitYear)
    tempRecN	<-rep(0,InitYear)
    tempRecS	<-rep(0,InitYear)
    tempCatchAtAgeN	<-matrix(ncol=MaxAge,nrow=InitYear)
    tempCatchAtAgeS	<-matrix(ncol=MaxAge,nrow=InitYear)

    HistoricalFs<-HistoricalF
    HistoricalFn<-HistoricalF

    for (j in 2:InitYear)
    {
      for (i in 2:(MaxAge-1))
      {
        tempNn[j,i]		<-tempNn[j-1,i-1]*exp(-HistoricalFn[j]*vulnN[1,i-1])*exp(-NatMn[j])
        tempNs[j,i]		<-tempNs[j-1,i-1]*exp(-HistoricalFs[j]*vulnS[1,i-1])*exp(-NatMs[j])

      }
      tempNn[j,MaxAge]	<- (tempNn[j-1,(MaxAge-1)])*exp(-HistoricalFn[j]*vulnN[1,MaxAge])*exp(-NatMn[j])+ tempNn[j-1,MaxAge]*exp(-HistoricalFn[j]*vulnN[1,MaxAge])*exp(-NatMn[j])
      tempNs[j,MaxAge]	<-(tempNs[j-1,(MaxAge-1)])*exp(-HistoricalFs[j]*vulnS[1,MaxAge])*exp(-NatMs[j])+ tempNs[j-1,MaxAge]*exp(-HistoricalFs[j]*vulnS[1,MaxAge])*exp(-NatMs[j])

      EggsN			<-sum(tempNn[j-1,]*matureN[1,]*WeightAtAgeN[1,])
      EggsS			<-sum(tempNs[j-1,]*matureS[1,]*WeightAtAgeS[1,])

      tempNn[j,1]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1,j],recType="BH",NatMin=NatMn[j],
                                 vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,], MaxAge = MaxAge)
      tempNs[j,1]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1,j],recType="BH",NatMin=NatMs[j],
                                 vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,], MaxAge = MaxAge)
      tempRecN[j]		<-tempNn[j,1]
      tempRecS[j]		<-tempNs[j,1]

      tempCatchAtAgeN[j,]	<-((vulnN[1,]*HistoricalFn[j])/(vulnN[1,]*HistoricalFn[j]+NatMn[j])) * (1-exp(-(vulnN[1,]*HistoricalFn[j]+NatMn[j]))) * tempNn[j-1,]
      tempCatchN[j]		<-sum(tempCatchAtAgeN[j,]*WeightAtAgeN[1,])

      tempCatchAtAgeS[j,]	<-((vulnS[1,]*HistoricalFs[j])/(vulnS[1,]*HistoricalFs[j]+NatMs[j])) * (1-exp(-(vulnS[1,]*HistoricalFs[j]+NatMs[j]))) * tempNs[j-1,]
      tempCatchS[j]		<-sum(tempCatchAtAgeS[j,]*WeightAtAgeS[1,])

    }
  }

  #===============================================================
  # BEGIN SIMULATION OF ASSESSMENT AND HARVEST
  #===============================================================
  #==tempNn and tempNs from above are the starting points
  projNn	<-array(dim=c(SimYear,MaxAge,Nsim))
  projNs	<-array(dim=c(SimYear,MaxAge,Nsim))

  projCatchAtAgeN	<-array(dim=c(SimYear,MaxAge,Nsim))
  projCatchAtAgeS	<-array(dim=c(SimYear,MaxAge,Nsim))

  projCatchN	<-matrix(nrow=Nsim,ncol=SimYear)
  projCatchS	<-matrix(nrow=Nsim,ncol=SimYear)
  projBn	<-matrix(nrow=Nsim,ncol=SimYear)
  projBs	<-matrix(nrow=Nsim,ncol=SimYear)
  projSSBn	<-matrix(nrow=Nsim,ncol=SimYear)
  projSSBs	<-matrix(nrow=Nsim,ncol=SimYear)
  projExpBn	<-matrix(nrow=Nsim,ncol=SimYear)
  projExpBs	<-matrix(nrow=Nsim,ncol=SimYear)
  projSurvN	<-matrix(nrow=Nsim,ncol=SimYear)
  projSurvS	<-matrix(nrow=Nsim,ncol=SimYear)
  projRecN	<-matrix(nrow=Nsim,ncol=SimYear)
  projRecS	<-matrix(nrow=Nsim,ncol=SimYear)
  projFmortN	<-matrix(nrow=Nsim,ncol=SimYear)
  projFmortS	<-matrix(nrow=Nsim,ncol=SimYear)

  projCatLenFreqN	<-array(dim=c(SimYear,MaxAge,Nsim))
  projCatLenFreqS	<-array(dim=c(SimYear,MaxAge,Nsim))
  projSurvLenFreqN	<-array(dim=c(SimYear,MaxAge,Nsim))
  projSurvLenFreqS	<-array(dim=c(SimYear,MaxAge,Nsim))

  #==============================================================
  # calculate the catch proportion at (length) for assessment
  #==============================================================
  tempCatAtLenN<-matrix(nrow=MaxAge,ncol=MaxAge)
  tempCatAtLenS<-matrix(nrow=MaxAge,ncol=MaxAge)

  for(x in 1:Nsim)
    for(y in 2:InitYear)
    {
      #==make length frequencies for catch==
      for(w in 1:nrow(tempCatAtLenN))
      {
        probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
        ProbN<-probtemp/sum(probtemp)
        tempCatAtLenN[w,]<-tempCatchAtAgeN[y,w]*ProbN
        probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
        ProbS<-probtemp/sum(probtemp)
        tempCatAtLenS[w,]<-tempCatchAtAgeS[y,w]*ProbS
      }
      projCatLenFreqN[y,,x]<-apply(tempCatAtLenN,2,sum)
      projCatLenFreqS[y,,x]<-apply(tempCatAtLenS,2,sum)
    }


  #==============================================================
  # calculate the survey proportion at length for assessment
  #==============================================================
  tempSurvAtLenN<-matrix(nrow=MaxAge,ncol=MaxAge)
  tempSurvAtLenS<-matrix(nrow=MaxAge,ncol=MaxAge)
  for(x in 1:Nsim)
    for(y in 2:InitYear)
    {
      #==make length frequencies for s==
      for(w in 1:nrow(tempCatAtLenN))
      {
        probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn)
        ProbN<-probtemp/sum(probtemp)
        tempSurvAtLenN[w,] <- tempNn[y,w]*survSelN[y,w]*ProbN
        probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs)
        ProbS<-probtemp/sum(probtemp)
        tempSurvAtLenS[w,]<-tempNs[y,w]*survSelS[y,w]*ProbS
      }
      projSurvLenFreqN[y,,x]<-apply(tempSurvAtLenN,2,sum)
      projSurvLenFreqS[y,,x]<-apply(tempSurvAtLenS,2,sum)
    }

  #==transfer historical time series from above
  for(x in 1:Nsim)
  {
    projNn[1:InitYear,,x]			<-tempNn
    projNs[1:InitYear,,x]			<-tempNs
    projCatchN[x,1:InitYear]		<-tempCatchN
    projCatchS[x,1:InitYear]		<-tempCatchS
    projRecN[x,1:InitYear]			<-tempRecN
    projRecS[x,1:InitYear]			<-tempRecS
    projFmortN[x,1:InitYear]		<-HistoricalF
    projFmortS[x,1:InitYear]		<-HistoricalF
    projSurvN[x,1:InitYear]			<-apply(tempNn*survSelN[1:InitYear,]*WeightAtAgeN[1:InitYear,],1,sum)
    projSurvS[x,1:InitYear]			<-apply(tempNs*survSelS[1:InitYear,]*WeightAtAgeS[1:InitYear,],1,sum)
  }

  for(x in 1:InitYear)
  {
    projBn[,x] <- sum(projNn[x,,1]*WeightAtAgeN[1,])
    projBs[,x]	<-sum(projNs[x,,1]*WeightAtAgeS[1,])
    projSSBn[,x]	<- sum(projNn[x,,1]*matureN[1,]*WeightAtAgeN[1,])
    projSSBs[,x]	<-sum(projNs[x,,1]*matureS[1,]*WeightAtAgeS[1,])
    projExpBn[,x]	<-sum(projNn[x,,1]*vulnN[1,]*WeightAtAgeN[1,])
    projExpBs[,x]	<-sum(projNs[x,,1]*vulnS[1,]*WeightAtAgeS[1,])
  }

  CatchErrorN		<-matrix(rnorm(Nsim*SimYear,0,CatchCVn),nrow=Nsim,ncol=SimYear)
  CatchErrorS		<-matrix(rnorm(Nsim*SimYear,0,CatchCVs),nrow=Nsim,ncol=SimYear)
  CatchAssessN	<-projCatchN*exp(CatchErrorN)
  CatchAssessS	<-projCatchS*exp(CatchErrorS)

  CPUEErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
  CPUEErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
  CPUEAssessN		<-projExpBn*exp(CPUEErrorN)
  CPUEAssessS		<-projExpBs*exp(CPUEErrorS)

  SurvErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
  SurvErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
  SurvAssessN		<-projSurvN*exp(SurvErrorN)
  SurvAssessS		<-projSurvS*exp(SurvErrorS)

  #==storage for management and assessment quantities
  FMSY<-matrix(nrow=Nsim,ncol=SimYear)
  BMSY<-matrix(nrow=Nsim,ncol=SimYear)
  TAC<-matrix(nrow=Nsim,ncol=SimYear)
  CurBio<-matrix(nrow=Nsim,ncol=SimYear)


  # Tidy Up Data ------------------------------------------------------------

  # Deal with Numbers at age

  total_numbers_n <- adply(projNn, c(3)) %>% #convert array to dataframe with index for iteration
    dplyr::mutate(year = 1:SimYear) %>%
    gather('age','total_numbers',2:(MaxAge+1)) %>%
    dplyr::rename(iteration = X1) %>%
    dplyr::mutate(region = 'north')

  total_numbers_s <- adply(projNs, c(3)) %>% #convert array to dataframe with index for iteration
    dplyr::mutate(year = 1:SimYear) %>%
    gather('age','total_numbers',2:(MaxAge+1)) %>%
    dplyr::rename(iteration = X1) %>%
    dplyr::mutate(region = 'south')

  total_numbers <- rbind(total_numbers_n,total_numbers_s)

#   numbers <- left_join(total_numbers, survey_numbers, by = c('iteration','year','age','region')) %>%
    numbers <- total_numbers %>%
    gather('number_type','numbers_at_age', contains('_numbers')) #%>%

  length_at_age_n <- as.data.frame(LenAtAgeN) %>%
    dplyr::mutate(year = 1:SimYear, region = 'north') %>%
    gather('age','mean_length', contains('V')) %>%
    dplyr::mutate(age = gsub('V','', age))

  length_at_age_s <- as.data.frame(LenAtAgeS) %>%
    dplyr::mutate(year = 1:SimYear, region = 'south') %>%
    gather('age','mean_length', contains('V')) %>%
    dplyr::mutate(age = gsub('V','', age))

  length_at_age <- rbind(length_at_age_n, length_at_age_s)

  age_structure <- numbers %>%
    left_join(length_at_age, by = c('year','region','age')) %>%
    dplyr::mutate(age = as.numeric(age)) %>%
    arrange(age)

  # Deal with length frequencies


  length_mid_key <- data_frame(length_bin_index = 1:length(LengthBinsMid), LengthBinsMid = LengthBinsMid)

  survey_lenfreq_n <- adply(projSurvLenFreqN, c(3)) %>% #convert array to dataframe with index for iteration
    dplyr::mutate(year = 1:SimYear) %>%
    gather('length_bin_index','survey_numbers',2:(length(LengthBinsMid)+1)) %>%
    dplyr::rename(iteration = X1) %>%
    dplyr::mutate(region = 'north', length_bin_index = as.numeric(length_bin_index)) %>%
      left_join(length_mid_key, by = 'length_bin_index') %>%
    dplyr::select(-length_bin_index)

  survey_lenfreq_s <- adply(projSurvLenFreqS, c(3)) %>% #convert array to dataframe with index for iteration
    dplyr::mutate(year = 1:SimYear) %>%
    gather('length_bin_index','survey_numbers',2:(length(LengthBinsMid)+1)) %>%
    dplyr::rename(iteration = X1) %>%
    dplyr::mutate(region = 'south', length_bin_index = as.numeric(length_bin_index)) %>%
    left_join(length_mid_key, by = 'length_bin_index') %>%
    dplyr::select(-length_bin_index)


  survey_length_frequencies <- rbind(survey_lenfreq_n, survey_lenfreq_s)

  catch_lenfreq_n <- adply(projCatLenFreqN, c(3)) %>% #convert array to dataframe with index for iteration
    dplyr::mutate(year = 1:SimYear) %>%
    gather('length_bin_index','catch_numbers',2:(length(LengthBinsMid)+1)) %>%
    dplyr::rename(iteration = X1) %>%
    dplyr::mutate(region = 'north', length_bin_index = as.numeric(length_bin_index)) %>%
    left_join(length_mid_key, by = 'length_bin_index') %>%
    dplyr::select(-length_bin_index)


  catch_lenfreq_s <- adply(projCatLenFreqS, c(3)) %>% #convert array to dataframe with index for iteration
    dplyr::mutate(year = 1:SimYear) %>%
    gather('length_bin_index','catch_numbers',2:(length(LengthBinsMid)+1)) %>%
    dplyr::rename(iteration = X1) %>%
    dplyr::mutate(region = 'south', length_bin_index = as.numeric(length_bin_index)) %>%
    left_join(length_mid_key, by = 'length_bin_index') %>%
    dplyr::select(-length_bin_index)


  catch_length_frequencies <- rbind(catch_lenfreq_n, catch_lenfreq_s)

  length_frequencies <- left_join(catch_length_frequencies, survey_length_frequencies,
                                  by = c('iteration', 'year','region','LengthBinsMid')) %>%
    gather('lenfreq_type','numbers', contains('_numbers')) %>%
    dplyr::mutate(lenfreq_type = gsub("\\_.*","",lenfreq_type))

  # Deal with Biomass

  biomass_n <- data_frame(year = 1:SimYear, region = 'north', total_biomass = as.vector(projBn))

  biomass_s <- data_frame(year = 1:SimYear, region = 'south', total_biomass = as.vector(projBs))

  total_biomass <- rbind(biomass_n, biomass_s)

  ss_biomass_n <- data_frame(year = 1:SimYear, region = 'north', ss_biomass = as.vector(projSSBn))

  ss_biomass_s <- data_frame(year = 1:SimYear, region = 'south', ss_biomass = as.vector(projSSBs))

  ss_biomass <- rbind(ss_biomass_n, ss_biomass_s)

  exp_biomass_n <- data_frame(year = 1:SimYear, region = 'north', exp_biomass = as.vector(projExpBn))

  exp_biomass_s <- data_frame(year = 1:SimYear, region = 'south', exp_biomass = as.vector(projExpBs))

  exp_biomass <- rbind(exp_biomass_n, exp_biomass_s)

  survey_biomass_n <- data_frame(year = 1:SimYear, region = 'north', survey_biomass = as.vector(SurvAssessN))

  survey_biomass_s <- data_frame(year = 1:SimYear, region = 'south', survey_biomass = as.vector(SurvAssessS))

  survey_biomass <- rbind(survey_biomass_n, survey_biomass_s)

  biomass <- left_join(total_biomass, ss_biomass, by = c('year', 'region')) %>%
    left_join(exp_biomass, by = c('year','region')) %>%
    left_join(survey_biomass, by = c('year', 'region')) %>%
    gather('biomass_type','biomass', contains('_biomass')) %>%
    dplyr::mutate(biomass_type = gsub("\\_.*","",biomass_type))

  # Deal with catch

  catch_n <- data_frame(year = 1:SimYear, region = 'north', total_catch = as.vector(projCatchN))

  catch_s <- data_frame(year = 1:SimYear, region = 'south', total_catch = as.vector(projCatchS))

  total_catch <- rbind(catch_n, catch_s)

  surv_catch_n <- data_frame(year = 1:SimYear, region = 'north', survey_catch = as.vector(CatchAssessN))

  surv_catch_s <- data_frame(year = 1:SimYear, region = 'south', survey_catch = as.vector(CatchAssessS))

  surv_catch <- rbind(surv_catch_n, surv_catch_s)

  catch <- left_join(total_catch, surv_catch, by = c('year', 'region')) %>%
    gather('catch_type','catch', contains('_catch')) %>%
    dplyr::mutate(catch_type = gsub("\\_.*","",catch_type))

  # Deal with CPUE

  cpue_n <- data_frame(year = 1:SimYear, region = 'north', cpue = as.vector(CPUEAssessN))

  cpue_s <- data_frame(year = 1:SimYear, region = 'south', cpue = as.vector(CPUEAssessS))

  cpue <- rbind(cpue_n, cpue_s)


# Make Froese Indicators --------------------------------------------------

# Calculate Pmat: proportion of catch that is mature

  length_at_age_key<- data_frame(age = 1:MaxAge, mean_length = LenAtAgeS[1,])

  calc_lopt <- data_frame(age = 1:MaxAge, virgin_bio = VirInitN*WeightAtAgeN[1,] ) %>%
    left_join(length_at_age_key, by = 'age') %>%
    mutate(max_bio= virgin_bio == max(virgin_bio, na.rm = T)) %>%
    filter(max_bio == T)

  lopt <- calc_lopt$mean_length

  length_50_mat <- LinfN[1]*(1-exp(-VonKn[1]*(mean(c(mat50n,mat50s))-t0n[1])))


froese_indicators <- catch_length_frequencies %>%
    group_by(iteration,year,LengthBinsMid) %>%
    summarise(catch_at_length = sum(catch_numbers, na.rm = T)) %>%
    ungroup() %>%
    group_by(iteration,year) %>%
      mutate(total_catch = sum(catch_at_length, na.rm = T) ,
             is_mature = LengthBinsMid >= length_50_mat,
             is_opt =  LengthBinsMid >= (0.9 * lopt) & LengthBinsMid <= (1.1 * lopt),
             is_mega =  LengthBinsMid >= (1.1 * lopt) & LengthBinsMid <= mean(c(LinfN, LinfS))) %>%
  ungroup() %>%
  mutate( pl_mat = (catch_at_length / total_catch) * is_mature,
          pl_opt = (catch_at_length / total_catch) * is_opt,
          pl_mega = (catch_at_length / total_catch) * is_mega) %>%
#   filter(LengthBinsMid > length_50_mat & total_catch >0)
    group_by(year) %>%
  summarise(p_mat = sum(pl_mat, na.rm = T), p_opt = sum(pl_opt, na.rm = T), p_mega = sum(pl_mega, na.rm = T)) %>%
  ungroup() %>%
  mutate(p_obj = p_mat + p_opt + p_mega )

cope_punt_plot <- froese_indicators %>%
  gather('Indicator','Value',contains('p_')) %>%
  ggplot(aes(year,Value, fill = Indicator)) +
  geom_line(aes(color = Indicator), size = 1, alpha = 0.75) +
  geom_point(shape = 21, size = 2) +
  plot_theme #+
#   scale_fill_economist() +
#   scale_color_economist()


# Make Plots --------------------------------------------------------------


  length_ogive <- length_at_age %>%
    mutate(age = as.numeric(age)) %>%
    group_by(age) %>%
    summarise(mean_length = mean(mean_length, na.rm = T))

  maturity_ogive <- as.data.frame(matureS) %>%
    mutate(year = 1:SimYear) %>%
    gather('age','maturity',1:MaxAge) %>%
    mutate(age = gsub('V','',age)) %>%
    mutate(age = as.numeric(age)) %>%
    group_by(age) %>%
    summarise(mean_maturity = mean(maturity))

  weight_ogive <- as.data.frame(WeightAtAgeN) %>%
    mutate(year = 1:SimYear) %>%
    gather('age','weight',1:MaxAge) %>%
    mutate(age = gsub('V','',age)) %>%
    mutate(age = as.numeric(age)) %>%
    group_by(age) %>%
    summarise(mean_weight = mean(weight))

  movement_ogive <- as.data.frame(MovementS) %>%
    mutate(year = 1:SimYear) %>%
    gather('age','movement',1:MaxAge) %>%
    mutate(age = gsub('V','',age)) %>%
    mutate(age = as.numeric(age)) %>%
    group_by(age) %>%
    summarise(mean_movement = mean(movement))

  selectivity_ogive <- as.data.frame(vulnN) %>%
    mutate(year = 1:SimYear) %>%
    gather('age','selectivity',1:MaxAge) %>%
    mutate(age = gsub('V','',age)) %>%
    mutate(age = as.numeric(age)) %>%
    group_by(age) %>%
    summarise(mean_selectivity= mean(selectivity))

  life_plot <- length_ogive %>%
    left_join(maturity_ogive, by = 'age') %>%
    left_join(weight_ogive, by = 'age') %>%
    left_join(movement_ogive, by = 'age') %>%
    left_join(selectivity_ogive, by = 'age') %>%
    gather('metric','value', contains('mean_')) %>%
    ggplot(aes(age,value)) +
    geom_line(size = 2) +
    facet_wrap(~metric, scales = 'free_y') +
    plot_theme

  ssb<-seq(1,VirBioN,VirBioN/100)
  record<-ssb
  for(x in 1:length(record))
  {
    record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1],recType="BH",NatMin=NatMn[1],
                           vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,], MaxAge = MaxAge)
  }

  recruitment_plot <- data_frame(ssb = ssb, recruits = record) %>%
    ggplot(aes(ssb, recruits)) +
    geom_line(size = 2) +
    plot_theme


#   plot(record~ssb)

catch_plot <- catch %>%
    group_by(year,catch_type) %>%
    summarise(total_catch = sum(catch)) %>%
    ggplot(aes(year,total_catch, fill = catch_type)) +
    geom_point(shape = 21, size = 2) +
    xlab('Year') +
    ylab('Catch') +
  plot_theme

  cpue_plot <- cpue %>%
    group_by(year) %>%
    summarise(mean_cpue = mean(cpue)) %>%
    ggplot(aes(year,mean_cpue)) +
    geom_point(shape = 21, size = 2) +
    xlab('Year') +
    ylab('CPUE') +
    plot_theme

  biomass_plot <- biomass %>%
    group_by(year,biomass_type) %>%
    summarise(total_biomass = sum(biomass))  %>%
    ggplot(aes(year,total_biomass, fill = biomass_type)) +
    geom_point(shape = 21, alpha = 0.75, size = 2) +
    xlab('Year') +
    ylab('Biomass') +
    plot_theme
#   arg <- NA
#
#   for (i in 1:30){
#
#     arg[i] <- sum(age_structure$numbers_at_age[age_structure$age == i], na.rm = T)
#   }

  length_plot <- length_frequencies %>%
    ungroup() %>%
    group_by(iteration,year,LengthBinsMid,lenfreq_type) %>%
    summarize(total_numbers = sum(numbers, na.rm = T)) %>%
    ggplot(aes(LengthBinsMid,total_numbers, fill = lenfreq_type)) +
    geom_density(stat = 'identity', alpha = 0.6, color = 'black') +
    facet_wrap(~year) +
    xlab('Length Bin (cm)') +
    ylab('Numbers') +
    plot_theme


  age_structure_plot <- age_structure %>%
    ungroup() %>%
    group_by(iteration, year, age, number_type) %>%
    summarise(total_numbers = sum(numbers_at_age, na.rm = T)) %>%
    ggplot(aes(age,total_numbers, fill = number_type)) +
    geom_density(stat = 'identity',alpha = 0.6, color = 'black') +
    facet_wrap(~year) +
    xlab('Age (years)') +
    ylab('Numbers') +
    plot_theme

  local_files <- ls()

  plot_files <- local_files[grep('_plot', local_files, fixed = T)]

  plot_list <- list()

  for (i in 1:length(plot_files))
  {
    eval(parse( text = paste('plot_list$',plot_files[i],' <- ',plot_files[i], sep = '')))
  }


return(list(catch_data = catch, cpue_data = cpue, age_data = age_structure, length_data = length_frequencies, biomass_data = biomass,
            plots = plot_list))
}