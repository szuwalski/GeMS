#' Read in CTL csv file
#'
#' @param input Name of csv file
#'
#' @return List of inputs from CTL file
#'
#' @export
#' @examples
#' \dontrun{
#' ReadCTLfile("Cod_Base_CTL")
#' 
#' ReadCTLfile("~/GeMS/Examples/Cod_1_Production/Cod_Base_CTL")
#' }
#' 

ReadCTLfile<-function(input)
{
 OM						<-NULL
 CTLfilename 			<- paste0(input,".csv")
 SearchPool				<-as.matrix(read.csv(CTLfilename))
 OM$Nsim				<-TakeOut("Nsim",SearchPool)
 OM$SimYear				<-TakeOut("SimYear",SearchPool)
 OM$InitYear			<-TakeOut("InitYear",SearchPool)
 OM$AssessmentType		<-TakeOut("AssessmentType",SearchPool)
 OM$FisheryIndepenDat	<-TakeOut("FisheryIndepenDat",SearchPool)
 OM$LifeHistoryPlots	<-TakeOut("LifeHistoryPlots",SearchPool)
 OM$EstimationPlots		<-TakeOut("EstimationPlots",SearchPool)
 OM$AssessmentData		<-TakeOut("AssessmentData",SearchPool)
 OM$PlotYieldCurve		<-TakeOut("PlotYieldCurve",SearchPool)

 OM$CatchCVn			<-TakeOut("CatchCVn",SearchPool)
 OM$CatchCVs			<-TakeOut("CatchCVs",SearchPool)
 OM$CPUEcvN				<-TakeOut("CPUEcvN",SearchPool)
 OM$CPUEcvS				<-TakeOut("CPUEcvS",SearchPool)
 OM$IndexCVn			<-TakeOut("IndexCVn",SearchPool)
 OM$IndexCVs			<-TakeOut("IndexCVs",SearchPool)
 OM$LenSampleN			<-TakeOut("LenSampleN",SearchPool)
 OM$LenSampleS			<-TakeOut("LenSampleS",SearchPool)
 OM$GrowthSDn			<-TakeOut("GrowthSDn",SearchPool)
 OM$GrowthSDs			<-TakeOut("GrowthSDs",SearchPool)
 OM$LengthBinN			<-TakeOut("LengthBinN",SearchPool)

 OM$MaxAge				<-TakeOut("MaxAge",SearchPool)
 OM$NatMn				<-TakeOut("NatMn",SearchPool)
 OM$VonKn				<-TakeOut("VonKn",SearchPool)
 OM$LinfN				<-TakeOut("LinfN",SearchPool)
 OM$t0n					<-TakeOut("t0n",SearchPool)
 OM$mat50n				<-TakeOut("mat50n",SearchPool)
 OM$mat95n				<-TakeOut("mat95n",SearchPool)
 OM$alphaN				<-TakeOut("alphaN",SearchPool)
 OM$betaN				<-TakeOut("betaN",SearchPool)
 OM$MaxMovingN			<-TakeOut("MaxMovingN",SearchPool)
 OM$Move50n				<-TakeOut("Move50n",SearchPool)
 OM$Move95n				<-TakeOut("Move95n",SearchPool)
 OM$steepnessN			<-TakeOut("steepnessN",SearchPool)
 OM$RzeroN				<-TakeOut("RzeroN",SearchPool)
 OM$sigmaRn				<-TakeOut("sigmaRn",SearchPool)

 OM$sel50n				<-TakeOut("sel50n",SearchPool)
 OM$sel95n				<-TakeOut("sel95n",SearchPool)
 OM$surv50n				<-TakeOut("surv50n",SearchPool)
 OM$surv95n				<-TakeOut("surv95n",SearchPool)

 OM$NatMs				<-TakeOut("NatMs",SearchPool)
 OM$VonKs				<-TakeOut("VonKs",SearchPool)
 OM$LinfS				<-TakeOut("LinfS",SearchPool)
 OM$t0s					<-TakeOut("t0s",SearchPool)
 OM$mat50s				<-TakeOut("mat50s",SearchPool)
 OM$mat95s				<-TakeOut("mat95s",SearchPool)
 OM$alphaS				<-TakeOut("alphaS",SearchPool)
 OM$betaS				<-TakeOut("betaS",SearchPool)
 OM$MaxMovingS			<-TakeOut("MaxMovingS",SearchPool)
 OM$Move50s				<-TakeOut("Move50s",SearchPool)
 OM$Move95s				<-TakeOut("Move95s",SearchPool)
 OM$steepnessS			<-TakeOut("steepnessS",SearchPool)
 OM$RzeroS				<-TakeOut("RzeroS",SearchPool)
 OM$sigmaRs				<-TakeOut("sigmaRs",SearchPool)

 OM$sel50s				<-TakeOut("sel50s",SearchPool)
 OM$sel95s				<-TakeOut("sel95s",SearchPool)
 OM$surv50s				<-TakeOut("surv50s",SearchPool)
 OM$surv95s				<-TakeOut("surv95s",SearchPool)

 OM$HistoricalFn		<-TakeOut(SearchTerm="HistoricalFn",SearchPool=SearchPool)
 OM$HistoricalFs		<-TakeOut(SearchTerm="HistoricalFs",SearchPool=SearchPool)
 OM$PastFsdN			<-TakeOut("PastFsdN",SearchPool)
 OM$PastFsdS			<-TakeOut("PastFsdS",SearchPool)

 OM$HarvestControlN		<-TakeOut("HarvestControlN",SearchPool)
 OM$ConstantCatchN		<-TakeOut("ConstantCatchN",SearchPool)
 OM$ConstantFn			<-TakeOut("ConstantFn",SearchPool)
 OM$HCalphaN			<-TakeOut("HCalphaN",SearchPool)
 OM$HCbetaN				<-TakeOut("HCbetaN",SearchPool)

 OM$HarvestControlS		<-TakeOut("HarvestControlS",SearchPool)
 OM$ConstantCatchS		<-TakeOut("ConstantCatchS",SearchPool)
 OM$ConstantFs			<-TakeOut("ConstantFs",SearchPool)
 OM$HCalphaS			<-TakeOut("HCalphaS",SearchPool)
 OM$HCbetaS				<-TakeOut("HCbetaS",SearchPool)

 OM$start_assessment	<- TakeOut("start_assessment",SearchPool)
 OM$SmalNum				<-TakeOut("SmallNum",SearchPool)
 OM$InitSmooth			<-TakeOut("InitSmooth",SearchPool)
 OM$FmortPen			<-TakeOut("FmortPen",SearchPool)
 OM$RecruitPen			<-TakeOut("RecruitPen",SearchPool)
 OM$Mpenalty			<-TakeOut("Mpenalty",SearchPool)
 OM$EstM				<-TakeOut("EstM",SearchPool)
 OM$TimeVaryM			<-TakeOut("TimeVaryM",SearchPool)

 OM$EstGrowthK			<-TakeOut("EstGrowthK",SearchPool)
 OM$TimeVaryGrowthK		<-TakeOut("TimeVaryGrowthK",SearchPool)
 OM$EstLinf				<-TakeOut("EstLinf",SearchPool)
 OM$TimeVaryLinf		<-TakeOut("TimeVaryLinf",SearchPool)
 OM$Growthpenalty		<-TakeOut("Growthpenalty",SearchPool)
 OM$TimeVarySel50		<-TakeOut("TimeVarySel50",SearchPool)
 OM$TimeVarySel95		<-TakeOut("TimeVarySel95",SearchPool)
 OM$SelPenalty			<-TakeOut("SelPenalty",SearchPool)
 OM$ProjectTimeVary		<-TakeOut("ProjectTimeVary",SearchPool)
 OM$InitValSig			<-TakeOut("InitValSig",SearchPool)

 OM$InitBzeroMod		<-TakeOut("InitBzeroMod",SearchPool)
 OM$InitGrowthRate		<-TakeOut("InitGrowthRate",SearchPool)
 OM$InitShape			<-TakeOut("InitShape",SearchPool)
 OM$estInit         	<-TakeOut("estInit",SearchPool)
 OM$InitBioProd     	<-TakeOut("InitBioProd",SearchPool)
 
 OM$TwoPop				<-TakeOut("TwoPop",SearchPool)

 list(OM=OM)
}