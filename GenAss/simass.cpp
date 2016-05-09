 #include <math.h>
 #include <admodel.h>
  #include <time.h>
 
 ofstream CheckFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
 dvariable gammln(dvariable& xx)
 {
   RETURN_ARRAYS_INCREMENT();	
        dvariable x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
               0.1208650973866179e-2,-0.5395239384953e-5}; 
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) 
          {
            ser += cof[j]/(y+1.);
          }
         dvariable value_=-tmp+log(2.5066282746310005*ser/x);
    RETURN_ARRAYS_DECREMENT();         
         return(value_); 
 }
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <simass.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_post = new ofstream("eval.csv");;
 call_no = 0;
  styr.allocate("styr");
  endyr.allocate("endyr");
  maxAge.allocate("maxAge");
  Ages.allocate(1,maxAge,"Ages");
  cpueYr.allocate("cpueYr");
  cpueYears.allocate(1,cpueYr,"cpueYears");
  survYr.allocate("survYr");
  survYears.allocate(1,survYr,"survYears");
  catchLenYr.allocate("catchLenYr");
  catchLenYears.allocate(1,survYr,"catchLenYears");
  catchSampN.allocate(styr,endyr,"catchSampN");
  survLenYr.allocate("survLenYr");
  survLenYears.allocate(1,survYr,"survLenYears");
  survSampN.allocate(styr,endyr,"survSampN");
cout<<"SurvLenYear"<<survLenYears<<endl;
  survBiomass.allocate(1,survYr,"survBiomass");
  survCV.allocate(styr,endyr,"survCV");
  cpueIndex.allocate(1,cpueYr,"cpueIndex");
  cpueCV.allocate(styr,endyr,"cpueCV");
  catchBiomass.allocate(1,survYr,"catchBiomass");
  catchCV.allocate(styr,endyr,"catchCV");
  LengthBinsMid.allocate(1,maxAge,"LengthBinsMid");
cout<<"catchBiomass"<<catchBiomass<<endl; 
  catchLenNum.allocate(styr,endyr,1,maxAge,"catchLenNum");
  survLenNum.allocate(styr,endyr,1,maxAge,"survLenNum");
 ad_comm::change_datafile_name("SimAss.ctl");
  weightPars.allocate(1,2,"weightPars");
  projectTimeVary.allocate("projectTimeVary");
  EstGrowthK.allocate("EstGrowthK");
  TimeVaryGrowthK.allocate("TimeVaryGrowthK");
  EstLinf.allocate("EstLinf");
  TimeVaryLinf.allocate("TimeVaryLinf");
  lengthParsIn.allocate(1,3,"lengthParsIn");
  EstM.allocate("EstM");
  NatMin.allocate("NatMin");
  TimeVaryM.allocate("TimeVaryM");
  TimeVarySel50.allocate("TimeVarySel50");
  TimeVarySel95.allocate("TimeVarySel95");
  maturity.allocate(1,2,"maturity");
  steepness.allocate("steepness");
  smallNum.allocate("smallNum");
  InitSmoothWeight.allocate("InitSmoothWeight");
  FisheryIndependentData.allocate("FisheryIndependentData");
  GrowthSD.allocate("GrowthSD");
  FmortPen.allocate("FmortPen");
  RecruitPen.allocate("RecruitPen");
  Mpenalty.allocate("Mpenalty");
  Growthpenalty.allocate("Growthpenalty");
  SelPenalty.allocate("SelPenalty");
  HarvestControl.allocate("HarvestControl");
  ConstantCatch.allocate("ConstantCatch");
  ConstantF.allocate("ConstantF");
  HCalpha.allocate("HCalpha");
  HCbeta.allocate("HCbeta");
  LengthBinN.allocate("LengthBinN");
cout<<"steep"<<steepness<<endl;
 Nproj = 125;
 logMaxAge = log(maxAge);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  srv_sel50.allocate(0.01,maxAge,2,"srv_sel50");
  srv_sel95.allocate(1,maxAge,2,"srv_sel95");
  stNatLen.allocate(1,maxAge,-10,20,1,"stNatLen");
  log_avg_fmort_dir.allocate(-5,5,1,"log_avg_fmort_dir");
  fmort_dir_dev.allocate(styr,endyr,-20,20,1,"fmort_dir_dev");
  mean_log_rec.allocate(-5,25,1,"mean_log_rec");
  rec_dev.allocate(styr,endyr,-20,20,1,"rec_dev");
  log_avg_NatM.allocate(-5,3,EstM,"log_avg_NatM");
  m_dev.allocate(styr,endyr,-20,20,TimeVaryM,"m_dev");
  log_avg_GrowthK.allocate(-5,2,EstGrowthK,"log_avg_GrowthK");
  log_avg_Linf.allocate(0,8,EstLinf,"log_avg_Linf");
  GrowthK_dev.allocate(styr,endyr,-20,20,TimeVaryGrowthK,"GrowthK_dev");
  Linf_dev.allocate(styr,endyr,-20,20,TimeVaryLinf,"Linf_dev");
  SelPars50.allocate(0.01,maxAge,2,"SelPars50");
  SelPars95.allocate(1,maxAge,2,"SelPars95");
  SelPars_dev50.allocate(styr,endyr,-20,20,TimeVarySel50,"SelPars_dev50");
  SelPars_dev95.allocate(styr,endyr,-20,20,TimeVarySel95,"SelPars_dev95");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  catchLenFreq.allocate(styr,endyr,1,maxAge,"catchLenFreq");
  #ifndef NO_AD_INITIALIZE
    catchLenFreq.initialize();
  #endif
  survLenFreq.allocate(styr,endyr,1,maxAge,"survLenFreq");
  #ifndef NO_AD_INITIALIZE
    survLenFreq.initialize();
  #endif
  sel_fish.allocate(styr,endyr+Nproj,1,maxAge,"sel_fish");
  #ifndef NO_AD_INITIALIZE
    sel_fish.initialize();
  #endif
  sel_srv.allocate(1,maxAge,"sel_srv");
  #ifndef NO_AD_INITIALIZE
    sel_srv.initialize();
  #endif
  MatAtAge.allocate(1,maxAge,"MatAtAge");
  #ifndef NO_AD_INITIALIZE
    MatAtAge.initialize();
  #endif
  LengthAtAge.allocate(styr,endyr+Nproj,1,maxAge,"LengthAtAge");
  #ifndef NO_AD_INITIALIZE
    LengthAtAge.initialize();
  #endif
  WeightAtAge.allocate(styr,endyr+Nproj,1,maxAge,"WeightAtAge");
  #ifndef NO_AD_INITIALIZE
    WeightAtAge.initialize();
  #endif
  NatM.allocate(styr,endyr+Nproj,"NatM");
  #ifndef NO_AD_INITIALIZE
    NatM.initialize();
  #endif
  fmort_dir.allocate(styr,endyr+Nproj,"fmort_dir");
  #ifndef NO_AD_INITIALIZE
    fmort_dir.initialize();
  #endif
  F_dir.allocate(styr,endyr+Nproj,1,maxAge,"F_dir");
  #ifndef NO_AD_INITIALIZE
    F_dir.initialize();
  #endif
  F_trawl.allocate(styr,endyr+Nproj,1,maxAge,"F_trawl");
  #ifndef NO_AD_INITIALIZE
    F_trawl.initialize();
  #endif
  Baranov.allocate(1,maxAge,"Baranov");
  #ifndef NO_AD_INITIALIZE
    Baranov.initialize();
  #endif
  NatAge.allocate(styr,endyr+Nproj,1,maxAge,"NatAge");
  #ifndef NO_AD_INITIALIZE
    NatAge.initialize();
  #endif
  Spbio.allocate(styr,endyr+Nproj,"Spbio");
  predSurvBio.allocate(styr,endyr+Nproj,"predSurvBio");
  #ifndef NO_AD_INITIALIZE
    predSurvBio.initialize();
  #endif
  predCpueBio.allocate(styr,endyr+Nproj,"predCpueBio");
  #ifndef NO_AD_INITIALIZE
    predCpueBio.initialize();
  #endif
  predCpueBioLen.allocate(styr,endyr+Nproj,1,maxAge,"predCpueBioLen");
  #ifndef NO_AD_INITIALIZE
    predCpueBioLen.initialize();
  #endif
  predSurvNatAge.allocate(styr,endyr+Nproj,1,maxAge,"predSurvNatAge");
  #ifndef NO_AD_INITIALIZE
    predSurvNatAge.initialize();
  #endif
  predSurvAgeFreq.allocate(styr,endyr+Nproj,1,maxAge,"predSurvAgeFreq");
  #ifndef NO_AD_INITIALIZE
    predSurvAgeFreq.initialize();
  #endif
  predSurvLenFreq.allocate(styr,endyr+Nproj,1,maxAge,"predSurvLenFreq");
  #ifndef NO_AD_INITIALIZE
    predSurvLenFreq.initialize();
  #endif
  predCatchAtAge.allocate(styr,endyr+Nproj,1,maxAge,"predCatchAtAge");
  #ifndef NO_AD_INITIALIZE
    predCatchAtAge.initialize();
  #endif
  predCatchBio.allocate(styr,endyr+Nproj,"predCatchBio");
  predCatchAgeFreq.allocate(styr,endyr+Nproj,1,maxAge,"predCatchAgeFreq");
  #ifndef NO_AD_INITIALIZE
    predCatchAgeFreq.initialize();
  #endif
  predCatchLenFreq.allocate(styr,endyr+Nproj,1,maxAge,"predCatchLenFreq");
  #ifndef NO_AD_INITIALIZE
    predCatchLenFreq.initialize();
  #endif
  CpueBio_like.allocate("CpueBio_like");
  #ifndef NO_AD_INITIALIZE
  CpueBio_like.initialize();
  #endif
  SurvBio_like.allocate("SurvBio_like");
  #ifndef NO_AD_INITIALIZE
  SurvBio_like.initialize();
  #endif
  SurvLen_like.allocate("SurvLen_like");
  #ifndef NO_AD_INITIALIZE
  SurvLen_like.initialize();
  #endif
  Catch_like.allocate("Catch_like");
  #ifndef NO_AD_INITIALIZE
  Catch_like.initialize();
  #endif
  CatchLen_like.allocate("CatchLen_like");
  #ifndef NO_AD_INITIALIZE
  CatchLen_like.initialize();
  #endif
  initsmo_penal.allocate("initsmo_penal");
  #ifndef NO_AD_INITIALIZE
  initsmo_penal.initialize();
  #endif
  Rec_pen.allocate("Rec_pen");
  #ifndef NO_AD_INITIALIZE
  Rec_pen.initialize();
  #endif
  F_pen.allocate("F_pen");
  #ifndef NO_AD_INITIALIZE
  F_pen.initialize();
  #endif
  M_pen.allocate("M_pen");
  #ifndef NO_AD_INITIALIZE
  M_pen.initialize();
  #endif
  Growth_pen1.allocate("Growth_pen1");
  #ifndef NO_AD_INITIALIZE
  Growth_pen1.initialize();
  #endif
  Growth_pen2.allocate("Growth_pen2");
  #ifndef NO_AD_INITIALIZE
  Growth_pen2.initialize();
  #endif
  Sel50_pen.allocate("Sel50_pen");
  #ifndef NO_AD_INITIALIZE
  Sel50_pen.initialize();
  #endif
  Sel95_pen.allocate("Sel95_pen");
  #ifndef NO_AD_INITIALIZE
  Sel95_pen.initialize();
  #endif
  F35.allocate("F35");
  #ifndef NO_AD_INITIALIZE
  F35.initialize();
  #endif
  SBPRF35.allocate("SBPRF35");
  #ifndef NO_AD_INITIALIZE
  SBPRF35.initialize();
  #endif
  Bmsy.allocate("Bmsy");
  #ifndef NO_AD_INITIALIZE
  Bmsy.initialize();
  #endif
  FOFL.allocate("FOFL");
  #ifndef NO_AD_INITIALIZE
  FOFL.initialize();
  #endif
  OFL.allocate("OFL");
  #ifndef NO_AD_INITIALIZE
  OFL.initialize();
  #endif
  FutRec.allocate("FutRec");
  #ifndef NO_AD_INITIALIZE
  FutRec.initialize();
  #endif
  FutMort.allocate("FutMort");
  #ifndef NO_AD_INITIALIZE
  FutMort.initialize();
  #endif
  Bzero.allocate("Bzero");
  #ifndef NO_AD_INITIALIZE
  Bzero.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 int i,j;
 dvariable nn;
  // for(i=1;i<=maxAge;i++)
   // WeightAtAge(i) = weightPars(1) * pow(LengthBinsMid(i),weightPars(2));
   
  
    for(i=styr;i<=endyr;i++)
    {
	   nn =0;
      for(j=1;j<=LengthBinN;j++)
       nn +=  catchLenNum(i,j);
      for(j=1;j<=LengthBinN;j++)
	  catchLenFreq(i,j) = catchLenNum(i,j)/nn;
    }
  for(i=styr;i<=endyr;i++)
    {
	   nn =0;
      for(j=1;j<=LengthBinN;j++)
       nn +=  survLenNum(i,j);
      for(j=1;j<=LengthBinN;j++)
	  survLenFreq(i,j) = survLenNum(i,j)/nn;
    }
}

void model_parameters::initializationfunction(void)
{
}

void model_parameters::userfunction(void)
{
  f =0.0;
  ofstream& post= *pad_post;
  getselectivity();
   // cout<<"Select"<<endl;
  getmortality();
  // cout<<"Mort"<<endl;
  getmaturity();
  // cout<<"Mat"<<endl;
  get_num_at_len();
  // cout<<"Numbers"<<endl;
  evaluate_the_objective_function();
    if (mceval_phase())
   {
    // Find_F35();
    // Find_OFL();
    // post << Bmsy << " " << F35 << " " << FOFL << " " << OFL << " ";   
    // post << srv_sel50 << " " << srv_sel95 << " " << log_avg_fmort_dir << " " << fmort_dir_dev << " " << log_avg_fmort_trawl << " ";
    // post << mean_log_rec << " " << rec_dev << " " <<endl;
   }
}

void model_parameters::getselectivity(void)
{
  ofstream& post= *pad_post;
 int  i,j,k;
 dvar_vector tempSelProj(1,maxAge);
    for (j=1;j<=maxAge;j++)
      sel_srv(j)		   =  1./(1.+mfexp(-1.*log(19.)*(Ages(j)-srv_sel50)/(srv_sel95-srv_sel50)));
 //calculate fishery selectivity by year
  for(i=styr;i<=endyr;i++)
    for (j=1;j<=maxAge;j++)
	 sel_fish(i,j)       =  1./(1.+mfexp(-1.*log(19.)*(Ages(j)-(SelPars50+SelPars_dev50(i)))/((SelPars95+SelPars_dev95(i))-(SelPars50+SelPars_dev50(i)))));
	 // sel_fish(i,j)       =  1./(1.+mfexp(-1.*log(19.)*(Ages(j)-mfexp(log_avg_SelPars50+SelPars_dev50(i)))/(mfexp(log_avg_SelPars95+SelPars_dev95(i))-mfexp(log_avg_SelPars50+SelPars_dev50(i)))));
 if(TimeVarySel50>0  | TimeVarySel95>0)
	{
       for(k=1;k<=maxAge;k++)	
	    for(j=endyr-projectTimeVary+1;j<=endyr;j++)	
			tempSelProj(k)+=sel_fish(j,k);
        tempSelProj = tempSelProj/(projectTimeVary);
	}
 if(TimeVarySel50<0  & TimeVarySel50<0)
   tempSelProj = sel_fish(endyr);
  for(i=endyr+1;i<=endyr+Nproj;i++)
   sel_fish(i) = tempSelProj;
}

void model_parameters::getmaturity(void)
{
  ofstream& post= *pad_post;
 int j,i,k;
  dvar_vector tempLenAtAge(1,maxAge);
  dvar_vector tempWgtAtAge(1,maxAge);
  //calculate length at age
  for(i=styr;i<=endyr;i++)
    for (j=1;j<=maxAge;j++)
   {
	  if(EstGrowthK>0)
          LengthAtAge(i,j) = lengthParsIn(2)*(1-mfexp(-1*mfexp(log_avg_GrowthK+GrowthK_dev(i))*(Ages(j)-lengthParsIn(3))));		
	  if(EstLinf>0)
        LengthAtAge(i,j) = mfexp(log_avg_Linf+Linf_dev(i))*(1-mfexp(-1*lengthParsIn(1)*(Ages(j)-lengthParsIn(3))));		
	  if(EstGrowthK<=0 & EstLinf<=0)
	    LengthAtAge(i,j) = lengthParsIn(2)*(1-mfexp(-1*lengthParsIn(1)*(Ages(j)-lengthParsIn(3))));
	//  LengthAtAge(i,j) = Len1 + (Len2 - Len1)* ( (1-mfexp(-1*mfexp(log_avg_GrowthK+GrowthK_dev(i))*(Ages(j)-lengthParsIn(3))))/(1-mfexp(-1*mfexp(log_avg_GrowthK+GrowthK_dev(i))*(MaxAge-lengthParsIn(3))))
	 WeightAtAge(i,j) = weightPars(1) * pow(LengthAtAge(i,j),weightPars(2));
   }
  //==taking the average of the calculated values rather than the average of the parameters may introduce wonkiness
 tempLenAtAge.initialize();
 tempWgtAtAge.initialize(); 
 if(TimeVaryGrowthK>0  | TimeVaryLinf>0)
	{
       for(k=1;k<=maxAge;k++)	
	    for(j=endyr-projectTimeVary+1;j<=endyr;j++)	
			tempLenAtAge(k)+=LengthAtAge(j,k);
           tempLenAtAge = tempLenAtAge/(projectTimeVary);
       for(k=1;k<=maxAge;k++)	
		  for(j=endyr-projectTimeVary+1;j<=endyr;j++)	
			tempWgtAtAge(k)+=WeightAtAge(j,k);
           tempWgtAtAge = tempWgtAtAge/(projectTimeVary);		    
	}
 if(TimeVaryGrowthK<0  & TimeVaryLinf<0)
	{
	tempLenAtAge = LengthAtAge(endyr);
	tempWgtAtAge = WeightAtAge(endyr);
	}
	for(i=endyr+1;i<=endyr+Nproj;i++)
	{
	  LengthAtAge(i) = tempLenAtAge;
	  WeightAtAge(i) = tempWgtAtAge;
	}	  
  // logistic maturity curves
   for (j=1;j<=maxAge;j++)
    MatAtAge(j)       =  1./(1.+mfexp(-1.*log(19) *  (Ages(j)-maturity(1))/(maturity(2)-maturity(1))));
  // ==========================================================================
}

void model_parameters::getmortality(void)
{
  ofstream& post= *pad_post;
 int i,j;
 dvariable tempMproj;
 for (i=styr;i<=endyr;i++)
 {
	if(EstM>0)
		NatM(i) = mfexp(log_avg_NatM + m_dev(i));
	if(EstM<=0)
		NatM(i) = NatMin;
   fmort_dir(i)     	= mfexp(log_avg_fmort_dir    + fmort_dir_dev(i));
   F_dir(i) 			= sel_fish(i)*fmort_dir(i);
  }
   if(TimeVaryM>0 )
	{
	    for(j=endyr-projectTimeVary+1;j<=endyr;j++)	
			tempMproj+=NatM(j);
        tempMproj = tempMproj/(projectTimeVary);
	}
   if(TimeVaryM<0 )
   tempMproj = NatM(endyr);
  for(j=endyr+1;j<=(endyr+Nproj);j++)
	  NatM(j) = tempMproj;
 // ==========================================================================
}

void model_parameters::get_num_at_len(void)
{
  ofstream& post= *pad_post;
 NatAge.initialize();
 Spbio.initialize();
 predSurvBio.initialize();
 predCpueBio.initialize();
 predCatchAtAge.initialize();
 predCatchBio.initialize();
 NatAge(styr) = mfexp(stNatLen);
 for(ipass=styr;ipass<=endyr;ipass++) 
    get_num_at_len_yr();
	// ==========================================================================
}

void model_parameters::get_num_at_len_yr(void)
{
  ofstream& post= *pad_post;
 int i,j,k,sex;
 dvar_vector predSurvBioLen(1,maxAge);
 dvar_vector expRate(1,maxAge);
 dvariable temp;
 dvar_vector tempN(1,maxAge);
 dvar_vector dummyN(1,maxAge);
 dvar_matrix dummyAllProbs(1,LengthBinN,1,LengthBinN);
 dvar_vector nn(1,LengthBinN);
 dvariable pi;
 pi = 3.14159265;
  i = ipass;
  //==cpue index==
   dummyN				    = elem_prod(sel_fish(i),NatAge(i));
   predCpueBioLen(i)	= elem_prod(dummyN,WeightAtAge(i));
   for(j=1;j<=maxAge;j++)
    predCpueBio(i)		+= predCpueBioLen(i,j);
  // cout<<"predCpueBio"<<predCpueBio<<endl;
  //==survey==
  predSurvNatAge(i) 	= elem_prod(sel_srv,NatAge(i));
  predSurvBioLen		= elem_prod(predSurvNatAge(i),WeightAtAge(i));
   for(j=1;j<=maxAge;j++)
    predSurvBio(i)		+= predSurvBioLen(j);
  // cout<<"predSurvBio"<<predSurvBio<<endl;
   // cout<<"LengthAtAge"<<LengthAtAge<<endl;
   // cout<<"LengthBinsMid"<<LengthBinsMid<<endl;
   // cout<<"sel_srv"<<sel_srv<<endl;
   // cout<<"NatAge(i)"<<NatAge(i)<<endl;
   // cout<<"WeightAtAge"<<WeightAtAge<<endl;
  //==survey length frequencies==
	//find probability of length at age given variability
	for(k=1;k<=LengthBinN;k++)
    for(j=1;j<=LengthBinN;j++)
 	 dummyAllProbs(k,j) = (1/(GrowthSD*sqrt(2*pi)))*mfexp(-1*(square(LengthBinsMid(j)-LengthAtAge(i,k)))/(2*GrowthSD*GrowthSD));
	  // cout<<"dummyAllProbs"<<dummyAllProbs<<endl;
   //make a matrix for the N at length for each age given survey data
    nn = rowsum(dummyAllProbs);
 	for(k=1;k<=LengthBinN;k++)
	  dummyAllProbs(k) = (dummyAllProbs(k)/nn(k))*predSurvNatAge(i,k);
	  // cout<<"dummyAllProbs"<<dummyAllProbs<<endl;
  //make N at length for entire population	  
	  predSurvLenFreq(i) = colsum(dummyAllProbs);
	  // cout<<"predSurvLenFreq"<<predSurvLenFreq<<endl;
	  predSurvLenFreq(i) = predSurvLenFreq(i)/sum(predSurvLenFreq(i));
	  // cout<<"predSurvLenFreq"<<predSurvLenFreq<<endl;
 //==spawning biomass==
   dummyN					= elem_prod(NatAge(i),MatAtAge);
   tempN						= elem_prod(dummyN,WeightAtAge(i));
   Spbio(i) = 0;
   for(j=1;j<=maxAge;j++)
     Spbio(i)				+= tempN(j);
     // cout<<"Spbio"<<Spbio<<endl;
  // ==fishery==
  //==catch at length==
  Baranov = 0;
  for(j=1;j<=maxAge;j++)
  	Baranov(j) = (F_dir(i,j)/(NatM(i)+F_dir(i,j)))*(1-mfexp(-1.0*(F_dir(i,j)+NatM(i))));
    predCatchAtAge(i)	=  elem_prod(NatAge(i),Baranov);
  // cout<<"Selfish"<<sel_fish<<endl;
  // cout<<"predCatchAtAge(i)"<<predCatchAtAge(i)<<endl;
 //find probability of length at age for catch given variability
	for(k=1;k<=LengthBinN;k++)
    for(j=1;j<=LengthBinN;j++)
 	 dummyAllProbs(k,j) = (1/(GrowthSD*sqrt(2*pi)))*mfexp(-1*(square(LengthBinsMid(j)-LengthAtAge(i,k)))/(2*GrowthSD*GrowthSD));
  // cout<<"dummyAllProbs"<<dummyAllProbs<<endl;
   //make a matrix for the N at length for each age given survey data
    nn = rowsum(dummyAllProbs);
 	for(k=1;k<=LengthBinN;k++)
	  dummyAllProbs(k) = (dummyAllProbs(k)/nn(k))*predCatchAtAge(i,k);
  // cout<<"dummyAllProbs2"<<dummyAllProbs<<endl;
  //make N at length for entire population	  
	predCatchLenFreq(i) = colsum(dummyAllProbs);
  // cout<<"predCatchLenFreq(i)"<<predCatchLenFreq(i)<<endl;
    predCatchLenFreq(i) = predCatchLenFreq(i)/sum(predCatchLenFreq(i));
  // cout<<"predCatchLenFreq(i)"<<predCatchLenFreq(i)<<endl;
 //==aggregate catch biomass==	
    tempN					=	elem_prod(predCatchAtAge(i),WeightAtAge(i));	   
	predCatchBio(i) = 0;
  for(j=1;j<=maxAge;j++)
   predCatchBio(i)		+= tempN(j);
        // cout<<"predCatchBio"<<predCatchBio<<endl;
  // ==some die....but everyone else has a birthday!
  for(j =2;j<=(maxAge-1);j++)
   NatAge(i+1,j)			  = NatAge(i,j-1)*mfexp(-1*(NatM(i)+F_dir(i,j-1)));
  NatAge(i+1,maxAge)  = NatAge(i,maxAge)*mfexp(-1*(NatM(i)+F_dir(i,maxAge))) + NatAge(i,maxAge-1)*mfexp(-1*(NatM(i)+F_dir(i,maxAge-1)));
  if(i <= endyr)
   NatAge(i+1,1) 			= mfexp(mean_log_rec + rec_dev(i));
  else
   NatAge(i+1,1) 			= FutRec;
          // cout<<"NatAge"<<NatAge<<endl;
}

void model_parameters::evaluate_the_objective_function(void)
{
  ofstream& post= *pad_post;
 int i,j;
 dvariable tempVal;
 CpueBio_like.initialize();
 SurvBio_like.initialize();
 SurvLen_like.initialize();
 Catch_like.initialize();
 CatchLen_like.initialize();
 initsmo_penal.initialize();
 F_pen.initialize();
 Rec_pen.initialize();
 M_pen.initialize();
 Growth_pen1.initialize(); 
 Growth_pen2.initialize();
 Sel50_pen.initialize(); 
 Sel95_pen.initialize(); 
  //CPUE biomass
   for(i=styr+1;i<=endyr;i++)
	CpueBio_like +=  square(log(cpueIndex(i) + smallNum) - log(predCpueBio(i) + smallNum))  / ((log(cpueCV(i)*cpueCV(i)+ 1)  ));
 //Survey biomass
	for(i=styr;i<=endyr;i++)
	 SurvBio_like +=  square(log(survBiomass(i) + smallNum) - log(predSurvBio(i) + smallNum))  / ((log(survCV(i)*survCV(i)+ 1)  ));
    // cout<<"SurvBio"<<SurvBio_like<<endl;
  //catch biomass
   for(i=styr+1;i<=endyr;i++)
  	Catch_like +=  square(log(catchBiomass(i) + smallNum) - log(predCatchBio(i) + smallNum))  / ((log(catchCV(i)*catchCV(i)+ 1)  ));
    // cout<<"Catch_like"<<Catch_like<<endl;	
   // survey length frequencies
	for(i=styr+1;i<=endyr;i++)
	 for(j=1;j<=maxAge;j++)
	 {
      if(survLenFreq(i,j)>0.001)
	   SurvLen_like -= survSampN(i) *survLenFreq(i,j) *log(predSurvLenFreq(i,j) +smallNum);
	 }
    // cout<<"SurvLen_like"<<SurvLen_like<<endl;	
   //catch length frequencies
	for(i=styr+1;i<=endyr;i++)
	 for(j=1;j<=maxAge;j++)
	 {
      if(catchLenFreq(i,j)>0.001)
	   CatchLen_like -= catchSampN(i) *catchLenFreq(i,j) *log(predCatchLenFreq(i,j)+smallNum );		//same sample size for survey and catch right now...
      }  
    // cout<<"CatchLen_like"<<CatchLen_like<<endl;	
  initsmo_penal = InitSmoothWeight*(norm2(first_difference(stNatLen)));
  Rec_pen = RecruitPen*(norm2(first_difference(rec_dev)));
  F_pen = FmortPen*(norm2(first_difference(fmort_dir_dev)));
  M_pen = Mpenalty*(norm2(first_difference(m_dev)));
  Growth_pen1 = Growthpenalty*(norm2(first_difference(GrowthK_dev)));
  Growth_pen2= Growthpenalty*(norm2(first_difference(Linf_dev)));
  Sel50_pen = SelPenalty*(norm2(first_difference(SelPars_dev50)));
  Sel95_pen = SelPenalty*(norm2(first_difference(SelPars_dev95)));
  f = 0;
 if(FisheryIndependentData==1)
 {
  f += SurvBio_like;
  f += SurvLen_like;
 }
 f += CpueBio_like;
 f += Catch_like; 
 f += CatchLen_like;
 f += initsmo_penal;
 f += Rec_pen;
 f += F_pen;
 f += M_pen;
 f += Growth_pen1;
 f += Growth_pen2;
 f += Sel50_pen;
 f += Sel95_pen;
 cout << current_phase() << " " << call_no << " " << f << endl;
 // ==========================================================================
}

void model_parameters::get_fut_mortality(void)
{
  ofstream& post= *pad_post;
 int i,j;
  for (i=ipass;i<=endyr+100;i++)
  {
   if (IsB0 == 0)
   	fmort_dir(i) = 0; 
   else 
	fmort_dir(i) = FutMort; 
    F_dir(i) = sel_fish(i)*fmort_dir(i);
   }
}

void model_parameters::Find_F35(void)
{
  ofstream& post= *pad_post;
  dvariable Ratio,Target,B0;
  dvariable Btest;
  int icnt,i;
  dvar_vector tempN(1,maxAge);
  // Find B0 for given recruitment
  IsB0 = 0;
  FutRec = mfexp(mean_log_rec);
  FutRec = 100000;
  FutMort = 0;
  ipass = endyr+1;
  get_fut_mortality();
  for (ipass=endyr+1;ipass<=endyr+100;ipass++) get_num_at_len_yr(); 
 // put biomass in here soon
  B0 = Spbio(endyr +99);
  cout<<"B0"<<B0<<endl;
  Target = 0.35;
  IsB0 = 1;
  FutMort = 0.3;
  for (icnt=1;icnt<=20;icnt++)
   {
    ipass = endyr+1;
    get_fut_mortality();
    for (ipass=endyr+1;ipass<=endyr+100;ipass++) get_num_at_len_yr(); 
     Btest = Spbio(endyr +99);
   Ratio = Btest/B0;
    cout << FutMort << " " << Ratio << endl;
    FutMort = FutMort * Ratio / Target;
   }
  Bzero = B0;
  F35 = FutMort; 
  SBPRF35 = Btest/FutRec;
   cout << FutMort << " " << Ratio << endl; 
}

void model_parameters::Find_OFL(void)
{
  ofstream& post= *pad_post;
  dvariable Fmsy,Rbar,nn,alpha,beta;
  int BMSY_Yr1, BMSY_Yr2,ii,Iyr,kk,jj;
  if(HarvestControl==1)
  {
	FOFL = F35;
   get_fut_mortality();
   get_num_at_len_yr();
   OFL = 0;
   for(ii = 1;ii<=maxAge;ii++)
    OFL += predCatchAtAge(ipass,ii) *WeightAtAge(ipass,ii);
  }
   if(HarvestControl==2)
   OFL = ConstantCatch;
  if(HarvestControl==3)
  {
  FOFL = ConstantF;
   get_fut_mortality();
   get_num_at_len_yr();
   OFL = 0;
   for(ii = 1;ii<=maxAge;ii++)
    OFL += predCatchAtAge(ipass,ii) *WeightAtAge(ipass,ii);
  }
  if(HarvestControl==4)
  {
  BMSY_Yr1 = 1;BMSY_Yr2 = endyr-8;
  alpha = 0.05;
  beta = 0.25;
  // Define Fmsy
  Fmsy = F35;
 // Find Rbar (Dynamic or not)
  Rbar = 0; nn= 0;
  for (Iyr=BMSY_Yr1;Iyr<=BMSY_Yr2;Iyr++)
   {
    Rbar += mfexp(mean_log_rec+rec_dev(Iyr));
    nn += 1;
   }
  Rbar = Rbar / nn;
 // Specify the BMSY proxy
  Bmsy = SBPRF35 * Rbar;
  // cout << "Rbar" << Rbar << endl;
  // cout << "Bmsy" << Bmsy << endl;
  // cout<<"FMSY"<<Fmsy<<endl;
  // year for projection
  ipass = endyr+1;
  // Define future recruitment 
  if (ipass > endyr) FutRec = Rbar;
  // Find FOFL
  FutMort = Fmsy;
  get_fut_mortality();
  get_num_at_len_yr();
      // cout<<"Mspbio"<<Spbio(ipass)<<endl;
  if (Spbio(ipass) < Bmsy)
   {
    FutMort = 0;
    get_fut_mortality();
    get_num_at_len_yr();
    if (Spbio(ipass) > beta*Bmsy)
     {
      FutMort = Fmsy/2;
      get_fut_mortality();
      get_num_at_len_yr();
      for (ii=1;ii<=10;ii++)
       {
        FutMort = Fmsy*(Spbio(ipass)/Bmsy-alpha)/(1-alpha);
        get_fut_mortality();
        get_num_at_len_yr();
       }
     }
   }
   FOFL = FutMort;
   // cout<<"FOFL"<<FOFL<<endl;
   get_fut_mortality();
   get_num_at_len_yr();
   OFL = 0;
   for(ii = 1;ii<=maxAge;ii++)
    OFL += predCatchAtAge(ipass,ii) *WeightAtAge(ipass,ii);
  }
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
 int i,j;
  Find_F35();
  Find_OFL();
  report<<"SBPRF35"<<endl;
  report<<SBPRF35<<endl; 
  report<<"Curr Bio"<<endl;
  report<<Spbio(endyr)<<endl;
  report<<"B35"<<endl;
  report<<Bmsy<<endl;
  report<<"F35"<<endl;
  report<<F35<<endl;
  report<<"FOFL"<<endl;
  report<<FOFL<<endl;
  report<<"OFL"<<endl;
  report<<OFL<<endl;
  report<<"Bzero"<<endl;
  report<<Bzero<<endl;
 report << "#pred survey bio"<<endl;
 report << predSurvBio << endl;
 report << "#obs survey bio"<<endl;
 report << survBiomass << endl;
 report << "#pred catch bio"<<endl;
 report << predCatchBio << endl;
 report << "#obs catch bio"<<endl;
 report << catchBiomass << endl;
 report << "#pred cpue bio"<<endl;
 report << predCpueBio << endl;
 report << "#obs cpue bio"<<endl;
 report << cpueIndex << endl;
 report << "#pred surv len freq" <<endl;
    for(j=styr;j<=endyr;j++)
     report <<predSurvLenFreq(j) << endl; 
 report << "#obs surv len freq" <<endl;
    for(j=styr;j<=endyr;j++)
     report <<survLenFreq(j) << endl;  	 
 report << "#pred catch len freq" <<endl;
    for(j=styr;j<=endyr;j++)
     report <<predCatchLenFreq(j) << endl; 
 report << "#obs catch len freq" <<endl;
    for(j=styr;j<=endyr;j++)
     report <<catchLenFreq(j) << endl; 	 
  report<<"#spawning biomass"<< endl;
  report<<Spbio<<endl;
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{5000,5000,5000,5000,5000,5000,5000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{.1,.001,.001}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::final_calcs()
{
 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 if(FisheryIndependentData==1)
 {
 cout<<"Survey Bio"<<SurvBio_like<<endl;
 cout<<"Survey len"<<SurvLen_like<<endl;
 }
 cout<<"CPUE"<<CpueBio_like<<endl;
 cout<<"Catch"<<Catch_like<<endl;
 cout<<"Catch len"<<CatchLen_like<<endl;
 cout<<"initsmo_penal"<<initsmo_penal<<endl;
 cout<<"F_pen"<<F_pen<<endl; 
 cout<<"Rec_pen"<<Rec_pen<<endl; 
 cout<<"M_pen"<<M_pen<<endl; 
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_post;
  pad_post = NULL;
}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 4000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  // the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(900);
  
  time(&start);
  CheckFile.open("Check.Out");
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
