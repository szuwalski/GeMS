DATA_SECTION
//Generic assessment model

  !!CLASS ofstream post("eval.csv");
  int call_no;
  !! call_no = 0;

  // .DAT file
  init_int styr                                                    			 // start year of the model
  init_int endyr															//end year of the model
  init_int maxAge														//maximum age class
  init_vector Ages(1,maxAge)
  init_int cpueYr
  init_ivector cpueYears(1,cpueYr)
  init_int survYr
  init_ivector survYears(1,survYr)
  
  init_int catchLenYr
  init_ivector catchLenYears(1,survYr)
  init_number catchSampN
  init_int survLenYr
  init_ivector survLenYears(1,survYr)
  init_number survSampN
 !!cout<<"SurvLenYear"<<survLenYears<<endl;
 
  init_vector survBiomass(1,survYr)
  init_number survCV
  init_vector cpueIndex(1,cpueYr)
  init_number cpueCV
  init_vector catchBiomass(1,survYr)
  init_number catchCV
  init_vector LengthBinsMid(1,maxAge)
  !!cout<<"catchBiomass"<<catchBiomass<<endl; 
  
  init_matrix catchLenNum(styr,endyr,1,maxAge)	
  init_matrix survLenNum(styr,endyr,1,maxAge)			//insert when functionality exists


  // .CTL file
  !! ad_comm::change_datafile_name("SimAss.ctl");
  init_vector weightPars(1,2)
  init_vector lengthPars(1,3)
  init_number EstM
  init_number NatMin
  init_vector maturity(1,2)
  init_number steepness
  init_number smallNum
  init_number InitSmoothWeight
  init_number FisheryIndependentData
  init_number GrowthSD
  init_number FmortPen
  init_number RecruitPen
  init_number HarvestControl
  init_number ConstantCatch
  init_number ConstantF
  init_number HCalpha
  init_number HCbeta
  init_number LengthBinN
  !!cout<<"steep"<<steepness<<endl;
  
  int ipass;
  int Nproj;
  int IsB0;                                                         // Set to 0 for B0
  !! Nproj = 125;

// =======================================================================
PARAMETER_SECTION

	init_bounded_number fish_sel50(0.01,maxAge,2)
	init_bounded_number fish_sel95(1,maxAge,2)
	init_bounded_number srv_sel50(0.01,maxAge,2)
	init_bounded_number srv_sel95(1,maxAge,2)
	init_bounded_vector stNatLen(1,maxAge,-10,20,1)		//using this instead of recruitments because the number of years that would be required for long lived species would be a lot
	
	init_bounded_number log_avg_fmort_dir(-5,5,1)
	init_bounded_dev_vector fmort_dir_dev(styr,endyr,-20,20,1)

	init_bounded_number mean_log_rec(-5,25,1)
	init_bounded_dev_vector rec_dev(styr,endyr,-20,20,1)
	
	init_bounded_number NatM(0.001,0.8,EstM)
   
   // init_bounded_number dummy(1,2)
  
   //end of estimated parameters///
   objective_function_value f
 
  //storage for calculation of length frequencies
  matrix catchLenFreq(styr,endyr,1,maxAge)	
  matrix survLenFreq(styr,endyr,1,maxAge)	
  
  //Biological processes
  vector sel_fish(1,maxAge)											    // directed selectivity
  vector sel_srv(1,maxAge)            	                        		// Survey selectivity (do I need different eras?)
  vector MatAtAge(1,maxAge)
  vector LengthAtAge(1,maxAge)
   
  //fishing mortality and survival by area
  vector fmort_dir(styr,endyr+Nproj)								//Directed fully selected F 
  matrix F_dir(styr,endyr+Nproj,1,maxAge)
  matrix F_trawl(styr,endyr+Nproj,1,maxAge)
   vector Baranov(1,maxAge)
   
  //Numbers at length by sex (ladies first)
  matrix NatAge(styr,endyr+Nproj,1,maxAge)
  sdreport_vector Spbio(styr,endyr+Nproj)

 //predicted total numbers in the survey
   vector predSurvBio(styr,endyr+Nproj)
   vector predCpueBio(styr,endyr+Nproj)
   matrix predCpueBioLen(styr,endyr+Nproj,1,maxAge)
   matrix predSurvNatAge(styr,endyr+Nproj,1,maxAge)
   matrix predSurvAgeFreq(styr,endyr+Nproj,1,maxAge)
   matrix predSurvLenFreq(styr,endyr+Nproj,1,maxAge)   
   
  //predicted catches, discards and bycatches at length
  matrix predCatchAtAge(styr,endyr+Nproj,1,maxAge)
  sdreport_vector predCatchBio(styr,endyr+Nproj);  
  matrix predCatchAgeFreq(styr,endyr+Nproj,1,maxAge)
  matrix predCatchLenFreq(styr,endyr+Nproj,1,maxAge)
  
  //likelihood components
  number CpueBio_like
  number SurvBio_like
  number SurvLen_like
  number Catch_like
  number CatchLen_like
  number initsmo_penal
  number Rec_pen
  number F_pen
 
  //weight at length
  vector WeightAtAge(1,maxAge)
  
  //reference points
  number F35;                                                       // F35
  number SBPRF35;                                                   // SBPR(F=F35);
  number Bmsy;                                                      // Bmsy
  number FOFL;														//FOFL
  number OFL;
  number FutRec;                                                     // Future recruitment
  number FutMort;                                                    // Future target
  number Bzero;

//  ========================================================================
PRELIMINARY_CALCS_SECTION
 int i,j;
 dvariable nn;
  // for(i=1;i<=maxAge;i++)
   // WeightAtAge(i) = weightPars(1) * pow(LengthBinsMid(i),weightPars(2));
   
  
//make length frequencies from observed
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

// =============================================================
INITIALIZATION_SECTION
 fish_sel50 5.916511
 fish_sel95 7.831472
 srv_sel50 2.683375
 srv_sel95 4.365736
 NatM NatMin
// ==========================================================================
PROCEDURE_SECTION
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
  
// ==========================================================================
FUNCTION getselectivity
 int j;

  // logistic selectivity curves in terms of age
    for (j=1;j<=maxAge;j++)
     { 
	  sel_fish(j)       =  1./(1.+mfexp(-1.*log(19.)*(Ages(j)-fish_sel50)/(fish_sel95-fish_sel50)));
      sel_srv(j)		   =  1./(1.+mfexp(-1.*log(19.)*(Ages(j)-srv_sel50)/(srv_sel95-srv_sel50)));
	 } 

// ==========================================================================
FUNCTION getmaturity
 int j;

  // logistic maturity curves
   for (j=1;j<=maxAge;j++)
   {
    MatAtAge(j)       =  1./(1.+mfexp(-1.*log(19) *  (Ages(j)-maturity(1))/(maturity(2)-maturity(1))));
    LengthAtAge(j) = lengthPars(2)*(1-mfexp(-1*lengthPars(1)*(Ages(j)-lengthPars(3))));
	WeightAtAge(j) = weightPars(1) * pow(LengthAtAge(j),weightPars(2));
   }

  // ==========================================================================
FUNCTION getmortality
 int i,j;
 for (i=styr;i<=endyr;i++)
 {
   fmort_dir(i)     	= mfexp(log_avg_fmort_dir    + fmort_dir_dev(i));
   F_dir(i) 			= sel_fish*fmort_dir(i);
  }

 // ==========================================================================
FUNCTION get_num_at_len
//to estimate an initial recruitment in the start year?
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
FUNCTION get_num_at_len_yr
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
   dummyN				    = elem_prod(sel_fish,NatAge(i));
   predCpueBioLen(i)	= elem_prod(dummyN,WeightAtAge);
   for(j=1;j<=maxAge;j++)
    predCpueBio(i)		+= predCpueBioLen(i,j);
  // cout<<"predCpueBio"<<predCpueBio<<endl;
  
  //==survey==
  predSurvNatAge(i) 	= elem_prod(sel_srv,NatAge(i));
  predSurvBioLen		= elem_prod(predSurvNatAge(i),WeightAtAge);
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
 	 dummyAllProbs(k,j) = (1/(GrowthSD*sqrt(2*pi)))*mfexp(-1*(square(LengthBinsMid(j)-LengthAtAge(k)))/(2*GrowthSD*GrowthSD));
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
   tempN						= elem_prod(dummyN,WeightAtAge);
   Spbio(i) = 0;
   for(j=1;j<=maxAge;j++)
     Spbio(i)				+= tempN(j);
     // cout<<"Spbio"<<Spbio<<endl;
	 
  // ==fishery==
  //==catch at length==
  Baranov = 0;
  for(j=1;j<=maxAge;j++)
  	Baranov(j) = (F_dir(i,j)/(NatM+F_dir(i,j)))*(1-mfexp(-1.0*(F_dir(i,j)+NatM)));
    predCatchAtAge(i)	=  elem_prod(NatAge(i),Baranov);
  // cout<<"Selfish"<<sel_fish<<endl;
  // cout<<"predCatchAtAge(i)"<<predCatchAtAge(i)<<endl;

 //find probability of length at age for catch given variability
	for(k=1;k<=LengthBinN;k++)
    for(j=1;j<=LengthBinN;j++)
 	 dummyAllProbs(k,j) = (1/(GrowthSD*sqrt(2*pi)))*mfexp(-1*(square(LengthBinsMid(j)-LengthAtAge(k)))/(2*GrowthSD*GrowthSD));
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
    tempN					=	elem_prod(predCatchAtAge(i),WeightAtAge);	   
	predCatchBio(i) = 0;
  for(j=1;j<=maxAge;j++)
   predCatchBio(i)		+= tempN(j);
        // cout<<"predCatchBio"<<predCatchBio<<endl;
		
  // ==some die....but everyone else has a birthday!
  for(j =2;j<=(maxAge-1);j++)
   NatAge(i+1,j)			  = NatAge(i,j-1)*mfexp(-1*(NatM+F_dir(i,j-1)));
  NatAge(i+1,maxAge)  = NatAge(i,maxAge)*mfexp(-1*(NatM+F_dir(i,maxAge))) + NatAge(i,maxAge-1)*mfexp(-1*(NatM+F_dir(i,maxAge-1)));

  if(i <= endyr)
   NatAge(i+1,1) 			= mfexp(mean_log_rec + rec_dev(i));
  else
   NatAge(i+1,1) 			= FutRec;
          // cout<<"NatAge"<<NatAge<<endl;
// ==========================================================================                                                             
FUNCTION evaluate_the_objective_function
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

  //CPUE biomass
   for(i=styr+1;i<=endyr;i++)
	CpueBio_like +=  square(log(cpueIndex(i) + smallNum) - log(predCpueBio(i) + smallNum))  / ((log(survCV*survCV+ 1)  ));
 
 //Survey biomass
	for(i=styr;i<=endyr;i++)
	 SurvBio_like +=  square(log(survBiomass(i) + smallNum) - log(predSurvBio(i) + smallNum))  / ((log(survCV*survCV+ 1)  ));

    // cout<<"SurvBio"<<SurvBio_like<<endl;
 
  //catch biomass
   for(i=styr+1;i<=endyr;i++)
  	Catch_like +=  square(log(catchBiomass(i) + smallNum) - log(predCatchBio(i) + smallNum))  / ((log(catchCV*catchCV+ 1)  ));
    // cout<<"Catch_like"<<Catch_like<<endl;	
	
   // survey length frequencies
	for(i=styr+1;i<=endyr;i++)
	 for(j=1;j<=maxAge;j++)
	 {
      if(survLenFreq(i,j)>0.001)
	   SurvLen_like -= survSampN *survLenFreq(i,j) *log(predSurvLenFreq(i,j) +smallNum);
	 }
    // cout<<"SurvLen_like"<<SurvLen_like<<endl;	
	
   //catch length frequencies
	for(i=styr+1;i<=endyr;i++)
	 for(j=1;j<=maxAge;j++)
	 {
      if(catchLenFreq(i,j)>0.001)
	   CatchLen_like -= survSampN *catchLenFreq(i,j) *log(predCatchLenFreq(i,j)+smallNum );		//same sample size for survey and catch right now...
      }  
    // cout<<"CatchLen_like"<<CatchLen_like<<endl;	
	
  initsmo_penal = InitSmoothWeight*(norm2(first_difference(stNatLen)));
  Rec_pen = RecruitPen*(norm2(first_difference(rec_dev)));
  F_pen = FmortPen*(norm2(first_difference(fmort_dir_dev)));
  
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

 cout << current_phase() << " " << call_no << " " << f << endl;
 
 // ==========================================================================

FUNCTION get_fut_mortality
 int i,j;

  for (i=ipass;i<=endyr+100;i++)
  {
   if (IsB0 == 0)
   	fmort_dir(i) = 0; 
   else 
	fmort_dir(i) = FutMort; 

    F_dir(i) = sel_fish*fmort_dir(i);
   }
 
// ==========================================================================

FUNCTION Find_F35
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
  
// Find F35%  
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

// ==========================================================================

FUNCTION Find_OFL
  dvariable Fmsy,Rbar,nn,alpha,beta;
  int BMSY_Yr1, BMSY_Yr2,ii,Iyr,kk,jj;
 
  if(HarvestControl==1)
  {
	FOFL = F35;
   get_fut_mortality();
   get_num_at_len_yr();
   OFL = 0;
   for(ii = 1;ii<=maxAge;ii++)
    OFL += predCatchAtAge(ipass,ii) *WeightAtAge(ii);
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
    OFL += predCatchAtAge(ipass,ii) *WeightAtAge(ii);
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
    OFL += predCatchAtAge(ipass,ii) *WeightAtAge(ii);
  }
  
// ==========================================================

REPORT_SECTION
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

 
//===============================================================================

GLOBALS_SECTION
 #include <math.h>
 #include <admodel.h>
  #include <time.h>
 
 ofstream CheckFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

// EVALUATE THE GAMMA FUNCTION
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

// ===============================================================================
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 5000,5000,5000,5000,5000,5000,5000
//  convergence_criteria 1,1,1,1,.01,.001,1e-3,1e-3
 convergence_criteria .1,.001,.001

TOP_OF_MAIN_SECTION
  arrmblsize = 4000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  // the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(900);
  
  time(&start);
  CheckFile.open("Check.Out");

FINAL_SECTION
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
 
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;

