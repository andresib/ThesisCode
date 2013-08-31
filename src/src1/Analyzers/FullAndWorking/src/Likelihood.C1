#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>


float byHandc(double x){
  //  double p=exp(-(discriminant - 0.25)/(2*(0.082)*(0.082)));
  float p=185.268 + 23746.5*x - 294140*x**2 + 1478110*x**3 - 3991330*x**4 + 6283470*x**5 - 5783640*x**6 +2889290*x**7 - 605695*x**8; 
  return p;
}
float byHandb(double x){
  //  double p=exp(-(discriminant - 0.25)/(2*(0.082)*(0.082)));
  float p=62.3237 + 2392.79*x - 22118.5*x**2 + 72492.7*x**3 - 77981.6*x**4 - 86901*x**5 + 294793*x**6 - 266675*x**7 + 84802.5*x**8;
  return p;
}

void Likelihood()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111111);
  gStyle->SetPalette(60);
  gStyle->SetOptFit(1);
  bool bTaggingLikelihood;
  bool madWeightLikelihood_Gen;
  bool madWeightLikelihood_RecoWithoutTF;
  bool madWeightLikelihood_Reco;
  bool bTagging_MWLikelihood;
  bool madWeightLikelihoodRecoD;
  bool bTag_madWeightLikelihood_RecoD;
  bool MC;
  MC=true; 
  bTaggingLikelihood=false;
  madWeightLikelihood_Gen=true;
  madWeightLikelihood_RecoWithoutTF=false;
  madWeightLikelihood_Reco=true;
  madWeightLikelihood_RecoD=true;
  bTag_madWeightLikelihood_RecoD=true;
  bTagging_MWLikelihoodGen=false;
  bTagging_MWLikelihoodReco=false;
  bTagging_MWLikelihoodRecoWithoutTF=false;
  
  TFile histos("LikelihoodBtag.root","update");

  TH1F *Reco_Matching=new TH1F("1: Matched Events, 3: Right Matching" , "1: Matched Events, 3: Right Matching" , 5 , 0 ,5 );

  TH1F *MatchingWithPreSelection=new TH1F("MatchingWithPreSelection" , "MatchingWithPreSelection" , 101 , 0 ,100 );
  TH1F *Btag_MatchingWithPreSelection=new TH1F("Btag_MatchingWithPreSelection" , "Btag_MatchingWithPreSelection" , 101 , 0 , 100 );
  TH1F *bTaggingLikelihood_Total=new TH1F("bTagging likelihood (all)" , "bTagging likelihood (all)" , 300 , -0.1 , 0.5 );
  TH1F *bTaggingLikelihood_Real =new TH1F("bTagging likelihood (real)" , "bTagging likelihood (real)" , 300 , -0.1 , 0.5 );
  TH1F *bTaggingLikelihood_Best =new TH1F("bTagging likelihood (best)" , "bTagging likelihood (best)" , 300 , -0.1 , 0.5 );
  TH1F *bTaggingLikelihood_NumberOfJetsIdentified =new TH1F("bTagging Likelihood: Number of jets that were identified" , "bTagging Likelihood: Number of jets that were identified" , 11 , -1 , 10 );
  TH1F *bTaggingLikelihood_Difference1_2 =new TH1F("bTagging Likelihood: Absolute difference between the two different highest likelihooks" , "bTagging Likelihood: Absolute difference between the two highest likelihooks" , 5000 , 0 , 1 );
  TH1F *bTaggingLikelihood_DifferenceReal_1 =new TH1F("bTagging Likelihood: Difference between the real likelihood and the highest one" , "bTagging Likelihood: Difference between the real likelihood and the highest one" , 5000 , -1 , 0 ); 
  
  
  
  
  TH1F *MWLikelihood_Gen_Total=new TH1F("MW likelihood (all) - Gen" , "MW likelihood (all) - Gen" , 3000 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_Gen_Real =new TH1F("MW likelihood (real) - Gen" , "MW likelihood (real) - Gen" , 3000 ,1e-33 , 1e-24 );
  TH1F *MWLikelihood_Gen_Best =new TH1F("MW likelihood (best) - Gen" , "MW likelihood (best) - Gen" , 251 , 0 , 250 );
  TH1F *MWLikelihood_Gen_NumberOfJetsIdentified =new TH1F("MW Likelihood - Gen: Number of jets that were identified" , "MW Likelihood - Gen: Number of jets that were identified" , 11 , -1 , 10 );
  TH1F *MWLikelihood_Gen_JetsIdentified =new TH1F("MW Likelihood - Gen: Jets that were identified" , "MW Likelihood - Gen: Jets that were identified" , 11 , -1 , 10 );
  TH1F *MWLikelihood_Gen_Difference1_2 =new TH1F("MW Likelihood - Gen: Absolute difference between the two different highest likelihooks" , "MW Likelihood - Gen: Absolute difference between the two highest likelihooks" , 3000 , 0 , 1e-24 );
  TH1F *MWLikelihood_Gen_DifferenceReal_1 =new TH1F("MW Likelihood - Gen: Difference between the real likelihood and the highest one" , "MW Likelihood - Gen: Difference between the real likelihood and the highest one" , 3000 , 0 , 1e-24 ); 
  
  
  
  
  
  TH1F *MWLikelihood_RecoWithoutTF_Total=new TH1F("MW likelihood (all) - Reco Without TF" , "MW likelihood (all) - Reco Without TF" , 300 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_RecoWithoutTF_Real =new TH1F("MW likelihood (real) - Reco Without TF" , "MW likelihood (real) - Reco Without TF" , 300 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_RecoWithoutTF_Best =new TH1F("MW likelihood (best) - Reco Without TF" , "MW likelihood (best) - Reco Without TF" , 300 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_RecoWithoutTF_NumberOfJetsIdentified =new TH1F("MW Likelihood - Reco Without TF: Number of jets that were identified" , "MW Likelihood - Reco Without TF: Number of jets that were identified" , 11 , -1 , 10 );
  TH1F *MWLikelihood_RecoWithoutTF_JetsIdentified =new TH1F("MW Likelihood - Reco Without TF: Jets that were identified" , "MW Likelihood - Reco Without TF: Jets that were identified" , 11 , -1 , 10 );
  TH1F *MWLikelihood_RecoWithoutTF_Difference1_2 =new TH1F("MW Likelihood - Reco Without TF: Absolute difference between the two different highest likelihooks" , "MW Likelihood - Reco Without TF: Absolute difference between the two highest likelihooks" , 300 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_RecoWithoutTF_DifferenceReal_1 =new TH1F("MW Likelihood - Reco Without TF: Difference between the real likelihood and the highest one" , "MW Likelihood - Reco Without TF: Difference between the real likelihood and the highest one" , 300 , 1e-33 , 1e-24 ); 
  
  
  
  
  
  
  TH1F *MWLikelihood_Reco_Total=new TH1F("MW likelihood (all) - Reco" , "MW likelihood (all) - Reco" , 3000 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_Reco_Real =new TH1F("MW likelihood (real) - Reco" , "MW likelihood (real) - Reco" , 3000 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_Reco_Best =new TH1F("MW likelihood (best) - Reco" , "MW likelihood (best) - Reco" , 251 , 0 , 250 );
  TH1F *MWLikelihood_Reco_NumberOfJetsIdentified =new TH1F("MW Likelihood - Reco: Number of jets that were identified" , "MW Likelihood - Reco: Number of jets that were identified" , 11 , -1 , 10 );
  TH1F *MWLikelihood_Reco_JetsIdentified =new TH1F("MW Likelihood - Reco: Jets that were identified" , "MW Likelihood - Reco: Jets that were identified" , 11 , -1 , 10 );
  TH1F *MWLikelihood_Reco_Difference1_2 =new TH1F("MW Likelihood - Reco: Absolute difference between the two different highest likelihooks" , "MW Likelihood - Reco: Absolute difference between the two highest likelihooks" , 3000 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_Reco_DifferenceReal_1 =new TH1F("MW Likelihood - Reco: Difference between the real likelihood and the highest one" , "MW Likelihood - Reco: Difference between the real likelihood and the highest one" , 3000 , 1e-33 , 1e-24 ); 
  

  TH1F *MWLikelihood_RecoD_Total=new TH1F("MW likelihood (all) - RecoD" , "MW likelihood (all) - RecoD" , 3000 , 1e-33 , 1e-24 );
  TH1F *MWLikelihood_RecoD_Best =new TH1F("MW likelihood (best) - RecoD" , "MW likelihood (best) - RecoD" , 251 , 0 ,250 );
  TH1F *MWLikelihood_RecoD_Difference1_2 =new TH1F("MW Likelihood - RecoD: Absolute difference between the two different highest likelihooks" , "MW Likelihood - RecoD: Absolute difference between the two highest likelihooks" , 3000 , 1e-33 , 1e-24 );



  TH1F *Btag_MWLikelihood_RecoD_Total=new TH1F("Btag_MW likelihood (all) - RecoD" , "Btag_MW likelihood (all) - RecoD" , 3000 , 1e-33 , 1e-24 );
  TH1F *Btag_MWLikelihood_RecoD_Best =new TH1F("Btag_MW likelihood (best) - RecoD" , "Btag_MW likelihood (best) - RecoD" , 251 , 0 , 250 );
  TH1F *Btag_MWLikelihood_RecoD_Difference1_2 =new TH1F("Btag_MW Likelihood - RecoD: Absolute difference between the two different highest likelihooks" , "Btag_MW Likelihood - RecoD: Absolute difference between the two highest likelihooks" , 3000 , 1e-33 , 1e-24 );




  TH1F *BTaggingMWLikelihood_Gen_Total=new TH1F("BTaggingMW likelihood (all) - Gen" , "BTaggingMW likelihood (all) - Gen" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Gen_Real =new TH1F("BTaggingMW likelihood (real) - Gen" , "BTaggingMW likelihood (real) - Gen" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Gen_Best =new TH1F("BTaggingMW likelihood (best) - Gen" , "BTaggingMW likelihood (best) - Gen" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Gen_NumberOfJetsIdentified =new TH1F("BTaggingMW Likelihood - Gen: Number of jets that were identified" , "BTaggingMW Likelihood - Gen: Number of jets that were identified" , 11 , -1, 10 );
  TH1F *BTaggingMWLikelihood_Gen_JetsIdentified =new TH1F("BTaggingMW Likelihood - Gen: Jets that were identified" , "BTaggingMW Likelihood - Gen: Jets thatwere identified" , 11 , -1 , 10 );
  TH1F *BTaggingMWLikelihood_Gen_Difference1_2 =new TH1F("BTaggingMW Likelihood - Gen: Absolute difference between the two different highest likelihooks" , "BTaggingMW Likelihood - Gen: Absolute difference between the two highest likelihooks" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Gen_DifferenceReal_1 =new TH1F("BTaggingMW Likelihood - Gen: Difference between the real likelihood and the highest one" , "BTaggingMW Likelihood - Gen: Difference between the real likelihood and the highest one" , 300 , 1e-33 , 1e-24 );







  TH1F *BTaggingMWLikelihood_RecoWithoutTF_Total=new TH1F("BTaggingMW likelihood (all) - Reco Without TF" , "BTaggingMW likelihood (all) - Reco Without TF" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_Real =new TH1F("BTaggingMW likelihood (real) - Reco Without TF" , "BTaggingMW likelihood (real) - Reco Without TF" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_Best =new TH1F("BTaggingMW likelihood (best) - Reco Without TF" , "BTaggingMW likelihood (best) - Reco Without TF" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_NumberOfJetsIdentified =new TH1F("BTaggingMW Likelihood - Reco Without TF: Number of jets that were identified" , "BTaggingMW Likelihood - Reco Without TF: Number of jets that were identified" , 11 , -1, 10 );
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_JetsIdentified =new TH1F("BTaggingMW Likelihood - Reco Without TF: Jets that were identified" , "BTaggingMW Likelihood - Reco Without TF: Jets that were identified" , 11 , -1 , 10 );
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_Difference1_2 =new TH1F("BTaggingMW Likelihood - Reco Without TF: Absolute difference between the two different highest likelihooks" , "BTaggingMW Likelihood - Reco Without TF: Absolute difference between the two highest likelihooks" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_DifferenceReal_1 =new TH1F("BTaggingMW Likelihood - Reco Without TF: Difference between the real likelihood and the highest one" , "BTaggingMW Likelihood - Reco Without TF: Difference between the real likelihood and the highest one" , 300 , 1e-33 , 1e-24 );


  TH1F *BTaggingMWLikelihood_Reco_Total=new TH1F("BTaggingMW likelihood (all) - Reco" , "BTaggingMW likelihood (all) - Reco" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Reco_Real =new TH1F("BTaggingMW likelihood (real) - Reco" , "BTaggingMW likelihood (real) - Reco" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Reco_Best =new TH1F("BTaggingMW likelihood (best) - Reco" , "BTaggingMW likelihood (best) - Reco" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Reco_NumberOfJetsIdentified =new TH1F("BTaggingMW Likelihood - Reco: Number of jets that were identified" , "BTaggingMW Likelihood - Reco: Number of jets that were identified" , 11 , -1, 10 );
  TH1F *BTaggingMWLikelihood_Reco_JetsIdentified =new TH1F("BTaggingMW Likelihood - Reco: Jets that were identified" , "BTaggingMW Likelihood - Reco: Jets that were identified" , 11 , -1 , 10 );
  TH1F *BTaggingMWLikelihood_Reco_Difference1_2 =new TH1F("BTaggingMW Likelihood - Reco: Absolute difference between the two different highest likelihooks" , "BTaggingMW Likelihood - Reco: Absolute difference between the two highest likelihooks" , 300 , 1e-33 , 1e-24 );
  TH1F *BTaggingMWLikelihood_Reco_DifferenceReal_1 =new TH1F("BTaggingMW Likelihood - Reco: Difference between the real likelihood and the highest one" , "BTaggingMW Likelihood - Reco: Difference between the real likelihood and the highest one" , 300 , 1e-33 , 1e-24 );

  TH1F *BTaggingMWLikelihood_Gen_Inv_Mass_t_ =new TH1F("BTaggingMW Likelihood - Gen: Inv. Mass t_" , "BTaggingMW Likelihood - Gen: Inv. Mass t_" , 3000 , 100 , 400 );
  TH1F *BTaggingMWLikelihood_Reco_Inv_Mass_t_ =new TH1F("BTaggingMW Likelihood - Reco: Inv. Mass t_" , "BTaggingMW Likelihood - Reco: Inv. Mass t_" , 3000 , 100 , 400 );
  TH1F *BTaggingMWLikelihood_Gen_Inv_Mass_t =new TH1F("BTaggingMW Likelihood - Gen: Inv. Mass t" , "BTaggingMW Likelihood - Gen: Inv. Mass t" , 3000 , 100 , 400 );
  TH1F *BTaggingMWLikelihood_Reco_Inv_Mass_t =new TH1F("BTaggingMW Likelihood - Reco: Inv. Mass t" , "BTaggingMW Likelihood - Reco: Inv. Mass t" , 3000 , 100 , 400 );

  TH1F *MWLikelihood_Gen_Inv_Mass_t_ =new TH1F("MW Likelihood - Gen: Inv. Mass t_" , "MW Likelihood - Gen: Inv. Mass t_" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_Reco_Inv_Mass_t_ =new TH1F("MW Likelihood - Reco: Inv. Mass t_" , "MW Likelihood - Reco: Inv. Mass t_" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_RecoD_Inv_Mass_t_ =new TH1F("MW Likelihood - RecoD: Inv. Mass t_" , "MW Likelihood - RecoD: Inv. Mass t_" ,100 , 0 , 1000 );  
  TH1F *Btag_MWLikelihood_RecoD_Inv_Mass_t_ =new TH1F("Btag_MW Likelihood - RecoD: Inv. Mass t_" , "Btag_MW Likelihood - RecoD: Inv. Mass t_" ,100 , 0 , 1000 );  
  

  TH1F *MWLikelihood_Gen_Inv_Mass_t__b =new TH1F("MW Likelihood - Gen: Inv. Mass t_ - b" , "MW Likelihood - Gen: Inv. Mass t_ - b" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_Reco_Inv_Mass_t__b =new TH1F("MW Likelihood - Reco: Inv. Mass t_ - b" , "MW Likelihood - Reco: Inv. Mass t_ - b" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_RecoD_Inv_Mass_t__b =new TH1F("MW Likelihood - RecoD: Inv. Mass t_ - b" , "MW Likelihood - RecoD: Inv. Mass t_ - b" ,100 , 0 , 1000 );  
  TH1F *Btag_MWLikelihood_RecoD_Inv_Mass_t__b =new TH1F("Btag_MW Likelihood - RecoD: Inv. Mass t_ - b" , "Btag_MW Likelihood - RecoD: Inv. Mass t_ - b" ,100 , 0 , 1000 );  
  




  TH1F *MWLikelihood_Gen_Inv_Mass_t =new TH1F("MW Likelihood - Gen: Inv. Mass t" , "MW Likelihood - Gen: Inv. Mass t" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_Reco_Inv_Mass_t =new TH1F("MW Likelihood - Reco: Inv. Mass t" , "MW Likelihood - Reco: Inv. Mass t" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_RecoD_Inv_Mass_t =new TH1F("MW Likelihood - RecoD: Inv. Mass t" , "MW Likelihood - RecoD: Inv. Mass t" , 100 , 0 , 1000 );
  TH1F *Btag_MWLikelihood_RecoD_Inv_Mass_t =new TH1F("Btag_MW Likelihood - RecoD: Inv. Mass t" , "Btag_MW Likelihood - RecoD: Inv. Mass t" , 100 , 0 , 1000 );


  TH1F *MWLikelihood_Gen_Inv_Mass_t_b =new TH1F("MW Likelihood - Gen: Inv. Mass t - b" , "MW Likelihood - Gen: Inv. Mass t - b" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_Reco_Inv_Mass_t_b =new TH1F("MW Likelihood - Reco: Inv. Mass t - b" , "MW Likelihood - Reco: Inv. Mass t - b" , 100 , 0 , 1000 );
  TH1F *MWLikelihood_RecoD_Inv_Mass_t_b =new TH1F("MW Likelihood - RecoD: Inv. Mass t - b" , "MW Likelihood - RecoD: Inv. Mass t - b" , 100 , 0 , 1000 );
  TH1F *Btag_MWLikelihood_RecoD_Inv_Mass_t_b =new TH1F("Btag_MW Likelihood - RecoD: Inv. Mass t - b" , "Btag_MW Likelihood - RecoD: Inv. Mass t - b" , 100 , 0 , 1000 );


  TH1F *MWLikelihood_WOReco_Inv_Mass_t_ =new TH1F("MW Likelihood - Reco Without TF: Inv. Mass t_" , "MW Likelihood - Reco Without TF: Inv. Mass t_" , 3000 , 100 , 400 );
  TH1F *MWLikelihood_WOReco_Inv_Mass_t =new TH1F("MW Likelihood - Reco Without TF: Inv. Mass t" , "MW Likelihood - Reco Without TF: Inv. Mass t" , 3000 , 100 , 400);
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_Inv_Mass_t_ =new TH1F("BTaggingMW Likelihood - Reco Without TF: Inv. Mass t_" , "BTaggingMW Likelihood - Reco Without TF: Inv. Mass t_" , 3000 , 100 , 400 );
  TH1F *BTaggingMWLikelihood_RecoWithoutTF_Inv_Mass_t =new TH1F("BTaggingMW Likelihood - Reco Without TF: Inv. Mass t" , "BTaggingMW Likelihood - Reco Without TF: Inv. Mass t" , 3000 , 100 , 400 );
 





  if (madWeightLikelihood_RecoD)
    {
      ifstream list;
      list.open("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/MadWeight/Results/background/top/recoD/list.txt");
      string idD; 
      list>>idD;
      cout<<idD<<endl;
      bool noMore = false;
      while (!list.eof())
	{
   	  float  error, likelihood,likelihood1, likelihood2, realLikelihood;
          int numberOfPermutation1, numberOfPermutation2;
	  
	  
          likelihood1=0;
          likelihood2=0;
	  

          ifstream mwlikelihood;
          const char * file= ("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/MadWeight/Results/background/top/recoD/"+idD).c_str();
          cout<<file<<endl;
	  mwlikelihood.open(file);
          //string line;
	  //getline(mwlikelihood,line);
          //cout<<line<<endl;
	  string tmp1;
          mwlikelihood>>tmp1;
          cout<<"JP"<<tmp1<<endl;
	  mwlikelihood>>likelihood>>error;
          cout<<"likelihood: "<<likelihood<<endl;
          ifstream LHCO;
          LHCO.open(("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/LHCO/background/top/recoD/"+idD).c_str());
          double tmp_;
	  double jeta[4];
	  double jphi[4];
	  double jpt[4];
	  double jm[4];
	  double eeta,ephi,ept,em;
	  double meta,mphi,mpt,mm;
          double tmp_B;
	  double jetaB[4];
	  double jphiB[4];
	  double jptB[4];
	  double jmB[4];
	  double eetaB,ephiB,eptB,emB;
	  double metaB,mphiB,mptB,mmB;
          LHCO>>tmp_>>tmp_;
          cout<<tmp_<<endl;
          LHCO>>tmp_;
          LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  
	  const char init = tmp1[0];
          int combination = 0;
          int bestCombination = 0;
          while (init=='1')
            {
              combination++;
              MWLikelihood_RecoD_Total->Fill(likelihood);
	      if (likelihood>likelihood1)
		{
		  likelihood2 = likelihood1;
		  likelihood1 = likelihood;
                  bestCombination = combination;
                  for (int k=0 ; k<4 ; k++)
                  {
                    jetaB[k]=jeta[k];
                    jphiB[k]=jphi[k];
                    jptB[k]=jpt[k];
                    jmB[k]=jm[k];
                    cout<<"entro"<<endl;
                  } 
                  eetaB=eeta;
                  ephiB=ephi;
                  eptB=ept;
                  emB=em;
                  metaB=meta;
                  mphiB=mphi;
                  mptB=mpt;
                  mmB=mm;
		}
              LHCO>>tmp_>>tmp_;
              cout<<tmp_<<endl;
              LHCO>>tmp_;
              LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      
              mwlikelihood>>tmp1>>likelihood>>error;
              const char init = tmp1[0];
	    }
	  mwlikelihood.close();
	  MatchingWithPreSelection->Fill(bestCombination);
          likelihood1 =  -log(likelihood1);
	  MWLikelihood_RecoD_Best->Fill(likelihood1);
	  float _difference1_2 = likelihood1 - likelihood2;
	  MWLikelihood_RecoWithoutTF_Difference1_2->Fill(_difference1_2);
	  LHCO.close();
	  vector<TLorentzVector> j(4);
	  TLorentzVector e;
	  TLorentzVector m;
          j[0].SetPtEtaPhiM(jptB[0],jetaB[0],jphiB[0],jmB[0]);
	  j[1].SetPtEtaPhiM(jptB[1],jetaB[1],jphiB[1],jmB[1]);
	  j[2].SetPtEtaPhiM(jptB[2],jetaB[2],jphiB[2],jmB[2]);
	  j[3].SetPtEtaPhiM(jptB[3],jetaB[3],jphiB[3],jmB[3]);
	  e.SetPtEtaPhiM(eptB,eetaB,ephiB,emB);
	  m.SetPtEtaPhiM(mptB,metaB,mphiB,mmB);
	  double invmasst, invmasst_;
	  TLorentzVector t = j[2] + j[3] + j[0];
          invmasst = t.M();
	  invmasst_=(e+m+j[1]).M(); 
          //cout<<"Mass: "<<invmasst<<endl;

          MWLikelihood_RecoD_Inv_Mass_t->Fill(invmasst);  
	  MWLikelihood_RecoD_Inv_Mass_t_->Fill(invmasst_);
	  
          TLorentzVector t_b = j[2] + j[3];

          double invmasst_b = (t_b).M();


          double invmasst__b = (e+m).M(); 
          MWLikelihood_RecoD_Inv_Mass_t_b->Fill(invmasst_b);  
	  MWLikelihood_RecoD_Inv_Mass_t__b->Fill(invmasst__b);

          cout<<"RecoD:"<<invmasst_<<endl;
	  list>>idD;
          cout<<idD<<endl; 
	}
      
      }


  if (bTag_madWeightLikelihood_RecoD)
    {
      ifstream list;
      list.open("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/MadWeight/Results/background/top/recoD/list.txt");
      string idD; 
      list>>idD;
      cout<<idD<<endl;
      bool noMore = false;
      while (!list.eof())
	{
   	  float  error, likelihood,likelihood1, likelihood2, realLikelihood;
          int numberOfPermutation1, numberOfPermutation2;
	  
	  
          likelihood1=0;
          likelihood2=0;
	  

          ifstream mwlikelihood;
          const char * file= ("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/MadWeight/Results/background/top/recoD/"+idD).c_str();
          cout<<file<<endl;
	  mwlikelihood.open(file);
          //string line;
	  //getline(mwlikelihood,line);
          //cout<<line<<endl;
	  string tmp1;
          mwlikelihood>>tmp1;
          cout<<"JP"<<tmp1<<endl;
	  mwlikelihood>>likelihood>>error;
          cout<<"likelihood: "<<likelihood<<endl;
          ifstream LHCO;
          LHCO.open(("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/LHCO/background/top/recoD/"+idD).c_str());
          double tmp_;
	  double jeta[4];
	  double jphi[4];
	  double jpt[4];
	  double jm[4];
	  double eeta,ephi,ept,em;
	  double meta,mphi,mpt,mm;
          double tmp_B;
	  double jetaB[4];
	  double jphiB[4];
	  double jptB[4];
	  double jmB[4];
	  double eetaB,ephiB,eptB,emB;
	  double metaB,mphiB,mptB,mmB;
          LHCO>>tmp_>>tmp_;
          cout<<tmp_<<endl;
          LHCO>>tmp_;
          LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	  ifstream discriminantsD;
          const char * file1= ("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/crab/discriminantsDCommon.txt");
	  discriminantsD.open(file1);
          const char init = tmp1[0];
          double discriminants1[4];
          double discriminantsB[4]; 
          double likelihoodBtag;
          double eventperm;
          int combination = 0;
          int bestCombination = 0;
          while (init=='1')
            {
              combination ++; 
              discriminantsD>>eventperm>>discriminants1[0]>>discriminants1[1]>>discriminants1[2]>>discriminants1[3]; 
	      likelihoodBtag=byHandb(discriminants1[0])*byHandb(discriminants1[1])*byHandc(discriminants1[2])*byHandc(discriminants1[3]);//bb~jj~ 
              likelihood = likelihood * likelihoodBtag;

	      Btag_MWLikelihood_RecoD_Total->Fill(likelihood);
	      if (likelihood>likelihood1)
		{
		  likelihood2 = likelihood1;
		  likelihood1 = likelihood;
                  bestCombination = combination; 
                  for (int k=0 ; k<4 ; k++)
                  {
                    jetaB[k]=jeta[k];
                    jphiB[k]=jphi[k];
                    jptB[k]=jpt[k];
                    jmB[k]=jm[k];
                    discriminantsB[k]=discriminants1[k];
                    cout<<"entro"<<endl;
                  } 
                  eetaB=eeta;
                  ephiB=ephi;
                  eptB=ept;
                  emB=em;
                  metaB=meta;
                  mphiB=mphi;
                  mptB=mpt;
                  mmB=mm;
		}
              LHCO>>tmp_>>tmp_;
              cout<<tmp_<<endl;
              LHCO>>tmp_;
              LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      
              mwlikelihood>>tmp1>>likelihood>>error;
              const char init = tmp1[0];
	    }
	  mwlikelihood.close();
          discriminantsD.close();
	  Btag_MatchingWithPreSelection->Fill(bestCombination);
	  likelihood1 = -log(likelihood1);
	  Btag_MWLikelihood_RecoD_Best->Fill(likelihood1);
	  float _difference1_2 = likelihood1 - likelihood2;
	  Btag_MWLikelihood_RecoD_Difference1_2->Fill(_difference1_2);
	  LHCO.close();
	  vector<TLorentzVector> j(4);
	  TLorentzVector e;
	  TLorentzVector m;
          j[0].SetPtEtaPhiM(jptB[0],jetaB[0],jphiB[0],jmB[0]);
	  j[1].SetPtEtaPhiM(jptB[1],jetaB[1],jphiB[1],jmB[1]);
	  j[2].SetPtEtaPhiM(jptB[2],jetaB[2],jphiB[2],jmB[2]);
	  j[3].SetPtEtaPhiM(jptB[3],jetaB[3],jphiB[3],jmB[3]);
	  e.SetPtEtaPhiM(eptB,eetaB,ephiB,emB);
	  m.SetPtEtaPhiM(mptB,metaB,mphiB,mmB);
	  double invmasst, invmasst_;
	  TLorentzVector t = j[2] + j[3] + j[0];
          invmasst = t.M();
	  invmasst_=(e+m+j[1]).M(); 
          //cout<<"Mass: "<<invmasst<<endl;
          Btag_MWLikelihood_RecoD_Inv_Mass_t->Fill(invmasst);  
	  Btag_MWLikelihood_RecoD_Inv_Mass_t_->Fill(invmasst_);


	  TLorentzVector t_b = j[2] + j[3];

          double invmasst_b = (t_b).M();



          double invmasst__b = (e+m).M();
          Btag_MWLikelihood_RecoD_Inv_Mass_t_b->Fill(invmasst_b);
          Btag_MWLikelihood_RecoD_Inv_Mass_t__b->Fill(invmasst__b);




          cout<<"RecoD:"<<invmasst_<<endl;
	  list>>idD;
          cout<<idD<<endl; 
	}
      
      }







  if (MC)
    {
      double Discriminants[4];
      string id;
  
      ifstream discriminants;
      // discriminants.open("BTAGPerformance/discriminantsCommon.txt");
      // discriminants>>id>>Discriminants[0]>>Discriminants[1]>>Discriminants[2]>>Discriminants[3];
      discriminants.open("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/MadWeight/common.txt");
      discriminants>>id;
      while (!discriminants.eof())
	{
	  cout<<"ID: "<<id<<endl;

	  const int numberOfJets1 = 4; //later it could be 5 to include ISR    
	  const int numberOfPermutations = 24; //later it can be improved depending on the numberOfJets1!                                                                                       
      
	  int permutationMap[numberOfPermutations][numberOfJets1]= {{1,2,3,4},
								    {2,1,3,4},
								    {3,1,2,4},
								    {1,3,2,4},
								    {2,3,1,4},
								    {3,2,1,4},
								    {4,2,1,3},
								    {2,4,1,3},
								    {1,4,2,3},
								    {4,1,2,3},
								    {2,1,4,3},
								    {1,2,4,3},
								    {1,3,4,2},
								    {3,1,4,2},
								    {4,1,3,2},
								    {1,4,3,2},
								    {3,4,1,2},
								    {4,3,1,2},
								    {4,3,2,1},
								    {3,4,2,1},
								    {2,4,3,1},
								    {4,2,3,1},
								    {3,2,4,1},
								    {2,3,4,1}}; 
      
      
      
	  if (bTaggingLikelihood)
	    {
	      float likelihood,likelihood1, likelihood2, realLikelihood;
	      int numberOfPermutation1, numberOfPermutation2;
	  
	  
	      likelihood1=0;
	      likelihood2=0;
	  
	      numberOfPermutation1=0;
	      numberOfPermutation2=0;
	      int numberOfPermutation = 0;
	      for (int numberOfPermutation=0;numberOfPermutation<24;numberOfPermutation++)
		{
		  likelihood=byHandb(Discriminants[permutationMap[numberOfPermutation][0]-1])*byHandb(Discriminants[permutationMap[numberOfPermutation][1]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][2]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][3]-1]);//bb~jj~  
		  bTaggingLikelihood_Total->Fill(likelihood);

		  /////////////
		  //cout<<permutationMap[numberOfPermutation][0]-1<<" "<<permutationMap[numberOfPermutation][1]-1<<" "<<permutationMap[numberOfPermutation][2]-1<<" "<<permutationMap[numberOfPermutation][3]-1<<endl;
	      
	      
		  //cout<<Discriminants[permutationMap[numberOfPermutation][0]-1]<<" "<<Discriminants[permutationMap[numberOfPermutation][1]-1]<<" "<<Discriminants[permutationMap[numberOfPermutation][2]-1]<<" "<<Discriminants[permutationMap[numberOfPermutation][3]-1]<<endl;
																											 
		  //	      cout<<byHand(Discriminants[permutationMap[numberOfPermutation][0]-1])<<" "<<byHand(Discriminants[permutationMap[numberOfPermutation][1]-1])<<" "<<byHand(Discriminants[permutationMap[numberOfPermutation][2]-1]))<<" "<<byHand(Discriminants[permutationMap[numberOfPermutation][3]-1]))<<endl;
											       
											       


		  ///////////


		  //cout<<"Btagging likelihood: "<<likelihood<<endl;	
		  if (likelihood>likelihood1)
		    {
		      likelihood2 = likelihood1;
		      likelihood1 = likelihood;
		      numberOfPermutation2=numberOfPermutation1;
		      numberOfPermutation1=numberOfPermutation;
		    }
		}
	      realLikelihood=byHandb(Discriminants[0])*byHandb(Discriminants[1])*byHandc(Discriminants[2])*byHandc(Discriminants[3]);
	      bTaggingLikelihood_Real->Fill(realLikelihood);
	      bTaggingLikelihood_Best->Fill(likelihood1);
	      cout<<"JP:  btag:  real="<<realLikelihood<<" best="<<likelihood1<<endl;
	      float _differenceReal_1 = realLikelihood - likelihood1;
	      bTaggingLikelihood_DifferenceReal_1->Fill(_differenceReal_1);
	      float _difference1_2 = likelihood1 - likelihood2;
	      bTaggingLikelihood_Difference1_2->Fill(_difference1_2);
	      int _numberOfJetsIdentified = 0;
	      if ( permutationMap[numberOfPermutation1][0] == 1 || permutationMap[numberOfPermutation1][0] == 2  )
		{
		  _numberOfJetsIdentified++;
		}  
	      if ( permutationMap[numberOfPermutation1][1] == 1 || permutationMap[numberOfPermutation1][1] == 2  )
		{
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][2] == 3 || permutationMap[numberOfPermutation1][2] == 4  )
		{
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][3] == 3 || permutationMap[numberOfPermutation1][3] == 4  )
		{
		  _numberOfJetsIdentified++;
		}
	      bTaggingLikelihood_NumberOfJetsIdentified->Fill(_numberOfJetsIdentified);
	  
	    }
      
      
      
      
      
      
      
	  if (madWeightLikelihood_Gen)
	    {
	      float likelihood,likelihood1, likelihood2, realLikelihood, error;
	      int numberOfPermutation1, numberOfPermutation2;


	      likelihood1=0;
	      likelihood2=0;
	      numberOfPermutation1=0;
	      numberOfPermutation2=0;
	      int numberOfPermutation = 0;
	      ifstream mwlikelihood;
	      mwlikelihood.open(("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/MadWeight/Results/background/top/gen/"+id).c_str());
	      string tmp;
	      mwlikelihood>>tmp;
              cout<<tmp<<endl;
              mwlikelihood>>likelihood>>error;
	      //cout<<"Likelihood:"<<likelihood<<endl;

	      for (int numberOfPermutation=0;numberOfPermutation<24;numberOfPermutation++)
		{
		  /*if (tmp[0]!=1)
		    {
		    cout<<"el archivo de pesos con id: "<<id<<" esta corrupto."<<endl;
		    exit(0);
		    }*/
		  if (numberOfPermutation==0)
		    {
		      realLikelihood=likelihood;
		    }
		  //cout<<"JP2: "<<tmp<<endl;
		  //cout<<"JP1: "<<likelihood<<endl;

		  MWLikelihood_Gen_Total->Fill(likelihood);
		  if (likelihood>likelihood1)
		    {
		      likelihood2 = likelihood1;
		      likelihood1 = likelihood;
		      numberOfPermutation2=numberOfPermutation1;
		      numberOfPermutation1=numberOfPermutation;
		    }
		  mwlikelihood>>tmp>>likelihood>>error;
		}
	      mwlikelihood.close();
	      cout<<"JP:  gen:  real="<<realLikelihood<<" best="<<likelihood1<<endl;
	      //cout<<"Real JP: "<<realLikelihood<<endl;
	      MWLikelihood_Gen_Real->Fill(realLikelihood);
	      likelihood1 = -log(likelihood1);
	      MWLikelihood_Gen_Best->Fill(likelihood1);
	      float _differenceReal_1 = realLikelihood - likelihood1;
	      MWLikelihood_Gen_DifferenceReal_1->Fill(_differenceReal_1);
	      float _difference1_2 = likelihood1 - likelihood2;
	      MWLikelihood_Gen_Difference1_2->Fill(_difference1_2);
	      int _numberOfJetsIdentified = 0;
	      if ( permutationMap[numberOfPermutation1][0] == 1 )
		{
		  MWLikelihood_Gen_JetsIdentified->Fill(0);
		  _numberOfJetsIdentified++;
		}  
	      if ( permutationMap[numberOfPermutation1][1] == 2 )
		{
		  MWLikelihood_Gen_JetsIdentified->Fill(1); 	      
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][2] == 3 || permutationMap[numberOfPermutation1][2] == 4  )
		{
		  MWLikelihood_Gen_JetsIdentified->Fill(2);
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][3] == 4 || permutationMap[numberOfPermutation1][3] == 3  )
		{
		  MWLikelihood_Gen_JetsIdentified->Fill(3);
		  _numberOfJetsIdentified++;
		}
	      MWLikelihood_Gen_NumberOfJetsIdentified->Fill(_numberOfJetsIdentified);
	      ///New in 2013
	      ifstream LHCO;
	      LHCO.open(("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/LHCO/background/top/gen/"+id).c_str());
	      double tmp_;
	      double jeta[4];
	      double jphi[4];
	      double jpt[4];
	      double jm[4];
	      double eeta,ephi,ept,em;
	      double meta,mphi,mpt,mm;
	      LHCO>>tmp_;
              cout<<tmp_<<endl;
	      LHCO>>tmp_;
	      cout<<"Event: "<<tmp_<<endl;
	      LHCO>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO.close();
	      cout<<"mpt="<<mpt<<" ephi="<<ephi<<" beta="<<jeta[0]<<" b_eta="<<jeta[1]<<"q eta="<<jeta[2]<<" q_eta="<<jeta[3]<<" eeta="<<eeta<<endl;
	      vector<TLorentzVector> j(4);
	      TLorentzVector e;
	      TLorentzVector m;
	      j[0].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][0]-1],jeta[permutationMap[numberOfPermutation1][0]-1],jphi[permutationMap[numberOfPermutation1][0]-1],jm[permutationMap[numberOfPermutation1][0]-1]);
	      j[1].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][1]-1],jeta[permutationMap[numberOfPermutation1][1]-1],jphi[permutationMap[numberOfPermutation1][1]-1],jm[permutationMap[numberOfPermutation1][1]-1]);
	      j[2].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][2]-1],jeta[permutationMap[numberOfPermutation1][2]-1],jphi[permutationMap[numberOfPermutation1][2]-1],jm[permutationMap[numberOfPermutation1][2]-1]);
	      j[3].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][3]-1],jeta[permutationMap[numberOfPermutation1][3]-1],jphi[permutationMap[numberOfPermutation1][3]-1],jm[permutationMap[numberOfPermutation1][3]-1]);
	      e.SetPtEtaPhiM(ept,eeta,ephi,em);
	      m.SetPtEtaPhiM(mpt,meta,mphi,mm);
	      double invmasst, invmasst_;
	      TLorentzVector t = j[2] + j[3] + j[0];
	      invmasst = t.M();
	      invmasst_=(e+m+j[1]).M(); 
	      //cout<<"MassJP: "<<invmasst<<endl;
	      MWLikelihood_Gen_Inv_Mass_t->Fill(invmasst);  
	      MWLikelihood_Gen_Inv_Mass_t_->Fill(invmasst_);
	      

             
	      TLorentzVector t_b = j[2] + j[3];

	      double invmasst_b = (t_b).M();



	      double invmasst__b = (e+m).M();

              MWLikelihood_Gen_Inv_Mass_t_b->Fill(invmasst_b);  
	      MWLikelihood_Gen_Inv_Mass_t__b->Fill(invmasst__b);


              if (_numberOfJetsIdentified==4)
		{
		  cout<<"MassJPGen: "<<invmasst<<endl;  
		}
	    }
      
      
      
      
      
      
          
      
      
	  if (madWeightLikelihood_Reco)
	    {
	     	     
	      float likelihood,likelihood1, likelihood2, realLikelihood, error;
	      int numberOfPermutation1, numberOfPermutation2;
	  
	  
	      likelihood1=0;
	      likelihood2=0;
	  
	      numberOfPermutation1=0;
	      numberOfPermutation2=0;
	      int numberOfPermutation = 0;
              cout<<"el id es:"<<id<<endl;
	      ifstream mwlikelihood;
	      mwlikelihood.open(("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/MadWeight/Results/background/top/reco/"+id).c_str());
	      string tmp;
	      mwlikelihood>>tmp;
              cout<<tmp<<endl;
              mwlikelihood>>likelihood>>error;

	      for (int numberOfPermutation=0;numberOfPermutation<24;numberOfPermutation++)
		{
		  /*if (tmp!=1) 
		    {	
		    cout<<"el archivo de pesos con id: "<<id<<" esta corrupto."<<endl;
		    exit(0);
		    }*/ 
		  if (numberOfPermutation==0)
		    {
		      realLikelihood=likelihood;
		    }
		  //cout<<likelihood<<endl;
	      
		  MWLikelihood_Reco_Total->Fill(likelihood);
		  if (likelihood>likelihood1)
		    {
		      likelihood2 = likelihood1;
		      likelihood1 = likelihood;
		      numberOfPermutation2=numberOfPermutation1;
		      numberOfPermutation1=numberOfPermutation;
		    }
		  mwlikelihood>>tmp>>likelihood>>error;
		}
	      mwlikelihood.close();
	      Reco_Matching->Fill(1);
              if ((numberOfPermutation1==0) || (numberOfPermutation1==11))
		{
		  Reco_Matching->Fill(3);
		}
              MWLikelihood_Reco_Real->Fill(realLikelihood);
	      likelihood1 = -log(likelihood1);
	      MWLikelihood_Reco_Best->Fill(likelihood1);
	      cout<<"JP:  reco:  real="<<realLikelihood<<" best="<<likelihood1<<endl;
	      float _differenceReal_1 = realLikelihood - likelihood1;
	      MWLikelihood_Reco_DifferenceReal_1->Fill(_differenceReal_1);
	      float _difference1_2 = likelihood1 - likelihood2;
	      MWLikelihood_Reco_Difference1_2->Fill(_difference1_2);
	      int _numberOfJetsIdentified = 0;
	      if ( permutationMap[numberOfPermutation1][0] == 1 )
		{
		  MWLikelihood_Reco_JetsIdentified->Fill(0);
		  _numberOfJetsIdentified++;
		}  
	      if ( permutationMap[numberOfPermutation1][1] == 2 )
		{
		  MWLikelihood_Reco_JetsIdentified->Fill(1); 	      
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][2] == 3 || permutationMap[numberOfPermutation1][2] == 4 )
		{
		  MWLikelihood_Reco_JetsIdentified->Fill(2);
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][3] == 4  || permutationMap[numberOfPermutation1][3] == 3   )
		{
		  MWLikelihood_Reco_JetsIdentified->Fill(3);
		  _numberOfJetsIdentified++;
		}
	      MWLikelihood_Reco_NumberOfJetsIdentified->Fill(_numberOfJetsIdentified);
	      ///New in 2013
	      ifstream LHCO;
	      LHCO.open(("/afs/cern.ch/work/j/jgomezca/thesis/CMSSW_5_3_9/LHCO/background/top/reco/"+id).c_str());
	      double tmp_;
	      double jeta[4];
	      double jphi[4];
	      double jpt[4];
	      double jm[4];
	      double eeta,ephi,ept,em;
	      double meta,mphi,mpt,mm;
	      LHCO>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO.close();
	      cout<<"mpt="<<mpt<<" ephi="<<ephi<<" beta="<<jeta[0]<<" b_eta="<<jeta[1]<<" qeta="<<jeta[2]<<" q_eta="<<jeta[3]<<" eeta="<<eeta<<endl;
	      vector<TLorentzVector> j(4);
	      TLorentzVector e;
	      TLorentzVector m;
	      j[0].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][0]-1],jeta[permutationMap[numberOfPermutation1][0]-1],jphi[permutationMap[numberOfPermutation1][0]-1],jm[permutationMap[numberOfPermutation1][0]-1]);
	      j[1].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][1]-1],jeta[permutationMap[numberOfPermutation1][1]-1],jphi[permutationMap[numberOfPermutation1][1]-1],jm[permutationMap[numberOfPermutation1][1]-1]);
	      j[2].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][2]-1],jeta[permutationMap[numberOfPermutation1][2]-1],jphi[permutationMap[numberOfPermutation1][2]-1],jm[permutationMap[numberOfPermutation1][2]-1]);
	      j[3].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][3]-1],jeta[permutationMap[numberOfPermutation1][3]-1],jphi[permutationMap[numberOfPermutation1][3]-1],jm[permutationMap[numberOfPermutation1][3]-1]);
	      e.SetPtEtaPhiM(ept,eeta,ephi,em);
	      m.SetPtEtaPhiM(mpt,meta,mphi,mm);
	      double invmasst, invmasst_;
	      TLorentzVector t = j[2] + j[3] + j[0];
	      invmasst = t.M();
	      invmasst_=(e+m+j[1]).M(); 
	      //cout<<"Mass: "<<invmasst<<endl;
	      MWLikelihood_Reco_Inv_Mass_t->Fill(invmasst);  
	      MWLikelihood_Reco_Inv_Mass_t_->Fill(invmasst_);
	      
	      TLorentzVector t_b = j[2] + j[3];

	      double invmasst_b = (t_b).M();



	      double invmasst__b = (e+m).M();
	      MWLikelihood_Reco_Inv_Mass_t_b->Fill(invmasst_b);
	      MWLikelihood_Reco_Inv_Mass_t__b->Fill(invmasst__b);
 


              if (_numberOfJetsIdentified==4)
		{
		  cout<<"MassJPReco: "<<invmasst<<endl;
		}

	    }
      
      






      
	  if (bTagging_MWLikelihoodRecoWithoutTF)
	    {
	      float blikelihood,likelihood,likelihood1, likelihood2, realLikelihood;
	      int numberOfPermutation1, numberOfPermutation2;


	      likelihood1=0;
	      likelihood2=0;

	      numberOfPermutation1=0;
	      numberOfPermutation2=0;
	      int numberOfPermutation = 0;
	      ifstream mwlikelihood;
	      mwlikelihood.open(("BTAGPerformance/ResultsWOReco/"+id+"0/Analysis/Events/fermi/fermi_weights.out").c_str());
	      string tmp;
	      mwlikelihood>>tmp>>likelihood>>error;

	      for (int numberOfPermutation=0;numberOfPermutation<24;numberOfPermutation++)
		{
		  /*if (tmp!=1)
		    {
		    cout<<"el archivo de pesos con id: "<<id<<" esta corrupto."<<endl;
		    exit(0);
		    }*/
              
		  //cout<<likelihood<<endl;
		  blikelihood=byHandb(Discriminants[permutationMap[numberOfPermutation][0]-1])*byHandb(Discriminants[permutationMap[numberOfPermutation][1]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][2]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][3]-1]);//bb~jj~  
		  likelihood=likelihood*blikelihood;	      
		  if (numberOfPermutation==0)
		    {
		      realLikelihood=likelihood;
		    }	
	      
	      
		  BTaggingMWLikelihood_RecoWithoutTF_Total->Fill(likelihood);
		  if (likelihood>likelihood1)
		    {
		      likelihood2 = likelihood1;
		      likelihood1 = likelihood;
		      numberOfPermutation2=numberOfPermutation1;
		      numberOfPermutation1=numberOfPermutation;
		    }
		  mwlikelihood>>tmp>>likelihood>>error;
		}
	      mwlikelihood.close();
	      BTaggingMWLikelihood_RecoWithoutTF_Real->Fill(realLikelihood);
	      BTaggingMWLikelihood_RecoWithoutTF_Best->Fill(likelihood1);
	      cout<<"JP:  btagrecoWO:  real="<<realLikelihood<<" best="<<likelihood1<<" solo btag="<<blikelihood<<endl;
	      float _differenceReal_1 = realLikelihood - likelihood1;
	      BTaggingMWLikelihood_RecoWithoutTF_DifferenceReal_1->Fill(_differenceReal_1);
	      float _difference1_2 = likelihood1 - likelihood2;
	      BTaggingMWLikelihood_RecoWithoutTF_Difference1_2->Fill(_difference1_2);
	      int _numberOfJetsIdentified = 0;
	      if ( permutationMap[numberOfPermutation1][0] == 1 )
		{
		  BTaggingMWLikelihood_RecoWithoutTF_JetsIdentified->Fill(0);
		  _numberOfJetsIdentified++;
		}  
	      if ( permutationMap[numberOfPermutation1][1] == 2  )
		{
		  BTaggingMWLikelihood_RecoWithoutTF_JetsIdentified->Fill(1); 	      
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][2] == 3 || permutationMap[numberOfPermutation1][2] == 4  )
		{
		  BTaggingMWLikelihood_RecoWithoutTF_JetsIdentified->Fill(2);
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][3] == 4 || permutationMap[numberOfPermutation1][3] == 3  )
		{
		  BTaggingMWLikelihood_RecoWithoutTF_JetsIdentified->Fill(3);
		  _numberOfJetsIdentified++;
		}
	      BTaggingMWLikelihood_RecoWithoutTF_NumberOfJetsIdentified->Fill(_numberOfJetsIdentified);
	      ///New in 2013
	      ifstream LHCO;
	      LHCO.open(("BTAGPerformance/ResultsWOReco/"+id+"0/Analysis/Events/input.lhco").c_str());
	      double tmp_;
	      double jeta[4];
	      double jphi[4];
	      double jpt[4];
	      double jm[4];
	      double eeta,ephi,ept,em;
	      double meta,mphi,mpt,mm;
	      LHCO>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO.close();
	      cout<<"mpt="<<mpt<<" ephi="<<ephi<<" beta="<<jeta[0]<<" b_eta="<<jeta[1]<<" qeta="<<jeta[2]<<" q_eta="<<jeta[3]<<" eeta="<<eeta<<endl;
	      vector<TLorentzVector> j(4);
	      TLorentzVector e;
	      TLorentzVector m;
	      j[0].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][0]-1],jeta[permutationMap[numberOfPermutation1][0]-1],jphi[permutationMap[numberOfPermutation1][0]-1],jm[permutationMap[numberOfPermutation1][0]-1]);
	      j[1].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][1]-1],jeta[permutationMap[numberOfPermutation1][1]-1],jphi[permutationMap[numberOfPermutation1][1]-1],jm[permutationMap[numberOfPermutation1][1]-1]);
	      j[2].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][2]-1],jeta[permutationMap[numberOfPermutation1][2]-1],jphi[permutationMap[numberOfPermutation1][2]-1],jm[permutationMap[numberOfPermutation1][2]-1]);
	      j[3].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][3]-1],jeta[permutationMap[numberOfPermutation1][3]-1],jphi[permutationMap[numberOfPermutation1][3]-1],jm[permutationMap[numberOfPermutation1][3]-1]);
	      e.SetPtEtaPhiM(ept,eeta,ephi,em);
	      m.SetPtEtaPhiM(mpt,meta,mphi,mm);
	      double invmasst, invmasst_;
	      TLorentzVector t = j[2] + j[3] + j[0];
	      invmasst = t.M();
	      invmasst_=(e+m+j[1]).M(); 
	      //cout<<"Mass: "<<invmasst<<endl;
	      BTaggingMWLikelihood_RecoWithoutTF_Inv_Mass_t->Fill(invmasst);  
	      BTaggingMWLikelihood_RecoWithoutTF_Inv_Mass_t_->Fill(invmasst_);
	      if (_numberOfJetsIdentified==4)
		{
		  cout<<"MassJPBtagRecoWO: "<<invmasst<<endl;
		}

	    }
      
      
 
	  if (bTagging_MWLikelihoodReco)
	    {
	      float likelihood,blikelihood,likelihood1, likelihood2, realLikelihood;
	      int numberOfPermutation1, numberOfPermutation2;


	      likelihood1=0;
	      likelihood2=0;

	      numberOfPermutation1=0;
	      numberOfPermutation2=0;
	      int numberOfPermutation = 0;
	      ifstream mwlikelihood;
	      mwlikelihood.open(("BTAGPerformance/ResultsReco/"+id+"0/Analysis/Events/fermi/fermi_weights.out").c_str());
	      string tmp;
	      mwlikelihood>>tmp>>likelihood>>error;

	      for (int numberOfPermutation=0;numberOfPermutation<24;numberOfPermutation++)
		{
		  /*if (tmp!=1)
		    {
		    cout<<"el archivo de pesos con id: "<<id<<" esta corrupto."<<endl;
		    exit(0);
		    }*/
              
		  //cout<<likelihood<<endl;
		  blikelihood=byHandb(Discriminants[permutationMap[numberOfPermutation][0]-1])*byHandb(Discriminants[permutationMap[numberOfPermutation][1]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][2]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][3]-1]);//bb~jj~  
		  likelihood=likelihood*blikelihood;
		  if (numberOfPermutation==0)
		    {
		      realLikelihood=likelihood;
		    }	      
		  BTaggingMWLikelihood_Reco_Total->Fill(likelihood);
		  if (likelihood>likelihood1)
		    {
		      likelihood2 = likelihood1;
		      likelihood1 = likelihood;
		      numberOfPermutation2=numberOfPermutation1;
		      numberOfPermutation1=numberOfPermutation;
		    }
		  mwlikelihood>>tmp>>likelihood>>error;
		}
	      mwlikelihood.close();
	      BTaggingMWLikelihood_Reco_Real->Fill(realLikelihood);
	      BTaggingMWLikelihood_Reco_Best->Fill(likelihood1);
	      cout<<"JP:  btagreco:  real="<<realLikelihood<<" best="<<likelihood1<<" solo btag="<<blikelihood<<endl;
	      float _differenceReal_1 = realLikelihood - likelihood1;
	      BTaggingMWLikelihood_Reco_DifferenceReal_1->Fill(_differenceReal_1);
	      float _difference1_2 = likelihood1 - likelihood2;
	      BTaggingMWLikelihood_Reco_Difference1_2->Fill(_difference1_2);
	      int _numberOfJetsIdentified = 0;
	      if ( permutationMap[numberOfPermutation1][0] == 1 )
		{
		  BTaggingMWLikelihood_Reco_JetsIdentified->Fill(0);
		  _numberOfJetsIdentified++;
		}  
	      if ( permutationMap[numberOfPermutation1][1] == 2  )
		{
		  BTaggingMWLikelihood_Reco_JetsIdentified->Fill(1); 	      
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][2] == 3 || permutationMap[numberOfPermutation1][2] == 4  )
		{
		  BTaggingMWLikelihood_Reco_JetsIdentified->Fill(2);
		  _numberOfJetsIdentified++;
		}
	      if ( permutationMap[numberOfPermutation1][3] == 4 || permutationMap[numberOfPermutation1][3] == 3 )
		{
		  BTaggingMWLikelihood_Reco_JetsIdentified->Fill(3);
		  _numberOfJetsIdentified++;
		}
	      BTaggingMWLikelihood_Reco_NumberOfJetsIdentified->Fill(_numberOfJetsIdentified);
	      ///New in 2013
	      ifstream LHCO;
	      LHCO.open(("BTAGPerformance/ResultsReco/"+id+"0/Analysis/Events/input.lhco").c_str());
	      double tmp_;
	      double jeta[4];
	      double jphi[4];
	      double jpt[4];
	      double jm[4];
	      double eeta,ephi,ept,em;
	      double meta,mphi,mpt,mm;
	      LHCO>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	      LHCO.close();
	      cout<<"mpt="<<mpt<<" ephi="<<ephi<<" beta="<<jeta[0]<<" b_eta="<<jeta[1]<<" qeta="<<jeta[2]<<" q_eta="<<jeta[3]<<" eeta="<<eeta<<endl;
	      vector<TLorentzVector> j(4);
	      TLorentzVector e;
	      TLorentzVector m;
	      j[0].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][0]-1],jeta[permutationMap[numberOfPermutation1][0]-1],jphi[permutationMap[numberOfPermutation1][0]-1],jm[permutationMap[numberOfPermutation1][0]-1]);
	      j[1].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][1]-1],jeta[permutationMap[numberOfPermutation1][1]-1],jphi[permutationMap[numberOfPermutation1][1]-1],jm[permutationMap[numberOfPermutation1][1]-1]);
	      j[2].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][2]-1],jeta[permutationMap[numberOfPermutation1][2]-1],jphi[permutationMap[numberOfPermutation1][2]-1],jm[permutationMap[numberOfPermutation1][2]-1]);
	      j[3].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][3]-1],jeta[permutationMap[numberOfPermutation1][3]-1],jphi[permutationMap[numberOfPermutation1][3]-1],jm[permutationMap[numberOfPermutation1][3]-1]);
	      e.SetPtEtaPhiM(ept,eeta,ephi,em);
	      m.SetPtEtaPhiM(mpt,meta,mphi,mm);
	      double invmasst, invmasst_;
	      TLorentzVector t = j[2] + j[3] + j[0];
	      invmasst = t.M();
	      invmasst_=(e+m+j[1]).M(); 
	      //cout<<"Mass: "<<invmasst<<endl;
	      BTaggingMWLikelihood_Reco_Inv_Mass_t->Fill(invmasst);  
	      BTaggingMWLikelihood_Reco_Inv_Mass_t_->Fill(invmasst_);
	      if (_numberOfJetsIdentified==4)
		{
		  cout<<"MassJPBtagReco: "<<invmasst<<endl;
		}

	    }
      
      
 

	  /*
	    if (bTagging_MWLikelihoodGen)
	    {
	    float likelihood,blikelihood,likelihood1, likelihood2, realLikelihood;
	    int numberOfPermutation1, numberOfPermutation2;


	    likelihood1=0;
	    likelihood2=0;

	    numberOfPermutation1=0;
	    numberOfPermutation2=0;
	    int numberOfPermutation = 0;
	    ifstream mwlikelihood;
	    mwlikelihood.open(("BTAGPerformance/ResultsGen/"+id+"0/Analysis/Events/fermi/fermi_weights.out").c_str());
	    string tmp;
	    mwlikelihood>>tmp>>likelihood>>error;

	    for (int numberOfPermutation=0;numberOfPermutation<24;numberOfPermutation++)
            {
            
              
	    cout<<likelihood<<endl;
	    blikelihood=byHandb(Discriminants[permutationMap[numberOfPermutation][0]-1])*byHandb(Discriminants[permutationMap[numberOfPermutation][1]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][2]-1])*byHandc(Discriminants[permutationMap[numberOfPermutation][3]-1]);//bb~jj~  
	    likelihood=likelihood*blikelihood;
	    if (numberOfPermutation==0)
	    {
	    realLikelihood=likelihood;
	    }
	    BTaggingMWLikelihood_Gen_Total->Fill(likelihood);
	    if (likelihood>likelihood1)
	    {
	    likelihood2 = likelihood1;
	    likelihood1 = likelihood;
	    numberOfPermutation2=numberOfPermutation1;
	    numberOfPermutation1=numberOfPermutation;
	    }
	    mwlikelihood>>tmp>>likelihood>>error;
	    }
	    mwlikelihood.close();
	    BTaggingMWLikelihood_Gen_Real->Fill(realLikelihood);
	    BTaggingMWLikelihood_Gen_Best->Fill(likelihood1);
	    float _differenceReal_1 = realLikelihood - likelihood1;
	    BTaggingMWLikelihood_Gen_DifferenceReal_1->Fill(_differenceReal_1);
	    float _difference1_2 = likelihood1 - likelihood2;
	    BTaggingMWLikelihood_Gen_Difference1_2->Fill(_difference1_2);
	    int _numberOfJetsIdentified = 0;
	    if ( permutationMap[numberOfPermutation1][0] == 1 )
	    {
	    BTaggingMWLikelihood_Gen_JetsIdentified->Fill(0);
	    _numberOfJetsIdentified++;
	    }  
	    if ( permutationMap[numberOfPermutation1][1] == 2  || permutationMap[numberOfPermutation1][1] == 3  )
	    {
	    BTaggingMWLikelihood_Gen_JetsIdentified->Fill(1); 	      
	    _numberOfJetsIdentified++;
	    }
	    if ( permutationMap[numberOfPermutation1][2] == 2  || permutationMap[numberOfPermutation1][2] == 3  )
	    {
	    BTaggingMWLikelihood_Gen_JetsIdentified->Fill(2);
	    _numberOfJetsIdentified++;
	    }
	    if ( permutationMap[numberOfPermutation1][3] == 4  )
	    {
	    BTaggingMWLikelihood_Gen_JetsIdentified->Fill(3);
	    _numberOfJetsIdentified++;
	    }
	    BTaggingMWLikelihood_Gen_NumberOfJetsIdentified->Fill(_numberOfJetsIdentified);
	    ///New in 2013
	    ifstream LHCO;
	    LHCO.open(("BTAGPerformance/ResultsGen/"+id+"0/Analysis/Events/input.lhco").c_str());
	    double tmp_;
	    double jeta[4];
	    double jphi[4];
	    double jpt[4];
	    double jm[4];
	    double eeta,ephi,ept,em;
	    double meta,mphi,mpt,mm;
	    LHCO>>tmp_>>tmp_>>tmp_;
	    LHCO>>tmp_>>tmp_>>jeta[0]>>jphi[0]>>jpt[0]>>jm[0]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	    LHCO>>tmp_>>tmp_>>jeta[1]>>jphi[1]>>jpt[1]>>jm[1]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	    LHCO>>tmp_>>tmp_>>jeta[2]>>jphi[2]>>jpt[2]>>jm[2]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	    LHCO>>tmp_>>tmp_>>jeta[3]>>jphi[3]>>jpt[3]>>jm[3]>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	    LHCO>>tmp_>>tmp_>>eeta>>ephi>>ept>>em>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	    LHCO>>tmp_>>tmp_>>meta>>mphi>>mpt>>mm>>tmp_>>tmp_>>tmp_>>tmp_>>tmp_;
	    LHCO.close();
	    cout<<"mpt="<<mpt<<"ephi="<<ephi<<"beta="<<jeta[0]<<"b_eta="<<jeta[1]<<"qeta="<<jeta[2]<<"q_eta="<<jeta[3]<<"eeta="<<eeta<<endl;
	    vector<TLorentzVector> j(4);
	    TLorentzVector e;
	    TLorentzVector m;
	    j[0].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][0]-1],jeta[permutationMap[numberOfPermutation1][0]-1],jphi[permutationMap[numberOfPermutation1][0]-1],jm[permutationMap[numberOfPermutation1][0]-1]);
	    j[1].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][1]-1],jeta[permutationMap[numberOfPermutation1][1]-1],jphi[permutationMap[numberOfPermutation1][1]-1],jm[permutationMap[numberOfPermutation1][1]-1]);
	    j[2].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][2]-1],jeta[permutationMap[numberOfPermutation1][2]-1],jphi[permutationMap[numberOfPermutation1][2]-1],jm[permutationMap[numberOfPermutation1][2]-1]);
	    j[3].SetPtEtaPhiM(jpt[permutationMap[numberOfPermutation1][3]-1],jeta[permutationMap[numberOfPermutation1][3]-1],jphi[permutationMap[numberOfPermutation1][3]-1],jm[permutationMap[numberOfPermutation1][3]-1]);
	    e.SetPtEtaPhiM(ept,eeta,ephi,em);
	    m.SetPtEtaPhiM(mpt,meta,mphi,mm);
	    double invmasst, invmasst_;
	    TLorentzVector t = j[2] + j[3] + j[0];
	    invmasst = t.M();
	    invmasst_=(e+m+j[1]).M(); 
	    cout<<"Mass: "<<invmasst<<endl;
	    BTaggingMWLikelihood_Gen_Inv_Mass_t->Fill(invmasst);  
	    BTaggingMWLikelihood_Gen_Inv_Mass_t_->Fill(invmasst_);
	  
	    }
      
	  */




     
      
      
	  //discriminants>>id>>Discriminants[0]>>Discriminants[1]>>Discriminants[2]>>Discriminants[3];
	  discriminants>>id;
	}
    }
      histos.Write();
      histos.Close();
      exit(0);
    }
