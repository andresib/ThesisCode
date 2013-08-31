// -*- C++ -*-
//
// Package:    BTAGPerformance
// Class:      BTAGPerformance
// 
/**\class BTAGPerformance BTAGPerformance.cc bTagAnalyzer/BTAGPerformance/src/BTAGPerformance.cc
   
Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Juan Pablo Gomez Cardona,42 R-023,+41227662349,
//         Created:  Sun Sep  23 17:11:03 CET 2012
// $Id$
//
//


// system include files
#include <memory>


// user include files

//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/Common/interface/TriggerNames.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "JMTucker/CMSSWIntro/interface/GenUtilities.h"
//#include "JMTucker/CMSSWIntro/interface/TriggerHelper.h"
//#include "JMTucker/CMSSWIntro/interface/Utilities.h"







#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"


#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <DataFormats/JetReco/interface/GenJetCollection.h>



#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>


#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//#include "DataFormats/Math/interface/PtEtaPhiM4D.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>



#include <stdio.h>
#include <stdlib.h>
//
// class declaration
//

class FullAndWorking : public edm::EDAnalyzer {
public:
  explicit FullAndWorking(const edm::ParameterSet&);
  ~FullAndWorking();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  double selector,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14;
  std::string bdisc_name;
  bool  top, neutralino, RFijo, LHCOWithMC, LHCOWithData, MC;
  bool lastCombination, lastCombination2, lastCombination3;
  std::ofstream   filet;
  std::ofstream   filet_;
  std::ofstream   file5;  	
  std::ofstream   file6;
  std::ofstream   file7;  	
  std::ofstream   file8;
  std::ofstream   file9;
  /*0404
  //edm::InputTag  bTagAlgo;
  //edm::InputTag  matchingAlgo;
  std::ofstream   file;
  std::ofstream   file1;
  std::ofstream   file10;
  std::ofstream   prueba;
  std::ofstream   file11;
  int numberOfEventsJP, semiNumberOfBJetsJP,numberOf1JetsJP,numberOf2JetsJP;
  0404*/
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  bool nextComb (int j, std::vector<int>& comb, std::vector<int>& oth, int n, int r); 
  double fit_cl(double j); 
  double fit_b(double j);
  struct JetRefCompare :
    public std::binary_function<edm::RefToBase<reco::Jet>, edm::RefToBase<reco::Jet>, bool> {
    inline bool operator () (const edm::RefToBase<reco::Jet> &j1,
			     const edm::RefToBase<reco::Jet> &j2) const
    { return j1.id() < j2.id() || (j1.id() == j2.id() && j1.key() < j2.key()); }
  };
  
  typedef std::map<edm::RefToBase<reco::Jet>, unsigned int, JetRefCompare> FlavourMap;
  
  // ----------member data ---------------------------
  TH1D * SAnal_Inv_Mass_W_Data0;        
  TH1D * SAnal_Inv_Mass_W__Data0;        
  TH1D * SAnal_Inv_Mass_t_Data0;        
  TH1D * SAnal_Inv_Mass_t__Data0;        
  TH1D * SAnal_Inv_Mass_t_bData0;        
  TH1D * SAnal_Inv_Mass_t__bData0; 




  TH1D * jets;
  TH1D * in;
  TH1D * inb;
  TH1D * TotalEvents;

  TH1D * h_NumberOfJets;
  TH1D * h_NumberOfJets_300;
  TH1D * h_NumberOfJets__300;
    
  TH1D * MatchedSelection;
  TH1D * Reco_Inv_Mass_top;
  TH1D * Reco_Inv_Mass_top_;
  
  TH1D * Reco_Inv_Mass_b;
  TH1D * Reco_Inv_Mass_b_;
  
  TH1D * Reco_Inv_Mass_W;
  TH1D * Reco_Inv_Mass_W_;
  
  TH1D * Gen_Inv_Mass_top;
  TH1D * Gen_Inv_Mass_top_;
  
  TH1D * Gen_Inv_Mass_b;
  TH1D * Gen_Inv_Mass_b_;
  
  TH1D * Gen_Inv_Mass_W;
  TH1D * Gen_Inv_Mass_W_;
  
  TH1D * discriminantbTotal;
  TH1D * discriminantclTotal;

  TH1D * Difference_InPhi_MET_Gen_Reco;  
  TH1D * Difference_InEta_MET_Gen_Reco;  
  TH1D * Difference_InPt_MET_Gen_Reco; 
  TH1D * Difference_InMass_MET_Gen_Reco; 
  TH1D * Difference_InPhi_Muon_Gen_Reco;  
  TH1D * Difference_InEta_Muon_Gen_Reco;  
  TH1D * Difference_InPt_Muon_Gen_Reco; 
  TH1D * Difference_InMass_Muon_Gen_Reco; 

  TH1D * Difference_InEnergy_bJets_Gen_Reco;
  TH1D * Difference_InPt_bJets_Gen_Reco;
  TH1D * Difference_InEta_bJets_Gen_Reco;
  TH1D * Difference_InPhi_bJets_Gen_Reco;  
  TH1D * Difference_InEnergy_bJets_Gen_Reco_20;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_40;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_60;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_80;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_100;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_120;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_140;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_160;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_180;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_200;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_220;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_240;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_260;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_280;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_300;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_320;

  TH1D * Difference_InEnergy_clJets_Gen_Reco;
  TH1D * Difference_InPt_clJets_Gen_Reco;
  TH1D * Difference_InEta_clJets_Gen_Reco;
  TH1D * Difference_InPhi_clJets_Gen_Reco; 
  TH1D * Difference_InEnergy_clJets_Gen_Reco_20;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_40;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_60;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_80;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_100;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_120;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_140;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_160;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_180;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_200;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_220;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_240;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_260;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_280;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_300;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_320;


  TH1D * Sel_Inv_Mass_b_Data0;
  TH1D * Sel_Inv_Mass_b__Data0;
  TH1D * Sel_Inv_Mass_W_Data0;
  TH1D * Sel_Inv_Mass_W__Data0;
  TH1D * Sel_Inv_Mass_t_Data0;
  TH1D * Sel_Inv_Mass_t__Data0;
  TH1D * Sel_Inv_Mass_t_Data1;
  TH1D * Sel_Inv_Mass_t__Data1;
  TH1D * Sel_Inv_Mass_b_Data1;
  TH1D * Sel_Inv_Mass_b__Data1;
  TH1D * Sel_Inv_Mass_W_Data1;
  TH1D * Sel_Inv_Mass_W__Data1;

  TH1D * deltaPhiJetLepbLepS;
  TH1D * deltaPhiJetMETbLepS;
  TH1D * deltaPhiJetLepMETbLepS;
  TH1D * deltaPhiJetLepMETbHadS;

  TH1D * deltaPhiJetQbHadS;
  TH1D * deltaPhiJetQ_bHadS;
  TH1D * deltaPhiJetQQ_bHadS;


  TH1D *  anti_bTagging;
  TH1D *  deltaR_bJetLep;
  TH1D *  et_proportion;
  TH1D *  transverse_momentum;
  TH1D *  mt;


  
/*0404

  TH1D * discriminant2_inclusive;
  TH1D * discriminantb2_inclusive;
  TH1D * discriminantb_2_inclusive;
0404*/  



/*0404  
  TH1D * discriminantbTotalS;
  TH1D * discriminantclTotalS;

  TH1D * bTagDiscS;
0404*/

/*0404  
  TH1D * Difference_InPhi_MET_Gen_RecoS;  
  TH1D * Difference_InEta_MET_Gen_RecoS;  
  TH1D * Difference_InPt_MET_Gen_RecoS; 
  TH1D * Difference_InMass_MET_Gen_RecoS; 
  TH1D * Difference_InPhi_Muon_Gen_RecoS;  
  TH1D * Difference_InEta_Muon_Gen_RecoS;  
  TH1D * Difference_InPt_Muon_Gen_RecoS; 
  TH1D * Difference_InMass_Muon_Gen_RecoS; 
0404*/  

/*0404
  TH1D * Difference_InEnergy_bJets_Gen_RecoS;
  TH1D * Difference_InPt_bJets_Gen_RecoS;
  TH1D * Difference_InEta_bJets_Gen_RecoS;
  TH1D * Difference_InPhi_bJets_Gen_RecoS;  
  TH1D * Difference_InEnergy_bJets_Gen_Reco_20S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_40S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_60S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_80S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_100S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_120S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_140S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_160S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_180S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_200S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_220S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_240S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_260S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_280S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_300S;
  TH1D * Difference_InEnergy_bJets_Gen_Reco_320S;
0404*/
  /* en caso de ser necesario hacer una funcion de transferencia tanto en eta como en energia, tomado de Top Analysis Group Tool.  La idea es que apenas se definen, falta hacer el constructor y llenarlas, ademas faltan las mismas para cl.
     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0087;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0087;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0087;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0174;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0174;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0174;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0261;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0261;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0261;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0348;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0348;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0348;


     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0435;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0435;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0435;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0522;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0522;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0522;
 

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0609;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0609;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0609;
 

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0696;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0696;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0696;


     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0783;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0783;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0783;
 

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0870;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0870;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0870;


     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta0957;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta0957;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta0957;


     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1044;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1044;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1044;


     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1131;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1131;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1131;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1218;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1218;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1218;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1305;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1305;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1305;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1392;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1392;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1392;
  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1479;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1479;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1479;
  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1566;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1566;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1566;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1653;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1653;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1653;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1740;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1740;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1740;

     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1830;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1830;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1830;
  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta1930;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta1930;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta1930;
   
     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta2043;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta2043;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta2043;
 
     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta2172;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta2172;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta2172;
 
     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta2322;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta2322;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta2322;


     TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta2500;  
     TH1D * Difference_InEnergy_bJets_Gen_Reco_20_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_40_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_60_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_80_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_100_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_120_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_140_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_160_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_180_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_200_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_220_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_240_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_260_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_280_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_300_Eta2500;
     TH1D * Difference_InEnergy_bJets_Gen_Reco_320_Eta2500;
  */


  //Inicio: estan que vienen a continuacion es para hacer la funcion de transferencia en eta, tal como la hacen en el grupo de tops, podria ser el paso a seguir si no funciona solo con energia,  y si con esta tampoco, seria incluir la de MET y leptons y si aun asi no, seria incluir la mezcla de energia y eta, es decir, lo que esta comentado arriba
  //falta hacerlo pero la idea es hacer:
  // TH1D * Difference_InEnergy_bJets_Gen_Reco_Eta2500;
  //para cada uno de los eta     
  //Inicio: estan que vienen a continuacion es para hacer la funcion de transferencia en eta, tal como la hacen en el grupo de tops, podria ser el paso a seguir si no funciona solo con energia,  y si con esta tampoco, seria incluir la de MET y leptons y si aun asi no, seria incluir la mezcla de energia y eta, es decir, lo que esta comentado arriba

/*0404  
  TH1D * Difference_InEnergy_clJets_Gen_RecoS;
  TH1D * Difference_InPt_clJets_Gen_RecoS;
  TH1D * Difference_InEta_clJets_Gen_RecoS;
  TH1D * Difference_InPhi_clJets_Gen_RecoS; 
  TH1D * Difference_InEnergy_clJets_Gen_Reco_20S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_40S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_60S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_80S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_100S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_120S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_140S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_160S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_180S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_200S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_220S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_240S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_260S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_280S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_300S;
  TH1D * Difference_InEnergy_clJets_Gen_Reco_320S;
  
  
  TH1D * numberOfJetsInPhi;
  TH1D * inSel;
0404*/  
/*0404  
  TH1D * Scaled_Gen_Inv_Mass_top;
  TH1D * Scaled_Gen_Inv_Mass_top_;
  
  TH1D * Scaled_Gen_Inv_Mass_b;
  TH1D * Scaled_Gen_Inv_Mass_b_;
  
  TH1D * Scaled_Gen_Inv_Mass_W;
  TH1D * Scaled_Gen_Inv_Mass_W_;
  
  
  TH2D * METvsNumberOfJets;
  TH1D * gProportionMETJets;
  TH1D * ProportionMETJets;

 
  TH1D * Sel_Inv_Mass_b;
  TH1D * Sel_Inv_Mass_b_;
  TH1D * Sel_Inv_Mass_W;
0404*/

/*0404  
  TH1D * deltaR_q;
  TH1D * deltaR_q_;
  TH1D * deltaR_b;
  TH1D * deltaR_b_;
  
  TH1D * sel_deltaR_q;
  TH1D * sel_deltaR_q_;
  TH1D * sel_deltaR_b;
  TH1D * sel_deltaR_b_;
  TH2D * discriminantJPvsN;
  TH2D * deltaPhiJetLep_vs_deltaPhiJetMET_bLep;
  TH2D * deltaPhiJetQ_vs_deltaPhiJetQ__bHad;

  TH1D * deltaPhiJetLepbLep;
  TH1D * deltaPhiJetMETbLep;
  TH1D * deltaPhiJetLepMETbLep;

  TH1D * deltaPhiJetQbHad;
  TH1D * deltaPhiJetQ_bHad;
  TH1D * deltaPhiJetQQ_bHad;





  TH1D * discriminantJPSelection;
  TH1D * discriminantJPMatching;
  TH2D * discriminantJPvsN_Data;
  TH1D * discriminantJPSelection_Data;
0404*/
  /*
    TH1D * ElDiscriminante1_Data0;
    TH1D * ElDiscriminante2_Data0;
    TH1D * ElDiscriminante1_Data1;
    TH1D * ElDiscriminante2_Data1;
  */
/*0404
  TH1D * ElDiscriminante1;
  TH1D * ElDiscriminante2;
0404*/
  /*
    TH1D * Discriminant_lep_W_10_10;
    TH1D * Discriminant_lep_t_10_10;
    TH1D * Discriminant_had_W_10_10;
    TH1D * Discriminant_had_t_10_10;

    TH1D * Discriminant_lep_W_30_30;
    TH1D * Discriminant_lep_t_30_30;
    TH1D * Discriminant_had_W_30_30;
    TH1D * Discriminant_had_t_30_30;

    TH1D * Discriminant_lep_W_50_50;
    TH1D * Discriminant_lep_t_50_50;
    TH1D * Discriminant_had_W_50_50;
    TH1D * Discriminant_had_t_50_50;

    TH1D * Discriminant_lep_W_70_70;
    TH1D * Discriminant_lep_t_70_70;
    TH1D * Discriminant_had_W_70_70;
    TH1D * Discriminant_had_t_70_70;

    TH1D * Discriminant_lep_W_90_90;
    TH1D * Discriminant_lep_t_90_90;
    TH1D * Discriminant_had_W_90_90;
    TH1D * Discriminant_had_t_90_90;


    TH1D * Discriminant_lep_W_10_90;
    TH1D * Discriminant_lep_t_10_90;
    TH1D * Discriminant_had_W_10_90;
    TH1D * Discriminant_had_t_10_90;

    TH1D * Discriminant_lep_W_30_70;
    TH1D * Discriminant_lep_t_30_70;
    TH1D * Discriminant_had_W_30_70;
    TH1D * Discriminant_had_t_30_70;

    TH1D * Discriminant_lep_W_70_30;
    TH1D * Discriminant_lep_t_70_30;
    TH1D * Discriminant_had_W_70_30;
    TH1D * Discriminant_had_t_70_30;

    TH1D * Discriminant_lep_W_90_10;
    TH1D * Discriminant_lep_t_90_10;
    TH1D * Discriminant_had_W_90_10;
    TH1D * Discriminant_had_t_90_10;

    TH1D * Discriminant_total_t_10_10;
    TH1D * Discriminant_total_W_10_10;
    TH1D * Discriminant_total_t_30_30;
    TH1D * Discriminant_total_W_30_30;
    TH1D * Discriminant_total_t_50_50;
    TH1D * Discriminant_total_W_50_50;
    TH1D * Discriminant_total_t_70_70;
    TH1D * Discriminant_total_W_70_70;
    TH1D * Discriminant_total_t_90_90;
    TH1D * Discriminant_total_W_90_90;
    TH1D * Discriminant_total_t_10_90;
    TH1D * Discriminant_total_W_10_90;
    TH1D * Discriminant_total_t_30_70;
    TH1D * Discriminant_total_W_30_70;
    TH1D * Discriminant_total_t_70_30;
    TH1D * Discriminant_total_W_70_30;
    TH1D * Discriminant_total_t_90_10;
    TH1D * Discriminant_total_W_90_10;
  */


/*0404
  TH1D * Discriminant_lep_W_90_90;
  TH1D * Discriminant_lep_t_90_90;
  TH1D * Discriminant_had_W_90_90;
  TH1D * Discriminant_had_t_90_90;
  TH1D * Discriminant_total_t_90_90;
  TH1D * Discriminant_total_W_90_90;

  TH1D * Discriminant_lep_W_85_85;
  TH1D * Discriminant_lep_t_85_85;
  TH1D * Discriminant_had_W_85_85;
  TH1D * Discriminant_had_t_85_85;
  TH1D * Discriminant_total_t_85_85;
  TH1D * Discriminant_total_W_85_85;

  TH1D * Discriminant_lep_W_85_95;
  TH1D * Discriminant_lep_t_85_95;
  TH1D * Discriminant_had_W_85_95;
  TH1D * Discriminant_had_t_85_95;
  TH1D * Discriminant_total_t_85_95;
  TH1D * Discriminant_total_W_85_95;

  TH1D * Discriminant_lep_W_95_85;
  TH1D * Discriminant_lep_t_95_85;
  TH1D * Discriminant_had_W_95_85;
  TH1D * Discriminant_had_t_95_85;
  TH1D * Discriminant_total_t_95_85;
  TH1D * Discriminant_total_W_95_85;

  TH1D * Discriminant_lep_W_95_95;
  TH1D * Discriminant_lep_t_95_95;
  TH1D * Discriminant_had_W_95_95;
  TH1D * Discriminant_had_t_95_95;
  TH1D * Discriminant_total_t_95_95;
  TH1D * Discriminant_total_W_95_95;
0404*/

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FullAndWorking::FullAndWorking(const edm::ParameterSet& iConfig)
  
{
  edm::Service<TFileService> fs;
  top=iConfig.getUntrackedParameter<bool>("top",true);
  RFijo=iConfig.getUntrackedParameter<bool>("RFijo",true);
  LHCOWithMC=iConfig.getUntrackedParameter<bool>("LHCOWithMC",false);
  MC=iConfig.getUntrackedParameter<bool>("MC",true);
  neutralino=iConfig.getUntrackedParameter<bool>("neutralino",true); 
  bdisc_name=iConfig.getParameter<std::string>("bdisc_name");
  selector=iConfig.getUntrackedParameter<double>("selector");
  c1=iConfig.getUntrackedParameter<double>("c1");
  c2=iConfig.getUntrackedParameter<double>("c2");
  c3=iConfig.getUntrackedParameter<double>("c3");
  c4=iConfig.getUntrackedParameter<double>("c4");
  c5=iConfig.getUntrackedParameter<double>("c5");
  c6=iConfig.getUntrackedParameter<double>("c6");
  c7=iConfig.getUntrackedParameter<double>("c7");
  c8=iConfig.getUntrackedParameter<double>("c8");
  c9=iConfig.getUntrackedParameter<double>("c9");
  c10=iConfig.getUntrackedParameter<double>("c10");
  c11=iConfig.getUntrackedParameter<double>("c11");
  c12=iConfig.getUntrackedParameter<double>("c12");
  c13=iConfig.getUntrackedParameter<double>("c13");
  c14=iConfig.getUntrackedParameter<double>("c14");
  
  in= fs->make<TH1D>("h_Matched_Jets_are_in_the_set of_--_Jets_with_biggest_pt" , "h_Matched_Jets_are_in_the_set_of_--_Jets_with_biggest_pt" , 15 , -5 , 10 );
  inb= fs->make<TH1D>("h_Matched_bJets_are_in_the_set_of_--_Jets_with_biggest_discriminator" , "h_Matched_Jets_are_in_set_of_--_Jets_with_biggest_discriminator" , 15 , -5 , 10 );
  TotalEvents= fs->make<TH1D>("h_Number_of_Total_Events" , "h_Number_of_Total_Events" , 3 , -1 , 1 );
  h_NumberOfJets= fs->make<TH1D>("NumberOfJets" , "NumberOfJets" , 100 , -1 , 98 );
  h_NumberOfJets_300= fs->make<TH1D>("NumberOfJets_300" , "NumberOfJets_300" , 100 , -1 , 98 );
  h_NumberOfJets__300= fs->make<TH1D>("NumberOfJets__300" , "NumberOfJets__300" , 100 , -1 , 98 );
  
  MatchedSelection = fs->make<TH1D>("h_MatchedSelection" , "h_MatchedSelection" , 10020 , -20 , 10000 );

  discriminantbTotal = fs->make<TH1D>("h_Discriminant_b_Total" , "h_Discriminant_b_Total" , 300 , -1 , 2 );
  discriminantclTotal = fs->make<TH1D>("h_Discriminant_cl_Total" , "h_Discriminant_cl_Total" , 300 , -1 , 2 );
  Difference_InPhi_Muon_Gen_Reco = fs->make<TH1D>("h_Difference_InPhi_Muon_Gen_Reco" , "h_Difference_InPhi_Muon_Gen_Reco" , 400 , -6.5 , 6.5 );   
  Difference_InPt_Muon_Gen_Reco = fs->make<TH1D>("h_Difference_InPt_Muon_Gen_Reco" , "h_Difference_InPt_Muon_Gen_Reco" , 4000 , -200 , 200 );   
  Difference_InEta_Muon_Gen_Reco = fs->make<TH1D>("h_Difference_InEta_Muon_Gen_Reco" , "h_Difference_InEta_Muon_Gen_Reco" , 400 , -6.5 , 6.5 );
  Difference_InMass_Muon_Gen_Reco = fs->make<TH1D>("h_Difference_InMass_Muon_Gen_Reco" , "h_Difference_InMass_Muon_Gen_Reco" , 400 , -6.5 , 6.5 );
  Difference_InPhi_MET_Gen_Reco = fs->make<TH1D>("h_Difference_InPhi_MET_Gen_Reco" , "h_Difference_InPhi_MET_Gen_Reco" , 400 , -6.5 , 6.5 );   
  Difference_InPt_MET_Gen_Reco = fs->make<TH1D>("h_Difference_InPt_MET_Gen_Reco" , "h_Difference_InPt_MET_Gen_Reco" , 4000 , -200 , 200 );   
  Difference_InEta_MET_Gen_Reco = fs->make<TH1D>("h_Difference_InEta_MET_Gen_Reco" , "h_Difference_InEta_MET_Gen_Reco" , 400 , -6.5 , 6.5 );   
  Difference_InMass_MET_Gen_Reco = fs->make<TH1D>("h_Difference_InMass_MET_Gen_Reco" , "h_Difference_InMass_MET_Gen_Reco" , 400 , -6.5 , 6.5 );   
  
  Difference_InEnergy_bJets_Gen_Reco = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco" , "h_Difference_InEnergy_bJets_Gen_Reco" , 4000 , -200 , 200 );
  Difference_InPt_bJets_Gen_Reco = fs->make<TH1D>("h_Difference_InPt_bJets_Gen_Reco" , "h_Difference_InPt_bJets_Gen_Reco" , 4000 , -200 , 200 );
  Difference_InEta_bJets_Gen_Reco = fs->make<TH1D>("h_Difference_InEta_bJets_Gen_Reco" , "h_Difference_InEta_bJets_Gen_Reco" , 400 , -20 , 20 );
  Difference_InPhi_bJets_Gen_Reco = fs->make<TH1D>("h_Difference_InPhi_bJets_Gen_Reco" , "h_Difference_InPhi_bJets_Gen_Reco" , 400 , -20 , 20 );
  Difference_InEnergy_bJets_Gen_Reco_20 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_20" , "h_Difference_InEnergy_bJets_Gen_Reco_20" , 4000 , -200 , 200 );   
  Difference_InEnergy_bJets_Gen_Reco_40 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_40" , "h_Difference_InEnergy_bJets_Gen_Reco_40" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_60 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_60" , "h_Difference_InEnergy_bJets_Gen_Reco_60" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_80 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_80" , "h_Difference_InEnergy_bJets_Gen_Reco_80" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_100 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_100" , "h_Difference_InEnergy_bJets_Gen_Reco_100" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_120 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_120" , "h_Difference_InEnergy_bJets_Gen_Reco_120" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_140 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_140" , "h_Difference_InEnergy_bJets_Gen_Reco_140" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_160 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_160" , "h_Difference_InEnergy_bJets_Gen_Reco_160" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_180 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_180" , "h_Difference_InEnergy_bJets_Gen_Reco_180" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_200 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_200" , "h_Difference_InEnergy_bJets_Gen_Reco_200" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_220 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_220" , "h_Difference_InEnergy_bJets_Gen_Reco_220" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_240 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_240" , "h_Difference_InEnergy_bJets_Gen_Reco_240" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_260 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_260" , "h_Difference_InEnergy_bJets_Gen_Reco_260" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_280 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_280" , "h_Difference_InEnergy_bJets_Gen_Reco_280" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_300 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_300" , "h_Difference_InEnergy_bJets_Gen_Reco_300" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_320 = fs->make<TH1D>("h_Difference_InEnergy_bJets_Gen_Reco_320" , "h_Difference_InEnergy_bJets_Gen_Reco_320" , 4000 , -200 , 200 );
  
  Difference_InEnergy_clJets_Gen_Reco = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco" , "h_Difference_InEnergy_clJets_Gen_Reco" , 4000 , -200 , 200 );
  Difference_InPt_clJets_Gen_Reco = fs->make<TH1D>("h_Difference_InPt_clJets_Gen_Reco" , "h_Difference_InPt_clJets_Gen_Reco" , 4000 , -200 , 200 );
  Difference_InEta_clJets_Gen_Reco = fs->make<TH1D>("h_Difference_InEta_clJets_Gen_Reco" , "h_Difference_InEta_clJets_Gen_Reco" , 400 , -20 , 20 );
  Difference_InPhi_clJets_Gen_Reco = fs->make<TH1D>("h_Difference_InPhi_clJets_Gen_Reco" , "h_Difference_InPhi_clJets_Gen_Reco" , 400 , -20 , 20 );
  Difference_InEnergy_clJets_Gen_Reco_20 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_20" , "h_Difference_InEnergy_clJets_Gen_Reco_20" , 4000 , -200 , 200 );   
  Difference_InEnergy_clJets_Gen_Reco_40 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_40" , "h_Difference_InEnergy_clJets_Gen_Reco_40" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_60 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_60" , "h_Difference_InEnergy_clJets_Gen_Reco_60" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_80 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_80" , "h_Difference_InEnergy_clJets_Gen_Reco_80" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_100 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_100" , "h_Difference_InEnergy_clJets_Gen_Reco_100" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_120 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_120" , "h_Difference_InEnergy_clJets_Gen_Reco_120" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_140 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_140" , "h_Difference_InEnergy_clJets_Gen_Reco_140" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_160 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_160" , "h_Difference_InEnergy_clJets_Gen_Reco_160" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_180 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_180" , "h_Difference_InEnergy_clJets_Gen_Reco_180" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_200 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_200" , "h_Difference_InEnergy_clJets_Gen_Reco_200" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_220 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_220" , "h_Difference_InEnergy_clJets_Gen_Reco_220" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_240 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_240" , "h_Difference_InEnergy_clJets_Gen_Reco_240" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_260 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_260" , "h_Difference_InEnergy_clJets_Gen_Reco_260" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_280 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_280" , "h_Difference_InEnergy_clJets_Gen_Reco_280" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_300 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_300" , "h_Difference_InEnergy_clJets_Gen_Reco_300" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_320 = fs->make<TH1D>("h_Difference_InEnergy_clJets_Gen_Reco_320" , "h_Difference_InEnergy_clJets_Gen_Reco_320" , 4000 , -200 , 200 );
  
  Reco_Inv_Mass_top = fs->make<TH1D>("h_Reco_Inv_Mass_top" , "h_Reco_Inv_Mass_top" , 150 , 100 , 500 );
  Reco_Inv_Mass_top_ = fs->make<TH1D>("h_Reco_Inv_Mass_top_" , "h_Reco_Inv_Mass_top_" , 150 , 100 , 500 );
  Reco_Inv_Mass_b = fs->make<TH1D>("h_Reco_Inv_Mass_b" , "h_Reco_Inv_Mass_b" , 40 , 0 , 25 );
  Reco_Inv_Mass_b_ = fs->make<TH1D>("h_Reco_Inv_Mass_b_" , "h_Reco_Inv_Mass_b_" , 40 , 0 , 25 );
  Reco_Inv_Mass_W = fs->make<TH1D>("h_Reco_Inv_Mass_W" , "h_Reco_Inv_Mass_W" , 40 , 10 , 150 );
  Reco_Inv_Mass_W_ = fs->make<TH1D>("h_Reco_Inv_Mass_W_" , "h_Reco_Inv_Mass_W_" , 40 , 10 , 150 );

  anti_bTagging = fs->make<TH1D>("h_anti_bTagging" , "h_anti_bTagging" , 2000 , 0 , 10000 );
  deltaR_bJetLep = fs->make<TH1D>("h_deltaR_bJetLep" , "h_deltaR_bJetLep" , 50 ,0 , 100 ); 
  et_proportion = fs->make<TH1D>("h_et_proportion" , "h_et_proportion" , 20 ,0 , 2 );
  transverse_momentum = fs->make<TH1D>("h_transverse_momentum" , "h_transverse_momentum" , 2000 ,-100 , 100 );
  mt = fs->make<TH1D>("h_mt" , "h_mt" , 20 ,0 , 300 );
  
  Gen_Inv_Mass_top = fs->make<TH1D>("h_Gen_Inv_Mass_top" , "h_Gen_Inv_Mass_top" , 40 , 100 , 250 );
  Gen_Inv_Mass_top_ = fs->make<TH1D>("h_Gen_Inv_Mass_top_" , "h_Gen_Inv_Mass_top_" , 40 , 100 , 250 );
  Gen_Inv_Mass_b = fs->make<TH1D>("h_Gen_Inv_Mass_b" , "h_Gen_Inv_Mass_b" , 40 , 0 , 25 );
  Gen_Inv_Mass_b_ = fs->make<TH1D>("h_Gen_Inv_Mass_b_" , "h_Gen_Inv_Mass_b_" , 40 , 0 , 25 );
  Gen_Inv_Mass_W = fs->make<TH1D>("h_Gen_Inv_Mass_W" , "h_Gen_Inv_Mass_W" , 40 , 10 , 150 );
  Gen_Inv_Mass_W_ = fs->make<TH1D>("h_Gen_Inv_Mass_W_" , "h_Gen_Inv_Mass_W_" , 40 , 10 , 150 );


   
  Sel_Inv_Mass_b_Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_b_Data0" , "h_Sel_Inv_Mass_b_Data0" , 50 , 0 , 25 );
  Sel_Inv_Mass_b__Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_b__Data0" , "h_Sel_Inv_Mass_b__Data0" , 50 , 0 , 25 );
  Sel_Inv_Mass_W_Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_W_Data0" , "h_Sel_Inv_Mass_W_Data0" , 50 , 0 , 200 );
  Sel_Inv_Mass_W__Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_W__Data0" , "h_Sel_Inv_Mass_W__Data0" , 50 , 0 , 200 );
  Sel_Inv_Mass_t_Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_t_Data0" , "h_Sel_Inv_Mass_t_Data0" , 50 , 0, 400 );
  Sel_Inv_Mass_t__Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_t__Data0" , "h_Sel_Inv_Mass_t__Data0" , 50 , 0 , 400 );




  SAnal_Inv_Mass_W_Data0 = fs->make<TH1D>("h_SAnal_Inv_Mass_W_Data0" , "h_SAnal_Inv_Mass_W_Data0" , 100 , 0 , 1000 );
  SAnal_Inv_Mass_W__Data0 = fs->make<TH1D>("h_SAnal_Inv_Mass_W__Data0" , "h_SAnal_Inv_Mass_W__Data0" , 100 , 0 , 1000 );
  SAnal_Inv_Mass_t_Data0 = fs->make<TH1D>("h_SAnal_Inv_Mass_t_Data0" , "h_SAnal_Inv_Mass_t_Data0" , 100 , 0 , 1000 );
  SAnal_Inv_Mass_t__Data0 = fs->make<TH1D>("h_SAnal_Inv_Mass_t__Data0" , "h_SAnal_Inv_Mass_t__Data0" , 100 , 0 , 1000 );
  SAnal_Inv_Mass_t_bData0 = fs->make<TH1D>("h_SAnal_Inv_Mass_t_bData0" , "h_SAnal_Inv_Mass_t_bData0" , 100 , 0, 1000 );
  SAnal_Inv_Mass_t__bData0 = fs->make<TH1D>("h_SAnal_Inv_Mass_t__bData0" , "h_SAnal_Inv_Mass_t__bData0" , 100 , 0 , 1000 );
  

 

  deltaPhiJetLepbLepS= fs->make<TH1D>("h_deltaPhiJetLep_bLepS" , "h_deltaPhiJetLep_bLepS" , 30 , 0 , 10);
  deltaPhiJetMETbLepS= fs->make<TH1D>("h_deltaPhiJetMET_bLepS" , "h_deltaPhiJetMET_bLepS" , 30 , 0 , 10);
  deltaPhiJetLepMETbLepS= fs->make<TH1D>("h_deltaPhiJetLepMET_bLepS" , "h_deltaPhiJetLepMET_bLepS" , 30 , 0 , 10);
  deltaPhiJetLepMETbHadS= fs->make<TH1D>("h_deltaPhiJetLepMET_bHadS" , "h_deltaPhiJetLepMET_bHadS" , 30 , 0 , 10);

  deltaPhiJetQbHadS= fs->make<TH1D>("h_deltaPhiJetQ_bHadS" , "h_deltaPhiJetQ_bHadS" , 30 , 0 , 10);
  deltaPhiJetQ_bHadS= fs->make<TH1D>("h_deltaPhiJetQ__bHadS" , "h_deltaPhiJetQ__bHadS" , 30 , 0 , 10);
  deltaPhiJetQQ_bHadS= fs->make<TH1D>("h_deltaPhiJetQQ__bHadS" , "h_deltaPhiJetQQ__bHadS" , 30 , 0 , 10);


  jets = fs->make<TH1D>("h_jets" , "h_jets" , 20 , -0 , 20 );

/*0404
  //now do what ever initialization is needed
  //  bTagAlgo=iConfig.getParameter<edm::InputTag>("bTagAlgo");
  //matchingAlgo=iConfig.getParameter<edm::InputTag>("matchingAlgo");

  Sel_Inv_Mass_b = fs->make<TH1D>("Sel_Inv_Mass_b" , "Sel_Inv_Mass_b" , 40 , 0 , 25 );
  Sel_Inv_Mass_b_ = fs->make<TH1D>("Sel_Inv_Mass_b_" , "Sel_Inv_Mass_b_" , 40 , 0 , 25 );
  Sel_Inv_Mass_W = fs->make<TH1D>("Sel_Inv_Mass_W" , "Sel_Inv_Mass_W" , 40 , 10 , 150 );
  
  discriminant2_inclusive = fs->make<TH1D>("Discriminant2_inclusive" , "Discriminant2_inclusive" , 300 , -1 , 2 );
  discriminantb2_inclusive = fs->make<TH1D>("Discriminant_inclusive b-Jets2" , "Discriminant_inclusive b-Jets2" , 300 , -1 , 2 );
  discriminantb_2_inclusive = fs->make<TH1D>("Discriminant_inclusive _b-Jets2" , "Discriminant_inclusive _b-Jets2" , 300 , -1 , 2 );
  

  bTagDiscS= fs->make<TH1D>("bTagDiscS" , "bTagDiscS " , 400 , -1 , 2000 );



  
  ////JPPPP
  discriminantbTotalS = fs->make<TH1D>("Discriminant b TotalS" , "Discriminant b TotalS" , 300 , -1 , 2 );
  discriminantclTotalS = fs->make<TH1D>("Discriminant cl TotalS" , "Discriminant cl TotalS" , 300 , -1 , 2 );
  
  Difference_InPhi_Muon_Gen_RecoS = fs->make<TH1D>("Difference_InPhi_Muon_Gen_RecoS" , "Difference_InPhi_Muon_Gen_RecoS" , 400 , -6.5 , 6.5 );   
  Difference_InPt_Muon_Gen_RecoS = fs->make<TH1D>("Difference_InPt_Muon_Gen_RecoS" , "Difference_InPt_Muon_Gen_RecoS" , 4000 , -200 , 200 );   
  Difference_InEta_Muon_Gen_RecoS = fs->make<TH1D>("Difference_InEta_Muon_Gen_RecoS" , "Difference_InEta_Muon_Gen_RecoS" , 400 , -6.5 , 6.5 );
  Difference_InMass_Muon_Gen_RecoS = fs->make<TH1D>("Difference_InMass_Muon_Gen_RecoS" , "Difference_InMass_Muon_Gen_RecoS" , 400 , -6.5 , 6.5 );
  Difference_InPhi_MET_Gen_RecoS = fs->make<TH1D>("Difference_InPhi_MET_Gen_RecoS" , "Difference_InPhi_MET_Gen_RecoS" , 400 , -6.5 , 6.5 );   
  Difference_InPt_MET_Gen_RecoS = fs->make<TH1D>("Difference_InPt_MET_Gen_RecoS" , "Difference_InPt_MET_Gen_RecoS" , 4000 , -200 , 200 );   
  Difference_InEta_MET_Gen_RecoS = fs->make<TH1D>("Difference_InEta_MET_Gen_RecoS" , "Difference_InEta_MET_Gen_RecoS" , 400 , -6.5 , 6.5 );   
  Difference_InMass_MET_Gen_RecoS = fs->make<TH1D>("Difference_InMass_MET_Gen_RecoS" , "Difference_InMass_MET_Gen_RecoS" , 400 , -6.5 , 6.5 );   
  
  Difference_InEnergy_bJets_Gen_RecoS = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_RecoS" , "Difference_InEnergy_bJets_Gen_RecoS" , 4000 , -200 , 200 );
  Difference_InPt_bJets_Gen_RecoS = fs->make<TH1D>("Difference_InPt_bJets_Gen_RecoS" , "Difference_InPt_bJets_Gen_RecoS" , 4000 , -200 , 200 );
  Difference_InEta_bJets_Gen_RecoS = fs->make<TH1D>("Difference_InEta_bJets_Gen_RecoS" , "Difference_InEta_bJets_Gen_RecoS" , 400 , -20 , 20 );
  Difference_InPhi_bJets_Gen_RecoS = fs->make<TH1D>("Difference_InPhi_bJets_Gen_RecoS" , "Difference_InPhi_bJets_Gen_RecoS" , 400 , -20 , 20 );
  Difference_InEnergy_bJets_Gen_Reco_20S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_20S" , "Difference_InEnergy_bJets_Gen_Reco_20S" , 4000 , -200 , 200 );   
  Difference_InEnergy_bJets_Gen_Reco_40S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_40S" , "Difference_InEnergy_bJets_Gen_Reco_40S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_60S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_60S" , "Difference_InEnergy_bJets_Gen_Reco_60S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_80S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_80S" , "Difference_InEnergy_bJets_Gen_Reco_80S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_100S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_100S" , "Difference_InEnergy_bJets_Gen_Reco_100S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_120S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_120S" , "Difference_InEnergy_bJets_Gen_Reco_120S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_140S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_140S" , "Difference_InEnergy_bJets_Gen_Reco_140S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_160S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_160S" , "Difference_InEnergy_bJets_Gen_Reco_160S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_180S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_180S" , "Difference_InEnergy_bJets_Gen_Reco_180S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_200S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_200S" , "Difference_InEnergy_bJets_Gen_Reco_200S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_220S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_220S" , "Difference_InEnergy_bJets_Gen_Reco_220S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_240S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_240S" , "Difference_InEnergy_bJets_Gen_Reco_240S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_260S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_260S" , "Difference_InEnergy_bJets_Gen_Reco_260S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_280S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_280S" , "Difference_InEnergy_bJets_Gen_Reco_280S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_300S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_300S" , "Difference_InEnergy_bJets_Gen_Reco_300S" , 4000 , -200 , 200 );
  Difference_InEnergy_bJets_Gen_Reco_320S = fs->make<TH1D>("Difference_InEnergy_bJets_Gen_Reco_320S" , "Difference_InEnergy_bJets_Gen_Reco_320S" , 4000 , -200 , 200 );
  
  Difference_InEnergy_clJets_Gen_RecoS = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_RecoS" , "Difference_InEnergy_clJets_Gen_RecoS" , 4000 , -200 , 200 );
  Difference_InPt_clJets_Gen_RecoS = fs->make<TH1D>("Difference_InPt_clJets_Gen_RecoS" , "Difference_InPt_clJets_Gen_RecoS" , 4000 , -200 , 200 );
  Difference_InEta_clJets_Gen_RecoS = fs->make<TH1D>("Difference_InEta_clJets_Gen_RecoS" , "Difference_InEta_clJets_Gen_RecoS" , 400 , -20 , 20 );
  Difference_InPhi_clJets_Gen_RecoS = fs->make<TH1D>("Difference_InPhi_clJets_Gen_RecoS" , "Difference_InPhi_clJets_Gen_RecoS" , 400 , -20 , 20 );
  Difference_InEnergy_clJets_Gen_Reco_20S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_20S" , "Difference_InEnergy_clJets_Gen_Reco_20S" , 4000 , -200 , 200 );   
  Difference_InEnergy_clJets_Gen_Reco_40S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_40S" , "Difference_InEnergy_clJets_Gen_Reco_40S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_60S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_60S" , "Difference_InEnergy_clJets_Gen_Reco_60S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_80S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_80S" , "Difference_InEnergy_clJets_Gen_Reco_80S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_100S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_100S" , "Difference_InEnergy_clJets_Gen_Reco_100S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_120S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_120S" , "Difference_InEnergy_clJets_Gen_Reco_120S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_140S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_140S" , "Difference_InEnergy_clJets_Gen_Reco_140S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_160S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_160S" , "Difference_InEnergy_clJets_Gen_Reco_160S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_180S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_180S" , "Difference_InEnergy_clJets_Gen_Reco_180S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_200S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_200S" , "Difference_InEnergy_clJets_Gen_Reco_200S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_220S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_220S" , "Difference_InEnergy_clJets_Gen_Reco_220S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_240S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_240S" , "Difference_InEnergy_clJets_Gen_Reco_240S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_260S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_260S" , "Difference_InEnergy_clJets_Gen_Reco_260S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_280S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_280S" , "Difference_InEnergy_clJets_Gen_Reco_280S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_300S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_300S" , "Difference_InEnergy_clJets_Gen_Reco_300S" , 4000 , -200 , 200 );
  Difference_InEnergy_clJets_Gen_Reco_320S = fs->make<TH1D>("Difference_InEnergy_clJets_Gen_Reco_320S" , "Difference_InEnergy_clJets_Gen_Reco_320S" , 4000 , -200 , 200 );
  
  ////JPPPP







  numberOfJetsInPhi= fs->make<TH1D>("number Of Jets In Phi" , "number Of Jets In Phi" , 11 , -1 , 10 );
  inSel= fs->make<TH1D>("-- Matched Jets are in the selection criteria used", "-- Matched Jets are in the selection criteria used" , 15 , -5 , 10 );
  
  
  
  Scaled_Gen_Inv_Mass_top = fs->make<TH1D>("Scaled_Gen_Inv_Mass_top" , "Scaled_Gen_Inv_Mass_top" , 40 , 100 , 250 );
  Scaled_Gen_Inv_Mass_top_ = fs->make<TH1D>("Scaled_Gen_Inv_Mass_top_" , "Scaled_Gen_Inv_Mass_top_" , 40 , 100 , 250 );
  Scaled_Gen_Inv_Mass_b = fs->make<TH1D>("Scaled_Gen_Inv_Mass_b" , "Scaled_Gen_Inv_Mass_b" , 40 , 0 , 25 );
  Scaled_Gen_Inv_Mass_b_ = fs->make<TH1D>("Scaled_Gen_Inv_Mass_b_" , "Scaled_Gen_Inv_Mass_b_" , 40 , 0 , 25 );
  Scaled_Gen_Inv_Mass_W = fs->make<TH1D>("Scaled_Gen_Inv_Mass_W" , "Scaled_Gen_Inv_Mass_W" , 40 , 10 , 150 );
  Scaled_Gen_Inv_Mass_W_ = fs->make<TH1D>("Scaled_Gen_Inv_Mass_W_" , "Scaled_Gen_Inv_Mass_W_" , 40 , 10 , 150 );
  
  
  
  gProportionMETJets = fs->make<TH1D>("gProportionMETJets" , "gProportionMETJets" , 40 , 0 , 1 );
  ProportionMETJets = fs->make<TH1D>("ProportionMETJets" , "ProportionMETJets" , 40 , 0 , 1 );
  METvsNumberOfJets = fs->make<TH2D>("METvsNumberOfJets" , "METvsNumberOfJets" , 105 , -5 , 100, 30, -5, 25 );
  
  

  Sel_Inv_Mass_b_Data1 = fs->make<TH1D>("Sel_Inv_Mass_b_Data1" , "Sel_Inv_Mass_b_Data1" , 40 , 0 , 25 );
  Sel_Inv_Mass_b__Data1 = fs->make<TH1D>("Sel_Inv_Mass_b__Data1" , "Sel_Inv_Mass_b__Data1" , 40 , 0 , 25 );
  Sel_Inv_Mass_W_Data1 = fs->make<TH1D>("Sel_Inv_Mass_W_Data1" , "Sel_Inv_Mass_W_Data1" , 40 , 10 , 150 );
  Sel_Inv_Mass_W__Data1 = fs->make<TH1D>("Sel_Inv_Mass_W__Data1" , "Sel_Inv_Mass_W__Data1" , 40 , 10 , 150 );
  Sel_Inv_Mass_t_Data1 = fs->make<TH1D>("Sel_Inv_Mass_t_Data1" , "Sel_Inv_Mass_t_Data1" , 9000 , 100, 1000 );
  Sel_Inv_Mass_t__Data1 = fs->make<TH1D>("Sel_Inv_Mass_t__Data1" , "Sel_Inv_Mass_t__Data1" , 9000 , 100 , 1000 );


  sel_deltaR_b = fs->make<TH1D>("sel_deltaR_b" , "sel_deltaR_b" , 70 , -20 , 50 );
  sel_deltaR_b_ = fs->make<TH1D>("sel_deltaR_b_" , "sel_deltaR_b_" , 70 , -20 , 50 );
  sel_deltaR_q = fs->make<TH1D>("sel_deltaR_q" , "sel_deltaR_q" , 70 , -20 , 50 );
  sel_deltaR_q_ = fs->make<TH1D>("sel_deltaR_q_" , "sel_deltaR_q_" , 70 , -20 , 50 );
  
  deltaR_b = fs->make<TH1D>("deltaR_b" , "deltaR_b" , 70 , -20 , 50 );
  deltaR_b_ = fs->make<TH1D>("deltaR_b_" , "deltaR_b_" , 70 , -20 , 50 );
  deltaR_q = fs->make<TH1D>("deltaR_q" , "deltaR_q" , 40 , -70 , 50 );
  deltaR_q_ = fs->make<TH1D>("deltaR_q_" , "deltaR_q_" , 70 , -20 , 50 );
  
  discriminantJPSelection = fs->make<TH1D>("discriminantJPSelection" , "discriminantJPSelection" , 2000, -1000, 1000 );
  discriminantJPMatching = fs->make<TH1D>("discriminantJPMatching" , "discriminantJPMatching" , 2000, -1000, 1000 );
  discriminantJPvsN = fs->make<TH2D>("discriminantJPvsN" , "discriminantJPvsN" , 110 , -5 , 105, 2000, -1000, 1000 );
  
  discriminantJPSelection_Data = fs->make<TH1D>("discriminantJPSelection_Data" , "discriminantJPSelection_Data" , 2000, -1000, 1000 );
  discriminantJPvsN_Data = fs->make<TH2D>("discriminantJPvsN_Data" , "discriminantJPvsN_Data" , 110 , -5 , 105, 2000, -1000, 1000 );
  0404*/
  /*
    ElDiscriminante1_Data0 = fs->make<TH1D>("ElDiscriminante1_Data0" , "ElDiscriminante1_Data0" , 1100 , -5 , 105);
    ElDiscriminante2_Data0 = fs->make<TH1D>("ElDiscriminante2_Data0" , "ElDiscriminante2_Data0" , 1100 , -5 , 105);

    ElDiscriminante1_Data1 = fs->make<TH1D>("ElDiscriminante1_Data1" , "ElDiscriminante1_Data1" , 1100 , -5 , 105);
    ElDiscriminante2_Data1 = fs->make<TH1D>("ElDiscriminante2_Data1" , "ElDiscriminante2_Data1" , 1100 , -5 , 105);
  */
/*0404
  ElDiscriminante1 = fs->make<TH1D>("ElDiscriminante1" , "ElDiscriminante1" , 1100 , -5 , 105);
  ElDiscriminante2 = fs->make<TH1D>("ElDiscriminante2" , "ElDiscriminante2" , 1100 , -5 , 105);
  deltaPhiJetLep_vs_deltaPhiJetMET_bLep= fs->make<TH2D>("deltaPhiJetLep_vs_deltaPhiJetMET_bLep" , "deltaPhiJetLep_vs_deltaPhiJetMET_bLep" , 300 , 0 , 3, 300, 0, 3 );
  deltaPhiJetQ_vs_deltaPhiJetQ__bHad= fs->make<TH2D>("deltaPhiJetQ_vs_deltaPhiJetQ__bHad" , "deltaPhiJetQ_vs_deltaPhiJetQ__bHad" , 300 , 0 , 3, 300, 0, 3 );
  deltaPhiJetLepbLep= fs->make<TH1D>("deltaPhiJetLep_bLep" , "deltaPhiJetLep_bLep" , 30 , 0 , 10);
  deltaPhiJetMETbLep= fs->make<TH1D>("deltaPhiJetMET_bLep" , "deltaPhiJetMET_bLep" , 30 , 0 , 10);
  deltaPhiJetLepMETbLep= fs->make<TH1D>("deltaPhiJetLepMET_bLep" , "deltaPhiJetLepMET_bLep" , 30 , 0 , 10);

  deltaPhiJetQbHad= fs->make<TH1D>("deltaPhiJetQ_bHad" , "deltaPhiJetQ_bHad" , 30 , 0 , 10);
  deltaPhiJetQ_bHad= fs->make<TH1D>("deltaPhiJetQ__bHad" , "deltaPhiJetQ__bHad" , 30 , 0 , 10);
  deltaPhiJetQQ_bHad= fs->make<TH1D>("deltaPhiJetQQ__bHad" , "deltaPhiJetQQ__bHad" , 30 , 0 , 10);
0404*/

  /*
    Discriminant_had_W_10_10 = fs->make<TH1D>("Discriminant_had_W_10_10" , "Discriminant_had_W_10_10" , 60 , -1 , 2);
    Discriminant_lep_W_10_10 = fs->make<TH1D>("Discriminant_lep_W_10_10" , "Discriminant_lep_W_10_10" , 60 , -1 , 2);
    Discriminant_total_W_10_10 = fs->make<TH1D>("Discriminant_total_W_10_10" , "Discriminant_total_W_10_10" , 60 , -1 , 2);
    Discriminant_had_t_10_10 = fs->make<TH1D>("Discriminant_had_t_10_10" , "Discriminant_had_t_10_10" , 60 , -1 , 2);
    Discriminant_lep_t_10_10 = fs->make<TH1D>("Discriminant_lep_t_10_10" , "Discriminant_lep_t_10_10" , 60 , -1 , 2);
    Discriminant_total_t_10_10 = fs->make<TH1D>("Discriminant_total_t_10_10" , "Discriminant_total_t_10_10" , 60 , -1 , 2);

    Discriminant_had_W_30_30 = fs->make<TH1D>("Discriminant_had_W_30_30" , "Discriminant_had_W_30_30" , 60 , -1 , 2);
    Discriminant_lep_W_30_30 = fs->make<TH1D>("Discriminant_lep_W_30_30" , "Discriminant_lep_W_30_30" , 60 , -1 , 2);
    Discriminant_total_W_30_30 = fs->make<TH1D>("Discriminant_total_W_30_30" , "Discriminant_total_W_30_30" , 60 , -1 , 2);
    Discriminant_had_t_30_30 = fs->make<TH1D>("Discriminant_had_t_30_30" , "Discriminant_had_t_30_30" , 60 , -1 , 2);
    Discriminant_lep_t_30_30 = fs->make<TH1D>("Discriminant_lep_t_30_30" , "Discriminant_lep_t_30_30" , 60 , -1 , 2);
    Discriminant_total_t_30_30 = fs->make<TH1D>("Discriminant_total_t_30_30" , "Discriminant_total_t_30_30" , 60 , -1 , 2);

    Discriminant_had_W_50_50 = fs->make<TH1D>("Discriminant_had_W_50_50" , "Discriminant_had_W_50_50" , 60 , -1 , 2);
    Discriminant_lep_W_50_50 = fs->make<TH1D>("Discriminant_lep_W_50_50" , "Discriminant_lep_W_50_50" , 60 , -1 , 2);
    Discriminant_total_W_50_50 = fs->make<TH1D>("Discriminant_total_W_50_50" , "Discriminant_total_W_50_50" , 60 , -1 , 2);
    Discriminant_had_t_50_50 = fs->make<TH1D>("Discriminant_had_t_50_50" , "Discriminant_had_t_50_50" , 60 , -1 , 2);
    Discriminant_lep_t_50_50 = fs->make<TH1D>("Discriminant_lep_t_50_50" , "Discriminant_lep_t_50_50" , 60 , -1 , 2);
    Discriminant_total_t_50_50 = fs->make<TH1D>("Discriminant_total_t_50_50" , "Discriminant_total_t_50_50" , 60 , -1 , 2);

    Discriminant_had_W_70_70 = fs->make<TH1D>("Discriminant_had_W_70_70" , "Discriminant_had_W_70_70" , 60 , -1 , 2);
    Discriminant_lep_W_70_70 = fs->make<TH1D>("Discriminant_lep_W_70_70" , "Discriminant_lep_W_70_70" , 60 , -1 , 2);
    Discriminant_total_W_70_70 = fs->make<TH1D>("Discriminant_total_W_70_70" , "Discriminant_total_W_70_70" , 60 , -1 , 2);
    Discriminant_had_t_70_70 = fs->make<TH1D>("Discriminant_had_t_70_70" , "Discriminant_had_t_70_70" , 60 , -1 , 2);
    Discriminant_lep_t_70_70 = fs->make<TH1D>("Discriminant_lep_t_70_70" , "Discriminant_lep_t_70_70" , 60 , -1 , 2);
    Discriminant_total_t_70_70 = fs->make<TH1D>("Discriminant_total_t_70_70" , "Discriminant_total_t_70_70" , 60 , -1 , 2);

    Discriminant_had_W_90_90 = fs->make<TH1D>("Discriminant_had_W_90_90" , "Discriminant_had_W_90_90" , 60 , -1 , 2);
    Discriminant_lep_W_90_90 = fs->make<TH1D>("Discriminant_lep_W_90_90" , "Discriminant_lep_W_90_90" , 60 , -1 , 2);
    Discriminant_total_W_90_90 = fs->make<TH1D>("Discriminant_total_W_90_90" , "Discriminant_total_W_90_90" , 60 , -1 , 2);
    Discriminant_had_t_90_90 = fs->make<TH1D>("Discriminant_had_t_90_90" , "Discriminant_had_t_90_90" , 60 , -1 , 2);
    Discriminant_lep_t_90_90 = fs->make<TH1D>("Discriminant_lep_t_90_90" , "Discriminant_lep_t_90_90" , 60 , -1 , 2);
    Discriminant_total_t_90_90 = fs->make<TH1D>("Discriminant_total_t_90_90" , "Discriminant_total_t_90_90" , 60 , -1 , 2);

    Discriminant_had_W_10_90 = fs->make<TH1D>("Discriminant_had_W_10_90" , "Discriminant_had_W_10_90" , 60 , -1 , 2);
    Discriminant_lep_W_10_90 = fs->make<TH1D>("Discriminant_lep_W_10_90" , "Discriminant_lep_W_10_90" , 60 , -1 , 2);
    Discriminant_total_W_10_90 = fs->make<TH1D>("Discriminant_total_W_10_90" , "Discriminant_total_W_10_90" , 60 , -1 , 2);
    Discriminant_had_t_10_90 = fs->make<TH1D>("Discriminant_had_t_10_90" , "Discriminant_had_t_10_90" , 60 , -1 , 2);
    Discriminant_lep_t_10_90 = fs->make<TH1D>("Discriminant_lep_t_10_90" , "Discriminant_lep_t_10_90" , 60 , -1 , 2);
    Discriminant_total_t_10_90 = fs->make<TH1D>("Discriminant_total_t_10_90" , "Discriminant_total_t_10_90" , 60 , -1 , 2);

    Discriminant_had_W_30_70 = fs->make<TH1D>("Discriminant_had_W_30_70" , "Discriminant_had_W_30_70" , 60 , -1 , 2);
    Discriminant_lep_W_30_70 = fs->make<TH1D>("Discriminant_lep_W_30_70" , "Discriminant_lep_W_30_70" , 60 , -1 , 2);
    Discriminant_total_W_30_70 = fs->make<TH1D>("Discriminant_total_W_30_70" , "Discriminant_total_W_30_70" , 60 , -1 , 2);
    Discriminant_had_t_30_70 = fs->make<TH1D>("Discriminant_had_t_30_70" , "Discriminant_had_t_30_70" , 60 , -1 , 2);
    Discriminant_lep_t_30_70 = fs->make<TH1D>("Discriminant_lep_t_30_70" , "Discriminant_lep_t_30_70" , 60 , -1 , 2);
    Discriminant_total_t_30_70 = fs->make<TH1D>("Discriminant_total_t_30_70" , "Discriminant_total_t_30_70" , 60 , -1 , 2);

    Discriminant_had_W_70_30 = fs->make<TH1D>("Discriminant_had_W_70_30" , "Discriminant_had_W_70_30" , 60 , -1 , 2);
    Discriminant_lep_W_70_30 = fs->make<TH1D>("Discriminant_lep_W_70_30" , "Discriminant_lep_W_70_30" , 60 , -1 , 2);
    Discriminant_total_W_70_30 = fs->make<TH1D>("Discriminant_total_W_70_30" , "Discriminant_total_W_70_30" , 60 , -1 , 2);
    Discriminant_had_t_70_30 = fs->make<TH1D>("Discriminant_had_t_70_30" , "Discriminant_had_t_70_30" , 60 , -1 , 2);
    Discriminant_lep_t_70_30 = fs->make<TH1D>("Discriminant_lep_t_70_30" , "Discriminant_lep_t_70_30" , 60 , -1 , 2);
    Discriminant_total_t_70_30 = fs->make<TH1D>("Discriminant_total_t_70_30" , "Discriminant_total_t_70_30" , 60 , -1 , 2);

    Discriminant_had_W_90_10 = fs->make<TH1D>("Discriminant_had_W_90_10" , "Discriminant_had_W_90_10" , 60 , -1 , 2);
    Discriminant_lep_W_90_10 = fs->make<TH1D>("Discriminant_lep_W_90_10" , "Discriminant_lep_W_90_10" , 60 , -1 , 2);
    Discriminant_total_W_90_10 = fs->make<TH1D>("Discriminant_total_W_90_10" , "Discriminant_total_W_90_10" , 60 , -1 , 2);
    Discriminant_had_t_90_10 = fs->make<TH1D>("Discriminant_had_t_90_10" , "Discriminant_had_t_90_10" , 60 , -1 , 2);
    Discriminant_lep_t_90_10 = fs->make<TH1D>("Discriminant_lep_t_90_10" , "Discriminant_lep_t_90_10" , 60 , -1 , 2);
    Discriminant_total_t_90_10 = fs->make<TH1D>("Discriminant_total_t_90_10" , "Discriminant_total_t_90_10" , 60 , -1 , 2);
  */

/*0404
  Discriminant_had_W_90_90 = fs->make<TH1D>("Discriminant_had_W_90_90" , "Discriminant_had_W_90_90" , 60 , -1 , 2);
  Discriminant_lep_W_90_90 = fs->make<TH1D>("Discriminant_lep_W_90_90" , "Discriminant_lep_W_90_90" , 60 , -1 , 2);
  Discriminant_total_W_90_90 = fs->make<TH1D>("Discriminant_total_W_90_90" , "Discriminant_total_W_90_90" , 60 , -1 , 2);
  Discriminant_had_t_90_90 = fs->make<TH1D>("Discriminant_had_t_90_90" , "Discriminant_had_t_90_90" , 60 , -1 , 2);
  Discriminant_lep_t_90_90 = fs->make<TH1D>("Discriminant_lep_t_90_90" , "Discriminant_lep_t_90_90" , 60 , -1 , 2);
  Discriminant_total_t_90_90 = fs->make<TH1D>("Discriminant_total_t_90_90" , "Discriminant_total_t_90_90" , 60 , -1 , 2);


  Discriminant_had_W_85_85 = fs->make<TH1D>("Discriminant_had_W_85_85" , "Discriminant_had_W_85_85" , 60 , -1 , 2);
  Discriminant_lep_W_85_85 = fs->make<TH1D>("Discriminant_lep_W_85_85" , "Discriminant_lep_W_85_85" , 60 , -1 , 2);
  Discriminant_total_W_85_85 = fs->make<TH1D>("Discriminant_total_W_85_85" , "Discriminant_total_W_85_85" , 60 , -1 , 2);
  Discriminant_had_t_85_85 = fs->make<TH1D>("Discriminant_had_t_85_85" , "Discriminant_had_t_85_85" , 60 , -1 , 2);
  Discriminant_lep_t_85_85 = fs->make<TH1D>("Discriminant_lep_t_85_85" , "Discriminant_lep_t_85_85" , 60 , -1 , 2);
  Discriminant_total_t_85_85 = fs->make<TH1D>("Discriminant_total_t_85_85" , "Discriminant_total_t_85_85" , 60 , -1 , 2);


  Discriminant_had_W_95_95 = fs->make<TH1D>("Discriminant_had_W_95_95" , "Discriminant_had_W_95_95" , 60 , -1 , 2);
  Discriminant_lep_W_95_95 = fs->make<TH1D>("Discriminant_lep_W_95_95" , "Discriminant_lep_W_95_95" , 60 , -1 , 2);
  Discriminant_total_W_95_95 = fs->make<TH1D>("Discriminant_total_W_95_95" , "Discriminant_total_W_95_95" , 60 , -1 , 2);
  Discriminant_had_t_95_95 = fs->make<TH1D>("Discriminant_had_t_95_95" , "Discriminant_had_t_95_95" , 60 , -1 , 2);
  Discriminant_lep_t_95_95 = fs->make<TH1D>("Discriminant_lep_t_95_95" , "Discriminant_lep_t_95_95" , 60 , -1 , 2);
  Discriminant_total_t_95_95 = fs->make<TH1D>("Discriminant_total_t_95_95" , "Discriminant_total_t_95_95" , 60 , -1 , 2);


  Discriminant_had_W_95_85 = fs->make<TH1D>("Discriminant_had_W_95_85" , "Discriminant_had_W_95_85" , 60 , -1 , 2);
  Discriminant_lep_W_95_85 = fs->make<TH1D>("Discriminant_lep_W_95_85" , "Discriminant_lep_W_95_85" , 60 , -1 , 2);
  Discriminant_total_W_95_85 = fs->make<TH1D>("Discriminant_total_W_95_85" , "Discriminant_total_W_95_85" , 60 , -1 , 2);
  Discriminant_had_t_95_85 = fs->make<TH1D>("Discriminant_had_t_95_85" , "Discriminant_had_t_95_85" , 60 , -1 , 2);
  Discriminant_lep_t_95_85 = fs->make<TH1D>("Discriminant_lep_t_95_85" , "Discriminant_lep_t_95_85" , 60 , -1 , 2);
  Discriminant_total_t_95_85 = fs->make<TH1D>("Discriminant_total_t_95_85" , "Discriminant_total_t_95_85" , 60 , -1 , 2);


  Discriminant_had_W_85_95 = fs->make<TH1D>("Discriminant_had_W_85_95" , "Discriminant_had_W_85_95" , 60 , -1 , 2);
  Discriminant_lep_W_85_95 = fs->make<TH1D>("Discriminant_lep_W_85_95" , "Discriminant_lep_W_85_95" , 60 , -1 , 2);
  Discriminant_total_W_85_95 = fs->make<TH1D>("Discriminant_total_W_85_95" , "Discriminant_total_W_85_95" , 60 , -1 , 2);
  Discriminant_had_t_85_95 = fs->make<TH1D>("Discriminant_had_t_85_95" , "Discriminant_had_t_85_95" , 60 , -1 , 2);
  Discriminant_lep_t_85_95 = fs->make<TH1D>("Discriminant_lep_t_85_95" , "Discriminant_lep_t_85_95" , 60 , -1 , 2);
  Discriminant_total_t_85_95 = fs->make<TH1D>("Discriminant_total_t_85_95" , "Discriminant_total_t_85_95" , 60 , -1 , 2);
0404*/
}



FullAndWorking::~FullAndWorking()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
FullAndWorking::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //system("./prueba.sh"); 
  std::cout<<"JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"<<std::endl; 
  using namespace edm;
  using namespace std;
  double Total_Energy_Jets;
  ostringstream id;
  id << iEvent.id().event();   // Convert value into a string.
  string s_id = id.str();      // Get the created string from the output stream.
  /*0404
  if (top)
    {
      file<<"top"<<endl;
    }
  else
    {
      file<<"stop"<<endl;
    }
  
  
  //bool goodEvent=false;
  0404*/
  double discriminants[4];
  double genToLHCO[9][5];//[5,-5,2,-1,11,met,1000022,1000022,-12][eta,phi,pt,mass,btag] if top met=-12 else met=(-12)+(1000022)+(1000022)	
  double recoToLHCO[6][5];//[5,-5,2,-1,11,-12][eta,phi,pt,mass,btag]
  double genEnergy[9];//[5,-5,2,-1,11,met,1000022,1000022,-12]  if top met=-12 else met=(-12)+(1000022)+(1000022)
  double recoEnergy[6];//[5,-5,2,-1,11,-12]
  jets->Fill(1);
  double genMuonPhi=0;
  //file<<"Number of the Event :"<<iEvent.id()<<endl;
  //bool MC =true;
  double InvMass;
  double InvMass_;
  
  
  if (MC)
    {
      Handle<reco::GenParticleCollection> gParticles;
      iEvent.getByLabel("genParticles", gParticles);
      const reco::GenParticleCollection* gPARTICLES = gParticles.product();
      //cout<<"Number of Particles: "<<gPARTICLES->size()<<endl;
      bool t,t_,sw,sw_,w,w_,b,b_,q,q_,m_,vm_,nsw,nsw_;
      sw=false;
      sw_=false;
      w=false;
      w_=false;
      b=false;
      b_=false;
      q=false;
      q_=false;
      m_=false;
      vm_=false;
      nsw=false;
      nsw_=false;
      t=false;
      t_=false;	
      unsigned int stDaughters, st_Daughters,tDaughters, t_Daughters, SWDaughters,SW_Daughters,WDaughters,W_Daughters;
      tDaughters=0;
      t_Daughters=0;
      WDaughters=0;
      W_Daughters=0;
      SWDaughters=0;
      SW_Daughters=0;
      stDaughters=0;
      st_Daughters=0;
      bool GoodDecay = false;
      int numberOfParticlesInTheProcess=gPARTICLES->size();
      for (int j = 0; j<5; j++)
	{
	  genToLHCO[5][j]=0;
	  genToLHCO[6][j]=0;
	  genToLHCO[7][j]=0;
	}
      for(unsigned int i =0; i<gPARTICLES->size();i++)
	{
	  if (!top)
	    {
	      if (!neutralino)
		{
		  const reco::GenParticle& gcandPARTICLE = (*gPARTICLES)[i];
		  //cout<<"pdgId: "<<gcandPARTICLE.pdgId()<<endl;
		  if (gcandPARTICLE.pdgId() == 1000006)
		    {
		      InvMass = gcandPARTICLE.mass();
		      //file<<"st found"<<endl;
		      if (gcandPARTICLE.numberOfDaughters())
			{
			  stDaughters = gcandPARTICLE.numberOfDaughters();
			  for ( unsigned istDaughter=0; istDaughter<stDaughters; ++istDaughter )
			    {
			      const reco::Candidate* stdaught = gcandPARTICLE.daughter(istDaughter);
			      //file<<"st's daughter :"<<stdaught->pdgId()<<endl;
			      if ( (stdaught->pdgId()) == 1000024 )
				{
				  sw=true;
				  SWDaughters = stdaught->numberOfDaughters();
				  //file<<"SW found"<<endl;
				  for ( unsigned iSWDaughter=0; iSWDaughter<SWDaughters; ++iSWDaughter )
				    {
				      const reco::Candidate* swdaught = stdaught->daughter(iSWDaughter);
				      //file<<"SW's duaghters pdgId: "<<swdaught->pdgId()<<endl;
				      if ( (swdaught->pdgId())== 24 )
					{
					  w=true;
					  WDaughters = swdaught->numberOfDaughters();
					  //file<<"W found"<<endl;
					  for ( unsigned iWDaughter=0; iWDaughter<WDaughters; ++iWDaughter )
					    {
					      const reco::Candidate* wdaught = swdaught->daughter(iWDaughter);
					      //file<<"W's duaghters pdgId: "<<wdaught->pdgId()<<endl;
					      if ( (wdaught->pdgId())< 5 && (wdaught->pdgId()) > 0 )
						{
						  q=true;
						  //file<<"q found"<<endl;
						  genToLHCO[2][0]=wdaught->eta();
						  genToLHCO[2][1]=wdaught->phi();
						  genToLHCO[2][2]=wdaught->pt();
						  genToLHCO[2][3]=wdaught->mass();//0;
						  genToLHCO[2][4]=0;
						  genEnergy[2]=wdaught->energy();
						}
					      if ( (wdaught->pdgId()) > -5 && (wdaught->pdgId()) < 0 )
						{
						  q_=true;
						  //file<<"q_ found"<<endl;
						  genToLHCO[3][0]=wdaught->eta();
						  genToLHCO[3][1]=wdaught->phi();
						  genToLHCO[3][2]=wdaught->pt();
						  genToLHCO[3][3]=wdaught->mass();//0;
						  genToLHCO[3][4]=0;
						  genEnergy[3]=wdaught->energy();
						}
					    }
					}
				      if ( (swdaught->pdgId()) == 1000022 )
					{
					  nsw=true;
					  //file<<"nsw found"<<endl;
					  genToLHCO[6][0]=swdaught->eta();
					  genToLHCO[6][1]=swdaught->phi();
					  genToLHCO[6][2]=swdaught->pt();
					  genToLHCO[6][3]=swdaught->mass();//0;
					  genToLHCO[6][4]=0;
					  genEnergy[6]=swdaught->energy();
				      
					}
				    }
				}
			      if ( (stdaught->pdgId()) == 5 )
				{
				  b = true;
				  //file<<"b found"<<endl;
				  genToLHCO[0][0]=stdaught->eta();
				  genToLHCO[0][1]=stdaught->phi();
				  genToLHCO[0][2]=stdaught->pt();
				  genToLHCO[0][3]=stdaught->mass();//4.7;
				  genToLHCO[0][4]=2;
				  genEnergy[0]=stdaught->energy();
				}
			      /*
				if ( (stdaught->pdgId()) < 4 && (stdaught->pdgId()) > 0 )
				{
				q=true;
				file<<"q found"<<endl;
				genToLHCO[2][0]=stdaught->eta();
				genToLHCO[2][1]=stdaught->phi();
				genToLHCO[2][2]=stdaught->pt();
				genToLHCO[2][3]=stdaught->mass();//0;
				genToLHCO[2][4]=0;
				}
				if ( (stdaught->pdgId()) > -4 && (stdaught->pdgId()) < 0 )
				{
				q_=true;
				file<<"q_ found"<<endl;
				genToLHCO[3][0]=stdaught->eta();
				genToLHCO[3][1]=stdaught->phi();
				genToLHCO[3][2]=stdaught->pt();
				genToLHCO[3][3]=stdaught->mass();//0;
				genToLHCO[3][4]=0;
				}
			      */
			  
			    }
			}
		    }
	      
		  if (gcandPARTICLE.pdgId() == -1000006)
		    {
		      InvMass_ = gcandPARTICLE.mass();
		      //file<<"st_ found"<<endl;
		      if (gcandPARTICLE.numberOfDaughters())
			{
			  st_Daughters = gcandPARTICLE.numberOfDaughters();
			  for ( unsigned ist_Daughter=0; ist_Daughter<st_Daughters; ++ist_Daughter )
			    {
			      const reco::Candidate* st_daught = gcandPARTICLE.daughter(ist_Daughter);
			      //file<<"st_'s daughter :"<<st_daught->pdgId()<<endl;
			      if ( (st_daught->pdgId()) == -1000024 )
				{
				  sw_=true;
				  SW_Daughters = st_daught->numberOfDaughters();
				  //file<<"SW_ found"<<endl;
				  for ( unsigned iSW_Daughter=0; iSW_Daughter<SW_Daughters; ++iSW_Daughter )
				    {
				      const reco::Candidate* sw_daught = st_daught->daughter(iSW_Daughter);
				      //file<<"SW's duaghters pdgId: "<<sw_daught->pdgId()<<endl;
				      if ( (sw_daught->pdgId())== -24 )
					{
					  w_=true;
					  W_Daughters = sw_daught->numberOfDaughters();
					  //file<<"W_ found"<<endl;
					  for ( unsigned iW_Daughter=0; iW_Daughter<W_Daughters; ++iW_Daughter )
					    {
					      const reco::Candidate* w_daught = sw_daught->daughter(iW_Daughter);
					      //file<<"W_'s duaghters pdgId: "<<w_daught->pdgId()<<endl;
					      if ( (w_daught->pdgId())== 13 )
						{
						  m_=true;
						  //file<<"m found"<<endl;
						  genMuonPhi=w_daught->phi();
						  genToLHCO[4][0]=w_daught->eta();
						  genToLHCO[4][1]=w_daught->phi();
						  genToLHCO[4][2]=w_daught->pt();
						  genToLHCO[4][3]=w_daught->mass();//0;
						  genToLHCO[4][4]=0;
						  genEnergy[4]=w_daught->energy();
						}
					      if ( (w_daught->pdgId())== -14 )
						{
						  vm_=true;
						  //file<<"vm found"<<endl;
						  genToLHCO[8][0]=w_daught->eta();
						  genToLHCO[8][1]=w_daught->phi();
						  genToLHCO[8][2]=w_daught->pt();
						  genToLHCO[8][3]=w_daught->mass();//0;
						  genToLHCO[8][4]=0;
						  genEnergy[8]=w_daught->energy();
						} 
					  
					  
					    }
					}
				      if ( (sw_daught->pdgId()) == 1000022 )
					{
					  nsw_=true;
					  //file<<"nsw_ found"<<endl;
					  genToLHCO[7][0]=sw_daught->eta();
					  genToLHCO[7][1]=sw_daught->phi();
					  genToLHCO[7][2]=sw_daught->pt();
					  genToLHCO[7][3]=sw_daught->mass();//0;
					  genToLHCO[7][4]=0;
					  genEnergy[7]=sw_daught->energy();
					}
				    }
				}
			      if ( (st_daught->pdgId()) == -5 )
				{
				  b_ = true;
				  //file<<"b found"<<endl;
				  genToLHCO[1][0]=st_daught->eta();
				  genToLHCO[1][1]=st_daught->phi();
				  genToLHCO[1][2]=st_daught->pt();
				  genToLHCO[1][3]=st_daught->mass();//4.7;
				  genToLHCO[1][4]=2;
				  genEnergy[1]=st_daught->energy();
				}
			  
			      /*
				if ( (st_daught->pdgId())== 13 )
				{
				m_=true;
				file<<"e found"<<endl;
				genMuonPhi=st_daught->phi();
				genToLHCO[4][0]=st_daught->eta();
				genToLHCO[4][1]=st_daught->phi();
				genToLHCO[4][2]=st_daught->pt();
				genToLHCO[4][3]=st_daught->mass();//0;
				genToLHCO[4][4]=0;
				genEnergy[4]=st_daught->energy();
				}
				if ( (st_daught->pdgId())== -14 )
				{
				vm_=true;
				file<<"ve found"<<endl;
				genToLHCO[8][0]=st_daught->eta();
				genToLHCO[8][1]=st_daught->phi();
				genToLHCO[8][2]=st_daught->pt();
				genToLHCO[8][3]=st_daught->mass();//0;
				genToLHCO[8][4]=0;
				genEnergy[8]=st_daught->energy();
				}
			      */
			    }
			}
		    }
                }
	      if (neutralino)
                {
		  const reco::GenParticle& gcandPARTICLE = (*gPARTICLES)[i];
		  //cout<<"pdgId: "<<gcandPARTICLE.pdgId()<<endl;
		  if (gcandPARTICLE.pdgId() == 1000006)
		    {
		      InvMass = gcandPARTICLE.mass();
		      //file<<"st found"<<endl;
		      if (gcandPARTICLE.numberOfDaughters())
			{
			  stDaughters = gcandPARTICLE.numberOfDaughters();
			  for ( unsigned istDaughter=0; istDaughter<stDaughters; ++istDaughter )
			    {
			      const reco::Candidate* stdaught = gcandPARTICLE.daughter(istDaughter);
			      //file<<"st's daughter :"<<stdaught->pdgId()<<endl;
			      if ( (stdaught->pdgId()) == 6 )
				{
				  t=true;
				  tDaughters = stdaught->numberOfDaughters();
				  //file<<"t found"<<endl;
				  for ( unsigned itDaughter=0; itDaughter<tDaughters; ++itDaughter )
				    {
				      const reco::Candidate* tdaught = stdaught->daughter(itDaughter);
				      //file<<"t's duaghters pdgId: "<<tdaught->pdgId()<<endl;
				      if ( (tdaught->pdgId())== 24 )
					{
					  w=true;
					  WDaughters = tdaught->numberOfDaughters();
					  //file<<"W found"<<endl;
					  for ( unsigned iWDaughter=0; iWDaughter<WDaughters; ++iWDaughter )
					    {
					      const reco::Candidate* wdaught = tdaught->daughter(iWDaughter);
					      //file<<"W's duaghters pdgId: "<<wdaught->pdgId()<<endl;
					      if ( (wdaught->pdgId())< 5 && (wdaught->pdgId()) > 0 )
						{
						  q=true;
						  //file<<"q found"<<endl;
						  genToLHCO[2][0]=wdaught->eta();
						  genToLHCO[2][1]=wdaught->phi();
						  genToLHCO[2][2]=wdaught->pt();
						  genToLHCO[2][3]=wdaught->mass();//0;
						  genToLHCO[2][4]=0;
						  genEnergy[2]=wdaught->energy();
						}
					      if ( (wdaught->pdgId()) > -5 && (wdaught->pdgId()) < 0 )
						{
						  q_=true;
						  //file<<"q_ found"<<endl;
						  genToLHCO[3][0]=wdaught->eta();
						  genToLHCO[3][1]=wdaught->phi();
						  genToLHCO[3][2]=wdaught->pt();
						  genToLHCO[3][3]=wdaught->mass();//0;
						  genToLHCO[3][4]=0;
						  genEnergy[3]=wdaught->energy();
						}
					    }
					}
				      if ( (tdaught->pdgId()) == 5 )
					{
					  b=true;
					  //file<<"b found"<<endl;
					  genToLHCO[0][0]=tdaught->eta();
					  genToLHCO[0][1]=tdaught->phi();
					  genToLHCO[0][2]=tdaught->pt();
					  genToLHCO[0][3]=tdaught->mass();//0;
					  genToLHCO[0][4]=0;
					  genEnergy[0]=tdaught->energy();
				      
					}
				    }
				}
			      if ( (stdaught->pdgId()) == 1000022 )
				{
				  nsw = true;
				  //file<<"nsw found"<<endl;
				  genToLHCO[6][0]=stdaught->eta();
				  genToLHCO[6][1]=stdaught->phi();
				  genToLHCO[6][2]=stdaught->pt();
				  genToLHCO[6][3]=stdaught->mass();//4.7;
				  genToLHCO[6][4]=2;
				  genEnergy[6]=stdaught->energy();
				}
			      /*
				if ( (stdaught->pdgId()) < 4 && (stdaught->pdgId()) > 0 )
				{
				q=true;
				file<<"q found"<<endl;
				genToLHCO[2][0]=stdaught->eta();
				genToLHCO[2][1]=stdaught->phi();
				genToLHCO[2][2]=stdaught->pt();
				genToLHCO[2][3]=0;
				genToLHCO[2][4]=0;
				}
				if ( (stdaught->pdgId()) > -4 && (stdaught->pdgId()) < 0 )
				{
				q_=true;
				file<<"q_ found"<<endl;
				genToLHCO[3][0]=stdaught->eta();
				genToLHCO[3][1]=stdaught->phi();
				genToLHCO[3][2]=stdaught->pt();
				genToLHCO[3][3]=0;
				genToLHCO[3][4]=0;
				}
			      */
			  
			    }
			}
		    }
	      
		  if (gcandPARTICLE.pdgId() == -1000006)
		    {
		      InvMass_ = gcandPARTICLE.mass();
		      //file<<"st_ found"<<endl;
		      if (gcandPARTICLE.numberOfDaughters())
			{
			  st_Daughters = gcandPARTICLE.numberOfDaughters();
			  for ( unsigned ist_Daughter=0; ist_Daughter<st_Daughters; ++ist_Daughter )
			    {
			      const reco::Candidate* st_daught = gcandPARTICLE.daughter(ist_Daughter);
			      //file<<"st_'s daughter :"<<st_daught->pdgId()<<endl;
			      if ( (st_daught->pdgId()) == -6 )
				{
				  t_=true;
				  t_Daughters = st_daught->numberOfDaughters();
				  //file<<"t_ found"<<endl;
				  for ( unsigned it_Daughter=0; it_Daughter<t_Daughters; ++it_Daughter )
				    {
				      const reco::Candidate* t_daught = st_daught->daughter(it_Daughter);
				      //file<<"t_'s duaghters pdgId: "<<t_daught->pdgId()<<endl;
				      if ( (t_daught->pdgId())== -24 )
					{
					  w_=true;
					  W_Daughters = t_daught->numberOfDaughters();
					  //file<<"W_ found"<<endl;
					  for ( unsigned iW_Daughter=0; iW_Daughter<W_Daughters; ++iW_Daughter )
					    {
					      const reco::Candidate* w_daught = t_daught->daughter(iW_Daughter);
					      //file<<"W_'s duaghters pdgId: "<<w_daught->pdgId()<<endl;
					      if ( (w_daught->pdgId())== 13 )
						{
						  m_=true;
						  //file<<"m found"<<endl;
						  genMuonPhi=w_daught->phi();
						  genToLHCO[4][0]=w_daught->eta();
						  genToLHCO[4][1]=w_daught->phi();
						  genToLHCO[4][2]=w_daught->pt();
						  genToLHCO[4][3]=w_daught->mass();//0;
						  genToLHCO[4][4]=0;
						  genEnergy[4]=w_daught->energy();
						}
					      if ( (w_daught->pdgId())== -14 )
						{
						  vm_=true;
						  //file<<"vm found"<<endl;
						  genToLHCO[8][0]=w_daught->eta();
						  genToLHCO[8][1]=w_daught->phi();
						  genToLHCO[8][2]=w_daught->pt();
						  genToLHCO[8][3]=w_daught->mass();//0;
						  genToLHCO[8][4]=0;
						  genEnergy[8]=w_daught->energy();
						} 
					  
					  
					    }
					}
				      if ( (t_daught->pdgId()) == -5 )
					{
					  b_=true;
					  //file<<"nsw_ found"<<endl;
					  genToLHCO[1][0]=t_daught->eta();
					  genToLHCO[1][1]=t_daught->phi();
					  genToLHCO[1][2]=t_daught->pt();
					  genToLHCO[1][3]=t_daught->mass();//0;
					  genToLHCO[1][4]=0;
					  genEnergy[1]=t_daught->energy();
					}
				    }
				}
			      if ( (st_daught->pdgId()) == 1000022 )
				{
				  nsw_ = true;
				  //file<<"b found"<<endl;
				  genToLHCO[7][0]=st_daught->eta();
				  genToLHCO[7][1]=st_daught->phi();
				  genToLHCO[7][2]=st_daught->pt();
				  genToLHCO[7][3]=st_daught->mass();//4.7;
				  genToLHCO[7][4]=2;
				  genEnergy[7]=st_daught->energy();
				}
			  
			      /*
				if ( (st_daught->pdgId())== 13 )
				{
				m_=true;
				file<<"e found"<<endl;
				genMuonPhi=st_daught->phi();
				genToLHCO[4][0]=st_daught->eta();
				genToLHCO[4][1]=st_daught->phi();
				genToLHCO[4][2]=st_daught->pt();
				genToLHCO[4][3]=0;
				genToLHCO[4][4]=0;
				genEnergy[4]=st_daught->energy();
				}
				if ( (st_daught->pdgId())== -14 )
				{
				vm_=true;
				file<<"ve found"<<endl;
				genToLHCO[8][0]=st_daught->eta();
				genToLHCO[8][1]=st_daught->phi();
				genToLHCO[8][2]=st_daught->pt();
				genToLHCO[8][3]=0;
				genToLHCO[8][4]=0;
				genEnergy[8]=st_daught->energy();
				}
			      */
			    }
			}
		    }
                }		
	    }
	  if (top)
	    {
	      const reco::GenParticle& gcandPARTICLE = (*gPARTICLES)[i];
	      //cout<<"pdgId: "<<gcandPARTICLE.pdgId()<<endl;
	      if (gcandPARTICLE.pdgId() == 6)
		{
		  InvMass = gcandPARTICLE.mass();
		  //file<<"t found"<<endl;
		  if (gcandPARTICLE.numberOfDaughters())
		    {
		      tDaughters = gcandPARTICLE.numberOfDaughters();
		      for ( unsigned iDaughter=0; iDaughter<tDaughters; ++iDaughter )
			{
			  const reco::Candidate* daught = gcandPARTICLE.daughter(iDaughter);
			  //file<<"t's dughter :"<<daught->pdgId()<<endl;
			  if ( (daught->pdgId()) == 24 )
			    {
			      w=true;
			      WDaughters = daught->numberOfDaughters();
			      //file<<"W found"<<endl;
			      for ( unsigned iWDaughter=0; iWDaughter<WDaughters; ++iWDaughter )
				{
				  const reco::Candidate* decay = daught->daughter(iWDaughter);
				  int decayId = decay->pdgId();
				  //file<<"W_'s duaghters pdgId: "<<decayId<<endl;
				  if ( decayId < 5 && decayId > 0 )
				    {
				      q=true;
				      //file<<"q found"<<endl;
				      genToLHCO[2][0]=decay->eta();
				      genToLHCO[2][1]=decay->phi();
				      genToLHCO[2][2]=decay->pt();
				      genToLHCO[2][3]=decay->mass();//0;
				      genToLHCO[2][4]=0;
				      genEnergy[2]=decay->energy();
				    }
				  if ( decayId > -5 && decayId < 0 )
				    {
				      q_=true;
				      //file<<"q_ found"<<endl;
				      genToLHCO[3][0]=decay->eta();
				      genToLHCO[3][1]=decay->phi();
				      genToLHCO[3][2]=decay->pt();
				      genToLHCO[3][3]=decay->mass();//0;
				      genToLHCO[3][4]=0;
				      genEnergy[3]=decay->energy();
				    }
				}
			    }
			  if ( (daught->pdgId()) == 5 )
			    {
			      b = true;
			      //file<<"b found"<<endl;
			      genToLHCO[0][0]=daught->eta();
			      genToLHCO[0][1]=daught->phi();
			      genToLHCO[0][2]=daught->pt();
			      genToLHCO[0][3]=daught->mass();//4.7;
			      genToLHCO[0][4]=2;
			      genEnergy[0]=daught->energy();
			    }
			  /*
			    if ( (daught->pdgId()) < 4 && (daught->pdgId()) > 0 )
			    {
			    q=true;
			    file<<"q found"<<endl;
			    genToLHCO[2][0]=daught->eta();
			    genToLHCO[2][1]=daught->phi();
			    genToLHCO[2][2]=daught->pt();
			    genToLHCO[2][3]=daught->mass();//0;
			    genToLHCO[2][4]=0;
			    }
			    if ( (daught->pdgId()) > -4 && (daught->pdgId()) < 0 )
			    {
			    q_=true;
			    file<<"q_ found"<<endl;
			    genToLHCO[3][0]=daught->eta();
			    genToLHCO[3][1]=daught->phi();
			    genToLHCO[3][2]=daught->pt();
			    genToLHCO[3][3]=daught->mass();//0;
			    genToLHCO[3][4]=0;
			    }
			  
			  */
			}
		    }
		}
	      
	      if (gcandPARTICLE.pdgId() == -6)
		{
		  InvMass_ = gcandPARTICLE.mass();
		  //file<<"t_ found"<<endl;
		  if (gcandPARTICLE.numberOfDaughters())
		    {
		      t_Daughters = gcandPARTICLE.numberOfDaughters();
		      for ( unsigned iDaughter=0; iDaughter<t_Daughters; ++iDaughter )
			{
			  const reco::Candidate* daught = gcandPARTICLE.daughter(iDaughter);
			  //file<<"t_'s dughter :"<<daught->pdgId()<<endl;
			  if ( (daught->pdgId()) == -24 )
			    {
			      w_=true;
			      //file<<"W_ found"<<endl;
			      W_Daughters = daught->numberOfDaughters();
			      for ( unsigned iW_Daughter=0; iW_Daughter<W_Daughters; ++iW_Daughter )
				{
				  const reco::Candidate* W_decay = daught->daughter(iW_Daughter);
				  int decayId = W_decay->pdgId();
				  //file<<"W_'s duaghters pdgId: "<<decayId<<endl;
				  if ( decayId == 13 )
				    {
				      m_=true;
				      //file<<"m found"<<endl;
				      genMuonPhi=W_decay->phi();
				      genToLHCO[4][0]=W_decay->eta();
				      genToLHCO[4][1]=W_decay->phi();
				      genToLHCO[4][2]=W_decay->pt();
				      genToLHCO[4][3]=W_decay->mass();//0;
				      genToLHCO[4][4]=0;
				      genEnergy[4]=W_decay->energy();
				    }
				  if ( decayId == -14 )
				    {
				      vm_=true;
				      //file<<"vm found"<<endl;
				      genToLHCO[5][0]=W_decay->eta();
				      genToLHCO[5][1]=W_decay->phi();
				      genToLHCO[5][2]=W_decay->pt();
				      genToLHCO[5][3]=W_decay->mass();//0;
				      genToLHCO[5][4]=0;
				      genEnergy[5]=W_decay->energy();
				    }
				  
				}
			    }
			  if ( (daught->pdgId()) == -5 )
			    {
			      b_ = true;
			      //file<<"b_ found"<<endl;
			      genToLHCO[1][0]=daught->eta();
			      genToLHCO[1][1]=daught->phi();
			      genToLHCO[1][2]=daught->pt();
			      genToLHCO[1][3]=daught->mass();//4.7;
			      genToLHCO[1][4]=2; 
			      genEnergy[1]=daught->energy();
			    }
                          /*
			    if ( (daught->pdgId()) == 13 )
			    {
			    m_=true;
			    file<<"m found"<<endl;
			    genMuonPhi=daught->phi();
			    genToLHCO[4][0]=daught->eta();
			    genToLHCO[4][1]=daught->phi();
			    genToLHCO[4][2]=daught->pt();
			    genToLHCO[4][3]=daught->mass();//0;
			    genToLHCO[4][4]=0;
			    genEnergy[4]=daught->energy();
			    }
			    if ( (daught->pdgId()) == -14 )
			    {
			    vm_=true;
			    file<<"vm found"<<endl;
			    genToLHCO[5][0]=daught->eta();
			    genToLHCO[5][1]=daught->phi();
			    genToLHCO[5][2]=daught->pt();
			    genToLHCO[5][3]=daught->mass();//0;
			    genToLHCO[5][4]=0;
			    genEnergy[5]=daught->energy();
			    }
			  */
			}
		    }
		}
	    }
	}
      
      if (b && q && q_ && b_ && m_ && vm_ && (top || (nsw && nsw_)) )// && tDaughters == 2 && t_Daughters == 2 && WDaughters == 3 && W_Daughters == 3 )//quasi-inclusive process                
	{
	 /*0404 
         if (!(tDaughters == 2 && t_Daughters == 2 && WDaughters == 3 && W_Daughters == 3))
	    {
	      file<<"JP:The decay has another particles"<<endl;
	    }
	  if (!w)
	    {
	      file<<"JP:The decay has not a W"<<endl;
	    }
	  if (!w_)
	    {
	      file<<"JP:The decay has not a W_"<<endl;
	    }
	  
	  
	  
	  if (!(stDaughters == 2 && st_Daughters == 2 && SWDaughters == 3 && SW_Daughters == 3)) //this part is not updated so if you want to use it please update it
	    {
	      file<<"JP:The decay has another particles"<<endl;
	    }
	  if (!w)
	    {
	      file<<"JP:The decay has not a W"<<endl;
	    }
	  if (!w_)
	    {
	      file<<"JP:The decay has not a W_"<<endl;
	    }
	  
	  
	  file<<"GOOOOODJJJPPP"<<endl;
	  file1<<"Number of the Event :"<<iEvent.id()<<endl;
	  file1<<"Number of particles :"<<numberOfParticlesInTheProcess<<endl;      
	  for(unsigned int i =0; i<gPARTICLES->size();i++)
	    {
	      const reco::GenParticle& gcandPARTICLE = (*gPARTICLES)[i];
	      file1<<"pdgId: "<<gcandPARTICLE.pdgId()<<endl;
	    }
	  0404*/
	  
	  jets->Fill(3);
	  GoodDecay=true;
	  
	  if (!top)//adding nsw, nsw_ and ve
	    {
	      math::PtEtaPhiMLorentzVectorD neutrino (genToLHCO[8][2],genToLHCO[8][0],genToLHCO[8][1],genToLHCO[8][3]);
	      math::PtEtaPhiMLorentzVectorD neutralino (genToLHCO[6][2],genToLHCO[6][0],genToLHCO[6][1],genToLHCO[6][3]);
	      math::PtEtaPhiMLorentzVectorD neutralino_ (genToLHCO[7][2],genToLHCO[7][0],genToLHCO[7][1],genToLHCO[7][3]);
	      math::PtEtaPhiMLorentzVectorD met = neutrino + neutralino + neutralino_;
	      genToLHCO[5][0]=met.eta();
	      genToLHCO[5][1]=met.phi();
	      genToLHCO[5][2]=met.pt();
	      genToLHCO[5][3]=0;
	      genToLHCO[5][4]=0;
	      genEnergy[5]=met.energy();
	      
	    }
	  
	  
	  
	  //////////////////////////////////////////////////////////////////////////////////
	  //Get reco information
	  
	  
	  // para mover cuando se quiera hacer seleccion de jets en caso de datos reales 05032013
	  Handle<pat::METCollection> pfMet;
	  iEvent.getByLabel("patMETsPF", pfMet);
	  const pat::METCollection* PFMET = pfMet.product();
	  const pat::MET& pfMET = (*PFMET)[0];
	  recoToLHCO[5][0]=pfMET.eta();
	  recoToLHCO[5][1]=pfMET.phi();
	  recoToLHCO[5][2]=pfMET.pt();
	  recoToLHCO[5][3]=pfMET.mass();//0;
	  recoToLHCO[5][4]=0;
	  recoEnergy[5]=pfMET.energy();
	  //file<<"pfmetJJJPPP"<<pfMET.energy()<<endl;
	  cout<<"paso1"<<endl;
	  
	  
	  Handle<pat::MuonCollection> pfMuon;
	  iEvent.getByLabel("semilepMuonsPF", pfMuon);
	  const pat::MuonCollection* PFMuon = pfMuon.product();
	  const pat::Muon& pfMUON = (*PFMuon)[0];
	  recoToLHCO[4][0]=pfMUON.eta();
	  recoToLHCO[4][1]=pfMUON.phi();
	  recoToLHCO[4][2]=pfMUON.pt();
	  recoToLHCO[4][3]=pfMUON.mass();//0;
	  recoToLHCO[4][4]=0;
	  recoEnergy[4]=pfMUON.energy();
	  //file<<"pfMuonJJJPPP"<<pfMUON.energy()<<endl;
	  cout<<"paso2"<<endl;
	  
	  
	  Handle<pat::JetCollection> pfJets;
	  iEvent.getByLabel("selectedPatJetsPF", pfJets);
	  const pat::JetCollection* PFJet = pfJets.product();
	  int numberOfJets = PFJet->size();
	  bool negativeDiscriminant= false;
	  const int numberOfCombinations = 100; 
	  double W_Combinations[numberOfCombinations][numberOfJets];//en realidad es numberOfJets-1 ya que se necesitan al menos 2 para b1 y b2 pero se gasta la primera casilla escribiendo el deltaWHad
	  double b1_Combinations[numberOfCombinations][numberOfJets];//en realidad es numberOfJets-2 ya que se necesitan al menos 2 para el W y uno para b2 pero se gasta la primera casilla escribiendo el deltab1, pero por facilidad
	  double b2_Combinations[numberOfCombinations][numberOfJets];//en realidad es numberOfJets-2 ya que se necesitan al menos 2 para el W y uno para b1 pero se gasta la primera casilla escribiendo el deltab2, pero por facilidad
	  double deltaWHad, deltab1, deltab2;
	  bool endOfComparison;
	  for (int j=0; j<numberOfCombinations;j++)
	    {
	      for (int k=0; k<numberOfJets;k++)
		{
                  W_Combinations[j][k]=1000;
                  b1_Combinations[j][k]=1000;
                  b2_Combinations[j][k]=1000;
		} 	
	    }	
          int a, b, c, d, e, f;	  
          if (RFijo)
	    {
	      a=1;
	      b=2;
	    }
          else
	    {
	      a=0;
	      b=numberOfJets-2;
	    }

 	  //for (int jetNumber=0; jetNumber<numberOfJets-2; jetNumber++)  //R variable
	  //for (int jetNumber=1;jetNumber<2;jetNumber++)  //R fijo 
          for (int jetNumber=a; jetNumber<b; jetNumber++)  
	    {	
	      int N = numberOfJets;
	      int R1 = jetNumber;
	      const int N1 = R1+1;
	      //prueba<<"NW="<<N<<" CW="<<N1<<endl; 
	      cout<<"NW="<<N<<" CW="<<N1<<endl;
	      vector<int> combinationWHad(N1);
	      vector<int> othersWHad(N-N1);
	      for (int j=0; j<N1; j++)
		{
		  combinationWHad[j]=j + 1;
		}
	      for (int j=0; j<N-N1; j++)
		{
		  othersWHad[j]=N1 + j + 1;
		}
	      
	      lastCombination = false;
	      while (!lastCombination)
		{
		  //                                                                                                                                                                                   
		  for (int j=0; j<N1;j++)
		    {
		      //prueba<<"combWHad "<<combinationWHad[j]<<" ";
		      cout<<"combWHad "<<combinationWHad[j]<<" ";
		    }
		  //prueba<<endl;
		  cout<<endl;
		  for(int j=0; j<N-N1;j++)
		    {
		      //prueba<<"othersWHad "<<othersWHad[j]<<" ";
		      cout<<"othersWHad "<<othersWHad[j]<<" ";
		    }
		  //prueba<<endl;
		  cout<<endl;
		  //                                                                                                                                                                                    
		  math::PtEtaPhiMLorentzVectorD WHad (0,0,0,0);
		  for (int j = 0; j < N1; j++)
		    {
		      const pat::Jet& pfJETtmp = (*PFJet)[combinationWHad[j]-1];
                      if (j<2)
			{
                          discriminants[j+2]=pfJETtmp.bDiscriminator(bdisc_name); 
			}
	  	      math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		      WHad = jetTmp + WHad;
		    }
		  
		  //prueba<<"massWHad="<<WHad.M()<<endl;
		  cout<<"massWHad="<<WHad.M()<<endl;
		  deltaWHad=abs(WHad.M()-80.385);//73.29);
		  
		  //////
		  //cout<<"paso1"<<endl;
		  int M = N - N1;
		  if (RFijo)
		    {
		      c=0;
		      d=1;
		    }
		  else
		    {
		      c=0;
		      d=M-1;
		    }
		  //for (int l=0; l<M-1; l++) //R Variable
		  //for (int l=0;l<1;l++) //R fijo
		  for (int l=c; l<d; l++)    
		    {
		      //cout<<"paso2"<<endl;
		      int R2 = l;
		      const int N2 = R2+1;
		      //prueba<<"Nb1="<<M<<" Cb1="<<N2<<endl;
		      cout<<"Nb1="<<M<<" Cb1="<<N2<<endl;
	              vector<int> combinationb1(N2);
		      vector<int> othersb1(M-N2);
		      for (int k=0; k<N2; k++)
			{
			  combinationb1[k]=k + 1;
			}
		      for (int k=0; k<M-N2; k++)
			{
			  othersb1[k]=N2 + k + 1;
			}
		      
		      lastCombination2= false;
		      while (!lastCombination2)
			{
			  //cout<<"paso3"<<endl;
			  //                                                                                                                                                                    
			  for (int k=0; k<N2;k++)
			    {
			      //prueba<<"combb1 "<<othersWHad[combinationb1[k]-1]<<" ";
			      cout<<"combb1 "<<othersWHad[combinationb1[k]-1]<<" ";
			    }
			  //prueba<<endl;
			  cout<<endl;
			  for(int k=0; k<M-N2;k++)
			    {
			      //prueba<<"othersb1 "<<othersWHad[othersb1[k]-1]<<" ";
			      cout<<"othersb1 "<<othersWHad[othersb1[k]-1]<<" ";
			    }
			  //prueba<<endl;
			  cout<<endl;
			  //   
			  math::PtEtaPhiMLorentzVectorD b1 (0,0,0,0);
			  for (int k = 0; k < N2; k++)
			    {
			      const pat::Jet& pfJETtmp = (*PFJet)[othersWHad[combinationb1[k]-1]-1];
			      if (k<1)
				{
				  discriminants[k]=pfJETtmp.bDiscriminator(bdisc_name);
				}
			      math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
			      b1 = jetTmp + b1;
			    }
			  
			  //prueba<<"massb1="<<b1.M()<<endl;
			  cout<<"massb1="<<b1.M()<<endl;
			  
			  
			  deltab1=abs(b1.M()-4.7);//12.46);
                          
			  
			  //////
			  //cout<<"paso1a"<<endl;
			  int M1 = M - N2;
			  int tmp = M1 -1;  
			  //cout<<"M="<<M<<" N2="<<N2<<" M1="<<M1<<" M1-N3="<<tmp<<endl;
			  if (RFijo)
			    {
			      e=0;
			      f=1;
			    }
			  else
			    {
			      e=0;
			      f=M1;
			    }
			  //for (int l1=0; l1<M1; l1++)//R variable
			  //for (int l1=0;l1<1;l1++)//R fijo 
			  for (int l1=e; l1<f; l1++)  
			    {
			      //cout<<"paso2a"<<endl;
			      int R3 = l1;
			      const int N3 = R3+1;
			      const int tmp1=M1-N3;
			      //cout<<"M1-N3="<<tmp1<<endl;
			      //prueba<<"Nb2="<<M1<<" Cb2="<<N3<<endl;
			      cout<<"Nb2="<<M1<<" Cb2="<<N3<<endl;
			      vector<int> combinationb2(N3);
			      cout<<"paso2b"<<endl;
			      vector<int> othersb2(M1-N3);
			      //cout<<"paso2c"<<endl;
			      for (int k1=0; k1<N3; k1++)
				{
				  //cout<<"paso2d"<<endl;
				  combinationb2[k1]=k1 + 1;
				}
			      //cout<<"paso2e"<<endl;
			      for (int k1=0; k1<M1-N3; k1++)
				{
				  //cout<<"paso2f"<<endl;
				  othersb2[k1]=N3 + k1 + 1;
				}
			      //cout<<"paso2g"<<endl;
			      lastCombination3= false;
			      while (!lastCombination3)
				{
				  //cout<<"paso3a"<<endl;
				  //                                                                                                                                                       
				  for (int k1=0; k1<N3;k1++)
				    {
				      //cout<<"paso3b"<<endl;
				      //prueba<<"combb2 "<<othersWHad[othersb1[combinationb2[k1]-1]-1]<<" ";
				      cout<<"combb2 "<<othersWHad[othersb1[combinationb2[k1]-1]-1]<<" ";
				    }
				  //cout<<"paso3c"<<endl;
				  //prueba<<endl;
				  cout<<endl;
				  for(int k1=0; k1<M1-N3;k1++)
				    {
				      //cout<<"paso3d"<<endl;
				      //prueba<<"othersb2 "<<othersWHad[othersb1[othersb2[k1]-1]-1]<<" ";
				      cout<<"othersb2 "<<othersWHad[othersb1[othersb2[k1]-1]-1]<<" ";
				    }
				  //cout<<"paso3e"<<endl;
				  //prueba<<endl;
				  cout<<endl;
				  // 
				  math::PtEtaPhiMLorentzVectorD b2 (0,0,0,0);
				  //cout<<"paso3f"<<endl;
				  for (int k1 = 0; k1 < N3; k1++)
				    {
				      //cout<<"paso3g"<<endl;
				      const pat::Jet& pfJETtmp = (*PFJet)[othersWHad[othersb1[combinationb2[k1]-1]-1]-1];
				      if (k1<1)
					{
					  discriminants[k1+1]=pfJETtmp.bDiscriminator(bdisc_name);
					}
				      math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
				      b2 = jetTmp + b2;
				    }
				  //prueba<<"massbLep="<<b2.M()<<endl;
				  cout<<"massbLep="<<b2.M()<<endl;
				  //cout<<"paso3h"<<endl;
                                  double deltaRJetLep_bLep=pow(b2.eta()-recoToLHCO[4][0],2)+pow(b2.phi()-recoToLHCO[4][1],2);
                                 
				
				  double deltaPhiJetLep_bLep=abs(b2.phi()-recoToLHCO[4][1]);
                                  if (deltaPhiJetLep_bLep>3.141592)
                                    {
                                      deltaPhiJetLep_bLep=deltaPhiJetLep_bLep-3.1416;
                                    }
                                  double deltaPhiJetMET_bLep=abs(b2.phi()-recoToLHCO[5][1]);
                                  if (deltaPhiJetMET_bLep>3.141592)
                                    {
                                      deltaPhiJetMET_bLep=deltaPhiJetMET_bLep-3.1416;
                                    }
                                  double deltaPhi_bLep= deltaPhiJetLep_bLep + deltaPhiJetMET_bLep;
                                  

				  double deltaPhiJetLep_bHad_=abs(b1.phi()-recoToLHCO[4][1]);
                                  if (deltaPhiJetLep_bHad_>3.141592)
                                    {
                                      deltaPhiJetLep_bHad_=deltaPhiJetLep_bHad_-3.1416;
                                    }
                                  double deltaPhiJetMET_bHad_=abs(b1.phi()-recoToLHCO[5][1]);
                                  if (deltaPhiJetMET_bHad_>3.141592)
                                    {
                                      deltaPhiJetMET_bHad_=deltaPhiJetMET_bHad_-3.1416;
                                    }
                                  double deltaPhi_bHad_= deltaPhiJetLep_bHad_ + deltaPhiJetMET_bHad_;

				  math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[4][2],recoToLHCO[4][0],recoToLHCO[4][1],recoToLHCO[4][3]);
				  math::PtEtaPhiMLorentzVectorD met (recoToLHCO[5][2],recoToLHCO[5][0],recoToLHCO[5][1],recoToLHCO[5][3]);
				  math::PtEtaPhiMLorentzVectorD WLep = ele + met ;
				  math::PtEtaPhiMLorentzVectorD tLep = WLep + b2;
				  math::PtEtaPhiMLorentzVectorD tHad = WHad + b1;
                                  double deltatHad=abs(tHad.M()-172.5);
				  double deltatLep=abs(tLep.M()-172.5);

                    
				  double deltaPhiJetQ_bHad=abs(b1.phi()-recoToLHCO[2][1]);
                                  if (deltaPhiJetQ_bHad>3.141592)
                                    {
                                      deltaPhiJetQ_bHad=deltaPhiJetQ_bHad-3.1416;
                                    }
                                  double deltaPhiJetQ__bHad=abs(b1.phi()-recoToLHCO[3][1]);
                                  if (deltaPhiJetQ__bHad>3.141592)
                                    {
                                      deltaPhiJetQ__bHad=deltaPhiJetQ__bHad-3.1416;
                                    }
                                  double deltaPhi_bHad= deltaPhiJetQ_bHad + deltaPhiJetQ__bHad;
				  double deltaPhi_WbHad= abs(b1.phi() -WHad.phi());
				  if (deltaPhi_WbHad>3.141592)
                                    {
                                      deltaPhi_WbHad=deltaPhi_WbHad-3.1416;
                                    }
				  //deltaPhi_bHad = deltaPhi_bHad - 2.648;
                                  //deltaPhi_bLep = deltaPhi_bLep - 3.352;
                                  //deltaPhi_WbHad = deltaPhi_WbHad - 3.582;
                                  //prueba<<"massbLep="<<b2.M()<<endl;
				  cout<<"massbLep="<<b2.M()<<endl;
				  //cout<<"paso3h"<<endl;
				  deltab2=abs(b2.M()-4.7);//10.77);
				  //double discriminantJP = sqrt(pow(deltatHad,2.0)/30 + pow(deltaWHad,2.0)/11.75 + pow(deltab1,2.0)/3.36 + pow(deltab2,2.0)/4.1);//Solo neutralinos
				  //	  double discriminantJP = (pow(deltaWHad,2.0)/11.76 + pow(deltab1,2.0)/3.36 + pow(deltab2,2.0)/4.15)*(pow(deltaPhi_bLep,2.0) +1.26*pow(deltaPhi_bHad,2.0));//el que mejor va
                                  //double discriminantJP =  (( pow(deltaWHad,2.0)/10.29 + pow(deltab1,2.0)/2.587 + pow(deltab2,2.0)/3.926 )+100*(0.15*pow(deltaPhi_bLep,2.0)/0.6071 +pow(deltaPhi_bHad,2.0)/1.156 ))/100000;La que estaba cuadrando
				  //double discriminantJP = ( (pow(deltaWHad,2.0)/11.76 + pow(deltab1,2.0)/3.36 + pow(deltab2,2.0)/4.15 + 1)*(pow(deltaPhi_bLep,2.0) +1*pow(2*deltaPhi_WbHad,2.0) + 1) - 1 )/1000000;  

				  double bTagDisc=fit_b(discriminants[0])*fit_b(discriminants[1])*fit_cl(discriminants[2])*fit_cl(discriminants[3]);
				  double AntibTagDisc=fit_cl(discriminants[0])*fit_cl(discriminants[1])*fit_b(discriminants[2])*fit_b(discriminants[3]);
				  double bTagDiscriminant=10000-bTagDisc;
				  double discriminantJP;
				  // double discriminantJP =(3*(1*(100*(bTagDiscriminant + 0.875*(deltaPhi_bLep/deltaPhi_bHad_)) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/100000;//OJO LA QUE SE DEBE UTILIZAR PARA ANALISIS ALTERNO.
				  //double discriminantJP =(1*(1*(100*(10*bTagDiscriminant + 0.875*((0.1+deltaPhi_bLep)/(0.1+deltaPhi_bHad_))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) +2*pow(deltaWHad,2.0))/100000;
				  if (selector==1)
				    {
				      discriminantJP =(3*(1*(300*(bTagDiscriminant+ 1875*((3.241592)/(0.1+deltaPhi_bHad_))  + 18750*((deltaPhi_bHad)/(3.141592)) + 0.875*((deltaPhi_bLep)/(3.141592))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/1000000;//MATCHING EN 6
				    }
				  if (selector==2)
				    {
				      discriminantJP =(0.1*(1*(300*(10*bTagDiscriminant+ 187500000*((3.241592)/(0.1+deltaPhi_bHad_))  + 18750000000*((deltaPhi_bHad)/(3.141592)) + 87.5*((deltaPhi_bLep)/(3.141592))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/100000000000000000; 
				    }
				  if (selector==3)
				    {
				      discriminantJP =(3*(1*(100*(bTagDiscriminant + 87.5000*((3.2416)/(0.1+deltaPhi_bHad_))  + 10*(1875000000*((deltaPhi_bHad)/(3.1416))+ 0.875*(deltaPhi_bLep/3.1416))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/1000000000000000000;
				    }
				  if (selector==4)
                                    {
                                      discriminantJP =100/(bTagDisc+100) +  3.2416/(0.1+deltaPhi_bHad_)  + (deltaPhi_bHad)/(3.1416) + (deltaPhi_bLep)/(3.1416) + pow(deltab2-4,2.0)/4.0 + pow(deltab1-4,2.0)/4.0 + pow(deltaWHad,2.0)/20;
				    }
				  if (selector==5)
				    {
				      discriminantJP = 4*100/(bTagDisc+100) + 2* 3.2416/(0.1+deltaPhi_bHad_)  +8*((deltaPhi_bHad)/(3.1416)+ (deltaPhi_bLep)/(3.1416)) + pow(deltab2-4,2.0)/4.0 + pow(deltab1-4,2.0)/4.0 + pow(deltaWHad,2.0)/20;
				    }
				  if (selector==7)
                                    {
                                      discriminantJP = (c1*c10/(bTagDisc+c10) + c2* c11/(c11+deltaPhi_bHad_)  +c3*(deltaPhi_bHad)/(3.1416)+ c4*(deltaPhi_bLep)/(3.1416) + c5*pow(deltab2-c9,2.0)/c13 +c6*pow(deltab1-c8,2.0)/c14 +c7*pow(deltaWHad,2.0)/20)/c12;
                                    }

				
                                  
                                  if (selector==9)
                                    {
                                      discriminantJP = pow(AntibTagDisc,c1)*pow(deltaRJetLep_bLep,c2)*pow(deltaWHad,2)/c3;
                                    }
				  if (selector==10)
                                    {
                                      discriminantJP = pow(AntibTagDisc,c1)*pow(deltaRJetLep_bLep,c2)*(pow(deltaWHad,2)+pow(deltatHad,2))/c3;
                                    }



				  //double bTagDiscriminant=exp(-fit_b(discriminants[0])*fit_b(discriminants[1])*fit_cl(discriminants[2])*fit_cl(discriminants[3])/1000);
                                  //double discriminantJP =(bTagDiscriminant + 0.5*(deltaPhi_bLep/deltaPhi_bHad_));


				  //Double  discriminantJP = (pow(deltaWHad,2.0)/11.76 + 10*1.26*pow(deltaPhi_bHad,2.0)*pow(deltab1,2.0)/3.36 + 10*pow(deltaPhi_bLep,2.0)*pow(deltab2,2.0)/4.15)/100000.0;

				  //*(pow(deltaPhi_bLep,1.0) + pow(deltaPhi_bHad,1.0)); 
				  //double discriminantJP = sqrt(pow(deltatHad,2.0)/17.55 + pow(deltatLep,2.0)/45.11 + pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316);//*(pow(deltaPhi_bLep,1.0) + pow(deltaPhi_bHad,1.0));
				  //double discriminantJP = sqrt(pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)*(pow(deltaPhi_bLep,1.0) + pow(deltaPhi_bHad,1.0));
				  //double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)*(pow(deltaPhi_bLep,2.0) + pow(deltaPhi_bHad,2.0));
				  //double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)*(pow(deltaPhi_bLep,2.0)/pow(deltaPhi_bHad_,2.0)); 
				  
				  //double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + (deltaPhi_bHad/3.2)*pow(deltab1,2.0)/4.175 + (deltaPhi_bLep/3.2)*pow(deltab2,2.0)/4.316);   
              
     

				  //          double discriminantJP=pow(deltaWHad,2.0)/13.24+pow(deltab1,2.0)/4.175+pow(deltab2,2.0)/4.316+pow(deltaPhi_bLep,2.0)/1;
				  //          double discriminantJP =  pow(deltaWHad,2.0)/12.57 + pow(deltab1,2.0)/4.434 + (deltaPhi_bLep/3.2)*pow(deltab2,2.0)/4.412 ;
				  //        double discriminantJP=(deltaWHad/13.24+1)*(deltaPhi_bHad/3.2*deltab1/4.175+1)*(deltaPhi_bLep/3.2*deltab2/4.316+1);
		      
				 
				  ////
				  endOfComparison = false;
				  int j=0;
				  while (!endOfComparison)
				    {
				      //cout<<"j="<<j<<endl;
				      //prueba<<"discriminantJP= "<<discriminantJP<<" WComb[j][0]= "<<W_Combinations[j][0]<<endl;
				      cout<<"discriminantJP= "<<discriminantJP<<" WComb[j][0]= "<<W_Combinations[j][0]<<" j="<<j<<endl;
				      if (discriminantJP<W_Combinations[j][0])
					{
					  //prueba<<"Clasifica Combinacion"<<endl;
					  cout<<"Clasifica Combinacion"<<endl;
					  endOfComparison=true;
					  for (int k=numberOfCombinations-2;k>j-1;k--)
					    {
					      //cout<<"k="<<k<<endl;
					      for (int l=0; l<numberOfJets;l++)
						{
						  //cout<<"l="<<l<<endl;
						  W_Combinations[k+1][l]=W_Combinations[k][l];
						  b1_Combinations[k+1][l]=b1_Combinations[k][l];
						  b2_Combinations[k+1][l]=b2_Combinations[k][l];
						}
					    }
					  for (int l=1; l<N1+1;l++)
					    {
					      //cout<<"lW="<<l<<endl;
					      W_Combinations[j][l]=combinationWHad[l-1];
					      //prueba<<" JPW "<<combinationWHad[l-1];
					    }
					  //prueba<<endl;
				  
					  for (int l=1; l<N2+1;l++)
					    {
					      //cout<<"lb1="<<l<<endl;
					      b1_Combinations[j][l]=othersWHad[combinationb1[l-1]-1];//combinationb1[l-1];
					      //prueba<<" JPb1 "<<combinationb1[l-1];
					    }
					  //prueba<<endl;
					  
					  for (int l=1; l<N3+1;l++)
					    {
					      //cout<<"lb2="<<l<<endl;
					      b2_Combinations[j][l]=othersWHad[othersb1[combinationb2[l-1]-1]-1];//combinationb2[l-1];
					      //prueba<<" JPb2 "<<combinationb2[l-1];
					    }
					  //prueba<<endl;
				  
					  W_Combinations[j][0]=discriminantJP;
					  for(int l=N1+1; l<numberOfJets;l++)
					    {
					      //cout<<"lWo="<<l<<endl;
					      W_Combinations[j][l]=1000;
					    }
				  
					  for(int l=N2+1; l<numberOfJets;l++)
					    {
					      //cout<<"lb1o="<<l<<endl;
					      b1_Combinations[j][l]=1000;
					    }
				  
					  for(int l=N3+1; l<numberOfJets;l++)
					    {
					      //cout<<"lb2o="<<l<<endl;
					      b2_Combinations[j][l]=1000;
					    }
				  
					}
				      if (j<99)
					{
					  j++;
					}
				      else
					{
					  endOfComparison = true;
					}
				    }
				  ////
			    
				  cout<<"paso5a"<<endl;
				  lastCombination3 = nextComb(R3, combinationb2, othersb2, M1, R3);
				  cout<<"paso6a"<<lastCombination3<<endl;
				}
			    }
			  
			  //////
			  
			  
			  //cout<<"paso5"<<endl;
			  cout<<"lastCombBefore="<<lastCombination2<<" R2="<<R2<<endl;
			  lastCombination2 = nextComb(R2, combinationb1, othersb1, M, R2);
			  cout<<"lastCombAfter="<<lastCombination2<<" R2="<<R2<<endl;
			}
		    }
		  
		  lastCombination = nextComb(R1, combinationWHad, othersWHad, N, R1);  
		  
		}
	    }
	  double mass[3];
	  double minDiscriminantJP=1000;
	  double sel_deltaR[4];
	  bool all_assigned;  
	  for (int j=0; j<numberOfCombinations;j++)
	    {
	      double deltaR_temp[4];
	      double deltaR[4];
	      bool q_assigned = false;
	      bool q__assigned = false;
	      bool b_assigned = false;
	      bool b__assigned = false;
	      bool endOfLine=false;
	      int l=0;
	      //prueba<<"combinacion W j="<<j<<endl;
	      math::PtEtaPhiMLorentzVectorD WHad (0,0,0,0);
	      while (!endOfLine)
		{
		  //prueba<<W_Combinations[j][l]<<" ";
                  if (l>0)
		    {
		      const pat::Jet& pfJETtmp = (*PFJet)[W_Combinations[j][l]-1];
		      
                      if (l<3)
			{
                          recoToLHCO[l+1][0]=pfJETtmp.eta();
			  recoToLHCO[l+1][1]=pfJETtmp.phi();
			  recoToLHCO[l+1][2]=pfJETtmp.pt();
			  recoToLHCO[l+1][3]=pfJETtmp.mass();//0;
			  recoToLHCO[l+1][4]=0;
			  recoEnergy[l+1]=pfJETtmp.energy();
                          discriminants[l+1]=pfJETtmp.bDiscriminator(bdisc_name);
		        }



                      deltaR_temp[2]=sqrt(pow(pfJETtmp.eta()-genToLHCO[2][0],2.0)+pow(pfJETtmp.phi()-genToLHCO[2][1],2.0)+pow(pfJETtmp.pt()-genToLHCO[2][2],2.0)+pow(pfJETtmp.mass()-genToLHCO[2][3],2.0));
                      deltaR_temp[3]=sqrt(pow(pfJETtmp.eta()-genToLHCO[3][0],2.0)+pow(pfJETtmp.phi()-genToLHCO[3][1],2.0)+pow(pfJETtmp.pt()-genToLHCO[3][2],2.0)+pow(pfJETtmp.mass()-genToLHCO[3][3],2.0));
                      if (deltaR_temp[2]<deltaR_temp[3])
			{
			  deltaR[2]=deltaR_temp[2];
			  q_assigned=true; 
			}
                      else
			{
			  deltaR[3]=deltaR_temp[3];
			  q__assigned=true;
			}
		      math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		      WHad = jetTmp + WHad;
		    } 
		  l++;
		  if (W_Combinations[j][l]==1000)
		    {
		      endOfLine=true;
		    }
		}
	      //prueba<<endl;
	      //prueba<<"combinacion b1 j="<<j<<endl;
              math::PtEtaPhiMLorentzVectorD b1 (0,0,0,0);
	      endOfLine=false;
	      l=0;
	      while (!endOfLine)
		{
		  //prueba<<b1_Combinations[j][l]<<" ";
                  if (l>0)
		    {
		      const pat::Jet& pfJETtmp = (*PFJet)[b1_Combinations[j][l]-1];
                      if (l==1)
			{
			  recoToLHCO[0][0]=pfJETtmp.eta();
			  recoToLHCO[0][1]=pfJETtmp.phi();
			  recoToLHCO[0][2]=pfJETtmp.pt();
			  recoToLHCO[0][3]=pfJETtmp.mass();//4.7;
			  recoToLHCO[0][4]=2;
			  recoEnergy[0]=pfJETtmp.energy();
			  discriminants[0]=pfJETtmp.bDiscriminator(bdisc_name);
			}
			  
                      deltaR_temp[0]=sqrt(pow(pfJETtmp.eta()-genToLHCO[0][0],2.0)+pow(pfJETtmp.phi()-genToLHCO[0][1],2.0)+pow(pfJETtmp.pt()-genToLHCO[0][2],2.0)+pow(pfJETtmp.mass()-genToLHCO[0][3],2.0));
                      deltaR_temp[1]=sqrt(pow(pfJETtmp.eta()-genToLHCO[1][0],2.0)+pow(pfJETtmp.phi()-genToLHCO[1][1],2.0)+pow(pfJETtmp.pt()-genToLHCO[1][2],2.0)+pow(pfJETtmp.mass()-genToLHCO[1][3],2.0));
                      if (deltaR_temp[0]<deltaR_temp[1])
			{
			  deltaR[0]=deltaR_temp[0];
			  b_assigned=true; 
			}
                      else
			{
			  deltaR[1]=deltaR_temp[1];
			  b__assigned=true;
			}
		      math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		      b1 = jetTmp + b1;
		    } 
		  l++;
		  if (b1_Combinations[j][l]==1000)
		    {
		      endOfLine=true;
		    }
		}
	      //prueba<<endl;
	      //prueba<<"combinacion b2 j="<<j<<endl;
              math::PtEtaPhiMLorentzVectorD b2 (0,0,0,0);
	      endOfLine=false;
	      l=0;
	      while (!endOfLine)
		{
		  //prueba<<b2_Combinations[j][l]<<" ";
		  if (l>0)
		    {
		      const pat::Jet& pfJETtmp = (*PFJet)[b2_Combinations[j][l]-1];
                      if (l==1)
			{
			  recoToLHCO[1][0]=pfJETtmp.eta();
			  recoToLHCO[1][1]=pfJETtmp.phi();
			  recoToLHCO[1][2]=pfJETtmp.pt();
			  recoToLHCO[1][3]=pfJETtmp.mass();//4.7;
			  recoToLHCO[1][4]=2;
			  recoEnergy[1]=pfJETtmp.energy();
			  discriminants[1]=pfJETtmp.bDiscriminator(bdisc_name); 
			}

                      deltaR_temp[0]=sqrt(pow(pfJETtmp.eta()-genToLHCO[0][0],2.0)+pow(pfJETtmp.phi()-genToLHCO[0][1],2.0)+pow(pfJETtmp.pt()-genToLHCO[0][2],2.0)+pow(pfJETtmp.mass()-genToLHCO[0][3],2.0));
                      deltaR_temp[1]=sqrt(pow(pfJETtmp.eta()-genToLHCO[1][0],2.0)+pow(pfJETtmp.phi()-genToLHCO[1][1],2.0)+pow(pfJETtmp.pt()-genToLHCO[1][2],2.0)+pow(pfJETtmp.mass()-genToLHCO[1][3],2.0));
                      if (deltaR_temp[0]<deltaR_temp[1])
			{
			  deltaR[0]=deltaR_temp[0];
			  b_assigned=true; 
			}
                      else
			{
			  deltaR[1]=deltaR_temp[1];
			  b__assigned=true;
			  if (deltaR[1]<0)
			    {
			      //prueba<<"JJJJPPPP";
			    }
			}
		      math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		      b2 = jetTmp + b2;
		    } 
                  l++;
		  if (b2_Combinations[j][l]==1000)
		    {
		      endOfLine=true;
		    }
		}
	      //prueba<<endl;
	      
              //double discriminantJP = (abs(WHad.M()-80.385) + 1000)*(abs(bHad.M()-4.7) + 1000)*(abs(bLep.M()-4.7) + 1000);
              //double discriminantJP = sqrt(pow(WHad.M()-73.29,2.0)+pow(bHad.M()-12.46,2.0)+pow(bLep.M()-10.77,2.0));
              /*0404
	      double discriminantJP = sqrt(pow(WHad.M()-80.385,2.0)+pow(b1.M()-4.7,2.0) +pow(b2.M()-4.7,2.0));
	      
	      if (j==0)
		{
		  discriminantJPSelection->Fill(discriminantJP);
		}
	      discriminantJPvsN->Fill(j,discriminantJP);
              if (discriminantJP < minDiscriminantJP)
		{
		  minDiscriminantJP=discriminantJP;
		  mass[0]=WHad.M();
		  mass[1]=b1.M();
		  mass[2]=b2.M();
		  for (int l=0; l<4; l++)
		    {
		      sel_deltaR[l]=deltaR[l];
		    }
		  if ( b_assigned && b__assigned && q_assigned && q__assigned)
		    {
		      all_assigned = true;
		    }
		  else
		    {
		      all_assigned = false;
		    }  
		}
              
              if (j==0)
		{
	    
		  ///JPPP
		  for (int i=0;i<2;i++)
		    {
		      //Difference2_InEnergy_bJets_inclusive_vs_PartonEnergy->Fill(genEnergy[i],recoEnergy[i]-genEnergy[i]);
		      Difference_InEnergy_bJets_Gen_RecoS->Fill(recoEnergy[i]-genEnergy[i]);
		      Difference_InPt_bJets_Gen_RecoS->Fill(recoToLHCO[i][2]-genToLHCO[i][2]);
		      Difference_InEta_bJets_Gen_RecoS->Fill(recoToLHCO[i][0]-genToLHCO[i][0]);
		      Difference_InPhi_bJets_Gen_RecoS->Fill(recoToLHCO[i][1]-genToLHCO[i][1]);
		      if (genEnergy[i]<=20)
			{ 		  
			  Difference_InEnergy_bJets_Gen_Reco_20S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>20)&&(genEnergy[i]<40))
			{ 		  
			  Difference_InEnergy_bJets_Gen_Reco_40S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>40)&&(genEnergy[i]<60))
			{ 		  
			  Difference_InEnergy_bJets_Gen_Reco_60S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>60)&&(genEnergy[i]<80))
			{ 		  
			  Difference_InEnergy_bJets_Gen_Reco_80S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>80)&&(genEnergy[i]<100))
			{ 		  
			  Difference_InEnergy_bJets_Gen_Reco_100S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>100)&&(genEnergy[i]<120))
			{ 		  
			  Difference_InEnergy_bJets_Gen_Reco_120S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>120)&&(genEnergy[i]<140))
			{
			  Difference_InEnergy_bJets_Gen_Reco_140S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>140)&&(genEnergy[i]<160))
			{
			  Difference_InEnergy_bJets_Gen_Reco_160S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>160)&&(genEnergy[i]<180))
			{
			  Difference_InEnergy_bJets_Gen_Reco_180S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>180)&&(genEnergy[i]<200))
			{
			  Difference_InEnergy_bJets_Gen_Reco_200S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>200)&&(genEnergy[i]<220))
			{
			  Difference_InEnergy_bJets_Gen_Reco_220S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>220)&&(genEnergy[i]<240))
			{
			  Difference_InEnergy_bJets_Gen_Reco_240S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>240)&&(genEnergy[i]<260))
			{
			  Difference_InEnergy_bJets_Gen_Reco_260S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>260)&&(genEnergy[i]<280))
			{
			  Difference_InEnergy_bJets_Gen_Reco_280S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>280)&&(genEnergy[i]<300))
			{
			  Difference_InEnergy_bJets_Gen_Reco_300S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>300)&&(genEnergy[i]<320))
			{
			  Difference_InEnergy_bJets_Gen_Reco_320S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		  
		  
		    } 
	      
		  for (int i=2;i<4;i++)
		    {
		      //Difference2_InEnergy_lJets_inclusive_vs_PartonEnergy->Fill(genEnergy[i],recoEnergy[i]-genEnergy[i]);
		      Difference_InEnergy_clJets_Gen_RecoS->Fill(recoEnergy[i]-genEnergy[i]);
		      Difference_InPt_clJets_Gen_RecoS->Fill(recoToLHCO[i][2]-genToLHCO[i][2]);
		      Difference_InEta_clJets_Gen_RecoS->Fill(recoToLHCO[i][0]-genToLHCO[i][0]);
		      Difference_InPhi_clJets_Gen_RecoS->Fill(recoToLHCO[i][1]-genToLHCO[i][1]);
		      if (genEnergy[i]<=20)
			{ 		  
			  Difference_InEnergy_clJets_Gen_Reco_20S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>20)&&(genEnergy[i]<40))
			{ 		  
			  Difference_InEnergy_clJets_Gen_Reco_40S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>40)&&(genEnergy[i]<60))
			{ 		  
			  Difference_InEnergy_clJets_Gen_Reco_60S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>60)&&(genEnergy[i]<80))
			{ 		  
			  Difference_InEnergy_clJets_Gen_Reco_80S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>80)&&(genEnergy[i]<100))
			{ 		  
			  Difference_InEnergy_clJets_Gen_Reco_100S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>100)&&(genEnergy[i]<120))
			{ 		  
			  Difference_InEnergy_clJets_Gen_Reco_120S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>120)&&(genEnergy[i]<140))
			{
			  Difference_InEnergy_clJets_Gen_Reco_140S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>140)&&(genEnergy[i]<160))
			{
			  Difference_InEnergy_clJets_Gen_Reco_160S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>160)&&(genEnergy[i]<180))
			{
			  Difference_InEnergy_clJets_Gen_Reco_180S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>180)&&(genEnergy[i]<200))
			{
			  Difference_InEnergy_clJets_Gen_Reco_200S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>200)&&(genEnergy[i]<220))
			{
			  Difference_InEnergy_clJets_Gen_Reco_220S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>220)&&(genEnergy[i]<240))
			{
			  Difference_InEnergy_clJets_Gen_Reco_240S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>240)&&(genEnergy[i]<260))
			{
			  Difference_InEnergy_clJets_Gen_Reco_260S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>260)&&(genEnergy[i]<280))
			{
			  Difference_InEnergy_clJets_Gen_Reco_280S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>280)&&(genEnergy[i]<300))
			{
			  Difference_InEnergy_clJets_Gen_Reco_300S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		      if ((genEnergy[i]>300)&&(genEnergy[i]<320))
			{
			  Difference_InEnergy_clJets_Gen_Reco_320S->Fill(recoEnergy[i]-genEnergy[i]);
			}
		  
		  
		    }
		  Difference_InEta_MET_Gen_RecoS->Fill(recoToLHCO[5][0]-genToLHCO[5][0]);
		  Difference_InPhi_MET_Gen_RecoS->Fill(recoToLHCO[5][1]-genToLHCO[5][1]);
		  Difference_InPt_MET_Gen_RecoS->Fill(recoToLHCO[5][2]-genToLHCO[5][2]);
		  Difference_InMass_MET_Gen_RecoS->Fill(recoToLHCO[5][3]-genToLHCO[5][3]);
		  Difference_InEta_Muon_Gen_RecoS->Fill(recoToLHCO[4][0]-genToLHCO[4][0]);
		  Difference_InPhi_Muon_Gen_RecoS->Fill(recoToLHCO[4][1]-genToLHCO[4][1]);
		  Difference_InPt_Muon_Gen_RecoS->Fill(recoToLHCO[4][2]-genToLHCO[4][2]);
		  Difference_InMass_Muon_Gen_RecoS->Fill(recoToLHCO[4][3]-genToLHCO[4][3]);
		}	      
	      ///JPPP
              0404*/
	    }
          
          /*0404
          Sel_Inv_Mass_W -> Fill(mass[0]);
          Sel_Inv_Mass_b -> Fill(mass[1]);
          Sel_Inv_Mass_b_ -> Fill(mass[2]);
	  
          sel_deltaR_b->Fill(sel_deltaR[0]);
          sel_deltaR_b_->Fill(sel_deltaR[1]);
          sel_deltaR_q->Fill(sel_deltaR[2]);
          sel_deltaR_q_->Fill(sel_deltaR[3]);
          0404*/
	  /*
	    if (all_assigned)
	    {
	    sel_deltaR_b->Fill(-10);
	    sel_deltaR_b_->Fill(-10);
	    sel_deltaR_q->Fill(-10);
	    sel_deltaR_q_->Fill(-10);
	    }
	  */	  
	  //fin para mover cuando se quiera hacer seleccion de jets en caso de datos reales 05032013    
	  
	  
	  
	  
	  /*
	  //ejemplo para quitar
	  const pat::Jet& pfJET = (*PFJet)[0];
	  int parton = pfJET.genParton()->pdgId();
	  file<<"partonJJJPPP"<<parton<<endl;
	  double disc = pfJET.bDiscriminator(bdisc_name);
	  //const std::vector< std::pair< std::string, float > > & discri = pfJET.getPairDiscri(); 
	  file<<"discJJJPPP  "<<disc<<endl;  
	  const reco::GenParticle* genParton= pfJET.genParton();
	  int parton1 = genParton->pdgId();
          file<<"parton1JJJPPP"<<parton1<<endl;
	  cout<<"paso3"<<endl;
          //fin ejemplo
	  */
	  int realConfiguration[4];
	  for (int i=0; i<4 ; i++)
	    {
	      realConfiguration[i]=-1;
	    }
	  
	  //file11<<"gen       Matchgen       reco  "<<endl;
          //int numberOfJets = PFJet->size();
	  
	  /* //////0503
	     N = numberOfJets;
	     R = 2;
	     const int N1 = 3;
	     vector<int> combinacion(N1);
	     for (int jp=0; jp<N1; jp++)
	     {
	     combinacion[jp]=jp + 1;
	     }
	     for (int prueba1 =0; prueba1 < 100; prueba1++)
	     {
	     nextComb(R, combinacion);
	     prueba<<"JP: comb "<<prueba1<<" N= "<<N<<endl;
	     for (int jp=0; jp<N1; jp++)
	     {
	     prueba<<combinacion[jp]<<" ";
	     }
	     prueba<<endl;
             
	     }
	     //////0503end
	     */
	  
	  
	  
	  //bool negativeDiscriminant= false;
	  // for (pat::JetCollection::iterator iter=
          Total_Energy_Jets = 0;
	  for (int jetNumber=0; jetNumber<numberOfJets; jetNumber++) 
	    {
	      
	      const pat::Jet& pfJET = (*PFJet)[jetNumber];
              Total_Energy_Jets = Total_Energy_Jets + pfJET.energy(); 
	      //const reco::GenParticle& genParton =*(pfJET.genParton());
	      const reco::GenParticle* genParton= pfJET.genParton();
	      cout<<"paso3a"<<endl;
              cout<<(float)genToLHCO[jetNumber][0]<<endl;
	      //double eta=(((*PFJet)[jetNumber]).genParton())->eta();
	      //cout<<eta<<endl;
	      //cout<<(float)((pfJET.genParton())->eta())<<endl;
	      if (genParton)
		{
		  cout<<genParton->eta()<<endl;
		  
		  for (int i=0; i<4; i++)
		    {
		      cout<<"paso3a1"<<endl;
		      //file10<<"la diferencia entre el gen y el genMatching es: "<<abs((float)genToLHCO[i][0]-(float)(genParton->eta()))<<endl;
		      cout<<"paso3a2"<<endl;
		      if (abs((float)genToLHCO[i][0]-(float)(genParton->eta()))==0)
			{
			  cout<<"paso3b"<<endl;
			  //file11<<"jet: "<<i<<" eta: "<<genToLHCO[i][0]<<"    "<<genParton->eta()<<"     "<<pfJET.eta()<<endl;
			  //file11<<"jet: "<<i<<" phi: "<<genToLHCO[i][1]<<"    "<<genParton->phi()<<"     "<<pfJET.phi()<<endl;
			  //file11<<"jet: "<<i<<" pt: "<<genToLHCO[i][2]<<"    "<<genParton->pt()<<"     "<<pfJET.pt()<<endl;
			  realConfiguration[i]=jetNumber;
			  recoToLHCO[i][0]=pfJET.eta();
			  recoToLHCO[i][1]=pfJET.phi();
			  recoToLHCO[i][2]=pfJET.pt();
			  recoToLHCO[i][3]=pfJET.mass();//4.7;
			  recoToLHCO[i][4]=2;
			  recoEnergy[i]=pfJET.energy();
			  if (pfJET.bDiscriminator(bdisc_name)<0)
			    {
			      negativeDiscriminant = true;
			    }
			  discriminants[i]=pfJET.bDiscriminator(bdisc_name);
			}
		    }
		  for (int i=2; i<4; i++)
		    {
		      //recoToLHCO[i][3]=0;
		      recoToLHCO[i][4]=0;
		    }
		}
	    }
	  
	  //////////////////////////////////////////////////////////Nueva
	  
	  cout<<"Paso4"<<endl;
	  //////////////////////////////
	  bool realGood = true;
	  for (int i=0; i<4; i++)
	    {
	      if ( realConfiguration[i]==-1 || negativeDiscriminant )//mirar bien que significa un discriminante negativo
		{
		  realGood=false;
		}
	    }
	  cout<<"Paso4a"<<endl;
	  if (realGood)
	    {
	      jets->Fill(5);
              cout<<"Paso4b"<<endl;
	      //making comparison between realConfiguration and combinations got from selection
	      bool Matched=false;
              bool MatchedByType=false;
              MatchedSelection->Fill(-10);
	      // solo para comparar seleccion con real
              double deltaR[4];
              for (int j=0; j<4;j++)
		{
		  cout<<"Paso4c"<<endl;
		  //prueba<<"real: "<<realConfiguration[j];
		  //para ver lo de los bjets en phi respecto al lepton
		  if (j<2)
		    {
		      //file11<<"diff en phi with lepton of reco "<<j<<"="<<recoToLHCO[j][1]-recoToLHCO[4][1]<<endl;
		    } 
		  //end para ver lo de los bjets en phi respecto al lepton
		  deltaR[j]= sqrt( pow(recoToLHCO[j][0]-genToLHCO[j][0],2.0 )+pow(recoToLHCO[j][1]-genToLHCO[j][1],2.0)+pow(recoToLHCO[j][2]-genToLHCO[j][2],2.0)+pow(recoToLHCO[j][3]-genToLHCO[j][3],2.0));   
		  cout<<"Paso4d"<<endl;  
		}
              //prueba<<endl;
              /*0404
              deltaR_b->Fill(deltaR[0]);
              deltaR_b_->Fill(deltaR[1]);
              deltaR_q->Fill(deltaR[2]);
              deltaR_q_->Fill(deltaR[3]);
              0404*/
	      cout<<"Paso5"<<endl;
	      // end solo para comparar seleccion con real
              for (int j=0; j<numberOfCombinations;j++)
		{
		  cout<<"Paso6"<<endl;
		  int bMatched=0;
                  int clMatched=0;
                  int matched=0;
                  bool endOfLine=false;
		  int l=1;
		  while (!endOfLine)
		    {
		      cout<<"Paso7"<<endl;
		      for (int i=0; i<4; i++)
			{
			  cout<<"Paso8"<<endl;
			  if (realConfiguration[i]==W_Combinations[j][l]-1)
			    {
			      cout<<"Paso9"<<endl;
			      matched++;
			      if (i>1)
				{
				  cout<<"Paso10"<<endl;
				  clMatched++;
				} 
			    } 
			} 
                      l++;
		      cout<<"Paso11"<<endl;
		      if (W_Combinations[j][l]==1000)
			{
			  cout<<"Paso12"<<endl;
			  endOfLine=true;
			}
		    }
		  
                  endOfLine=false;
		  l=1;
		  while (!endOfLine)
		    {
		      for (int i=0; i<4; i++)
			{
			  if (realConfiguration[i]==b1_Combinations[j][l]-1)
			    {
			      matched++;
			      if (i<2)
				{
				  bMatched++;
				} 
			    } 
			}
		      l++;
		      if (b1_Combinations[j][l]==1000)
			{
			  endOfLine=true;
			}
		    }
		  
		  endOfLine=false;
		  l=1;
		  while (!endOfLine)
		    {
		      for (int i=0; i<4; i++)
			{
			  if (realConfiguration[i]==b2_Combinations[j][l]-1)
			    {
			      matched++;
			      if (i<2)
				{
				  bMatched++;
				} 
			    } 
			}
		      l++;
		      if (b2_Combinations[j][l]==1000)
			{
			  endOfLine=true;
			}
		    }
		  cout<<"Paso30"<<endl;
                  if ((matched==4)&&(!Matched))
		    {
                      Matched=true;
		    }
                  if ((clMatched==2) &&(bMatched==2)&&(!MatchedByType))
		    {
                      MatchedByType=true;
		    }
                  if (Matched)
		    {
                      MatchedSelection->Fill((j+1)*10);
		    }
		  cout<<"Paso35"<<endl;
                  if (MatchedByType)
		    {
                      MatchedSelection->Fill((j+1)*10 + 1);
		    }
		  cout<<"Paso40"<<endl;
		}
	      //end making comparison between realConfiguration and combinations got from selection
	      
	      /////making plot of inv mass reco
              math::PtEtaPhiMLorentzVectorD jet3 (recoToLHCO[2][2],recoToLHCO[2][0],recoToLHCO[2][1],recoToLHCO[2][3]);
	      math::PtEtaPhiMLorentzVectorD jet4 (recoToLHCO[3][2],recoToLHCO[3][0],recoToLHCO[3][1],recoToLHCO[3][3]);
	      math::PtEtaPhiMLorentzVectorD jet1 (recoToLHCO[0][2],recoToLHCO[0][0],recoToLHCO[0][1],recoToLHCO[0][3]);
	      math::PtEtaPhiMLorentzVectorD W = jet3 + jet4;
	      math::PtEtaPhiMLorentzVectorD top = W + jet1;
	      math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[4][2],recoToLHCO[4][0],recoToLHCO[4][1],recoToLHCO[4][3]);
	      math::PtEtaPhiMLorentzVectorD met (recoToLHCO[5][2],recoToLHCO[5][0],recoToLHCO[5][1],recoToLHCO[5][3]);
	      math::PtEtaPhiMLorentzVectorD jet2 (recoToLHCO[1][2],recoToLHCO[1][0],recoToLHCO[1][1],recoToLHCO[1][3]);
	      math::PtEtaPhiMLorentzVectorD W_ = ele + met ;
	      math::PtEtaPhiMLorentzVectorD top_ = W_ + jet2;
	      Reco_Inv_Mass_top->Fill(top.M());
	      Reco_Inv_Mass_top_->Fill(top_.M());
	      Reco_Inv_Mass_b_->Fill(jet2.M());
	      Reco_Inv_Mass_b->Fill(jet1.M());
	      Reco_Inv_Mass_W_->Fill(W_.M());
	      Reco_Inv_Mass_W->Fill(W.M());

              double Anti_bTagging = fit_cl(discriminants[0])*fit_cl(discriminants[1])*fit_b(discriminants[2])*fit_b(discriminants[3]);
              anti_bTagging->Fill(Anti_bTagging);
              double DeltaR_bJetLep = sqrt(pow(jet2.Phi()-ele.Phi(),2)+pow(jet2.Eta()-ele.Eta(),2));
              deltaR_bJetLep->Fill(DeltaR_bJetLep);
              double Et_proportion = (jet1.energy() + jet2.energy() + jet3.energy() + jet4.energy())/Total_Energy_Jets;
	      et_proportion->Fill(Et_proportion);
              double Transverse_momentum = (jet1+jet2+jet3+jet4+met+ele).Pt();
	      transverse_momentum->Fill(Transverse_momentum);
              double Mt=sqrt(2*met.Pt()*ele.Pt()*(1-cos(met.Phi()-ele.Phi())));
              mt->Fill(Mt);

	      /*0404
              double discriminantJP = sqrt(pow(W.M()-80.385,2.0)+pow(jet1.M()-4.7,2.0)+pow(jet2.M()-4.7,2.0));
              discriminantJPMatching->Fill(discriminantJP); 
	      double elDiscriminante1=exp(-abs(top.M()-172.5))*abs(top_.M()-172.5); 
	      ElDiscriminante1 -> Fill(elDiscriminante1);
	      double elDiscriminante2=exp(-abs(W.M()-80.385))*abs(W_.M()-80.385); 
	      ElDiscriminante2 -> Fill(elDiscriminante2);
              0404*/
	      /////making plot of inv mass gen                                                                                                                                            
	      math::PtEtaPhiMLorentzVectorD gjet3 (genToLHCO[2][2],genToLHCO[2][0],genToLHCO[2][1],genToLHCO[2][3]);
	      math::PtEtaPhiMLorentzVectorD gjet4 (genToLHCO[3][2],genToLHCO[3][0],genToLHCO[3][1],genToLHCO[3][3]);
	      math::PtEtaPhiMLorentzVectorD gjet1 (genToLHCO[0][2],genToLHCO[0][0],genToLHCO[0][1],genToLHCO[0][3]);
	      math::PtEtaPhiMLorentzVectorD gW = gjet3 + gjet4;
	      math::PtEtaPhiMLorentzVectorD gtop = gW + gjet1;
	      math::PtEtaPhiMLorentzVectorD gele (genToLHCO[4][2],genToLHCO[4][0],genToLHCO[4][1],genToLHCO[4][3]);
	      math::PtEtaPhiMLorentzVectorD gmet (genToLHCO[5][2],genToLHCO[5][0],genToLHCO[5][1],genToLHCO[5][3]);
	      math::PtEtaPhiMLorentzVectorD gjet2 (genToLHCO[1][2],genToLHCO[1][0],genToLHCO[1][1],genToLHCO[1][3]);
	      math::PtEtaPhiMLorentzVectorD gW_ = gele + gmet ;
	      math::PtEtaPhiMLorentzVectorD gtop_ = gW_ + gjet2;
	      Gen_Inv_Mass_top->Fill(gtop.M());
	      Gen_Inv_Mass_top_->Fill(gtop_.M());
	      Gen_Inv_Mass_b_->Fill(gjet2.M());
	      Gen_Inv_Mass_b->Fill(gjet1.M());
	      Gen_Inv_Mass_W_->Fill(gW_.M());
	      Gen_Inv_Mass_W->Fill(gW.M());
	      
	      /*0404
	      double METJets =recoEnergy[5]/(recoEnergy[0]+recoEnergy[1]+recoEnergy[2]+recoEnergy[3]);
	      ProportionMETJets->Fill(METJets);
	      double gMETJets =genEnergy[5]/(genEnergy[0]+genEnergy[1]+genEnergy[2]+genEnergy[3]); 
	      gProportionMETJets->Fill(gMETJets);
              METvsNumberOfJets->Fill(recoEnergy[5],numberOfJets);
	      0404*/
	      
	      
	      
	      
	      
	      
	      //file10<<"JJJJJJJJJJJJJJJJJJJJ"<<endl;
	      //file7<<s_id<<" "<<discriminants[0]<<" "<<discriminants[1]<<" "<<discriminants[2]<<" "<<discriminants[3]<<" "<<InvMass<<" "<<InvMass_<<endl; 
              /*0404
              double deltaPhiJetLep_bLep=abs(recoToLHCO[1][1]-recoToLHCO[4][1]);
	      double deltaPhiJetq_bHad=abs(recoToLHCO[0][1]-recoToLHCO[2][1]);
	      double deltaPhiJetMET_bLep=abs(recoToLHCO[1][1]-recoToLHCO[5][1]);
	      double deltaPhiJetq__bHad=abs(recoToLHCO[0][1]-recoToLHCO[3][1]);
              if (deltaPhiJetLep_bLep>3.141592)
		{
		  deltaPhiJetLep_bLep=deltaPhiJetLep_bLep-3.1416;
		}
              if (deltaPhiJetMET_bLep>3.141592)
		{
		  deltaPhiJetMET_bLep=deltaPhiJetMET_bLep-3.1416;
		}
              if (deltaPhiJetq_bHad>3.141592)
		{
		  deltaPhiJetq_bHad=deltaPhiJetq_bHad-3.1416;
		}
              if (deltaPhiJetq__bHad>3.141592)
		{
		  deltaPhiJetq__bHad=deltaPhiJetq__bHad-3.1416;
		}

	      deltaPhiJetLep_vs_deltaPhiJetMET_bLep->Fill(deltaPhiJetLep_bLep, deltaPhiJetMET_bLep);
	      deltaPhiJetQ_vs_deltaPhiJetQ__bHad->Fill(deltaPhiJetq_bHad, deltaPhiJetq__bHad);
              
              deltaPhiJetLepbLep->Fill(deltaPhiJetLep_bLep);
              deltaPhiJetMETbLep->Fill(deltaPhiJetMET_bLep);
              deltaPhiJetLepMETbLep->Fill(deltaPhiJetLep_bLep + deltaPhiJetMET_bLep);	
	      
              deltaPhiJetQbHad->Fill(deltaPhiJetq_bHad);
              deltaPhiJetQ_bHad->Fill(deltaPhiJetq__bHad);
              deltaPhiJetQQ_bHad->Fill(deltaPhiJetq_bHad + deltaPhiJetq__bHad);
	
	      for (int i=0;i<4;i++)
		{
		  discriminant2_inclusive->Fill( discriminants[i] );  
		}
	      
	      discriminantb2_inclusive->Fill( discriminants[0] );
	      discriminantb_2_inclusive->Fill( discriminants[1] );
              0404*/
	      discriminantbTotal->Fill( discriminants[0] );
	      discriminantbTotal->Fill( discriminants[1] );
	      discriminantclTotal->Fill( discriminants[2] );
	      discriminantclTotal->Fill( discriminants[3] );
	      
	      //goodEvent=true;
	      for (int i=0;i<2;i++)
		{
		  //Difference2_InEnergy_bJets_inclusive_vs_PartonEnergy->Fill(genEnergy[i],recoEnergy[i]-genEnergy[i]);
		  Difference_InEnergy_bJets_Gen_Reco->Fill(recoEnergy[i]-genEnergy[i]);
                  Difference_InPt_bJets_Gen_Reco->Fill(recoToLHCO[i][2]-genToLHCO[i][2]);
                  Difference_InEta_bJets_Gen_Reco->Fill(recoToLHCO[i][0]-genToLHCO[i][0]);
                  Difference_InPhi_bJets_Gen_Reco->Fill(recoToLHCO[i][1]-genToLHCO[i][1]);
		  if (genEnergy[i]<=20)
		    { 		  
		      Difference_InEnergy_bJets_Gen_Reco_20->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>20)&&(genEnergy[i]<40))
		    { 		  
		      Difference_InEnergy_bJets_Gen_Reco_40->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>40)&&(genEnergy[i]<60))
		    { 		  
		      Difference_InEnergy_bJets_Gen_Reco_60->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>60)&&(genEnergy[i]<80))
		    { 		  
		      Difference_InEnergy_bJets_Gen_Reco_80->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>80)&&(genEnergy[i]<100))
		    { 		  
		      Difference_InEnergy_bJets_Gen_Reco_100->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>100)&&(genEnergy[i]<120))
		    { 		  
		      Difference_InEnergy_bJets_Gen_Reco_120->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>120)&&(genEnergy[i]<140))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_140->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>140)&&(genEnergy[i]<160))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_160->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>160)&&(genEnergy[i]<180))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_180->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>180)&&(genEnergy[i]<200))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_200->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>200)&&(genEnergy[i]<220))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_220->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>220)&&(genEnergy[i]<240))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_240->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>240)&&(genEnergy[i]<260))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_260->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>260)&&(genEnergy[i]<280))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_280->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>280)&&(genEnergy[i]<300))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_300->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>300)&&(genEnergy[i]<320))
		    {
		      Difference_InEnergy_bJets_Gen_Reco_320->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  
		  
		} 
	      
	      for (int i=2;i<4;i++)
		{
		  //Difference2_InEnergy_lJets_inclusive_vs_PartonEnergy->Fill(genEnergy[i],recoEnergy[i]-genEnergy[i]);
		  Difference_InEnergy_clJets_Gen_Reco->Fill(recoEnergy[i]-genEnergy[i]);
                  Difference_InPt_clJets_Gen_Reco->Fill(recoToLHCO[i][2]-genToLHCO[i][2]);
                  Difference_InEta_clJets_Gen_Reco->Fill(recoToLHCO[i][0]-genToLHCO[i][0]);
                  Difference_InPhi_clJets_Gen_Reco->Fill(recoToLHCO[i][1]-genToLHCO[i][1]);
		  if (genEnergy[i]<=20)
		    { 		  
		      Difference_InEnergy_clJets_Gen_Reco_20->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>20)&&(genEnergy[i]<40))
		    { 		  
		      Difference_InEnergy_clJets_Gen_Reco_40->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>40)&&(genEnergy[i]<60))
		    { 		  
		      Difference_InEnergy_clJets_Gen_Reco_60->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>60)&&(genEnergy[i]<80))
		    { 		  
		      Difference_InEnergy_clJets_Gen_Reco_80->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>80)&&(genEnergy[i]<100))
		    { 		  
		      Difference_InEnergy_clJets_Gen_Reco_100->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>100)&&(genEnergy[i]<120))
		    { 		  
		      Difference_InEnergy_clJets_Gen_Reco_120->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>120)&&(genEnergy[i]<140))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_140->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>140)&&(genEnergy[i]<160))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_160->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>160)&&(genEnergy[i]<180))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_180->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>180)&&(genEnergy[i]<200))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_200->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>200)&&(genEnergy[i]<220))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_220->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>220)&&(genEnergy[i]<240))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_240->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>240)&&(genEnergy[i]<260))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_260->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>260)&&(genEnergy[i]<280))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_280->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>280)&&(genEnergy[i]<300))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_300->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>300)&&(genEnergy[i]<320))
		    {
		      Difference_InEnergy_clJets_Gen_Reco_320->Fill(recoEnergy[i]-genEnergy[i]);
		    }
		  
		  
		}
	      Difference_InEta_MET_Gen_Reco->Fill(recoToLHCO[5][0]-genToLHCO[5][0]);
	      Difference_InPhi_MET_Gen_Reco->Fill(recoToLHCO[5][1]-genToLHCO[5][1]);
	      Difference_InPt_MET_Gen_Reco->Fill(recoToLHCO[5][2]-genToLHCO[5][2]);
	      Difference_InMass_MET_Gen_Reco->Fill(recoToLHCO[5][3]-genToLHCO[5][3]);
	      Difference_InEta_Muon_Gen_Reco->Fill(recoToLHCO[4][0]-genToLHCO[4][0]);
	      Difference_InPhi_Muon_Gen_Reco->Fill(recoToLHCO[4][1]-genToLHCO[4][1]);
	      Difference_InPt_Muon_Gen_Reco->Fill(recoToLHCO[4][2]-genToLHCO[4][2]);
	      Difference_InMass_Muon_Gen_Reco->Fill(recoToLHCO[4][3]-genToLHCO[4][3]);
	      ////////////////////////////////////////////////////////Begin parenthesis
	      
	      
	      const int s= 20; //this is the number of biggest jets to have into account                                                                                                                                       
	      double selectedJets_Pt[s];
	      double selectedJets[s];
	      int j;
	      for (int i=0; i<s;i++)
		{
		  selectedJets_Pt[i]=1000000;//initial values in the array
		  
		  selectedJets[i]=-1;                                                                                                                            
		}
	      for (int i=0; i<s;i++)
		{
		  double maxvalue = -1000000;//initial value bigger than any one else
		  int jetSelected=-1;
		  for (int k = 0; k != numberOfJets; ++k)
		    {
		      const pat::Jet& aJet = (*PFJet)[k];
		      if (i==0)
			{
			  j=0;
			}
		      else
			{
			  j=i-1;
			}
		      if ((aJet.pt()<selectedJets_Pt[j])&&(aJet.pt() > maxvalue))
			{
			  maxvalue=aJet.pt();
			  jetSelected=k;
			}
		    }
		  selectedJets_Pt[i]=maxvalue;
		  //prueba<<maxvalue<<endl;
		  selectedJets[i]=jetSelected;
		}
	      
	      
	      ////////////////////////////////////////////////////////              
	      ////////////////////////////////////////////////////////              
	      
	      
	      
	      bool in4=true;
	      bool in5=true;
	      bool in6=true;
	      bool in7=true;
	      bool in8=true;
	      bool in9=true;
	      bool in10=true;
	      bool in11=true;
	      bool in12=true;
	      bool in13=true;
	      bool in14=true;
	      bool in15=true;
	      
	      
	      for (int i=0;i<4;i++)
		{
		  
		  // prueba<<recoToLHCO[i][1]<<endl;
		  
		  if (recoToLHCO[i][2]<selectedJets_Pt[3])
		    {
		      in4=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[4])
		    {
		      in5=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[5])
		    {
		      in6=false;
		    }
		  
		  if (recoToLHCO[i][2]<selectedJets_Pt[6])
		    {
		      in7=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[7])
		    {
		      in8=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[8])
		    {
		      in9=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[9])
		    {
		      in10=false;
		    }
		  
		  if (recoToLHCO[i][2]<selectedJets_Pt[10])
		    {
		      in11=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[11])
		    {
		      in12=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[12])
		    {
		      in13=false;
		    }
		  if (recoToLHCO[i][2]<selectedJets_Pt[13])
		    {
		      in14=false;
		    }
		  
		  if (recoToLHCO[i][2]<selectedJets_Pt[14])
		    {
		      in15=false;
		    }
		  
		  /*0404
		  
		  if (recoToLHCO[i][0]>=0)
		    {
		      numberOfJetsInPhi->Fill(1);
		    }
		  else
		    {
		      numberOfJetsInPhi->Fill(3);
		    }
		    0404*/
		}
	      if(in4)
		{
		  in->Fill(4);
		}
	      if(in5)
		{
		  in->Fill(5);
		}
	      if(in6)
		{
		  in->Fill(6);
		}
	      if(in7)
		{
		  in->Fill(7);
		}
              if(in8)
		{
		  in->Fill(8);
		}
	      if(in9)
		{
		  in->Fill(9);
		}
	      if(in10)
		{
		  in->Fill(10);
		}
	      if(in11)
		{
		  in->Fill(11);
		}
              if(in12)
		{
		  in->Fill(12);
		}
	      if(in13)
		{
		  in->Fill(13);
		}
	      if(in14)
		{
		  in->Fill(14);
		}
	      if(in15)
		{
		  in->Fill(15);
		}
	      in->Fill(-1);
	      
	      ////////////////////////////////////////////////////////
	      
	      
	      double selectedJets_disc[s];
	      double disc;
	      for (int i=0; i<s;i++)
		{
		  selectedJets_disc[i]=1000000;//initial values in the array
		  
		  selectedJets[i]=-1;                                                                                                                            
		}
	      for (int i=0; i<s;i++)
		{
		  double maxvalue = -1000000;//initial value bigger than any one else
		  int jetSelected=-1;
		  for (int k = 0; k != numberOfJets; ++k)
		    {
		      const pat::Jet& aJet = (*PFJet)[k];
		      disc = aJet.bDiscriminator(bdisc_name);
		      if (i==0)
			{
			  j=0;
			}
		      else
			{
			  j=i-1;
			}
		      if ((disc < selectedJets_disc[j]) && (disc > maxvalue))
			{
			  maxvalue=disc;
			  jetSelected=k;
			}
		    }
		  selectedJets_disc[i]=maxvalue;
		  //prueba<<maxvalue<<endl;
		  selectedJets[i]=jetSelected;
		}
	      
	      
	      //////////////////////////////////////////////////////////////////////
	      ////////////////////////////////////////////////////////              
	      
	      
	      
	      bool inb2=true;
	      bool inb3=true;
	      bool inb4=true;
	      bool inb5=true;
	      for (int i=0;i<2;i++)
		{
		  
		  // prueba<<recoToLHCO[i][1]<<endl;
		  
		  if (discriminants[i]<selectedJets_disc[1])
		    {
		      inb2=false;
		    }
		  if (discriminants[i]<selectedJets_disc[2])
		    {
		      inb3=false;
		    }
		  if (discriminants[i]<selectedJets_disc[3])
		    {
		      inb4=false;
		    }
		  if (discriminants[i]<selectedJets_disc[4])
		    {
		      inb5=false;
		    }
		}
	      if(inb2)
		{
		  inb->Fill(2);
		}
	      if(inb3)
		{
		  inb->Fill(3);
		}
	      if(inb4)
		{
		  inb->Fill(4);
		}
	      if(inb5)
		{
		  inb->Fill(5);
		}
	      inb->Fill(-1);
	      
	      
	      
	      
	      /*0404
	      
	      ///selection of jets with biggest pt after selection with discriminant 
	      
	      
	      const int s1= 3; //this is the number of biggest jets to have into account                                                                                                                                       
	      //double selectedJets_Pt[s];
	      //int j;
	      for (int i=0; i<s1;i++)
		{
		  selectedJets_Pt[i]=1000000;//initial values in the array
		}
	      for (int i=3; i-3<s1; i++)
		{
		  selectedJets[i]=-1;                                                                                                                            
		}
	      for (int i=0; i<s1;i++)
		{
		  double maxvalue = -1000000;//initial value bigger than any one else
		  int jetSelected=-1;
		  bool notselectedyet = true;
		  for (int k = 0; k != numberOfJets; ++k)
		    {
		      const pat::Jet& aJet = (*PFJet)[k];
		      if (i==0)
			{
			  j=0;
			}
		      else
			{
			  j=i-1;
			}
		      for (int l=0; l<3; l++)
			{
			  if (selectedJets[l]==k)
			    {
			      notselectedyet = false;
			    }
			} 
		      if (notselectedyet && (aJet.pt()<selectedJets_Pt[j+3])&&(aJet.pt() > maxvalue))
			{
			  maxvalue=aJet.pt();
			  jetSelected=k;
			}
		    }
		  selectedJets_Pt[i]=maxvalue;
		  //prueba<<maxvalue<<endl;
		  selectedJets[i+3]=jetSelected;
		}
	      
	      ////////////////////////////////////////////////////////              
	      
	      
	      
	      double inSelection=0;
	      for (int i=0;i<4;i++)
		{
		  for (int j=0;j<6;j++)
		    {
		      if (realConfiguration[i]==selectedJets[j])
			{
			  inSelection++;
			}
		    }
		}
	      inSel->Fill(inSelection);
	      inSel->Fill(-1);
	      
	      
	      
	      0404*/
	      
	      
	      
	      ////////////////////////////////////////////////////////(End parenthesis)
	    
  


  if (MC && LHCOWithMC)
    {
      //jets->Fill(7);
      const int numberOfJets1 = 4; //later it could be 5 to include ISR    
      const int numberOfPermutations = 24; //later it can be improved as numberOfJets!                                                                                       
      file7<<s_id<<" "<<discriminants[0]<<" "<<discriminants[1]<<" "<<discriminants[2]<<" "<<discriminants[3]<<" "<<InvMass<<" "<<InvMass_<<endl; 
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
      
      
      
      for (int numberOfPermutation=0;numberOfPermutation<24;numberOfPermutation++)
	{
	  
	  ///////////Writing the recoLHCO file
	  file8<<"0"<<" "<<s_id<<numberOfPermutation<<" "<<"6"<<endl;
	  for (int j=0;j<4;j++)
	    {
	      file8<<j+1
		   <<" "
		   <<"4"
		   <<" "
		   <<recoToLHCO[permutationMap[numberOfPermutation][j]-1][0]
		   <<" "
		   <<recoToLHCO[permutationMap[numberOfPermutation][j]-1][1]
		   <<" "
		   <<recoToLHCO[permutationMap[numberOfPermutation][j]-1][2]
		   <<" "
		   <<recoToLHCO[j][3]
		   <<" "
		   <<"1"
		   <<" "
		   <<recoToLHCO[j][4]
		   <<" "
		   <<"0"
		   <<" "
		   <<"0"
		   <<" "
		   <<"0"
		   <<endl;
	    }
	  file8<<"5"
	       <<" "
	       <<"1"
	       <<" "
	       <<recoToLHCO[4][0]
	       <<" "
	       <<recoToLHCO[4][1]
	       <<" "
	       <<recoToLHCO[4][2]
	       <<" "
	       <<recoToLHCO[4][3]
	       <<" "
	       <<"-1"
	       <<" "
	       <<recoToLHCO[4][4]
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<endl;
	  file8<<"6"
	       <<" "
	       <<"6"
	       <<" "
	       <<recoToLHCO[5][0]
	       <<" "
	       <<recoToLHCO[5][1]
	       <<" "
	       <<recoToLHCO[5][2]
	       <<" "
	       <<recoToLHCO[5][3]
	       <<" "
	       <<"1"
	       <<" "
	       <<recoToLHCO[5][4]
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<endl;
	  ///////////
	  if (MC)
	    {
	      ///////////Writing the genLHCO file
	      file9<<"0"<<" "<<s_id<<numberOfPermutation<<" "<<"6"<<endl;
	      for (int j=0;j<4;j++)
		{
		  file9<<j+1
		       <<" "
		       <<"4"
		       <<" "
		       <<genToLHCO[permutationMap[numberOfPermutation][j]-1][0]
		       <<" "
		       <<genToLHCO[permutationMap[numberOfPermutation][j]-1][1]
		       <<" "
		       <<genToLHCO[permutationMap[numberOfPermutation][j]-1][2]
		       <<" "
		       <<genToLHCO[j][3]
		       <<" "
		       <<"1"
		       <<" "
		       <<genToLHCO[j][4]
		       <<" "
		       <<"0"
		       <<" "
		       <<"0"
		       <<" "
		       <<"0"
		       <<endl;
		}
	      file9<<"5"
		   <<" "
		   <<"1"
		   <<" "
		   <<genToLHCO[4][0]
		   <<" "
		   <<genToLHCO[4][1]
		   <<" "
		   <<genToLHCO[4][2]
		   <<" "
		   <<genToLHCO[4][3]
		   <<" "
		   <<"-1"
		   <<" "
		   <<genToLHCO[4][4]
		   <<" "
		   <<"0"
		   <<" "
		   <<"0"
		   <<" "
		   <<"0"
		   <<endl;
	      file9<<"6"
		   <<" "
		   <<"6"
		   <<" "
		   <<genToLHCO[5][0]
		   <<" "
		   <<genToLHCO[5][1]
		   <<" "
		   <<genToLHCO[5][2]
		   <<" "
		   <<genToLHCO[5][3]
		   <<" "
		   <<"1"
		   <<" "
		   <<genToLHCO[5][4]
		   <<" "
		   <<"0"
		   <<" "
		   <<"0"
		   <<" "
		   <<"0"
		   <<endl;
	      ///////////
	      
	      
	    }
	}
    }
}
	} 
    }
  
  

  
  // para mover cuando se quiera hacer seleccion de jets en caso de datos reales 05032013
  if (LHCOWithData)
    {
      bool lepton, MET, Jets;
      lepton = false;
      MET =  false;
      Jets = false;
      int leptonS,METS, JetsS;
      leptonS = 0;
      METS = 0;
      JetsS = 0;
      int leptonCharge;
      Handle<pat::METCollection> pfMet;
      iEvent.getByLabel("patMETsPF", pfMet);
      //if(pfMet.isValid())
      //{
      //  cout<<"FUNCIONO MET"<<endl;
      
      //MET = true;
      const pat::METCollection* PFMET = pfMet.product();
      METS = PFMET->size();
      const pat::MET& pfMET = (*PFMET)[0];
      recoToLHCO[5][0]=pfMET.eta();
      recoToLHCO[5][1]=pfMET.phi();
      recoToLHCO[5][2]=pfMET.pt();
      recoToLHCO[5][3]=pfMET.mass();//0;
      recoToLHCO[5][4]=0;
      recoEnergy[5]=pfMET.energy();
      //file<<"pfmetJJJPPP"<<pfMET.energy()<<endl;
      cout<<"paso1"<<endl;
      //}  
       	  
      Handle<pat::MuonCollection> pfMuon;
      iEvent.getByLabel("semilepMuonsPF", pfMuon);
      //if(pfMuon.isValid())
      //{
      //  cout<<"FUNCIONO MUON"<<endl;
      
      //lepton = true;
      const pat::MuonCollection* PFMuon = pfMuon.product();
      leptonS = PFMuon->size();
      if (leptonS ==1)
      {
      const pat::Muon& pfMUON = (*PFMuon)[0];
      recoToLHCO[4][0]=pfMUON.eta();
      recoToLHCO[4][1]=pfMUON.phi();
      recoToLHCO[4][2]=pfMUON.pt();
      recoToLHCO[4][3]=pfMUON.mass();//0;
      recoToLHCO[4][4]=0;
      recoEnergy[4]=pfMUON.energy();
      leptonCharge = pfMUON.charge();
      //file<<"pfMuonJJJPPP"<<pfMUON.energy()<<endl;
      cout<<"paso2"<<endl;
      }
      //}	 

      Handle<pat::ElectronCollection> pfElectron;
      iEvent.getByLabel("semilepElectronsPF", pfElectron);
      //if(pfElectron.isValid())
      //{
      //  cout<<"FUNCIONO ELECTRON"<<endl;
      
      //lepton = true;
      const pat::ElectronCollection* PFElectron = pfElectron.product();
      leptonS = leptonS + PFElectron->size();
      if ((leptonS == 1) && (PFElectron->size()==1))
      {
      const pat::Electron& pfELECTRON = (*PFElectron)[0];
      recoToLHCO[4][0]=pfELECTRON.eta();
      recoToLHCO[4][1]=pfELECTRON.phi();
      recoToLHCO[4][2]=pfELECTRON.pt();
      recoToLHCO[4][3]=pfELECTRON.mass();//0;
      recoToLHCO[4][4]=0;
      recoEnergy[4]=pfELECTRON.energy();
      leptonCharge = pfELECTRON.charge();
      //file<<"pfMuonJJJPPP"<<pfMUON.energy()<<endl;
      cout<<"paso3"<<endl;
      }
      if (leptonS>1)
      {
        cout<<"WrongJP: more than 1 lepton"<<endl;   
      }
      //}	  
      cout<<"paso4"<<endl;	  
      Handle<pat::JetCollection> pfJets;
      iEvent.getByLabel("selectedPatJetsPF", pfJets);
      //if(pfJets.isValid())
      //{
      //  cout<<"FUNCIONO Jets"<<endl;
      //Jets = true;
      const pat::JetCollection* PFJet = pfJets.product();
      JetsS = PFJet->size();
      //}
      cout<<"paso5"<<endl;
      cout<<"JP0404 MET="<<METS<<" leptons="<<leptonS<<" Jets="<<JetsS<<endl;
      if( (METS==1) && (leptonS==1) && (JetsS>3) ) //&& MET && lepton && Jets 
      {
        cout<<"FUNCIONO TODO"<<endl;

////2606
  bool InvMassTopWithStandardAnalysis = true;
  if (InvMassTopWithStandardAnalysis)
  {

	      ////////////////////////////////////////////////////////Begin parenthesis
	      
	      
	      const int s= 4; //this is the number of biggest jets to have into account                                                                                                                                       
	      double selectedJets_Pt[s];
	      double selectedJets[s];
	      int j;
	      for (int i=0; i<s;i++)
		{
		  selectedJets_Pt[i]=1000000;//initial values in the array
		  
		  selectedJets[i]=-1;                                                                                                                            
		}
	      for (int i=0; i<s;i++)
		{
		  double maxvalue = -1000000;//initial value bigger than any one else
		  int jetSelected=-1;
		  for (int k = 0; k != JetsS; ++k)
		    {
		      const pat::Jet& aJet = (*PFJet)[k];
		      if (i==0)
			{
			  j=0;
			}
		      else
			{
			  j=i-1;
			}
		      if ((aJet.pt()<selectedJets_Pt[j])&&(aJet.pt() > maxvalue))
			{
			  maxvalue=aJet.pt();
			  jetSelected=k;
			}
		    }
		  selectedJets_Pt[i]=maxvalue;
		  //prueba<<maxvalue<<endl;
		  selectedJets[i]=jetSelected;
		}
	      
	      
	      ////////////////////////////////////////////////////////              

              ////////////////////////////////////////////////////////
	      
	      
	      double selectedJets_disc[s];
              double selectedJetsd[s];
	      double disc;
	      for (int i=0; i<s;i++)
		{
		  selectedJets_disc[i]=1000000;//initial values in the array
		  
		  selectedJetsd[i]=-1;                                                                                                                            
		}
	      for (int i=0; i<s;i++)
		{
		  double maxvalue = -1000000;//initial value bigger than any one else
		  int jetSelected=-1;
		  for (int k = 0; k != 4; ++k)
		    {
		      const pat::Jet& aJet = (*PFJet)[selectedJets[k]];
		      disc = aJet.bDiscriminator(bdisc_name);
		      if (i==0)
			{
			  j=0;
			}
		      else
			{
			  j=i-1;
			}
		      if ((disc < selectedJets_disc[j]) && (disc > maxvalue))
			{
			  maxvalue=disc;
			  jetSelected=selectedJets[k];
			}
		    }
		  selectedJets_disc[i]=maxvalue;
		  //prueba<<maxvalue<<endl;
		  selectedJetsd[i]=jetSelected;
		}
	      
	      
	      //////////////////////////////////////////////////////////////////////
	      


              const pat::Jet& pfJETb1 = (*PFJet)[selectedJetsd[0]];
	      math::PtEtaPhiMLorentzVectorD jetb1 (pfJETb1.pt(),pfJETb1.eta(),pfJETb1.phi(),pfJETb1.mass());  
	      const pat::Jet& pfJETb2 = (*PFJet)[selectedJetsd[1]];
	      math::PtEtaPhiMLorentzVectorD jetb2 (pfJETb2.pt(),pfJETb2.eta(),pfJETb2.phi(),pfJETb2.mass());
              const pat::Jet& pfJETq1 = (*PFJet)[selectedJetsd[2]];
	      math::PtEtaPhiMLorentzVectorD jetq1 (pfJETq1.pt(),pfJETq1.eta(),pfJETq1.phi(),pfJETq1.mass());
              const pat::Jet& pfJETq2 = (*PFJet)[selectedJetsd[3]];
	      math::PtEtaPhiMLorentzVectorD jetq2 (pfJETq2.pt(),pfJETq2.eta(),pfJETq2.phi(),pfJETq2.mass());  
              math::PtEtaPhiMLorentzVectorD W = jetq1 + jetq2;
	      math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[4][2],recoToLHCO[4][0],recoToLHCO[4][1],recoToLHCO[4][3]);
	      math::PtEtaPhiMLorentzVectorD met (recoToLHCO[5][2],recoToLHCO[5][0],recoToLHCO[5][1],recoToLHCO[5][3]);
	      math::PtEtaPhiMLorentzVectorD W_ = ele + met ;
	      


	      math::PtEtaPhiMLorentzVectorD top_1 = W_ + jetb2;
	      math::PtEtaPhiMLorentzVectorD top1 = W + jetb1;
              math::PtEtaPhiMLorentzVectorD top_2 = W_ + jetb1;
	      math::PtEtaPhiMLorentzVectorD top2 = W + jetb2;

              
	      double deltaW_= (W_.M()-80.385);
	      double deltaW= (W.M()-80.385);
	      double deltaTop1 = (top1.M()-172.5);
              double deltaTop_1 = (top_1.M()-172.5);
	      double deltaTop2 = (top2.M()-172.5);
              double deltaTop_2 = (top_2.M()-172.5);


              math::PtEtaPhiMLorentzVectorD top_; 
	      math::PtEtaPhiMLorentzVectorD top;
              math::PtEtaPhiMLorentzVectorD b_; 
	      math::PtEtaPhiMLorentzVectorD b;
            
              if (deltaTop_1 + deltaTop1 < deltaTop_2 + deltaTop2)
              {
                 top_ = top_1;
                 top = top1;
                 b_ = jetb2;
                 b = jetb1;
              }
              else
              {
                 top_ = top_2;
                 top = top2;
                 b_ = jetb1;
                 b = jetb2;
              }

	      math::PtEtaPhiMLorentzVectorD t_b = top - b;
	      math::PtEtaPhiMLorentzVectorD t__b = top_ - b_;
 
              SAnal_Inv_Mass_W_Data0 -> Fill(W.M());        
	      SAnal_Inv_Mass_W__Data0 -> Fill(W_.M());        
              SAnal_Inv_Mass_t_Data0 -> Fill(top.M());        
	      SAnal_Inv_Mass_t__Data0 -> Fill(top_.M());        
	      SAnal_Inv_Mass_t_bData0 -> Fill(t_b.M());        
	      SAnal_Inv_Mass_t__bData0 -> Fill(t__b.M());        
			 

  }

////2606END



















      //const pat::JetCollection* PFJet = pfJets.product();
      bool negativeDiscriminant= false;
      const int numberOfCombinations = 100; 
      int numberOfJets = JetsS;
      double W_Combinations[numberOfCombinations][numberOfJets];//en realidad es numberOfJets-1 ya que se necesitan al menos 2 para b1 y b2 pero se gasta la primera casilla escribiendo el deltaWHad
      double b1_Combinations[numberOfCombinations][numberOfJets];//en realidad es numberOfJets-2 ya que se necesitan al menos 2 para el W y uno para b2 pero se gasta la primera casilla escribiendo el deltab1, pero por facilidad
      double b2_Combinations[numberOfCombinations][numberOfJets];//en realidad es numberOfJets-2 ya que se necesitan al menos 2 para el W y uno para b1 pero se gasta la primera casilla escribiendo el deltab2, pero por facilidad
      double deltaWHad, deltab1, deltab2;
      bool endOfComparison;
      for (int j=0; j<numberOfCombinations;j++)
	{
	  for (int k=0; k<numberOfJets;k++)
	    {
	      W_Combinations[j][k]=1000;
	      b1_Combinations[j][k]=1000;
	      b2_Combinations[j][k]=1000;
	    } 	
	}	
      int a, b, c, d, e, f;	  
      if (RFijo)
	{
	  a=1;
	  b=2;
	}
      else
	{
	  a=0;
	  b=numberOfJets-2;
	}

      //for (int jetNumber=0; jetNumber<numberOfJets-2; jetNumber++)  //R variable
      //for (int jetNumber=1;jetNumber<2;jetNumber++)  //R fijo 
      for (int jetNumber=a; jetNumber<b; jetNumber++)  
	{	
	  int N = numberOfJets;
	  int R1 = jetNumber;
	  const int N1 = R1+1;
	  //prueba<<"NW="<<N<<" CW="<<N1<<endl; 
	  cout<<"NW="<<N<<" CW="<<N1<<endl;
	  vector<int> combinationWHad(N1);
	  vector<int> othersWHad(N-N1);
	  for (int j=0; j<N1; j++)
	    {
	      combinationWHad[j]=j + 1;
	    }
	  for (int j=0; j<N-N1; j++)
	    {
	      othersWHad[j]=N1 + j + 1;
	    }
	      
	  lastCombination = false;
	  while (!lastCombination)
	    {
	      //                                                                                                                                                                                   
	      for (int j=0; j<N1;j++)
		{
		  //prueba<<"combWHad "<<combinationWHad[j]<<" ";
		  cout<<"combWHad "<<combinationWHad[j]<<" ";
		}
	      //prueba<<endl;
	      cout<<endl;
	      for(int j=0; j<N-N1;j++)
		{
		  //prueba<<"othersWHad "<<othersWHad[j]<<" ";
		  cout<<"othersWHad "<<othersWHad[j]<<" ";
		}
	      //prueba<<endl;
	      cout<<endl;
	      //                                                                                                                                                                                    
	      math::PtEtaPhiMLorentzVectorD WHad (0,0,0,0);
	      for (int j = 0; j < N1; j++)
		{
		  const pat::Jet& pfJETtmp = (*PFJet)[combinationWHad[j]-1];
		  if (j<2)
		    {
		      discriminants[j+2]=pfJETtmp.bDiscriminator(bdisc_name);
		    }

		  math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		  WHad = jetTmp + WHad;
		}
		  
	      //prueba<<"massWHad="<<WHad.M()<<endl;
	      cout<<"massWHad="<<WHad.M()<<endl;
	      deltaWHad=abs(WHad.M()-80.385);//73.29);
		  
	      //////
	      //cout<<"paso1"<<endl;
	      int M = N - N1;
	      if (RFijo)
		{
		  c=0;
		  d=1;
		}
	      else
		{
		  c=0;
		  d=M-1;
		}
	      //for (int l=0; l<M-1; l++) //R Variable
	      //for (int l=0;l<1;l++) //R fijo
	      for (int l=c; l<d; l++)    
		{
		  //cout<<"paso2"<<endl;
		  int R2 = l;
		  const int N2 = R2+1;
		  //prueba<<"Nb1="<<M<<" Cb1="<<N2<<endl;
		  cout<<"Nb1="<<M<<" Cb1="<<N2<<endl;
		  vector<int> combinationb1(N2);
		  vector<int> othersb1(M-N2);
		  for (int k=0; k<N2; k++)
		    {
		      combinationb1[k]=k + 1;
		    }
		  for (int k=0; k<M-N2; k++)
		    {
		      othersb1[k]=N2 + k + 1;
		    }
		      
		  lastCombination2= false;
		  while (!lastCombination2)
		    {
		      //cout<<"paso3"<<endl;
		      //                                                                                                                                                                    
		      for (int k=0; k<N2;k++)
			{
			  //prueba<<"combb1 "<<othersWHad[combinationb1[k]-1]<<" ";
			  cout<<"combb1 "<<othersWHad[combinationb1[k]-1]<<" ";
			}
		      //prueba<<endl;
		      cout<<endl;
		      for(int k=0; k<M-N2;k++)
			{
			  //prueba<<"othersb1 "<<othersWHad[othersb1[k]-1]<<" ";
			  cout<<"othersb1 "<<othersWHad[othersb1[k]-1]<<" ";
			}
		      //prueba<<endl;
		      cout<<endl;
		      //   
		      math::PtEtaPhiMLorentzVectorD b1 (0,0,0,0);
		      for (int k = 0; k < N2; k++)
			{
			  const pat::Jet& pfJETtmp = (*PFJet)[othersWHad[combinationb1[k]-1]-1];
			  if (k<1)
			    {
			      discriminants[k]=pfJETtmp.bDiscriminator(bdisc_name);
			    }

			  math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
			  b1 = jetTmp + b1;
			}
			  
		      //prueba<<"massb1="<<b1.M()<<endl;
		      cout<<"massb1="<<b1.M()<<endl;
			  
			  
		      deltab1=abs(b1.M()-4.7);//12.46);
                          
			  
		      //////
		      //cout<<"paso1a"<<endl;
		      int M1 = M - N2;
		      int tmp = M1 -1;  
		      //cout<<"M="<<M<<" N2="<<N2<<" M1="<<M1<<" M1-N3="<<tmp<<endl;
		      if (RFijo)
			{
			  e=0;
			  f=1;
			}
		      else
			{
			  e=0;
			  f=M1;
			}
		      //for (int l1=0; l1<M1; l1++)//R variable
		      //for (int l1=0;l1<1;l1++)//R fijo 
		      for (int l1=e; l1<f; l1++)  
			{
			  //cout<<"paso2a"<<endl;
			  int R3 = l1;
			  const int N3 = R3+1;
			  const int tmp1=M1-N3;
			  //cout<<"M1-N3="<<tmp1<<endl;
			  //prueba<<"Nb2="<<M1<<" Cb2="<<N3<<endl;
			  cout<<"Nb2="<<M1<<" Cb2="<<N3<<endl;
			  vector<int> combinationb2(N3);
			  cout<<"paso2b"<<endl;
			  vector<int> othersb2(M1-N3);
			  //cout<<"paso2c"<<endl;
			  for (int k1=0; k1<N3; k1++)
			    {
			      //cout<<"paso2d"<<endl;
			      combinationb2[k1]=k1 + 1;
			    }
			  //cout<<"paso2e"<<endl;
			  for (int k1=0; k1<M1-N3; k1++)
			    {
			      //cout<<"paso2f"<<endl;
			      othersb2[k1]=N3 + k1 + 1;
			    }
			  //cout<<"paso2g"<<endl;
			  lastCombination3= false;
			  while (!lastCombination3)
			    {
			      //cout<<"paso3a"<<endl;
			      //                                                                                                                                                       
			      for (int k1=0; k1<N3;k1++)
				{
				  //cout<<"paso3b"<<endl;
				  //prueba<<"combb2 "<<othersWHad[othersb1[combinationb2[k1]-1]-1]<<" ";
				  cout<<"combb2 "<<othersWHad[othersb1[combinationb2[k1]-1]-1]<<" ";
				}
			      //cout<<"paso3c"<<endl;
			      //prueba<<endl;
			      cout<<endl;
			      for(int k1=0; k1<M1-N3;k1++)
				{
				  //cout<<"paso3d"<<endl;
				  //prueba<<"othersb2 "<<othersWHad[othersb1[othersb2[k1]-1]-1]<<" ";
				  cout<<"othersb2 "<<othersWHad[othersb1[othersb2[k1]-1]-1]<<" ";
				}
			      //cout<<"paso3e"<<endl;
			      //prueba<<endl;
			      cout<<endl;
			      // 
			      math::PtEtaPhiMLorentzVectorD b2 (0,0,0,0);
			      //cout<<"paso3f"<<endl;
			      for (int k1 = 0; k1 < N3; k1++)
				{
				  //cout<<"paso3g"<<endl;
				  const pat::Jet& pfJETtmp = (*PFJet)[othersWHad[othersb1[combinationb2[k1]-1]-1]-1];
				  if (k1<1)
				    {
				      discriminants[k1+1]=pfJETtmp.bDiscriminator(bdisc_name);
				    }

				  math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
				  b2 = jetTmp + b2;
				}
			      double deltaRJetLep_bLep=pow(b2.eta()-recoToLHCO[4][0],2)+pow(b2.phi()-recoToLHCO[4][1],2);

			      double deltaPhiJetLep_bLep=abs(b2.phi()-recoToLHCO[4][1]);
			      if (deltaPhiJetLep_bLep>3.141592)
				{
				  deltaPhiJetLep_bLep=deltaPhiJetLep_bLep-3.1416;
				}
			      double deltaPhiJetMET_bLep=abs(b2.phi()-recoToLHCO[5][1]);
			      if (deltaPhiJetMET_bLep>3.141592)
				{
				  deltaPhiJetMET_bLep=deltaPhiJetMET_bLep-3.1416;
				}
			      double deltaPhi_bLep= deltaPhiJetLep_bLep + deltaPhiJetMET_bLep;


			      double deltaPhiJetLep_bHad_=abs(b1.phi()-recoToLHCO[4][1]);
			      if (deltaPhiJetLep_bHad_>3.141592)
				{
				  deltaPhiJetLep_bHad_=deltaPhiJetLep_bHad_-3.1416;
				}
			      double deltaPhiJetMET_bHad_=abs(b1.phi()-recoToLHCO[5][1]);
			      if (deltaPhiJetMET_bHad_>3.141592)
				{
				  deltaPhiJetMET_bHad_=deltaPhiJetMET_bHad_-3.1416;
				}
			      double deltaPhi_bHad_= deltaPhiJetLep_bHad_ + deltaPhiJetMET_bHad_;



			      double deltaPhiJetQ_bHad=abs(b1.phi()-recoToLHCO[2][1]);
			      if (deltaPhiJetQ_bHad>3.141592)
				{
				  deltaPhiJetQ_bHad=deltaPhiJetQ_bHad-3.1416;
				}
			      double deltaPhiJetQ__bHad=abs(b1.phi()-recoToLHCO[3][1]);
			      if (deltaPhiJetQ__bHad>3.141592)
				{
				  deltaPhiJetQ__bHad=deltaPhiJetQ__bHad-3.1416;
				}
			      double deltaPhi_bHad= deltaPhiJetQ_bHad + deltaPhiJetQ__bHad;

			      math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[4][2],recoToLHCO[4][0],recoToLHCO[4][1],recoToLHCO[4][3]);
			      math::PtEtaPhiMLorentzVectorD met (recoToLHCO[5][2],recoToLHCO[5][0],recoToLHCO[5][1],recoToLHCO[5][3]);
			      math::PtEtaPhiMLorentzVectorD WLep = ele + met ;
			      math::PtEtaPhiMLorentzVectorD tLep = WLep + b2;
			      math::PtEtaPhiMLorentzVectorD tHad = WHad + b1;
			      double deltatHad=abs(tHad.M()-172.5);
			      double deltatLep=abs(tLep.M()-172.5);

			      //prueba<<"massbLep="<<b2.M()<<endl;
			      cout<<"massbLep="<<b2.M()<<endl;
			      //cout<<"paso3h"<<endl;
			      deltab2=abs(b2.M()-4.7);//10.77);
                
			      double deltaPhi_WbHad= abs(b1.phi() -WHad.phi());
			      if (deltaPhi_WbHad>3.141592)
				{
				  deltaPhi_WbHad=deltaPhi_WbHad-3.1416;
				}

			      //deltaPhi_WbHad = deltaPhi_WbHad - 3.556;
			      //deltaPhi_bHad = deltaPhi_bHad - 2.648;
			      //deltaPhi_bLep = deltaPhi_bLep - 3.267;
			      //double discriminantJP = sqrt(pow(deltatHad,2.0)/30 + pow(deltaWHad,2.0)/11.75 + pow(deltab1,2.0)/3.36 + pow(deltab2,2.0)/4.1);
			      //double discriminantJP = (pow(deltaWHad,2.0)/11.76 + pow(deltab1,2.0)/3.36 + pow(deltab2,2.0)/4.15)*(pow(deltaPhi_bLep,2.0) +1.26*pow(deltaPhi_bHad,2.0));
			      //double discriminantJP = (( pow(deltaWHad,2.0)/10.29 + pow(deltab1,2.0)/2.587 + pow(deltab2,2.0)/3.926 )+100*(0.15*pow(deltaPhi_bLep,2.0)/0.6071 + pow(deltaPhi_bHad,2.0)/1.156 ))/100000;
			      //double discriminantJP = (pow(deltaPhi_bLep,2.0)/0.5943 + pow(deltaPhi_WbHad,2.0)/1.18);
			      //double discriminantJP = (pow(deltaWHad,2.0)/11.76 + pow(deltab1,2.0)/3.36 + pow(deltab2,2.0)/4.15)*(pow(deltaPhi_bLep,2.0) +1.26*pow(deltaPhi_WbHad,2.0));
			      //double discriminantJP = ( (pow(deltaWHad,2.0)/11.76 + pow(deltab1,2.0)/3.36 + pow(deltab2,2.0)/4.15 + 1)*(pow(deltaPhi_bLep,2.0) +1*pow(2*deltaPhi_WbHad,2.0) + 1) - 1 )/1000000;

			      double bTagDisc=fit_b(discriminants[0])*fit_b(discriminants[1])*fit_cl(discriminants[2])*fit_cl(discriminants[3]);
			      double AntibTagDisc=fit_cl(discriminants[0])*fit_cl(discriminants[1])*fit_b(discriminants[2])*fit_b(discriminants[3]);
			      double bTagDiscriminant=10000-fit_b(discriminants[0])*fit_b(discriminants[1])*fit_cl(discriminants[2])*fit_cl(discriminants[3]);
			      double discriminantJP;
			      if (selector==1)
				{
				  discriminantJP =(3*(1*(300*(bTagDiscriminant+ 1875*((3.241592)/(0.1+deltaPhi_bHad_))  + 18750*((deltaPhi_bHad)/(3.141592)) + 0.875*((deltaPhi_bLep)/(3.141592))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/1000000;//MATCHING EN 6                                                                                            
				}
			      if (selector==2)
				{
				  discriminantJP =(0.1*(1*(300*(10*bTagDiscriminant+ 187500000*((3.241592)/(0.1+deltaPhi_bHad_))  + 18750000000*((deltaPhi_bHad)/(3.141592)) + 87.5*((deltaPhi_bLep)/(3.141592))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/100000000000000000;
				}
			      if (selector==3)
				{
				  discriminantJP =(3*(1*(100*(bTagDiscriminant + 87.5000*((3.2416)/(0.1+deltaPhi_bHad_))  + 10*(1875000000*((deltaPhi_bHad)/(3.1416))+ 0.875*(deltaPhi_bLep/3.1416))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/1000000000000000000;
				}
			      if (selector==4)
				{
				  discriminantJP =100/(bTagDisc+100) +  3.2416/(0.1+deltaPhi_bHad_)  + (deltaPhi_bHad)/(3.1416) + (deltaPhi_bLep)/(3.1416) + pow(deltab2-4,2.0)/4.0 + pow(deltab1-4,2.0)/4.0 + pow(deltaWHad,2.0)/20;
				}
			      if (selector==5)
				{
				  discriminantJP = 4*100/(bTagDisc+100) + 2* 3.2416/(0.1+deltaPhi_bHad_)  +8*((deltaPhi_bHad)/(3.1416)+ (deltaPhi_bLep)/(3.1416)) + pow(deltab2-4,2.0)/4.0 + pow(deltab1-4,2.0)/4.0 + pow(deltaWHad,2.0)/20;
				}

			      if (selector==7)
				{
				  discriminantJP = (c1*c10/(bTagDisc+c10) + c2* c11/(c11+deltaPhi_bHad_)  +c3*(deltaPhi_bHad)/(3.1416)+ c4*(deltaPhi_bLep)/(3.1416) + c5*pow(deltab2-c9,2.0)/c13 +c6*pow(deltab1-c8,2.0)/c14 +c7*pow(deltaWHad,2.0)/20)/c12;
				}
			      if (selector==9)
				{
				  discriminantJP = pow(AntibTagDisc,c1)*pow(deltaRJetLep_bLep,c2)*pow(deltaWHad,2)/c3;
				}
			      if (selector==10)
				{
				  discriminantJP = pow(AntibTagDisc,c1)*pow(deltaRJetLep_bLep,c2)*(pow(deltaWHad,2)+pow(deltatHad,2))/c3;
				}


			      //double discriminantJP =(1*(1*(1000*(10*bTagDiscriminant + 0.875*((0.1+deltaPhi_bLep)/(0.1+deltaPhi_bHad_))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) +2*pow(deltaWHad,2.0))/100000;
			      //double discriminantJP =(3*(1*(300*(bTagDiscriminant + 1875*((3.241592)/(0.1+deltaPhi_bHad_))+ 18750*((deltaPhi_bHad)/(3.141592)) + 0.875*((deltaPhi_bLep)/(3.141592))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/1000000;
			      //				  double discriminantJP =(0.1*(1*(300*(10*bTagDiscriminant+ 187500000*((3.241592)/(0.1+deltaPhi_bHad_))  + 18750000000*((deltaPhi_bHad)/(3.141592)) + 87.5*((deltaPhi_bLep)/(3.141592))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/100000000000000000;
			      // double discriminantJP =(3*(1*(100*(bTagDiscriminant + 1875000*((3.2416)/(0.1+deltaPhi_bHad_))  + 1875000000*((deltaPhi_bHad)/(3.1416))+ 0.875*(deltaPhi_bLep/3.1416)) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/10000000000000000;
			      //double discriminantJP =(3*(1*(100*(bTagDiscriminant + 87.5000*((3.2416)/(0.1+deltaPhi_bHad_))  + 10*(1875000000*((deltaPhi_bHad)/(3.1416))+ 0.875*(deltaPhi_bLep/3.1416))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/1000000000000000000;

			      //				  double discriminantJP =(3*(1*(300*(bTagDiscriminant + 1875000*((3.2416)/(0.1+deltaPhi_bHad_))  + 1*(1875000000*((deltaPhi_bHad)/(3.1416))+ 0.875*(deltaPhi_bLep/3.1416))) + pow(deltab2,2.0)) + pow(deltab1,2.0)) + pow(deltaWHad,2.0))/1000000000000000000;
			      //double bTagDiscriminant=exp(-fit_b(discriminants[0])*fit_b(discriminants[1])*fit_cl(discriminants[2])*fit_cl(discriminants[3])/1000);
			      //double discriminantJP =(bTagDiscriminant + 0.5*(deltaPhi_bLep/deltaPhi_bHad_));



			      //double discriminantJP = (pow(deltaWHad,2.0)/11.76 + 10*1.26*pow(deltaPhi_bHad,2.0)*pow(deltab1,2.0)/3.36 + 10*pow(deltaPhi_bLep,2.0)*pow(deltab2,2.0)/4.15)/100000.0;
				  

			      //double discriminantJP = sqrt(pow(deltatHad,2.0)/61.22 + pow(deltatLep,2.0)/84.87 + pow(deltaWHad,2.0)/15.06 + pow(deltab1,2.0)/3.561 + pow(deltab2,2.0)/4.388);//*(pow(deltaPhi_bLep,1.0) + pow(deltaPhi_bHad,1.0)); 
			      //double discriminantJP = sqrt(pow(deltatHad,2.0)/17.55 + pow(deltatLep,2.0)/45.11 + pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316);//*(pow(deltaPhi_bLep,1.0) + pow(deltaPhi_bHad,1.0));
			      //double discriminantJP =  sqrt(pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)*(pow(deltaPhi_bLep,1.0) + pow(deltaPhi_bHad,1.0));

			      //double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)*(pow(deltaPhi_bLep,2.0)+pow(deltaPhi_bHad,2.0));
			      //double discriminantJP=(deltaWHad/13.24+1)*(deltaPhi_bHad/3.2*deltab1/4.175+1)*(deltaPhi_bLep/3.2*deltab2/4.316+1);
			      //				  double discriminantJP=(deltaWHad/13.24+1)*(deltab1/4.175+1)*(deltab2/4.316+1)*(deltaPhi_bLep/1.1+1)*(deltaPhi_bHad/1.1+1);
			      //				  double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)*(pow(deltaPhi_bLep-2.7,2.0)/1.2 + pow(deltaPhi_bHad-2.7,2.0)/1.2);
			      //				  double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)+(pow(deltaPhi_bLep,2.0)/1.2 + pow(deltaPhi_bHad,2.0)/1.2);                         	                          
			      //				  double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + (deltaPhi_bLep/3.2)*pow(deltab2,2.0)/4.316); 
			      //				  double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + pow(deltab1,2.0)/4.175 + pow(deltab2,2.0)/4.316)*(pow(deltaPhi_bLep,2.0) + pow(deltaPhi_bHad,2.0));
			      //double discriminantJP=(deltaWHad+1)*(deltab1+1)*(deltab2+1)*(deltaPhi_bLep+1);
			      //				  double discriminantJP =  (pow(deltaWHad,2.0)/13.24 + (deltaPhi_bHad/3.2)*pow(deltab1,2.0)/4.175 + (deltaPhi_bLep/3.2)*pow(deltab2,2.0)/4.316);
				 
			      ////
			      endOfComparison = false;
			      int j=0;
			      while (!endOfComparison)
				{
				  //cout<<"j="<<j<<endl;
				  //prueba<<"discriminantJP= "<<discriminantJP<<" WComb[j][0]= "<<W_Combinations[j][0]<<endl;
				  cout<<"discriminantJP= "<<discriminantJP<<" WComb[j][0]= "<<W_Combinations[j][0]<<" j="<<j<<endl;
				  if (discriminantJP<W_Combinations[j][0])
				    {
				      //prueba<<"Clasifica Combinacion"<<endl;
				      cout<<"Clasifica Combinacion"<<endl;
				      endOfComparison=true;
				      for (int k=numberOfCombinations-2;k>j-1;k--)
					{
					  //cout<<"k="<<k<<endl;
					  for (int l=0; l<numberOfJets;l++)
					    {
					      //cout<<"l="<<l<<endl;
					      W_Combinations[k+1][l]=W_Combinations[k][l];
					      b1_Combinations[k+1][l]=b1_Combinations[k][l];
					      b2_Combinations[k+1][l]=b2_Combinations[k][l];
					    }
					}
				      for (int l=1; l<N1+1;l++)
					{
					  //cout<<"lW="<<l<<endl;
					  W_Combinations[j][l]=combinationWHad[l-1];
					  //prueba<<" JPW "<<combinationWHad[l-1];
					}
				      //prueba<<endl;
				  
				      for (int l=1; l<N2+1;l++)
					{
					  //cout<<"lb1="<<l<<endl;
					  b1_Combinations[j][l]=othersWHad[combinationb1[l-1]-1];//combinationb1[l-1];
					  //prueba<<" JPb1 "<<combinationb1[l-1];
					}
				      //prueba<<endl;
					  
				      for (int l=1; l<N3+1;l++)
					{
					  //cout<<"lb2="<<l<<endl;
					  b2_Combinations[j][l]=othersWHad[othersb1[combinationb2[l-1]-1]-1];//combinationb2[l-1];
					  //prueba<<" JPb2 "<<combinationb2[l-1];
					}
				      //prueba<<endl;
				  
				      W_Combinations[j][0]=discriminantJP;
				      for(int l=N1+1; l<numberOfJets;l++)
					{
					  //cout<<"lWo="<<l<<endl;
					  W_Combinations[j][l]=1000;
					}
				  
				      for(int l=N2+1; l<numberOfJets;l++)
					{
					  //cout<<"lb1o="<<l<<endl;
					  b1_Combinations[j][l]=1000;
					}
				  
				      for(int l=N3+1; l<numberOfJets;l++)
					{
					  //cout<<"lb2o="<<l<<endl;
					  b2_Combinations[j][l]=1000;
					}
				  
				    }
				  if (j<99)
				    {
				      j++;
				    }
				  else
				    {
				      endOfComparison = true;
				    }
				}
			      ////
			    
			      cout<<"paso5a"<<endl;
			      lastCombination3 = nextComb(R3, combinationb2, othersb2, M1, R3);
			      cout<<"paso6a"<<lastCombination3<<endl;
			    }
			}
			  
		      //////
			  
			  
		      //cout<<"paso5"<<endl;
		      cout<<"lastCombBefore="<<lastCombination2<<" R2="<<R2<<endl;
		      lastCombination2 = nextComb(R2, combinationb1, othersb1, M, R2);
		      cout<<"lastCombAfter="<<lastCombination2<<" R2="<<R2<<endl;
		    }
		}
		  
	      lastCombination = nextComb(R1, combinationWHad, othersWHad, N, R1);  
		  
	    }
	}





      bool ThereIsACombination = false;
      double mass[3];
      double minDiscriminantJP=1000;
      for (int j=0; j<numberOfCombinations;j++)
	{
	  ThereIsACombination = false;   
	  //prueba<<"combinacion b1 j="<<j<<endl;
	  math::PtEtaPhiMLorentzVectorD b1 (0,0,0,0);
	  bool endOfLine=false;
	  int l=0;
	  while (!endOfLine)
	    {
	      //prueba<<b1_Combinations[j][l]<<" ";
	      if (l>0)
		{
		  const pat::Jet& pfJETtmp = (*PFJet)[b1_Combinations[j][l]-1];
		  math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		  b1 = jetTmp + b1;
		  if (l==1)
		    {
		      ThereIsACombination = true;
                      recoToLHCO[0][0]=pfJETtmp.eta();
		      recoToLHCO[0][1]=pfJETtmp.phi();
		      recoToLHCO[0][2]=pfJETtmp.pt();
		      recoToLHCO[0][3]=pfJETtmp.mass();//4.7;
		      recoToLHCO[0][4]=2;
		      recoEnergy[0]=pfJETtmp.energy();
		      discriminants[0]=pfJETtmp.bDiscriminator(bdisc_name);
		      file6<<"0"<<" "<<s_id<<j<<" "<<"6"<<endl;
		      file6<<1
			   <<" "
			   <<"4"
			   <<" "
			   <<recoToLHCO[0][0]
			   <<" "
			   <<recoToLHCO[0][1]
			   <<" "
			   <<recoToLHCO[0][2]
			   <<" "
			   <<recoToLHCO[0][3]
			   <<" "
			   <<"1"
			   <<" "
			   <<recoToLHCO[0][4]
			   <<" "
			   <<"0"
			   <<" "
			   <<"0"
			   <<" "
			   <<"0"
			   <<endl; 
		    }
		} 
	      l++;
	      if (b1_Combinations[j][l]==1000)
		{
		  endOfLine=true;
		}
	    }
	  //prueba<<endl;
	  //prueba<<"combinacion b2 j="<<j<<endl;
	  math::PtEtaPhiMLorentzVectorD b2 (0,0,0,0);
	  endOfLine=false;
	  l=0;
	  while (!endOfLine)
	    {
	      //prueba<<b2_Combinations[j][l]<<" ";
	      if (l>0)
		{
		  const pat::Jet& pfJETtmp = (*PFJet)[b2_Combinations[j][l]-1];
		  math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		  b2 = jetTmp + b2;
		  if (l==1)
		    {
		      recoToLHCO[1][0]=pfJETtmp.eta();
		      recoToLHCO[1][1]=pfJETtmp.phi();
		      recoToLHCO[1][2]=pfJETtmp.pt();
		      recoToLHCO[1][3]=pfJETtmp.mass();//4.7;
		      recoToLHCO[1][4]=2;
		      recoEnergy[1]=pfJETtmp.energy();
		      discriminants[1]=pfJETtmp.bDiscriminator(bdisc_name); 
		      file6<<2
			   <<" "
			   <<"4"
			   <<" "
			   <<recoToLHCO[1][0]
			   <<" "
			   <<recoToLHCO[1][1]
			   <<" "
			   <<recoToLHCO[1][2]
			   <<" "
			   <<recoToLHCO[1][3]
			   <<" "
			   <<"1"
			   <<" "
			   <<recoToLHCO[1][4]
			   <<" "
			   <<"0"
			   <<" "
			   <<"0"
			   <<" "
			   <<"0"
			   <<endl; 
			      
			      
		    }
		} 
	      l++;
	      if (b2_Combinations[j][l]==1000)
		{
		  endOfLine=true;
		}
	    }
	  //prueba<<endl;
	  endOfLine=false;
	  l=0;
	  //prueba<<"combinacion W j="<<j<<endl;
	  math::PtEtaPhiMLorentzVectorD WHad (0,0,0,0);
	  while (!endOfLine)
	    {  
	      //prueba<<W_Combinations[j][l]<<" ";
	      if (l>0)
		{
		  const pat::Jet& pfJETtmp = (*PFJet)[W_Combinations[j][l]-1];
		  math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
		  WHad = jetTmp + WHad;
		  if (l<3)
		    {
		      recoToLHCO[l+1][0]=pfJETtmp.eta();
		      recoToLHCO[l+1][1]=pfJETtmp.phi();
		      recoToLHCO[l+1][2]=pfJETtmp.pt();
		      recoToLHCO[l+1][3]=pfJETtmp.mass();//0;
		      recoToLHCO[l+1][4]=0;
		      recoEnergy[l+1]=pfJETtmp.energy();
		      discriminants[l+1]=pfJETtmp.bDiscriminator(bdisc_name);
			  
		      file6<<l+2
			   <<" "
			   <<"4"
			   <<" "
			   <<recoToLHCO[l+1][0]
			   <<" "
			   <<recoToLHCO[l+1][1]
			   <<" "
			   <<recoToLHCO[l+1][2]
			   <<" "
			   <<recoToLHCO[l+1][3]
			   <<" "
			   <<"1"
			   <<" "
			   <<recoToLHCO[l+1][4]
			   <<" "
			   <<"0"
			   <<" "
			   <<"0"
			   <<" "
			   <<"0"
			   <<endl; 
		    }
		} 
	      l++;
	      if (W_Combinations[j][l]==1000)
		{
		  endOfLine=true;
		}
	    }
	  //prueba<<endl;
	      
	  //double discriminantJP = (abs(WHad.M()-80.385) + 1000)*(abs(bHad.M()-4.7) + 1000)*(abs(bLep.M()-4.7) + 1000);
	  //double discriminantJP = sqrt(pow(WHad.M()-73.29,2.0)+pow(bHad.M()-12.46,2.0)+pow(bLep.M()-10.77,2.0));
	  double discriminantJP = sqrt(pow(WHad.M()-80.385,2.0)+pow(b1.M()-4.7,2.0) +pow(b2.M()-4.7,2.0));
	  /*0404
	  if (j==0)
	    {
	      discriminantJPSelection_Data->Fill(discriminantJP);
	    }
	  discriminantJPvsN_Data->Fill(j,discriminantJP);
	  if (discriminantJP < minDiscriminantJP)// es lo mismo que decir j=0
	    {
	      minDiscriminantJP=discriminantJP;
	      mass[0]=WHad.M();
	      mass[1]=b1.M();
	      mass[2]=b2.M();
	    }
          0404*/
          if (ThereIsACombination)
          {
	  file6<<"5"
	       <<" "
	       <<"1"
	       <<" "
	       <<recoToLHCO[4][0]
	       <<" "
	       <<recoToLHCO[4][1]
	       <<" "
	       <<recoToLHCO[4][2]
	       <<" "
	       <<recoToLHCO[4][3]
	       <<" "
	       <<leptonCharge
	       <<" "
	       <<recoToLHCO[4][4]
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<endl; 
	      
	  file6<<"6"
	       <<" "
	       <<"6"
	       <<" "
	       <<recoToLHCO[5][0]
	       <<" "
	       <<recoToLHCO[5][1]
	       <<" "
	       <<recoToLHCO[5][2]
	       <<" "
	       <<recoToLHCO[5][3]
	       <<" "
	       <<"1"
	       <<" "
	       <<recoToLHCO[5][4]
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<" "
	       <<"0"
	       <<endl;
	      
	  file5<<s_id<<j<<" "<<discriminants[0]<<" "<<discriminants[1]<<" "<<discriminants[2]<<" "<<discriminants[3]<<endl;
	      
	  if (j<1)
	    { 
	      /////making plot of inv mass reco
	      math::PtEtaPhiMLorentzVectorD jet3 (recoToLHCO[2][2],recoToLHCO[2][0],recoToLHCO[2][1],recoToLHCO[2][3]);
	      math::PtEtaPhiMLorentzVectorD jet4 (recoToLHCO[3][2],recoToLHCO[3][0],recoToLHCO[3][1],recoToLHCO[3][3]);
	      math::PtEtaPhiMLorentzVectorD jet1 (recoToLHCO[0][2],recoToLHCO[0][0],recoToLHCO[0][1],recoToLHCO[0][3]);
	      math::PtEtaPhiMLorentzVectorD W = jet3 + jet4;
	      math::PtEtaPhiMLorentzVectorD top = W + jet1;
	      math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[4][2],recoToLHCO[4][0],recoToLHCO[4][1],recoToLHCO[4][3]);
	      math::PtEtaPhiMLorentzVectorD met (recoToLHCO[5][2],recoToLHCO[5][0],recoToLHCO[5][1],recoToLHCO[5][3]);
	      math::PtEtaPhiMLorentzVectorD jet2 (recoToLHCO[1][2],recoToLHCO[1][0],recoToLHCO[1][1],recoToLHCO[1][3]);
	      math::PtEtaPhiMLorentzVectorD W_ = ele + met ;
	      math::PtEtaPhiMLorentzVectorD top_ = W_ + jet2;
	      /*
		double deltaTop = abs(top.M()-172.5);
		double deltaTop_ = abs(top_.M()-172.5);
		double deltaW = abs(W.M()-80.385);
		double deltaW_= abs(W_.M()-80.385);
		if (deltaTop<45)
		{
		deltaTop=0;
		} 
		if (deltaTop_<40)
		{
		deltaTop=0;
		}
		if (deltaW<25)
		{
		deltaW=0;
		} 
		if (deltaW_<35)
		{
		deltaW_=0;
		}
	      */
	      double deltaW_= (W_.M()-80.385);
	      double deltaTop_ = (top_.M()-172.5);
	      double deltaW= (W.M()-80.385);
	      double deltaTop = (top.M()-172.5);
	      if (j==0)
		{
		  /*
		  //
		  double discriminant_lep_W_10_10=pow(((deltaW_+1)/3.5),50)/(pow(((deltaW_+1)/3.5),50)+10);
		  Discriminant_lep_W_10_10 -> Fill(discriminant_lep_W_10_10);
		  double discriminant_lep_t_10_10=pow(((deltaTop_+1)/3.6),50)/(pow(((deltaTop_+1)/3.6),50)+10);
		  Discriminant_lep_t_10_10 -> Fill(discriminant_lep_t_10_10);

		  double discriminant_had_W_10_10=10/(pow(((deltaW+2.5)/1.7),50)+10);
		  Discriminant_had_W_10_10 -> Fill(discriminant_had_W_10_10);
		  double discriminant_had_t_10_10=10/(pow(((deltaTop-12)/2.6),50)+10);
		  Discriminant_had_t_10_10 -> Fill(discriminant_had_t_10_10);

		  double discriminant_total_W_10_10=discriminant_lep_W_10_10 * discriminant_had_W_10_10;
		  Discriminant_total_W_10_10 -> Fill(discriminant_total_W_10_10);
		  double discriminant_total_t_10_10=discriminant_lep_t_10_10 * discriminant_had_t_10_10;
		  Discriminant_total_t_10_10 -> Fill(discriminant_total_t_10_10);
		  //
		  //
		  double discriminant_lep_W_30_30=pow(((deltaW_+1)/10.5),50)/(pow(((deltaW_+1)/10.5),50)+10);
		  Discriminant_lep_W_30_30 -> Fill(discriminant_lep_W_30_30);
		  double discriminant_lep_t_30_30=pow(((deltaTop_+1)/10.8),50)/(pow(((deltaTop_+1)/10.8),50)+10);
		  Discriminant_lep_t_30_30 -> Fill(discriminant_lep_t_30_30);

		  double discriminant_had_W_30_30=10/(pow(((deltaW+2.5)/5.1),50)+10);
		  Discriminant_had_W_30_30 -> Fill(discriminant_had_W_30_30);
		  double discriminant_had_t_30_30=10/(pow(((deltaTop-12)/7.8),50)+10);
		  Discriminant_had_t_30_30 -> Fill(discriminant_had_t_30_30);

		  double discriminant_total_W_30_30=discriminant_lep_W_30_30 * discriminant_had_W_30_30;
		  Discriminant_total_W_30_30 -> Fill(discriminant_total_W_30_30);
		  double discriminant_total_t_30_30=discriminant_lep_t_30_30 * discriminant_had_t_30_30;
		  Discriminant_total_t_30_30 -> Fill(discriminant_total_t_30_30);
		  //		      
		  //
		  double discriminant_lep_W_50_50=pow(((deltaW_+1)/17.5),50)/(pow(((deltaW_+1)/17.5),50)+10);
		  Discriminant_lep_W_50_50 -> Fill(discriminant_lep_W_50_50);
		  double discriminant_lep_t_50_50=pow(((deltaTop_+1)/18),50)/(pow(((deltaTop_+1)/18),50)+10);
		  Discriminant_lep_t_50_50 -> Fill(discriminant_lep_t_50_50);

		  double discriminant_had_W_50_50=10/(pow(((deltaW+2.5)/8.5),50)+10);
		  Discriminant_had_W_50_50 -> Fill(discriminant_had_W_50_50);
		  double discriminant_had_t_50_50=10/(pow(((deltaTop-12)/13),50)+10);
		  Discriminant_had_t_50_50 -> Fill(discriminant_had_t_50_50);

		  double discriminant_total_W_50_50=discriminant_lep_W_50_50 * discriminant_had_W_50_50;
		  Discriminant_total_W_50_50 -> Fill(discriminant_total_W_50_50);
		  double discriminant_total_t_50_50=discriminant_lep_t_50_50 * discriminant_had_t_50_50;
		  Discriminant_total_t_50_50 -> Fill(discriminant_total_t_50_50);
		  //		      
		  //
		  double discriminant_lep_W_70_70=pow(((deltaW_+1)/24.5),50)/(pow(((deltaW_+1)/24.5),50)+10);
		  Discriminant_lep_W_70_70 -> Fill(discriminant_lep_W_70_70);
		  double discriminant_lep_t_70_70=pow(((deltaTop_+1)/25.2),50)/(pow(((deltaTop_+1)/25.2),50)+10);
		  Discriminant_lep_t_70_70 -> Fill(discriminant_lep_t_70_70);

		  double discriminant_had_W_70_70=10/(pow(((deltaW+2.5)/11.9),50)+10);
		  Discriminant_had_W_70_70 -> Fill(discriminant_had_W_70_70);
		  double discriminant_had_t_70_70=10/(pow(((deltaTop-12)/18.2),50)+10);
		  Discriminant_had_t_70_70 -> Fill(discriminant_had_t_70_70);

		  double discriminant_total_W_70_70=discriminant_lep_W_70_70 * discriminant_had_W_70_70;
		  Discriminant_total_W_70_70 -> Fill(discriminant_total_W_70_70);
		  double discriminant_total_t_70_70=discriminant_lep_t_70_70 * discriminant_had_t_70_70;
		  Discriminant_total_t_70_70 -> Fill(discriminant_total_t_70_70);
		  //		      
		  //
		  double discriminant_lep_W_90_90=pow(((deltaW_+1)/42),50)/(pow(((deltaW_+1)/42),50)+10);
		  Discriminant_lep_W_90_90 -> Fill(discriminant_lep_W_90_90);
		  double discriminant_lep_t_90_90=pow(((deltaTop_+1)/43.2),50)/(pow(((deltaTop_+1)/43.2),50)+10);
		  Discriminant_lep_t_90_90 -> Fill(discriminant_lep_t_90_90);

		  double discriminant_had_W_90_90=10/(pow(((deltaW+2.5)/20.4),50)+10);
		  Discriminant_had_W_90_90 -> Fill(discriminant_had_W_90_90);
		  double discriminant_had_t_90_90=10/(pow(((deltaTop-12)/31.5),50)+10);
		  Discriminant_had_t_90_90 -> Fill(discriminant_had_t_90_90);

		  double discriminant_total_W_90_90=discriminant_lep_W_90_90 * discriminant_had_W_90_90;
		  Discriminant_total_W_90_90 -> Fill(discriminant_total_W_90_90);
		  double discriminant_total_t_90_90=discriminant_lep_t_90_90 * discriminant_had_t_90_90;
		  Discriminant_total_t_90_90 -> Fill(discriminant_total_t_90_90);
		  //		      
		  //
		  double discriminant_lep_W_10_90=pow(((deltaW_+1)/3.5),50)/(pow(((deltaW_+1)/3.5),50)+10);
		  Discriminant_lep_W_10_90 -> Fill(discriminant_lep_W_10_90);
		  double discriminant_lep_t_10_90=pow(((deltaTop_+1)/3.6),50)/(pow(((deltaTop_+1)/3.6),50)+10);
		  Discriminant_lep_t_10_90 -> Fill(discriminant_lep_t_10_90);

		  double discriminant_had_W_10_90=10/(pow(((deltaW+2.5)/42),50)+10);
		  Discriminant_had_W_10_90 -> Fill(discriminant_had_W_10_90);
		  double discriminant_had_t_10_90=10/(pow(((deltaTop-12)/43.2),50)+10);
		  Discriminant_had_t_10_90 -> Fill(discriminant_had_t_10_90);

		  double discriminant_total_W_10_90=discriminant_lep_W_10_90 * discriminant_had_W_10_90;
		  Discriminant_total_W_10_90 -> Fill(discriminant_total_W_10_90);
		  double discriminant_total_t_10_90=discriminant_lep_t_10_90 * discriminant_had_t_10_90;
		  Discriminant_total_t_10_90 -> Fill(discriminant_total_t_10_90);
		  //		      
		  //
		  double discriminant_lep_W_30_70=pow(((deltaW_+1)/10.5),50)/(pow(((deltaW_+1)/10.5),50)+10);
		  Discriminant_lep_W_30_70 -> Fill(discriminant_lep_W_30_70);
		  double discriminant_lep_t_30_70=pow(((deltaTop_+1)/10.8),50)/(pow(((deltaTop_+1)/10.8),50)+10);
		  Discriminant_lep_t_30_70 -> Fill(discriminant_lep_t_30_70);

		  double discriminant_had_W_30_70=10/(pow(((deltaW+2.5)/11.9),50)+10);
		  Discriminant_had_W_30_70 -> Fill(discriminant_had_W_30_70);
		  double discriminant_had_t_30_70=10/(pow(((deltaTop-12)/18.2),50)+10);
		  Discriminant_had_t_30_70 -> Fill(discriminant_had_t_30_70);

		  double discriminant_total_W_30_70=discriminant_lep_W_30_70 * discriminant_had_W_30_70;
		  Discriminant_total_W_30_70 -> Fill(discriminant_total_W_30_70);
		  double discriminant_total_t_30_70=discriminant_lep_t_30_70 * discriminant_had_t_30_70;
		  Discriminant_total_t_30_70 -> Fill(discriminant_total_t_30_70);
		  //		      
		  //
		  double discriminant_lep_W_70_30=pow(((deltaW_+1)/24.5),50)/(pow(((deltaW_+1)/24.5),50)+10);
		  Discriminant_lep_W_70_30 -> Fill(discriminant_lep_W_70_30);
		  double discriminant_lep_t_70_30=pow(((deltaTop_+1)/25.2),50)/(pow(((deltaTop_+1)/25.2),50)+10);
		  Discriminant_lep_t_70_30 -> Fill(discriminant_lep_t_70_30);

		  double discriminant_had_W_70_30=10/(pow(((deltaW+2.5)/5.1),50)+10);
		  Discriminant_had_W_70_30 -> Fill(discriminant_had_W_70_30);
		  double discriminant_had_t_70_30=10/(pow(((deltaTop-12)/7.8),50)+10);
		  Discriminant_had_t_70_30 -> Fill(discriminant_had_t_70_30);

		  double discriminant_total_W_70_30=discriminant_lep_W_70_30 * discriminant_had_W_70_30;
		  Discriminant_total_W_70_30 -> Fill(discriminant_total_W_70_30);
		  double discriminant_total_t_70_30=discriminant_lep_t_70_30 * discriminant_had_t_70_30;
		  Discriminant_total_t_70_30 -> Fill(discriminant_total_t_70_30);
		  //		      
		  //
		  double discriminant_lep_W_90_10=pow(((deltaW_+1)/42),50)/(pow(((deltaW_+1)/42),50)+10);
		  Discriminant_lep_W_90_10 -> Fill(discriminant_lep_W_90_10);
		  double discriminant_lep_t_90_10=pow(((deltaTop_+1)/43.2),50)/(pow(((deltaTop_+1)/43.2),50)+10);
		  Discriminant_lep_t_90_10 -> Fill(discriminant_lep_t_90_10);

		  double discriminant_had_W_90_10=10/(pow(((deltaW+2.5)/1.7),50)+10);
		  Discriminant_had_W_90_10 -> Fill(discriminant_had_W_90_10);
		  double discriminant_had_t_90_10=10/(pow(((deltaTop-12)/2.6),50)+10);
		  Discriminant_had_t_90_10 -> Fill(discriminant_had_t_90_10);

		  double discriminant_total_W_90_10=discriminant_lep_W_90_10 * discriminant_had_W_90_10;
		  Discriminant_total_W_90_10 -> Fill(discriminant_total_W_90_10);
		  double discriminant_total_t_90_10=discriminant_lep_t_90_10 * discriminant_had_t_90_10;
		  Discriminant_total_t_90_10 -> Fill(discriminant_total_t_90_10);
		  //		      
		  */

                  /*0404 
		  //
		  double discriminant_lep_W_90_90=pow(((deltaW_+0.5)/36),50)/(pow(((deltaW_+0.5)/36),50)+10);
		  Discriminant_lep_W_90_90 -> Fill(discriminant_lep_W_90_90);
		  double discriminant_lep_t_90_90=pow(((deltaTop_-5)/39.5),50)/(pow(((deltaTop_-5)/39.5),50)+10);
		  Discriminant_lep_t_90_90 -> Fill(discriminant_lep_t_90_90);

		  double discriminant_had_W_90_90=10/(pow(((deltaW-6)/15.8),50)+10);
		  Discriminant_had_W_90_90 -> Fill(discriminant_had_W_90_90);
		  double discriminant_had_t_90_90=10/(pow(((deltaTop-16)/25),50)+10);
		  Discriminant_had_t_90_90 -> Fill(discriminant_had_t_90_90);

		  double discriminant_total_W_90_90=discriminant_lep_W_90_90 * discriminant_had_W_90_90;
		  Discriminant_total_W_90_90 -> Fill(discriminant_total_W_90_90);
		  double discriminant_total_t_90_90=discriminant_lep_t_90_90 * discriminant_had_t_90_90;
		  Discriminant_total_t_90_90 -> Fill(discriminant_total_t_90_90);
		  //


		  //
		  double discriminant_lep_W_85_85=pow(((deltaW_+0.5)/30.5),50)/(pow(((deltaW_+0.5)/30.5),50)+10);
		  Discriminant_lep_W_85_85 -> Fill(discriminant_lep_W_85_85);
		  double discriminant_lep_t_85_85=pow(((deltaTop_-5)/32.9),50)/(pow(((deltaTop_-5)/32.9),50)+10);
		  Discriminant_lep_t_85_85 -> Fill(discriminant_lep_t_85_85);

		  double discriminant_had_W_85_85=10/(pow(((deltaW-6)/13.25),50)+10);
		  Discriminant_had_W_85_85 -> Fill(discriminant_had_W_85_85);
		  double discriminant_had_t_85_85=10/(pow(((deltaTop-16)/20.7),50)+10);
		  Discriminant_had_t_85_85 -> Fill(discriminant_had_t_85_85);

		  double discriminant_total_W_85_85=discriminant_lep_W_85_85 * discriminant_had_W_85_85;
		  Discriminant_total_W_85_85 -> Fill(discriminant_total_W_85_85);
		  double discriminant_total_t_85_85=discriminant_lep_t_85_85 * discriminant_had_t_85_85;
		  Discriminant_total_t_85_85 -> Fill(discriminant_total_t_85_85);
		  //

		  //
		  double discriminant_lep_W_95_95=pow(((deltaW_+0.5)/46),50)/(pow(((deltaW_+0.5)/46),50)+10);
		  Discriminant_lep_W_95_95 -> Fill(discriminant_lep_W_95_95);
		  double discriminant_lep_t_95_95=pow(((deltaTop_-5)/49.5),50)/(pow(((deltaTop_-5)/49.5),50)+10);
		  Discriminant_lep_t_95_95 -> Fill(discriminant_lep_t_95_95);

		  double discriminant_had_W_95_95=10/(pow(((deltaW-6)/20),50)+10);
		  Discriminant_had_W_95_95 -> Fill(discriminant_had_W_95_95);
		  double discriminant_had_t_95_95=10/(pow(((deltaTop-16)/31),50)+10);
		  Discriminant_had_t_95_95 -> Fill(discriminant_had_t_95_95);

		  double discriminant_total_W_95_95=discriminant_lep_W_95_95 * discriminant_had_W_95_95;
		  Discriminant_total_W_95_95 -> Fill(discriminant_total_W_95_95);
		  double discriminant_total_t_95_95=discriminant_lep_t_95_95 * discriminant_had_t_95_95;
		  Discriminant_total_t_95_95 -> Fill(discriminant_total_t_95_95);
		  //
		  //
		  double discriminant_lep_W_85_95=pow(((deltaW_+0.5)/30.5),50)/(pow(((deltaW_+0.5)/30.5),50)+10);
		  Discriminant_lep_W_85_95 -> Fill(discriminant_lep_W_85_95);
		  double discriminant_lep_t_85_95=pow(((deltaTop_-5)/32.9),50)/(pow(((deltaTop_-5)/32.9),50)+10);
		  Discriminant_lep_t_85_95 -> Fill(discriminant_lep_t_85_95);

		  double discriminant_had_W_85_95=10/(pow(((deltaW-6)/20),50)+10);
		  Discriminant_had_W_85_95 -> Fill(discriminant_had_W_85_95);
		  double discriminant_had_t_85_95=10/(pow(((deltaTop-16)/31),50)+10);
		  Discriminant_had_t_85_95 -> Fill(discriminant_had_t_85_95);

		  double discriminant_total_W_85_95=discriminant_lep_W_85_95 * discriminant_had_W_85_95;
		  Discriminant_total_W_85_95 -> Fill(discriminant_total_W_85_95);
		  double discriminant_total_t_85_95=discriminant_lep_t_85_95 * discriminant_had_t_85_95;
		  Discriminant_total_t_85_95 -> Fill(discriminant_total_t_85_95);
		  //
		  //
		  double discriminant_lep_W_95_85=pow(((deltaW_+0.5)/46),50)/(pow(((deltaW_+0.5)/46),50)+10);
		  Discriminant_lep_W_95_85 -> Fill(discriminant_lep_W_95_85);
		  double discriminant_lep_t_95_85=pow(((deltaTop_-5)/49.5),50)/(pow(((deltaTop_-5)/49.5),50)+10);
		  Discriminant_lep_t_95_85 -> Fill(discriminant_lep_t_95_85);
		    
		  double discriminant_had_W_95_85=10/(pow(((deltaW-6)/13.25),50)+10);
		  Discriminant_had_W_95_85 -> Fill(discriminant_had_W_95_85);
		  double discriminant_had_t_95_85=10/(pow(((deltaTop-16)/20.7),50)+10);
		  Discriminant_had_t_95_85 -> Fill(discriminant_had_t_95_85);

		  double discriminant_total_W_95_85=discriminant_lep_W_95_85 * discriminant_had_W_95_85;
		  Discriminant_total_W_95_85 -> Fill(discriminant_total_W_95_85);
		  double discriminant_total_t_95_85=discriminant_lep_t_95_85 * discriminant_had_t_95_85;
		  Discriminant_total_t_95_85 -> Fill(discriminant_total_t_95_85);
		  //

                  0404*/

		  double deltaPhiJetLep_bLep=abs(recoToLHCO[1][1]-recoToLHCO[4][1]);
		  double deltaPhiJetq_bHad=abs(recoToLHCO[0][1]-recoToLHCO[2][1]);
		  double deltaPhiJetMET_bLep=abs(recoToLHCO[1][1]-recoToLHCO[5][1]);
		  double deltaPhiJetq__bHad=abs(recoToLHCO[0][1]-recoToLHCO[3][1]);
		  if (deltaPhiJetLep_bLep>3.141592)
		    {
		      deltaPhiJetLep_bLep=deltaPhiJetLep_bLep-3.1416;
		    }
		  if (deltaPhiJetMET_bLep>3.141592)
		    {
		      deltaPhiJetMET_bLep=deltaPhiJetMET_bLep-3.1416;
		    }
		  if (deltaPhiJetq_bHad>3.141592)
		    {
		      deltaPhiJetq_bHad=deltaPhiJetq_bHad-3.1416;
		    }
		  if (deltaPhiJetq__bHad>3.141592)
		    {
		      deltaPhiJetq__bHad=deltaPhiJetq__bHad-3.1416;
		    }




		  deltaPhiJetLepbLepS->Fill(deltaPhiJetLep_bLep);
		  deltaPhiJetMETbLepS->Fill(deltaPhiJetMET_bLep);
		  deltaPhiJetLepMETbLepS->Fill(deltaPhiJetLep_bLep + deltaPhiJetMET_bLep);	
	      
		  deltaPhiJetQbHadS->Fill(deltaPhiJetq_bHad);
		  deltaPhiJetQ_bHadS->Fill(deltaPhiJetq__bHad);
		  deltaPhiJetQQ_bHadS->Fill(deltaPhiJetq_bHad + deltaPhiJetq__bHad);	


		  double deltaPhiJetMET_bHad=abs(recoToLHCO[0][1]-recoToLHCO[5][1]);
		  double deltaPhiJetLep_bHad=abs(recoToLHCO[0][1]-recoToLHCO[4][1]);
		  if (deltaPhiJetLep_bHad>3.141592)
		    {
		      deltaPhiJetLep_bHad=deltaPhiJetLep_bHad-3.1416;
		    }
		  if (deltaPhiJetMET_bHad>3.141592)
		    {
		      deltaPhiJetMET_bHad=deltaPhiJetMET_bHad-3.1416;
		    }

		  deltaPhiJetLepMETbHadS->Fill(deltaPhiJetLep_bHad + deltaPhiJetMET_bHad);
		  /*0404
		  double bTagDisc=fit_b(discriminants[0])*fit_b(discriminants[1])*fit_cl(discriminants[2])*fit_cl(discriminants[3]);
		  bTagDiscS->Fill(bTagDisc);
		  discriminantbTotalS->Fill( discriminants[0] );
		  discriminantbTotalS->Fill( discriminants[1] );
		  discriminantclTotalS->Fill( discriminants[2] );
		  discriminantclTotalS->Fill( discriminants[3] );
                  0404*/
		  Sel_Inv_Mass_W_Data0 -> Fill(W.M());        
		  Sel_Inv_Mass_W__Data0 -> Fill(W_.M());        
		  Sel_Inv_Mass_t_Data0 -> Fill(top.M());        
		  Sel_Inv_Mass_t__Data0 -> Fill(top_.M());        
		  Sel_Inv_Mass_b_Data0 -> Fill(jet1.M());
		  Sel_Inv_Mass_b__Data0 -> Fill(jet2.M());
                  jets->Fill(7);
                  h_NumberOfJets->Fill(numberOfJets);
                  if (top_.M() > 210)
		    {
		      filet_<<s_id<<" "<<top_.M()<<endl;
                      cout<<"NOJ_"<<numberOfJets<<endl; 
                      h_NumberOfJets__300->Fill(numberOfJets);
		    }
		  if (top.M() > 210)
                    {
		      filet<<s_id<<" "<<top.M()<<endl;
                      cout<<"NOJ"<<numberOfJets<<endl;
		      h_NumberOfJets_300->Fill(numberOfJets);
                    }

		}
              }
	      /*
		if (j==1)
		{
		double elDiscriminante1=exp(-deltaTop)*deltaTop_; 
		ElDiscriminante1_Data1 -> Fill(elDiscriminante1);
		double elDiscriminante2=exp(-deltaW)*deltaW_; 
		ElDiscriminante2_Data1 -> Fill(elDiscriminante2);
		Sel_Inv_Mass_W_Data1 -> Fill(W.M());        
		Sel_Inv_Mass_W__Data1 -> Fill(W_.M());        
		Sel_Inv_Mass_t_Data1 -> Fill(top.M());        
		Sel_Inv_Mass_t__Data1 -> Fill(top_.M());        
		Sel_Inv_Mass_b_Data1 -> Fill(jet1.M());
		Sel_Inv_Mass_b__Data1 -> Fill(jet2.M());
		}
	      */
	    }
	}
      }
    }
  //fin para mover cuando se quiera hacer seleccion de jets en caso de datos reales 05032013    
	  
	  
	  
	  
	  
	  
	  
} 

// ------------ method called once each job just before starting event loop  ------------
double 
FullAndWorking::fit_cl(double x){
  float p=1.53 + 107.14*x - 1354.27*pow(x,2) + 6865.04*pow(x,3) - 18705.7*pow(x,4) + 29732.6*pow(x,5) - 27643.8*pow(x,6) + 13958.1*pow(x,7) - 2960.73*pow(x,8);
  if (p<0)
    {
      p=0.1;
    }
  return p;
}
double 
FullAndWorking::fit_b(double x){
  float p=0.41 + 13.7*x - 154.99*pow(x,2) + 727.69*pow(x,3) - 1865.42*pow(x,4) + 2827.49*pow(x,5) - 2506.64*pow(x,6) + 1165.72*pow(x,7) - 201.49*pow(x,8);
  return p;
}

bool 
FullAndWorking::nextComb(int j, std::vector<int>& comb, std::vector<int>& others, int n, int r)
{
  int l=0;
  int k=0;
  bool lastComb;
  if (comb[0]!=n-r)
    {
      lastComb = false;
      if (comb[j] < n - r + j)
	{
	  std::cout<<"n= "<<n<<std::endl;
	  std::cout<<"r= "<<r<<std::endl;
	  comb[j] = comb[j] + 1;
	  std::cout<<"j1= "<<j<<" comb[j]= "<<comb[j]<<std::endl;
	} 
      else
	{
	  nextComb(j-1, comb, others, n, r);
          comb[j]= comb[j-1] + 1;
	  std::cout<<"j2= "<<j<<" comb[j]= "<<comb[j]<<std::endl;
	}
    }
  else
    {
      lastComb = true;
      std::cout<<"lastComb"<<lastComb<<std::endl;
    }
  for (int i=1; i<n+1; i++)
    {
      if (k<r+1)
	{
	  if (comb[k]!=i)
	    {
	      others[l] = i;
	      std::cout<<"l="<<l<<" others[l]="<<i<<std::endl;
	      l++;
	    }
	  else
	  {
	    k++;
	  }
	}
      else
	{
	  others[l] = i;
	  std::cout<<"l="<<l<<" others[l]="<<i<<std::endl;
	  l++;
	}

    }
  return lastComb; 
}

void 
FullAndWorking::beginJob()
{
  filet.open("top210.txt");
  filet_.open("top_210.txt");
  if (top)
    { 
      file5.open("DiscriminantsD.txt");
      file6.open("recoLHCOD.txt");
      file7.open("Discriminants.txt");
      file8.open("recoLHCO.txt");
      file9.open("genLHCO.txt");
      /*0404
      file.open("Events.txt");
      file1.open("Filter.txt");
      file10.open("Matching.txt");
      prueba.open("Prueba.txt");
      file11.open("MatchingDetail.txt");  
      0404*/
    }
  else
    {
      if (neutralino)
	{
          file5.open("snDiscriminantsD.txt");
          file6.open("snrecoLHCOD.txt");
	  file7.open("snDiscriminants.txt");
	  file8.open("snrecoLHCO.txt");
	  file9.open("sngenLHCO.txt");
	  /*0404
          file.open("snEvents.txt");
	  file1.open("snFilter.txt");
	  file10.open("snMatching.txt");
	  prueba.open("snPrueba.txt");
	  file11.open("snMatchingDetail.txt");
          0404*/ 
	}
      else
	{
          file5.open("scDiscriminantsD.txt");
          file6.open("screcoLHCOD.txt");
	  file7.open("scDiscriminants.txt");
	  file8.open("screcoLHCO.txt");
	  file9.open("scgenLHCO.txt");
	  /*0404
          file.open("scEvents.txt");
	  file1.open("scFilter.txt");
	  file10.open("scMatching.txt");
	  prueba.open("scPrueba.txt");
	  file11.open("scMatchingDetail.txt");
          0404*/
	}
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FullAndWorking::endJob() 
{
  filet.close();
  filet_.close();
  file5.close();
  file6.close();
  file7.close();
  file8.close();
  file9.close();
  ifstream numberOfEvents;
  numberOfEvents.open("TotalNumberofEvents.txt");
  int events;
  numberOfEvents>>events;
  for (int i =0; i<events; i++)
    {
       TotalEvents->Fill(0);
    } 
  /*0404
  file.close();
  file1.close();
  file10.close();
  prueba.close();
  file11.close();
  0404*/
  //system("./prueba.sh");
}

// ------------ method called when starting to processes a run  ------------
void 
FullAndWorking::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
FullAndWorking::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FullAndWorking::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FullAndWorking::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FullAndWorking::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FullAndWorking);	
