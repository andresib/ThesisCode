#include <memory>


// user include files

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/Common/interface/TriggerNames.h"

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

class PreselectionStudy : public edm::EDAnalyzer {
public:
  explicit PreselectionStudy(const edm::ParameterSet&);
  ~PreselectionStudy();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  std::ofstream   file11;
  std::string bdisc_name;
  std::string outputDiscriminants;
  std::string outputLHCO;
  std::string postfix;
  bool lastCombination, lastCombination2, lastCombination3;
  bool RFijo;
  std::ofstream   file5;
  std::ofstream   file6;


  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  double byHandb(double x);
  double byHandcl(double x);
  bool nextComb (int j, std::vector<int>& comb, std::vector<int>& oth, int n, int r); 
};
  

  
// ----------member data ---------------------------


TH1D * hW;
TH1D * hW_;
TH1D * ht;
TH1D * ht_; 

TH1D * h1_AntibTagDisc;
TH1D * h1_bTagDisc;

  
TH1D *  deltaR_bJetLep;
TH1D *  et_proportion;
TH1D *  transverse_momentum;
TH1D *  mt;


TH1D * Matching;  
TH1D * MatchedSelection;


////


TH1D * hW_p;
TH1D * hW__p;
TH1D * ht_p;
TH1D * ht__p;

TH1D * h1_AntibTagDisc_p;
TH1D * h1_bTagDisc_p;


TH1D *  deltaR_bJetLep_p;
TH1D *  et_proportion_p;
TH1D *  transverse_momentum_p;
TH1D *  mt_p;


TH1D * Matching_p;
TH1D * MatchedSelection_p;



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PreselectionStudy::PreselectionStudy(const edm::ParameterSet& iConfig)
//This code is used to analyze the number of events in which the matching algorithm doesn't work (between reco and gen), as well as, when it works, in this case these results are used to analyze the criteria used by the top group to select the right permutation. This analysis is done dividing the possible answers in three categories, the first one in which the criteria gives the same jets in the same order (matched), the second when the jets are the right ones but the order is not (wrong) and, the third in which the jets are not the right ones (unmatched).
{
  edm::Service<TFileService> fs;
  
  bdisc_name=iConfig.getParameter<std::string>("bdisc_name");
  RFijo=iConfig.getUntrackedParameter<bool>("RFijo",true);
  outputLHCO=iConfig.getParameter<std::string>("outputLHCO");
  outputDiscriminants=iConfig.getParameter<std::string>("outputDiscriminants");
  postfix=iConfig.getParameter<std::string>("postfix");
  
  
  hW = fs->make<TH1D>("hW" , "hW" ,  100 , 0 , 1000 );                
  hW_ = fs->make<TH1D>("hW_" , "hW_" ,  100 , 0 , 1000 );                
  ht = fs->make<TH1D>("ht" , "ht" ,  100 , 0 , 1000 );                
  ht_ = fs->make<TH1D>("ht_" , "ht_" ,  100 , 0 , 1000 );                


  h1_AntibTagDisc = fs->make<TH1D>("h1_AntibTagDisc" , "AntibTagDisc" , 20000 , -1000 , 1000);
  h1_bTagDisc = fs->make<TH1D>("h1_bTagDisc" , "bTagDisc" , 20000 , -1000 , 1000);


  deltaR_bJetLep = fs->make<TH1D>("deltaR_bJetLep" , "deltaR_bJetLep" , 50 ,0 , 100 ); 
  et_proportion = fs->make<TH1D>("et_proportion" , "et_proportion" , 20 ,0 , 2 );
  transverse_momentum = fs->make<TH1D>("transverse_momentum" , "transverse_momentum" , 2000 ,-100 , 100 );
  mt = fs->make<TH1D>("h_mt" , "h_mt" , 20 ,0 , 300 );

  Matching= fs->make<TH1D>("Matching" , "Events, cuts, MC, matching reco" , 5 , 0 , 5);

  MatchedSelection = fs->make<TH1D>("MatchedSelection" , "MatchedSelection" , 10020 , -20 , 10000 );


  ////

  hW_p = fs->make<TH1D>("hW_p" , "hW" ,  100 , 0 , 1000 );
  hW__p = fs->make<TH1D>("hW__p" , "hW_" ,  100 , 0 , 1000 );
  ht_p = fs->make<TH1D>("ht_p" , "ht" ,  100 , 0 , 1000 );
  ht__p = fs->make<TH1D>("ht__p" , "ht_" ,  100 , 0 , 1000 );


  h1_AntibTagDisc_p = fs->make<TH1D>("h1_AntibTagDisc_p" , "AntibTagDisc" , 20000 , -1000 , 1000);
  h1_bTagDisc_p = fs->make<TH1D>("h1_bTagDisc_p" , "bTagDisc" , 20000 , -1000 , 1000);


  deltaR_bJetLep_p = fs->make<TH1D>("deltaR_bJetLep_p" , "deltaR_bJetLep" , 50 ,0 , 100 );
  et_proportion_p = fs->make<TH1D>("et_proportion_p" , "et_proportion" , 20 ,0 , 2 );
  transverse_momentum_p = fs->make<TH1D>("transverse_momentum_p" , "transverse_momentum" , 2000 ,-100 , 100 );
  mt_p = fs->make<TH1D>("h_mt_p" , "h_mt" , 20 ,0 , 300 );

  Matching_p = fs->make<TH1D>("Matching_p" , "Events, cuts, MC, matching reco" , 5 , 0 , 5);

  MatchedSelection_p = fs->make<TH1D>("MatchedSelection_p" , "MatchedSelection" , 10020 , -20 , 10000 );






}



PreselectionStudy::~PreselectionStudy()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PreselectionStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  ostringstream id;
  id << iEvent.id().event();   // Convert value into a string.
  string s_id = id.str();      // Get the created string from the output stream.

  double discriminants[4];
  double Total_Energy_Jets;
  
  Matching->Fill(1); //number of Events


  //selecting only events with MET, 1 lepton and more than 3 jets 
  //  bool lepton, MET, Jets;
  //lepton = false;
  //MET =  false;
  //Jets = false;
  int leptonS,METS, JetsS;
  leptonS = 0;
  METS = 0;
  JetsS = 0;
  int leptonCharge;

  Handle<pat::JetCollection> pfJets;
  iEvent.getByLabel("selectedPatJetsPF", pfJets);
  //if(pfJets.isValid())
  //{
  //  cout<<"FUNCIONO Jets"<<endl;
  //Jets = true;
  const pat::JetCollection* PFJet = pfJets.product();
  JetsS = PFJet->size();
  //}

  //arrays used to store the information required to write the LHCO files
  double genToLHCO[6][5];//[5,-5,2,-1,11,-12][eta,phi,pt,mass,btag]
  double recoToLHCO[JetsS+2][5];//[jets,11,-12][eta,phi,pt,mass,btag]


  Handle<pat::METCollection> pfMet;
  iEvent.getByLabel("patMETsPF", pfMet);
  //if(pfMet.isValid())
  //{
  //  cout<<"FUNCIONO MET"<<endl;
      
  //MET = true;
  const pat::METCollection* PFMET = pfMet.product();
  METS = PFMET->size();
  const pat::MET& pfMET = (*PFMET)[0];
  recoToLHCO[JetsS+1][0]=pfMET.eta();
  recoToLHCO[JetsS+1][1]=pfMET.phi();
  recoToLHCO[JetsS+1][2]=pfMET.pt();
  recoToLHCO[JetsS+1][3]=pfMET.mass();//0;
  recoToLHCO[JetsS+1][4]=0;
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
      recoToLHCO[JetsS][0]=pfMUON.eta();
      recoToLHCO[JetsS][1]=pfMUON.phi();
      recoToLHCO[JetsS][2]=pfMUON.pt();
      recoToLHCO[JetsS][3]=pfMUON.mass();
      recoToLHCO[JetsS][4]=0;
      leptonCharge = pfMUON.charge();
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
      recoToLHCO[JetsS][0]=pfELECTRON.eta();
      recoToLHCO[JetsS][1]=pfELECTRON.phi();
      recoToLHCO[JetsS][2]=pfELECTRON.pt();
      recoToLHCO[JetsS][3]=pfELECTRON.mass();
      recoToLHCO[JetsS][4]=0;
      leptonCharge = pfELECTRON.charge();
      cout<<"paso3"<<endl;
    }
  if (leptonS>1)
    {
      cout<<"WrongJP: more than 1 lepton"<<endl;   
    }
  //}	  

  
  cout<<"MET="<<METS<<" leptons="<<leptonS<<" Jets="<<JetsS<<endl;
  bool compliant = false; 
  if( (METS==1) && (leptonS==1) && (JetsS>3) )  
    {
      for (int j=0; j<JetsS ; j++)
        {
          const pat::Jet& pfJETtmp = (*PFJet)[j];
          recoToLHCO[j][0]=pfJETtmp.eta();
          recoToLHCO[j][1]=pfJETtmp.phi();
          recoToLHCO[j][2]=pfJETtmp.pt();
          recoToLHCO[j][3]=pfJETtmp.mass();
          recoToLHCO[j][4]=0;
        }    
      compliant = true;
      Matching->Fill(2);//number of events with more that 3 jets, MET and just 1 lepton (eletron or muon) 
    }











  // Now we are going to read the Gen level information and store it to be written in LHCO files
  bool t,t_,w,w_,b,b_,q,q_,m_,vm_, rare;
  w=false;
  w_=false;
  b=false;
  b_=false;
  q=false;
  q_=false;
  m_=false;
  vm_=false;
  t=false;
  t_=false;

  int realConfiguration[4]; //array used to stored the number of the jet associated with each one of the partons
  bool matchingReco = false; //flag used to know if there is a matching between reco and gen level for all the jets coming from the tt~ decay. 
  bool tHad; //flag use to know if the hadronic branch is coming from the t decay





  if (compliant)
    {
      ////////////////////////////////////////////////////////////////////MC  
      Handle<reco::GenParticleCollection> gParticles;
      iEvent.getByLabel("genParticles", gParticles);
      const reco::GenParticleCollection* gPARTICLES = gParticles.product();
      //cout<<"Number of Particles: "<<gPARTICLES->size()<<endl;
      unsigned int tDaughters, t_Daughters, WDaughters,W_Daughters;
      tDaughters=0;
      t_Daughters=0;
      WDaughters=0;
      W_Daughters=0;
      rare = false; //(variable used to know if one of the final states is generated more than one time (q q_ b b_ m l), in this case the event is rejected at MC level  
      //bool GoodDecay = false;
      //int numberOfParticlesInTheProcess=gPARTICLES->size();
      
      for(unsigned int i =0; i<gPARTICLES->size();i++)
	{
	  const reco::GenParticle& gcandPARTICLE = (*gPARTICLES)[i];
	  //cout<<"pdgId: "<<gcandPARTICLE.pdgId()<<endl;
	  if (gcandPARTICLE.pdgId() == 6)
	    {
	      
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
                                  if (!q)
				    {
				      q=true;
                                      tHad=true;
				    }
                                  else
				    {
				      rare=true;
				    } 
				  //file<<"q found"<<endl;
				  genToLHCO[2][0]=decay->eta();
				  genToLHCO[2][1]=decay->phi();
				  genToLHCO[2][2]=decay->pt();
				  genToLHCO[2][3]=decay->mass();//0;
				  genToLHCO[2][4]=0;
				  
				}
			      if ( decayId > -5 && decayId < 0 )
				{
                                  if (!q_)
				    {
				      q_=true;
				    }
                                  else
				    {
				      rare=true;
				    }				  
				  //file<<"q_ found"<<endl;
				  genToLHCO[3][0]=decay->eta();
				  genToLHCO[3][1]=decay->phi();
				  genToLHCO[3][2]=decay->pt();
				  genToLHCO[3][3]=decay->mass();//0;
				  genToLHCO[3][4]=0;
				  
				}
			      if ( decayId == -13 )
				{
				  if (!m_)
				    {
				      m_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
                                  //file<<"m found"<<endl;
		
				  genToLHCO[4][0]=decay->eta();
				  genToLHCO[4][1]=decay->phi();
				  genToLHCO[4][2]=decay->pt();
				  genToLHCO[4][3]=decay->mass();//0;
				  genToLHCO[4][4]=0;
		
				}
			      if ( decayId == 14 )
				{
				  if (!vm_)
				    {
				      vm_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
				  //file<<"vm found"<<endl;
				  genToLHCO[5][0]=decay->eta();
				  genToLHCO[5][1]=decay->phi();
				  genToLHCO[5][2]=decay->pt();
				  genToLHCO[5][3]=decay->mass();//0;
				  genToLHCO[5][4]=0;
		
				}
			      if ( decayId == -11 )
				{
				  if (!m_)
				    {
				      m_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
                                  //file<<"m found"<<endl;
		
				  genToLHCO[4][0]=decay->eta();
				  genToLHCO[4][1]=decay->phi();
				  genToLHCO[4][2]=decay->pt();
				  genToLHCO[4][3]=decay->mass();//0;
				  genToLHCO[4][4]=0;
		
				}
			      if ( decayId == 12 )
				{
				  if (!vm_)
				    {
				      vm_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
				  //file<<"vm found"<<endl;
				  genToLHCO[5][0]=decay->eta();
				  genToLHCO[5][1]=decay->phi();
				  genToLHCO[5][2]=decay->pt();
				  genToLHCO[5][3]=decay->mass();//0;
				  genToLHCO[5][4]=0;
		
				}

			    }
			}
		      if ( (daught->pdgId()) == 5 )
			{
			  if (!b)
			    {
			      b=true;
			    }
			  else
			    {
			      rare=true;
			    }
			  //file<<"b found"<<endl;
			  genToLHCO[0][0]=daught->eta();
			  genToLHCO[0][1]=daught->phi();
			  genToLHCO[0][2]=daught->pt();
			  genToLHCO[0][3]=daught->mass();//4.7;
			  genToLHCO[0][4]=2;
			 
			}
		      /* this part was written to analyze the cases in which the final states are obtained by direct decay 
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
				  if (!m_)
				    {
				      m_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
                                  //file<<"m found"<<endl;
		
				  genToLHCO[4][0]=W_decay->eta();
				  genToLHCO[4][1]=W_decay->phi();
				  genToLHCO[4][2]=W_decay->pt();
				  genToLHCO[4][3]=W_decay->mass();//0;
				  genToLHCO[4][4]=0;
		
				}
			      if ( decayId == -14 )
				{
				  if (!vm_)
				    {
				      vm_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
				  //file<<"vm found"<<endl;
				  genToLHCO[5][0]=W_decay->eta();
				  genToLHCO[5][1]=W_decay->phi();
				  genToLHCO[5][2]=W_decay->pt();
				  genToLHCO[5][3]=W_decay->mass();//0;
				  genToLHCO[5][4]=0;
		
				}
			      if ( decayId == 11 )
				{
				  if (!m_)
				    {
				      m_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
                                  //file<<"m found"<<endl;
		
				  genToLHCO[4][0]=W_decay->eta();
				  genToLHCO[4][1]=W_decay->phi();
				  genToLHCO[4][2]=W_decay->pt();
				  genToLHCO[4][3]=W_decay->mass();//0;
				  genToLHCO[4][4]=0;
		
				}
			      if ( decayId == -12 )
				{
				  if (!vm_)
				    {
				      vm_=true;
				    }
                                  else
				    {
				      rare=true;
				    }
				  //file<<"vm found"<<endl;
				  genToLHCO[5][0]=W_decay->eta();
				  genToLHCO[5][1]=W_decay->phi();
				  genToLHCO[5][2]=W_decay->pt();
				  genToLHCO[5][3]=W_decay->mass();//0;
				  genToLHCO[5][4]=0;
		
				}
			      if ( decayId < 5 && decayId > 0 )
				{
                                  if (!q)
				    {
				      q=true;
				    }
                                  else
				    {
				      rare=true;
				    } 
				  //file<<"q found"<<endl;
				  genToLHCO[2][0]=W_decay->eta();
				  genToLHCO[2][1]=W_decay->phi();
				  genToLHCO[2][2]=W_decay->pt();
				  genToLHCO[2][3]=W_decay->mass();//0;
				  genToLHCO[2][4]=0;
				  
				}
			      if ( decayId > -5 && decayId < 0 )
				{
                                  if (!q_)
				    {
				      q_=true;
				    }
                                  else
				    {
				      rare=true;
				    }				  
				  //file<<"q_ found"<<endl;
				  genToLHCO[3][0]=W_decay->eta();
				  genToLHCO[3][1]=W_decay->phi();
				  genToLHCO[3][2]=W_decay->pt();
				  genToLHCO[3][3]=W_decay->mass();//0;
				  genToLHCO[3][4]=0;
				  
				}
				  
			    }
			}
		      if ( (daught->pdgId()) == -5 )
			{
			  if (!b_)
			    {
			      b_=true;
			    }
			  else
			    {
			      rare=true;
			    }
			  //file<<"b_ found"<<endl;
			  genToLHCO[1][0]=daught->eta();
			  genToLHCO[1][1]=daught->phi();
			  genToLHCO[1][2]=daught->pt();
			  genToLHCO[1][3]=daught->mass();//4.7;
			  genToLHCO[1][4]=2; 
		
			}
		      /*this part was written to analyze the cases in which the final states are obtained by direct decay 
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
    

      
      if (b && q && q_ && b_ && m_ && vm_ && !rare)                
	{
	  Matching->Fill(3);  //number of events with semileptonic decay     
	 	  

	  for (int i=0; i<4 ; i++)
	    {
	      realConfiguration[i]=-1;
	    }
	  
          //Now it is going to be done the matching between the gen level (partons) and the reco level (jets) 
	  Total_Energy_Jets = 0;
	  for (int jetNumber=0; jetNumber<JetsS; jetNumber++) 
	    {
	      
	      const pat::Jet& pfJET = (*PFJet)[jetNumber];
              Total_Energy_Jets = Total_Energy_Jets + pfJET.energy();
	      const reco::GenParticle* genParton= pfJET.genParton();
	      if (genParton)
		{
		  cout<<genParton->eta()<<endl;
		  for (int i=0; i<4; i++)
		    {
		      //file10<<"the difference between gen and genMatching is: "<<abs((float)genToLHCO[i][0]-(float)(genParton->eta()))<<endl;
		      if (abs((float)genToLHCO[i][0]-(float)(genParton->eta()))==0)
			{
			  //file11<<"jet: "<<i<<" eta: "<<genToLHCO[i][0]<<"    "<<genParton->eta()<<"     "<<pfJET.eta()<<endl;  //file11 is used to study the matching 
			  //file11<<"jet: "<<i<<" phi: "<<genToLHCO[i][1]<<"    "<<genParton->phi()<<"     "<<pfJET.phi()<<endl;
			  //file11<<"jet: "<<i<<" pt: "<<genToLHCO[i][2]<<"    "<<genParton->pt()<<"     "<<pfJET.pt()<<endl;
			  realConfiguration[i]=jetNumber;
                          discriminants[i]=pfJET.bDiscriminator(bdisc_name);
			}
		    }
	      
		}
	    }
	  
	 	 
	  //////////////////////////////
	  int numberOfMatched=0;
	  for (int i=0; i<4; i++)
	    {
	      if ( realConfiguration[i]!=-1 )
		{
		  numberOfMatched++;
		}
	    }
          if (numberOfMatched == 4)
            {
              matchingReco=true;
            } 

	  
	  if (matchingReco)
	    {

	      Matching->Fill(4); //number of events with the 4 main partons matched with some jet at reco level.     
	      ///////////Writing the recoLHCO. If there is a matching for all the 4 partons, in this file is written all the jets(4 or more), MET and the lepton main information.
	


	      ////In this part we are filling the histograms that later are used to determine the sigma for W, W_, t and t_ 
    

	   
	      math::PtEtaPhiMLorentzVectorD jetb1 (recoToLHCO[realConfiguration[0]][2],recoToLHCO[realConfiguration[0]][0],recoToLHCO[realConfiguration[0]][1],recoToLHCO[realConfiguration[0]][3]); 
	      cout<<"pasoJPPPPa"<<endl;
	      math::PtEtaPhiMLorentzVectorD jetb2 (recoToLHCO[realConfiguration[1]][2],recoToLHCO[realConfiguration[1]][0],recoToLHCO[realConfiguration[1]][1],recoToLHCO[realConfiguration[1]][3]);
	      cout<<"pasoJPPPPb"<<endl;
	      math::PtEtaPhiMLorentzVectorD jetq1 (recoToLHCO[realConfiguration[2]][2],recoToLHCO[realConfiguration[2]][0],recoToLHCO[realConfiguration[2]][1],recoToLHCO[realConfiguration[2]][3]);
	      cout<<"pasoJPPPPc"<<endl;
	      math::PtEtaPhiMLorentzVectorD jetq2 (recoToLHCO[realConfiguration[3]][2],recoToLHCO[realConfiguration[3]][0],recoToLHCO[realConfiguration[3]][1],recoToLHCO[realConfiguration[3]][3]);
	      cout<<"pasoJPPPPd"<<endl;
	      math::PtEtaPhiMLorentzVectorD W = jetq1 + jetq2;
	      cout<<"pasoJPPPP1"<<endl;
	      math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[JetsS][2],recoToLHCO[JetsS][0],recoToLHCO[JetsS][1],recoToLHCO[JetsS][3]);
	      math::PtEtaPhiMLorentzVectorD met (recoToLHCO[JetsS+1][2],recoToLHCO[JetsS+1][0],recoToLHCO[JetsS+1][1],recoToLHCO[JetsS+1][3]);
	      math::PtEtaPhiMLorentzVectorD W_ = ele + met ;
	      cout<<"pasoJPPPP2"<<endl;
              if (W_.M()<30)
		{
		  //file11<<recoToLHCO[JetsS][2]<<" "<<recoToLHCO[JetsS][0]<<" "<<recoToLHCO[JetsS][1]<<" "<<recoToLHCO[JetsS][3]<<endl;
		  //file11<<recoToLHCO[JetsS + 1][2]<<" "<<recoToLHCO[JetsS + 1][0]<<" "<<recoToLHCO[JetsS + 1][1]<<" "<<recoToLHCO[JetsS + 1][3]<<endl;
		}
	      hW->Fill(W.M());
	      hW_->Fill(W_.M());
              math::PtEtaPhiMLorentzVectorD t; 
	      math::PtEtaPhiMLorentzVectorD t_;
              if (tHad)
		{
	          t = W + jetb1;
	          t_ = W_ + jetb2;
		}
              else
		{
	          t = W_+ jetb1;
	          t_ = W + jetb2;
		}



              ht->Fill(t.M());
	      ht_->Fill(t_.M());

	      ////
	    

	      double AntibTagDisc=byHandcl(discriminants[0])*byHandcl(discriminants[1])*byHandb(discriminants[2])*byHandb(discriminants[3]);
              file11<<"ab "<<AntibTagDisc<<endl;
	      h1_AntibTagDisc->Fill(AntibTagDisc/10000);            
	      double bTagDisc=byHandb(discriminants[0])*byHandb(discriminants[1])*byHandcl(discriminants[2])*byHandcl(discriminants[3]);
	      h1_bTagDisc->Fill(bTagDisc/10000000);
	      file11<<"b "<<bTagDisc<<endl;
              double DeltaR_bJetLep = sqrt(pow(jetb2.Phi()-ele.Phi(),2)+pow(jetb2.Eta()-ele.Eta(),2));
              deltaR_bJetLep->Fill(DeltaR_bJetLep);
              double Et_proportion = (jetb1.energy() + jetb2.energy() + W.energy())/Total_Energy_Jets;
	      et_proportion->Fill(Et_proportion);
              double Transverse_momentum = (jetb1+jetb2+W+met+ele).Pt();
	      transverse_momentum->Fill(Transverse_momentum);
              double Mt=sqrt(2*met.Pt()*ele.Pt()*(1-cos(met.Phi()-ele.Phi())));
              mt->Fill(Mt);
    
            
	      ///////
	      int numberOfJets = PFJet->size();  
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
                                 
				      /*
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
				      */
				      math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[JetsS][2],recoToLHCO[JetsS][0],recoToLHCO[JetsS][1],recoToLHCO[JetsS][3]);
				      math::PtEtaPhiMLorentzVectorD met (recoToLHCO[JetsS+1][2],recoToLHCO[JetsS+1][0],recoToLHCO[JetsS+1][1],recoToLHCO[JetsS+1][3]);
				      math::PtEtaPhiMLorentzVectorD WLep = ele + met ;
				      math::PtEtaPhiMLorentzVectorD tLep = WLep + b2;
				      math::PtEtaPhiMLorentzVectorD tHad = WHad + b1;
				      double deltatHad=abs(tHad.M()-172.5);
				      double deltatLep=abs(tLep.M()-172.5);

				      /*
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
					*/
				      double bTagDisc=byHandb(discriminants[0])*byHandb(discriminants[1])*byHandcl(discriminants[2])*byHandcl(discriminants[3]);
				      double AntibTagDisc=byHandcl(discriminants[0])*byHandcl(discriminants[1])*byHandb(discriminants[2])*byHandb(discriminants[3]);
				      double bTagDiscriminant=10000-bTagDisc;
				      double discriminantJP;

				     
				      double Et_proportion = (b1.energy() + b2.energy() + WLep.energy())/Total_Energy_Jets;
				      double Transverse_momentum = (b1+b2+W+met+ele).Pt();
				      double Mt=sqrt(2*met.Pt()*ele.Pt()*(1-cos(met.Phi()-ele.Phi())));
				     
				      ////
				      hW_p->Fill(WHad.M());
                                      hW__p->Fill(WLep.M());
				      ht_p->Fill(tHad.M());
				      ht__p->Fill(tLep.M());
				      h1_AntibTagDisc_p->Fill(AntibTagDisc/10000);
				      h1_bTagDisc_p->Fill(bTagDisc/10000000);
				      et_proportion_p->Fill(Et_proportion);
				      transverse_momentum_p->Fill(Transverse_momentum);
				      mt_p->Fill(Mt);



				      ////







				  
				      /*
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
				      */
                                      discriminantJP = pow((AntibTagDisc-10.05)/24.39,2)+pow((deltaRJetLep_bLep-2.603)/1.286,2)+(pow(deltaWHad/13.79,2)+pow(deltatHad/70.45,2))+pow((Et_proportion-0.8459)/0.1916,2)+pow((Transverse_momentum-39.24)/23.54,2)+pow((Mt-63.18)/30.48,2);

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




	      /////estudio 100 combinaciones
	      // end solo para comparar seleccion con real
              bool Matched=false;
              bool MatchedByType=false;
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
	      

	      ////

















	      ///////
	      int matching=1;  //(1:Unmatched, 2:Wrong,  3:Matched)


	      ///////




              ///////
	      bool ThereIsACombination = false;
	      double mass[3];
	      double minDiscriminantJP=1000;
	      for (int j=0; j<numberOfCombinations;j++)
		{

		  ////////////////
		  int numberOfMatchings=0;
		  int numberOfRightMatchings=0; 
		  for (int k=0; k<4 ; k++)
		    { 
		      if (b1_Combinations[j][1]-1==realConfiguration[k])
			{
			  numberOfMatchings++;
			  if (k==0)
			    {
			      numberOfRightMatchings++;
			    }
			}
		      if (b2_Combinations[j][1]-1==realConfiguration[k])
			{
			  numberOfMatchings++;
			  if (k==1)
			    {
			      numberOfRightMatchings++;
			    }
			}
		      if (W_Combinations[j][1]-1==realConfiguration[k])
			{
			  numberOfMatchings++;
			  if ((k==2)||(k==3))
			    {
			      numberOfRightMatchings++;
			    }
			}
		      if (W_Combinations[j][2]-1==realConfiguration[k])
			{
			  numberOfMatchings++;
			  if ((k==2)||(k==3))
			    {
			      numberOfRightMatchings++;
			    }
			}
		    }
		  if (numberOfMatchings == 4)
		    {
		      matching=2; 
		      if (numberOfRightMatchings == 4)
			{
			  matching=3;
			}
		    }
		  /*
		    b1_Combinations[j][1]-1==realConfiguration[0] //es el jet b1
		    b2_Combinations[j][1]-1==realConfiguration[1] //es el jet b2 
		    W_Combinations[j][1]-1==realConfiguration[2 o 3]  //son q y q~
		    W_Combinations[j][2]-1==realConfiguration[2 o 3]
		  */
		  ///////////////



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
			      discriminants[0]=pfJETtmp.bDiscriminator(bdisc_name);
			      file6<<"0"<<" "<<s_id<<j<<" "<<matching<<endl;
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
	      
		  double discriminantJP = sqrt(pow(WHad.M()-80.385,2.0)+pow(b1.M()-4.7,2.0) +pow(b2.M()-4.7,2.0));
	  
		  if (ThereIsACombination)
		    {
		      file6<<"5"
			   <<" "
			   <<"1"
			   <<" "
			   <<recoToLHCO[JetsS][0]
			   <<" "
			   <<recoToLHCO[JetsS][1]
			   <<" "
			   <<recoToLHCO[JetsS][2]
			   <<" "
			   <<recoToLHCO[JetsS][3]
			   <<" "
			   <<leptonCharge
			   <<" "
			   <<recoToLHCO[JetsS][4]
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
			   <<recoToLHCO[JetsS+1][0]
			   <<" "
			   <<recoToLHCO[JetsS+1][1]
			   <<" "
			   <<recoToLHCO[JetsS+1][2]
			   <<" "
			   <<recoToLHCO[JetsS+1][3]
			   <<" "
			   <<"1"
			   <<" "
			   <<recoToLHCO[JetsS+1][4]
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
			  math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[JetsS][2],recoToLHCO[JetsS][0],recoToLHCO[JetsS][1],recoToLHCO[JetsS][3]);
			  math::PtEtaPhiMLorentzVectorD met (recoToLHCO[JetsS+1][2],recoToLHCO[JetsS+1][0],recoToLHCO[JetsS+1][1],recoToLHCO[JetsS+1][3]);
			  math::PtEtaPhiMLorentzVectorD jet2 (recoToLHCO[1][2],recoToLHCO[1][0],recoToLHCO[1][1],recoToLHCO[1][3]);
			  math::PtEtaPhiMLorentzVectorD W_ = ele + met ;
			  math::PtEtaPhiMLorentzVectorD top_ = W_ + jet2;
	      
			  double deltaW_= (W_.M()-80.385);
			  double deltaTop_ = (top_.M()-172.5);
			  double deltaW= (W.M()-80.385);
			  double deltaTop = (top.M()-172.5);
			  /*
			    if (j==0)
			    {
		  
			    Sel_Inv_Mass_W_Data0 -> Fill(W.M());        
			    Sel_Inv_Mass_W__Data0 -> Fill(W_.M());        
			    Sel_Inv_Mass_t_Data0 -> Fill(top.M());        
			    Sel_Inv_Mass_t__Data0 -> Fill(top_.M());        
			    Sel_Inv_Mass_b_Data0 -> Fill(jet1.M());
			    Sel_Inv_Mass_b__Data0 -> Fill(jet2.M());
          
			    }
			  */
			}
	      
		    }
		}
	      ///////
	    }
	}
    }
}

double 
PreselectionStudy::byHandcl(double x){
  float p=185.268 + 23746.5*x - 294140*pow(x,2) + 1478110*pow(x,3) - 3991330*pow(x,4) + 6283470*pow(x,5) - 5783640*pow(x,6) +2889290*pow(x,7) - 605695*pow(x,8);
  return p;
}
double 
PreselectionStudy::byHandb(double x){
  float p=62.3237 + 2392.79*x - 22118.5*pow(x,2) + 72492.7*pow(x,3) - 77981.6*pow(x,4) - 86901*pow(x,5) + 294793*pow(x,6) - 266675*pow(x,7) + 84802.5*pow(x,8);
  return p;
}	
bool 
PreselectionStudy::nextComb(int j, std::vector<int>& comb, std::vector<int>& others, int n, int r)
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
PreselectionStudy::beginJob()
{
  file5.open(outputDiscriminants+"/Discriminants_MC_"+postfix+".txt");  
  file6.open(outputLHCO+"/LHCO_MC_"+postfix+".txt");
  file11.open("PruebaMatching.txt");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PreselectionStudy::endJob() 
{
  file5.close();  
  file6.close();
  file11.close();
 
   
}

// ------------ method called when beginning the processing of a run  ------------                                                                                                                
void
PreselectionStudy::beginRun(edm::Run const&, edm::EventSetup const&)
{
}



// ------------ method called when ending the processing of a run  ------------
void 
PreselectionStudy::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PreselectionStudy::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PreselectionStudy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PreselectionStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PreselectionStudy);	
