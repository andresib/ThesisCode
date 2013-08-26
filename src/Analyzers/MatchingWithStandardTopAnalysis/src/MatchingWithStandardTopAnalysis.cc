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

class MatchingWithStandardTopAnalysis : public edm::EDAnalyzer {
public:
  explicit MatchingWithStandardTopAnalysis(const edm::ParameterSet&);
  ~MatchingWithStandardTopAnalysis();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  std::ofstream   file11;
  std::string bdisc_name;
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
};
  

  
// ----------member data ---------------------------


TH1D * Matching;  
TH1D * MatchingStandardTopAnalysis;  

TH1D * W_Wrong;
TH1D * W__Wrong;        
TH1D * t_Wrong;        
TH1D * t__Wrong;        

TH1D * W_Matched;        
TH1D * W__Matched;        
TH1D * t_Matched;        
TH1D * t__Matched;        

TH1D * W_Unmatched;        
TH1D * W__Unmatched;        
TH1D * t_Unmatched;        
TH1D * t__Unmatched;        



TH1D * hW;
TH1D * hW_;
TH1D * ht;
TH1D * ht_; 

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MatchingWithStandardTopAnalysis::MatchingWithStandardTopAnalysis(const edm::ParameterSet& iConfig)
//This code is used to analyze the number of events in which the matching algorithm doesn't work (between reco and gen), as well as, when it works, in this case these results are used to analyze the criteria used by the top group to select the right permutation. This analysis is done dividing the possible answers in three categories, the first one in which the criteria gives the same jets in the same order (matched), the second when the jets are the right ones but the order is not (wrong) and, the third in which the jets are not the right ones (unmatched).
{
  edm::Service<TFileService> fs;
  bdisc_name=iConfig.getParameter<std::string>("bdisc_name");

  Matching= fs->make<TH1D>("Matching" , "Events, cuts, MC, matching reco" , 5 , 0 , 5);
  MatchingStandardTopAnalysis= fs->make<TH1D>("MatchingStandardTopAnalysis" , "Unmatched, matched and wrong" , 5 , 0 , 5);  
  //Invariant mass of W, W_, t and t_ for the three cases: Wrong, Matched and Unmatched. 
  W_Wrong = fs->make<TH1D>("W_Wrong" , "W_Wrong" ,  100 , 0 , 1000 );        
  W__Wrong = fs->make<TH1D>("W__Wrong" , "W__Wrong" ,  100 , 0 , 1000 );        
  t_Wrong = fs->make<TH1D>("t_Wrong" , "t_Wrong" ,  100 , 0 , 1000 );        
  t__Wrong = fs->make<TH1D>("t__Wrong" , "t__Wrong" ,  100 , 0 , 1000 );        

  W_Matched = fs->make<TH1D>("W_Matched" , "W_Matched" ,  100 , 0 , 1000 );                
  W__Matched = fs->make<TH1D>("W__Matched" , "W__Matched" ,  100 , 0 , 1000 );                
  t_Matched = fs->make<TH1D>("t_Matched" , "t_Matched" ,  100 , 0 , 1000 );                
  t__Matched = fs->make<TH1D>("t__Matched" , "t__Matched" ,  100 , 0 , 1000 );                

  W_Unmatched = fs->make<TH1D>("W_Unmatched" , "W_Unmatched" ,  100 , 0 , 1000 );                
  W__Unmatched = fs->make<TH1D>("W__Unmatched" , "W__Unmatched" ,  100 , 0 , 1000 );                
  t_Unmatched = fs->make<TH1D>("t_Unmatched" , "t_Unmatched" ,  100 , 0 , 1000 );                
  t__Unmatched = fs->make<TH1D>("t__Unmatched" , "t__Unmatched" ,  100 , 0 , 1000 );                

  hW = fs->make<TH1D>("hW" , "hW" ,  100 , 0 , 1000 );                
  hW_ = fs->make<TH1D>("hW_" , "hW_" ,  100 , 0 , 1000 );                
  ht = fs->make<TH1D>("ht" , "ht" ,  100 , 0 , 1000 );                
  ht_ = fs->make<TH1D>("ht_" , "ht_" ,  100 , 0 , 1000 );                


}



MatchingWithStandardTopAnalysis::~MatchingWithStandardTopAnalysis()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchingWithStandardTopAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  ostringstream id;
  id << iEvent.id().event();   // Convert value into a string.
  string s_id = id.str();      // Get the created string from the output stream.


  
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
	  for (int jetNumber=0; jetNumber<JetsS; jetNumber++) 
	    {
	      
	      const pat::Jet& pfJET = (*PFJet)[jetNumber];
	      const reco::GenParticle* genParton= pfJET.genParton();
	      if (genParton)
		{
		  cout<<genParton->eta()<<endl;
		  for (int i=0; i<4; i++)
		    {
		      //file10<<"the difference between gen and genMatching is: "<<abs((float)genToLHCO[i][0]-(float)(genParton->eta()))<<endl;
		      if (abs((float)genToLHCO[i][0]-(float)(genParton->eta()))==0)
			{
			  file11<<"jet: "<<i<<" eta: "<<genToLHCO[i][0]<<"    "<<genParton->eta()<<"     "<<pfJET.eta()<<endl;  //file11 is used to study the matching 
			  file11<<"jet: "<<i<<" phi: "<<genToLHCO[i][1]<<"    "<<genParton->phi()<<"     "<<pfJET.phi()<<endl;
			  file11<<"jet: "<<i<<" pt: "<<genToLHCO[i][2]<<"    "<<genParton->pt()<<"     "<<pfJET.pt()<<endl;
			  realConfiguration[i]=jetNumber;
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
              /*
	      if (W_.M()<30)
		{
		  file11<<recoToLHCO[JetsS][2]<<" "<<recoToLHCO[JetsS][0]<<" "<<recoToLHCO[JetsS][1]<<" "<<recoToLHCO[JetsS][3]<<endl;
		  file11<<recoToLHCO[JetsS + 1][2]<<" "<<recoToLHCO[JetsS + 1][0]<<" "<<recoToLHCO[JetsS + 1][1]<<" "<<recoToLHCO[JetsS + 1][3]<<endl;
		}
	      */
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
	    }

	

    


	  //now we are going to calculate the invariant mass of t and t~ following the guides given by the top group.  Here is only had in mind the permutations and not the resolution, it is meant that there are not kinematical fits applied in order to improve the resolution.
	  bool InvMassTopWithStandardAnalysis = true;
	  if (InvMassTopWithStandardAnalysis && matchingReco)//it is asked matchingReco to be true in order to assure that the algorithm of matching worked properly and thus, to be able to make an analysis of the matching got when it is used this criteria to select the permutation 
	    {
	      //now we are going to select the 4 jets with the highest pt.  Note:  this in fact is not important because in PAT they are organized by pt, but I  realized it when I had written the code
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
		  double maxvalue = -1000000;//initial value that is lower than any one else
		  int jetSelected=-1;
		  for (int k = 0; k != JetsS; ++k)
		    {
		      const pat::Jet& aJet = (*PFJet)[k];
		      cout<<"Jet "<<k<<" pt="<<aJet.pt()<<endl;
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
		  cout<<"Selected Jet "<<i<<" pt="<<maxvalue<<endl;
		  //prueba<<maxvalue<<endl;
		  selectedJets[i]=jetSelected;
		  cout<<"PT sel jet "<<i<<" jet number: "<<jetSelected;
		}
	      
	      
	      ////////////////////////////////////////////////////////              

	      ////////////////////////////////////////////////////////
	      //now we are going to select from the last 4 jets,  the 2 with the highest b-discriminant. 
	      
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
		  double maxvalue = -1000000;//initial value lower than any one else
		  int jetSelected=-1;
		  for (int k = 0; k != 4; ++k)
		    {
		      const pat::Jet& aJet = (*PFJet)[selectedJets[k]];
		      disc = aJet.bDiscriminator(bdisc_name);
		      cout<<"Jet "<<k<<" disc="<<disc<<endl;
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
		      if (disc == -1)
			{
			  maxvalue=disc;
			  jetSelected=selectedJets[k];
			}

		    }
		  selectedJets_disc[i]=maxvalue;
		  cout<<"Selected Jet "<<i<<" disc="<<maxvalue<<endl;
		  //prueba<<minvalue<<endl;
		  selectedJetsd[i]=jetSelected;
		  cout<<"bTag sel "<<i<<" jet number "<<jetSelected;
		}
	      
	      
	      int numberOfMatchings = 0;
	      bool MatchingTop = false;//flaf use to decide if the jets selected with the standar top analysis are the matched jets to partons 
	      for (int i=0; i<4; i++)
		{
		  for (int j=0; j<4; j++)
		    {
		      if (selectedJetsd[i] == realConfiguration[j])
			{
			  numberOfMatchings++;
			}
		    }
		}
	  
	      if (numberOfMatchings == 4)
		{
		  MatchingTop = true; 
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



	      math::PtEtaPhiMLorentzVectorD top_1 = W_ + jetb1;
	      math::PtEtaPhiMLorentzVectorD top1 = W + jetb1;
	      math::PtEtaPhiMLorentzVectorD top_2 = W_ + jetb2;
	      math::PtEtaPhiMLorentzVectorD top2 = W + jetb2;

              
	      double deltaW_= (W_.M()-80.385);
	      double deltaW= (W.M()-80.385);
	      double deltaTop1 = (top1.M()-172.5);
	      double deltaTop_1 = (top_1.M()-172.5);
	      double deltaTop2 = (top2.M()-172.5);
	      double deltaTop_2 = (top_2.M()-172.5);

          
	      //these are the sigma values obtained with ht, ht_, hw and hw_  
	      double chiW = deltaW/13.55;
	      double chiW_ = deltaW_/38.35;
	      double chiTop1 = deltaTop1/70.04;
	      double chiTop_1 = deltaTop_1/72.3;
	      double chiTop2 = deltaTop2/70.04;
	      double chiTop_2 = deltaTop_2/72.3;

	      double chiASquare = pow(chiW,2) + pow(chiW_,2) + pow(chiTop1,2) + pow(chiTop_2,2);  
	      double chiBSquare = pow(chiW,2) + pow(chiW_,2) + pow(chiTop2,2) + pow(chiTop_1,2);
	 

	      math::PtEtaPhiMLorentzVectorD top_; 
	      math::PtEtaPhiMLorentzVectorD top;
	      math::PtEtaPhiMLorentzVectorD b_; 
	      math::PtEtaPhiMLorentzVectorD b;
            

	      if (chiASquare < chiBSquare)
		{
		  top_ = top_2;
		  top = top1;
		  b_ = jetb2;
		  b = jetb1;
		}
	      else
		{
		  top_ = top_1;
		  top = top2;
		  b_ = jetb1;
		  b = jetb2;
		}
	      bool MatchingTopRight = false; //flag use to decide if the matching done with the standard top analysis is right in order 
	      const pat::Jet& realJETb1 = (*PFJet)[realConfiguration[0]];
	      math::PtEtaPhiMLorentzVectorD realjetb1 (realJETb1.pt(),realJETb1.eta(),realJETb1.phi(),realJETb1.mass());
	      const pat::Jet& realJETb2 = (*PFJet)[realConfiguration[1]];
	      math::PtEtaPhiMLorentzVectorD realjetb2 (realJETb2.pt(),realJETb2.eta(),realJETb2.phi(),realJETb2.mass());
	      if (MatchingTop)
		{
		  if ((b==realjetb1)&&(b_==realjetb2))
		    {
		      MatchingTopRight = true;
		    }
		}
	      

	      math::PtEtaPhiMLorentzVectorD t_b = top - b;
	      math::PtEtaPhiMLorentzVectorD t__b = top_ - b_;
	  
	      if (!MatchingTop)
		{
		  MatchingStandardTopAnalysis->Fill(1);
		  cout<<W.M()<<" "<<W_.M()<<" "<<top.M()<<" "<< top_.M()<<endl;              
		  W_Unmatched->Fill(W.M());        
		  W__Unmatched->Fill(W_.M());        
		  t_Unmatched->Fill(top.M());        
		  t__Unmatched->Fill(top_.M());        
	      
		}
	      else
		{
		  if (MatchingTopRight)
		    {
		      MatchingStandardTopAnalysis->Fill(2);
		      W_Matched -> Fill(W.M());        
		      W__Matched -> Fill(W_.M());        
		      t_Matched -> Fill(top.M());        
		      t__Matched -> Fill(top_.M());
		    }
		  else
		    {
		      MatchingStandardTopAnalysis->Fill(3);
		      W_Wrong -> Fill(W.M());        
		      W__Wrong -> Fill(W_.M());        
		      t_Wrong -> Fill(top.M());        
		      t__Wrong -> Fill(top_.M());
		    }   
		}
	      
	    }
 
	}
    }
}
	

void 
MatchingWithStandardTopAnalysis::beginJob()
{

  //  file11.open("PruebaMatching.txt");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MatchingWithStandardTopAnalysis::endJob() 
{
 
  // file11.close();
 
   
}

// ------------ method called when beginning the processing of a run  ------------                                                                                                                
void
MatchingWithStandardTopAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}



// ------------ method called when ending the processing of a run  ------------
void 
MatchingWithStandardTopAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MatchingWithStandardTopAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MatchingWithStandardTopAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MatchingWithStandardTopAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchingWithStandardTopAnalysis);	
