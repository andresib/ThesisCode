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

class FitBTag_TF : public edm::EDAnalyzer {
public:
  explicit FitBTag_TF(const edm::ParameterSet&);
  ~FitBTag_TF();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  //std::ofstream   file11;
  std::string s_bdisc_name;
  
  TH1D * h1_discriminantbTotal;
  TH1D * h1_discriminantclTotal;

  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_20;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_40;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_60;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_80;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_100;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_120;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_140;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_160;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_180;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_200;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_220;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_240;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_260;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_280;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_300;
  TH1D * h1_Difference_InEnergy_bJets_Gen_Reco_320;

  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_20;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_40;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_60;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_80;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_100;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_120;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_140;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_160;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_180;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_200;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_220;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_240;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_260;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_280;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_300;
  TH1D * h1_Difference_InEnergy_clJets_Gen_Reco_320;

  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
};
  

  
// ----------member data ---------------------------



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FitBTag_TF::FitBTag_TF(const edm::ParameterSet& iConfig)
//This code is used to analyze the number of events in which the matching algorithm doesn't work (between reco and gen), as well as, when it works, in this case these results are used to analyze the criteria used by the top group to select the right permutation. This analysis is done dividing the possible answers in three categories, the first one in which the criteria gives the same jets in the same order (matched), the second when the jets are the right ones but the order is not (wrong) and, the third in which the jets are not the right ones (unmatched).
{
  edm::Service<TFileService> fs;
  s_bdisc_name=iConfig.getParameter<std::string>("s_bdisc_name");

  h1_discriminantbTotal = fs->make<TH1D>("h1_Discriminant_b_Total" , "Discriminant_b_Total" , 300 , -1 , 2 );
  h1_discriminantclTotal = fs->make<TH1D>("h1_Discriminant_cl_Total" , "Discriminant_cl_Total" , 300 , -1 , 2 );

  h1_Difference_InEnergy_bJets_Gen_Reco_20 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_20" , "Difference_InEnergy_bJets_Gen_Reco_20" , 4000 , -200 , 200 );   
  h1_Difference_InEnergy_bJets_Gen_Reco_40 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_40" , "Difference_InEnergy_bJets_Gen_Reco_40" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_60 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_60" , "Difference_InEnergy_bJets_Gen_Reco_60" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_80 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_80" , "Difference_InEnergy_bJets_Gen_Reco_80" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_100 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_100" , "Difference_InEnergy_bJets_Gen_Reco_100" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_120 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_120" , "Difference_InEnergy_bJets_Gen_Reco_120" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_140 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_140" , "Difference_InEnergy_bJets_Gen_Reco_140" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_160 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_160" , "Difference_InEnergy_bJets_Gen_Reco_160" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_180 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_180" , "Difference_InEnergy_bJets_Gen_Reco_180" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_200 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_200" , "Difference_InEnergy_bJets_Gen_Reco_200" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_220 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_220" , "Difference_InEnergy_bJets_Gen_Reco_220" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_240 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_240" , "Difference_InEnergy_bJets_Gen_Reco_240" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_260 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_260" , "Difference_InEnergy_bJets_Gen_Reco_260" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_280 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_280" , "Difference_InEnergy_bJets_Gen_Reco_280" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_300 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_300" , "Difference_InEnergy_bJets_Gen_Reco_300" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_bJets_Gen_Reco_320 = fs->make<TH1D>("h1_Difference_InEnergy_bJets_Gen_Reco_320" , "Difference_InEnergy_bJets_Gen_Reco_320" , 4000 , -200 , 200 );

  h1_Difference_InEnergy_clJets_Gen_Reco_20 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_20" , "Difference_InEnergy_clJets_Gen_Reco_20" , 4000 , -200 , 200 );   
  h1_Difference_InEnergy_clJets_Gen_Reco_40 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_40" , "Difference_InEnergy_clJets_Gen_Reco_40" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_60 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_60" , "Difference_InEnergy_clJets_Gen_Reco_60" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_80 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_80" , "Difference_InEnergy_clJets_Gen_Reco_80" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_100 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_100" , "Difference_InEnergy_clJets_Gen_Reco_100" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_120 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_120" , "Difference_InEnergy_clJets_Gen_Reco_120" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_140 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_140" , "Difference_InEnergy_clJets_Gen_Reco_140" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_160 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_160" , "Difference_InEnergy_clJets_Gen_Reco_160" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_180 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_180" , "Difference_InEnergy_clJets_Gen_Reco_180" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_200 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_200" , "Difference_InEnergy_clJets_Gen_Reco_200" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_220 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_220" , "Difference_InEnergy_clJets_Gen_Reco_220" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_240 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_240" , "Difference_InEnergy_clJets_Gen_Reco_240" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_260 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_260" , "Difference_InEnergy_clJets_Gen_Reco_260" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_280 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_280" , "Difference_InEnergy_clJets_Gen_Reco_280" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_300 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_300" , "Difference_InEnergy_clJets_Gen_Reco_300" , 4000 , -200 , 200 );
  h1_Difference_InEnergy_clJets_Gen_Reco_320 = fs->make<TH1D>("h1_Difference_InEnergy_clJets_Gen_Reco_320" , "Difference_InEnergy_clJets_Gen_Reco_320" , 4000 , -200 , 200 );
}



FitBTag_TF::~FitBTag_TF()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
FitBTag_TF::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  ostringstream id;
  id << iEvent.id().event();   // Convert value into a string.
  string s_id = id.str();      // Get the created string from the output stream.

  //selecting only events with MET, 1 lepton and more than 3 jets 
  //bool lepton, MET, Jets;
  //lepton = false;
  //MET =  false;
  //Jets = false;
  int leptonS,METS, JetsS;
  leptonS = 0;
  METS = 0;
  JetsS = 0;
 

  double discriminants[4];

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
  double genEnergy[6];//[5,-5,2,-1,11,-12]
  double recoEnergy[JetsS];//[jets]


  Handle<pat::METCollection> pfMet;
  iEvent.getByLabel("patMETsPF", pfMet);
  //if(pfMet.isValid())
  //{
  //  cout<<"FUNCIONO MET"<<endl;
      
  //MET = true;
  const pat::METCollection* PFMET = pfMet.product();
  METS = PFMET->size();
  const pat::MET& pfMET = (*PFMET)[0];
  //}  
       	  
  Handle<pat::MuonCollection> pfMuon;
  iEvent.getByLabel("semilepMuonsPF", pfMuon);
  //if(pfMuon.isValid())
  //{
  //  cout<<"FUNCIONO MUON"<<endl;
  //lepton = true;
  const pat::MuonCollection* PFMuon = pfMuon.product();
  leptonS = PFMuon->size();
  //}	 

  Handle<pat::ElectronCollection> pfElectron;
  iEvent.getByLabel("semilepElectronsPF", pfElectron);
  //if(pfElectron.isValid())
  //{
  //  cout<<"FUNCIONO ELECTRON"<<endl;
      
  //lepton = true;
  const pat::ElectronCollection* PFElectron = pfElectron.product();
  leptonS = leptonS + PFElectron->size();
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
          recoEnergy[j]=pfJETtmp.energy();
        }    
      compliant = true;
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
      bool GoodDecay = false;
      int numberOfParticlesInTheProcess=gPARTICLES->size();
      
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
                                    }
                                  else
				    {
				      rare=true;
				    } 
				  //file<<"q found"<<endl;
				  genEnergy[2]=decay->energy();
				  
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
				  genEnergy[3]=decay->energy();
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
				  genEnergy[4]=decay->energy();
		
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
				  genEnergy[5]=decay->energy();

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
				  genEnergy[4]=decay->energy();
		
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
				  genEnergy[5]=decay->energy();

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
			  genEnergy[0]=daught->energy();

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
                                  genEnergy[4]=W_decay->energy();

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
				  genEnergy[5]=W_decay->energy();

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
                		  genEnergy[4]=W_decay->energy();

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
				  genEnergy[5]=W_decay->energy();

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
				  genEnergy[2]=W_decay->energy();

				  
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
				  genEnergy[3]=W_decay->energy();

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
			  genEnergy[1]=daught->energy();
 
		
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
		      if (abs((float)genEnergy[i]-(float)(genParton->energy()))==0)
			{
			  //file11<<"jet: "<<i<<" eta: "<<genToLHCO[i][0]<<"    "<<genParton->eta()<<"     "<<pfJET.eta()<<endl;  //file11 is used to study the matching 
			  //file11<<"jet: "<<i<<" phi: "<<genToLHCO[i][1]<<"    "<<genParton->phi()<<"     "<<pfJET.phi()<<endl;
			  //file11<<"jet: "<<i<<" pt: "<<genToLHCO[i][2]<<"    "<<genParton->pt()<<"     "<<pfJET.pt()<<endl;
			  realConfiguration[i]=jetNumber;
                          discriminants[i]=pfJET.bDiscriminator(s_bdisc_name);
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

              h1_discriminantbTotal->Fill( discriminants[0] );
	      h1_discriminantbTotal->Fill( discriminants[1] );
	      h1_discriminantclTotal->Fill( discriminants[2] );
	      h1_discriminantclTotal->Fill( discriminants[3] );


	      for (int i=0;i<2;i++)
		{
		  if (genEnergy[i]<=20)
		    { 		  
		      h1_Difference_InEnergy_bJets_Gen_Reco_20->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>20)&&(genEnergy[i]<40))
		    { 		  
		      h1_Difference_InEnergy_bJets_Gen_Reco_40->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>40)&&(genEnergy[i]<60))
		    { 		  
		      h1_Difference_InEnergy_bJets_Gen_Reco_60->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>60)&&(genEnergy[i]<80))
		    { 		  
		      h1_Difference_InEnergy_bJets_Gen_Reco_80->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>80)&&(genEnergy[i]<100))
		    { 		  
		      h1_Difference_InEnergy_bJets_Gen_Reco_100->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>100)&&(genEnergy[i]<120))
		    { 		  
		      h1_Difference_InEnergy_bJets_Gen_Reco_120->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>120)&&(genEnergy[i]<140))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_140->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>140)&&(genEnergy[i]<160))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_160->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>160)&&(genEnergy[i]<180))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_180->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>180)&&(genEnergy[i]<200))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_200->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>200)&&(genEnergy[i]<220))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_220->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>220)&&(genEnergy[i]<240))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_240->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>240)&&(genEnergy[i]<260))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_260->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>260)&&(genEnergy[i]<280))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_280->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>280)&&(genEnergy[i]<300))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_300->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>300)&&(genEnergy[i]<320))
		    {
		      h1_Difference_InEnergy_bJets_Gen_Reco_320->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  
		  
		} 
	      
	      for (int i=2;i<4;i++)
		{
		  if (genEnergy[i]<=20)
		    { 		  
		      h1_Difference_InEnergy_clJets_Gen_Reco_20->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>20)&&(genEnergy[i]<40))
		    { 		  
		      h1_Difference_InEnergy_clJets_Gen_Reco_40->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>40)&&(genEnergy[i]<60))
		    { 		  
		      h1_Difference_InEnergy_clJets_Gen_Reco_60->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>60)&&(genEnergy[i]<80))
		    { 		  
		      h1_Difference_InEnergy_clJets_Gen_Reco_80->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>80)&&(genEnergy[i]<100))
		    { 		  
		      h1_Difference_InEnergy_clJets_Gen_Reco_100->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>100)&&(genEnergy[i]<120))
		    { 		  
		      h1_Difference_InEnergy_clJets_Gen_Reco_120->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>120)&&(genEnergy[i]<140))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_140->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>140)&&(genEnergy[i]<160))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_160->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>160)&&(genEnergy[i]<180))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_180->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>180)&&(genEnergy[i]<200))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_200->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>200)&&(genEnergy[i]<220))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_220->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>220)&&(genEnergy[i]<240))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_240->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>240)&&(genEnergy[i]<260))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_260->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>260)&&(genEnergy[i]<280))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_280->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>280)&&(genEnergy[i]<300))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_300->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  if ((genEnergy[i]>300)&&(genEnergy[i]<320))
		    {
		      h1_Difference_InEnergy_clJets_Gen_Reco_320->Fill(recoEnergy[realConfiguration[i]]-genEnergy[i]);
		    }
		  
		  
		}
	    }
	}
    }
}
	

void 
FitBTag_TF::beginJob()
{
//  file11.open("PruebaMatching.txt");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FitBTag_TF::endJob() 
{
 // file11.close();
}

// ------------ method called when beginning the processing of a run  ------------                                                                                                                
void
FitBTag_TF::beginRun(edm::Run const&, edm::EventSetup const&)
{
}



// ------------ method called when ending the processing of a run  ------------
void 
FitBTag_TF::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FitBTag_TF::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FitBTag_TF::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FitBTag_TF::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FitBTag_TF);	
