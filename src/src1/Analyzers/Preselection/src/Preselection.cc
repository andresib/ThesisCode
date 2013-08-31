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

class Preselection : public edm::EDAnalyzer {
public:
  explicit Preselection(const edm::ParameterSet&);
  ~Preselection();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  std::string bdisc_name;
  std::string outputDiscriminants;
  std::string outputLHCO;
  std::string postfix;
  bool  RFijo;
  bool lastCombination, lastCombination2, lastCombination3;
  std::ofstream   file5;
  std::ofstream   file6;
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  bool nextComb (int j, std::vector<int>& comb, std::vector<int>& oth, int n, int r); 
  double byHandb(double x);
  double byHandcl(double x);
  
  // ----------member data ---------------------------
  
  TH1D * Sel_Inv_Mass_b_Data0;
  TH1D * Sel_Inv_Mass_b__Data0;
  TH1D * Sel_Inv_Mass_W_Data0;
  TH1D * Sel_Inv_Mass_W__Data0;
  TH1D * Sel_Inv_Mass_t_Data0;
  TH1D * Sel_Inv_Mass_t__Data0;
  

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
Preselection::Preselection(const edm::ParameterSet& iConfig)
  
{
  edm::Service<TFileService> fs;
  RFijo=iConfig.getUntrackedParameter<bool>("RFijo",true);
  bdisc_name=iConfig.getParameter<std::string>("bdisc_name");
  outputLHCO=iConfig.getParameter<std::string>("outputLHCO");
  outputDiscriminants=iConfig.getParameter<std::string>("outputDiscriminants");
  postfix=iConfig.getParameter<std::string>("postfix");
  /* 
     Sel_Inv_Mass_b_Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_b_Data0" , "h_Sel_Inv_Mass_b_Data0" , 50 , 0 , 25 );
     Sel_Inv_Mass_b__Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_b__Data0" , "h_Sel_Inv_Mass_b__Data0" , 50 , 0 , 25 );
     Sel_Inv_Mass_W_Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_W_Data0" , "h_Sel_Inv_Mass_W_Data0" , 50 , 0 , 200 );
     Sel_Inv_Mass_W__Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_W__Data0" , "h_Sel_Inv_Mass_W__Data0" , 50 , 0 , 200 );
     Sel_Inv_Mass_t_Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_t_Data0" , "h_Sel_Inv_Mass_t_Data0" , 50 , 0, 400 );
     Sel_Inv_Mass_t__Data0 = fs->make<TH1D>("h_Sel_Inv_Mass_t__Data0" , "h_Sel_Inv_Mass_t__Data0" , 50 , 0 , 400 );
  */



 

}



Preselection::~Preselection()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Preselection::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  ostringstream id;
  id << iEvent.id().event();   // Convert value into a string.
  string s_id = id.str();      // Get the created string from the output stream.

  double discriminants[4];
  
 

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
  //  //cout<<"FUNCIONO Jets"<<endl;
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
  //  //cout<<"FUNCIONO MET"<<endl;
      
  //MET = true;
  const pat::METCollection* PFMET = pfMet.product();
  METS = PFMET->size();
  const pat::MET& pfMET = (*PFMET)[0];
  recoToLHCO[5][0]=pfMET.eta();
  recoToLHCO[5][1]=pfMET.phi();
  recoToLHCO[5][2]=pfMET.pt();
  recoToLHCO[5][3]=pfMET.mass();//0;
  recoToLHCO[5][4]=0;
  //}  
       	  
  Handle<pat::MuonCollection> pfMuon;
  iEvent.getByLabel("semilepMuonsPF", pfMuon);
  //if(pfMuon.isValid())
  //{
  //  //cout<<"FUNCIONO MUON"<<endl;
  //lepton = true;
  const pat::MuonCollection* PFMuon = pfMuon.product();
  leptonS = PFMuon->size();
  if (leptonS ==1)
    {
      const pat::Muon& pfMUON = (*PFMuon)[0];
      recoToLHCO[4][0]=pfMUON.eta();
      recoToLHCO[4][1]=pfMUON.phi();
      recoToLHCO[4][2]=pfMUON.pt();
      recoToLHCO[4][3]=pfMUON.mass();
      recoToLHCO[4][4]=0;
      leptonCharge = pfMUON.charge();
      //cout<<"paso2"<<endl;
    }
  //}	 

  Handle<pat::ElectronCollection> pfElectron;
  iEvent.getByLabel("semilepElectronsPF", pfElectron);
  //if(pfElectron.isValid())
  //{
  //  //cout<<"FUNCIONO ELECTRON"<<endl;
      
  //lepton = true;
  const pat::ElectronCollection* PFElectron = pfElectron.product();
  leptonS = leptonS + PFElectron->size();
  if ((leptonS == 1) && (PFElectron->size()==1))
    {
      const pat::Electron& pfELECTRON = (*PFElectron)[0];
      recoToLHCO[4][0]=pfELECTRON.eta();
      recoToLHCO[4][1]=pfELECTRON.phi();
      recoToLHCO[4][2]=pfELECTRON.pt();
      recoToLHCO[4][3]=pfELECTRON.mass();
      recoToLHCO[4][4]=0;
      leptonCharge = pfELECTRON.charge();
      //cout<<"paso3"<<endl;
    }
  if (leptonS>1)
    {
      //cout<<"WrongJP: more than 1 lepton"<<endl;   
    }
  //}	  

  
  //cout<<"MET="<<METS<<" leptons="<<leptonS<<" Jets="<<JetsS<<endl;
  if( (METS==1) && (leptonS==1) && (JetsS>3) )  
    {
      
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
	  //cout<<"NW="<<N<<" CW="<<N1<<endl;
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
		  //cout<<"combWHad "<<combinationWHad[j]<<" ";
		}
	      //prueba<<endl;
	      //cout<<endl;
	      for(int j=0; j<N-N1;j++)
		{
		  //prueba<<"othersWHad "<<othersWHad[j]<<" ";
		  //cout<<"othersWHad "<<othersWHad[j]<<" ";
		}
	      //prueba<<endl;
	      //cout<<endl;
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
	      //cout<<"massWHad="<<WHad.M()<<endl;
	      deltaWHad=abs(WHad.M()-80.385);//73.29);
		  
	      //////
	      ////cout<<"paso1"<<endl;
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
		  ////cout<<"paso2"<<endl;
		  int R2 = l;
		  const int N2 = R2+1;
		  //prueba<<"Nb1="<<M<<" Cb1="<<N2<<endl;
		  //cout<<"Nb1="<<M<<" Cb1="<<N2<<endl;
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
		      ////cout<<"paso3"<<endl;
		      //                                                                                                                                                                    
		      for (int k=0; k<N2;k++)
			{
			  //prueba<<"combb1 "<<othersWHad[combinationb1[k]-1]<<" ";
			  //cout<<"combb1 "<<othersWHad[combinationb1[k]-1]<<" ";
			}
		      //prueba<<endl;
		      //cout<<endl;
		      for(int k=0; k<M-N2;k++)
			{
			  //prueba<<"othersb1 "<<othersWHad[othersb1[k]-1]<<" ";
			  //cout<<"othersb1 "<<othersWHad[othersb1[k]-1]<<" ";
			}
		      //prueba<<endl;
		      //cout<<endl;
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
		      //cout<<"massb1="<<b1.M()<<endl;
			  
			  
		      deltab1=abs(b1.M()-4.7);//12.46);
                          
			  
		      //////
		      ////cout<<"paso1a"<<endl;
		      int M1 = M - N2;
		      int tmp = M1 -1;  
		      ////cout<<"M="<<M<<" N2="<<N2<<" M1="<<M1<<" M1-N3="<<tmp<<endl;
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
			  ////cout<<"paso2a"<<endl;
			  int R3 = l1;
			  const int N3 = R3+1;
			  const int tmp1=M1-N3;
			  ////cout<<"M1-N3="<<tmp1<<endl;
			  //prueba<<"Nb2="<<M1<<" Cb2="<<N3<<endl;
			  //cout<<"Nb2="<<M1<<" Cb2="<<N3<<endl;
			  vector<int> combinationb2(N3);
			  //cout<<"paso2b"<<endl;
			  vector<int> othersb2(M1-N3);
			  ////cout<<"paso2c"<<endl;
			  for (int k1=0; k1<N3; k1++)
			    {
			      ////cout<<"paso2d"<<endl;
			      combinationb2[k1]=k1 + 1;
			    }
			  ////cout<<"paso2e"<<endl;
			  for (int k1=0; k1<M1-N3; k1++)
			    {
			      ////cout<<"paso2f"<<endl;
			      othersb2[k1]=N3 + k1 + 1;
			    }
			  ////cout<<"paso2g"<<endl;
			  lastCombination3= false;
			  while (!lastCombination3)
			    {
			      ////cout<<"paso3a"<<endl;
			      //                                                                                                                                                       
			      for (int k1=0; k1<N3;k1++)
				{
				  ////cout<<"paso3b"<<endl;
				  //prueba<<"combb2 "<<othersWHad[othersb1[combinationb2[k1]-1]-1]<<" ";
				  //cout<<"combb2 "<<othersWHad[othersb1[combinationb2[k1]-1]-1]<<" ";
				}
			      ////cout<<"paso3c"<<endl;
			      //prueba<<endl;
			      //cout<<endl;
			      for(int k1=0; k1<M1-N3;k1++)
				{
				  ////cout<<"paso3d"<<endl;
				  //prueba<<"othersb2 "<<othersWHad[othersb1[othersb2[k1]-1]-1]<<" ";
				  //cout<<"othersb2 "<<othersWHad[othersb1[othersb2[k1]-1]-1]<<" ";
				}
			      ////cout<<"paso3e"<<endl;
			      //prueba<<endl;
			      //cout<<endl;
			      // 
			      math::PtEtaPhiMLorentzVectorD b2 (0,0,0,0);
			      ////cout<<"paso3f"<<endl;
			      for (int k1 = 0; k1 < N3; k1++)
				{
				  ////cout<<"paso3g"<<endl;
				  const pat::Jet& pfJETtmp = (*PFJet)[othersWHad[othersb1[combinationb2[k1]-1]-1]-1];
				  if (k1<1)
				    {
				      discriminants[k1+1]=pfJETtmp.bDiscriminator(bdisc_name);
				    }

				  math::PtEtaPhiMLorentzVectorD jetTmp (pfJETtmp.pt(),pfJETtmp.eta(),pfJETtmp.phi(),pfJETtmp.mass());  
				  b2 = jetTmp + b2;
				}
			      double deltaRJetLep_bLep=pow(b2.eta()-recoToLHCO[4][0],2)+pow(b2.phi()-recoToLHCO[4][1],2);


			      math::PtEtaPhiMLorentzVectorD ele (recoToLHCO[4][2],recoToLHCO[4][0],recoToLHCO[4][1],recoToLHCO[4][3]);
			      math::PtEtaPhiMLorentzVectorD met (recoToLHCO[5][2],recoToLHCO[5][0],recoToLHCO[5][1],recoToLHCO[5][3]);
			      math::PtEtaPhiMLorentzVectorD WLep = ele + met ;
			      math::PtEtaPhiMLorentzVectorD tLep = WLep + b2;
			      math::PtEtaPhiMLorentzVectorD tHad = WHad + b1;
			      double deltatHad=abs(tHad.M()-172.5);
			      double deltatLep=abs(tLep.M()-172.5);

			      //prueba<<"massbLep="<<b2.M()<<endl;
			      //cout<<"massbLep="<<b2.M()<<endl;
			      ////cout<<"paso3h"<<endl;
			      

			      double bTagDisc=byHandb(discriminants[0])*byHandb(discriminants[1])*byHandcl(discriminants[2])*byHandcl(discriminants[3]);
			      double AntibTagDisc=byHandcl(discriminants[0])*byHandcl(discriminants[1])*byHandb(discriminants[2])*byHandb(discriminants[3]);
			      double bTagDiscriminant=10000-bTagDisc;
			      double discriminantJP;
			      
			      discriminantJP = pow(AntibTagDisc,2)*pow(deltaRJetLep_bLep,2)*(pow(deltaWHad,2)*pow(deltatHad,2))/100000000000000000000.0;


			      endOfComparison = false;
			      int j=0;
			      while (!endOfComparison)
				{
				  ////cout<<"j="<<j<<endl;
				  //prueba<<"discriminantJP= "<<discriminantJP<<" WComb[j][0]= "<<W_Combinations[j][0]<<endl;
				  //cout<<"discriminantJP= "<<discriminantJP<<" WComb[j][0]= "<<W_Combinations[j][0]<<" j="<<j<<endl;
				  if (discriminantJP<W_Combinations[j][0])
				    {
				      //prueba<<"Clasifica Combinacion"<<endl;
				      //cout<<"Clasifica Combinacion"<<endl;
				      endOfComparison=true;
				      for (int k=numberOfCombinations-2;k>j-1;k--)
					{
					  ////cout<<"k="<<k<<endl;
					  for (int l=0; l<numberOfJets;l++)
					    {
					      ////cout<<"l="<<l<<endl;
					      W_Combinations[k+1][l]=W_Combinations[k][l];
					      b1_Combinations[k+1][l]=b1_Combinations[k][l];
					      b2_Combinations[k+1][l]=b2_Combinations[k][l];
					    }
					}
				      for (int l=1; l<N1+1;l++)
					{
					  ////cout<<"lW="<<l<<endl;
					  W_Combinations[j][l]=combinationWHad[l-1];
					  //prueba<<" JPW "<<combinationWHad[l-1];
					}
				      //prueba<<endl;
				  
				      for (int l=1; l<N2+1;l++)
					{
					  ////cout<<"lb1="<<l<<endl;
					  b1_Combinations[j][l]=othersWHad[combinationb1[l-1]-1];//combinationb1[l-1];
					  //prueba<<" JPb1 "<<combinationb1[l-1];
					}
				      //prueba<<endl;
					  
				      for (int l=1; l<N3+1;l++)
					{
					  ////cout<<"lb2="<<l<<endl;
					  b2_Combinations[j][l]=othersWHad[othersb1[combinationb2[l-1]-1]-1];//combinationb2[l-1];
					  //prueba<<" JPb2 "<<combinationb2[l-1];
					}
				      //prueba<<endl;
				  
				      W_Combinations[j][0]=discriminantJP;
				      for(int l=N1+1; l<numberOfJets;l++)
					{
					  ////cout<<"lWo="<<l<<endl;
					  W_Combinations[j][l]=1000;
					}
				  
				      for(int l=N2+1; l<numberOfJets;l++)
					{
					  ////cout<<"lb1o="<<l<<endl;
					  b1_Combinations[j][l]=1000;
					}
				  
				      for(int l=N3+1; l<numberOfJets;l++)
					{
					  ////cout<<"lb2o="<<l<<endl;
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
			    
			      //cout<<"paso5a"<<endl;
			      lastCombination3 = nextComb(R3, combinationb2, othersb2, M1, R3);
			      //cout<<"paso6a"<<lastCombination3<<endl;
			    }
			}
			  
		      //////
			  
			  
		      ////cout<<"paso5"<<endl;
		      //cout<<"lastCombBefore="<<lastCombination2<<" R2="<<R2<<endl;
		      lastCombination2 = nextComb(R2, combinationb1, othersb1, M, R2);
		      //cout<<"lastCombAfter="<<lastCombination2<<" R2="<<R2<<endl;
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
    }
}
//fin para mover cuando se quiera hacer seleccion de jets en caso de datos reales 05032013    
  
	  
	  
	  
	  
	  
	  
 

// ------------ method called once each job just before starting event loop  ------------
double
Preselection::byHandcl(double x){
  float p=185.268 + 23746.5*x - 294140*pow(x,2) + 1478110*pow(x,3) - 3991330*pow(x,4) + 6283470*pow(x,5) - 5783640*pow(x,6) +2889290*pow(x,7) - 605695*pow(x,8);
  return p;
}
double 
Preselection::byHandb(double x){
  float p=62.3237 + 2392.79*x - 22118.5*pow(x,2) + 72492.7*pow(x,3) - 77981.6*pow(x,4) - 86901*pow(x,5) + 294793*pow(x,6) - 266675*pow(x,7) + 84802.5*pow(x,8);
  return p;
}

bool 
Preselection::nextComb(int j, std::vector<int>& comb, std::vector<int>& others, int n, int r)
{
  int l=0;
  int k=0;
  bool lastComb;
  if (comb[0]!=n-r)
    {
      lastComb = false;
      if (comb[j] < n - r + j)
	{
	  //std::cout<<"n= "<<n<<std::endl;
	  //std::cout<<"r= "<<r<<std::endl;
	  comb[j] = comb[j] + 1;
	  //std::cout<<"j1= "<<j<<" comb[j]= "<<comb[j]<<std::endl;
	} 
      else
	{
	  nextComb(j-1, comb, others, n, r);
          comb[j]= comb[j-1] + 1;
	  //std::cout<<"j2= "<<j<<" comb[j]= "<<comb[j]<<std::endl;
	}
    }
  else
    {
      lastComb = true;
      //std::cout<<"lastComb"<<lastComb<<std::endl;
    }
  for (int i=1; i<n+1; i++)
    {
      if (k<r+1)
	{
	  if (comb[k]!=i)
	    {
	      others[l] = i;
	      //std::cout<<"l="<<l<<" others[l]="<<i<<std::endl;
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
	  //std::cout<<"l="<<l<<" others[l]="<<i<<std::endl;
	  l++;
	}

    }
  return lastComb; 
}

void 
Preselection::beginJob()
{
  file5.open(outputDiscriminants+"/Discriminants_data_"+postfix+".txt");  
  file6.open(outputLHCO+"/LHCO_data_"+postfix+".txt");
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Preselection::endJob() 
{
  file5.close();  
  file6.close();
}

// ------------ method called when starting to processes a run  ------------
void 
Preselection::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Preselection::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Preselection::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Preselection::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Preselection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Preselection);	
