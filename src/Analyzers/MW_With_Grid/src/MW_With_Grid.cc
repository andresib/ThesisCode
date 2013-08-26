// -*- C++ -*-
//
// Package:    MW_With_Grid
// Class:      MW_With_Grid
// 
/**\class MW_With_Grid MW_With_Grid.cc Analyzers/MW_With_Grid/src/MW_With_Grid.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Juan Pablo Gomez Cardona,42 R-021,,
//         Created:  Wed Aug  7 13:37:03 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>



//
// class declaration
//
using namespace std;
class MW_With_Grid : public edm::EDAnalyzer {
public:
  explicit MW_With_Grid(const edm::ParameterSet&);
  ~MW_With_Grid();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  std::string s_id; 
  std::string s_prefix1;
  std::string s_postfix1;
  std::ofstream   os_event_bestPermutation;
  std::ofstream   os_eventsWithoutLHCO;
  std::string s_lhco_source_address;
  bool b_runningLocally;



  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  void run();
  void selectingBestPerm(std::string weights);
  // ----------member data ---------------------------
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
MW_With_Grid::MW_With_Grid(const edm::ParameterSet& iConfig)

{
  s_lhco_source_address=iConfig.getParameter<std::string>("s_lhco_source_address");
  b_runningLocally=iConfig.getParameter<bool>("b_runningLocally");
  //now do what ever initialization is needed

}


MW_With_Grid::~MW_With_Grid()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MW_With_Grid::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ostringstream oss_id;
  oss_id << iEvent.id().event();   // Convert value into a string.                                                                                                                            
  s_id = oss_id.str();      // Get the created string from the output stream.
  run();
}

void 
MW_With_Grid::run()
{
   if (b_runningLocally)
    {
      s_prefix1 = "data/";
      s_postfix1 = "Locally.sh";
    }
  else
    {
      s_prefix1 = "src/data/";
      s_postfix1 = ".sh";
    }
  string s_fileAddress = s_prefix1 + s_lhco_source_address + s_id;
  //Copying input file into Madweight
  const char * cc_command1 = ("cp " + s_fileAddress + " " + s_prefix1 + "Events/input.lhco").c_str();
  system(cc_command1);
  //Running MadWeight
  const char * cc_command2 = ("./"+s_prefix1+"runningMadweight" + s_postfix1).c_str();
  system(cc_command2);
  selectingBestPerm(s_prefix1 + "Events/fermi/fermi_weights.out");

}

void
MW_With_Grid::selectingBestPerm(std::string weights)
{
  ifstream is_Weights;
  is_Weights.open(weights);
  float f_tmp1; 
  float f_weight,f_error,f_weight1, f_error1;
  is_Weights>>f_tmp1;
  cout<<"JP"<<f_tmp1<<endl;
  is_Weights>>f_weight>>f_error;
  f_weight1 = 0.;
  int i_permutation = 0;
  int i_bestPermutation = 0;
  while (f_tmp1/2<1) //Madweight doesn't work if it only has one mass for the top, so I had to put 2 masses but only the first one is the real.  In this step it is selected only the permutations for the first mass   
  {
      i_permutation ++;
      cout<<i_permutation<<endl;
      if (f_weight>f_weight1)
	{
	  f_weight1 = f_weight;
          f_error1 = f_error;
          i_bestPermutation = i_permutation;
	}
      is_Weights>>f_tmp1>>f_weight>>f_error;
    }
  os_event_bestPermutation<<s_id<<" "<<i_bestPermutation<<" "<<f_weight1<<" "<<f_error1<<std::endl;

}







// ------------ method called once each job just before starting event loop  ------------
void 
MW_With_Grid::beginJob()
{
  os_event_bestPermutation.open("Event_BestPermutation_Weight_Error.txt");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MW_With_Grid::endJob() 
{
  os_event_bestPermutation.close();
}

// ------------ method called when starting to processes a run  ------------
void 
MW_With_Grid::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MW_With_Grid::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MW_With_Grid::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MW_With_Grid::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MW_With_Grid::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MW_With_Grid);
