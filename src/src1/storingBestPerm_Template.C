#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <vector>


void storingBestPerm(){
  int s_id = -input-;
  ofstream   os_event_bestPermutation;
   ifstream is_Weights;
  is_Weights.open("Events/fermi/fermi_weights.out");
  float f_tmp1;
  float f_weight,f_error,f_weight1, f_error1;
  is_Weights>>f_tmp1;
  cout<<"JP"<<f_tmp1<<endl;
  is_Weights>>f_weight>>f_error;
  f_weight1 = 0.;
  int i_permutation = 0;
  int i_bestPermutation = 0;
  while ((f_tmp1/2<1)&&(!is_Weights.eof())) //In case there is more that one mass.  This is not important now but I prefer leave it.                                                                          
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
  if (f_weight1 !=0)
    {
      os_event_bestPermutation.open("../../resultsLocal/mw_-output-/-input-JP");  //solo para probar 
    }
  else
    {
      os_event_bestPermutation.open("../../resultsLocal/mw_-output-/-input-");
    }
  os_event_bestPermutation<<s_id<<" "<<i_bestPermutation<<" "<<f_weight1<<" "<<f_error1<<std::endl;
  os_event_bestPermutation.close();
  exit(0);
}


