#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <vector>


void LHCO-name-(){
  ifstream input;
  input.open("-input-");
  string line;
  int count=0;
  bool alreadyCharged = false;
  string tmp1,event_perm,tmp2;
  while ((!input.eof())  )//&&(count<5))
    {
      ofstream output;
      string event1;
      if (!alreadyCharged)
	{
	  input>>tmp1>>event_perm>>tmp2;
	}
      //cout<<tmp1<<" "<<event_perm<<" "<<tmp2<<endl;
      if (tmp1.length()!=0)
	{
	  int length=event_perm.length();
	  //cout<<event_perm<<endl;
	  for (int i=0;i<length-1;i++)
	    {
	      event1+=event_perm[i];
	      //cout<<event1<<endl;
	    }
	  cout<<"Event: "<<event1<<endl;
	  ////
	  output.open(("-output-"+event1+".lhco").c_str());
	  output<<tmp1<<" "<<event_perm<<" "<<tmp2<<endl;
	  getline(input,line);
	  for (int i=0;i<6;i++) 
	    {
	      getline(input,line);
	      output<<line<<endl;
	      //cout<<line<<endl;
	    }
	  //cout<<line<<endl;   
	  bool sameEvent=true;
	  int count1=0;
	  while ( (sameEvent) && (!input.eof())    )//&& (count1<3) )
	    {
	      count1++;
	      bool stillSame=false;
	      input>>tmp1;
              if (input.eof()) break;
              input>>event_perm>>tmp2;
	      //cout<<"JP: "<<tmp1<<" "<<event_perm<<" "<<tmp2<<endl;
	      if (tmp1.length()!=0)
		{
		  string event1a;
		  int length=event_perm.length();
		  for (int i=0;i<length;i++)
		    {
		      event1a+=event_perm[i];
		      //cout<<event1a<<endl;
		      //cout<<event1<<endl;
		      if (event1a==event1)
			{
			  stillSame=true;
			  //cout<<"IGUAL"<<endl;
			}
		    }
		  cout<<stillSame<<endl;
		  if (stillSame)
		    {
		      output<<tmp1<<" "<<event_perm<<" "<<tmp2<<endl;
		      getline(input,line);
		      for (int i=0;i<6;i++)
			{
			  getline(input,line);
			  output<<line<<endl;
			  //cout<<line<<endl;
			}
		      //cout<<line<<endl;
		    }
		  else
		    {
		      sameEvent=false;
		      //cout<<"DIFERENTE"<<endl;
		      alreadyCharged=true;
		    }
		  
		}
	    }
	  ////      
	  
	  
	  /*
	    output.open(("/afs/cern.ch/work/j/jgomezca/thesis/mientras/CMSSW_5_3_9/src/LHCO/recoD/"+event1+".lhco").c_str());
	    output<<tmp1<<" "<<event_perm<<" "<<tmp2;
	    for (int i=1;i<169;i++)
	    {
	    getline(input,line);
	    output<<line<<endl;
	    }
	    cout<<line<<endl;
	  */
	  output.close();
	}
      count++;
    }
  exit(0);
}


