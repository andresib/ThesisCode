#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <vector>


void AddingFilePrefix_-name-(){
  ifstream input;
  input.open("-input-/list.txt");
  ofstream output;
  output.open("-input-/listWithPrefix.txt");
  string line;
  getline(input,line);
  output<<"file:"<<line<<endl;
  while ((!input.eof()))
    {
      getline(input,line);
      if (line.length() > 1) 
	{
	  output<<"file:"<<line<<endl;
	}
    }
  output.close();
  exit(0);
}


