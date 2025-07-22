#include <stdlib.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>

#include "CutID.C"
//#include "tdrstyle.C"
//#include "Variables.C"
using namespace std;
  void run_ParseEAs(){
    CutID t; 
  double Fin[5][3]= {0};

cout<<"start"<<endl;
     ifstream myfile; 
             myfile.open("../EAS/EA_1_60_TruePV_test.txt");
//             myfile.open("../EAS/Effective_areas_1_60_TruePv_test.txt");

  cout<<"------------------------------------------------------------"<<endl;
  cout<<"                   Reading the EAS "<<endl;
  cout<<"-----------------------------------------------------------"<<endl;

  
  int i = 0; 
  if(myfile.is_open()){   // blocked this line for try
 while(!myfile.eof()) {         //blocked this line as well
    //charged ecal hcal
      double iC,iN,iP;      
      myfile>>iC>>iN>>iP;
      Fin[i][0] = iC; 
      Fin[i][1] = iN; 
      Fin[i][2] = iP; 




      cout<<i<<" "<<iC<<" "<<iN<<" "<<iP<<endl;
            i++;
      if( i > 4) goto conti; 

       }
  }

 conti:
    
      cout<<"----------------------Barrel--------------"<<endl;
  t.CutBasedID(1,Fin);
  cout<<"----------------------End Cap--------------"<<endl;
  t.CutBasedID(2,Fin);
  

  
 cout<<"..................End of Programme............"<<endl;



  myfile.close();
//  return 0;
 cout<<"done"<<endl; 

}

