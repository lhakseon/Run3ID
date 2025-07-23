#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TH2.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;

void Reader_neu(){
  
  TMVA::Tools::Instance();  
  
  //ofstream myfileM; 
  //myfileM.open("MediumR.txt");
  
    ofstream myfileL; 
    // myfileM.open("Medium_default_flatpt_factory.txt");
    //  myfileL.open("very_Loose_neu_.txt");
    //  myfileL.open("Loose_neu_.txt");
    //    myfileL.open("Medium_neu_.txt");
 //   myfileL.open("Medium.txt");
//    myfileL.open("myloose.txt");
//    myfileL.open("mytight.txt");
 //   myfileL.open("mymedium.txt");
    myfileL.open("mycustom.txt");
  
 double mymarking = 0.9;




         TString methodName = "Cut_Loose_r";
        TString weightfile =  "./dataset/weights/TMVAClassification_TruePV_Cut_Loose_r.weights.xml";

  
  TMVA::Reader *reader = new TMVA::Reader( "!Color" );
   
  float Sieie,ToE,isoC,pt,isoecal,isohcal;
      reader->AddVariable("Sieie",&Sieie);
      reader->AddVariable("ToE",&ToE);
      reader->AddVariable("(isoC>0)? isoC-0.0 : 0.0", &isoC);
      reader->AddVariable("(isoecal-(0.000583374*Ppt) > 0 ) ? isoecal-(0.000583374*Ppt) : 0.0",&isoecal);
      reader->AddVariable("(isohcal-(0.0112465*Ppt+(1.46174e-05)*Ppt*Ppt) > 0 ) ? isohcal-(0.0112465*Ppt+(1.46174e-05)*Ppt*Ppt) : 0.0",&isohcal);
  

  reader->AddSpectator("Ppt",&pt);
  reader->BookMVA(methodName,weightfile); 


  TMVA::MethodCuts* mcuts = dynamic_cast<TMVA::MethodCuts*>(reader->FindCutsMVA(methodName) );
  std::vector<Double_t> cutsMin;
  std::vector<Double_t> cutsMax;
 
//double SEF = mymarking;
   //       double SEF =  0.80;   
  //          double SEF =  0.91;   
//  //    //very loose   
//     double SEF =  0.89;   //loose   
//        double SEF =  0.818;   //custMedium   
//        double SEF =  0.82;   //Medium   
 //     double SEF =  0.73;   //Tight   
//      double SEF =  0.82;   
//
/*std::vector<double> sef_values = {
    0.99, 0.98, 0.97, 0.96, 0.95,
    0.94, 0.93, 0.92, 0.91, 0.90,
    0.89, 0.88, 0.87, 0.86, 0.85,
    0.84, 0.83, 0.82, 0.81, 0.80,
    0.79, 0.78, 0.77, 0.76, 0.75,
    0.74, 0.73, 0.72, 0.71, 0.70
};*/
 std::vector<double> sef_values;
   for (double val = 0.99; val >= 0.70; val -= 0.01)
        sef_values.push_back(val);

for (double SEF : sef_values) { 
  if(mcuts)mcuts->GetCuts(SEF, cutsMin, cutsMax ); 
//  myfileL<<" "<<cutsMax[0]<<" "<<cutsMax[1]<<" "<<cutsMax[2]<<" "<<cutsMax[3]<<" "<<cutsMax[4]<<" "<<endl;  
  myfileL<<SEF<<" "<<cutsMax[0]<<" "<<cutsMax[1]<<" "<<cutsMax[2]<<" "<<cutsMax[3]<<" "<<cutsMax[4]<<" "<<"\n"<<endl;  
}


//myfileL<<" "<<cutsMax[0]<<" "<<cutsMax[1]<<" "<<cutsMax[2]<<" "<<cutsMax[3]<<endl;  
  
 
  
  /*     double SEF =  0.81;   
  if(mcuts)mcuts->GetCuts(SEF, cutsMin, cutsMax ); 
  myfileL<<" "<<cutsMax[0]<<" "<<cutsMax[1]<<" "<<cutsMax[2]<<" "<<cutsMax[3]<<" "<<cutsMax[4]<<endl;  
  */


  
  /*    double SEF =  0.69;   
  if(mcuts)mcuts->GetCuts(SEF, cutsMin, cutsMax ); 
  myfileL<<" "<<cutsMax[0]<<" "<<cutsMax[1]<<" "<<cutsMax[2]<<" "<<cutsMax[3]<<" "<<cutsMax[4]<<" "<<endl;  
  */
  delete reader;





   myfileL.close(); 
  //  myfileM.close(); 
   // myfileT.close(); 
  cout<<"DONE READING CUTS"<<endl;
  
}
