#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/DataLoader.h"
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
using namespace TMVA;


void Reg(){
  
  TMVA::Tools::Instance();
  std::cout << "==> Start TMVARegression" << std::endl;
      ifstream myfile("99per_bar.txt");



  ostringstream xcS,xcH,xcC,xcecal,xchcal;  
    double xS,xH,xC,xecal,xhcal;
 
  if(myfile.is_open()){
    while(!myfile.eof()){
            myfile>>xS>>xH>>xC>>xecal>>xhcal;
    }
  }


  xcS<<xS;
  xcH<<xH;
  xcC<<xC;
  xcecal<<xecal;
  xchcal<<xhcal;
  
    cout<< xS<< " "<< xH<< " " << xC<< " "<< xecal<< " " <<xhcal<< endl; 
  //Output file 
   TString outfileName ("Loose_w_BAR_scaled_TruePV.root");
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  //Declaring the factory
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_TruePV", outputFile,      "!V:!Silent:Color:DrawProgressBar" );
  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");
  dataloader->AddVariable("Sieie", 'F');
  dataloader->AddVariable("ToE", 'F');
  dataloader->AddVariable("(isoC>0)? isoC-0.0 : 0.0", 'F');
  dataloader->AddVariable("(isoecal-(0.000583374*Ppt) > 0 ) ? isoecal-(0.000583374*Ppt) : 0.0", 'F');
  dataloader->AddVariable("(isohcal-(0.0112465*Ppt+(1.46174e-05)*Ppt*Ppt) > 0 ) ? isohcal-(0.0112465*Ppt+(1.46174e-05)*Ppt*Ppt) : 0.0", 'F');
  dataloader->AddSpectator("Ppt", 'F');



  TFile* input = TFile::Open("../../Isopt/Mergedrun3barrel.root");
  
  // --- Register the regression tree
  TTree *signal = (TTree*)input->Get("t_S");
  TTree *background = (TTree*)input->Get("t_B");
  signal->ls();
  dataloader->AddSignalTree(signal, 1.0);
  dataloader->AddBackgroundTree(background, 1.0);
  dataloader->SetSignalWeightExpression("weighT");
  dataloader->SetBackgroundWeightExpression("weighT");



  TCut mycuts =" Ppt > 200 ";
  TCut mycutb =" Ppt > 200 "; 
  dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "");


  TString methodName = "Cut_Loose_r";
  TString methodOptions ="!H:!V:FitMethod=GA:EffMethod=EffSEl:PopSize=1200:Steps=100"; 
  methodOptions +=":VarProp[0]=FMin:VarProp[1]=FMin:VarProp[2]=FMin:VarProp[3]=FMin:VarProp[4]=FMin";
   

   //commenting this in april  
  methodOptions +=":CutRangeMax[0]="+xcS.str(); 
  methodOptions +=":CutRangeMax[1]="+xcH.str();   
  methodOptions +=":CutRangeMax[2]="+xcC.str();
  methodOptions +=":CutRangeMax[3]="+xcecal.str();
  methodOptions +=":CutRangeMax[4]="+xchcal.str();
  factory->BookMethod(dataloader, TMVA::Types::kCuts, methodName, methodOptions);  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();    
   
   // Save the output
   outputFile->Close();



   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;      
   delete factory;
   delete dataloader;
}
