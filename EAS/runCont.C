#include "cstdlib"
#include "TSystem.h"
#include "TableWriter.C"
#include "ContourBuilderpv.C"

void runCont(){
  //-----------------Macro Describtion----------------------------------------------------------
  //
  // This macro wraps the whole procedure in a single shell. 
  // Initialy it caluclates the effective areas -macros (AreaCalc.C + Fitter.C) 
  // and then feed them in a Tree Maker, that makes a Tree with the variables of interest 
  // separated to background/signal and weight based on the pt, for possible MVA techinques 
  // and to derive the cuts. 
  // This file to run, needs the full construction of the CUTID class and also the extra macros
  // for each of the specified process. 
  //
  //--------------------------------------------------------------------------------------------
  double Fin[5][3]= {0};
  double FinEr[5][3]= {0};
  double in,ein,ip,eip,ic,eic;

  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Creating the bins and retrieving the effective area for each"<<endl;
  cout<<"-----------------------------------------------------------"<<endl;








  double emin = 0.0;
  double emax = 0.5;
  ContourBuilderpv(1,emin,emax,in,ein,ip,eip,ic,eic);
  Fin[0][0] = ic;
  Fin[0][1] = in;
  Fin[0][2] = ip;

  FinEr[0][0] = eic;
  FinEr[0][1] = ein;
  FinEr[0][2] = eip;



  emin = 0.5; 
  emax = 1.0;
  ContourBuilderpv(2,emin,emax,in,ein,ip,eip,ic,eic);
  Fin[1][0] = ic;
  Fin[1][1] = in;
  Fin[1][2] = ip;

  FinEr[1][0] = eic;
  FinEr[1][1] = ein;
  FinEr[1][2] = eip;

  emin =1.0; 
  emax =1.4442;
  ContourBuilderpv(3,emin,emax,in,ein,ip,eip,ic,eic);
  Fin[2][0] = ic;
  Fin[2][1] = in;
  Fin[2][2] = ip;
  FinEr[2][0] = eic;
  FinEr[2][1] = ein;
  FinEr[2][2] = eip;
  

 
  emin = 1.566;
  emax = 2.0;
  ContourBuilderpv(4,emin,emax,in,ein,ip,eip,ic,eic);
  Fin[3][0] = ic;
  Fin[3][1] = in;
  Fin[3][2] = ip;
  FinEr[3][0] = eic;
  FinEr[3][1] = ein;
  FinEr[3][2] = eip;


  emin = 2.0;
  emax = 2.5;
  ContourBuilderpv(5,emin,emax,in,ein,ip,eip,ic,eic); 
  Fin[4][0] = ic;
  Fin[4][1] = in;
  Fin[4][2] = ip;
  FinEr[4][0] = eic;
  FinEr[4][1] = ein;
  FinEr[4][2] = eip;


  for(int j  = 0; j < 4; j++)cout<<endl;
  cout<<"Derived the Effective areas --- "<<endl;
  for(int j  = 0; j < 4; j++)cout<<endl;

  TableWriter(Fin,FinEr);


  cout<<"WRITTING THE EFFECTIVE AREAS"<<endl;
  for(int j  = 0; j < 4; j++)cout<<endl;

  cout<<"----------------------------------------------------"<<endl;
  cout<<"Creating the Trees for the Cut Optimization"<<endl;
  for(int j  = 0; j < 4; j++)cout<<endl;

  
  cout<<"----------------------------------------------------"<<endl;



}
