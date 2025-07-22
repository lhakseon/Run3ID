#include <iostream>
#include <fstream>

void TableWriter(double EAB[5][3],double EABerr[5][3]){

  ofstream myfile; 
  
  myfile.open("Effective_areas_1_60_TruePv_test.txt");
  

  myfile<<"| *bin* | *EA worstchargedhadrons* | *ECAL Isolation* | *HCAL Isolation* |"<<endl;
  myfile<<"| 0.0<abs(&eta;)<0.5|"<< EAB[0][0] <<" "<< EABerr[0][0] <<"|"<< EAB[0][1] <<" "<< EABerr[0][1] <<"|"<< EAB[0][2] <<" "<< EABerr[0][2] <<"|"<<endl;
  myfile<<"| 0.5<abs(&eta;)<1.0|"<< EAB[1][0] <<" "<< EABerr[1][0] <<"|"<< EAB[1][1] <<" "<< EABerr[1][1] <<"|"<< EAB[1][2] <<" "<< EABerr[1][2] <<"|"<<endl;
  myfile<<"| 1.0<abs(&eta;)<1.4442 |"<<  EAB[2][0] <<" "<<  EABerr[2][0] <<" | "<< EAB[2][1] <<" "<< EABerr[2][1] <<"|"<< EAB[2][2] <<" "<< EABerr[2][2] <<"|"<<endl;
  myfile<<"| 1.566<abs(&eta;)<2.0|"<< EAB[3][0]   <<" "<< EABerr[3][0] <<" | "<< EAB[3][1] <<" "<< EABerr[3][1]<<" | "<< EAB[3][2] <<" "<< EABerr[3][2]<<" | "<<endl;
  myfile<<"| 2.0 <abs(&eta;) < 2.5 |"<<EAB[4][0]  <<" "<<EABerr[4][0]<<" | "<< EAB[4][1] <<" "<< EABerr[4][1]<<" | "<< EAB[4][2] <<" "<< EABerr[4][2] <<" | "<<endl;
  //myfile<<"| 2.0 <abs(&eta;) < 2.5 |"<<EAB[5][0]  <<" "<<EABerr[5][0]<<" | "<< EAB[5][1] <<" "<< EABerr[5][1]<<" | "<< EAB[5][2] <<" "<< EABerr[5][2] <<" | "<<endl;
  // myfile<<"| 2.3 <abs(&eta;) < 2.4 |"<<EAB[6][0]  <<" "<<EABerr[6][0]<<" | "<< EAB[6][1] <<" "<< EABerr[6][1]<<" | "<< EAB[6][2] <<" "<<EABerr[6][2] <<" | "<<endl;
  // myfile<<"| abs(&eta;) > 2.4 |"<<EAB[7][0]       <<" "<<EABerr[7][0]<<" | "<< EAB[7][1] <<" "<< EABerr[7][1]<<" | "<< EAB[7][2] <<" "<<EABerr[7][2] <<" | "<<endl;

  myfile.close();


  ofstream myfile2;  
  myfile2.open("EA_1_60_TruePV_test.txt");
  myfile2<<EAB[0][0]<<" "<<EAB[0][1]<<" "<<EAB[0][2]<<endl;
  myfile2<<EAB[1][0]<<" "<<EAB[1][1]<<" "<<EAB[1][2]<<endl;
  myfile2<<EAB[2][0]<<" "<<EAB[2][1]<<" "<<EAB[2][2]<<endl;
  myfile2<<EAB[3][0]<<" "<<EAB[3][1]<<" "<<EAB[3][2]<<endl;
  myfile2<<EAB[4][0]<<" "<<EAB[4][1]<<" "<<EAB[4][2]<<endl;
  //    myfile2<<EAB[5][0]<<" "<<EAB[5][1]<<" "<<EAB[5][2]<<endl;
  //myfile2<<EAB[6][0]<<" "<<EAB[6][1]<<" "<<EAB[6][2]<<endl;
  //myfile2<<EAB[7][0]<<" "<<EAB[7][1]<<" "<<EAB[7][2]<<endl;
    
myfile2.close();




}
