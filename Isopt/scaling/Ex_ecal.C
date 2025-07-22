#include <TMath.h>
#include <cstdlib>
#include <TRandom.h>
#include <TGraph.h>
#include "vector"
#include <new>
#include "sstream"
#include <string>
#include <fstream>
#include <iostream>
void ErrCalc(TH1F *HIST,int binxn,double perc,double & X_val, double & errXL,double & errXH){
  

  TH1F *h1 = (TH1F*)HIST->Clone();
  int arsize = h1->GetXaxis()->GetNbins(); 

  double *eff; 
  double *eff_err; 
  double *cutV; 
  double *cutV_err; 
  
  eff      = new double[arsize];
  cutV     = new double[arsize];
  eff_err  = new double[arsize];
  cutV_err = new double[arsize];
  h1->Draw();
  int tot = h1->GetEntries();
  float tot_nz = tot - h1->GetBinContent(1); // copied from the barrel one not the original  
 int integ = 0; 
  for(int i  = 1; i < (h1->GetXaxis()->GetNbins() + 1); i++){ 
    double xCut = h1->GetXaxis()->GetBinLowEdge(i); 
    integ += h1->GetBinContent(i); 


    if(integ != 0 && tot != 0 ){ 
      eff[i -1] = (integ*1.0/tot);
            eff_err[i -1] = ((integ*1.0/tot)*sqrt(pow(sqrt(tot)/tot,2) + pow(sqrt(integ)/integ,2) )); 
 }
else{
      eff_err[i -1] = 0.0; 
      eff[i -1] = 0.0;
    }
      
    cutV[i -1] = (xCut);
    cutV_err[i - 1] = 0;
  }

    
  gStyle->SetOptStat(1);

  //Now draw the resulting curve 
  TCanvas *c1 = new TCanvas("c1","The Eff- cut value plot",1200,600);
    
  c1->Divide(2,1);

  c1->cd(1);
  double x[2] = {h1->GetXaxis()->GetBinCenter(1),h1->GetXaxis()->GetBinCenter(arsize)};
  double y[2] = {perc,perc};

  TGraph *lineP =  new TGraph(2,x,y); 
  lineP->SetLineColor(kRed);

  TGraphErrors *efC = new TGraphErrors(arsize ,cutV,eff,cutV_err,eff_err);
  efC->SetMarkerStyle(20);
  efC->SetMarkerSize(0.5);
    
  //   efC->Draw("AP");
  efC->GetXaxis()->SetTitle("Cut Value");
  efC->GetYaxis()->SetTitle("Efficiency");
    



  TMultiGraph *GRPhs = new TMultiGraph(); 
  GRPhs->Add(lineP,"l");
  GRPhs->Add(efC,"p");

  GRPhs->Draw("APL");

  c1->cd(2);

  h1->Draw();

  //  c1->SaveAs(graphname);   //it is block but here in barrel file we are opening it

  double * err_up; 
  double * err_down; 

  err_up = new double [arsize];
  err_down = new double [arsize];
        
  for(int i = 0;  i < arsize ; i++){
    
    if( fabs(eff[i]) > 0.00001 && fabs(eff_err[i]) > 0.0001){ // this loop is blocked also in ghe barrel one

      err_up[i] = eff[i] + eff_err[i];
      err_down[i] = eff[i] - eff_err[i];
      

    }
    else{        // this one is also block in barral one 

      
      err_up[i] = 0; 
      err_down[i] = 0; 


    }
  }
  


  //Extrapolation method to find the CutValue errors
  int up = -99; 
  int down = -99;  
  for( int i  = 0; i < arsize ; i++){
    if( err_up [i] > perc){ up = i; 
      break;
    }
  }
  for( int i  = 0; i < arsize ; i++){
    if( err_down [i] > perc){
      down = i; 
      break;
    }
  }
    
  double Usl   = 0; 
  double Ustom = 0; 
  double err_1 = 0; 
  double Dsl   = 0; 
  double Dstom = 0; 
  double err_2 = 0; 

  if(up != -99  ){
    double Usl   = (err_up[up] - err_up[up-1])/(cutV[up] - cutV[up-1]); 
    double Ustom = err_up[up] - Usl*cutV[up]; 
   
    err_1 = (perc - Ustom)/Usl; 
  }else{
    err_1 = 0; 
    Usl   = 0; 
    Ustom = 0; 
  }

  if(down != -99){
    double Dsl   = (err_down[down] - err_down[down-1])/(cutV[down] - cutV[down-1]); 
    double Dstom = err_down[down] - Dsl*cutV[down]; 
    err_2 = (perc - Dstom)/Dsl; 
  }else{
    err_2 = 0; 
    Dsl   = 0; 
    Dstom = 0; 
   
  }


  int xc = 0; 
   
  for(int i = 0; i < arsize; i++){
    if( eff[i] > perc ){ 
      xc = i; 
      break;
    }
  }
    
   

  if(fabs(eff[xc] - perc) > fabs(eff[xc - 1] - perc) ){
    X_val = cutV[xc -1 ]; 
   
  }else{
    X_val = cutV[xc ];
   
      }

  if(up == -99){
    errXL = 0 ; 


  }else{
    errXL = X_val - err_1;
  }

  if(down == -99){
      
    int dif = 0; 
    for(int i = 0; i < arsize; i++){
      if(eff[i] >  0.98){
	dif = i; 
	break; 
      }
    }
    errXH = 200 - X_val;


  }else{
    errXH = err_2 - X_val;
  }
    


  //Deleting all the dynamical arrays

    
  delete[] eff;
  delete[] cutV;

  delete[] eff_err;
  delete[] cutV_err;

  delete[] err_up; 
  delete[] err_down; 
    

}






void Ex_ecal(){

  
   TFile  *f1 =  new TFile ("../Mergedrun3barrel.root");
  float genPt,ppt,peta,Sie_ie,iso_P,iso_C,iso_N,iso_Ecal,iso_Hcal,to_e,weighT;
  int nvtx;
  gStyle->SetOptStat(0);
    TTree *t_S = (TTree*)f1->Get("t_S");
    TTree *t_B = (TTree*)f1->Get("t_B");

  //Signal Tree                                                                         
     
     // set Branch Address 


  t_S->SetBranchAddress("Sieie",&Sie_ie);
  t_S->SetBranchAddress("isoP",&iso_P);
  t_S->SetBranchAddress("isoC",&iso_C);
  t_S->SetBranchAddress("isoN",&iso_N);
  t_S->SetBranchAddress("isoecal",&iso_Ecal);
  t_S->SetBranchAddress("isohcal",&iso_Hcal);
  t_S->SetBranchAddress("ToE",&to_e);
  t_S->SetBranchAddress("weighT",&weighT);

  t_S->SetBranchAddress("Nvtx",&nvtx);
  t_S->SetBranchAddress("Peta",&peta);
  t_S->SetBranchAddress("Ppt",&ppt);
  t_S->SetBranchAddress("genPt",&genPt);

  //Background Tree                  
  t_B->SetBranchAddress("Sieie",&Sie_ie);
  t_B->SetBranchAddress("isoP",&iso_P);
  t_B->SetBranchAddress("isoC",&iso_C);
  t_B->SetBranchAddress("isoN",&iso_N);
  t_B->SetBranchAddress("isoecal",&iso_Ecal);
  t_B->SetBranchAddress("isohcal",&iso_Hcal);
  t_B->SetBranchAddress("ToE",&to_e);
  t_B->SetBranchAddress("weighT",&weighT);

  t_B->SetBranchAddress("Nvtx",&nvtx);
  t_B->SetBranchAddress("Peta",&peta);
  t_B->SetBranchAddress("Ppt",&ppt);
  
 TH2F *isoPptS = new TH2F("isoPptS","ECAL Isolation vs Pt",250,0,1000,1000,0,100);


  for(int i  = 0; i < t_S->GetEntries();i++){
    t_S->GetEntry(i);

    if(ppt < 200) continue;
    if(iso_Ecal == 0) continue;
    isoPptS->Fill(ppt,iso_Ecal);
  }

  cout<<"Builded the 2d HISTOGRAM"<<endl;



  TH2F *his2 =(TH2F*) isoPptS->Clone();



  int dim = his2->GetXaxis()->GetNbins(); 
  cout<<"dim"<<dim<<endl;  
  double * cutV; 
  double * errVH; 
  double * errVL; 
  double * binc; 
  double * bincerH; 
  double * bincerL; 
  cutV   = new double[dim];
  errVH   = new double[dim];
  errVL   = new double[dim];
  binc   = new double[dim]; 
  bincerH = new double[dim];
  bincerL = new double[dim];
  
  for(int i  = 1; i < dim ; i++){
      double xval = 0; 
      double errXH = 0; 
      double errXL = 0; 
      TH1D * r22 = his2->ProjectionY(" ",i,i+1," ");
      TH1F *h1 = (TH1F*)r22->Clone();
      cout<<"h1 get entries"<<h1->GetEntries()<<endl; 
      ErrCalc(h1,i,0.900,xval,errXL,errXH);
    cout<<"bin :"<<i<<" "<<xval<<"-"<<errXL<<"+ " << errXH<<endl;
    cutV[i-1]   = xval; 
    errVL[i-1]   = errXL; 
    errVH[i-1]   = errXH;
    binc[i-1]   = his2->GetXaxis()->GetBinCenter(i);
    bincerL[i-1] = 0;
    bincerH[i-1] = 0;
  }
  TGraphAsymmErrors * IsoptScaling = new TGraphAsymmErrors(dim,binc,cutV,bincerL,bincerH,errVL,errVH);
  TGraphAsymmErrors * IsoptScaling2 = new TGraphAsymmErrors(dim,binc,cutV,bincerL,bincerH,errVL,errVH);
  TGraphAsymmErrors * IsoptScalingLin = new TGraphAsymmErrors(dim,binc,cutV,bincerL,bincerH,errVL,errVH);

  TF1 *fn1 = new TF1("fn1","exp([0]*x + [1])",200,1000);
  TF1 *fn2 = new TF1("fn2","[1]*x + [2]*x*x + [0]",200,1000);
  TF1 *fnlin = new TF1("fnlin","[1]*x + [0]",200,1000);


  IsoptScaling->Fit("fn1","R");
  IsoptScaling2->Fit("fn2","R");
  IsoptScalingLin->Fit("fnlin","R");


  std::ofstream outfile("scale_para.txt", std::ios::app);
  if (outfile.is_open()) {
      double p1 = fnlin->GetParameter(1);
      outfile << "double ecal_scale = " << p1 << ";\n";
      outfile.close();
      std::cout << "Fitting parameters saved to scale_para.txt" << std::endl;
  }
  

  gStyle->SetOptFit(1);
  
  TCanvas *c3 = new TCanvas("c3","Iso Pt",1200,600);
  c3->Divide(2,1);
  
  c3->cd(1);
  IsoptScaling->SetMarkerStyle(24); 
  IsoptScaling->SetMarkerSize(0.4);
  IsoptScaling->Draw("AP");
  IsoptScaling->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  IsoptScaling->GetYaxis()->SetTitle("ECAL Isolation");

    IsoptScaling->GetYaxis()->SetRangeUser(0,100.0);
  IsoptScaling->GetXaxis()->SetRangeUser(0,1000);
  c3->cd(2); 
  his2->Draw("colz");
  his2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  his2->GetYaxis()->SetTitle("ECAL Isolation");
  his2->GetYaxis()->SetRangeUser(0,100);
//  c3->SaveAs("exponentialneu60.png");
//  c3->SaveAs("exponentialneu60.C");
  
  TCanvas *c4 = new TCanvas("c4","Iso Pt",1200,600);
  c4->Divide(2,1);
  
  c4->cd(1);
  IsoptScaling2->SetMarkerStyle(24); 
  IsoptScaling2->SetMarkerSize(0.4);
  IsoptScaling2->Draw("");
  IsoptScaling2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  IsoptScaling2->GetYaxis()->SetTitle("ECAL Isolation");
  IsoptScaling2->GetYaxis()->SetRangeUser(0,35.0);
  IsoptScaling2->GetYaxis()->SetRangeUser(0,100);
  c4->cd(2);
  his2->Draw("colz");
  his2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  his2->GetYaxis()->SetTitle("ECAL Isolation");
  his2->GetYaxis()->SetRangeUser(0,50);
//  c4->SaveAs("Quadratic1to60_neu.png");
//  c4->SaveAs("Quadratic1to60_neu.C");

  TCanvas *c6 = new TCanvas("c6","Iso Pt",1200,600);
  c6->Divide(2,1);  
  c6->cd(1);
  IsoptScalingLin->SetMarkerStyle(24); 
  IsoptScalingLin->SetMarkerSize(0.4);
  IsoptScalingLin->Draw("AP");
  IsoptScalingLin->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  IsoptScalingLin->GetYaxis()->SetTitle("ECAL Isolation ");

  c6->cd(2);
  his2->Draw("colz");
  his2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  his2->GetYaxis()->SetTitle("ECAL Isolation");
  c6->SaveAs("Straightline_Ecal.png");
  c6->SaveAs("Straightline_Ecal.C");

  TCanvas *c7 = new TCanvas("c7","Iso Pt",1200,600);

  isoPptS->Draw("colz");
  c7->SaveAs("isoPptS_neu.png");

}


