#include <cstdlib>
#include <TMath.h>
#include <stdio.h>
#include <iostream>
using namespace std;
void Scr99(){

TFile *f1 = new TFile("../../Isopt/Mergedrun3barrel.root");
  Float_t ToE,Sie,IsoC,IsoN,Ppt,isoecal,isohcal; 
  Double_t weighT;
  ofstream myfile; 
  myfile.open("99per_bar.txt");   
  TTree *t_S = (TTree*)f1->Get("t_S");
  t_S->SetBranchAddress("Ppt",&Ppt);
  t_S->SetBranchAddress("ToE",&ToE);
  t_S->SetBranchAddress("Sieie",&Sie);
  t_S->SetBranchAddress("isoC",&IsoC);
  t_S->SetBranchAddress("isoecal",&isoecal);
  t_S->SetBranchAddress("isohcal",&isohcal);

  weighT = 1;
  float Hbin_max = 200;
  float ecal_max = 200;
  float hcal_max = 200;
  int binss = 10000;
  TH1F *HPt = new TH1F("HPt","HPt ",10000,0,1000);
  TH1F *HS = new TH1F("HS","sieie ",10000,0,0.5);
  TH1F *HH = new TH1F("HH","h over e ",10000,0,0.5);
  TH1F *HC = new TH1F("HC","worst charge  iso ",binss,0,Hbin_max);
  TH1F *ecal = new TH1F("ecal","ecal iso ",binss,0,Hbin_max);
  TH1F *hcal = new TH1F("hcal","hcal iso ",binss,0,Hbin_max);

  cout<<" Total entries in the tree  " <<t_S->GetEntries()<<endl;

  double max_s,max_c,max_t; 
  max_s = 0; 
  max_c = 0; 
  max_t = 0; 
  double totS = 1;
  for(int i = 0; i < t_S->GetEntries(); i++){
    t_S->GetEntry(i);
    if(Ppt < 200 ) continue;
    HPt->Fill(Ppt);
    if(Sie > max_s)max_s = Sie; 
    if(ToE > max_t)max_t = ToE; 
    totS += weighT;
    HH->Fill(ToE,weighT);
    HS->Fill(Sie,weighT);
    double isoch = TMath::Max(IsoC-0.0 ,0.0);
    double isoHcal = TMath::Max(isohcal -(0.0112465*Ppt+1.46174e-05*Ppt*Ppt),0.0);
    double isoEcal = TMath::Max(isoecal -(0.000583374*Ppt),0.0);
    if(IsoC> max_c)max_c =IsoC; 
    HC->Fill(IsoC,weighT);
    ecal->Fill(isoEcal,weighT);
    hcal->Fill(isoHcal,weighT);
  }

  cout<<  t_S->GetEntries() <<  " " << totS<< endl;
  
  double xcsf=0; 
  double xchf=0; 
  double xccf=0; 
  double xcecalf=0; 
  double xchcalf=0; 
  
  int p1 = 0;
  int p2 = 0;
  int p3 = 0;
  int p4 = 0;
  int p5 = 0;
  int p6 = 0;
  int p7 = 0;

  double sie_b = 0; 
  double neu_b = 0; 
  double pho_b = 0; 
  double chg_b = 0; 
  double toe_b = 0; 


  cout<<std::fixed<<(HS->Integral())<<" "<<(HS->Integral(0,10001))<<" "<<(HS->GetEntries())<<endl;
  cout<<std::fixed<<(HH->Integral())<<" "<<(HH->Integral(0,10001))<<" "<<(HH->GetEntries())<<endl;
  cout<<std::fixed<<(HC->Integral())<<" "<<(HC->Integral(0,binss+1))<<" "<<(HC->GetEntries())<<endl;  
  cout<<std::fixed<<(ecal->Integral())<<" "<<(ecal->Integral(0,binss+1))<<" "<<(ecal->GetEntries())<<endl;
  cout<<std::fixed<<(hcal->Integral())<<" "<<(hcal->Integral(0,binss+1))<<" "<<(hcal->GetEntries())<<endl;

  double HSIntegral  = HS->Integral(0,10001);
  double HHIntegral  = HH->Integral(0,10001);
  double HCIntegral  = HC->Integral(0,binss+1);
  double ecalIntegral  = ecal->Integral(0,binss+1);
  double hcalIntegral  = hcal->Integral(0,binss+1);
  
  for(int i = 1 ; i <= 10000; i++){
    double xcs = (i*0.5)/10000; 
    double xch = (i*0.5)/10000; 
    if(1.0*HS->Integral(1,i)/(double)HSIntegral > 0.9999 && p1 == 0){
      p1 = 1;
      xcsf = xcs; 
    }
    if(1.0*HH->Integral(1,i)/(double)HHIntegral > 0.999 && p2 == 0){
      xchf = xch; 
      p2 =1;
    }
}



  for(int i = 1 ; i <= binss; i++){
    double xcc = (i*Hbin_max)/binss;
    double xcecal = (i*ecal_max)/binss;
    double xchcal = (i*hcal_max)/binss;


    if(1.0*HC->Integral(1,i)/(double)HCIntegral > 0.99 && p4 == 0){
      xccf = xcc; 
      p4 = 1; 
    }
    if(1.0*ecal->Integral(1,i)/(double)ecalIntegral > 0.99 && p6 == 0){
      p6 = 1;
      xcecalf = xcecal;
    }  
    if(1.0*hcal->Integral(1,i)/(double)hcalIntegral > 0.99 && p7 == 0){
      p7 = 1;
      xchcalf = xchcal;
    }  
  } 


//SigIetaIeta HoverE Charged Ecal Hcal
myfile<<" "<<xcsf<<"  "<<xchf<<" "<<xccf<<" "<<xcecalf<<"  "<<xchcalf<<endl;



  myfile.close();
  TFile *f = new TFile("Vars_bar.root","recreate");
  HS->Write();
  HH->Write();
  HC->Write();
  HPt->Write();
  ecal->Write();
  hcal->Write();
  f->Close();

}
