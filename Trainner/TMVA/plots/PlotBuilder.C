#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "input_params.h"
#include "TH1.h"
#include <sstream>
#include <string.h>

void PlotBuilder(){

  //Read the cuts: 


  ifstream tight;
  ifstream medium;
  ifstream loose;

  double isoCL,isoCM,isoCT;
  double isoecalL,isoecalM,isoecalT;
  double isoehalL,isohcalM,isohcalT;
  double isoPL,isoPM,isoPT;
  double isoNL,isoNM,isoNT;
  double sieL,sieM,sieT; 
  double hoeL,hoeM,hoeT;

Double_t bins_pt[] = { 0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,450,500,550,600,1000};
int nbinss = 24;
  TH1F *EffNVTX0 = new TH1F("EffNVTX0","Signal Eff vs NVTX 0",100,0,100);
  TH1F *EffPT0   = new TH1F("EffPT0","Signal Eff vs PT 0",nbinss,bins_pt);
  TH1F *EffETA0  = new TH1F("EffETA0","Signal Eff vs Eta 0",100,-5,5);
  TH1F *EffPHI0  = new TH1F("EffPHI0","Signal Eff vs PHI 0",100,-4,4);

  
  TH1F *EffNVTXM = new TH1F("EffNVTXM","Signal Eff vs NVTX M",100,0,100);
  TH1F *EffPTM   = new TH1F("EffPTM","Signal Eff vs PT M",nbinss,bins_pt);
  TH1F *EffETAM  = new TH1F("EffETAM","Signal Eff vs Eta M",100,-5,5);
  TH1F *EffPHIM  = new TH1F("EffPHIM","Signal Eff vs PHI M",100,-4,4);


  TH1F *EffNVTX0b = new TH1F("EffNVTX0b","Background Eff vs NVTX 0",100,0,100);
  TH1F *EffPT0b   = new TH1F("EffPT0b","Background Eff vs PT 0",nbinss,bins_pt);
  TH1F *EffETA0b  = new TH1F("EffETA0b","Background Eff vs Eta 0",100,-5,5);
  TH1F *EffPHI0b  = new TH1F("EffPHI0b","Background Eff vs PHI 0",100,-4,4);

  
  TH1F *EffNVTXMb = new TH1F("EffNVTXMb","Background Eff vs NVTX M",100,0,100);
  TH1F *EffPTMb   = new TH1F("EffPTMb","Background Eff vs PT M",nbinss,bins_pt);
  TH1F *EffETAMb  = new TH1F("EffETAMb","Background Eff vs Eta M",100,-5,5);
  TH1F *EffPHIMb  = new TH1F("EffPHIMb","Background Eff vs PHI M",100,-4,4);


  // Branch out Cuts

   

  TH1F *EffNVTXt = new TH1F("EffNVTXt","Signal Eff vs NVTX T s",100,0,100);
  TH1F *EffPTt  = new TH1F("EffPTt","Signal Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAt  = new TH1F("EffETAt","Signal Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXc = new TH1F("EffNVTXc","Signal Eff vs NVTX T s",100,0,100);
  TH1F *EffPTc  = new TH1F("EffPTc","Signal Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAc  = new TH1F("EffETAc","Signal Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXn = new TH1F("EffNVTXn","Signal Eff vs NVTX T s",100,0,100);
  TH1F *EffPTn  = new TH1F("EffPTn","Signal Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAn  = new TH1F("EffETAn","Signal Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXp = new TH1F("EffNVTXp","Signal Eff vs NVTX T s",100,0,100);
  TH1F *EffPTp  = new TH1F("EffPTp","Signal Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAp  = new TH1F("EffETAp","Signal Eff vs PT T s",100,-5,5);


  TH1F *EffNVTXbs = new TH1F("EffNVTXbs","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTbs  = new TH1F("EffPTbs","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAbs  = new TH1F("EffETAbs","Background Eff vs PT T s",100,-5,5);


  TH1F *EffNVTXbc = new TH1F("EffNVTXbc","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTbc  = new TH1F("EffPTbc","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAbc  = new TH1F("EffETAbc","Background Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXbn = new TH1F("EffNVTXbn","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTbn  = new TH1F("EffPTbn","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAbn  = new TH1F("EffETAbn","Background Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXbp = new TH1F("EffNVTXbp","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTbp  = new TH1F("EffPTbp","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAbp  = new TH1F("EffETAbp","Background Eff vs PT T s",100,-5,5);

  //PT
  TH1F *EffPTMs  = new TH1F("EffPTMs","Signal Eff vs PT s",nbinss,bins_pt);
  TH1F *EffPTMt  = new TH1F("EffPTMt","Signal Eff vs PT s",nbinss,bins_pt);
  TH1F *EffPTMc  = new TH1F("EffPTMc","Signal Eff vs PT s",nbinss,bins_pt);
  TH1F *EffPTMp  = new TH1F("EffPTMp","Signal Eff vs PT s",nbinss,bins_pt);
  TH1F *EffPTMn  = new TH1F("EffPTMn","Signal Eff vs PT s",nbinss,bins_pt);



  //NVTX
  TH1F *EffNVTXMs  = new TH1F("EffNVTXMs","Signal Eff vs NVTX s",100,0,100);
  TH1F *EffNVTXMt  = new TH1F("EffNVTXMt","Signal Eff vs NVTX s",100,0,100);
  TH1F *EffNVTXMc  = new TH1F("EffNVTXMc","Signal Eff vs NVTX s",100,0,100);
  TH1F *EffNVTXMp  = new TH1F("EffNVTXMp","Signal Eff vs NVTX s",100,0,100);
  TH1F *EffNVTXMn  = new TH1F("EffNVTXMn","Signal Eff vs NVTX s",100,0,100);



  //ETA
  TH1F *EffETAMs  = new TH1F("EffETAMs","Signal Eff vs ETA s",100,-5,5);
  TH1F *EffETAMt  = new TH1F("EffETAMt","Signal Eff vs ETA s",100,-5,5);
  TH1F *EffETAMc  = new TH1F("EffETAMc","Signal Eff vs ETA s",100,-5,5);
  TH1F *EffETAMp  = new TH1F("EffETAMp","Signal Eff vs ETA s",100,-5,5);
  TH1F *EffETAMn  = new TH1F("EffETAMn","Signal Eff vs ETA s",100,-5,5);




  //Medium 
  TH1F *EffNVTXMbs = new TH1F("EffNVTXMbs","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTMbs  = new TH1F("EffPTMbs","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAMbs  = new TH1F("EffETAMbs","Background Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXMbt = new TH1F("EffNVTXMbt","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTMbt  = new TH1F("EffPTMbt","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAMbt  = new TH1F("EffETAMbt","Background Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXMbc = new TH1F("EffNVTXMbc","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTMbc  = new TH1F("EffPTMbc","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAMbc  = new TH1F("EffETAMbc","Background Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXMbn = new TH1F("EffNVTXMbn","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTMbn  = new TH1F("EffPTMbn","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAMbn  = new TH1F("EffETAMbn","Background Eff vs PT T s",100,-5,5);

  TH1F *EffNVTXMbp = new TH1F("EffNVTXMbp","Background Eff vs NVTX T s",100,0,100);
  TH1F *EffPTMbp  = new TH1F("EffPTMbp","Background Eff vs PT T s",nbinss,bins_pt);
  TH1F *EffETAMbp  = new TH1F("EffETAMbp","Background Eff vs PT T s",100,-5,5);






  
  //Setting the Tree Branches

   TFile *finput = new TFile( "../../../Isopt/Mergedrun3barrel.root");
  float Ppt,Peta,Pphi,isoP,isoecal,isohcal,isoC,isoN,sieie,toe,weight;
  int   nvtx; 
  float myl_hoe,myl_sie,myl_c,myl_e,myl_h; 
  float mym_hoe,mym_sie,mym_c,mym_e,mym_h; 
  float myt_hoe,myt_sie,myt_c,myt_e,myt_h; 
  //  weight = 1.0; 

  finput->cd();
  TTree *t_S = (TTree*)finput->Get("t_S");
 TTree *t_B = (TTree*)finput->Get("t_B");

  //Signal Tree                                                                 
  t_S->SetBranchAddress("Sieie",&sieie);
  t_S->SetBranchAddress("isoP",&isoP);
  t_S->SetBranchAddress("isoecal",&isoecal);
  t_S->SetBranchAddress("isohcal",&isohcal);
  t_S->SetBranchAddress("isoC",&isoC);
  t_S->SetBranchAddress("isoN",&isoN);
  t_S->SetBranchAddress("ToE",&toe);
  t_S->SetBranchAddress("weighT",&weight);
  t_S->SetBranchAddress("Nvtx",&nvtx);
  t_S->SetBranchAddress("Peta",&Peta);
  t_S->SetBranchAddress("Ppt",&Ppt);

  //Background Tree                                                                 
  t_B->SetBranchAddress("Sieie",&sieie);
  t_B->SetBranchAddress("isoP",&isoP);
  t_B->SetBranchAddress("isoecal",&isoecal);
  t_B->SetBranchAddress("isohcal",&isohcal);
  t_B->SetBranchAddress("isoC",&isoC);
  t_B->SetBranchAddress("isoN",&isoN);
  t_B->SetBranchAddress("ToE",&toe);
  t_B->SetBranchAddress("weighT",&weight);
  t_B->SetBranchAddress("Nvtx",&nvtx);
  t_B->SetBranchAddress("Peta",&Peta);
  t_B->SetBranchAddress("Ppt",&Ppt);

  
  EffNVTX0->Sumw2();
  EffPT0->Sumw2();
  EffETA0->Sumw2();
  EffPHI0->Sumw2();
  
  
  EffNVTXM->Sumw2();
  EffPTM->Sumw2();
  EffETAM->Sumw2();
  EffPHIM->Sumw2();
  
  
  
  EffNVTX0b->Sumw2();
  EffPT0b->Sumw2();
  EffETA0b->Sumw2();
  EffPHI0b->Sumw2();
  
  
  EffNVTXMb->Sumw2();
  EffPTMb->Sumw2();
  EffETAMb->Sumw2();
  EffPHIMb->Sumw2();
  

  
  
  EffNVTXt->Sumw2();
  EffPTt->Sumw2();
  EffETAt->Sumw2();
  
  EffNVTXc->Sumw2();
  EffPTc->Sumw2();
  EffETAc->Sumw2();
  
  EffNVTXn->Sumw2();
  EffPTn->Sumw2();
  EffETAn->Sumw2();

  EffNVTXp->Sumw2();
  EffPTp->Sumw2();
  EffETAp->Sumw2();

  
  EffNVTXbs->Sumw2();
  EffPTbs->Sumw2();
  EffETAbs->Sumw2();


  EffNVTXbc->Sumw2();
  EffPTbc->Sumw2();
  EffETAbc->Sumw2();

  EffNVTXbn->Sumw2();
  EffPTbn->Sumw2();
  EffETAbn->Sumw2();

  EffNVTXbp->Sumw2();
  EffPTbp->Sumw2();
  EffETAbp->Sumw2();

  //PT    
  EffPTMs->Sumw2();
  EffPTMt->Sumw2();
  EffPTMc->Sumw2();
  EffPTMp->Sumw2();
  EffPTMn->Sumw2();
    

  //NVTX
  EffNVTXMs->Sumw2();
  EffNVTXMt->Sumw2();
  EffNVTXMc->Sumw2();
  EffNVTXMp->Sumw2();
  EffNVTXMn->Sumw2();
    

  //ETA
  EffETAMs->Sumw2();
  EffETAMt->Sumw2();
  EffETAMc->Sumw2();
  EffETAMp->Sumw2();
  EffETAMn->Sumw2();
    

  //Medium
  EffNVTXMbs->Sumw2();
  EffPTMbs->Sumw2();
  EffETAMbs->Sumw2();

  EffNVTXMbt->Sumw2();
  EffPTMbt->Sumw2();
  EffETAMbt->Sumw2();

  EffNVTXMbc->Sumw2();
  EffPTMbc->Sumw2();
  EffETAMbc->Sumw2();

  EffNVTXMbn->Sumw2();
  EffPTMbn->Sumw2();
  EffETAMbn->Sumw2();

  EffNVTXMbp->Sumw2();
  EffPTMbp->Sumw2();
  EffETAMbp->Sumw2();


double epscale = 0.000583374;
double hpscale = 0.0112465;
double hppscale = 1.46174e-5;
double medium_ept = epscale;
double medium_hpt = hpscale;
double medium_hptpt = hppscale;


Long64_t snentries = t_S->GetEntries();


for(int i = 0; i < t_S->GetEntries(); i++){
    t_S->GetEntry(i);
//  for(int i = 0; i < t_B->GetEntries(); i++){
          if (i % (snentries / 20) == 0) {
        std::cout << "Signal [Progress] " << i << " / " << snentries << " events processed..." << std::endl;
          }
        //    t_B->GetEntry(i);
    
    EffNVTX0->Fill(nvtx,weight);
    EffPT0->Fill(Ppt,weight);
    EffETA0->Fill(Peta,weight);

    //Medium Cut:
    if((sieie  < medium_sieie )&&
       (toe    < medium_toe )&&
       (isoC < medium_isoC)&&
       (isoecal < medium_ecal + (medium_ept*Ppt))&&
       (isohcal < medium_hcal + medium_hpt*Ppt+medium_hptpt*Ppt*Ppt)
       ){
      EffNVTXM->Fill(nvtx,weight);
      EffPTM->Fill(Ppt,weight);
      EffETAM->Fill(Peta,weight);
    }
    
  //Branch out Cuts

    //Branch out Cuts Medium
    if(sieie  < medium_sieie){
      EffPTMs->Fill(Ppt,weight);
      EffNVTXMs->Fill(nvtx,weight);
      EffETAMs->Fill(Peta,weight);
    }
    if(toe    <medium_toe ){
      EffPTMt->Fill(Ppt,weight);
      EffNVTXMt->Fill(nvtx,weight);
      EffETAMt->Fill(Peta,weight);
    }    
    if(isoecal < medium_ecal + (medium_ept*Ppt)){    
      EffPTMp->Fill(Ppt,weight);
      EffNVTXMp->Fill(nvtx,weight);
      EffETAMp->Fill(Peta,weight);
    }    
    if(isoC   < medium_isoC ){
      EffPTMc->Fill(Ppt,weight);
      EffNVTXMc->Fill(nvtx,weight);
      EffETAMc->Fill(Peta,weight);
    }
    if(isohcal < medium_hcal + medium_hpt*Ppt+medium_hptpt*Ppt*Ppt){
      EffPTMn->Fill(Ppt,weight);
      EffNVTXMn->Fill(nvtx,weight);
      EffETAMn->Fill(Peta,weight);
    }


    



  }  


  //plots for Background 

Long64_t nentries = t_B->GetEntries();

  for(int i = 0; i < t_B->GetEntries(); i++){
          if (i % (nentries / 20) == 0) {
        std::cout << "Background [Progress] " << i << " / " << nentries << " events processed..." << std::endl;
          }
    t_B->GetEntry(i);
    EffNVTX0b->Fill(nvtx,weight);
    EffPT0b->Fill(Ppt,weight);
    EffETA0b->Fill(Peta,weight);
  
    
       //Medium Cut:
    if((sieie  < medium_sieie )&&
       (toe    < medium_toe )&&
       (isoC < medium_isoC)&&
       (isoecal < medium_ecal + (medium_ept*Ppt))&&
       (isohcal < medium_hcal + medium_hpt*Ppt+medium_hptpt*Ppt*Ppt)
       ){
	 EffNVTXMb->Fill(nvtx,weight);
	 EffPTMb->Fill(Ppt,weight);
	 EffETAMb->Fill(Peta,weight);
       }
        //Branch out Cuts

	  //Medium
	  if(sieie  < medium_sieie){
	    EffNVTXMbs->Fill(nvtx,weight);
	    EffPTMbs->Fill(Ppt,weight);
	    EffETAMbs->Fill(Peta,weight);
	  }
	  if(toe    < medium_toe ){
	    EffNVTXMbt->Fill(nvtx,weight);
	    EffPTMbt->Fill(Ppt,weight);
	    EffETAMbt->Fill(Peta,weight);
	  }    
	  if(isoecal < medium_ecal + (medium_ept*Ppt)){    
	    EffNVTXMbp->Fill(nvtx,weight);
	    EffPTMbp->Fill(Ppt,weight);
	    EffETAMbp->Fill(Peta,weight);
	  }    
	  if(isoC   < medium_isoC ){
	    EffNVTXMbc->Fill(nvtx,weight);
	    EffPTMbc->Fill(Ppt,weight);
	    EffETAMbc->Fill(Peta,weight);
 
	  }
	  if(isohcal <  medium_hcal + medium_hpt*Ppt+medium_hptpt*Ppt*Ppt ){
	    EffNVTXMbn->Fill(nvtx,weight);
	    EffPTMbn->Fill(Ppt,weight);
	    EffETAMbn->Fill(Peta,weight);
      
	  }
	  }



       TFile *f1 = new TFile("Eff_updat_SR_TruePVID.root","recreate");
       f1->cd();
       EffNVTX0->Write();
       EffPT0->Write();
       EffETA0->Write();
       EffPHI0->Write();
  
  
       EffNVTXM->Write();
       EffPTM->Write();
       EffETAM->Write();
       EffPHIM->Write();
  


       EffNVTX0b->Write();
       EffPT0b->Write();
       EffETA0b->Write();
       EffPHI0b->Write();
  
  
       EffNVTXMb->Write();
       EffPTMb->Write();
       EffETAMb->Write();
       EffPHIMb->Write();
  




       EffNVTXt->Write();
       EffPTt->Write();
       EffETAt->Write();
    
       EffNVTXc->Write();
       EffPTc->Write();
       EffETAc->Write();


       EffNVTXn->Write();
       EffPTn->Write();
       EffETAn->Write();
  

       EffNVTXp->Write();
       EffPTp->Write();
       EffETAp->Write();

  
       EffNVTXbs->Write();
       EffPTbs->Write();
       EffETAbs->Write();

  
  
       EffNVTXbc->Write();
       EffPTbc->Write();
       EffETAbc->Write();

       EffNVTXbn->Write();
       EffPTbn->Write();
       EffETAbn->Write();

       EffNVTXbp->Write();
       EffPTbp->Write(); 
       EffETAbp->Write();


       //PT
       EffPTMs->Write();
       EffPTMt->Write();
       EffPTMc->Write();
       EffPTMp->Write();
       EffPTMn->Write();
    

       //NVTX
       EffNVTXMs->Write();
       EffNVTXMt->Write();
       EffNVTXMc->Write();
       EffNVTXMp->Write();
       EffNVTXMn->Write();
    

       //ETA
       EffETAMs->Write();
       EffETAMt->Write();
       EffETAMc->Write();
       EffETAMp->Write();
       EffETAMn->Write();
    
  
       //Medium
       EffNVTXMbs->Write();
       EffPTMbs->Write();
       EffETAMbs->Write();
  
       EffNVTXMbt->Write();
       EffPTMbt->Write();
       EffETAMbt->Write();
  
       EffNVTXMbc->Write();
       EffPTMbc->Write();
       EffETAMbc->Write();

       EffNVTXMbn->Write();
       EffPTMbn->Write();
       EffETAMbn->Write();

       EffNVTXMbp->Write();
       EffPTMbp->Write(); 
       EffETAMbp->Write();




       TFile *feta = new TFile("Eff1etaB_updat_SR_TruePVID.root","recreate");
       feta->cd();


       EffETA0->Write();
       EffETA0b->Write();





 
       // Plotter.C //
       TH1F *Metas = new TH1F("Metas","Medium Cut TruePV Efficiency Eta",100,-5,5);
       TH1F *Mpts = new TH1F("Mpts","Medium Cut TruePV Efficiency pt",nbinss,bins_pt);
       TH1F *Mnvtxs = new TH1F("Mnvtxs","Medium Cut TruePV Efficiency vertices",100,0,100);
       TH1F *Metab = new TH1F("Metab","Medium Cut b TruePV Efficiency Eta",100,-5,5);
       TH1F *Mptb = new TH1F("Mptb","Medium Cut b TruePV Efficiency pt",nbinss,bins_pt);
       TH1F *Mnvtxb = new TH1F("Mnvtxb","Medium Cut b TruePV Efficiency vertices",100,0,100);





       TH1F *SieaftM  = new TH1F("SieaftM","Sieie cut only",nbinss,bins_pt); 
       TH1F *ToEaftM  = new TH1F("ToEaftM","HoE cut only",nbinss,bins_pt); 
       TH1F *IsoPaftM = new TH1F("IsoPaftM","IsoEcal cut only",nbinss,bins_pt); 
       TH1F *IsoCaftM = new TH1F("IsoCaftM","IsoC cut only",nbinss,bins_pt); 
       TH1F *IsoNaftM = new TH1F("IsoNaftM","IsoHcal cut only",nbinss,bins_pt); 




       TH1F *SieaftbM  = new TH1F("Sieaftb","Sieie cut only",nbinss,bins_pt); 
       TH1F *ToEaftbM  = new TH1F("ToEaftb","HoE cut only",nbinss,bins_pt); 
       TH1F *IsoPaftbM = new TH1F("IsoPaftb","IsoEcal cut only",nbinss,bins_pt); 
       TH1F *IsoCaftbM = new TH1F("IsoCaftb","IsoC cut only",nbinss,bins_pt); 
       TH1F *IsoNaftbM = new TH1F("IsoNaftb","IsoHcal cut only",nbinss,bins_pt); 



  










       Metas->Divide(EffETAM,EffETA0,1.,1.,"B");
       Metab->Divide(EffETAMb,EffETA0b,1.,1.,"B");
       Mpts->Divide(EffPTM,EffPT0,1.,1.,"B");
       Mptb->Divide(EffPTMb,EffPT0b,1.,1.,"B");
       Mnvtxs->Divide(EffNVTXM,EffNVTX0,1.,1.,"B");
       Mnvtxb->Divide(EffNVTXMb,EffNVTX0b,1.,1.,"B");


       // the branch  out cuts 



       SieaftM->Divide(EffPTMs,EffPT0,1.,1.,"B"); 
       ToEaftM->Divide(EffPTMt,EffPT0,1.,1.,"B");
       IsoPaftM->Divide(EffPTMp,EffPT0,1.,1.,"B");
       IsoCaftM->Divide(EffPTMc,EffPT0,1.,1.,"B");
       IsoNaftM->Divide(EffPTMn,EffPT0,1.,1.,"B");



       SieaftbM->Divide(EffPTMbs,EffPT0b,1.,1.,"B"); 
       ToEaftbM->Divide(EffPTMbt,EffPT0b,1.,1.,"B");
       IsoPaftbM->Divide(EffPTMbp,EffPT0b,1.,1.,"B");
       IsoCaftbM->Divide(EffPTMbc,EffPT0b,1.,1.,"B");
       IsoNaftbM->Divide(EffPTMbn,EffPT0b,1.,1.,"B");




       TCanvas *c1  = new TCanvas("c1","Medium",600,600);
  c1->Divide(2,2);

  c1->cd(1);
  Mnvtxs->Draw();
  Mnvtxs->GetYaxis()->SetRangeUser(0,1.0);
 
  Mnvtxs->GetXaxis()->SetRangeUser(0,50);
 
  Mnvtxs->GetXaxis()->SetTitle("# Nvtx");
  Mnvtxs->SetLineColor(kRed);
  Mnvtxs->SetMarkerColor(kRed);
  Mnvtxs->SetMarkerSize(0.5);
  Mnvtxb->SetLineColor(kGreen);
  Mnvtxb->SetMarkerColor(kGreen);
  Mnvtxb->SetMarkerSize(0.5);
  Mnvtxb->Draw("same");
  
  c1->cd(2);
  Mpts->Draw();
  Mpts->GetYaxis()->SetRangeUser(0,1.0);
 
  Mpts->GetXaxis()->SetTitle("Pt GeVc^{-1}");
  Mpts->SetLineColor(kRed);
  Mpts->SetMarkerColor(kRed);
  Mpts->SetMarkerSize(0.5);
  Mptb->SetLineColor(kGreen);
  Mptb->SetMarkerColor(kGreen);
  Mptb->SetMarkerSize(0.5);
  Mptb->Draw("same");

  c1->cd(3);
  Metas->Draw();
  Metas->GetYaxis()->SetRangeUser(0,1.0);
  // Metas->GetXaxis()->SetRangeUser(1.5,3.0);
  Metas->GetXaxis()->SetRangeUser(-5.0,5.0);
 
  Metas->GetXaxis()->SetTitle("#eta");
  Metas->SetLineColor(kRed);
  Metas->SetMarkerColor(kRed);
  Metas->SetMarkerSize(0.5);
  Metab->SetLineColor(kGreen);
  Metab->SetMarkerColor(kGreen);
  Metab->SetMarkerSize(0.5);
  Metab->Draw("same");
//  c1->SaveAs("MediumEffBck_TruePV.png");
  c1->SaveAs(output_filename);






 



 
       TCanvas *c20 = new TCanvas("c20","Branch Out Cuts",900,600);
       c20->Divide(3,2);
  
       c20->cd(1);  

       SieaftM->Draw();
       SieaftM->SetMarkerSize(0.5);
 
       SieaftM->SetMarkerColor(kOrange -3);
       SieaftM->GetYaxis()->SetTitle("Only Sieie Cut TruePV Efficiency");
       SieaftM->GetXaxis()->SetTitle("Pt GeVc^{-1}");
 
       c20->cd(2);
  
       ToEaftM->Draw();
       ToEaftM->SetMarkerSize(0.5);
 
       ToEaftM->SetMarkerColor(kOrange -3);
       ToEaftM->GetYaxis()->SetTitle("Only H over E Cut TruePV Efficiency");
       ToEaftM->GetXaxis()->SetTitle("Pt GeVc^{-1}");

       c20->cd(3);

       IsoPaftM->Draw();
       IsoPaftM->SetMarkerSize(0.5);
 
       IsoPaftM->SetMarkerColor(kOrange -3);
       IsoPaftM->GetYaxis()->SetTitle("Only Iso Ecal Cut TruePV Efficiency");
       IsoPaftM->GetXaxis()->SetTitle("Pt GeVc^{-1}");

       c20->cd(4);

       IsoCaftM->Draw();
       IsoCaftM->SetMarkerSize(0.5);
       IsoCaftM->SetMarkerSize(0.5);
 
       IsoCaftM->SetMarkerColor(kOrange -3);
       IsoCaftM->GetYaxis()->SetTitle("Only Iso c Cut TruePV Efficiency");
       IsoCaftM->GetXaxis()->SetTitle("Pt GeVc^{-1}");
  
       c20->cd(5);
  
       IsoNaftM->Draw();
       IsoNaftM->SetMarkerSize(0.5);
 
       IsoNaftM->SetMarkerColor(kOrange -3);
       IsoNaftM->GetYaxis()->SetTitle("Only Iso Hcal Cut TruePV Efficiency");
       IsoNaftM->GetXaxis()->SetTitle("Pt GeVc^{-1}");

       c20->SaveAs("BranchOutCutsSignalWP_TruePV.png");

 

 









  
       TCanvas *c11 = new TCanvas("c11","Branch Out Cuts",900,600);
       c11->Divide(3,2);
  
       c11->cd(1);  

       SieaftM->Draw();
       SieaftbM->Draw("same");
       SieaftM->SetMarkerSize(0.5);
       SieaftbM->SetMarkerSize(0.5);
       SieaftbM->SetMarkerColor(kRed);
       SieaftM->GetYaxis()->SetTitle("Only Sieie Cut TruePV Efficiency");
       SieaftM->GetXaxis()->SetTitle("# Nvtx");
 
       c11->cd(2);
       ToEaftM->Draw();
       ToEaftbM->Draw("same");
       ToEaftM->SetMarkerSize(0.5);
       ToEaftbM->SetMarkerSize(0.5);
       ToEaftbM->SetMarkerColor(kRed);
       ToEaftM->GetYaxis()->SetTitle("Only HOE Cut TruePV Efficiency");
       ToEaftM->GetXaxis()->SetTitle("# Nvtx");

       c11->cd(3);
       IsoPaftM->Draw();
       IsoPaftbM->Draw("same");
       IsoPaftM->SetMarkerSize(0.5);
       IsoPaftbM->SetMarkerColor(kRed);
       IsoPaftbM->SetMarkerSize(0.5);
       IsoPaftM->GetYaxis()->SetTitle("Only iso Ecal Cut TruePV Efficiency");
       IsoPaftM->GetXaxis()->SetTitle("# Nvtx");

       c11->cd(4);
       IsoCaftM->Draw();
       IsoCaftbM->Draw("same");
       IsoCaftM->SetMarkerSize(0.5);
       IsoCaftbM->SetMarkerColor(kRed);
       IsoCaftbM->SetMarkerSize(0.5);
       IsoCaftM->GetYaxis()->SetTitle("Only iso c Cut TruePV Efficiency");
       IsoCaftM->GetXaxis()->SetTitle("# Nvtx");
  
       c11->cd(5);
       IsoNaftM->Draw();
       IsoNaftbM->Draw("same");
       IsoNaftM->SetMarkerSize(0.5);
       IsoNaftbM->SetMarkerSize(0.5);
       IsoNaftbM->SetMarkerColor(kRed);
       IsoNaftM->GetYaxis()->SetTitle("Only iso Hcal Cut TruePV Efficiency");
       IsoNaftM->GetXaxis()->SetTitle("# Nvtx");
  
       c11->SaveAs("BranchOutCutsNVTX_TruePV.png");




       TCanvas *c12 = new TCanvas("c12","Branch Out Cuts",900,600);
       c12->Divide(3,2);
  
       c12->cd(1);  

       SieaftM->Draw();
       SieaftbM->Draw("same");
       SieaftM->SetMarkerSize(0.5);
       SieaftbM->SetMarkerSize(0.5);
       SieaftbM->SetMarkerColor(kRed);
       SieaftM->GetYaxis()->SetTitle("Only Sieie Cut TruePV Efficiency");
       SieaftM->GetXaxis()->SetTitle("#eta");
 
       c12->cd(2);
       ToEaftM->Draw();
       ToEaftbM->Draw("same");
       ToEaftM->SetMarkerSize(0.5);
       ToEaftbM->SetMarkerSize(0.5);
       ToEaftbM->SetMarkerColor(kRed);
       ToEaftM->GetYaxis()->SetTitle("Only HOE Cut TruePV Efficiency");
       ToEaftM->GetXaxis()->SetTitle("#eta");

       c12->cd(3);
       IsoPaftM->Draw();
       IsoPaftbM->Draw("same");
       IsoPaftM->SetMarkerSize(0.5);
       IsoPaftbM->SetMarkerColor(kRed);
       IsoPaftbM->SetMarkerSize(0.5);
       IsoPaftM->GetYaxis()->SetTitle("Only iso ecal Cut TruePV Efficiency");
       IsoPaftM->GetXaxis()->SetTitle("#eta");

       c12->cd(4);
       IsoCaftM->Draw();
       IsoCaftbM->Draw("same");
       IsoCaftM->SetMarkerSize(0.5);
       IsoCaftbM->SetMarkerColor(kRed);
       IsoCaftbM->SetMarkerSize(0.5);
       IsoCaftM->GetYaxis()->SetTitle("Only iso c Cut TruePV Efficiency");
       IsoCaftM->GetXaxis()->SetTitle("#eta");
  
       c12->cd(5);
       IsoNaftM->Draw();
       IsoNaftbM->Draw("same");
       IsoNaftM->SetMarkerSize(0.5);
       IsoNaftbM->SetMarkerSize(0.5);
       IsoNaftbM->SetMarkerColor(kRed);
       IsoNaftM->GetYaxis()->SetTitle("Only iso hcaln Cut TruePV Efficiency");
       IsoNaftM->GetXaxis()->SetTitle("#eta");
  
       c12->SaveAs("BranchOutCutsEta_TruePV.png");




       TCanvas *cpt = new TCanvas("cpt","Pt Eff",500,500);
       cpt->cd();

       Mpts->Draw();
       Mpts->SetMarkerColor(kOrange -3);
       Mpts->SetLineColor(kOrange -3);
       Mpts->SetMarkerStyle(20);
       Mpts->SetMarkerSize(0.5);
 

  
       Mptb->SetLineColor(kAzure + 5);
       Mptb->SetMarkerColor(kAzure +5);
       Mptb->SetMarkerStyle(20);
       Mptb->SetMarkerSize(0.5);
       Mptb->Draw("esame");

 

       cpt->SaveAs("EfPT_TruePV.png");

       TCanvas *ceta = new TCanvas("ceta","Eta Eff",500,500);
       ceta->cd();
  
       Metas->SetMarkerColor(kOrange -3);
       Metas->SetLineColor(kOrange -3);
       Metas->SetMarkerStyle(20);
       Metas->SetMarkerSize(0.5);
       Metas->Draw();




       Metab->SetLineColor(kAzure + 5);
       Metab->SetMarkerColor(kAzure +5);
       Metab->SetMarkerStyle(20);
       Metab->SetMarkerSize(0.5);
       Metab->Draw("esame");


      
       ceta->SaveAs("EfETA_TruePV.png");


  
       TCanvas *cnvtx = new TCanvas("cnvtx","NVtx Eff",500,500);
       cnvtx->cd();

       Mnvtxs->Draw();
       Mnvtxs->SetMarkerColor(kOrange -3);
       Mnvtxs->SetLineColor(kOrange -3);
       Mnvtxs->SetMarkerStyle(20);
       Mnvtxs->SetMarkerSize(0.5);

  
       Mnvtxb->SetLineColor(kAzure + 5);
       Mnvtxb->SetMarkerColor(kAzure +5);
       Mnvtxb->SetMarkerStyle(20);
       Mnvtxb->SetMarkerSize(0.5);
       Mnvtxb->Draw("esame");


       cnvtx->Update();
       cnvtx->SaveAs("EfNVTX_TruePV.png");


       // sig-eff/bkg-eff plots
       TCanvas *crpt = new TCanvas("crpt","Pt SBR Eff",500,500);
       crpt->cd();

       TH1F *Mptr = (TH1F*)Mpts->Clone();
       Mptr->Divide(Mptb);
       Mptr->Draw();
       Mptr->SetMarkerColor(kOrange -3);
       Mptr->SetLineColor(kOrange -3);
       Mptr->SetMarkerStyle(20);
       Mptr->SetMarkerSize(0.5);

  
       auto lptr = new TLegend(0.55,0.15,0.9,0.3);
       lptr->SetBorderSize(0);
       lptr->SetFillColor(0);
       lptr->AddEntry(Mptr,"Medium","lp");
       lptr->Draw("same");     


       crpt->SaveAs("EfSBR_PT_TruePV.png");




       TCanvas *creta = new TCanvas("creta","Eta SBR Eff",500,500);
       creta->cd();

       TH1F *Metar = (TH1F*)Metas->Clone();
       Metar->Divide(Metab);
       Metar->Draw();
       Metar->SetMarkerColor(kOrange -3);
       Metar->SetLineColor(kOrange -3);
       Metar->SetMarkerStyle(20);
       Metar->SetMarkerSize(0.5);
 
  
       auto letar = new TLegend(0.55,0.15,0.9,0.3);
       letar->SetBorderSize(0);
       letar->SetFillColor(0);
       letar->AddEntry(Metar,"Medium","lp");
       letar->Draw("same");     


       creta->SaveAs("EfSBR_ETA_TruePV.png");




       TCanvas *crnvtx = new TCanvas("crnvtx","Nvtx SBR Eff",500,500);
       crnvtx->cd();

       TH1F *Mnvtxr = (TH1F*)Mnvtxs->Clone();
       Mnvtxr->Divide(Mnvtxb);
       Mnvtxr->Draw();
       Mnvtxr->SetMarkerColor(kOrange -3);
       Mnvtxr->SetLineColor(kOrange -3);
       Mnvtxr->SetMarkerStyle(20);
       Mnvtxr->SetMarkerSize(0.5);
 

  
       auto lnvtxr = new TLegend(0.55,0.15,0.9,0.3);
       lnvtxr->SetBorderSize(0);
       lnvtxr->SetFillColor(0);
       lnvtxr->AddEntry(Mnvtxr,"Medium","lp");
       lnvtxr->Draw("same");     


       crnvtx->SaveAs("EfSBR_NVTX_TruePV.png");





       TFile *f3 = new TFile("EtasB_updat_SR_TruePVID.root","recreate");
       f3->cd();
       Metas->Write();

       Metab->Write();


       TFile *fb = new TFile(testoutput,"recreate");
       fb->cd();
       Mpts->Write();
       Mptb->Write();
       Mnvtxs->Write();
       Mnvtxb->Write();
       }
