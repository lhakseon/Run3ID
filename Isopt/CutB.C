//#define CUTID_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//#include "CutID.h"
using namespace std ;
void CutID::CutBasedID(int reg,double EFFAr[5][3])
{

  TFile *f1;
        if(reg == 1) f1 =  new TFile("Mergedrun3barrel.root","recreate");
      if(reg == 2) f1 =  new TFile("Mergedrun3end.root","recreate");
  TTree *t_S = new TTree("t_S","CUT NEEDS THIS  ");
  TTree *t_B = new TTree("t_B","CUT NEEDS THIS ");
  
  TTree *t_R = new TTree("t_R","CUT NEEDS THIS ");


  TH1F *SieieSw = new TH1F("SieieSw","Sigma ie ie signal pt weighted",1000,0,1);
  TH1F *SieieBw = new TH1F("SieieBw","Sigma ie ie background pt weighted",1000,0,1);

  TH1F *TowSw = new TH1F("TowSw","H o E ie ie signal pt weighted",600,0,6);
  TH1F *TowBw = new TH1F("TowBw","H o E ie ie background pt weighted",600,0,6);

  TH1F *iSOPsw = new TH1F("iSOPsw","Iso gamma  signal pt weighted",300,0,30);
  TH1F *iSOPbw = new TH1F("iSOPbw","Iso gamma  background pt weighted",300,0,30);

  TH1F *iSOCsw = new TH1F("iSOCsw","Iso charge  signal pt weighted",300,0,30);
  TH1F *iSOCbw = new TH1F("iSOCbw","Iso charge  background pt weighted",300,0,30);

  TH1F *iSONsw = new TH1F("iSONsw","Iso neutral  signal pt weighted",300,0,30);
  TH1F *iSONbw = new TH1F("iSONbw","Iso neutral  background pt weighted",300,0,30);

  TH1F *iSOecalsw = new TH1F("iSOecalsw","Iso ecal  signal pt weighted",300,0,30);
  TH1F *iSOecalbw = new TH1F("iSOecalbw","Iso ecal  background pt weighted",300,0,30);

  TH1F *iSOhcalsw = new TH1F("iSOhcalsw","Iso hcal  signal pt weighted",300,0,30);
  TH1F *iSOhcalbw = new TH1F("iSOhcalbw","Iso hcal  background pt weighted",300,0,30);

  
  TH1F *pTs = new TH1F("pTs","Signal Pt ",960,40,1000);
  TH1F *pTb = new TH1F("pTb","Bckground Pt",960,40,1000);

  TH1F *etaS = new TH1F("etaS","Signal Eta   ",20,-5,5);
  TH1F *etaB = new TH1F("etaB","Bckground Eta",20,-5,5);  
  TH1F *etaSw = new TH1F("etaSw","Signal Eta weighted  ",200,-5,5);
  TH1F *etaBw = new TH1F("etaBw","Bckground Eta weighted",200,-5,5);

  TH1F *pTsw = new TH1F("pTsw","Signal Ptpt weighted ",960,40,1000);
  TH1F *pTbw = new TH1F("pTbw","Bckground Pt pt weighted",960,40,1000);
  TH2D *etaPts = new TH2D("etaPts","Eta vs pt",200,-5,5,960,40,1000);
  TH2D *etaPtb = new TH2D("etaPtb","Eta vs pt",200,-5,5,960,40,1000);
  TH2D *etaPtsw = new TH2D("etaPtsw","Eta vs pt weighted",200,-5,5,960,40,1000);
  TH2D *etaPtbw = new TH2D("etaPtbw","Eta vs pt weighted",200,-5,5,960,40,1000);
  
  float genPt,isoecal,isohcal,isoN,isoC,isoP,Peta,Ppt,Pphi,ToE,Sieie,weighT,weighTXS;

  int Pix,Nvtx;

  
  t_S->Branch("Peta",&Peta,"Peta/F");
  t_S->Branch("Pphi",&Pphi,"Pphi/F");
  t_S->Branch("Ppt",&Ppt,"Ppt/F");
  t_S->Branch("isoN",&isoN,"isoN/F");   
  t_S->Branch("isoP",&isoP,"isoP/F");
  t_S->Branch("isoC",&isoC,"isoC/F");
  t_S->Branch("isoecal",&isoecal,"isoecal/F");
  t_S->Branch("isohcal",&isohcal,"isohcal/F");

  

  t_S->Branch("Nvtx",&Nvtx,"Nvtx/I");
  t_S->Branch("genPt",&genPt,"genPt/F");
  t_S->Branch("Sieie",&Sieie,"Sieie/F");
  t_S->Branch("ToE",&ToE,"ToE/F");
  t_S->Branch("Pix",&Pix,"Pix/I");
  t_S->Branch("weighT",&weighT,"weightT/F");
  t_S->Branch("weighTXS",&weighTXS,"weightTXS/F");



  //---------Background Tree---------------

  t_B->Branch("Peta",&Peta,"Peta/F");
  t_B->Branch("Pphi",&Pphi,"Pphi/F");
  t_B->Branch("Ppt",&Ppt,"Ppt/F");
  t_B->Branch("isoN",&isoN,"isoN/F");
  t_B->Branch("isoP",&isoP,"isoP/F");
  t_B->Branch("isoC",&isoC,"isoC/F");
  t_B->Branch("isoecal",&isoecal,"isoecal/F");
  t_B->Branch("isohcal",&isohcal,"isohcal/F");


  

  t_B->Branch("Sieie",&Sieie,"Sieie/F");
  t_B->Branch("ToE",&ToE,"ToE/F");
  t_B->Branch("Pix",&Pix,"Pix/I");
  t_B->Branch("weighT",&weighT,"weightT/F");
  t_B->Branch("Nvtx",&Nvtx,"Nvtx/I");
  t_B->Branch("weighTXS",&weighTXS,"weightTXS/F");



 
  ///////////////////////////////////////////////////////////////////////////////

  //  ----------------  EVENT LOOP BEGINING --------------

  ///////////////////////////////////////////////////////////////////////////////

  // this first loop on the events is for deriving the weights of the trainning
  // so we are only interested on the pt distribution of back/ signal events
  

  if (fChain == 0) return;
  Long64_t nentries =fChain->GetEntriesFast();
  cout<<"nentries "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
        if(gedPhPt< 40) continue;
	if(fabs(MCVtz - Vtz) >= 0.1)continue;
        if(fabs( sqrt( pow(MCVtx-Vtx,2) +  pow(MCVty-Vty,2) ) ) >= 0.2 ) continue;
        if(reg == 1 && fabs(gedPhEta) > 1.4442 ) continue;
        if(reg == 2 && ( fabs(gedPhEta) < 1.566 || fabs(gedPhEta) > 2.5  )) continue;
        if(gedPhisPrompt){
      pTs->Fill(gedPhPt,gedPhweightXS);
      etaS->Fill(gedPhEta,gedPhweightXS);
      etaPts->Fill(gedPhEta,gedPhPt,gedPhweightXS);
      
    }
    
    if(gedPhisPrompt == 0){
      pTb->Fill(gedPhPt,gedPhweightXS);
      etaB->Fill(gedPhEta,gedPhweightXS);
      etaPtb->Fill(gedPhEta,gedPhPt,gedPhweightXS);
    }
}
  cout<<"end of first loop on events"<<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb =fChain->GetEntry(jentry);   nbytes += nb;
    if(gedPhPt < 40)continue;
    if(fabs(MCVtz - Vtz) >= 0.1)continue;
    if(fabs( sqrt( pow(MCVtx-Vtx,2) +  pow(MCVty-Vty,2) ) ) >= 0.2 ) continue;
    //---------------------------------------------------------------------
    //      CALC OF THE RHO SUB - ISOLATION 
    //---------------------------------------------------------------------
    double EAhcal = 0; 
    double EAecal = 0;
    double EAch = 0;	
    
    //Setting up the effective areas
    int biin;
    double eta = gedPhEta;
    // double eta = phoEta;

    //if(fabs(eta) < 1.0 )biin = 0;
    if(fabs(eta) > 0.0 && fabs(eta) < 0.5)biin = 0;
    if(fabs(eta) > 0.5 && fabs(eta) < 1.0)biin = 1;
    if(fabs(eta) > 1.0 && fabs(eta) < 1.4442)biin = 2;
    if(fabs(eta) > 1.566 && fabs(eta) < 2.0)biin = 3;
    if(fabs(eta) > 2.0 && fabs(eta) < 2.5)biin = 4;
    
    
    
    EAhcal  = EFFAr[biin][2];
    EAecal = EFFAr[biin][1];
    EAch  = EFFAr[biin][0];
        


    float ChgI = gedChgIso;
    float ecalI = gedecalIso;	
    float hcalI = gedhcalIso;	
        
    ChgI = max( ChgI - Rh*EAch,0.0);
    ecalI = max ( ecalI - Rh*EAecal,0.0);
    hcalI = max ( hcalI - Rh*EAhcal,0.0);


    //EOF THE ISOLATION CALCULATION

    if(reg == 1 && fabs(gedPhEta) > 1.4442 ) continue;
    if(reg == 2 && ( fabs(gedPhEta) < 1.566 || fabs(gedPhEta) > 2.5  )) continue;

    if(gedPhisPrompt){
      Peta=gedPhEta;
      Pphi=gedPhPhi;
      Ppt=gedPhPt;
      Pix = gedPhPixSeed;
      isoC = ChgI;
      isoecal = ecalI; 
      isohcal = hcalI; 
      Sieie =gedPhSieie ;
      ToE = gedPhTower;
      genPt = gedGenPt; 
      Nvtx = NVtx; 

      int binx = etaPts->FindBin(gedPhEta,gedPhPt);
      weighT = gedPhweightXS*( ( etaPts->GetBinContent(binx) == 0  ) ? 0.0  : 1./etaPts->GetBinContent(binx));
      weighTXS = gedPhweightXS;
      SieieSw->Fill(Sieie,weighT);
      TowSw->Fill(ToE,weighT);
      iSOPsw->Fill(isoP,weighT);
      iSONsw->Fill(isoN,weighT);
      iSOCsw->Fill(isoC,weighT);
      iSOecalsw->Fill(isoecal,weighT);
      iSOhcalsw->Fill(isohcal,weighT);
      pTsw->Fill(Ppt,weighT);
      etaSw->Fill(Peta,weighT);
      etaPtsw->Fill(Peta,Ppt,weighT);
      t_S->Fill();     
    
    
}
    if(gedPhisPrompt == 0){
      
      Peta = gedPhEta;
      Pphi = gedPhPhi;
      Ppt  = gedPhPt;
      Pix  = gedPhPixSeed;
      
      isoC = ChgI;
      isoecal = ecalI; 
      isohcal = hcalI; 
      
      Sieie = gedPhSieie ;
      ToE   = gedPhTower;
      Nvtx  = NVtx; 


      int binx = etaPtb->FindBin(gedPhEta,gedPhPt);
      weighT = gedPhweightXS*(( etaPtb->GetBinContent(binx) == 0  ) ? 0.0  : 1./etaPtb->GetBinContent(binx));
      weighTXS = gedPhweightXS;
      
      SieieBw->Fill(Sieie,weighT);
      TowBw->Fill(ToE,weighT);
      iSOPbw->Fill(isoP,weighT);
      iSONbw->Fill(isoN,weighT);
      iSOCbw->Fill(isoC,weighT);
      iSOecalbw->Fill(isoecal,weighT);
      iSOhcalbw->Fill(isohcal,weighT);
      pTbw->Fill(Ppt,weighT);
      etaBw->Fill(Peta,weighT);
      etaPtbw->Fill(Peta,Ppt,weighT);
      t_B->Fill();
					       
    }
    



  }//eof the event loop


  t_S->Write();
  t_B->Write();
  t_R->Write();


  SieieSw->Write();
  TowSw->Write();
  iSOPsw->Write();
  iSONsw->Write();
  iSOCsw->Write();  

  SieieBw->Write();
  TowBw->Write();
  iSOPbw->Write();
  iSONbw->Write();
  iSOCbw->Write();
  iSOecalsw->Write();
  iSOhcalsw->Write();
  iSOecalbw->Write();
  iSOhcalbw->Write();
  
  pTsw->Write();
  pTbw->Write();

  etaSw->Write();
  etaBw->Write();

  etaPts->Write();
  etaPtb->Write();

  etaPtsw->Write();
  etaPtbw->Write();

  t_S->Print();
  
  f1->Close();
}// loop function end







