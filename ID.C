#define ID_cxx
#include "ID.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <TMath.h>
#include <vector>
#include <TVector3.h>
#include <TProfile.h>
#include "TStopwatch.h"
#include <algorithm>
#include "TROOT.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;


int main(int argc, const char *argv[]){
  gROOT->ProcessLine("#include <vector>");

  string outfile = argv[2];
  long int nTotEvt = atof(argv[3]);
  //std::cout<<"r a valid value for maxEvents (parameter 3)."<<std          ::endl;                                                                                                                                                          
  long int nPrintEvt = atof(argv[4]);
  double evtWt = atof(argv[5]);
  ID m(argv[1],argv[2]);
  m.Loop(outfile,nTotEvt,nPrintEvt,evtWt);
  //return 0;                                                                                                                                                                                                                                

}

  void ID::Loop(string outputfile, int ntotEvent, int nPrintEvent, double evtWeight)
  {
    TStopwatch sw;
    sw.Start();

    cout<<"Output file : " << outputfile << endl;
    char *outfilename = const_cast<char*>(outputfile.c_str());
    fileName = new TFile(outfilename,"RECREATE");
    tree = new TTree("tree","EventsTree");
    TH1D *MatchPhos = new TH1D("MatchPhos","# of matched photons in event",5,0,5);

    
    float  Vtx,Vty,Vtz,MCVtx,MCVty,MCVtz,gedPhR9,gedPhEta,gedPhPt,gedPhPhi,gedPhoIso,gedChgIso,gedNeuIso,gedhcalIso,gedecalIso,gedPhTower,gedPhSieie,Rh,gedGenPt,gedPhweightXS;//,gedPhIDMVA;
    bool gedPhPixSeed;
     int NVtx,gedPhisPrompt,hascorrectPrimV,gedPhEleVeto;

    fileName->cd();

    //Float Branches                                                                                                                                                    

 tree->Branch("Vtx",&Vtx,"Vtx/F");
 tree->Branch("Vty",&Vty,"Vty/F");
 tree->Branch("Vtz",&Vtz,"Vtz/F");
 tree->Branch("MCVtx",&MCVtx,"MCVtx/F");
 tree->Branch("MCVty",&MCVty,"MCVty/F");
 tree->Branch("MCVtz",&MCVtz,"MCVtz/F");
 tree->Branch("gedPhR9",&gedPhR9,"gedPhR9/F");
 tree->Branch("Rh",&Rh,"Rh/F");
 tree->Branch("gedPhEta",&gedPhEta,"gedPhEta/F");
 tree->Branch("gedPhPhi",&gedPhPhi,"gedPhPhi/F");
 tree->Branch("gedPhPt",&gedPhPt,"gedPhPt/F");
 tree->Branch("gedPhSieie",&gedPhSieie,"gedPhSieie/F");
 tree->Branch("gedPhTower",&gedPhTower,"gedPhTower/F");
 tree->Branch("gedPhoIso",&gedPhoIso,"gedPhoIso/F");
 tree->Branch("gedChgIso",&gedChgIso,"gedChgIso/F");
 tree->Branch("gedNeuIso",&gedNeuIso,"gedNeuIso/F");


 tree->Branch("gedhcalIso",&gedhcalIso,"gedhcalIso/F");
 tree->Branch("gedecalIso",&gedecalIso,"gedecalIso/F");


 tree->Branch("gedGenPt",&gedGenPt,"gedGenPt/F");
 tree->Branch("gedPhweightXS",&gedPhweightXS,"gedPhweightXS/F");
 
 tree->Branch("NVtx",&NVtx,"NVtx/I");
 tree->Branch("gedPhPixSeed",&gedPhPixSeed,"gedPhPixSeed/I");
 tree->Branch("gedPhisPrompt",&gedPhisPrompt,"gedPhisPrompt/I");
 tree->Branch("hascorrectPrimV",&hascorrectPrimV,"hascorrectPrimV/I");
 tree->Branch("gedPhEleVeto",&gedPhEleVeto,"gedPhEleVeto/I");
 if (fChain == 0) return;
 Long64_t nentries = fChain->GetEntries();

 Long64_t nbytes = 0, nb = 0;
 double weightt1 = evtWeight;
 cout<<"weight of  this sample "<<weightt1<<endl;
 cout<<"total entries:"<<nentries<<endl;

 cout<<"nPrintevents"<<nPrintEvent<<" "<<ntotEvent<<endl;
 Long64_t nentriesToCheck = nentries;
 if (ntotEvent != -1LL && nentries > ntotEvent)nentriesToCheck = ntotEvent;

 for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++) {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;

   if (jentry % nPrintEvent == 0) std::cout << "  " << jentry  << "  Events Processed... " << std::endl;
   int mPHOS = 0;
   int mto1st = 0;

   Vtx= PV_x;
   Vty=PV_y;
   Vtz=PV_z;

   MCVtx=GenVtx_x;
   MCVty=GenVtx_y;
   MCVtz=GenVtx_z;
   if (nPhoton<0) continue;
   

   for(int ipho = 0; ipho < nPhoton; ipho++){
     double pEta = Photon_eta[ipho];
     double pPhi = Photon_phi[ipho];
     Rh = Rho_fixedGridRhoAll;

     NVtx = PV_npvs;
     gedPhPt = Photon_pt[ipho];
     gedPhEta = pEta;
     gedPhPhi = pPhi;
     gedPhSieie = Photon_sieie[ipho];
     gedPhTower = Photon_hoe[ipho];
     gedPhoIso = Photon_pfPhoIso03[ipho];
     gedChgIso = Photon_pfChargedIsoWorstVtx[ipho];
     gedNeuIso = Photon_pfChargedIsoWorstVtx[ipho];
     gedhcalIso = Photon_hcalPFClusterIso[ipho];
     gedecalIso = Photon_ecalPFClusterIso[ipho];


     gedPhPixSeed  = Photon_pixelSeed[ipho];
     gedPhweightXS = weightt1;

     int pass = 0;
     int HLTPhoIsPrescaled = 0;                  

     double genPt = -1;
     for(int imc = 0; imc < nGenPart; imc++){
         if(GenPart_pt[imc] < 15 ) continue;
         if(GenPart_status[imc] != 1)continue;
         if(GenPart_pdgId[imc] != 22)continue;
         pass++;
         double meta = GenPart_eta[imc];
         double mphi = GenPart_phi[imc];
         TVector3 mcphoton;
         TVector3 recoPHOTOn;
         mcphoton.SetPtEtaPhi(1.0,meta,mphi);
         recoPHOTOn.SetPtEtaPhi(1.0,pEta,pPhi);
         double DR = mcphoton.DrEtaPhi(recoPHOTOn);
         double dp = fabs(GenPart_pt[imc] - Photon_pt[ipho] )/GenPart_pt[imc];
         if(DR < 0.2 && dp < 0.2  ){
             mPHOS++;
             if(pass == 1 ){
                 genPt = GenPart_pt[imc];
                 HLTPhoIsPrescaled= 1;
             }
         }
     }//EOF MC Particles loop                                                                                                                                           

     gedPhisPrompt = HLTPhoIsPrescaled;
     gedGenPt = genPt;
     tree->Fill();
   }//EOF Photon Loop                                                                                                                                                   

   MatchPhos->Fill(mPHOS);

 }//EOF EVENT LOOP                                                                                                                                                      

 MatchPhos->Write();
                                                                       

 cout<<outputfile<< " created ...."<<endl;

 sw.Stop();
 std::cout << "RealTime : " << sw.RealTime() / 60.0 << " minutes" << std::endl;
 std::cout << "CPUTime  : " << sw.CpuTime()  / 60.0 << " minutes" << std::endl;

  }
  void ID::histograms(const char* file2)
  {
    fileName = new TFile(file2, "RECREATE");
    fileName->cd();
  }



