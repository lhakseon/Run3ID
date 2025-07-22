//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 22 04:00:35 2022 by ROOT version 6.14/09
// from TTree tree/EventsTree
// found on file: flat_pt.root
//////////////////////////////////////////////////////////

#ifndef CutID_h
#define CutID_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class CutID {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         Vtx;
   Float_t         Vty;
   Float_t         Vtz;
   Float_t         MCVtx;
   Float_t         MCVty;
   Float_t         MCVtz;
   Float_t         gedPhR9;
   Float_t         Rh;
   Float_t         gedPhEta;
   Float_t         gedPhPhi;
   Float_t         gedPhPt;
   Float_t         gedPhSieie;
   Float_t         gedPhTower;
   Float_t         gedPhoIso;
   Float_t         gedecalIso;
   Float_t         gedhcalIso;
   Float_t         gedChgIso;
   Float_t         gedNeuIso;
   Float_t         gedPhSCEta;
   Float_t         gedGenPt;
   Float_t         gedPhweightXS;
   Float_t         PhCITKChIso;
   Float_t         PhCITKNeuIso;
   Float_t         PhCITKPhoIso;
   Float_t         PhPUPPIChIso;
   Float_t         PhPUPPINeuIso;
   Float_t         PhPUPPIPhoIso;
   Float_t         gedPhIDMVA;
   Int_t           NVtx;
   Int_t           gedPhPixSeed;
   Int_t           gedPhisPrompt;
   Int_t           hascorrectPrimV;
   Int_t           gedPhEleVeto;

   // List of branches
   TBranch        *b_Vtx;   //!
   TBranch        *b_Vty;   //!
   TBranch        *b_Vtz;   //!
   TBranch        *b_MCVtx;   //!
   TBranch        *b_MCVty;   //!
   TBranch        *b_MCVtz;   //!
   TBranch        *b_gedPhR9;   //!
   TBranch        *b_Rh;   //!
   TBranch        *b_gedPhEta;   //!
   TBranch        *b_gedPhPhi;   //!
   TBranch        *b_gedPhPt;   //!
   TBranch        *b_gedPhSieie;   //!
   TBranch        *b_gedPhTower;   //!
   TBranch        *b_gedPhoIso;   //!
   TBranch        *b_gedhcalIso;   //!
   TBranch        *b_gedecalIso;   //!
   TBranch        *b_gedChgIso;   //!
   TBranch        *b_gedNeuIso;   //!
   TBranch        *b_gedPhSCEta;   //!
   TBranch        *b_gedGenPt;   //!
   TBranch        *b_gedPhweightXS;   //!
   TBranch        *b_PhCITKChIso;   //!
   TBranch        *b_PhCITKNeuIso;   //!
   TBranch        *b_PhCITKPhoIso;   //!
   TBranch        *b_PhPUPPIChIso;   //!
   TBranch        *b_PhPUPPINeuIso;   //!
   TBranch        *b_PhPUPPIPhoIso;   //!
   TBranch        *b_gedPhIDMVA;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_gedPhPixSeed;   //!
   TBranch        *b_gedPhisPrompt;   //!
   TBranch        *b_hascorrectPrimV;   //!
   TBranch        *b_gedPhEleVeto;   //!

   CutID(TTree *tree=0);
   virtual ~CutID();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     CutBasedID(int reg,double EFFAr[5][3]);
};

#endif

#ifdef CutID_cxx
CutID::CutID(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Flat_all.root");
      if (!f || !f->IsOpen()) {
       f = new TFile("../Flat_all.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

CutID::~CutID()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CutID::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CutID::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CutID::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Vtx", &Vtx, &b_Vtx);
   fChain->SetBranchAddress("Vty", &Vty, &b_Vty);
   fChain->SetBranchAddress("Vtz", &Vtz, &b_Vtz);
   fChain->SetBranchAddress("MCVtx", &MCVtx, &b_MCVtx);
   fChain->SetBranchAddress("MCVty", &MCVty, &b_MCVty);
   fChain->SetBranchAddress("MCVtz", &MCVtz, &b_MCVtz);
   fChain->SetBranchAddress("gedPhR9", &gedPhR9, &b_gedPhR9);
   fChain->SetBranchAddress("Rh", &Rh, &b_Rh);
   fChain->SetBranchAddress("gedPhEta", &gedPhEta, &b_gedPhEta);
   fChain->SetBranchAddress("gedPhPhi", &gedPhPhi, &b_gedPhPhi);
   fChain->SetBranchAddress("gedPhPt", &gedPhPt, &b_gedPhPt);
   fChain->SetBranchAddress("gedPhSieie", &gedPhSieie, &b_gedPhSieie);
   fChain->SetBranchAddress("gedPhTower", &gedPhTower, &b_gedPhTower);
   fChain->SetBranchAddress("gedPhoIso", &gedPhoIso, &b_gedPhoIso);
   fChain->SetBranchAddress("gedhcalIso", &gedhcalIso, &b_gedhcalIso);
   fChain->SetBranchAddress("gedecalIso", &gedecalIso, &b_gedecalIso);
   fChain->SetBranchAddress("gedChgIso", &gedChgIso, &b_gedChgIso);
   fChain->SetBranchAddress("gedNeuIso", &gedNeuIso, &b_gedNeuIso);
   fChain->SetBranchAddress("gedPhSCEta", &gedPhSCEta, &b_gedPhSCEta);
   fChain->SetBranchAddress("gedGenPt", &gedGenPt, &b_gedGenPt);
   fChain->SetBranchAddress("gedPhweightXS", &gedPhweightXS, &b_gedPhweightXS);
   fChain->SetBranchAddress("PhCITKChIso", &PhCITKChIso, &b_PhCITKChIso);
   fChain->SetBranchAddress("PhCITKNeuIso", &PhCITKNeuIso, &b_PhCITKNeuIso);
   fChain->SetBranchAddress("PhCITKPhoIso", &PhCITKPhoIso, &b_PhCITKPhoIso);
   fChain->SetBranchAddress("PhPUPPIChIso", &PhPUPPIChIso, &b_PhPUPPIChIso);
   fChain->SetBranchAddress("PhPUPPINeuIso", &PhPUPPINeuIso, &b_PhPUPPINeuIso);
   fChain->SetBranchAddress("PhPUPPIPhoIso", &PhPUPPIPhoIso, &b_PhPUPPIPhoIso);
   fChain->SetBranchAddress("gedPhIDMVA", &gedPhIDMVA, &b_gedPhIDMVA);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("gedPhPixSeed", &gedPhPixSeed, &b_gedPhPixSeed);
   fChain->SetBranchAddress("gedPhisPrompt", &gedPhisPrompt, &b_gedPhisPrompt);
   fChain->SetBranchAddress("hascorrectPrimV", &hascorrectPrimV, &b_hascorrectPrimV);
   fChain->SetBranchAddress("gedPhEleVeto", &gedPhEleVeto, &b_gedPhEleVeto);
   Notify();
}

Bool_t CutID::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CutID::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CutID::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CutID_cxx
