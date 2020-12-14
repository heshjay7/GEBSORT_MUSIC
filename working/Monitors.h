//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 23 15:50:10 2020 by ROOT version 6.15/01
// from TTree gen_tree/PSD Tree
// found on file: ../root_data/gen_run004.root
//////////////////////////////////////////////////////////

#ifndef Monitors_h
#define Monitors_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

class Monitors : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runID;
   Int_t           multi_l;
   Int_t           multi_r;
   Float_t         de_l[16];
   Float_t         de_r[16];
   Int_t           seg[16];
   Float_t         stp0;
   Float_t         stp17;
   Float_t         cath;
   Float_t         grid;
   Float_t         tac;
   ULong64_t       de_l_t[16];
   ULong64_t       de_r_t[16];
   ULong64_t       stp0_t;
   ULong64_t       stp17_t;
   ULong64_t       cath_t;
   ULong64_t       grid_t;
   ULong64_t       tac_t;

   // List of branches
   TBranch        *b_runID;   //!
   TBranch        *b_multi_l;   //!
   TBranch        *b_multi_r;   //!
   TBranch        *b_de_l;   //!
   TBranch        *b_de_r;   //!
   TBranch        *b_seg;   //!
   TBranch        *b_stp0;   //!
   TBranch        *b_stp17;   //!
   TBranch        *b_cath;   //!
   TBranch        *b_grid;   //!
   TBranch        *b_tac;   //!
   TBranch        *b_de_l_TimeStamp;   //!
   TBranch        *b_de_r_TimeStamp;   //!
   TBranch        *b_stp0_TimeStamp;   //!
   TBranch        *b_stp17_TimeStamp;   //!
   TBranch        *b_cath_TimeStamp;   //!
   TBranch        *b_grid_TimeStamp;   //!
   TBranch        *b_tac_TimeStamp;   //!

   Monitors(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Monitors() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Monitors,0);
};

#endif

#ifdef Monitors_cxx
void Monitors::Init(TTree *tree)
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("multi_l", &multi_l, &b_multi_l);
   fChain->SetBranchAddress("multi_r", &multi_r, &b_multi_r);
   fChain->SetBranchAddress("de_l", de_l, &b_de_l);
   fChain->SetBranchAddress("de_r", de_r, &b_de_r);
   fChain->SetBranchAddress("seg", seg, &b_seg);
   fChain->SetBranchAddress("stp0", &stp0, &b_stp0);
   fChain->SetBranchAddress("stp17", &stp17, &b_stp17);
   fChain->SetBranchAddress("cath", &cath, &b_cath);
   fChain->SetBranchAddress("grid", &grid, &b_grid);
   fChain->SetBranchAddress("tac", &tac, &b_tac);
   fChain->SetBranchAddress("de_l_t", de_l_t, &b_de_l_TimeStamp);
   fChain->SetBranchAddress("de_r_t", de_r_t, &b_de_r_TimeStamp);
   fChain->SetBranchAddress("stp0_t", &stp0_t, &b_stp0_TimeStamp);
   fChain->SetBranchAddress("stp17_t", &stp17_t, &b_stp17_TimeStamp);
   fChain->SetBranchAddress("cath_t", &cath_t, &b_cath_TimeStamp);
   fChain->SetBranchAddress("grid_t", &grid_t, &b_grid_TimeStamp);
   fChain->SetBranchAddress("tac_t", &tac_t, &b_tac_TimeStamp);
}

Bool_t Monitors::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Monitors_cxx
