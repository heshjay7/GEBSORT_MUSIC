//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May  7 16:10:02 2015 by ROOT version 5.32/00
// from TTree tree/tree
// found on file: run5113.root
//////////////////////////////////////////////////////////

#ifndef GeneralSortTraceProof_h
#define GeneralSortTraceProof_h

#include <TROOT.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TSelector.h>
#include <TF1.h>
#include <TGraph.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TProofOutputFile.h>

class GeneralSortTraceProof : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           RunNum;
   Int_t           NumHits;
   Short_t         id[200];   //[NumHits]
   UShort_t        tpe[200];   //[NumHits]
   UShort_t        tid[200];   //[NumHits]
   UShort_t        chan_id[200];   //[NumHits]
   UShort_t        board_id[200];   //[NumHits]
   UShort_t        geo_addr[200];   //[NumHits]
   UShort_t        packet_length[200];   //[NumHits]
   UShort_t        header_type[200];   //[NumHits]
   UShort_t        event_type[200];   //[NumHits]
   UShort_t        header_length[200];   //[NumHits]
   ULong64_t       event_timestamp[200];   //[NumHits]
   ULong64_t       last_disc_timestamp[200];   //[NumHits]
   ULong64_t       peak_timestamp[200];   //[NumHits]
   UShort_t        timestamp_match_flag[200];   //[NumHits]
   UShort_t        external_disc_flag[200];   //[NumHits]
   UShort_t        cfd_valid_flag[200];   //[NumHits]
   UShort_t        pileup_only_flag[200];   //[NumHits]
   UShort_t        offset_flag[200];   //[NumHits]
   UShort_t        sync_error_flag[200];   //[NumHits]
   UShort_t        general_error_flag[200];   //[NumHits]
   UShort_t        peak_valid_flag[200];   //[NumHits]
   UShort_t        pileup_flag[200];   //[NumHits]
   Int_t           sampled_baseline[200];   //[NumHits]
   Int_t           cfd_sample_0[200];   //[NumHits]
   Int_t           cfd_sample_1[200];   //[NumHits]
   Int_t           cfd_sample_2[200];   //[NumHits]
   Int_t           pre_rise_energy[200];   //[NumHits]
   Int_t           post_rise_energy[200];   //[NumHits]
   UShort_t        m1_begin_sample[200];   //[NumHits]
   UShort_t        m1_end_sample[200];   //[NumHits]
   UShort_t        m2_begin_sample[200];   //[NumHits]
   UShort_t        m2_end_sample[200];   //[NumHits]
   UShort_t        peak_sample[200];   //[NumHits]
   UShort_t        base_sample[200];   //[NumHits]
   Int_t           baseline[200];   //[NumHits]
   UShort_t        trace_length[200];   //[NumHits]
   Short_t         trace[200][1024];   //[NumHits]

   // List of branches
   TBranch        *b_RunNum;
   TBranch        *b_NumHits;   //!
   TBranch        *b_id;   //!
   TBranch        *b_tpe;   //!
   TBranch        *b_tid;   //!
   TBranch        *b_chan_id;   //!
   TBranch        *b_board_id;   //!
   TBranch        *b_geo_addr;   //!
   TBranch        *b_packet_length;   //!
   TBranch        *b_header_type;   //!
   TBranch        *b_event_type;   //!
   TBranch        *b_header_length;   //!
   TBranch        *b_event_timestamp;   //!
   TBranch        *b_last_disc_timestamp;   //!
   TBranch        *b_peak_timestamp;   //!
   TBranch        *b_timestamp_match_flag;   //!
   TBranch        *b_external_disc_flag;   //!
   TBranch        *b_cfd_valid_flag;   //!
   TBranch        *b_pileup_only_flag;   //!
   TBranch        *b_offset_flag;   //!
   TBranch        *b_sync_error_flag;   //!
   TBranch        *b_general_error_flag;   //!
   TBranch        *b_peak_valid_flag;   //!
   TBranch        *b_pileup_flag;   //!
   TBranch        *b_sampled_baseline;   //!
   TBranch        *b_cfd_sample_0;   //!
   TBranch        *b_cfd_sample_1;   //!
   TBranch        *b_cfd_sample_2;   //!
   TBranch        *b_pre_rise_energy;   //!
   TBranch        *b_post_rise_energy;   //!
   TBranch        *b_m1_begin_sample;   //!
   TBranch        *b_m1_end_sample;   //!
   TBranch        *b_m2_begin_sample;   //!
   TBranch        *b_m2_end_sample;   //!
   TBranch        *b_peak_sample;   //!
   TBranch        *b_base_sample;   //!
   TBranch        *b_baseline;   //!
   TBranch        *b_trace_length;   //!
   TBranch        *b_trace;   //!
   
   //======================= new tree, new file
   TString saveFileName;
   TFile *saveFile; //!
   TProofOutputFile * proofFile; //!
   TTree *newTree; //!

   //trace
   TClonesArray * arr ;//!
   TGraph * gTrace; //!
   TF1 * gFit; //!

   float te[30];    // energy from trace
   float te_r[30];  // rising time from frace
   float te_t[30];  // time
   
   float ttac[6];
   float ttac_t[6];
   float ttac_r[6];
   
   float trdt[8];
   float trdt_t[8];
   float trdt_r[8];
   
   float tezero[8];
   float tezero_t[8];
   float tezero_r[8];

   int runIDLast;
   TString fileNameTemp;
   
   //PSD struct
   typedef struct {
      Int_t   eventID;
      Int_t   runID;
      Float_t Energy[100];
      Float_t XF[100];
      Float_t XN[100];
      Float_t Ring[100];
      Float_t RDT[100];
      Float_t TAC[100];
      Float_t ELUM[32];
      Float_t EZERO[10];

      ULong64_t EnergyTimestamp[100];
      ULong64_t XFTimestamp[100];
      ULong64_t XNTimestamp[100];
      ULong64_t RingTimestamp[100];
      ULong64_t RDTTimestamp[100];
      ULong64_t TACTimestamp[100];
      ULong64_t ELUMTimestamp[32];
      ULong64_t EZEROTimestamp[10];
      
      Float_t x[100];
      
   } PSD;

   PSD psd; 
   
   //need to put NULL on pointer
 GeneralSortTraceProof(TTree * /*tree*/ = 0) : fChain(0), saveFile(0), proofFile(0), newTree(0), arr(0), gTrace(0), gFit(0) { }
   virtual ~GeneralSortTraceProof() { }
   
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

   ClassDef(GeneralSortTraceProof,0);

    
};

#endif

#ifdef GeneralSortTraceProof_cxx
void GeneralSortTraceProof::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // this will be called when a new file is loaded in a new tree? possibly
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("NumHits", &NumHits, &b_NumHits);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("tpe", tpe, &b_tpe);
   fChain->SetBranchAddress("tid", tid, &b_tid);
   fChain->SetBranchAddress("chan_id", chan_id, &b_chan_id);
   fChain->SetBranchAddress("board_id", board_id, &b_board_id);
   fChain->SetBranchAddress("geo_addr", geo_addr, &b_geo_addr);
   fChain->SetBranchAddress("packet_length", packet_length, &b_packet_length);
   fChain->SetBranchAddress("header_type", header_type, &b_header_type);
   fChain->SetBranchAddress("event_type", event_type, &b_event_type);
   fChain->SetBranchAddress("header_length", header_length, &b_header_length);
   fChain->SetBranchAddress("event_timestamp", event_timestamp, &b_event_timestamp);
   fChain->SetBranchAddress("last_disc_timestamp", last_disc_timestamp, &b_last_disc_timestamp);
   fChain->SetBranchAddress("peak_timestamp", peak_timestamp, &b_peak_timestamp);
   fChain->SetBranchAddress("timestamp_match_flag", timestamp_match_flag, &b_timestamp_match_flag);
   fChain->SetBranchAddress("external_disc_flag", external_disc_flag, &b_external_disc_flag);
   fChain->SetBranchAddress("cfd_valid_flag", cfd_valid_flag, &b_cfd_valid_flag);
   fChain->SetBranchAddress("pileup_only_flag", pileup_only_flag, &b_pileup_only_flag);
   fChain->SetBranchAddress("offset_flag", offset_flag, &b_offset_flag);
   fChain->SetBranchAddress("sync_error_flag", sync_error_flag, &b_sync_error_flag);
   fChain->SetBranchAddress("general_error_flag", general_error_flag, &b_general_error_flag);
   fChain->SetBranchAddress("peak_valid_flag", peak_valid_flag, &b_peak_valid_flag);
   fChain->SetBranchAddress("pileup_flag", pileup_flag, &b_pileup_flag);
   fChain->SetBranchAddress("sampled_baseline", sampled_baseline, &b_sampled_baseline);
   fChain->SetBranchAddress("cfd_sample_0", cfd_sample_0, &b_cfd_sample_0);
   fChain->SetBranchAddress("cfd_sample_1", cfd_sample_1, &b_cfd_sample_1);
   fChain->SetBranchAddress("cfd_sample_2", cfd_sample_2, &b_cfd_sample_2);
   fChain->SetBranchAddress("pre_rise_energy", pre_rise_energy, &b_pre_rise_energy);
   fChain->SetBranchAddress("post_rise_energy", post_rise_energy, &b_post_rise_energy);
   fChain->SetBranchAddress("m1_begin_sample", m1_begin_sample, &b_m1_begin_sample);
   fChain->SetBranchAddress("m1_end_sample", m1_end_sample, &b_m1_end_sample);
   fChain->SetBranchAddress("m2_begin_sample", m2_begin_sample, &b_m2_begin_sample);
   fChain->SetBranchAddress("m2_end_sample", m2_end_sample, &b_m2_end_sample);
   fChain->SetBranchAddress("peak_sample", peak_sample, &b_peak_sample);
   fChain->SetBranchAddress("base_sample", base_sample, &b_base_sample);
   fChain->SetBranchAddress("baseline", baseline, &b_baseline);
   fChain->SetBranchAddress("trace_length", trace_length, &b_trace_length);
   fChain->SetBranchAddress("trace", trace, &b_trace);

   int NumEntries = fChain->GetEntries();
   printf( "========== total Entry : %d\n", NumEntries);

   fileNameTemp = fChain->GetDirectory()->GetName();
   //remove any folder path to get the name;
   int found;
   do{
     found = fileNameTemp.First("/");
     fileNameTemp.Remove(0,found+1);
   }while( found >= 0 );
   
   found = fileNameTemp.First(".");
   fileNameTemp.Remove(found);
   //the fileNameTemp should be something like "xxx_run4563" now
   
   while( !fileNameTemp.IsDigit() ){
      fileNameTemp.Remove(0,1);
   }
   
   runIDLast = fileNameTemp.Atoi();
   
}

Bool_t GeneralSortTraceProof::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
  
   return kTRUE;
}

#endif // #ifdef GeneralSortTraceProof_cxx
