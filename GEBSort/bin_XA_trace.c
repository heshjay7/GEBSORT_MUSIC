#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "Rtypes.h"
#include "TROOT.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TKey.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TTree.h"
#include "gdecomp.h"
#include "GEBSort.h"
#include "GTMerge.h"
#include "bin_XA.h"

/* parameters */
extern DGSEVENT XAEvent[MAXCOINEV];
EvtList* EL = new EvtList(100);

TTree * tree;

extern int XAng;
extern PARS Pars;
int XAtlkup[NCHANNELS];
int XAtid[NCHANNELS];

/*-----------------------------------------------------*/
int sup_XA (){

  printf("------ runing sup_XA(), setup XA, define histogram, mapping \n");

  /* declarations */
  char str[256];
  int i, i1, i2, i7, i8;

  printf("====== creating tree\n");
  tree = new TTree("tree", "tree");

  tree->Branch("NumHits",&EL->NumHits,"NumHits/I");
  tree->Branch("id",EL->id,"id[NumHits]/s");
  tree->Branch("event_timestamp",EL->event_timestamp,"event_timestamp[NumHits]/l"); 
  tree->Branch("pre_rise_energy",EL->pre_rise_energy,"pre_rise_energy[NumHits]/I");
  tree->Branch("post_rise_energy",EL->post_rise_energy,"post_rise_energy[NumHits]/I"); 
  
  tree->Branch("chan_id",EL->chan_id,"chan_id[NumHits]/s");
  tree->Branch("board_id",EL->board_id,"board_id[NumHits]/s");
  tree->Branch("geo_addr",EL->geo_addr,"geo_addr[NumHits]/s");
  tree->Branch("packet_length",EL->packet_length,"packet_length[NumHits]/s");

  tree->Branch("header_type",EL->header_type,"header_type[NumHits]/s");
  tree->Branch("event_type",EL->event_type,"event_type[NumHits]/s");
  tree->Branch("header_length",EL->header_length,"header_length[NumHits]/s");

  tree->Branch("last_disc_timestamp",EL->last_disc_timestamp,"last_disc_timestamp[NumHits]/l");
  tree->Branch("peak_timestamp",EL->peak_timestamp,"peak_timestamp[NumHits]/l");

  tree->Branch("timestamp_match_flag",EL->timestamp_match_flag,"timestamp_match_flag[NumHits]/s");
  tree->Branch("external_disc_flag",EL->external_disc_flag,"external_disc_flag[NumHits]/s");
  tree->Branch("cfd_valid_flag",EL->cfd_valid_flag,"cfd_valid_flag[NumHits]/s");
  tree->Branch("pileup_only_flag",EL->pileup_only_flag,"pileup_only_flag[NumHits]/s");
  tree->Branch("offset_flag",EL->offset_flag,"offset_flag[NumHits]/s");
  tree->Branch("sync_error_flag",EL->sync_error_flag,"sync_error_flag[NumHits]/s");
  tree->Branch("general_error_flag",EL->general_error_flag,"general_error_flag[NumHits]/s");

  tree->Branch("peak_valid_flag",EL->peak_valid_flag,"peak_valid_flag[NumHits]/s");
  tree->Branch("pileup_flag",EL->pileup_flag,"pileup_flag[NumHits]/s");
  tree->Branch("sampled_baseline",EL->sampled_baseline,"sampled_baseline[NumHits]/I");
  tree->Branch("cfd_sample_0",EL->cfd_sample_0,"cfd_sample_0[NumHits]/I");
  tree->Branch("cfd_sample_1",EL->cfd_sample_1,"cfd_sample_1[NumHits]/I");
  tree->Branch("cfd_sample_2",EL->cfd_sample_2,"cfd_sample_2[NumHits]/I");

//  tree->Branch("m1_begin_sample",EL->m1_begin_sample,"m1_begin_sample[NumHits]/s");
//  tree->Branch("m1_end_sample",EL->m1_end_sample,"m1_end_sample[NumHits]/s");
//  tree->Branch("m2_begin_sample",EL->m2_begin_sample,"m2_begin_sample[NumHits]/s");
//  tree->Branch("m2_end_sample",EL->m2_end_sample,"m2_end_sample[NumHits]/s");

  tree->Branch("peak_sample",EL->peak_sample,"peak_sample[NumHits]/s");
  tree->Branch("base_sample",EL->base_sample,"base_sample[NumHits]/s");
  tree->Branch("baseline",EL->baseline,"baseline[NumHits]/I");
  tree->Branch("trace_length",EL->trace_length,"trace_length[NumHits]/s");
  
  //TODO make it controlable from outside
  tree->Branch("trace",EL->trace,"trace[NumHits][1024]/S");
  tree->GetBranch("trace")->SetCompressionSettings(205);
  
  return (0);
};


/*-----------------------------------------------------*/
int exit_XA (){

  /* declarations */
  printf ("\nbegin exit_XA\n");
  
  //tree->Print();
  
  /* done */
  printf ("done exit_XA\n");
    
  return (0);
};


/* ----------------------------------------------------------------- */
int bin_XA (GEB_EVENT * GEB_event) {

  /* declarations */
  int i, j, e, st;
  char str[128];

  /* prototypes */
  int GebTypeStr (int type, char str[]);
  int DGSEvDecompose_v3 (unsigned int *ev, int len, DGSEVENT * XAEvent, int tlkup[], int tid[]);
  if (Pars.CurEvNo <= Pars.NumToPrint) printf ("\nentered bin_XA:\n");

  for (i = 0; i < GEB_event->mult; i++){
    if (GEB_event->ptgd[i]->type == 24) GEB_event->ptgd[i]->type = GEB_TYPE_XA;
  }

  XAng = 0;

  for (i = 0; i < GEB_event->mult; i++){
    if (GEB_event->ptgd[i]->type == GEB_TYPE_XA){
      if (Pars.CurEvNo <= Pars.NumToPrint){
          GebTypeStr (GEB_event->ptgd[i]->type, str);
          printf ("bin_XA (header): %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str, GEB_event->ptgd[i]->timestamp);
      }

      st = DGSEvDecompose_v3 ((unsigned int *) GEB_event->ptinp[i], GEB_event->ptgd[i]->length / sizeof (unsigned int), &XAEvent[XAng], XAtlkup, XAtid);
      if (st != 0) return (0);
      XAng++;
    }
  }

  if (Pars.CurEvNo <= Pars.NumToPrint){
    printf ("XAng=%i\n", XAng);
    for (i = 0; i < XAng; i++){
      printf ("XAEvent[i].tid=%i ", XAEvent[i].tid);
      printf ("tpe=%i ", XAEvent[i].tpe);
      printf ("ts %lli\n", XAEvent[i].event_timestamp);
    };
  };
  
  EL->Reset();
  
  /* Loop tree event */
  EL->NumHits = XAng;
  for (i = 0; i < XAng; i++){

    EL->id[i] = XAEvent[i].id;
    EL->event_timestamp[i] = XAEvent[i].event_timestamp;
    EL->pre_rise_energy[i] = XAEvent[i].sum1;
    EL->post_rise_energy[i] = XAEvent[i].sum2;
    
    EL->chan_id[i] = XAEvent[i].chan_id;
    EL->board_id[i] = XAEvent[i].board_id;
    EL->geo_addr[i] = XAEvent[i].geo_addr;
    EL->packet_length[i] = XAEvent[i].packet_length;
    
    EL->header_type[i] = XAEvent[i].header_type;
    EL->event_type[i] = XAEvent[i].event_type;
    EL->header_length[i] = XAEvent[i].header_length;
  
    EL->last_disc_timestamp[i] = XAEvent[i].last_disc_timestamp;
    EL->peak_timestamp[i] = XAEvent[i].peak_timestamp;
    
    EL->timestamp_match_flag[i] = XAEvent[i].timestamp_match_flag;
    EL->external_disc_flag[i] = XAEvent[i].external_disc_flag;
    EL->cfd_valid_flag[i] = XAEvent[i].cfd_valid_flag;
    EL->pileup_only_flag[i] = XAEvent[i].pileup_only_flag;
    EL->offset_flag[i] = XAEvent[i].offset_flag;
    EL->sync_error_flag[i] = XAEvent[i].sync_error_flag;
    EL->general_error_flag[i] = XAEvent[i].general_error_flag;
  
    EL->peak_valid_flag[i] = XAEvent[i].peak_valid_flag;
    EL->pileup_flag[i] = XAEvent[i].pileup_flag;
    EL->sampled_baseline[i] = XAEvent[i].sampled_baseline;
    EL->cfd_sample_0[i] = XAEvent[i].cfd_sample_0;
    EL->cfd_sample_1[i] = XAEvent[i].cfd_sample_1;
    EL->cfd_sample_2[i] = XAEvent[i].cfd_sample_2;
  
  //EL->m1_begin_sample[i] = XAEvent[i].m1_begin_sample;
  //EL->m1_end_sample[i] = XAEvent[i].m1_end_sample;
  //EL->m2_begin_sample[i] = XAEvent[i].m2_begin_sample;
  //EL->m2_end_sample[i] = XAEvent[i].m2_end_sample;
  
    EL->peak_sample[i] = XAEvent[i].peak_sample;
    EL->base_sample[i] = XAEvent[i].base_sample;
    EL->baseline[i] = XAEvent[i].baseline;
    
    EL->trace_length[i] = XAEvent[i].traceLen;
    
    for(Int_t j=0; (j<EL->trace_length[i]&&j<1024); j++){
      EL->trace[i][j] = XAEvent[i].trace[j];
    }
    
  } 
      
  if (Pars.CurEvNo <= Pars.NumToPrint) printf ("exit bin_XA\n");
  return (0);
}
