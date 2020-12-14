#define GeneralSort_cxx

#include "GeneralSort.h"
#include <TH2.h>
#include <TMath.h>
#include <TMacro.h>
#include <TStyle.h>
#include <TBenchmark.h>

#define NUMPRINT 20 //>0
#define MAXNUMHITS 200 //Highest multiplicity
#define M 100.0 //M value for energy filter from digi setting, number of channel

ULong64_t MaxProcessedEntries=100000000;
ULong64_t NumEntries = 0;
int EffEntries;
ULong64_t ProcessedEntries = 0;

#include "../working/GeneralSortMapping.h"

//clock
TBenchmark gClock;
Bool_t shown = 0;

TFile *oFile;
TTree *gen_tree;

//PSD struct
typedef struct {
  Int_t runID;
  Float_t         de_l[16];
  Float_t         de_r[16];
  Int_t           seg[16];
  Float_t         cath;
  Float_t         grid;
  Float_t         stp0;
  Float_t         stp17;
  Float_t         tac;

  ULong64_t         de_l_TimeStamp[16];
  ULong64_t         de_r_TimeStamp[16];
  // ULong64_t         seg_TimeStamp[16];
  ULong64_t         cath_TimeStamp;
  ULong64_t         grid_TimeStamp;
  ULong64_t         stp0_TimeStamp;
  ULong64_t         stp17_TimeStamp;
  ULong64_t         tac_TimeStamp;
  
  Int_t       multi_l;
  Int_t       multi_r;


} MUSIC;

MUSIC music; 

TString saveFileName = "gen.root";

int runIDPresent = 0;
TString fileNum = "";

void GeneralSort::Begin(TTree * tree)
{
   
  TString option = GetOption();
  NumEntries = tree->GetEntries();
  EffEntries = TMath::Min(MaxProcessedEntries, NumEntries);

  saveFileName = tree->GetDirectory()->GetName();
  int findslat = saveFileName.Last('/');
  saveFileName.Remove(0, findslat+1);
  saveFileName = "../root_data/gen_" + saveFileName;

  printf("=============================================================\n");
  printf("=====================  GeneralSort.C  =======================\n");
  printf("=============================================================\n");
  printf("                    file : %s \n", tree->GetDirectory()->GetName());
  printf("          Number of Event: %llu\n", NumEntries);
  printf("Effective Number of Event: %d <= %llu\n", EffEntries, MaxProcessedEntries);  

  oFile = new TFile(saveFileName,"RECREATE");

  gen_tree = new TTree("gen_tree","MUSIC Tree");

  gen_tree->Branch("runID", &music.runID,"runID/I");
  gen_tree->Branch("multi_l", &music.multi_l,"multi_l/I");
  gen_tree->Branch("multi_r", &music.multi_r,"multi_r/I");

  gen_tree->Branch("de_l",  &music.de_l,     "de_l[16]/F");
  gen_tree->Branch("de_r",  &music.de_r,     "de_r[16]/F");
  gen_tree->Branch("seg",   &music.seg,      "seg[16]/I");
  gen_tree->Branch("stp0",  &music.stp0,  "stp0/F");
  gen_tree->Branch("stp17", &music.stp17, "stp17/F");
  gen_tree->Branch("cath",  &music.cath, "cath/F");
  gen_tree->Branch("grid",  &music.grid, "grid/F");
  gen_tree->Branch("tac",  &music.tac, "tac/F");

  gen_tree->Branch("de_l_t",  &music.de_l_TimeStamp,     "de_l_TimeStamp[16]/l");
  gen_tree->Branch("de_r_t",  &music.de_r_TimeStamp,     "de_r_TimeStamp[16]/l");
  gen_tree->Branch("stp0_t",  &music.stp0_TimeStamp,  "stp0_TimeStamp/l");
  gen_tree->Branch("stp17_t", &music.stp17_TimeStamp, "stp17_TimeStamp/l");
  gen_tree->Branch("cath_t",  &music.cath_TimeStamp, "cath_TimeStamp/l");
  gen_tree->Branch("grid_t",  &music.grid_TimeStamp, "grid_TimeStamp/l");
  gen_tree->Branch("tac_t",  &music.tac_TimeStamp, "tac_TimeStamp/l");
  
  printf("======= ID-MAP: \n");
  printf("%11s|", ""); 
  for(int i = 0 ; i < 10; i++ ) printf("%7d|", i);
  printf("\n");
  for(int i = 0 ; i < 96; i++ ) printf("-");
  for(int i = 0 ; i < 40; i ++){ //*************************
    if( (i) % 10 == 0 ) {
       printf("\n");
       if(((i+1)/10)/4+1 < 5) printf("%11s|", Form("VME%d-Dig%d", ((i+1)/10)/4+1, ((i+1)/10)%4+1)); 
    }
    if( 115 >= idDetMap[i] && idDetMap[i] >= 100){
      printf("\033[36m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); // Left,
    }else if( 215 >= idDetMap[i] && idDetMap[i] >= 200){
      printf("\033[91m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); // Right, 
    }else if( 304 >= idDetMap[i] && idDetMap[i] >= 300){
      printf("\033[92m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); // EZERO, 
    }/*else if( 450 >= idDetMap[i] && idDetMap[i] >= 400){
      printf("\033[93m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); // EZERO, 
    }else if(  99 >= idDetMap[i] && idDetMap[i] >= 0){    
      switch (idKindMap[i]) {
         case -1: printf("%7s|", ""); break;
         case  0: printf("\033[31m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); break; // RED
         case  1: printf("\033[32m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); break; // Green
         case  2: printf("\033[33m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); break; // Yellow
         case  3: printf("\033[34m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); break; // Blue
         case  4: printf("\033[35m%3d(%2d)\033[0m|", idDetMap[i], idKindMap[i]); break; // Magenta
         default: printf("%3d(%2d)|", idDetMap[i], idKindMap[i]); break; // no color
       }
     }*/
     else{
       printf("%7s|", "");
     }
  }
  printf("\n");
  gClock.Reset();
  gClock.Start("timer");

  printf("====================== started \n");

}

void GeneralSort::SlaveBegin(TTree * /*tree*/)
{
  
  TString option = GetOption();
  
}

Bool_t GeneralSort::Process(Long64_t entry)
{ 
  //=============================== Get Run Number
  if( entry == 0 ) {
      fileNum = fChain->GetDirectory()->GetName();
      
      //printf("----------------------- openning  %s \n", fileNum.Data());
      
      int findslat = fileNum.Last('/');
      fileNum.Remove(0, findslat+1);
      int found = fileNum.First(".");
      fileNum.Remove(found);
      //the fileNum should be something like "xxx_run4563" now
      while( !fileNum.IsDigit() ){
         fileNum.Remove(0,1);
      }
      runIDPresent = fileNum.Atoi();
   } 
   music.runID = runIDPresent; 
  
  
  
  ProcessedEntries++;
  if (ProcessedEntries<MaxProcessedEntries) {

    //=============================== Zero struct
    for (Int_t i=0;i<16;i++) {//num dets
      music.de_l[i]=TMath::QuietNaN();
      music.de_r[i]=TMath::QuietNaN();
      music.seg[i]=TMath::QuietNaN();

      music.de_l_TimeStamp[i]=0;
      music.de_r_TimeStamp[i]=0;
     }
     music.stp0=TMath::QuietNaN();
     music.stp17=TMath::QuietNaN();
     music.cath=TMath::QuietNaN();
     music.grid=TMath::QuietNaN();
     music.tac=TMath::QuietNaN();

     music.stp0_TimeStamp=0;
     music.stp17_TimeStamp=0;
     music.cath_TimeStamp=0;
     music.grid_TimeStamp=0;
     music.tac_TimeStamp=0;
     
     music.multi_l = 0;
     music.multi_r = 0;
     
    
    //=============================== Pull needed entries
    b_NumHits->GetEntry(entry);
    b_id->GetEntry(entry);
    b_pre_rise_energy->GetEntry(entry);
    b_post_rise_energy->GetEntry(entry);
    //   b_base_sample->GetEntry(entry);
    //    b_baseline->GetEntry(entry);
    b_event_timestamp->GetEntry(entry);

    //ID MUSIC Channels
    Int_t idKind = -1;
    Int_t idDet=-1; // Detector number
    
    
    //==============================================================
    /* --------------------- Loop over NumHits ------------------ */
    //==============================================================
    for (Int_t i=0;i<NumHits;i++) {
      Int_t idTemp = id[i] - idConst; // (idConst = 1010)
      idDet = idDetMap[idTemp];
      idKind = idKindMap[idTemp];
      
      
      ///=============================== MUSIC
      if ( 100 <= idDet && idDet < 116 && idKind < 2 ) {         
        switch(idKind)
          {
          case 0: // Left
            music.de_l[idDet-100] = (post_rise_energy[i]-pre_rise_energy[i])*1.0/M;
            music.de_l_TimeStamp[idDet-100] = event_timestamp[i];
            music.multi_l ++;
            break;
          case 1: // Right
            music.de_r[idDet-100] = ((float)(post_rise_energy[i])-(float)(pre_rise_energy[i]))/M;
            music.de_r_TimeStamp[idDet-100] = event_timestamp[i];
            music.multi_r ++;
            break;
          default:
            ;
            break;// 
          }
          music.seg[idDet-100] = idDet-100;
          
      }

      ///=============================== Strip 0
      if ( idDet == 200 ) {   
        music.stp0 = ((float)(post_rise_energy[i])-(float)(pre_rise_energy[i]))/M;
        music.stp0_TimeStamp = event_timestamp[i];
      }
       
      ///=============================== Strip 17
      if ( idDet == 201 ) {
        music.stp17 = ((float)(post_rise_energy[i])-(float)(pre_rise_energy[i]))/M*(-1);
        music.stp17_TimeStamp = event_timestamp[i];
      }

      ///=============================== Grid
      if ( idDet == 202 ) {
        music.grid = ((float)(post_rise_energy[i])-(float)(pre_rise_energy[i]))/M ;
        music.grid_TimeStamp = event_timestamp[i];
      }
      
      ///=============================== TAC
      /*if ( idDet == 203 ) {
        music.tac = ((float)(post_rise_energy[i]) -(float)(pre_rise_energy[i]))/M;
        music.tac_TimeStamp = event_timestamp[i];
      }*/

       ///=============================== CATHODE
      if ( idDet == 203 ) {
        music.cath = ((float)(post_rise_energy[i]) -(float)(pre_rise_energy[i]))/M;
        music.cath_TimeStamp = event_timestamp[i];
      }

      
    } // End NumHits Loop

    //Progress display
    /************************************************************************/
    oFile->cd(); //set focus on this file
    gen_tree->Fill();  

    gClock.Stop("timer");
    Double_t time = gClock.GetRealTime("timer");
    gClock.Start("timer");

    if ( !shown ) {
      if (fmod(time, 10) < 1 ){
         printf( "%10lld[%2d%%]|%3.0f min %5.2f sec | expect:%5.2f min\n", 
               entry, 
               TMath::Nint((entry+1)*100./EffEntries),
               TMath::Floor(time/60.), 
               time - TMath::Floor(time/60.)*60.,
               EffEntries*time/(entry+1.)/60.);
         shown = 1;
         gen_tree->Write();
      }
    }else{
      if (fmod(time, 10) > 9 ){
         shown = 0;
      }
    }

  }  
  return kTRUE;
}

void GeneralSort::SlaveTerminate()
{

}

void GeneralSort::Terminate()
{
  oFile->cd();
  gen_tree->Write();
  
  int savedEntries = gen_tree->GetEntries();
  oFile->Close();
  
  printf("=======================================================\n");
  printf(" Total processed entries : %3.1f k/%3.1f k [%4.1f%%] \n",EffEntries/1000.0, NumEntries/1000., EffEntries*100./NumEntries);
  gClock.Stop("timer");
  Double_t time = gClock.GetRealTime("timer");
  printf(" Total run time : %6.0f sec \n", time);
  printf(" Sorting rate   : %6.3f k/sec\n",EffEntries/time/1000.0);
  printf(" saved as \e[31m %s \e[m. events saved: %d\n", saveFileName.Data() , savedEntries); 
  printf("=======================================================\n");
  gROOT->ProcessLine(".q");
}
