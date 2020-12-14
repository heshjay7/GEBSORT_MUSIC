#define Monitors_cxx

#include "Monitors.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
#include <fstream>
#include <TMath.h>
#include <vector>
#include <TCanvas.h>
using namespace std;

//===================== User setting

int rawERange[3] = {400, 0, 250}; 

//==================== end

///======== histograms

///====== Raw data
TH1F* hcath;
TH1F* hgrid;
TH1F* htac;
TH1F* hstrip0;
TH1F* hstrip17;

TH1F*  hleft[16];
TH1F*  hright[16];

TH2F*  hL_seg;
TH2F*  hR_seg;

TH1F* hRate_l[16];
TH1F* hRate_r[16];

TH2F* hStp0Grid;
TH2F* hStp0Cathode;

TH1F * hTDiff;
TH1F * hMulti;

TH1F * hTDiffMax;
TH2F * hMultiTDiff;

///=========== anlayzed data
TH2F* hSum_seg;
TH2F* hSum_seg_g;

///======== 
TString canvasTitle;
int lastRunID;
bool contFlag;

int totEntries;
ULong64_t startTime, endTime;

///========= Calibration variables

void Monitors::Begin(TTree * tree)
{
  TString option = GetOption();
  
  canvasTitle.Form("Runs: ");
  lastRunID = -1;
  contFlag = false;
  
  ///==============================
  totEntries = tree->GetEntries();
  printf("==========================================\n");
  printf("total number of entries : %d \n", totEntries);
  
  ULong64_t time[16];
  tree->SetBranchAddress("de_l_t", time);
  
  bool timeFound = false;
  int checkedEvt = 0;
  
  for( int i = 0 ; i < totEntries; i++){
    tree->GetEntry(i);
    checkedEvt ++;
    for( int ch = 0; ch < 16 ; ch ++ ){
      startTime = time[ch];
      if( time[ch] != 0 ) {
        timeFound = true;
        break;
      }
    }
    if( timeFound == true) break;
  }
  
  timeFound = false;
  for( int i = totEntries - 1 ; i > 0; i--){
    tree->GetEntry(i);
    checkedEvt ++;
    for( int ch = 0; ch < 16 ; ch ++ ){
      endTime = time[ch];
      if( time[ch] != 0 ) {
        timeFound = true;
        break;
      }
    }
    if( timeFound == true) break;
  }
  
  int timeRange[2];
  timeRange[0] = startTime/1e8;
  timeRange[1] = endTime/1e8;
  
  double dTime = timeRange[1] - timeRange[0];
  
  timeRange[0] = timeRange[0] - 0.1 * dTime;
  timeRange[1] = timeRange[1] + 0.1 * dTime;

  printf("----------------------------\n");
  printf(" searched timeStamp in %d events \n", checkedEvt);
  printf("  start Time : %llu ch\n", startTime);
  printf("    end Time : %llu ch\n", endTime);
  printf("----------------------------\n");
  printf("  total time : %f sec = %f min \n", dTime, dTime/60.);
  
  //timeRange[0] = 400;
  //timeRange[1] = 800;
  

  ///======Definition of histograms; 
  hcath = new TH1F("hcath", "cath",rawERange[0],rawERange[1],1000);
  hgrid = new TH1F("hgrid", "grid",rawERange[0],-500,5000);
  htac  = new TH1F("htac",   "tac",rawERange[0],rawERange[1],rawERange[2]);
  hstrip0  = new TH1F("hstrip0",   "strip0",rawERange[0],-500,100);
  hstrip17 = new TH1F("hstrip17", "strip17",rawERange[0],-100,500);

  for(int i=0;i<16;i++){
    hleft[i] = new TH1F(Form("hleft%d",i), Form("left%d; abs(eng)",i), rawERange[0],rawERange[1],rawERange[2]);
    hright[i]= new TH1F(Form("hright%d",i),Form("right%d; abs(eng)",i),rawERange[0],rawERange[1],rawERange[2]);
    
    hRate_l[i] = new TH1F(Form("hRate_l%d",i), Form("rate id = %d", i), timeRange[1]-timeRange[0], timeRange[0], timeRange[1]);
    hRate_r[i] = new TH1F(Form("hRate_r%d",i), Form("rate id = %d", i), timeRange[1]-timeRange[0], timeRange[0], timeRange[1]);
  }
  
  hL_seg = new TH2F("hL_strip", "Left_vs_seg; id; abs(eng)",    16, 0, 16, rawERange[0],rawERange[1],rawERange[2]);
  hR_seg = new TH2F("hR_strip", "Right_vs_seg; id; abs(eng)", 16, 0, 16, rawERange[0],rawERange[1],rawERange[2]);
 
  hStp0Grid = new TH2F("hStp0Grid",  "Stp0 vs Grid; grid ; Stp0", rawERange[0], rawERange[1], rawERange[2], rawERange[0], rawERange[1], rawERange[2]); 
  hStp0Cathode = new TH2F("hStp0Cath","Stp0 vs Cath; Cath; Stp0", rawERange[0], rawERange[1], rawERange[2], rawERange[0], rawERange[1], rawERange[2]); 
  
  
  hTDiff    = new TH1F("hTDiff",    "Time diff in an event; [10 ns]", 3000, 0, 3000);
  hTDiffMax = new TH1F("hTDiffMax", "Max Time diff in an event; [10 ns]", 3000, 0, 3000);
  hMulti    = new TH1F("hMulti",    "Multiplicity", 40, 0, 40);
  
  hMultiTDiff = new TH2F("hMultiTDiff", "Multiplicity vs Max TimDiff; [10 ns]; Multiplicity", 3000, 0, 3000, 40, 0, 40);
  
  hSum_seg   = new  TH2F("hSum_seg",          "Sum_vs_seg; id; abs(eng)",    16, 0, 16, rawERange[0],rawERange[1],2*rawERange[2]);
  hSum_seg_g = new  TH2F("hSum_seg_g", "Sum_vs_seg(stp17 < -150); id; abs(eng)",    16, 0, 16, rawERange[0],rawERange[1],2*rawERange[2]);
  
  ///================================= load Calibration file
  
  
}

Bool_t Monitors::Process(Long64_t entry)
{
  b_runID->GetEntry(entry);
  b_multi_l->GetEntry(entry);
  b_multi_r->GetEntry(entry);
  b_de_l->GetEntry(entry);
  b_de_r->GetEntry(entry);
  b_de_l_TimeStamp->GetEntry(entry);
  b_de_r_TimeStamp->GetEntry(entry);
  
  b_seg->GetEntry(entry);
  b_stp0->GetEntry(entry);
  b_stp17->GetEntry(entry);
  b_cath->GetEntry(entry);
  b_grid->GetEntry(entry);
  b_tac->GetEntry(entry);
  
   /*********** forming canvas Title **********************************/ 
  if( entry == 0 ) {
     if( runID == lastRunID + 1 ) {
        int len = canvasTitle.Sizeof();
        if( contFlag == false) {
           canvasTitle.Remove(len-3);
           canvasTitle += " - ";
        }
        if( contFlag == true){
           canvasTitle.Remove(len-6);
        }
        contFlag = true;
     }
     if( runID > lastRunID + 1 ) contFlag = false;
     canvasTitle += Form("%03d, ", runID );
     lastRunID = runID;
  }

  /****************** Gate *****************************************/
  int multi = multi_l + multi_r;
      
  //if( multi < 16 ) return kTRUE;
  
  //if( grid < 2200 ) return kTRUE;
  
  /****************** Analysis *****************************************/
  
  hcath->Fill(cath);
  hgrid->Fill(grid);
  htac->Fill(tac);
  hstrip0->Fill(stp0);
  hstrip17->Fill(stp17);

  hStp0Grid->Fill(grid,stp0);
  hStp0Cathode->Fill(cath,stp0);
  
  ULong64_t timeRef = 0;
  int timeRefID = -1;
  int maxTimeDiff = 0;
  
  for(int i=0; i<16; i++) {
    
    hleft[i]->Fill(abs(de_l[i]));
    hright[i]->Fill(abs(de_r[i]));

    hL_seg->Fill(i,abs(de_l[i]));
    hR_seg->Fill(i,abs(de_r[i]));
    
    int tDiff = 0;
    
    if( de_l_t[i]>0 ) {
      
      hRate_l[i]->Fill(de_l_t[i]/1e8);
    
      if( timeRefID == -1  ) {
        timeRef = de_l_t[i]; timeRefID = i;
      }

      if( timeRefID != 1 ) {
        tDiff = de_l_t[i] > timeRef ? de_l_t[i] - timeRef : timeRef - de_l_t[i];
        hTDiff->Fill( tDiff );
        
        if( tDiff > maxTimeDiff ) maxTimeDiff = tDiff;
      }
    }
    
    if(  de_r_t[i]>0 ) {
      
      hRate_r[i]->Fill(de_r_t[i]/1e8);
    
      if( timeRefID == -1  ) {timeRef = de_r_t[i]; timeRefID = i;}

      
      if( timeRefID != 1 ) {
        tDiff = de_r_t[i] > timeRef ? de_r_t[i] - timeRef : timeRef - de_r_t[i];
        hTDiff->Fill( tDiff );
        
        if( tDiff > maxTimeDiff ) maxTimeDiff = tDiff;
      }
    
    }
  
    double scale = 1.0;
    
    hSum_seg->Fill( i, abs(de_l[i]) + scale * abs(de_r[i]));
    if( stp17 > -50) hSum_seg_g->Fill( i, abs(de_l[i]) + scale * abs(de_r[i]));
  }
  
  hMulti->Fill(multi);
  hTDiffMax->Fill(maxTimeDiff);
  hMultiTDiff->Fill(maxTimeDiff, multi);
  
  
  return kTRUE;
}


void Monitors::Terminate()
{
  printf("============ finishing.\n");

  gROOT->cd();
  gStyle->SetOptStat("neiou");
  
  int strLen = canvasTitle.Sizeof();
  canvasTitle.Remove(strLen-3);

  ///----------------------------------- Canvas - 1
  TCanvas * cSeg  = new TCanvas("cSeg","cSeg | " + canvasTitle,1000,600);
  cSeg->cd(); 
  cSeg->Divide(1,2);

  cSeg->cd(1);
  hL_seg->Draw("colz");

  cSeg->cd(2);
  hR_seg->Draw("colz");
  
  TCanvas * cRate = new TCanvas("cRate","cRate | " + canvasTitle, 0, 700 , 1000,600);
  cRate->cd(); 
  cRate->Divide(8,4);
  for( int i = 0 ; i < 32 ; i++){
    cRate->cd(i+1);
    if( i < 16 ) hRate_l[i]->Draw();
    if( i >= 16 ) hRate_r[i-16]->Draw();
  }

	///----------------------------------- Canvas - 2
/*
	TCanvas * cLeft = new TCanvas("hLeft","hLeft",1500,1500);
	cLeft->cd(); 
	cLeft->Divide(4,4);

	for(int i=0; i<16; i++){
		cLeft->cd(i+1);
		hleft[i]->SetTitle(Form("hLeft_%d",i+1));
		hleft[i]->Draw();
	}

	///----------------------------------- Canvas - 3

	TCanvas * cRight = new TCanvas("hRight","hRight",1500,1500);
	cRight->cd(); 
	cRight->Divide(4,4);

	for(int i=0; i<16; i++){
		cRight->cd(i+1);
		hright[i]->SetTitle(Form("hRight_%d",i+1));
		hright[i]->Draw();
	}

	///----------------------------------- Canvas - 4
*/
	TCanvas * cIndv = new TCanvas("hIndv","hIndv | " + canvasTitle, 1100, 0, 1000,1000);
	cIndv->cd(); 
	cIndv->Divide(2,3);

	cIndv->cd(1);	hstrip0->Draw("");

	cIndv->cd(2);	hstrip17->Draw("");

	cIndv->cd(3);	hgrid->Draw("");

	cIndv->cd(4);	hcath->Draw("");
  
	cIndv->cd(5); hStp0Grid->Draw("colz");

	cIndv->cd(6); hStp0Cathode->Draw("colz");
  
  ///----------------------------------- Canvas - 4
  TCanvas * cSum = new TCanvas("cSum","cSum | " + canvasTitle, 500, 0, 1000,1000);
	cSum->cd(); 
	cSum->Divide(2,2);
  
  cSum->cd(1); hMulti->Draw("colz");
  cSum->cd(2); hTDiff->Draw("colz");
  
  cSum->cd(3); hMultiTDiff->Draw("colz");
  cSum->cd(4); hTDiffMax->Draw("colz");

  
  //cSum->cd(1); hSum_seg->Draw("colz");
  //cSum->cd(2); hSum_seg_g->Draw("colz");

}

void Monitors::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
}
void Monitors::SlaveTerminate()
{

}
