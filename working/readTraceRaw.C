#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TClonesArray.h>
#include <TBenchmark.h>
#include <TMath.h>
#include <TF1.h>
#include <TMath.h>
#include <TLine.h>

  TGraph * gTrace;
  
  bool plotTrace = false;

void readTraceRaw(TString fileName){
   
/**///==============================================================   

   TFile * f1 = new TFile (fileName, "read");
   TTree * tree = (TTree *) f1->Get("tree");
   
   int totnumEntry = tree->GetEntries();
   printf( "========== total Entry : %d \n", totnumEntry);
   
   TCanvas * cRead = new TCanvas("cRead", "Read Trace", 0, 0, 1200, 600);
   cRead->Divide(1,1);
   for( int i = 1; i <= 2 ; i++){
      cRead->cd(i)->SetGrid();
   }
   cRead->SetGrid();
   
   
   TH1F * hppac = new TH1F("hppac", "hppac", 500, 0, 3000);
   
/**///==============================================================   

  Short_t         trace[200][1024];
  UShort_t        trace_length[200];
  Int_t           NumHits;
  tree->SetBranchAddress("trace_length", trace_length);
  tree->SetBranchAddress("trace", trace);
  tree->SetBranchAddress("NumHits", &NumHits);
   
	char s [1] ;
	char b [1] ;
	b[0] = 'q';
  
  
gTrace = new TGraph();
double base = 0;

   bool breakFlag = false;
   for( int ev = 0; ev < totnumEntry; ev++){
      tree->GetEntry(ev);
      double min = 1000000; 
      
      for(Int_t i=0;i<NumHits;i++) { 
          
        if ( trace_length[i] == 0 ) continue;
        
         gTrace->Clear();
        
         
         for ( long long int j = 0 ; j < TMath::Min(trace_length[i], (UShort_t)200); j++){
            if( trace[i][j] < 16000) {
              base = trace[i][j];
              gTrace->SetPoint(j, j, trace[i][j]);
            }else{
              gTrace->SetPoint(j, j, base);
            }
            
            if( 80 < j && j < 110){
              if( trace[i][j] < min ) min= trace[i][j];
            } 
             
          }
          
          //double min = gTrace->GetYaxis()->GetXmin();
          
          hppac->Fill((gTrace->Eval(50)) - min);

          if( plotTrace ){  
          gTrace->SetTitle(Form("min : %f ", min));
          
           cRead->cd(1);
           cRead->Clear();
            gTrace->GetYaxis()->SetRangeUser(6000, 8500);
            gTrace->SetMarkerStyle(3);
            gTrace->SetMarkerColor(2);
           gTrace->Draw("AP");
           cRead->Update();
           gSystem->ProcessEvents();
          
           gets(s);    
           if( s[0] == b[0] ) {
              breakFlag = true;
              break;
           }
          
          }
        
      }
      

   }
   
   cRead->cd(1);
   hppac->Draw();
   
   //gROOT->ProcessLine(".q");

}
