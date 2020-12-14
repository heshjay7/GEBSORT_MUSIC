void delg(int p1, int p2, int d1, int d2, int low, int high) {


 //gStyle->SetOptStat(0000000);
 //c1->Update();

dTgdssdr_clenergy->ProjectionY("p1",p1, p2);
dTgdssdr_clenergy->ProjectionY("p2",d1, d2);

TH1F *h1, *h2;


h1 = (TH1F*) gROOT->FindObject("p1");
h2 = (TH1F*) gROOT->FindObject("p2");

h1->Add(h2,-1);

 c1->Clear();
h1->GetXaxis()->SetRange(low,high);
h1->Draw();

}

 
