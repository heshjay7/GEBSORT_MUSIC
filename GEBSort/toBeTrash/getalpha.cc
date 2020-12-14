void getalpha(){

int low = 2000;
int high = 2700;





TCutG *cutg = new TCutG("cutA",12);
   cutg->SetVarX("e1t1log");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetLineColor(2);
   cutg->SetLineWidth(2);
   cutg->SetPoint(0,6265.39,3.80934);
   cutg->SetPoint(1,8083.51,3.87992);
   cutg->SetPoint(2,8008.77,6.44432);
   cutg->SetPoint(3,7462.52,7.02072);
   cutg->SetPoint(4,7401.37,9.40867);
   cutg->SetPoint(5,7173.09,9.70276);
   cutg->SetPoint(6,6955.67,9.42044);
   cutg->SetPoint(7,6368.66,10.6791);
   cutg->SetPoint(8,6167.55,10.7967);
   cutg->SetPoint(9,6242.28,4.95038);
   cutg->SetPoint(10,6242.28,4.95038);
   cutg->SetPoint(11,6265.39,3.80934);
   cutg->Draw("");

e1t1log->ProjectionX("pra",0,-1,"[cutA]");
//pra->Draw();
TCanvas *c1 = new TCanvas("c1","c1");
 c1->Clear();
 c1->Divide(1,2);
 c1->cd(1);
 e1t1log->Draw("colz");
 cutA->Draw("SAME");
 e1t1log->GetXaxis()->SetRange(low,high);
 e1t1log->GetYaxis()->SetRange(20,120);
 
 c1->cd(2);
 
 pra->Draw();
 pra->GetXaxis()->SetRange(low,high);

//decay_rate->Draw();
//pra->Draw();

	
}
