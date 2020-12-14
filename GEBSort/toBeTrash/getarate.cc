void getarate(){


TCanvas *c1 = new TCanvas("c1","c1");
 c1->Clear();
 c1->Divide(1,2);
 c1->cd(1);
 dssd_hitxy_phys_corrf->Draw("colz");
// dssd_hitxy_phys_corrf->GetXaxis()->SetRange(low,high);
// dssd_hitxy_phys_corrf->GetYaxis()->SetRange(20,120);
 
 c1->cd(2);
 recoil_rate->Draw();
 recoil_rate->GetXaxis()->SetRange(0,10000);

//decay_rate->Draw();
//pra->Draw();

	
}
