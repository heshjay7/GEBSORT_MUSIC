void int monitor(int low, int high) {


 gStyle->SetOptStat(0000000);
 c1->Update();

 c1->Clear();
 c1->Divide(1,7);

 c1->cd(1);
 decay_rate->GetXaxis()->SetRange(low,high);
 decay_rate->Draw();
 c1->cd(2);
 recoil_rate->GetXaxis()->SetRange(low,high);
 recoil_rate->Draw();
 c1->cd(3);
 emonl_rate->GetXaxis()->SetRange(low,high);
 emonl_rate->Draw();
 c1->cd(4);
 emonr_rate->GetXaxis()->SetRange(low,high);
 emonr_rate->Draw();
 c1->cd(5);
 FP_rate_left->GetXaxis()->SetRange(low,high);
 FP_rate_left->Draw();
 c1->cd(6);
 FP_rate_right->GetXaxis()->SetRange(low,high);
 FP_rate_right->Draw();
 c1->cd(7);
 corr_rate->GetXaxis()->SetRange(low,high);
 corr_rate->Draw();

}
