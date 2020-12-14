void int dssd() {


 gStyle->SetOptStat(0000000);
 c1->Update();

 c1->Clear();
 c1->Divide(2,4);

 c1->cd(1);
 //decay_rate->GetXaxis()->SetRange(low,high);
 dssd_en->Draw("colz");
 c1->cd(2);
 //cr->GetXaxis()->SetRange(low,high);
 dssd_en_raw->Draw("colz");
 c1->cd(3);
 //cx->GetXaxis()->SetRange(low,high);
 d_fr_e->Draw("colz");
 c1->cd(4);
 //ppacde->GetXaxis()->SetRange(low,high);
 d_ba_e->Draw("colz");
 c1->cd(5);
 //FP_rate_right->GetXaxis()->SetRange(low,high);
 r_fr_e->Draw("colz");
 c1->cd(6);
 //xde->GetXaxis()->SetRange(low,high);
 r_ba_e->Draw("colz");
 c1->cd(7);
 //xde->GetXaxis()->SetRange(low,high);
 sibox->Draw("colz");
 c1->cd(8);
 //xde->GetXaxis()->SetRange(low,high);
 e1t1log->Draw("colz");

}
