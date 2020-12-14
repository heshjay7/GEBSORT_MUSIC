void int ppac(int low, int high) {


 gStyle->SetOptStat(0000000);
 c1->Update();

 c1->Clear();
 c1->Divide(2,3);

 c1->cd(1);
 decay_rate->GetXaxis()->SetRange(low,high);
 cl->Draw();
 c1->cd(2);
 cr->GetXaxis()->SetRange(low,high);
 cr->Draw();
 c1->cd(3);
 cx->GetXaxis()->SetRange(low,high);
 cx->Draw();
 c1->cd(4);
 ppacde->GetXaxis()->SetRange(low,high);
 ppacde->Draw();
 c1->cd(5);
 //FP_rate_right->GetXaxis()->SetRange(low,high);
 clr->Draw("colz");
 c1->cd(6);
 //xde->GetXaxis()->SetRange(low,high);
 xde->Draw("colz");

}
