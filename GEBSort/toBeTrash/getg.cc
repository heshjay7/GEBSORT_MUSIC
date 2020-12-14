void int getg(int low, int high) {

 corr_gammas->ProjectionY("p1",695,710);
 corr_gammas->ProjectionY("p2",679,688);
 corr_gammas->ProjectionY("p3",665,672);
 corr_gammas->ProjectionY("p4",629,637);
 corr_gammas->ProjectionY("p5",594,603);
 corr_gammas->ProjectionY("p6",725,738);

 
 TCanvas *c1 = new TCanvas("c1","c1");
 c1->Clear();
 c1->Divide(1,6);
 c1->cd(1);
 p1->GetXaxis()->SetRange(low,high);
 p1->Draw();
 // decay_rate->Draw();
 c1->cd(2);
 p2->GetXaxis()->SetRange(low,high);
 p2->Draw();
 // recoil_rate->Draw();
 c1->cd(3);
 p3->GetXaxis()->SetRange(low,high);
 p3->Draw();
 //emonl_rate->Draw();
 c1->cd(4);
 p4->GetXaxis()->SetRange(low,high);
 p4->Draw();
 // left_rate->Draw();
 c1->cd(5);
 p5->GetXaxis()->SetRange(low,high);
 p5->Draw();
 c1->cd(6);
 p6->GetXaxis()->SetRange(low,high);
 p6->Draw();
}
