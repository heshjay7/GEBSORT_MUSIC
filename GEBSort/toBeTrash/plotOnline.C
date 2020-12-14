void plotOnline(int iRun=0){

  gStyle->SetOptStat(0);

  std::vector<int> runlist;

  TSystemDirectory *aDirectory = new TSystemDirectory("files","ROOT_FILES");

  double alphaX[] = {7467.37,7472.76,7629.18,7618.4,7467.37};
  double alphaY[] = {13.0591,8.83966,8.89241,13.0063,13.0591};
  TCutG *alphaCut = new TCutG("alphaCut",5,alphaX,alphaY);

  double isomerX[] = {73.0122,215.652,510.602,1231.05,1083.58,554.119,73.0122,73.0122};
  double isomerY[] = {9.64599,10.006,10.03,9.64599,5.99831,5.44636,6.0943,9.64599};
  TCutG *isomerCut = new TCutG("isomerCut",8,isomerX,isomerY);
  double recoilX1[] = {10786.9,12912.1,14963.9,19873.8,20093.6,10713.7,10786.9};
  double recoilY1[] = {-461.794,-461.066,-460.339,-460.339,-466.525,-468.345,-461.794};
  TCutG *recoilCut1 = new TCutG("recoilCut1",7,recoilX1,recoilY1);
  double recoilX2[] = {12032.7,13938,17382.2,20020.3,20093.6,17895.2,11959.4,12032.7};
  double recoilY2[] = {-515.536,-511.992,-510.135,-511.316,-514.186,-520.599,-520.937,-515.536};
  TCutG *recoilCut2 = new TCutG("recoilCut2",8,recoilX2,recoilY2);

  runlist.clear();
  TIter next(aDirectory->GetListOfFiles());
  TSystemFile *sfile = NULL;
  int runnr;
  while ((sfile = (TSystemFile*)next())){
    if(sscanf(sfile->GetName(),"run%d.root",&runnr))
    runlist.push_back(runnr);
  }
  for(int i=0;i<runlist.size()-1;i++){//For our purposes, we don't need to worry about using bubble sort
    if(runlist[i]>runlist[i+1]){
	int tmp = runlist[i+1]; 
	runlist[i+1]=runlist[i];
	runlist[i]=tmp;
	i=0;
    }
  }
  TFile *file = NULL;
  int lastRun = 0;
  lastRun = runlist.back();
 // if(iRun==0) iRun = lastRun; 

  file = new TFile(Form("ROOT_FILES/run%d.root",iRun));
std::cout << ((TH1D*)file->Get("decay_rate"))->GetNbinsX() << std::endl ;  std::cout << "Run #" << iRun << " ( analyzed = " << ((TH1D*)file->Get("decay_rate"))->GetNbinsX() / 60. << " min ) " << std::endl;
  TH1D *h1 = (TH1D*)file->Get("wheel_left");
  TH2D *h2 = (TH2D*)file->Get("dssd_en");
  TH2D *h3 = (TH2D*)file->Get("eftof");
  TH2D *h4 = (TH2D*)file->Get("e1t1log");
  h1->Draw();
  c1->Clear();
  c1->Modified();c1->Update();
  gSystem->ProcessEvents();
  c1->Divide(2,2);
  c1->cd(1);
  h1->Draw();
  h1->GetXaxis()->SetRangeUser(0,1100); 
  c1->cd(2);
  h2->Draw("colz");
  gPad->SetLogz();
  c1->cd(3);
  h3->Draw("colz");
  h3->GetXaxis()->SetRangeUser(0,30000); 
  h3->GetYaxis()->SetRangeUser(-480,-450); 
  //recoilCut1->Draw("same");
  gPad->SetLogz();
  c1->cd(4);
  h4->Draw("colz");
  h4->GetXaxis()->SetRangeUser(0,12000); 
  h4->GetYaxis()->SetRangeUser(6,14);
  //alphaCut->Draw("same");
  //isomerCut->Draw("same");
  gStyle->SetOptStat(1111111);
  c1->SetLogz(); 
  c1->Modified();c1->Update();
}
