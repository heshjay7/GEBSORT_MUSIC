{
bar = new TControlBar("vertical", "GRETINA GEB");

bar->AddButton("...............compile/load GSUtil...............",".L GSUtil.cc++");
bar->AddButton("load root file GTDATA/wsi.root","dload(\"GTDATA/wsi.root\")");
bar->AddButton("load root file GTDATA/nsi.root","dload(\"GTDATA/nsi.root\")");
bar->AddButton("load map file GTDATA/test.root","dload(\"GTDATA/test.root\")");
bar->AddButton("update spectra","update()", "Click Here to update spectra");
bar->AddButton("ls()","ls()");
bar->AddButton("mkcanvas","mkcanvas()");
bar->AddButton("display 2D dtbtev", "d2(\"dtbtev\")");
bar->AddButton("display 1D dtbtev", "pjx(\"dtbtev\",\"1D dtbtev\")");


bar->AddButton("update spectra","update()", "Click Here to update spectra");
bar->AddButton("display CCsum spectrum", "d1(\"CCsum\",0,1500)");
bar->AddButton("display hit pattern", "d1(\"hitpat\",0,130)");
bar->AddButton("display radius_all spectrum", "d1(\"radius_all\")");
bar->AddButton("display SMAP_allhits", "d2(\"SMAP_allhits\")");
bar->AddButton("display rate_mode2_min", "d1(\"rate_mode2_min\")");
bar->AddButton("display CCe matrix","d2(\"CCe\")");
/*bar->AddButton("findGT_cal_CC","findGT_cal_CC(\"207Bi\",\"CCenergy.cal\")", "Click Here to calibrate CC");*/

bar->AddButton("display fm spectrum", "d1(\"fm\")");
bar->AddButton("display sumTrackE FOM<0.5", "pjx(\"fomXe\",\"x\",0,0.5);d1(\"x\",0,1500)");
bar->AddButton("display sumTrackE FOM<0.6", "pjx(\"fomXe\",\"x\",0,0.6);d1(\"x\",0,1500)");
bar->AddButton("display sumTrackE FOM<0.7", "pjx(\"fomXe\",\"x\",0,0.7);d1(\"x\",0,1500)");
bar->AddButton("display sumTrackE FOM<0.8", "pjx(\"fomXe\",\"x\",0,0.8);d1(\"x\",0,1500)");
bar->AddButton("display sumTrackE FOM<1.0", "pjx(\"fomXe\",\"x\",0,1.0);d1(\"x\",0,1500)");
bar->AddButton("display sumTrackE FOM<2.0", "pjx(\"fomXe\",\"x\",0,2.0);d1(\"x\",0,1500)");
bar->AddButton("rebin","x->Rebin(2)");
bar->AddButton("smooth","x->Smooth(1)");

bar->AddButton("display SMAP_firsthits", "d2(\"SMAP_firsthits\")");
bar->AddButton("display rate_mode1", "d1(\"rate_mode1\")");

bar->AddButton("wrspe x", "wrspe(\"x\",\"x.spe\")");
bar->AddButton("curve.cc", ".x curve.cc");
bar->AddButton("quit", ".q");

bar->SetTextColor("blue");
bar->Show();
}
