{
 gROOT->Macro("../code/cm/logon.C+");
 int bin  = 8;
 int bin1 = 7;
 int binPt = 5;
 int binRap = 5;
 float cent[8]  ={22.1, 86.3, 130.0, 187.1, 261.4, 329.4, 381.3}; // with 40-50 and 50-100
 float cent1[8] ={14.2, 69.9, 130.0, 187.1, 261.4, 329.4, 381.3}; // with 40-60 and 60-100
 float cent2[8] ={17.8, 69.9, 130.0, 187.1, 261.4, 355.4};
 float pt [5] = {1.25, 3.75, 6.5, 10., 16.};
 float pte[5] = {1.25, 1.25, 1.5, 2., 4.};
 float rap[5] = {0.2 , 0.55, 0.85, 1.25, 1.95};
 float rape[5]= {0.2, 0.15, 0.15, 0.25, 0.45};
 //float cent[bin]={22.1, 86.3, 130.0, 187.1, 261.4, 329.4, 381.3};
 //float cent1[8]={14.2, 53.5, 86.3, 130.0, 187.1, 261.4, 329.4, 381.3};
 //float cent2[6]={22.1, 86.3, 130.0, 187.1, 261.4, 355.4};
 float centErr[bin]={6,6,6,6,6,6,6};
 float centnoErr[bin]={0,0,0,0};

 //alice shit pT>0 GeV/c 
 int binA=2;
 // float centAlice[binA]={72,308};
 // float centAliceErr[binA]={6,6};
 // float centAliceNoErr[binA]={0,0};
 // float ratioUpsAlice[binA]={0.634,0.341};
 // float ratioUpsAliceStat[binA]={0.120,0.074};
 // float ratioUpsAliceSyst[binA]={0.044,0.024};

 float centAlice[binA]={72,308};
 float centAliceErr[binA]={6,6};
 float centAliceNoErr[binA]={0,0};
 float ratioUpsAlice[binA]    ={0.512,0.351}; // 1S with 'alice' binning
 float ratioUpsAliceStat[binA]={0.036,0.026};
 float ratioUpsAliceSyst[binA]={0.028,0.029};

 float ratioUpsAlice2[binA]    ={0.203,0.056}; // 2S with 'alice' binning
 float ratioUpsAliceStat2[binA]={0.055,0.036};
 float ratioUpsAliceSyst2[binA]={0.035,0.050};
 //float ratioUpsAliceGlob[binA]={};

 //Raa with regit and pp2013

 float ratio1S[bin1]       ={/*1.047,0.981,*/1.134,0.617,0.573,0.429,0.356,0.364,0.324};
 float ratio1SstatErr[bin1]={/*0.202,0.166,*/0.181,0.098,0.062,0.044,0.038,0.040,0.039};
 float ratio1SsystErr[bin1]={/*0.186,0.167,*/0.154,0.084,0.037,0.057,0.041,0.050,0.035};
 float ratio2S[bin]       ={0.484,0.292,0.327,0.141,0.035,0.036,0.109};
 float ratio2SstatErr[bin]={0.337,0.122,0.105,0.068,0.055,0.063,0.059};
 float ratio2SsystErr[bin]={0.135,0.252,0.043,0.070,0.042,0.038,0.032};

 float ratioJpsi[bin]={0.610, 0.661, 0.506, 0.372, 0.220, 0.202};
 float ratioJpsistatErr[bin]={0.122, 0.131, 0.085, 0.058, 0.037, 0.030};
 float ratioJpsisystErr[bin]={0.097, 0.077, 0.048, 0.029, 0.016, 0.014};

 float raaPt1 [5] = {0.297,0.403,0.372,0.404,0.498};
 float raaPt1e [5] = {0.033,0.037,0.041,0.044,0.152};
 float raaPt1s [5] = {0.049,0.035,0.009,0.049,0.262};
 float raaRap1 [5] = {0.362,0.384,0.417,0.497,0.414};
 float raaRap1e[5] = {0.038,0.043,0.055,0.042,0.038};
 float raaRap1s[5] = {0.039,0.041,0.055,0.026,0.062};

 float raaPt2 [5] = {0.015,0.091,0.036,0.094,0.172};
 float raaPt2e [5] = {0.047,0.062,0.056,0.066,0.077};
 float raaPt2s [5] = {0.079,0.050,0.005,0.037,0.};
 float raaRap2 [5] = {0.071,0.128,0.133,0.09,0.247};
 float raaRap2e[5] = {0.053,0.060,0.076,0.056,0.069};
 float raaRap2s[5] = {0.032,0.033,0.038,0.082,0.092};

 float tot1SErr[bin];
 float tot2SErr[bin];

 for (int i=0; i<bin1; i++){
   tot1SErr[i] = sqrt(ratio1SstatErr[i]*ratio1SstatErr[i]+ratio1SsystErr[i]*ratio1SsystErr[i]);
   cout<<i<<" "<<tot1SErr[i]<<endl;
 }
 for (int i=0; i<bin; i++){
   tot2SErr[i] = sqrt(ratio2SstatErr[i]*ratio2SstatErr[i]+ratio2SsystErr[i]*ratio2SsystErr[i]);
   cout<<i<<" "<<tot2SErr[i]<<endl;
 }

 TCanvas *c1 = new TCanvas("c1","c1");

 c1->cd();
 TPad *p1 = new TPad("p1","p1",0.0,0.0,0.9,1.0);
 p1->SetBottomMargin(0.12);
 p1->SetTopMargin(0.03);
 p1->SetRightMargin(0.03);
 p1->SetLeftMargin(0.11);
 p1->Draw();
 p1->cd();

 TF1 *f4 = new TF1("f4","1",0,400);
 f4->SetLineWidth(1);
 f4->GetXaxis()->SetTitle("N_{part}");
 f4->GetYaxis()->SetTitle("R_{AA}");
 f4->GetYaxis()->SetTitleOffset(1.0);
 f4->GetYaxis()->SetRangeUser(0,1.5);
 f4->GetXaxis()->CenterTitle(kTRUE);
 f4->Draw();

 TGraphErrors *g = new TGraphErrors(bin1,cent1,ratio1S,centnoErr,ratio1SstatErr);
 g->SetMarkerColor(2);
 g->SetMarkerStyle(33);
 g->SetMarkerSize(2);

 TGraphErrors *g1Ssyst = new TGraphErrors(bin1,cent1,ratio1S,centnoErr,ratio1SsystErr);
 g1Ssyst->SetLineColor(kRed-9);
 g1Ssyst->SetLineWidth(18);
 g1Ssyst->SetMarkerSize(0);
 g1Ssyst->Draw("e");

 g->Draw("pe");
 TGraphErrors *g1circle = new TGraphErrors(bin1,cent1,ratio1S,centnoErr,ratio1SstatErr);
 g1circle->SetMarkerStyle(27);
 g1circle->SetMarkerSize(2);
 g1circle->SetLineColor(kBlack);
 g1circle->Draw("p");

 TGraphErrors *g2 = new TGraphErrors(bin1,cent,ratio2S,centnoErr,ratio2SstatErr);
 g2->SetMarkerColor(kTeal+3);
 g2->SetMarkerStyle(20);
 g2->SetMarkerSize(1.2);

 TGraphErrors *g2Ssyst = new TGraphErrors(bin1,cent,ratio2S,centnoErr,ratio2SsystErr);
 g2Ssyst->SetLineColor(kTeal+5);
 g2Ssyst->SetFillStyle(0);
 g2Ssyst->SetLineWidth(18);
 g2Ssyst->SetMarkerSize(0);
 g2Ssyst->Draw("2");

 g2->Draw("pe");
 TGraphErrors *g2circle = new TGraphErrors(bin1,cent,ratio2S,centnoErr,ratio2SstatErr);
 g2circle->SetMarkerStyle(24);
 g2circle->SetMarkerSize(1.2);
 g2circle->SetLineColor(kBlack);
 g2circle->Draw("p");
  

 //Alice shit
 TGraphErrors *gA = new TGraphErrors(binA,centAlice,ratioUpsAlice,centAliceNoErr,ratioUpsAliceStat);
 gA->SetMarkerColor(kGray);
 gA->SetMarkerStyle(33);
 gA->SetMarkerSize(2);
	
 TGraphErrors *gAsyst = new TGraphErrors(binA,centAlice,ratioUpsAlice,centAliceNoErr,ratioUpsAliceSyst);
 gAsyst->SetLineColor(kGray);
 gAsyst->SetLineWidth(25);
 // gAsyst->SetDrawStyle(3244);
 gAsyst->SetMarkerSize(0);
 // gAsyst->Draw("e");
 //gA->Draw("pe");
	
 TGraphErrors *gAcircle = new TGraphErrors(binA,centAlice,ratioUpsAlice,centAliceNoErr,ratioUpsAliceStat);
 gAcircle->SetMarkerStyle(33);
 gAcircle->SetMarkerSize(2);
 gAcircle->SetLineColor(kBlack);
 // gAcircle->Draw("p");

 TGraphErrors *gA2 = new TGraphErrors(binA,centAlice,ratioUpsAlice2,centAliceNoErr,ratioUpsAliceStat2);
 gA2->SetMarkerColor(kTeal);
 gA2->SetMarkerStyle(33);
 gA2->SetMarkerSize(2);
	
 TGraphErrors *gAsyst2 = new TGraphErrors(binA,centAlice,ratioUpsAlice2,centAliceNoErr,ratioUpsAliceSyst2);
 gAsyst2->SetLineColor(kTeal);
 gAsyst2->SetLineWidth(25);
 // gAsyst->SetDrawStyle(3244);
 gAsyst2->SetMarkerSize(0);
 // gAsyst2->Draw("e");
 //gA2->Draw("pe");
	
 TGraphErrors *gA2circle = new TGraphErrors(binA,centAlice,ratioUpsAlice2,centAliceNoErr,ratioUpsAliceStat2);
 gA2circle->SetMarkerStyle(33);
 gA2circle->SetMarkerSize(2);
 gA2circle->SetLineColor(kBlack);
 // gA2circle->Draw("p");
	


 float ppy[1]={1};
 float ppx[1]={377};
 float ppxEr[1]={0};
 float ppyErSystematic=sqrt(pow(0.06,2)+pow(0.023,2));
 cout<<"ppyErSystematic"<<ppyErSystematic<<endl;
 float ppyErSyst[1]={ppyErSystematic};
 // float ppyErtotal=sqrt(pow(ppyErSystematic,2)+pow(0.12,2));
 float ppyErtol[1]={0};
 // TBox *lumiY1S = new TBox(355.0,0.864,370.0,1.136);
 // lumiY1S->SetFillColor(kGreen-9);
 // lumiY1S->Draw("2");
 TGraphErrors *gpp = new TGraphErrors(1,ppx,ppy,ppxEr,ppyErSyst);
 gpp->SetMarkerSize(0);
 gpp->SetLineColor(kGray);
 gpp->SetLineWidth(10);
 //gpp->Draw("e");
 TGraphErrors *gpptol = new TGraphErrors(1,ppx,ppy,ppxEr,ppyErtol);
 gpptol->SetLineWidth(0);
 gpptol->SetLineColor(kRed-9);
 gpptol->SetMarkerColor(kRed);
 gpptol->SetMarkerStyle(21);
 gpptol->SetMarkerSize(1.1);
 //gpptol->Draw("p");
	


 float Y2Sppx[1]={392};
 float Y2SppyErSystematic=sqrt(pow(0.06,2)+pow(0.033,2));
 float Y2SppyErSyst[1]={Y2SppyErSystematic};
 float Y2SppyErtotal=sqrt(pow(Y2SppyErSystematic,2)+pow(0.202,2));
 float Y2SppyErtol[1]={0};

 // TBox *lumiY2S = new TBox(370.0,0.787,385.0,1.213);
 // lumiY2S->SetFillColor(kAzure-9);
 // /* lumiY2S->SetFillStyle(0); */
 // /* lumiY2S->SetLineWidth(2); */
 // /* lumiY2S->SetLineColor(kAzure+2); */
 // lumiY2S->Draw("2");
 TGraphErrors *Y2Sgpp = new TGraphErrors(1,Y2Sppx,ppy,ppxEr,Y2SppyErSyst);
 Y2Sgpp->SetMarkerSize(0);
 Y2Sgpp->SetLineColor(kGray);
 Y2Sgpp->SetLineWidth(10);
 //Y2Sgpp->Draw("e");
 TGraphErrors *Y2Sgpptol = new TGraphErrors(1,Y2Sppx,ppy,ppxEr,Y2SppyErtol);
 Y2Sgpptol->SetLineWidth(0);
 Y2Sgpptol->SetLineColor(kGray);
 Y2Sgpptol->SetMarkerColor(kGreen+2);
 Y2Sgpptol->SetMarkerStyle(20);
 Y2Sgpptol->SetMarkerSize(1.1);
 //Y2Sgpptol->Draw("p");

 /* TBox *lumiJpsi = new TBox(340.0,0.94,355.0,1.06); */
 /* lumiJpsi->SetFillColor(kBlue-7); */
 /* lumiJpsi->SetFillStyle(3244); */
 /* lumiJpsi->SetLineWidth(2); */
 /* lumiJpsi->SetLineColor(kBlue+2); */
 /* lumiJpsi->Draw("2"); */

 /* TBox *lumiJpsi_ = new TBox(340.0,0.94,355.0,1.06); */
 /* lumiJpsi_->SetFillStyle(0); */
 /* lumiJpsi_->SetLineWidth(1); */
 /* lumiJpsi_->SetLineColor(kBlue+2); */
 /* lumiJpsi_->Draw("2"); */
	
 // TBox *globAlice = new TBox(385.0,0.712,399.2,1.259);
 // globAlice->SetFillColor(kGray);
 // globAlice->Draw("2");
  
 f4->Draw("same");
 gPad->RedrawAxis();

  TLatex *l1CMS = new TLatex(20,1.39,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  l1CMS->SetTextFont(42);
  l1CMS->SetTextSize(0.032);
  // l1CMS->Draw();
  TLatex *lyCMS = new TLatex(200,1.39,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  lyCMS->SetTextSize(0.034);	
  lyCMS->Draw();
  TLatex *lyL= new TLatex(220,1.3,"L_{int} = 150 #mub^{-1}; |y| < 2.4");
  lyL->SetTextSize(0.029);
  lyL->Draw();

  // TLatex *l1A = new TLatex(210,1.4,"ALICE PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  // l1A->SetTextFont(42);
  // l1A->SetTextSize(0.027);
  // l1A->Draw();

  // TLatex *lyA = new TLatex(20,1.3,"ALICE Preliminary");
  // lyA->SetTextSize(0.027);
  // lyA->Draw();
  // TLatex *lPr = new TLatex(20,1.3,"L_{int} = 69 #mub^{-1}; 2.5 < y < 4");
  // lPr->SetTextSize(0.027);
  // lPr->Draw();
	
  // TLatex *mupt= new TLatex(20,1.24,"p_{T}^{#mu CMS} > 4 GeV/c; p_{T}^{#mu ALICE} > 0 GeV/c");
  // mupt->SetTextSize(0.027);
  // mupt->Draw();
  TLatex *centrality = new TLatex();
  centrality->SetTextSize(0.027);
  centrality->DrawLatex(368,0.57,"0-5%");
  centrality->DrawLatex(315,0.57,"5-10%");
  centrality->DrawLatex(245,0.57,"10-20%");
  centrality->DrawLatex(170,0.82,"20-30%");
  centrality->DrawLatex(120,0.82,"30-40%");
  centrality->DrawLatex(70,0.82,"40-60%");
  centrality->DrawLatex(10,0.75,"60-100%");
  // centrality->DrawLatex(17,0.62,"20-90%");
  // centrality->DrawLatex(260,0.325,"0-20%");

  TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
  legend->SetTextSize(0.029);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->AddEntry(g,"#varUpsilon(1S)","p");
  legend->AddEntry(g2, "#varUpsilon(2S)","p");
 
  legend->Draw();
  // TLegend *legendA = new TLegend(0.15,0.75,0.4,0.84);
  // legendA->SetFillStyle(0);
  // //legendA->SetFillColor(19);
  // legendA->SetBorderSize(0);
  // legendA->SetTextFont(42);
  // legendA->SetTextSize(0.026);
  // // legendA->AddEntry(gA, "ALICE Binning #varUpsilon(1S)","p");
  // legendA->AddEntry(gA2, "ALICE Binning #varUpsilon(2S)","p");
  // legendA->Draw();


 TLine *g3Slegend = new TLine(192,1.1,200,1.1);
 g3Slegend->SetLineWidth(2);
 //g3Slegend->Draw();
 TArrow *g3SlegendArrow = new TArrow(196,1.1,196,1.03,0.015,"|>");
 g3SlegendArrow->SetLineWidth(2);
 //g3SlegendArrow->Draw();



 c1->cd();
 TPad *p2 = new TPad("p2","p2",0.9,0.0,1.0,1.0);
 p2->SetBottomMargin(0.12);
 p2->SetTopMargin(0.03);
 p2->SetRightMargin(0.2);
 p2->SetLeftMargin(0.01);
 p2->SetTickx(0);
 p2->SetTicky(0);
 p2->Draw();
 p2->cd();

 //MB

 //pt > 4
 float centMB[1]={0.5};
 float centMB2S[1]={0.3};
 float centMBA[1]={0.4};
 float centMB3S[1]={0.7};
 float cent_limit_err[1]={0.1};
 // float raaMB1S[1]={0.564};
 // float raaMB1SstatErr[1]={0.077};
 // float raaMB1SsystErr[1]={0.071};
 // float raaMB2S[1]={0.120};
 // float raaMB2SstatErr[1]={0.038};
 // float raaMB2SsystErr[1]={0.020};
 // float raaMB3S[1]={0.033};
 // float raaMB3SstatErr[1]={0.035};//{0.035};
 // float raaMB3SsystErr[1]={0.006};

 float raaMB1S[1]={0.388};
 float raaMB1SstatErr[1]={0.018};
 float raaMB1SsystErr[1]={0.05};
 float raaMB2S[1]={0.106};
 float raaMB2SstatErr[1]={0.026};
 float raaMB2SsystErr[1]={0.021};
 float raaMB3S[1]={0.033};
 float raaMB3SstatErr[1]={0.035};//{0.035};
 float raaMB3SsystErr[1]={0.006};

 /// MB result!

 float raaMB1 [1]={0.388};
 float raaMB1e[1]={0.018};
 float raaMB1s[1]={0.042};

 float raaMB2 [1]={0.106};
 float raaMB2e[1]={0.026};
 float raaMB2s[1]={0.018};

 float raaMB3 [1]={0.046};
 float raaMB3e[1]={0.045};
 float raaMB3s[1]={0.028};

 float raaMB3Slimit[1]={0.1};
 float raaMBJpsi[1]={0.304};
 float raaMBJpsistatErr[1]={0.025};
 float raaMBJpsisystErr[1]={0.022};
 // float raaMBAlice[1]={0.439};
 // float raaMBAliceStat[1]={0.065};
 // float raaMBAliceSyst[1]={0.028};

 float raaMB1StotErr[0]=sqrt(pow(raaMB1e[0],2)+pow(raaMB1s[0],2));
 float raaMB2StotErr[0]=sqrt(pow(raaMB2e[0],2)+pow(raaMB2s[0],2));
 float raaMB3StotErr[0]=sqrt(pow(raaMB3e[0],2)+pow(raaMB3s[0],2));


 TF1 *f5 = new TF1("f4","1",0,1);
 f5->GetXaxis()->SetNdivisions(2);
 f5->SetLineWidth(1);
 f5->GetXaxis()->SetTitle("");
 f5->GetXaxis()->SetLabelColor(kWhite);
 f5->GetXaxis()->SetRangeUser(0,1);
 f5->GetYaxis()->SetTitle("R_{AA}");
 f5->GetYaxis()->SetTitleOffset(1.0);
 f5->GetYaxis()->SetRangeUser(0,1.5);
 f5->GetXaxis()->CenterTitle(kTRUE);
 f5->Draw();

 TGraphErrors *g1SMB = new TGraphErrors(1,centMB,raaMB1,centnoErr,raaMB1s);
 TGraphErrors *g1SMBtot = new TGraphErrors(1,centMB,raaMB1,centnoErr,raaMB1e);
 TGraphErrors *g1SMBcirc = new TGraphErrors(1,centMB,raaMB1,centnoErr,raaMB1s);

 g1SMBtot->SetMarkerStyle(33);
 g1SMBtot->SetMarkerColor(2);
 g1SMBtot->SetMarkerSize(2);
 g1SMBcirc->SetMarkerStyle(27);
 g1SMBcirc->SetMarkerColor(kBlack);
 g1SMBcirc->SetMarkerSize(2);
 g1SMB->SetLineColor(kRed-9);
 g1SMB->SetLineWidth(18);
 g1SMB->SetMarkerSize(0);
 g1SMB->Draw("e");
 g1SMBtot->Draw("pe");
 g1SMBcirc->Draw("p");

 float centMBErr[1]={0.15};
 TGraphErrors *g2SMB = new TGraphErrors(1,centMB2S,raaMB2,centnoErr,raaMB2s);
 TGraphErrors *g2SMBtot = new TGraphErrors(1,centMB2S,raaMB2,centnoErr,raaMB2e);
 TGraphErrors *g2SMBcirc = new TGraphErrors(1,centMB2S,raaMB2,centnoErr,raaMB2s);
 g2SMBtot->SetMarkerStyle(20);
 g2SMBtot->SetMarkerColor(kTeal+3);
 g2SMBtot->SetMarkerSize(1.5);
 g2SMBcirc->SetMarkerStyle(24);
 g2SMBcirc->SetMarkerColor(kBlack);
 g2SMBcirc->SetMarkerSize(1.5);
 g2SMB->SetLineColor(kTeal+5);
 g2SMB->SetFillStyle(0);
 g2SMB->SetLineWidth(18);
 g2SMB->SetMarkerSize(0);
 g2SMB->Draw("2");
 g2SMBtot->Draw("pe");
 g2SMBcirc->Draw("p");
 //ALICE 0-90% shit
 // TGraphErrors *gAMB = new TGraphErrors(1,centMBA,raaMBAlice,centnoErr,raaMBAliceSyst);
 // TGraphErrors *gAMBtot = new TGraphErrors(1,centMBA,raaMBAlice,centnoErr,raaMBAliceStat);
 // gAMBtot->SetMarkerStyle(33);
 // gAMBtot->SetMarkerColor(kBlack);
 // gAMBtot->SetMarkerSize(2);
 // gAMB->SetLineColor(kGray);
 // gAMB->SetFillStyle(0);
 // gAMB->SetLineWidth(18);
 // gAMB->SetMarkerSize(0);
 // gAMB->Draw("2");
 // gAMBtot->Draw("pe");

 /* TGraphErrors *gJpsiMB = new TGraphErrors(1,centMB,raaMBJpsi,centMBErr,raaMBJpsisystErr); */
 /* TGraphErrors *gJpsiMBtot = new TGraphErrors(1,centMB,raaMBJpsi,centnoErr,raaMBJpsistatErr); */
 /* gJpsiMBtot->SetMarkerStyle(33); */
 /* gJpsiMBtot->SetMarkerColor(kBlue+2); */
 /* gJpsiMBtot->SetMarkerSize(1.5); */
 /* gJpsiMB->SetFillStyle(3244); */
 /* gJpsiMB->SetFillColor(kBlue-7); */
 /* gJpsiMB->SetLineWidth(18); */
 /* gJpsiMB->SetMarkerSize(0); */
 /* gJpsiMB->Draw("2"); */

 /* TGraphErrors *gJpsiMB_ = new TGraphErrors(1,centMB,raaMBJpsi,centMBErr,raaMBJpsisystErr); */
 /* gJpsiMB_->SetFillStyle(0); */
 /* gJpsiMB_->SetLineColor(kBlue+1); */
 /* gJpsiMB_->SetLineWidth(1); */
 /* gJpsiMB_->SetMarkerSize(0); */
 /* gJpsiMB_->Draw("2"); */

 /* gJpsiMBtot->Draw("pe"); */


 // TGraphErrors *g3SMB = new TGraphErrors(1,centMB3S,raaMB3,centnoErr,raaMB3s);
 // TGraphErrors *g3SMBtot = new TGraphErrors(1,centMB3S,raaMB3,centnoErr,raaMB3e);
 // g3SMBtot->SetMarkerStyle(21);
 // g3SMBtot->SetMarkerColor(kBlue+2);
 // g3SMBtot->SetMarkerSize(1.3);
 // TGraphErrors *glimits = new TGraphErrors(1,centMB3S,raaMB3Slimit,cent_limit_err,centnoErr);
 // glimits->SetMarkerSize(0);
 // glimits->SetLineWidth(2);
 //glimits->Draw("D");
 // TArrow *g3SMBlimit = new TArrow(centMB3S[0],raaMB3Slimit[0],centMB3S[0],0,0.015,"|>");
 // g3SMBlimit->SetLineWidth(2);
 // g3SMBlimit->Draw();
 // g3SMB->Draw("e");
 //g3SMBtot->Draw("pe");

 TGraphErrors *g3SMB = new TGraphErrors(1,centMB,raaMB3,centnoErr,raaMB3s);
 TGraphErrors *g3SMBtot = new TGraphErrors(1,centMB,raaMB3,centnoErr,raaMB3e);

 g3SMBtot->SetMarkerStyle(21);
 g3SMBtot->SetMarkerColor(kBlue+2);
 g3SMBtot->SetMarkerSize(1.3);
 g3SMB->SetLineColor(kBlue-9);
 g3SMB->SetLineWidth(18);
 g3SMB->SetMarkerSize(0);
 g3SMB->Draw("e");
 g3SMBtot->Draw("pe");
 legend->AddEntry(g3SMBtot, "#varUpsilon(3S), 0-100% centrality","p");
 TLatex *MBcentrality = new TLatex(0.05,0.69,"0-100%");
 MBcentrality->SetTextSize(0.23);
 MBcentrality->Draw();
 TLatex *MBAcentrality = new TLatex(0.05,0.33,"0-90%");
 MBAcentrality->SetTextSize(0.23);
 // MBAcentrality->Draw();

 p2->Update();
 p1->Update();
 //	c1->SaveAs("Raa_JpsiY1SY2S_withMB.gif");
 c1->SaveAs("Raa_RegitY1SY2S_MB_v0.pdf");

 //// second canvas : raa vs Pt!
  TCanvas *cpt = new TCanvas("cpt","cpt");

 cpt->cd();
 TPad *ppt = new TPad("ppt","ppt",0.0,0.0,0.9,1.0);
 ppt->SetBottomMargin(0.12);
 ppt->SetTopMargin(0.03);
 ppt->SetRightMargin(0.03);
 ppt->SetLeftMargin(0.11);
 ppt->Draw();
 ppt->cd();

 TF1 *f4Pt = new TF1("f4Pt","1",0,22);
 f4Pt->SetLineWidth(1);
 f4Pt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} (GeV/c)");
 f4Pt->GetYaxis()->SetTitle("R_{AA}");
 f4Pt->GetYaxis()->SetTitleOffset(1.0);
 f4Pt->GetYaxis()->SetRangeUser(0,1.5);
 f4Pt->GetXaxis()->CenterTitle(kTRUE);
 f4Pt->Draw();

 TGraphErrors *gpt1 = new TGraphErrors(binPt,pt,raaPt1,pte,raaPt1e);
 gpt1->SetMarkerColor(2);
 gpt1->SetMarkerStyle(21);
 gpt1->SetMarkerSize(1.2);

 TGraphErrors *gpt1s = new TGraphErrors(binPt,pt,raaPt1,centnoErr,raaPt1s);
 gpt1s->SetLineColor(kRed-9);
 gpt1s->SetLineWidth(18);
 gpt1s->SetMarkerSize(0);
 gpt1s->Draw("e");
 
 gpt1->Draw("pe");

 TGraphErrors *gpt1circle = new TGraphErrors(binPt,pt,raaPt1,pte,raaPt1e);
 gpt1circle->SetMarkerStyle(25);
 gpt1circle->SetMarkerSize(1.2);
 gpt1circle->SetLineColor(kBlack);
 gpt1circle->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();

 // 2S pt dependence ??!!

 TGraphErrors *gpt2 = new TGraphErrors(binPt,pt,raaPt2,pte,raaPt2e);
 gpt2->SetMarkerColor(8);
 gpt2->SetMarkerStyle(21);
 gpt2->SetMarkerSize(1.2);

 TGraphErrors *gpt2s = new TGraphErrors(binPt,pt,raaPt2,centnoErr,raaPt2s);
 gpt2s->SetLineColor(kTeal+5);
 gpt2s->SetLineWidth(18);
 gpt2s->SetMarkerSize(0);
 gpt2s->Draw("e");
 
 gpt2->Draw("pe");

 TGraphErrors *gpt2circle = new TGraphErrors(binPt,pt,raaPt2,pte,raaPt2e);
 gpt2circle->SetMarkerStyle(25);
 gpt2circle->SetMarkerSize(1.2);
 gpt2circle->SetLineColor(kBlack);
 gpt2circle->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();

 TLatex *l1CMSpt = new TLatex(1,1.39,"#sqrt{s_{NN}} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.032);
 l1CMSpt->Draw();
 TLatex *lyLPt= new TLatex(1,1.3,"L_{int}^{PbPb} = 150 #mub^{-1}; |y| < 2.4");
 lyLPt->SetTextSize(0.029);
 lyLPt->Draw();
 ppt->Update();
 TLatex *lyLYPt2= new TLatex(1,1.21,"L_{int}^{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyLYPt2->SetTextSize(0.029);
 lyLYPt2->Draw();
 ppt->Update();
 TLegend *legendpt = new TLegend(0.182,0.8,0.5,0.7);
 legendpt->SetTextSize(0.029);
 legendpt->SetFillStyle(0);
 legendpt->SetFillColor(0);
 legendpt->SetBorderSize(0);
 legendpt->SetTextFont(42);
 legendpt->AddEntry(gpt1,"#varUpsilon(1S)","p");
 legendpt->AddEntry(gpt2,"#varUpsilon(2S)","p");
 legendpt->Draw();

 cpt->SaveAs("Raa_RegitY1S2S_pt.pdf");

 //// third canvas : raa vs rap!
  TCanvas *cy = new TCanvas("cy","cy");

 cy->cd();
 TPad *py = new TPad("py","py",0.0,0.0,0.9,1.0);
 py->SetBottomMargin(0.12);
 py->SetTopMargin(0.03);
 py->SetRightMargin(0.03);
 py->SetLeftMargin(0.11);
 py->Draw();
 py->cd();

 TF1 *f4Y = new TF1("f4Y","1",0,2.45);
 f4Y->SetLineWidth(1);
 f4Y->GetXaxis()->SetTitle("|y^{#Upsilon_{cand.}}|");
 f4Y->GetYaxis()->SetTitle("R_{AA}");
 f4Y->GetYaxis()->SetTitleOffset(1.0);
 f4Y->GetYaxis()->SetRangeUser(0,1.5);
 f4Y->GetXaxis()->CenterTitle(kTRUE);
 f4Y->Draw();

 TGraphErrors *gy1 = new TGraphErrors(binRap,rap,raaRap1,rape,raaRap1e);
 gy1->SetMarkerColor(2);
 gy1->SetMarkerStyle(21);
 gy1->SetMarkerSize(1.2);

 TGraphErrors *g1Ssyst = new TGraphErrors(binRap,rap,raaRap1,centnoErr,raaRap1s);
 g1Ssyst->SetLineColor(kRed-9);
 g1Ssyst->SetLineWidth(18);
 g1Ssyst->SetMarkerSize(0);
 g1Ssyst->Draw("e");

 gy1->Draw("pe");

 TGraphErrors *gy1circle = new TGraphErrors(bin,rap,raaRap1,rape,raaRap1e);
 gy1circle->SetMarkerStyle(25);
 gy1circle->SetMarkerSize(1.2);
 gy1circle->SetLineColor(kBlack);
 gy1circle->Draw("p");
 f4Y->Draw("same");
 gPad->RedrawAxis();

 TGraphErrors *gy2 = new TGraphErrors(binRap,rap,raaRap2,rape,raaRap2e);
 gy2->SetMarkerColor(kTeal+3);
 gy2->SetMarkerStyle(21);
 gy2->SetMarkerSize(1.2);

 TGraphErrors *g2Ssyst = new TGraphErrors(binRap,rap,raaRap2,centnoErr,raaRap2s);
 g2Ssyst->SetLineColor(kTeal+5);
 g2Ssyst->SetLineWidth(18);
 g2Ssyst->SetMarkerSize(0);
 g2Ssyst->Draw("e");

 gy2->Draw("pe");

 TGraphErrors *gy2circle = new TGraphErrors(bin,rap,raaRap2,rape,raaRap2e);
 gy2circle->SetMarkerStyle(25);
 gy2circle->SetMarkerSize(1.2);
 gy2circle->SetLineColor(kBlack);
 gy2circle->Draw("p");
 f4Y->Draw("same");
 gPad->RedrawAxis();

 TLatex *l1CMSy = new TLatex(0.1,1.39, "#sqrt{s_{NN}} = 2.76 TeV");
 l1CMSy->SetTextFont(42);
 l1CMSy->SetTextSize(0.032);
 l1CMSy->Draw();
 TLatex *lyLY= new TLatex(0.1,1.3,"L_{int}^{PbPb} = 150 #mub^{-1}; |y| < 2.4");
 lyLY->SetTextSize(0.029);
 lyLY->Draw();
 py->Update();
 TLatex *lyLY2= new TLatex(0.1,1.21,"L_{int}^{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyLY2->SetTextSize(0.029);
 lyLY2->Draw();
 py->Update();
 TLegend *legendy = new TLegend(0.182,0.8,0.5,0.7);
 legendy->SetTextSize(0.029);
 legendy->SetFillStyle(0);
 legendy->SetFillColor(0);
 legendy->SetBorderSize(0);
 legendy->SetTextFont(42);
 legendy->AddEntry(gy1,"#varUpsilon(1S)","p");
 legendy->AddEntry(gy2,"#varUpsilon(2S)","p");
 legendy->Draw();
 cy->SaveAs("Raa_RegitY1S2S_rap.pdf");


 // doubledifferential Npart/rapidity
 // !!!! first entry is |y|<1, second entry is 1<|y|<2.4 !!!!
 float raaRN1_c [2] = {0.364,0.407}; // central values
 float raaRN1e_c[2] = {0.036,0.047};
 float raaRN1s_c[2] = {0.033,0.029};

 float raaRN2_c [2] = {0.073,0.127};
 float raaRN2e_c[2] = {0.048,0.061};
 float raaRN2s_c[2] = {0.056,0.058};

 float raaRN3_c [2] = {-0.04,0.108};
 float raaRN3e_c[2] = {0.   ,0.105};
 float raaRN3s_c[2] = {0.   ,0.035};

 float raaRN1_p [2] = {0.423,0.596}; // non-central (20-100) values
 float raaRN1e_p[2] = {0.035,0.047};
 float raaRN1s_p[2] = {0.102,0.033};

 float raaRN2_p [2] = {0.168,0.296};
 float raaRN2e_p[2] = {0.057,0.076};
 float raaRN2s_p[2] = {0.160,0.027};

 float raaRN3_p [2] = {0.145,0.152};
 float raaRN3e_p[2] = {0.107,0.129};
 float raaRN3s_p[2] = {0.064,0.024};

 float rap_short_c[2] = {0.45,1.65};
 float rap_short_p[2] = {0.55,1.75};
 float rap_short_3[2] = {0.58,1.78};
 float rap_shorte[2]= {0.55,0.65};
 
/// fourth canvas : raa vs rap for two nPart bins!
 TCanvas *cyn = new TCanvas("cyn","cyn");

 cyn->cd();
 TPad *pyn = new TPad("pyn","pyn",0.0,0.0,0.9,1.0);
 pyn->SetBottomMargin(0.12);
 pyn->SetTopMargin(0.03);
 pyn->SetRightMargin(0.03);
 pyn->SetLeftMargin(0.11);
 pyn->Draw();
 pyn->cd();

 TF1 *f4Yn = new TF1("f4Y","1",0,2.45);
 f4Yn->SetLineWidth(1);
 f4Yn->GetXaxis()->SetTitle("|y^{#Upsilon_{cand.}}|");
 f4Yn->GetYaxis()->SetTitle("R_{AA}");
 f4Yn->GetYaxis()->SetTitleOffset(1.0);
 f4Yn->GetYaxis()->SetRangeUser(0,1.5);
 f4Yn->GetXaxis()->CenterTitle(kTRUE);
 f4Yn->Draw();

 TGraphErrors *gyc1 = new TGraphErrors(2,rap_short_c,raaRN1_c,rap_shorte,raaRN1e_c);
 gyc1->SetMarkerColor(2);
 gyc1->SetMarkerStyle(21);
 gyc1->SetMarkerSize(1.2);

 TGraphErrors *gc1Ssyst = new TGraphErrors(2,rap_short_c,raaRN1_c,centnoErr,raaRN1s_c);
 gc1Ssyst->SetLineColor(kRed-9);
 gc1Ssyst->SetLineWidth(18);
 gc1Ssyst->SetMarkerSize(0);
 gc1Ssyst->Draw("e");

 gyc1->Draw("pe");

 TGraphErrors *gcy1circle = new TGraphErrors(2,rap_short_c,raaRN1_c,rap_shorte,raaRN1e_c);
 gcy1circle->SetMarkerStyle(25);
 gcy1circle->SetMarkerSize(1.2);
 gcy1circle->SetLineColor(kBlack);
 gcy1circle->Draw("p");
 f4Yn->Draw("same");
 gPad->RedrawAxis();

 TGraphErrors *gyp1 = new TGraphErrors(2,rap_short_p,raaRN1_p,centnoErr,raaRN1e_p);
 gyp1->SetMarkerColor(kPink+8);
 gyp1->SetMarkerStyle(21);
 gyp1->SetMarkerSize(1.2);

 TGraphErrors *gp1Ssyst = new TGraphErrors(2,rap_short_p,raaRN1_p,centnoErr,raaRN1s_p);
 gp1Ssyst->SetLineColor(kPink+1);
 gp1Ssyst->SetLineWidth(18);
 gp1Ssyst->SetMarkerSize(0);
 gp1Ssyst->Draw("e");

 gyp1->Draw("pe");

 TGraphErrors *gpy1circle = new TGraphErrors(2,rap_short_p,raaRN1_p,centnoErr,raaRN1e_p);
 gpy1circle->SetMarkerStyle(25);
 gpy1circle->SetMarkerSize(1.2);
 gpy1circle->SetLineColor(kBlack);
 gpy1circle->Draw("p");
 f4Yn->Draw("same");
 gPad->RedrawAxis();


 TGraphErrors *gyc2 = new TGraphErrors(2,rap_short_c,raaRN2_c,rap_shorte,raaRN2e_c);
 gyc2->SetMarkerColor(kTeal+3);
 gyc2->SetMarkerStyle(20);
 gyc2->SetMarkerSize(1.2);

 TGraphErrors *gc2Ssyst = new TGraphErrors(2,rap_short_c,raaRN2_c,centnoErr,raaRN2s_c);
 gc2Ssyst->SetLineColor(kTeal+5);
 gc2Ssyst->SetLineWidth(18);
 gc2Ssyst->SetMarkerSize(0);
 gc2Ssyst->Draw("e");

 gyc2->Draw("pe");

 TGraphErrors *gcy2circle = new TGraphErrors(2,rap_short_c,raaRN2_c,rap_shorte,raaRN2e_c);
 gcy2circle->SetMarkerStyle(24);
 gcy2circle->SetMarkerSize(1.2);
 gcy2circle->SetLineColor(kBlack);
 gcy2circle->Draw("p");
 f4Yn->Draw("same");
 gPad->RedrawAxis();

 TGraphErrors *gyp2 = new TGraphErrors(2,rap_short_p,raaRN2_p,centnoErr,raaRN2e_p);
 gyp2->SetMarkerColor(kOrange+2);
 gyp2->SetMarkerStyle(20);
 gyp2->SetMarkerSize(1.2);

 TGraphErrors *gp2Ssyst = new TGraphErrors(2,rap_short_p,raaRN2_p,centnoErr,raaRN2s_p);
 gp2Ssyst->SetLineColor(kOrange-9);
 gp2Ssyst->SetLineWidth(18);
 gp2Ssyst->SetMarkerSize(0);
 gp2Ssyst->Draw("e");

 gyp2->Draw("pe");

 TGraphErrors *gpy2circle = new TGraphErrors(2,rap_short_p,raaRN2_p,centnoErr,raaRN2e_p);
 gpy2circle->SetMarkerStyle(24);
 gpy2circle->SetMarkerSize(1.2);
 gpy2circle->SetLineColor(kBlack);
 gpy2circle->Draw("p");
 f4Yn->Draw("same");
 gPad->RedrawAxis();

TGraphErrors *gyp3 = new TGraphErrors(2,rap_short_3,raaRN3_p,centnoErr,raaRN3e_p);
 gyp3->SetMarkerColor(kBlue+2);
 gyp3->SetMarkerStyle(33);
 gyp3->SetMarkerSize(1.5);

 TGraphErrors *gp3Ssyst = new TGraphErrors(2,rap_short_3,raaRN3_p,centnoErr,raaRN3s_p);
 gp3Ssyst->SetLineColor(kBlue-9);
 gp3Ssyst->SetLineWidth(18);
 gp3Ssyst->SetMarkerSize(0);
 //gp3Ssyst->Draw("e");

 //gyp3->Draw("pe");

 TGraphErrors *gpy3circle = new TGraphErrors(2,rap_short_3,raaRN3_p,centnoErr,raaRN3e_p);
 gpy3circle->SetMarkerStyle(27);
 gpy3circle->SetMarkerSize(1.5);
 gpy3circle->SetLineColor(kBlack);
 // gpy3circle->Draw("p");
 f4Yn->Draw("same");
 gPad->RedrawAxis();

TGraphErrors *gyc3 = new TGraphErrors(2,rap_short_c,raaRN3_c,centnoErr,raaRN3e_c);
 gyc3->SetMarkerColor(kViolet-3);
 gyc3->SetMarkerStyle(33);
 gyc3->SetMarkerSize(1.5);

 TGraphErrors *gc3Ssyst = new TGraphErrors(2,rap_short_c,raaRN3_c,centnoErr,raaRN3s_c);
 gc3Ssyst->SetLineColor(kViolet+7);
 gc3Ssyst->SetLineWidth(18);
 gc3Ssyst->SetMarkerSize(0);
 // gc3Ssyst->Draw("e");

 //gyc3->Draw("pe");

 TGraphErrors *gcy3circle = new TGraphErrors(2,rap_short_c,raaRN3_c,centnoErr,raaRN3e_c);
 gcy3circle->SetMarkerStyle(27);
 gcy3circle->SetMarkerSize(1.5);
 gcy3circle->SetLineColor(kBlack);
 // gcy3circle->Draw("p");
 f4Yn->Draw("same");
 gPad->RedrawAxis();

 l1CMSy->Draw();
 lyLY->Draw();
 pyn->Update();
 lyLY2->Draw();
 pyn->Update();

 TLegend *legendyn = new TLegend(0.5,0.7,0.9,0.9);
 legendyn->SetTextSize(0.029);
 legendyn->SetFillStyle(0);
 legendyn->SetFillColor(0);
 legendyn->SetBorderSize(0);
 legendyn->SetTextFont(42);
 legendyn->AddEntry(gyc1,"#varUpsilon(1S), 0-20%","p");
 legendyn->AddEntry(gyp1,"#varUpsilon(1S), 20-100%","p");
 legendyn->AddEntry(gyc2,"#varUpsilon(2S), 0-20%","p");
 legendyn->AddEntry(gyp2,"#varUpsilon(2S), 20-100%","p");
 // legendyn->AddEntry(gyc3,"#varUpsilon(3S), 0-20%","p");
 //legendyn->AddEntry(gyp3,"#varUpsilon(3S), 20-100%","p");
 legendyn->Draw();
 cyn->SaveAs("Raa_RegitY1S2S_rapNpart.pdf");






}
