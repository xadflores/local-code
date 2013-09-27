#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
//---misc---
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
//---Roofit includes
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooKeysPdf.h"
#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
//---ROOT includes---
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "RooHist.h"
#include "TPaveText.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TF1.h"
/*INPUTs*/
//default mass range for frame/fitting

double mass_l = 7.0;
double mass_h = 14.0;
double binw_ = 0.1;    //bin width of the histogram

const int nData = 7 ;
const char* sampleLegend[nData] = {"",
				   "pp #sqrt{s} = 7TeV",
				   "pp 2.76 TeV",
				   "PbPb Regit 2.76 TeV",
				   "PbPb ppReco TrkTrk",
				   "PbPb ppReco GlbGlb",
				   "PbPb hiReco GlbGlb"
};

const char* lumiLegend[nData]= {"",		       
				"L_{int} = ?nb^{-1}",
				"L_{int} = ?pb^{-1}",
				"L_{int} = 150 \mub^{-1}",
				"L_{int} = ? \mub^{-1}",
				"L_{int} = ? \mub^{-1}",
				"L_{int} = 150 \mub^{-1}",
};


using namespace std;
using namespace RooFit;
using namespace RooStats;

// incrementing 1 for each sample ; 0 pp276, 1 pp7 ; 2 PbPb; 3 pPb.

void fitpeaks(int bin, int area);
void rapbin(int bin, RooRealVar*,RooRealVar *);
void sidebandFit(bool, RooRealVar*, RooRealVar*, RooRealVar*, RooAbsPdf*, RooAbsPdf*, RooDataSet*, RooPlot*);
void suppression(RooRealVar*, RooFitResult*, RooAbsPdf*, RooAbsPdf*, int, RooDataSet*, RooPlot*);
void setUpLimits(float xMin, float xMax, RooAbsPdf*, RooDataSet*, RooRealVar*, float baseNll);
pair<double, double> ConfidencInterval(float, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf);

void fitUpsilonYields(){
  gROOT->Macro("../code/cm/logon.C+");
	for(int it=0;it<nBins;it++){
	  yAxisRap2S[it]=0.;
	  yAxisRapError2S[it]=0.;
	  yAxisRap2S_276[it]=0.;
	  yAxisRapError2S_276[it]=0.;
	  yAxisRap1S[it]=0.;
	  yAxisRapError1S[it]=0.;
	  yAxisRap3S[it]=0.;
	  yAxisRapError3S[it]=0.;
	  yAxisRap3S_276[it]=0.;
	  yAxisRapError3S_276[it]=0.;
	}   
	
	//different bins possible:
	for(int ah=3;ah<7;ah++){

	  		  if (MC) {
		    area=0;
		    /* finput = "../dimuonTree_MC_5TeV_completo.root";
		    figName_ = "masspeak_MC_1S2S"; 
		    if(ah==6)figs_="../output/1102/MC/";
		    if (TRKROT) figName_+="TrkRot_";
		    if (LS_constrain!=0) figName_+="LSconstr_";
		    figName_+=ah;
		    //bkgdModel=5;
		    bkgdModel=ah;
		    MC12=1;
		    fitpeaks(6,area);
		    MC12=0;*/
		  figs_="../output/2302/MC/";
		    finput = "../dimuonTree_MC_1S_5TeV_big.root";
		    figName_ = "masspeak_MC_1S_TrkRotNoLS"; 
		    figName_+=ah;
		    bkgdModel=ah;
		    MC1=1;
		    fitpeaks(6,area);
		    // MC1=0;
		    // finput = "../dimuonTree_MC_2S_5TeV_big.root";
		    // figName_ = "masspeak_MC_2S"; 
		    // MC2=1;
		    // fitpeaks(6,area);
		    // MC2=0;
		  }
	  if(pp276){
	    area=0;
	    finput = "../dimuonTree_upsiMiniTree_pp276tev_5p41_Run211739-211831_trigBit1_allTriggers0_pt4.root";
	   // finput = "../dimuonTree_upsiMiniTree_pp276tev_ppreco__Runa_trigBit1_allTriggers0_pt4.root";
	    figName_ = "masspeak_pp276_PULL";
	    figs_="../output/2302/pp276/";
	    if (TRKROT) figName_+="TrkRot_";
	    if (LS_constrain!=0) figName_+="LSconstr_";	  
	    figName_+=ah;
	    bkgdModel=ah;
	    //  bkgdModel=6;
	    fitpeaks(6,area);
	  }
	 if (pp7) {
	   area=1;
	   binw_=0.05;
	   
	   finput = "../dimuonTree_upsiMiniTree_pp7tev_glbglb_Runa_trigBit1_allTriggers0_pt4.root";
	   // dm3GlbGlb
	   mass_l = 8.5;
	   mass_h = 11.5;
	   figs_= "../output/0403/pp7d3/";
	   figName_ = "masspeak_pp7_ratiosd3_logo_LS";
	   TRKROT=1;
	   if (TRKROT) figName_+="TrkRot_";
	   if (LS_constrain!=0) figName_+="LSconstr_";
	   //bkgdModel=4;
	   bkgdModel=ah;
	   figName_+=ah;
	    fitpeaks(6,area);
	   
	   //dm0
	   finput = "../dimuonTree_upsiMiniTree_pp7tev_dimu0v1_Runa_trigBit1_allTriggers0_pt4.root";//trk-trk
	   mass_l = 8.5;
	   mass_h = 11.5;
	   figs_= "../output/0403/pp7d0/";
	   figName_ = "masspeak_pp7_ratiosd0_as_LS";		
	   TRKROT=1;
	   if (TRKROT) figName_+="TrkRot_";
	   if (LS_constrain!=0) figName_+="LSconstr_";
	   //bkgdModel=4;
	   bkgdModel=ah;
	   figName_+=ah;
	    fitpeaks(6,area);
	   //finput = "../dimuonTree_upsiMiniTree_pp7tev__Runa_trigBit1_allTriggers0_pt4.root";
	   
	   // dm3 trk-trk
	    mass_l = 8.5;
	   mass_h = 11.5;
	   figs_= "../output/0403/pp7d3/";
	   figName_ = "masspeak_pp7_ratiosd3trk_NominalA_Good";
	   TRKROT=1;
	   if (TRKROT) figName_+="TrkRot_";
	   if (LS_constrain!=0) figName_+="LSconstr_";
	   //bkgdModel=4;
	   bkgdModel=ah;
	   figName_+=ah;
	   // fitpeaks(6,area);
	 }
	 if (PbPb){ 
	   area=2;
	   mass_l = 7;
	   mass_h = 14;
	   binw_=0.14;
	   finput="../dimuonTree_upsiMiniTree_aa276tev_regitreco_glbglb_Runa_trigBit1_allTriggers0_pt4.root";
	   //   finput="../dimuonTree_upsiMiniTree_aa276tev_50100ppreco_glbglb_Runa_trigBit1_allTriggers0_pt4.root";
	   // finput = "../dimuonTree_upsiMiniTree_aa276tev_50100ppreco__Runa_trigBit1_allTriggers0_pt4.root"; 
	   figName_ = "masspeak_PbPbRegit_OurnominalANcuts_";
	   figs_="../output/0403/PbPbRegit/";
	   TRKROT=0;
	   if (TRKROT) figName_+="TrkRot_";
	   if (LS_constrain!=0) figName_+="LSconstr_";
	   figName_+=ah;
	   //bkgdModel=3; 
	   bkgdModel=ah; 
	   // fitpeaks(7,area);
	   // fitpeaks(8,area);
	   // fitpeaks(9,area);	
	   // fitpeaks(10,area);	
	   // fitpeaks(11,area);	
	   // fitpeaks(12,area);	
	   // fitpeaks(13,area);	
	   // fitpeaks(14,area);
	   
	   fitpeaks(15,area);	 
	      
	   // //   ----reprocessing HI2011 tree-------
	   //  finput="../dimuonTree_HI2011_fulldataset_trkRot.root";
	   // // figName_ = "masspeak_PbPbHI2011_";
	   // // figs_="../output/2602/2011/";
	   // // figName_+=ah;
	   // fitpeaks(7,area);
	   // fitpeaks(8,area);

	   //
	  

// GlbGlb pairs (do we recover the 2011 stats ?)
	   //----"zhen" tree------------------
	   	   finput="../dimuonTree_HI2011_fulldataset_trkRot.root";
	   figName_ = "masspeak_PbPb2011_OurnominalANcuts_";	  
	   figs_="../output/2802/2011/";	  
	   figName_+=ah;    			  
  	  // fitpeaks(9,area);	
  	  // fitpeaks(10,area);	
	  // fitpeaks(11,area);	
	  // fitpeaks(12,area);	
	  // fitpeaks(13,area);
  
	  // fitpeaks(14,area);	
	  fitpeaks(15,area);	 
	      

	   ///////centrality bins on regit sample : change the rapbins for each fit.
	   // 0-5% :     0-2
	   // 5-10 %:    2-4
	   // 10-20 %:   4-8
	   // 20-30 %:   8-12
	   // 30-40 %:   12-16
	   // 40-50 %:   16-20
	   // 50-100%:   20-40
	   
	 }
	 if (pA) {
	   area=3;
	   mass_l = 7;
	   mass_h = 14;
	   finput = "../dimuonTree_upsiMiniTree_pA5tev_18p4nb__Run210498-211256_trigBit1_allTriggers0_pt4.root";
	   figName_ = "masspeak_pPb_MCfixed"; 
	    figs_="../output/2302/pPb/";
	   if (TRKROT) figName_+="TrkRot_";
	   if (LS_constrain!=0) figName_+="LSconstr_";
	   figName_+=ah;
	   //bkgdModel=5;
	   bkgdModel=ah;
	   fitpeaks(6,area);
	 }
	 else if (!finput){
	   cout << "I don't have any file to eat !" << endl ;
	 }
	 cout << "bkg Model #" << bkgdModel << endl ;
		}/*
		 // for overlaying purposes :
		 finput = "../dimuonTree_upsiMiniTree_pp276tev_ppreco__Runa_trigBit1_allTriggers0_pt4.root";
		 figs_+="pp276/";
		 figName_ = "masspeak_pp276_"; 
		 if (TRKROT) figName_+="TrkRot_";
		 if (LS_constrain!=0) figName_+="LSconstr_";
		 fitpeaks(6);
		 mass_l=8.5;
		 mass_h=11.5;
		 finput = "../dimuonTree_upsiMiniTree_pp7tev__Runa_trigBit1_allTriggers0_pt4.root";
		 figs_="../output/0702/pp7/";
		 figName_ = "masspeak_pp7_"; 
		 TRKROT = 1;
		 if (TRKROT) figName_+="TrkRot_";
		 if (LS_constrain!=0) figName_+="LSconstr_";
		 fitpeaks(7);*/
		
}

void fitpeaks(int bin,int area){

  if(area==1) mass_l = 8.5, mass_h = 11.5;
  double mmin_ =  mass_l; 
  double mmax_ = mass_h;
  double fmin_ = mass_l;
  double fmax_ = mass_h;
  if(area==2){
    switch(bin){
    case 0 :
      cut_="";
      suffix_="_NegRap193_047";
      break;
    case 1 :
      cut_="upsRapidity>-0.94&&upsRapidity<-0.47";
      suffix_="_NegRap047_000";
      break;
    case 2 :
      cut_="upsRapidity>-0.47&&upsRapidity<0";
      suffix_="_PosRap000_047";
      break;
    case 3 :
      cut_="upsRapidity>0&&upsRapidity<1.46";
      suffix_="_PosRap047_193";
      break;
    case 4 :
      cut_="upsRapidity>-2.4&&upsRapidity<-0.47";
      suffix_="_Neg";
      break;
    case 5 :
      cut_="upsRapidity>-0.47&&upsRapidity<1.46";
      suffix_="_Pos";
      break;
    case 6 :
      cut_="upsRapidity>-2.4&&upsRapidity<1.46";
      suffix_="_All";
      break;
	   ///////centrality bins on regit sample : change the rapbins for each fit.
	   // 0-5% :     0-2
	   // 5-10 %:    2-4
	   // 10-20 %:   4-8
	   // 20-30 %:   8-12
	   // 30-40 %:   12-16
	   // 40-50 %:   16-20
	   // 50-100%:   20-40

      case 7: // now with vProb005 -------------------- //&&(abs(upsRapidity)<1) removed
      cut_="(Centrality>=0)&&(Centrality<2)";
      suffix_="_cntr0-5";
      break;
    case 8:
      cut_="(Centrality>=20)&&(Centrality<40)&&vProb>0.01";
      suffix_="_cntr50-100";
      break;
    case 9:
      cut_="Centrality>=0 && Centrality<2";
      suffix_="_cntr0-5";
      break;
    case 10:
      cut_="Centrality>=2 && Centrality<4";
      suffix_="_cntr5-10";
      break;
    case 11:
      cut_="Centrality>=4 && Centrality<8";
      suffix_="_cntr10-20";
      break;
    case 12:
      cut_="Centrality>=8 && Centrality<12";
      suffix_="_cntr20-30"; 
      break;
    case 13:
      cut_="Centrality>=12 && Centrality<16";
      suffix_="_cntr30-40"; 
      break;
    case 14:
      cut_="Centrality>=16 && Centrality<20";
      suffix_="_cntr40-50"; 
      break;
    case 15:
      cut_="Centrality>=20 && Centrality<50";
      suffix_="_cntr50-100"; 
      break;
    case 16:
      cut_="Centrality>=20 && Centrality<24";
      suffix_="_cntr50-60"; 
      break;
    default: break;
    }
  }
  else {
    switch(bin){
    case 0 :
      cut_="upsRapidity>-1.93&&upsRapidity<-0.47";
      suffix_="_NegRap193_047";
      break;
    case 1 :
      cut_="upsRapidity>-0.47&&upsRapidity<0.";
      suffix_="_NegRap047_000";
      break;
    case 2 :
      cut_="upsRapidity>0.&&upsRapidity<0.47";
      suffix_="_PosRap000_047";
      break;
    case 3 :
      cut_="upsRapidity>0.47&&upsRapidity<1.93";
      suffix_="_PosRap047_193";
      break;
    case 4 :
      cut_="upsRapidity>-1.93&&upsRapidity<0.";
      suffix_="_Neg";
      break;
    case 5 :
      cut_="upsRapidity>-0.&&upsRapidity<1.93";
      suffix_="_Pos";
      break;
    case 6 :
      cut_="upsRapidity>-1&&upsRapidity<1";
      suffix_="_Barrel";
      break;
    case 7 :
      cut_="upsRapidity>-1&&upsRapidity<1"; // deprecated
      suffix_="_OverlayPP";
    default: break;
    }
  }
  
  cout << "oniafitter processing"
       << "\n\tInput:  \t" << finput
       << "\n\tresults:\t" << figs_
       << endl;
  ofstream outfile(figs_+"fitresults-ppBarrel.out", ios_base::app);
  
  
  //read the data
  TFile f(finput,"read");
  gDirectory->Cd(finput+":/"+dirname_);
  TTree* theTree     = (TTree*)gROOT->FindObject("UpsilonTree");
  TTree* allsignTree     = (TTree*)gROOT->FindObject("UpsilonTree_allsign");
  if (PR_plot) {TRKROT = 1; PbPb=1;}
  if (TRKROT) TTree* trkRotTree = (TTree*)gROOT->FindObject("UpsilonTree_trkRot");
  cout << mmin_ << " " << mmax_ << " "<< bin << " 'bin' mmin et mmax"<<  endl ;
cout << mass_l << " " << mass_h << " " << bin << " 'bin' mass_l et mass_h" << endl ;
  RooRealVar* mass  = new RooRealVar("invariantMass","#mu#mu mass",mmin_,mmax_,"GeV/c^{2}");
  RooRealVar* upsPt  = new RooRealVar("upsPt","p_{T}(#Upsilon)",0,60,"GeV");
  RooRealVar* upsEta = new RooRealVar("upsEta",  "upsEta"  ,-7,7);
  RooRealVar* vProb = new RooRealVar("vProb",  "vProb"  ,0.01,1.00);
  RooRealVar* QQsign = new RooRealVar("QQsign",  "QQsign"  ,-1,5);
  RooRealVar* weight = new RooRealVar("weight",  "weight"  ,-2,2);
  if (area==2) RooRealVar* Centrality = new RooRealVar("Centrality","Centrality",0,40);
  RooRealVar* muPlusPt = new RooRealVar("muPlusPt","muPlusPt",muonpTcut,50);

  ////////CHECK THIS EACH TIME YOU WANT TO FIT STHG NEW !!!!!!!///////////////////


  RooRealVar* muPlusEta = new RooRealVar("muPlusEta","muPlusEta",-1,1);
  RooRealVar* muMinusPt = new RooRealVar("muMinusPt","muMinusPt",muonpTcut,50);
  RooRealVar* muMinusEta = new RooRealVar("muMinusEta","muMinusEta",-1,1);
  //updated rapidity to asymmetric collision
  if(area==3) { RooRealVar* upsRapidity = new RooRealVar("upsRapidity",  "upsRapidity"  ,-1.47,0.53);} 
  else { RooRealVar* upsRapidity = new RooRealVar("upsRapidity",  "upsRapidity"  ,-1,1);}
  //////////////THANK YOU////////////////////////////////////////////////////////

  //Declare war to data
  RooDataSet* data0, *data, *likesignData0, *likesignData, *TrkRotData0, *TrkRotData;
  //import data set
  if (area==2) data0 = new RooDataSet("data","data",theTree,RooArgSet(*mass,*upsRapidity,*vProb,*upsPt,*Centrality,*muMinusEta,*muPlusEta,*muPlusPt,*muMinusPt));//put QQsign out, doesn't change a thing
  else data0 = new RooDataSet("data","data",theTree,RooArgSet(*mass,*upsRapidity,*vProb,*upsPt,*muPlusEta,*muMinusEta,*muPlusPt,*muMinusPt,*QQsign));
  data0->Print();
  data = ( RooDataSet*) data0->reduce(Cut(cut_));
  data->Print();

  // import like-sign data set
  // too many variables...
  if (area==2) likesignData0 = new RooDataSet("likesignData","likesignData",allsignTree,RooArgSet(*mass,*Centrality,*vProb,*upsRapidity,*QQsign,*muMinusEta,*muPlusEta,*muPlusPt,*muMinusPt));
  else likesignData0 = new RooDataSet("likesignData","likesignData",allsignTree,RooArgSet(*mass,*upsRapidity,*vProb,*muPlusPt,*muMinusPt,*upsPt,*QQsign,*muMinusEta,*muPlusEta));
  likesignData0->Print();
  likesignData = ( RooDataSet*) likesignData0->reduce(Cut(cut_+"&&QQsign!=0"));
  likesignData->Print();

  
  //import track-rotation data set
  if (TRKROT) {
    if (area==2) TrkRotData0 = new RooDataSet("TrkRotData","TrkRotData",trkRotTree,RooArgSet(*mass,*upsRapidity,*muPlusPt,*muMinusPt,*vProb,*upsPt,*Centrality,*QQsign,*muMinusEta,*muPlusEta));
    else TrkRotData0 = new RooDataSet("TrkRotData","TrkRotData",trkRotTree,RooArgSet(*mass,*upsRapidity,*muPlusPt,*muMinusPt,*upsPt,*vProb,*QQsign,*muMinusEta,*muPlusEta));
    TrkRotData0->Print();
    if (PR_plot && RAA) TrkRotData = ( RooDataSet*) TrkRotData0->reduce(Cut(cut_+" && upsPt < 8.1"));
    else if (PR_plot && !RAA) TrkRotData = ( RooDataSet*) TrkRotData0->reduce(Cut(cut_+" && upsPt < 7.07"));
    else TrkRotData = ( RooDataSet*) TrkRotData0;//->reduce(Cut(cut_+" && QQsign != 0"));
    TrkRotData->Print();
  }
  //discuss this, read BPH now ! erase when done !
  // nothing to deal with states, just SB|DATA|SB at first but now i use R2 for plotting
  // my own mass range
  mass->setRange("R1",8.5,10.2);
  mass->setRange("R2",8.5,11.5);   
  mass->setRange("R3",10.2,11.5);
  
  const double M1S = 9.46;   //upsilon 1S pgd mass value
  const double M2S = 10.02;  //upsilon 2S pgd mass value
  const double M3S = 10.35;  //upsilon 3S pgd mass value

  RooRealVar *mean    = new RooRealVar("#mu_{#Upsilon(1S)}","#Upsilon mean",M1S,M1S-0.1,M1S+0.1);


  RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",
					    RooArgList(*mean));
  RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0*10.023/9.460",
					    RooArgList(*mean));
  RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@00*10.355/9.460",
					    RooArgList(*mean));

  RooRealVar *sigma1 = new RooRealVar("sigma1","Sigma_1",0.10,0.01,0.30);    //detector resolution
  RooRealVar *sigma2 = new RooRealVar("#sigma_{#Upsilon(1S)}","Sigma_1S",0.08,0.01,0.30); //Y(1S) resolution
  RooFormulaVar *reso1S = new RooFormulaVar("reso1S","@0",RooArgList(*sigma2));
  RooFormulaVar *reso2S = new RooFormulaVar("reso2S","@0*10.023/9.460",RooArgList(*sigma2));
  RooFormulaVar *reso3S = new RooFormulaVar("reso3S","@0*10.355/9.460",RooArgList(*sigma2));
  
  /// to describe final state radiation tail on the left of the peaks
  //RooRealVar *alpha  = new RooRealVar("alpha","tail shift",1.4,0,4);   // minbias fit value
  // for PbPb : alpha 3.9 +/- 0.9
  RooRealVar *alpha  = new RooRealVar("alpha","tail shift",0.982,0,4);   // MC value
  // RooRealVar *npow   = new RooRealVar("npow","power order",2.3,1,3);       // zhen's setting MC value
  RooRealVar *npow   = new RooRealVar("npow","power order",2.3,1,110);   //  my settings  // MC value  
// npow->setVal(1.57);
  // if (!fitMB)
   switch(area)
    {
    case 0:
      //    alpha->setVal(1.68);
      break;
    case 1:
      alpha->setVal(1.783);
      npow->setVal(1.57);     // pp 5 TeV MC values.
      //alpha->setVal(1.4);   // MC pp 7 TeV
      //npow->setVal(2.3); 
      break;
    case 2:
      //   alpha->setVal(1.39); // good MinBias Regit value (vProb > 0.01)
      //  npow->setVal(1.52);  // same. erf exp fit for both values
	  alpha->setVal(1.4); 
	  npow->setVal(2.3);  
      break;
    case 3:
      // alpha->setVal(1.3);
      //  npow->setVal(2.0);
      break;
    default:break; 
    }
   //alpha->setConstant(kTRUE); //FOR THE MOMENT
   //npow ->setConstant(kTRUE); //FOR THE MOMENT


  // relative fraction of the two peak components 
  RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.3,0.,1.);
  sigmaFraction->setVal(0);
  sigmaFraction->setConstant(kTRUE);
  
  // Upsilon 1S
  // RooCBShape  *gauss1S1 = new RooCBShape ("gauss1S1", "FSR cb 1s",
  //					  *mass,*mean1S,*sigma1,*alpha,*npow);
  RooCBShape  *gauss1S2 = new RooCBShape ("gauss1S2", "FSR cb 1s",
					  *mass,*mean1S,*reso1S,*alpha,*npow);
  //RooAddPdf *sig1S      = new RooAddPdf  ("sig1S","1S mass pdf",
  //					  RooArgList(*gauss1S1,*gauss1S2),*sigmaFraction);

  // mean->setVal(9.46030);
  //mean->setConstant(kTRUE);
  sigma1->setVal(0);
  sigma1->setConstant(kTRUE);
  
  // was going with if (!fitMB) statement
  // sigma2->setVal(width_);        //fix the resolution
  //   sigma2->setVal(0.0624);        //fix the resolution to minBias
  //sigma2->setVal(0.06229);        //fix the resolution to MC 5 TeV
  sigma2->setVal(0.062);         // MC pp 7TeV  (AN2011-455)
  // sigma2->setVal(0.092);            // nominal width fixing (MC) in AN2011-455
  //sigma2->setConstant(kTRUE);
 
  
  /// Upsilon 2S
  // RooCBShape  *gauss2S1 = new RooCBShape ("gauss2S1", "FSR cb 2s", 
  //					  *mass,*mean2S,*sigma1,*alpha,*npow); 
  RooCBShape  *gauss2S2 = new RooCBShape ("gauss2S2", "FSR cb 2s", 
					  *mass,*mean2S,*reso2S,*alpha,*npow); 
  //RooAddPdf *sig2S      = new RooAddPdf  ("sig2S","2S mass pdf",
  //					  RooArgList(*gauss2S1,*gauss2S2),*sigmaFraction);

  /// Upsilon 3S
// RooCBShape  *gauss3S1 = new RooCBShape ("gauss3S1", "FSR cb 3s", 
//					  *mass,*mean3S,*sigma1,*alpha,*npow); 
  RooCBShape  *gauss3S2 = new RooCBShape ("gauss3S2", "FSR cb 3s", 
					  *mass,*mean3S,*reso3S,*alpha,*npow); 
  //RooAddPdf *sig3S      = new RooAddPdf  ("sig3S","3S mass pdf",
  //                                        RooArgList(*gauss3S1,*gauss3S2),*sigmaFraction);

  // ohter bkgds
  RooRealVar *bkg_a1  = new RooRealVar("bkg_{a1}", "background a1", 0, -5, 5);
  RooRealVar *bkg_a2  = new RooRealVar("bkg_{a2}", "background a2", 0, -5, 5);
  RooRealVar *bkg_a3  = new RooRealVar("bkg_{a3}", "background a3", 0, -2, 2);
  RooAbsPdf  *bkgPdf  = new RooChebychev("bkg","background",
					 *mass, RooArgList(*bkg_a1,*bkg_a2,*bkg_a3));
  
  // only sideband region pdf, using RooPolynomial instead of RooChebychev for multiple ranges fit
  RooRealVar *SB_bkg_a1  = new RooRealVar("SB bkg_{a1}", "background a1", 0, -1, 1);
  RooRealVar *SB_bkg_a2  = new RooRealVar("SB bkg_{a2}", "background a2", 0, -1, 1);
  RooAbsPdf  *SB_bkgPdf  = new RooPolynomial("SB_bkg","side-band background",
					     *mass, RooArgList(*SB_bkg_a1,*SB_bkg_a2));
  //SB_bkg_a1->setVal(0);
  //SB_bkg_a1->setConstant(kTRUE);
  //SB_bkg_a2->setVal(0);
  //SB_bkg_a2->setConstant(kTRUE);
  
  /// Combined pdf
  // nt (PbPb) can be less 
  int nt = 100000;
  //bool fitfraction = true;
  RooRealVar *nbkgd = new RooRealVar("N_{bkg}","nbkgd",nt*0.75,0,10*nt); // pp 2.76 ?
  RooRealVar *SB_nbkgd = new RooRealVar("SB N_{bkg}","SB_nbkgd",nt*0.75,0,10*nt);
  RooRealVar *nsig1f  = new RooRealVar("N_{#Upsilon(1S)}","nsig1S",nt*0.25,0,10*nt);
   
  //use the YIELDs of 2S and 3S as free parameters
 
  //  RooRealVar *nsig2f  = new RooRealVar("N_{#Upsilon(2S)}","nsig2S",   nt*0.25,0,10*nt);
  //RooRealVar *nsig3f  = new RooRealVar("N_{#Upsilon(3S)}","nsig3S",   nt*0.25,0,10*nt);
  
  //use the RATIOs of 2S and 3S as free parameters
  
  RooRealVar *f2Svs1S = new RooRealVar("N_{2S}/N_{1S}","f2Svs1S",0.2,-0.1,1);  
  //RooRealVar *f23vs1S = new RooRealVar("N_{2S+3S}/N_{1S}","f23vs1S",0.45,-0.1,1);
    RooRealVar *f3Svs1S = new RooRealVar("N_{3S}/N_{1S}","f3Svs1S",0.10,0.,0.5);
  
   RooFormulaVar *nsig2f = new RooFormulaVar("nsig2S","@0*@1", RooArgList(*nsig1f,*f2Svs1S));
    RooFormulaVar *nsig3f = new RooFormulaVar("nsig3S","@0*@1", RooArgList(*nsig1f,*f3Svs1S));
    //RooFormulaVar *nsig3f = new RooFormulaVar("nsig3S","@0*@2-@0*@1", RooArgList(*nsig1f,*f2Svs1S,*f23vs1S));





   // f3Svs1S->setConstant(kTRUE);
  
  //force the ratio to the pp value
  f2Svs1S_pp->setConstant(kTRUE);
  f3Svs1S_pp->setConstant(kTRUE);
  RooFormulaVar *nsig2f_ = new RooFormulaVar("nsig2S_pp","@0*@1", RooArgList(*nsig1f,*f2Svs1S_pp)); 
  RooFormulaVar *nsig3f_ = new RooFormulaVar("nsig3S_pp","@0*@1", RooArgList(*nsig1f,*f3Svs1S_pp)); 
  
  //only sideband region pdf, using RooPolynomial instead of RooChebychev for multiple ranges fit
  RooAbsPdf  *SB_pdf = new RooAddPdf ("SB_pdf","sideband background pdf",
				      RooArgList(*SB_bkgPdf),
				      RooArgList(*SB_nbkgd));
  //only signal region pdf, using RooPolynomial instead of RooChebychev for multiple ranges fit
  RooAbsPdf  *S_pdf   = new RooAddPdf ("S_pdf","total signal+background pdf",
				       RooArgList(*gauss1S2,*gauss2S2,*gauss3S2,*SB_bkgPdf),
				       RooArgList(*nsig1f,*nsig2f,*nsig3f,*SB_nbkgd));
  //parameters for likesign
  RooRealVar m0shift("turnOn","turnOn",9.5,0,30) ; // turn On mass of error function
  RooRealVar width("width","width",2.36,0,100) ; // width of error fct
  RooRealVar par3("decay","decay",2, 0,30) ; // decay constant of exp
  RooGaussian* m0shift_constr;
  RooGaussian* width_constr;
  RooGaussian* par3_constr;

  RooRealVar *nLikesignbkgd = new RooRealVar("NLikesign_{bkg}","nlikesignbkgd",nt*0.75,0,10*nt);
  if (TRKROT) {
    nLikesignbkgd->setVal(TrkRotData->sumEntries());
    nLikesignbkgd->setError(sqrt(TrkRotData->sumEntries()));
    
  }
  else {
    nLikesignbkgd->setVal(likesignData->sumEntries());
    nLikesignbkgd->setError(sqrt(likesignData->sumEntries()));
  }
  nbkgd->setVal(data->sumEntries());
  nbkgd->setError(sqrt(data->sumEntries()));
  // nbkgd->setVal(data->sumEntries());
  // nbkgd->setError(sqrt(data->sumEntries()));

  // LIKE SIGN CONSTRAINTS : Use of plain gaussian PDF, N:Nlikesignbkgd, mean:the value of the LS Nbkgd, sigma: error on Nbkgd ( gaussian error sqrt(N) ?? )

  if (LS_constrain) {
    RooGaussian* nLikesignbkgd_constr = new RooGaussian("nLikesignbkgd_constr","nLikesignbkgd_constr",*nLikesignbkgd,RooConst(nLikesignbkgd->getVal()),RooConst(nLikesignbkgd->getError()));
  }
  else nLikesignbkgd->setConstant(kTRUE);
  
  // RESIDUAL BKGD : nBkgd - nLikesignbkgd
  cout << "nbkgd = " << nbkgd->getVal() << " ;  nLikesignbkgd = " << nLikesignbkgd->getVal() << endl;
  RooFormulaVar *nResidualbkgd = new RooFormulaVar("NResidual_{bkg}","@0-@1",RooArgList(*nbkgd,*nLikesignbkgd));
  cout << "nResidualbkgd = " << nResidualbkgd->getVal() << endl ;
  // BACKGROUND MODEL DEFINITION 
  
  switch (bkgdModel) {
  case 1 : 
    RooGenericPdf *LikeSignPdf = new  RooGenericPdf("Like-sign","likesign","exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",RooArgList(*mass,m0shift,width,par3));
     RooAbsPdf  *this_bkgd   = new RooAddPdf ("this_bkgd","Background1",
						    RooArgList(*LikeSignPdf),
						    RooArgList(*nLikesignbkgd));
     //TROUBLED AREA 
    if (TRKROT) {
      if(pp7){
	RooFitResult* fit_1st = this_bkgd->fitTo(*TrkRotData,Save()) ;
      }
      else {
	RooFitResult* fit_1st = LikeSignPdf->fitTo(*TrkRotData,NumCPU(8),Timer(kTRUE),Save()) ;
      }
    }
    else {
      if(pp7){
	RooFitResult* fit_1st = this_bkgd->fitTo(*likesignData,Save()) ;
      }
      else {
	RooFitResult* fit_1st = this_bkgd->fitTo(*likesignData,NumCPU(8),Timer(kTRUE),Save()) ;
      }
    }
     //use error function to fit the like-sign or TrkRot data, then : 
    //if LS_CONSTRAINT is on, use constrained values : for @1 @2 @3, build gaussian pdfs with parameter mean=@x, sigma=error(@x)

    // likesign data
    //LikeSignPdf.fitTo(*data) ;       //fit to unlikesign data, can always try!
    //fit_1st->Print();
    
    if (LS_constrain) {
      m0shift_constr = new RooGaussian("m0shift_constr","m0shift_constr",m0shift,RooConst(m0shift.getVal()),RooConst(m0shift.getError()));
      width_constr   = new RooGaussian("width_constr","width_constr",width,RooConst(width.getVal()),RooConst(width.getError()));
      par3_constr    = new RooGaussian("par3_constr","par3_constr",par3,RooConst(par3.getVal()),RooConst(par3.getError()));
    }
    else {
      m0shift.setConstant(kTRUE);
      width.setConstant(kTRUE);
      par3.setConstant(kTRUE);
    }
    bkg_a3->setVal(0);
    bkg_a3->setConstant(kTRUE);
    // make a linear combination of the (polynomial bkg)*(nBkgd-nLS) + (ErrorPdf)*(nLS)
     RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						    RooArgList(*bkgPdf,*LikeSignPdf),
						    RooArgList(*nResidualbkgd,*nLikesignbkgd));
    break;
  case 2 : 
 //use RooKeysPdf to smooth the like-sign, then fix the shape and normalization
 // Can also put some boolean args, what pick numbers for arguments?
    bkg_a3->setVal(0);
    bkg_a3->setConstant(kTRUE);
    if (TRKROT) RooKeysPdf *LikeSignPdf = new RooKeysPdf("Like-sign","likesign",*mass,*TrkRotData,3,1.5);
    else{ RooKeysPdf *LikeSignPdf = new RooKeysPdf("Like-sign","likesign",*mass,*likesignData,3,1.7);}
    RooAbsPdf  *pdf_combinedbkgd = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						    RooArgList(*bkgPdf,*LikeSignPdf),
						    RooArgList(*nResidualbkgd,*nLikesignbkgd));
    break;
  case 3 : //use error function to fit the whole bkg directly

    RooGenericPdf *LikeSignPdf = new RooGenericPdf("Like-sign","likesign","exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",RooArgList(*mass,m0shift,width,par3));
    RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						    RooArgList(*LikeSignPdf),
						    RooArgList(*nbkgd));
    break;
  case 4 :
    bkg_a3->setVal(0);
    bkg_a3->setConstant(kTRUE);
    RooAbsPdf *pdf_combinedbkgd = new RooAddPdf("pdf_combinedbkgd","total combined background pdf",
						RooArgList(*bkgPdf),
						RooArgList(*nbkgd));
  
    break;
  case 5 : //use lin.Comb :(nLS)*(erfExpPdf)+(nResidual)*(BkgdPdf)
  
    RooGenericPdf *LikeSignPdf = new RooGenericPdf("Like-sign","likesign","exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",RooArgList(*mass,m0shift,width,par3));
   bkg_a3->setVal(0);
   bkg_a3->setConstant(kTRUE);
    RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						    RooArgList(*LikeSignPdf,*bkgPdf),
						    RooArgList(*nbkgd,*nbkgd));
    break;
  case 6 : // pol3 chebychev
    RooGenericPdf *LikeSignPdf = new RooGenericPdf("Like-sign","likesign","exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",RooArgList(*mass,m0shift,width,par3));
    m0shift->setVal(0); m0shift->setConstant(kTRUE);
    par3->setVal(0);par3->setConstant(kTRUE);
    width->setVal(100);width->setConstant(kTRUE);		       
    nLikesignbkgd->setVal(0);
    nLikesignbkgd->setConstant(kTRUE);
    RooAbsPdf *pdf_combinedbkgd = new RooAddPdf("pdf_combinedbkgd","total combined background pdf",
						RooArgList(*bkgPdf,*LikeSignPdf),
						RooArgList(*nbkgd,*nLikesignbkgd));
	break;
  case 7 : // use erf*Exp+pol3
    RooGenericPdf *LikeSignPdf = new RooGenericPdf("Like-sign","likesign","exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",RooArgList(*mass,m0shift,width,par3));
    RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						    RooArgList(*LikeSignPdf,*bkgPdf),
						    RooArgList(*nLikesignbkgd,*nResidualbkgd));
  default :
    break;
  }

   
  //pdf with fixed ratio of the pp ratio
  // not used in fits right now.
  // RooAbsPdf  *pdf_pp   = new RooAddPdf ("pdf_pp","total signal+background pdf",
  //					RooArgList(*gauss1S2,*sig2S,*sig3S,*pdf_combinedbkgd),
  //					RooArgList(*nsig1f,*nsig2f_,*nsig3f_,*nbkgd));
  
  //the nominal fit with default pdf 
  if (LS_constrain) {
    RooAbsPdf  *pdf_unconstr   = new RooAddPdf ("pdf_unconstr","total signal+background pdf",
						RooArgList(*gauss1S2,*sig2S,*sig3S,*pdf_combinedbkgd),
						RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
    RooProdPdf *pdf  = new RooProdPdf ("pdf","total constr pdf",
				       RooArgSet(*pdf_unconstr,*m0shift_constr,*width_constr,*par3_constr,*nLikesignbkgd_constr));
    if(area==1){
      RooFitResult* fit_2nd = pdf->fitTo(*data,NumCPU(8),Timer(kTRUE),Range("R2"),Constrained(),Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    } 
    else {
      RooFitResult* fit_2nd = pdf->fitTo(*data,NumCPU(8),Timer(kTRUE),Constrained(),Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    }
  }
  else {
     RooAbsPdf  *pdf   = new RooAddPdf ("pdf","total signal+background pdf",
				       RooArgList(*gauss1S2,*gauss2S2,*gauss3S2,*pdf_combinedbkgd),
				       RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
    
    
    if(area==1){
      RooFitResult* fit_2nd = pdf->fitTo(*data,NumCPU(8),Timer(kTRUE),Range("R2"),Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    }
    else {
      RooFitResult* fit_2nd = pdf->fitTo(*data,NumCPU(8),Timer(kTRUE),Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    }
  }
  
  
   
  //draw the data and fit lines and pull and save plots AAAAAAAARRRGHHGRHHHMFPPNMF
  TCanvas c; c.cd();
 
  // if (area==1) int nbins = 2*ceil((mmax_-mmin_)/binw_); 
  int nbins = ceil((mmax_-mmin_)/binw_);
  RooPlot* frame = mass->frame(Bins(nbins),Range(mmin_,mmax_));
  TPad *pPad1 = new TPad("pPad1","This is pad1",0.05,0.35,0.95,0.97);
  pPad1->SetBottomMargin(0.00);
  pPad1->Draw();

  TPad *pPad2 = new TPad("pPad2","This is pad2",0.05,0.05,0.95,0.35);
  pPad2->SetTopMargin(0.0);
  pPad2->Draw();

  pPad1->cd();
  frame->Draw();
  
 
  data->plotOn(frame,Name("theData"),MarkerSize(0.8)); // plotting data
  pdf->plotOn(frame,Name("thePdf")); //plotting pdf
 
   RooArgSet * pars = pdf->getParameters(data);
  //RooArgSet * pars = LikeSignPdf->getParameters(likesignData);
  
  // if plotting ratios in a pdf ... gets old
  //if(yields)    rapbin_yields(bin, nsig1f, nsig2f,nsig3f);
  //  else if(ratios) rapbin(bin,f2Svs1S,f3Svs1S);
  
      //calculate chi2 in a mass range
  //----pull histogram
   RooHist *phPullm  = frame->pullHist(0,0,true); 
  phPullm->SetName("phPullm");
  double *ypull     = phPullm->GetY();
  
  TH1 *phData       = data->createHistogram("invariantMass",nbins); // create pull histo on data
  double Chi2       = 0;
  int nFullBinsPull = 0;
  for (int i=0; i < nbins+1; i++) 
    {
      if (phData->GetBinContent(i) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i],2);
    }

  //------------------------------------
  // for writing on canvas ---- new Chi2 version
  int nFitParam     = fit_2nd->floatParsFinal().getSize();
  double chi2FromRoo = frame->chiSquare(nFitParam); // reduced chi2
  int Dof           = nFullBinsPull - nFitParam; // nDof
  double myChi2 = chi2FromRoo*Dof; // true chi2
  double UnNormChi2 = Chi2; // same as it should be
  //-------------------------------------------
  //plot parameters
  float baseNll = fit_2nd->minNll();
  if(plotpars) { // not even sure we go in this condition...yes we do !!! fucker !
    paramOn_ = "_p";
    // may have to remove nll label here .. use made-up chi2 to pick best fit // 
    pdf->paramOn(frame,Layout(0.6,0.935,0.97)); // create your parameters' box on the frame of pdf draw

  }

  pPad2->cd();
  double mFrameMax = 0;
  RooPlot* prpFramePull = mass->frame(Title("Pull"),Bins(nbins),Range(mass_l,mass_h));
  prpFramePull->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  prpFramePull->GetXaxis()->CenterTitle(kTRUE);
  prpFramePull->GetXaxis()->SetLabelSize(0.085);
  // prpFramePull->GetXaxis()->SetTitleSize(0.08);
  prpFramePull->GetXaxis()->SetTitleOffset(1.3);
  // prpFramePull->GetYaxis()->SetMarkerSize(0.8);
  prpFramePull->GetYaxis()->SetLabelSize(0.085);
  prpFramePull->GetYaxis()->SetTitle("Pull");
  //  prpFramePull->GetYaxis()->SetTitleSize(0.01);
  prpFramePull->GetYaxis()->SetTitleOffset(1.3);
  prpFramePull->addPlotable(phPullm,"PX");  // take the pull histos to go on to this frame, in pPad2.
 
  if (prpFramePull->GetMinimum()*-1 > prpFramePull->GetMaximum()) mFrameMax = prpFramePull->GetMinimum()*-1;
  else mFrameMax = prpFramePull->GetMaximum();
  prpFramePull->SetMaximum(mFrameMax); 
  prpFramePull->SetMinimum(-1*mFrameMax); 
  prpFramePull->Draw();  // here, the pull should be drawn
  //-------------------------------------------

  //plot data
  pPad1->cd();
  data->plotOn(frame,Name("theData"),MarkerSize(0.8)); //plotting data (again?)
  // pdf->plotOn(frame,Components("bkg"),Name("theBkg"),LineStyle(5),LineColor(kGreen));
  pdf->plotOn(frame,Components("pdf_combinedbkgd"),LineStyle(kDashed)); // plotting the dashed background.
   if (plotLikeSign) {
     if (TRKROT) pdf->plotOn(frame,Components("Like-sign"),Name("theLikeSign"),LineStyle(9),LineColor(kMagenta)); // plot trkrot
     // else  pdf->plotOn(frame,Components("Like-sign"),Name("theLikeSign"),LineStyle(9),LineColor(kRed)); // plot the likesign pdf that we don't want now.
 }
 if (plotLikeSign) {
    if (TRKROT) TrkRotData->plotOn(frame,Name("theLikeSignData"),MarkerSize(0.8),MarkerColor(kMagenta),MarkerStyle(22));
    else likesignData->plotOn(frame,Name("theLikeSignData"),MarkerSize(0.8),MarkerColor(kRed),MarkerStyle(24));
    //LikeSignPdf->plotOn(frame,Name("theLikeSign"),VisualizeError(*fit_1st,1),FillColor(kOrange));
    //LikeSignPdf->plotOn(frame,Name("theLikeSign"),LineColor(kRed));

    // pdf->plotOn(frame,Name("thePdf"));


  }   
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetYaxis()->SetTitleOffset(1.3);
  if( area==1 ) frame->GetYaxis()->SetTitle("Entries / (0.05 GeV/c^{2})");
  else frame->GetYaxis()->SetTitle("Entries / (0.14 GeV/c^{2})");
  if (PR_plot && RAA) frame->GetYaxis()->SetRangeUser(0,1200);
  //frame->GetYaxis()->SetLabelSize(0.05);
  frame->Draw(); // now, the datapoints, pdf and pull should be drawn...
  //--------------------------------------------------------------------------------


  // cosmetics , legends
  
  //plot parameters ----- This is oooOLD
  if(!plotpars) {
    paramOn_ = ""; 
    TLatex latex1;
    latex1.SetNDC();
    if (PbPb) {
      switch (bin) {  
	if(bin<8){ latex1.DrawLatex(0.46,1.-0.05*3,"CMS pPb  #sqrt{s} = 5 TeV");
	  latex1.DrawLatex(0.5,1.-0.05*4.9,"L_{int} = 10.7 nb^{-1}"); }
	else{ latex.DrawLatex(0.46,1.-0.05*3,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV");
	  latex1.DrawLatex(0.5,1.-0.05*4.9,"L_{int} = 150 #mub^{-1}"); }
	
      case  0: latex1.DrawLatex(0.5,1.-0.05*6.2," -1.93 < y_{coll} < -0.47"); break;
      case  1: latex1.DrawLatex(0.5,1.-0.05*6.2," -0.47 < y_{coll} < 0"); break;
      case  2: latex1.DrawLatex(0.5,1.-0.05*6.2," 0 < y_{coll} < 0.47"); break;
      case  3: latex1.DrawLatex(0.5,1.-0.05*6.2," 0.47 < y_{coll} < 1.93"); break;
      case  4: latex1.DrawLatex(0.5,1.-0.05*6.2," -1.93 < y_{coll} < 0"); break;
      case  5: latex1.DrawLatex(0.5,1.-0.05*6.2," 0 < y_{coll} < 1.93"); break;
      case  6: latex1.DrawLatex(0.5,1.-0.05*6.2," |y_{coll}| < 1.93"); break;
      case  7: latex1.DrawLatex(0.5,1.-0.05*6.2,"|y_{coll}|<1.93, 60-100 %"); break;
      case  8: latex1.DrawLatex(0.5,1.-0.05*6.2,"|y_{coll}|<1.93, 0-10%"); break;
      default; break;
      }   
    }
    else {
      latex1.DrawLatex(0.46,1.-0.05*3,"CMS pp #sqrt{s_{NN}} = 7 TeV");
      latex1.DrawLatex(0.5,1.-0.05*4.2,"L_{int} = 621 fb^{-1}");
      latex1.DrawLatex(0.55,1.-0.05*5.4,"|y| < 1")
      // switch (bin) {  
      // case  0: latex1.DrawLatex(0.5,1.-0.05*6.2," -1.93 < y_{coll} < -0.47"); break;
      // case  1: latex1.DrawLatex(0.5,1.-0.05*6.2," -0.47 < y_{coll} < 0"); break;
      // case  2: latex1.DrawLatex(0.5,1.-0.05*6.2," 0 < y_{coll} < 0.47"); break;
      // case  3: latex1.DrawLatex(0.5,1.-0.05*6.2," 0.47 < y_{coll} < 1.93"); break;
      // case  4: latex1.DrawLatex(0.5,1.-0.05*6.2," -1.93 < y_{coll} < 0"); break;
      // case  5: latex1.DrawLatex(0.5,1.-0.05*6.2," 0 < y_{coll} < 1.93"); break;
      // case  6: latex1.DrawLatex(0.5,1.-0.05*6.2," |y_{coll}| < 1.93"); break;
      // default; break;
      // }
    }
   
     //  latex1.DrawLatex(0.5,1.-0.05*7.5,"p_{T}^{#mu} > 4 GeV/c");
    //latex1.DrawLatex(0.7,1.-0.05*8.8,"data");
    // latex1.DrawLatex(0.7,1.-0.05*9.8,"total fit");
    //latex1.DrawLatex(0.7,1.-0.05*10.8,"background");
    //TMarker M(11.35, 425, 20);
    //M.SetMarkerSize(1);
    //M.DrawMarker(11.35, 425);
    
    TLine L1; 
    L1.SetLineWidth(2);
    L1.SetLineColor(4);
    L1.DrawLineNDC(0.62, 1. - 0.05*9.55,
		   0.68, 1. - 0.05*9.55);
    
    TLine L2; 
    L2.SetLineWidth(2);
    L2.SetLineColor(4);
    L2.SetLineStyle(2);
    L2.DrawLineNDC(0.62, 1. - 0.05*10.55,
		   0.68, 1. - 0.05*10.55);
  } 
 TLatex latexb;
    latexb.SetNDC();
    //latexb.SetFont(32);
    latexb.DrawLatex(1-0.82,1.-0.05*3,"CMS pp");
    latexb.DrawLatex(1-0.82,1.-0.05*4,"#sqrt{s} = 7 TeV");
    latexb.DrawLatex(1-0.82,1.-0.05*5,"L_{int} =  5.4 pb^{-1}"); 
 
    latexb.DrawLatex(1-0.82,1.-0.05*6,"|y| < 1"); 
  pPad2->cd();
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.085);
  latex2.DrawLatex(0.7,1.-0.05*3.5,Form("#chi^{2}/ndf = %2.1f/%d",UnNormChi2,Dof));
  // now, the new Chi2 should be drawn....
 
  // c.SaveAs(figs_+figName_+paramOn_+suffix_+"2901.gif");
  cout << figName_ << endl;
  c.SaveAs(figs_+figName_+paramOn_+suffix_+".pdf");
    
  

  /*
    mass->setRange("R1S",8.8,9.7);
    mass->setRange("R2S",9.8,10.2);
    //pdf_combinedbkgd->fitTo(*data,Range("R1,R3"),Constrained(),Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    RooAbsReal* integral_1S = pdf_combinedbkgd->createIntegral(*mass,NormSet(*mass),Range("R1S")) ;
    cout << "1S bkgd integral = " << integral_1S->getVal() * (nbkgd->getVal()) << endl ;
    RooAbsReal* integral_2S = pdf_combinedbkgd->createIntegral(*mass,NormSet(*mass),Range("R2S")) ;
    cout << "2S bkgd integral = " << integral_2S->getVal() * (nbkgd->getVal()) << endl ;
    cout << "1S range count: " << data->sumEntries("invariantMass","R1S") <<endl;
    cout << "2S range count: " << data->sumEntries("invariantMass","R2S") <<endl;
    cout << "1S signal yield: " << data->sumEntries("invariantMass","R1S") - integral_1S->getVal() * (nbkgd->getVal()) << endl;
    cout << "2S signal yield: " << data->sumEntries("invariantMass","R2S") - integral_2S->getVal() * (nbkgd->getVal()) << endl;
  */

 
    if (PR_plot) {
      //overlay the pp ratio, calculate the significance of suppression
      m0shift.setConstant(kTRUE);
      width.setConstant(kTRUE);
      par3.setConstant(kTRUE);
      nbkgd->setConstant(kTRUE);
      sigma2->setConstant(kTRUE);
      alpha->setConstant(kTRUE);
      mean->setConstant(kTRUE);
      //  TrkRotData->plotOn(frame,Name("theTrkRotData"),MarkerSize(0.0),MarkerColor(kWhite),LineColor(kWhite),MarkerStyle(22));
      if (PbPb) suppression(nsig1f, fit_2nd, nfloatpars, pdf, pdf_pp, data, frame);
    }
    
    //setup the limits
 
    //setUpLimits(-0.08, 0.55, pdf, data, f23vs1S, baseNll);
    //setUpLimits(-0.1, 0.9, pdf, data, f2Svs1S, baseNll);
    
    //calculate the confidence interval with RooStats
    //ConfidencInterval(0.95, f3Svs1S, data, pdf);
    //ConfidencInterval(0.68, f2Svs1S, data, pdf);
    cout<<"!!!!! nFullBinsPull="<<nFullBinsPull<<"\tnFitParam="<<nFitParam<<endl;
 
    if(area!=2){
      outfile<<suffix_<<"+"<<figName_<<" "<<f2Svs1S->getVal()<<" "<<f2Svs1S->getError()<<" "<<f3Svs1S->getVal()<<" "<<f3Svs1S->getError()<<" "<<npow->getVal()<<" "<<npow->getError()<<" "<<alpha->getVal()<<" "<<alpha->getError()<<" "<<sigma2->getVal()<<" "<<sigma2->getError()<<" "<<2*nFitParam+2*baseNll<<" "<<fit_2nd->edm()<<" "<<myChi2<<" "<<UnNormChi2<<" "<<UnNormChi2/Dof<<" "<<TMath::Prob(UnNormChi2,Dof)<<" "<<Dof<<" "<<nFitParam<<" "<<baseNll<< endl;
    } 
    else if (area==2){
      outfile<<suffix_<<"+"<<figName_<<" "<<f2Svs1S->getVal()<<" "<<f2Svs1S->getError()<<" "<<f23vs1S->getVal()<<" "<<f23vs1S->getError()<<" "<<npow->getVal()<<" "<<npow->getError()<<" "<<alpha->getVal()<<" "<<alpha->getError()<<" "<<sigma2->getVal()<<" "<<sigma2->getError()<<" "<<2*nFitParam+2*baseNll<<" "<<fit_2nd->edm()<<" "<<UnNormChi2<<" "<<UnNormChi2/Dof<<" "<<TMath::Prob(UnNormChi2,Dof)<<" "<<Dof<<" "<<nFitParam<<" "<<baseNll<< endl;
    }


    // Print fit results 
    cout << endl << "figure name: "<< figName_ << endl << endl;
    cout << "the nominal fit with the default pdf "<< endl ;
    fit_2nd->Print() ;
    TString *fitResult ;
    fitResult=fit_2nd->Print();
    outfile<<fitResult<<endl;

    //   fit_2nd->Print(outfile); 
    // cout << "  free parameter = "<< nfloatpars << ", mychi2 = " << mychsq << ", ndof = " << myndof << endl << endl;
    cout << "nll = "<<baseNll<<endl;
    cout << "Chi2Prob= "<<TMath::Prob(UnNormChi2,Dof)<<endl; 
}
void rapbin(int bin, RooRealVar *f2Svs1S, RooRealVar *f3Svs1S){
  switch (bin)
    {
    case 0 :
      //RooArgSet *ratio2 = fit_2nd->getParameters(f2Svs1S);
      //     fit_2nd->Print();
      cout << " 2s = "<< f2Svs1S->getVal() << " +/- "<< f2Svs1S->getError()<< endl <<endl;
      x[0]=-1.2; ex[0]=0.73;
      yAxisRap2S[0]=f2Svs1S->getVal(); yAxisRapError2S[0]=f2Svs1S->getError();
      yAxisRap3S[0]=f3Svs1S->getVal(); yAxisRapError3S[0]=f3Svs1S->getError();
      break;
    case 1 :
      //RooArgSet *ratio2 = fit_2nd->getParameters(f2Svs1S);
      //     fit_2nd->Print();
      x[1]=-.235;ex[1]=.235;
      cout << " 2s = "<< f2Svs1S->getVal() << " +/- "<< f2Svs1S->getError()<< endl <<endl;
      yAxisRap2S[1]=f2Svs1S->getVal(); yAxisRapError2S[1]=f2Svs1S->getError();
      yAxisRap3S[1]=f3Svs1S->getVal(); yAxisRapError3S[1]=f3Svs1S->getError();
      break;
    case 2 :
      //RooArgSet *ratio2 = fit_2nd->getParameters(f2Svs1S);
      //     fit_2nd->Print();
      x[2]=0.235;ex[2]=.235;
      cout << " 2s = "<< f2Svs1S->getVal() << " +/- "<< f2Svs1S->getError()<< endl <<endl;
      yAxisRap2S[2]=f2Svs1S->getVal(); yAxisRapError2S[2]=f2Svs1S->getError();
      yAxisRap3S[2]=f3Svs1S->getVal(); yAxisRapError3S[2]=f3Svs1S->getError();
      break;
    case 3 :
      x[3]=1.2;ex[3]=.73;
      //RooArgSet *ratio2 = fit_2nd->getParameters(f2Svs1S);
      //     fit_2nd->Print();
      cout << " 2s = "<< f2Svs1S->getVal() << " +/- "<< f2Svs1S->getError()<< endl <<endl;
      yAxisRap2S[3]=f2Svs1S->getVal(); yAxisRapError2S[3]=f2Svs1S->getError();
      yAxisRap3S[3]=f3Svs1S->getVal(); yAxisRapError3S[3]=f3Svs1S->getError();
      /*	TGraphErrors *gr2 = new TGraphErrors(nBins,x,yAxisRap2S,ex,yAxisRapError2S);
	TGraphErrors *gr3 = new TGraphErrors(nBins,x,yAxisRap3S,ex,yAxisRapError3S);
	TF1 *Axis= new TF1("Axis","0",-2.4,2.4);
        TCanvas cRap; cRap.cd();
	gr2->SetName("Upsilon(2S) ratio");
	gr3->SetName("Upsilon(3S) ratio");
	gr2->SetTitle("Upsilon(2S) ratio");
	gr3->SetTitle("Upsilon(3S) ratio");
	gr2->SetMarkerColor(4);
	gr3->SetMarkerColor(3);
	gr2->SetMarkerStyle(21);
	gr3->SetMarkerStyle(21);
	Axis->GetXaxis()->SetTitle("Rapidity");
	Axis->GetYaxis()->SetTitle("2S/1S, 3S/1S ratios");
	Axis->Draw("");
	gr2->Draw("Psame");
	gr3->Draw("Psame");*/
     	    TLatex latex1;
    latex1.SetNDC();
    if (PbPb) {
      latex1.DrawLatex(0.46,1.-0.05*3,"CMS pPb  #sqrt{s} = 5 TeV");
      latex1.DrawLatex(0.5,1.-0.05*4.9,"L_{int} = 10.66 #nb^{-1}"); 
      switch (bin) {  
      case  0: latex1.DrawLatex(0.5,1.-0.05*6.2," -1.93 < y_{coll} < -0.47"); break;
      case  1: latex1.DrawLatex(0.5,1.-0.05*6.2," -0.47 < y_{coll} < 0"); break;
      case  2: latex1.DrawLatex(0.5,1.-0.05*6.2," 0 < y_{coll} < 0.47"); break;
      case  3: latex1.DrawLatex(0.5,1.-0.05*6.2," 0.47 < y_{coll} < 1.93"); break;
      case  4: latex1.DrawLatex(0.5,1.-0.05*6.2," -1.93 < y_{coll} < 0"); break;
      case  5: latex1.DrawLatex(0.5,1.-0.05*6.2," 0 < y_{coll} < 1.93"); break;
      case  6: latex1.DrawLatex(0.5,1.-0.05*6.2," |y_{coll}| < 1.93"); break;
      default; break;
      }   
    }
    else {
      latex1.DrawLatex(0.46,1.-0.05*3,"CMS pp  #sqrt{s_{NN}} = 7 TeV");
      latex1.DrawLatex(0.5,1.-0.05*4.9,"L_{int} = 5.1 fb^{-1}"); 
      latex1.DrawLatex(0.5,1.-0.05*6.2," |y_{coll}| < 1.93");
      default; break;
      }
    
      latex1.DrawLatex(0.5,1.-0.05*7.5,"p_{T}^{#mu} > 4 GeV/c");
      //latex1.DrawLatex(0.7,1.-0.05*8.8,"data");
      latex1.DrawLatex(0.7,1.-0.05*9.8,"3S/1S");
      latex1.DrawLatex(0.7,1.-0.05*10.8,"2S/1S");
      //TMarker M(11.35, 425, 20);
      //M.SetMarkerSize(1);
      //M.DrawMarker(11.35, 425);
      
      TLine L1; 
      L1.SetLineWidth(2);
      L1.SetLineColor(4);
      L1.DrawLineNDC(1-0.62, 1. - 0.05*9.55,
		     1-0.68, 1. - 0.05*9.55);
      
      TLine L2; 
      L2.SetLineWidth(2);
      L2.SetLineColor(kRed);
      L2.DrawLineNDC(1-0.62, 1. - 0.05*10.55,
		     1-0.68, 1. - 0.05*10.55);
	if(PbPb) cRap.SaveAs(figs_+"ratiosVSrapidity_pPbPLOPnb_test.pdf");
	else ;//cRap.SaveAs(figs_+"ratiosVSrapidity_pp.pdf");
	break; 
    case 4 :
      x[0]=-.965;ex[0]=0.965;
      //RooArgSet *ratio2 = fit_2nd->getParameters(f2Svs1S);
      //     fit_2nd->Print();
      cout << " 2s = "<< f2Svs1S->getVal() << " +/- "<< f2Svs1S->getError()<< endl <<endl;
      yAxisRap2S[0]=f2Svs1S->getVal(); yAxisRapError2S[0]=f2Svs1S->getError();
      yAxisRap3S[0]=f3Svs1S->getVal(); yAxisRapError3S[0]=f3Svs1S->getError();
      break; 
    case 5 :
      x[1]=.965;ex[1]=0.965;
      //RooArgSet *ratio2 = fit_2nd->getParameters(f2Svs1S);
      //     fit_2nd->Print();
      cout << " 2s = "<< f2Svs1S->getVal() << " +/- "<< f2Svs1S->getError()<< endl <<endl;
      yAxisRap2S[1]=f2Svs1S->getVal(); yAxisRapError2S[1]=f2Svs1S->getError();
      yAxisRap3S[1]=f3Svs1S->getVal(); yAxisRapError3S[1]=f3Svs1S->getError();
      TGraphErrors *gr2 = new TGraphErrors(nBins,x,yAxisRap2S,ex,yAxisRapError2S);
      TGraphErrors *gr3 = new TGraphErrors(nBins,x,yAxisRap3S,ex,yAxisRapError3S);
      TF1 *Axis= new TF1("Axis","0",-2.4,2.4);
      TCanvas cRap; cRap.cd();
      gr2->SetName("Upsilon(2S) ratio");
      gr3->SetName("Upsilon(3S) ratio");
      gr2->SetTitle("Upsilon(2S) ratio");
      gr3->SetTitle("Upsilon(3S) ratio");
      gr2->SetMarkerColor(4);
      gr3->SetMarkerColor(kRed);
      gr2->SetMarkerStyle(21);
      gr3->SetMarkerStyle(21);
      Axis->GetXaxis()->SetTitle("Rapidity");
      Axis->GetYaxis()->SetTitle("2S/1S, 3S/1S ratios");
      Axis->Draw("");
      gr2->Draw("Psame");
      gr3->Draw("Psame");
      TLatex latex1;
      latex1.SetNDC();
      if (PbPb) {
	latex1.DrawLatex(0.46,1.-0.05*3,"CMS pPb  #sqrt{s} = 5 TeV");
	latex1.DrawLatex(0.5,1.-0.05*4.9,"L_{int} = 10.66 #nb^{-1}"); 
	latex1.DrawLatex(0.5,1.-0.05*6.2," |y_{coll}| < 1.93");
      }
      else {
	latex1.DrawLatex(0.46,1.-0.05*3,"CMS pp  #sqrt{s} = 2.76 TeV");
	latex1.DrawLatex(0.5,1.-0.05*4.9,"L_{int} = 231 nb^{-1}"); 
	latex1.DrawLatex(0.5,1.-0.05*6.2," |y| < 1.93");
      }
      latex1.DrawLatex(0.5,1.-0.05*7.5,"p_{T}^{#mu} > 4 GeV/c");
      //latex1.DrawLatex(0.7,1.-0.05*8.8,"data");
      latex1.DrawLatex(0.7,1.-0.05*9.8,"2S/1S");
      latex1.DrawLatex(0.7,1.-0.05*10.8,"3S/1S");
      //TMarker M(11.35, 425, 20);
      //M.SetMarkerSize(1);
      //M.DrawMarker(11.35, 425);
      
      TLine L1; 
      L1.SetLineWidth(2);
      L1.SetLineColor(4);
      L1.DrawLineNDC(1-0.62, 1. - 0.05*9.55,
		     1-0.68, 1. - 0.05*9.55);
      
      TLine L2; 
      L2.SetLineWidth(2);
      L2.SetLineColor(kRed);
      L2.DrawLineNDC(1-0.62, 1. - 0.05*10.55,
		     1-0.68, 1. - 0.05*10.55);
      if(PbPb);//cRap.SaveAs(figs_+"ratiosVSfullRapidity_pPb10nbBF.pdf"); 
      else ;//cRap.SaveAs(figs_+"ratiosVSfullRapidity_pp276.pdf");
      break;
    case 6 :
      x_276[0]=0.;ex_276[0]=0;
      cout << "no bin" << endl ;
      yAxisRap2S_276[0]=f2Svs1S->getVal(); yAxisRapError2S_276[0]=f2Svs1S->getError();
      yAxisRap3S_276[0]=f3Svs1S->getVal(); yAxisRapError3S_276[0]=f3Svs1S->getError();   
      cout << "i'm not plotting " << endl;
    break;
    case 7 :
      //    X[0]=0.;ex[0]=1.93;
      cout << "no bin" << endl ;
      cout << "i'm plotting " << endl;
      TGraphErrors *gr2 = new TGraphErrors(nBins,x_276,yAxisRap2S_276,ex_276,yAxisRapError2S_276);
      TGraphErrors *gr3 = new TGraphErrors(nBins,x_276,yAxisRap3S_276,ex_276,yAxisRapError3S_276);
      yAxisRap2S[0]=f2Svs1S->getVal(); yAxisRapError2S[0]=f2Svs1S->getError();
      yAxisRap3S[0]=f3Svs1S->getVal(); yAxisRapError3S[0]=f3Svs1S->getError();
      gr2->SetName("Upsilon(2S) ratio");
      gr3->SetName("Upsilon(3S) ratio");
      gr2->SetTitle("Upsilon(2S) ratio");
      gr3->SetTitle("Upsilon(3S) ratio");
      gr2->SetMarkerColor(kRed);
      gr3->SetMarkerColor(4);
      gr2->SetMarkerStyle(21);
      gr3->SetMarkerStyle(21);
      TGraphErrors *gr4 = new TGraphErrors(nBins,x,yAxisRap2S,ex,yAxisRapError2S);
      TGraphErrors *gr5 = new TGraphErrors(nBins,x,yAxisRap3S,ex,yAxisRapError3S);
      cout << "Do Overlay !" << endl;
      cout << "Oh Yeah !"<< endl ;
      TCanvas cRap; cRap.cd();
      TF1 *Axis = new TF1("Axis","0",-1.93,1.93);
      gr4->SetName("Upsilon(2S) ratio");
      gr5->SetName("Upsilon(3S) ratio");
      gr4->SetTitle("Upsilon(2S) ratio");
      gr5->SetTitle("Upsilon(3S) ratio");
      gr4->SetMarkerColor(kRed+3);
      gr5->SetMarkerColor(38);
      gr4->SetMarkerStyle(21);
      gr5->SetMarkerStyle(21);
      Axis->GetXaxis()->SetTitle("Rapidity");
      Axis->GetYaxis()->SetTitle("Y(nS) ratios (n=2,3)");
      Axis->GetYaxis()->SetRangeUser(0,1);
      Axis->Draw("");
      cout << "i'm overlaying now " << endl;
      Axis->Draw("");
      gr2->Draw("Psame");
      gr3->Draw("Psame");
      gr4->Draw("Psame");
      gr5->Draw("Psame");
       
      TLatex latex1;
      latex1.SetNDC();
      if (PbPb) {
      	//latex1.DrawLatex(0.46,1.-0.05*3,"CMS pPb  #sqrt{s} = 5 TeV");
 	latex1.DrawLatex(0.2,1.-0.05*2.3,"CMS pPb  #sqrt{s} = 5 TeV, PbPb #sqrt{s_{NN}} = 2.76 TeV");
	latex1.DrawLatex(0.2,1.-0.05*3.2,"L_{int} = 10.66 nb^{-1}, 150 #mub^{-1}"); 
	latex1.DrawLatex(0.2,1.-0.05*4.1," |y_{coll}| < 1.93");
      }
      else {
	latex1.DrawLatex(0.2,1.-0.05*3,"CMS Preliminary");
	latex1.DrawLatex(0.2,1.-0.05*4,"pp #sqrt{s} = 7 TeV, 2.76 TeV");
	latex1.DrawLatex(0.6,1.-0.05*4,"L_{int} = 5 fb^{-1}, 231 nb^{-1}"); 
	latex1.DrawLatex(0.6,1.-0.05*5," |y_{coll}| < 1.93");
      }
      latex1.DrawLatex(0.2,1.-0.05*5.,"p_{T}^{#Upsilon} > 4 GeV/c");
      //latex1.DrawLatex(0.7,1.-0.05*8.8,"data");
      if(PbPb){
      latex1.DrawLatex(0.2,1.-0.05*6,"2S/1S - 5 TeV");
      latex1.DrawLatex(0.2,1.-0.05*6.95,"3S/1S - 5 TeV");
      }
      else {
      latex1.DrawLatex(0.2,1.-0.05*6,"2S/1S - 2.76 TeV");
      latex1.DrawLatex(0.2,1.-0.05*6.95,"3S/1S - 2.76 TeV");	
      }
      //TMarker M(11.35, 425, 20);
      //M.SetMarkerSize(1);
      //M.DrawMarker(11.35, 425);
      
      TLine L1; 
      L1.SetLineWidth(2);
      L1.SetLineColor(kRed);
      L1.DrawLineNDC(0.12, 1. - 0.05*5.70,
		     0.18, 1. - 0.05*5.70);
      
      TLine L2; 
      L2.SetLineWidth(2);
      L2.SetLineColor(4);
      L2.DrawLineNDC(0.12, 1. - 0.05*6.65,
		     0.18, 1. - 0.05*6.65);
      TLatex latex2;
      latex2.SetNDC();
      latex2.DrawLatex(0.2,1.-0.05*7.9,"2S/1S - 7 TeV");
      latex2.DrawLatex(0.2,1.-0.05*8.85,"3S/1S - 7 TeV");
      //TMarker M(11.35, 425, 20);
      //M.SetMarkerSize(1);
      //M.DrawMarker(11.35, 425);
      
      TLine L4; 
      L4.SetLineWidth(2);
      L4.SetLineColor(kRed+3);
      L4.DrawLineNDC(0.12, 1. - 0.05*7.6,
		     0.18, 1. - 0.05*7.6);
      
      TLine L5; 
      L5.SetLineWidth(2);
      L5.SetLineColor(38);
      L5.DrawLineNDC(.12, 1. - 0.05*8.55,
		     0.18, 1. - 0.05*8.55);
      if(PbPb)cRap.SaveAs(figs_+figName_+"ratiosVSfullRapidity_pPb20_PbPb100nbBF.pdf"); 
      cout << "file to be saved..." << endl;
      if(!PbPb) cRap.SaveAs(figs_+figName_+"ratiosVSfullRapidity_ppOverlay_test.pdf");
      break;
    default: break;
    }
    
}
void rapbin_yields(int bin, RooRealVar *nsig1f, RooRealVar *nsig2f, RooRealVar *nsig3f){
  switch (bin)
    {
    case 0 :
      //RooArgSet *ratio2 = fit_2nd->getParameters(nsig2f);
      //     fit_2nd->Print();
      cout << " 2s = "<< nsig2f->getVal() << " +/- "<< nsig2f->getError()<< endl <<endl;
      x[0]=-1.2; ex[0]=0.73;
      yAxisRap2S[0]=nsig2f->getVal(); yAxisRapError2S[0]=nsig2f->getError();
      yAxisRap1S[0]=nsig1f->getVal(); yAxisRapError1S[0]=nsig1f->getError();
      yAxisRap3S[0]=nsig3f->getVal(); yAxisRapError3S[0]=nsig3f->getError();
      break;
    case 1 :
      //RooArgSet *ratio2 = fit_2nd->getParameters(nsig2f);
      //     fit_2nd->Print();
      x[1]=-.235;ex[1]=.235;
      cout << " 2s = "<< nsig2f->getVal() << " +/- "<< nsig2f->getError()<< endl <<endl;
      yAxisRap2S[1]=nsig2f->getVal(); yAxisRapError2S[1]=nsig2f->getError();
      yAxisRap1S[1]=nsig1f->getVal(); yAxisRapError1S[1]=nsig1f->getError();
      yAxisRap3S[1]=nsig3f->getVal(); yAxisRapError3S[1]=nsig3f->getError();
      break;
    case 2 :
      //RooArgSet *ratio2 = fit_2nd->getParameters(nsig2f);
      //     fit_2nd->Print();
      x[2]=0.235;ex[2]=.235;
      cout << " 2s = "<< nsig2f->getVal() << " +/- "<< nsig2f->getError()<< endl <<endl;
      yAxisRap2S[2]=nsig2f->getVal(); yAxisRapError2S[2]=nsig2f->getError();
      yAxisRap1S[2]=nsig1f->getVal(); yAxisRapError1S[2]=nsig1f->getError();
      yAxisRap3S[2]=nsig3f->getVal(); yAxisRapError3S[2]=nsig3f->getError();
      break;
    case 3 :
      x[3]=1.2;ex[3]=.73;
      //RooArgSet *ratio2 = fit_2nd->getParameters(nsig2f);
      //     fit_2nd->Print();
      cout << " 2s = "<< nsig2f->getVal() << " +/- "<< nsig2f->getError()<< endl <<endl;
      yAxisRap2S[3]=nsig2f->getVal(); yAxisRapError2S[3]=nsig2f->getError();
      yAxisRap1S[3]=nsig1f->getVal(); yAxisRapError1S[3]=nsig1f->getError();
      yAxisRap3S[3]=nsig3f->getVal(); yAxisRapError3S[3]=nsig3f->getError();
      /*  TGraphErrors *gr1 = new TGraphErrors(nBins,x,yAxisRap1S,ex,yAxisRapError1S);
      TGraphErrors *gr2 = new TGraphErrors(nBins,x,yAxisRap2S,ex,yAxisRapError2S);
      TGraphErrors *gr3 = new TGraphErrors(nBins,x,yAxisRap3S,ex,yAxisRapError3S);
      TF1 *Axis= new TF1("Axis","0",-2.4,2.4);
      TCanvas cRap; cRap.cd();
      	gr1->SetName("Upsilon(1S) yield");
	gr2->SetName("Upsilon(2S) yield");
	gr3->SetName("Upsilon(3S) yield");
	gr1->SetTitle("Upsilon(1S) yield");
	gr2->SetTitle("Upsilon(2S) yield");
	gr3->SetTitle("Upsilon(3S) yield");
	gr1->SetMarkerColor(2);
	gr2->SetMarkerColor(4);
	gr3->SetMarkerColor(3);
	gr1->SetMarkerStyle(21);
	gr2->SetMarkerStyle(21);
	gr3->SetMarkerStyle(21);
	Axis->GetXaxis()->SetTitle("Rapidity");
	Axis->GetYaxis()->SetTitle("dN/dy");
	Axis->Draw("");
	gr1->Draw("Psame");
	gr2->Draw("Psame");
	gr3->Draw("Psame");
	if(PbPb) cRap.SaveAs(figs_+"yieldsVSrapidity_pPb.pdf");
	else cRap.SaveAs(figs_+"yieldsVSrapidity_pp.pdf");*/
	break; 
    case 4 :
      x[0]=-.965;ex[0]=.965;
      //RooArgSet *ratio2 = fit_2nd->getParameters(nsig2f);
      //     fit_2nd->Print();
      cout << " 2s = "<< nsig2f->getVal() << " +/- "<< nsig2f->getError()<<"bweuahaaaaaaa ! ! !" <<endl <<endl;
         cout << " 2s = "<< nsig2f->getVal() << " +/- "<< nsig2f->getError()<< endl <<endl;
      cout << " 3s = "<< nsig3f->getVal() << " +/- "<< nsig3f->getError()<< endl <<endl;
      yAxisRap1S[0]=nsig1f->getVal(); yAxisRapError1S[0]=nsig1f->getError();
      yAxisRap2S[0]=nsig2f->getVal(); yAxisRapError2S[0]=nsig2f->getError();
      yAxisRap3S[0]=nsig3f->getVal(); yAxisRapError3S[0]=nsig3f->getError();
      break; 
    case 5 :
      x[1]=.965;ex[1]=.965;
      //RooArgSet *ratio2 = fit_2nd->getParameters(nsig2f);
      //     fit_2nd->Print();
         cout << " 1s = "<< nsig1f->getVal() << " +/- "<< nsig1f->getError()<<"bweuahaaaaaaa ! ! !" << endl <<endl;
       cout << " 2s = "<< nsig2f->getVal() << " +/- "<< nsig2f->getError()<< endl <<endl;
      cout << " 3s = "<< nsig3f->getVal() << " +/- "<< nsig3f->getError()<< endl <<endl;
      yAxisRap1S[1]=nsig1f->getVal(); yAxisRapError1S[1]=nsig1f->getError();
      yAxisRap2S[1]=nsig2f->getVal(); yAxisRapError2S[1]=nsig2f->getError();
      yAxisRap3S[1]=nsig3f->getVal(); yAxisRapError3S[1]=nsig3f->getError();
      TGraphErrors *gr1 = new TGraphErrors(nBins,x,yAxisRap1S,ex,yAxisRapError1S);
      TGraphErrors *gr2 = new TGraphErrors(nBins,x,yAxisRap2S,ex,yAxisRapError2S);
      TGraphErrors *gr3 = new TGraphErrors(nBins,x,yAxisRap3S,ex,yAxisRapError3S);
      cout << "before TF1 def" << endl;
      TF1 *Axis= new TF1("Axis","0",-1.93,1.93);
      cout << "after TF1 def"<< endl ;
        TCanvas cRap; cRap.cd();
	gr1->SetName("Upsilon(1S) yield");
	gr2->SetName("Upsilon(2S) yield");
	gr3->SetName("Upsilon(3S) yield");
	gr2->SetTitle("Upsilon(2S) yield");
	gr3->SetTitle("Upsilon(3S) yield");
	gr1->SetMarkerColor(2);
	gr2->SetMarkerColor(4);
	gr3->SetMarkerColor(3);
	gr1->SetMarkerStyle(21);
	gr2->SetMarkerStyle(21);
	gr3->SetMarkerStyle(21);
	Axis->GetXaxis()->SetTitle("Rapidity");
	Axis->GetYaxis()->SetTitle("Y(nS) yields (n<4)");
	Axis->GetYaxis()->SetRangeUser(0,800);
	Axis->Draw("");
	gr1->Draw("Psame");
	cout << "i'm plotting " << endl;
	gr2->Draw("Psame");
	gr3->Draw("Psame");
	cRap.SaveAs(figs_+"yieldsVSrapidity_pp7_276.pdf"); 
      break;
      
    default: break;
    }
}

//calculate the confidence interval with RooStats
pair<double, double> ConfidencInterval(float CI, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf)  {  
	ProfileLikelihoodCalculator pl(*data,*pdf,*fnll);
	pl.SetConfidenceLevel(CI); // 95, 68, whatever
	LikelihoodInterval* inetrval = pl.GetInterval();
	LikelihoodIntervalPlot plot(interval);
	TCanvas c4; c4.cd(); 
	plot.Draw();
	TString intrvlName = fnll->GetTitle();
	c4.SaveAs(figs_+"nll_"+intrvlName+"_pt4.gif");
	c4.SaveAs(figs_+"nll_"+intrvlName+"_pt4.pdf");
	// print out the iterval on the Parameter of Interest
	cout <<endl<< CI <<"\% interval on " <<fnll->GetName()<<" is : ["<<
		interval->LowerLimit(*fnll) << ", "<<
		interval->UpperLimit(*fnll) << "] "<<endl;
	pair<double, double> CnfdncIntrvl;
	CnfdncIntrvl.first  = interval->LowerLimit(*fnll);
	CnfdncIntrvl.second = interval->UpperLimit(*fnll);
	return CnfdncIntrvl;
}


void suppression(RooRealVar *nsig1f, RooFitResult *fit_2nd, int nfloatpars, RooAbsPdf *pdf, RooAbsPdf *pdf_pp, RooDataSet *data, RooPlot *frame){
	//the layout fit with pp ratio pdf
	//nsig1f->setVal(1380);
	TCanvas c1; c1.cd();
	//pdf->plotOn(frame,Name("thePdf"));
	if (RAA) nsig1f->setVal(2464);//(2352);
	//nsig1f->setConstant(kTRUE);  //fix the 1S yield to the value from the default fit 
	//RooFitResult* fit_3rd = pdf_pp->fitTo(*data,Save(kTRUE),Extended(kTRUE));
	pdf_pp->plotOn(frame,Name("thePdf_pp"),LineStyle(kDashed),LineColor(kRed));
	//pdf_pp->plotOn(frame,Components("pdf_combinedbkgd"),LineColor(kRed),LineStyle(kDashed));
	/*
	   RooArgSet * pars_pp = pdf_pp->getParameters(data);
	   int nfloatpars_pp = pars_pp->selectByAttrib("Constant",kFALSE)->getSize();
	   float nll2 = fit_2nd->minNll();
	   float nll3 = fit_3rd->minNll();
	   double chi = 2.0*(nll3-nll2);
	   int deltaNDOF = nfloatpars-nfloatpars_pp;
	   double cl = 1.0 - TMath::Prob(chi,deltaNDOF)/2;
	   double Significance = -RooStats::PValueToSignificance(cl);
	   outfile<<endl<<"*************************************"<<endl;
	   outfile<<"the delta of S is         : "<<nll3-nll2<<endl;
	   outfile<<"the delta of ndof is      : "<<deltaNDOF<<endl;
	   outfile<<"the C.L. is               : "<<cl<<endl;
	   outfile<<"the significance level is : "<<Significance<<" (one side)"<<endl;
	   outfile<<"*************************************"<<endl<<endl;
	   outfile<<"the C.L. is : "<<cl<<",  the significance level is : "<<Significance<<" (one side)"<<endl<<endl;
	 */	
	//draw and save plots
	data->plotOn(frame,MarkerSize(0.8));
	frame->Draw();
	if(!plotpars) {
		if (RAA) float a=1160;
		else float a=730;
		float delta=a/62*5;
		TLatex latex1;
		latex1.DrawLatex(9.9,a, "CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV");
		a=a-delta/50.0;
		TLatex latex2;
		latex2.DrawLatex(9.9,a-1*delta,"Cent. 0-100%, |y| < 2.4");
		TLatex latex3;
		//^F^F^Flatex3.DrawLatex(9.9,a-2*delta,"p_{T} < 20 GeV/c");
		TLatex latex4;
		latex4.DrawLatex(9.9,a-2*delta, "L_{int} = 150 #mub^{-1}");
		TLatex latex5;
		latex5.DrawLatex(9.9,a-3*delta, "p_{T}^{#mu} > 4 GeV/c");

		if (!RAA) {
			TLatex latex9;
			//  latex9.DrawLatex(7.4,a, "Preliminary");

			latex9.DrawLatex(8,a, "data");  
			TMarker M(7.55,delta/5.0+a,20);
			M.Draw();
			//TLine L0(7.3,delta/4.0+a,7.8,delta/4.0+a);
			//L0.Draw();
			//TLine L1(7.55,a,7.55,a);
			//L1.Draw();

			TLatex latex10;
			latex10.DrawLatex(8,a-1*delta, "PbPb fit");
			/*	TLine L1(7.3,delta/5.0+a-1*delta,7.8,delta/5.0+a-1*delta);
			L1.SetLineColor(kBlue);
			L1.SetLineWidth(2);
			L1.Draw();*/

			TLatex latex11;
			latex11.DrawLatex(8,a-2*delta, "pp shape");
			TLine L2(7.3,delta/5.0+a-2*delta,7.8,delta/5.0+a-2*delta);
			L2.SetLineColor(kRed);
			L2.SetLineStyle(kDashed);
			L2.SetLineWidth(2);
			L2.Draw();

		}

		else {
			TLatex latex9;
			//  latex9.DrawLatex(7.4,a, "Preliminary");

			latex9.DrawLatex(8+3.5,a-5*delta, "data");
			TMarker M(7.55+3.5,delta/5.0+a-5*delta,20);
			M.Draw();
			//TLine L0(7.3,delta/4.0+a,7.8,delta/4.0+a);
			//L0.Draw();
			//TLine L1(7.55,a,7.55,a);
			//L1.Draw();

			TLatex latex10;
			latex10.DrawLatex(8+3.5,a-6*delta, "PbPb fit");
			TLine L1(7.3+3.5,delta/5.0+a-6*delta,7.8+3.5,delta/5.0+a-6*delta);
			L1.SetLineColor(kBlue);
			L1.SetLineWidth(2);
			L1.Draw();

			TLatex latex11;
			latex11.DrawLatex(8+3.5,a-7*delta, "pp shape");
			latex11.SetTextSize(0.035);
			latex11.DrawLatex(11.7,a-7.75*delta, "(R_{AA} scaled)");
			TLine L2(7.3+3.5,delta/5.0+a-7*delta,7.8+3.5,delta/5.0+a-7*delta);
			L2.SetLineColor(kRed);
			L2.SetLineStyle(kDashed);
			L2.SetLineWidth(2);
			L2.Draw();
		}

	}

	//	c1.SaveAs(figs_+"overlay1_"+figName_+paramOn_+suffix_+".gif");
	c1.SaveAs(figs_+"overlay1_"+figName_+paramOn_+suffix_+".pdf");
	nsig1f->setConstant(kFALSE);
}


void setUpLimits(float xMin, float xMax, RooAbsPdf *pdf, RooDataSet *data, RooRealVar *param, float baseNll){
	//setting up upper limits
	TCanvas c2; c2.cd();
	int totalBins = 100;
	float BinWidth = (xMax - xMin)/totalBins;
	TH1F *h1 = new TH1F("h1","h1",totalBins,xMin,xMax);
	h1->GetXaxis()->SetTitle(param->getTitle());
	//h1->GetXaxis()->SetRangeUser(-0.1,0.9);
	h1->GetYaxis()->SetTitle("Maximum likelihood");
	TH1F *h2 = new TH1F("h2","h2",totalBins,xMin,xMax);
	gStyle->SetOptStat(kFALSE);
	RooFitResult* fit_nll;
	double MinNll, L, cl;
	for (int i=1; i<=totalBins; i++) {
		param->setVal(xMin + BinWidth*(i-1));
		param->setConstant(kTRUE);
		fit_nll = pdf->fitTo(*data,Save(kTRUE));
		MinNll = fit_nll->minNll()-baseNll;
		L = TMath::Exp(-MinNll);
		cout<<"x = "<<param->getVal()<<", MinNll = "<<MinNll<<", L = "<<L<<endl;
		h1->SetBinContent(i,L);
	}
	for (int i=1; i<=totalBins; i++) {
		cout<<"bin "<<i<<" = "<< h1->GetBinContent(i)<<endl;
	}
	h1->Draw();
	cout<<endl<<"integral = "<<h1->Integral(1,totalBins)<<endl<<endl;
	h1->Scale(1.0/h1->Integral(1,totalBins));
	h1->Draw();
	c2.Update();

	//convoluted with systematic gaussian
	TH1F *h3 = new TH1F("h3","h3",totalBins,xMin,xMax);
	TH1F *h4 = new TH1F("h4","h4",totalBins,xMin,xMax);
	for (int i=1; i<=totalBins; i++) {
		float gausmean = xMin + BinWidth*(i-1);
		cout<<"gausmean = "<<gausmean;
		float gaussigma = 0.00001;
		TH1F *h_syst = new TH1F("h_syst","syst histogram",totalBins,xMin,xMax);
		for (Int_t k=0;k<10000;k++) {h_syst->Fill(gRandom->Gaus(gausmean,gaussigma));}
		double new_nll = 0;
		for (int j=1; j<=totalBins; j++){
			new_nll += (h_syst->GetBinContent(j) * h1->GetBinContent(j));
		}
		cout<<", new value after convolution in bin "<<i<<" = "<<new_nll<<endl;
		h3->SetBinContent(i,new_nll);
		delete h_syst;
	}
	cout<<endl<<"integral = "<<h3->Integral(1,totalBins)<<endl<<endl;
	h3->Scale(1.0/h3->Integral(1,totalBins));
	h3->SetLineColor(kMagenta);
	h3->SetLineWidth(2);
	h3->SetLineStyle(2);
	h3->Draw("same");

	//find out the upper limit
	c2.Update();
	float UpperLimit;
	double CI = 0.842; //0.842; //0.158;
	for (int i=1; i<=totalBins; i++) {
		cl = h1->Integral(1,i);
		cout<<"x = "<< xMin + BinWidth*(i-1) << ", y = " << h1->GetBinContent(i) <<", cl = "<<cl<<endl;
		if (cl<CI) UpperLimit = xMin + BinWidth*(i-1+1);
		h4->SetBinContent(i,cl);
	}   
	float rightmax = 1.1*h4->GetMaximum();
	float y_scale = gPad->GetUymax()/rightmax;
	cout<<gPad->GetUymax()<<" rightmax = "<<rightmax<<", scale = "<<y_scale<<endl;
	h4->SetLineWidth(2);
	h4->SetLineColor(kRed);
	h4->Scale(y_scale);
	h4->Draw("same");
	TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
	axis->SetLineColor(kRed); axis->SetLabelColor(kRed); 
	axis->SetTitle("confidence level"); axis->SetTitleColor(kRed);
	axis->Draw();

	TLine L1(xMin,y_scale*CI,xMax,y_scale*CI);
	L1.SetLineColor(kBlue);
	c2.SaveAs(figs_+"UpperLimits_"+paramName+"_hi_syst_Extend1.pdf");
	delete h1; delete h2; delete h3; delete h4;
	param->setConstant(kFALSE);
}


Double_t RooPlot::mychiSquare(const char* curvename, const char* histname, Int_t nFitParam, bool chisqRange, int startbin, int stopbin) const
{
	// Calculate and return reduced chi-squared of curve with given name with respect
	// to histogram with given name. If nFitParam is non-zero, it is used to reduce the
	// number of degrees of freedom for a chi^2 for a curve that was fitted to the
	// data with that number of floating parameters

	// Find curve object
	RooCurve* curve = (RooCurve*) findObject(curvename,RooCurve::Class()) ;
	if (!curve) {
		coutE(InputArguments) << "RooPlot::chiSquare(" << GetName() << ") cannot find curve" << endl ;
		return -1. ;
	}

	// Find histogram object
	RooHist* hist = (RooHist*) findObject(histname,RooHist::Class()) ;
	if (!hist) {
		coutE(InputArguments) << "RooPlot::chiSquare(" << GetName() << ") cannot find histogram" << endl ;
		return -1. ;
	}
	return curve->mychiSquare(*hist,nFitParam,chisqRange,startbin,stopbin) ;
}


Double_t RooCurve::mychiSquare(const RooHist& hist, Int_t nFitParam, bool chisqRange, int startbin, int stopbin) const 
{
	// Calculate the chi^2/NDOF of this curve with respect to the histogram
	// 'hist' accounting nFitParam floating parameters in case the curve
	// was the result of a fit

	if (chisqRange) {
		int i = startbin;
		int np = stopbin;
	}
	else {
		Int_t i=0;
		int np = hist.GetN() ;
	}

	Double_t x,y,eyl,eyh,exl,exh ;

	// Find starting and ending bin of histogram based on range of RooCurve
	Double_t xstart,xstop ;

#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
	GetPoint(0,xstart,y) ;
	GetPoint(GetN()-1,xstop,y) ;
#else
	const_cast<RooCurve*>(this)->GetPoint(0,xstart,y) ;
	const_cast<RooCurve*>(this)->GetPoint(GetN()-1,xstop,y) ;
#endif

	Int_t nbin(0) ;

	Double_t chisq(0) ;
	for (i ; i<np ; i++) {   

		// Retrieve histogram contents
		((RooHist&)hist).GetPoint(i,x,y) ;

		// Check if point is in range of curve
		if (x<xstart || x>xstop) continue ;

		nbin++ ;
		eyl = hist.GetEYlow()[i] ;
		eyh = hist.GetEYhigh()[i] ;
		exl = hist.GetEXlow()[i] ;
		exh = hist.GetEXhigh()[i] ;

		// Integrate function over this bin
		Double_t avg = average(x-exl,x+exh) ;

		// Add pull^2 to chisq
		if (y!=0) {      
			Double_t pull = (y>avg) ? ((y-avg)/eyl) : ((y-avg)/eyh) ;
			chisq += pull*pull ;
		}
	}

	// Return chisq/nDOF 
	return chisq / (nbin-nFitParam) ;
}

// static const unsigned int nBins = 1;
// bool pp7=1;
// bool pp276=0;
// bool pA=0;
// bool PbPb =0; //Input data sample.  1: PbPb data;   0: pp data
// bool MC=0;
// bool fitMB =0;        //1: fit the the Minbias sample(0-100%);   0: fit the each centrality bin
// double width_ = 0.0782;//new resolution for 2011 HI data
// bool plotpars = 1;     //1: plot parameters;   0: plot CMS label
// bool doMinos = 0;      //kFALSE;
// int bkgdModel = 4;     //Background Model.  1: LS erf*exp + polynomial; 2: LS RookeyPdf + polynomial; 3: erf*exp; 4: polynomial; 5: erf*exp+polynomial; 6: Implement simple shape
// bool plotLikeSign = 0; //1: plot likesign or trkRot data points and fit lines; 0: hide likesign or trkRot
// bool TRKROT =0;       //0: use likesign;   1: use track rotation
// bool LS_constrain = 0; //1: use constrain method
// bool PR_plot = 0;      //1: draw the PR plot
// bool RAA = 0;          //1: raa PR plot;   0: double ratio PR plot
// bool PbPb_23=1; // 1. enable 2S+3S/1S ratios
// RooRealVar *f2Svs1S_pp = new RooRealVar("N_{2S}/N_{1S}pp","Y(3S)/Y(1S) yields pp ratio",0.4,-1,5);
// RooRealVar *f3Svs1S_pp = new RooRealVar("N_{3S}/N_{1S}pp","Y(2S)/Y(1S) yields pp ratio",0.2,-1,5);

// double yAxisRap2S[nBins];
// double yAxisRap3S[nBins];
// double yAxisRap2S_276[nBins];
// double yAxisRap3S_276[nBins];
// double yAxisRap1S[nBins];
// double yAxisRapError1S[nBins];
// double yAxisRapError2S[nBins];
// double yAxisRapError3S[nBins];
// double yAxisRapError2S_276[nBins];
// double yAxisRapError3S_276[nBins];
// double x[nBins];
// double ex[nBins];
// double x_276[nBins];
// double ex_276[nBins];
