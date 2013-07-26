#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooKeysPdf.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

// Root stuff
#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"

#include "TMath.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"

// miscellaneous  
#include <fstream>
#include <iostream>

double mass_l =  7.0;
double mass_h = 14.0;
double binw   = 0.1;    //bin width of the histogram

const int nData   = 8;
const char* choseSampleLegend[nData] = {"",
					"pp #sqrt{s} = 7 TeV",
					"pp #sqrt{s} = 7 TeV",
					"PbPb #sqrt{s_{NN}} = 2.76 TeV (Regit)",
					"PbPb #sqrt{s_{NN}} = 2.76 TeV (pp reco TT)",
					"PbPb #sqrt{s_{NN}} = 2.76 TeV (pp reco GG)",
					"PbPb #sqrt{s_{NN}} = 2.76 TeV (HI reco GG)",
					"pp #sqrt{s} = 2.76 TeV"};

const char* choseSampleLumi[nData] = {"",
				      "L_{int} = 621 pb^{-1}",
				      "L_{int} = 966 pb^{-1}",
				      "L_{int} = 150 \mub^{-1}",
				      "L_{int} = xxx \mub^{-1}",
				      "L_{int} = xxx \mub^{-1}",
				      "L_{int} = 150 \mub^{-1}",
			      "L_{int} = 5.41 pb^{-1}"
};

// -------- make some local picks
bool doMinos      = 0;     //kFALSE;


using namespace RooFit;
using namespace RooStats;
pair<double, double> ConfidencInterval(float, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf);

void fitUpsilonYields(int choseSample    = 3, //Input data sample.  1: pp@7TeV data ; 2: pp@2.76TeV; 3: PbPb@276
		      int choseFitParams = 2, //0: (1s, 2s, 3s) 1: (1s, 2s/1s; 3s/1s); 2: (1S, (2s+3s)/1s)
		      int bkgdModel      = 3, //1:LS erf*exp + pol2; 2:LS RookeyPdf + pol2; 3:erf*exp; 4:pol2; 5:erf*exp+pol2 6:pol3
		      int fixFSR         = 1,
		      int fixSigma1      = 0, // 0 free; 1: fix
		      int centralityMin = 20,
		      int centralityMax = 40,
		      float muonEtaMin  = -1.,
		      float muonEtaMax  = 1, 
		      float dimuYMin    = -1., 
		      float dimuYMax    = 1.,
		      double muonpTcut  = 4, //single muon pT cut
		      bool plotBkg      = 0, //0: hide LS or trkRot; 1: plot LS or trkRot data points and fit lines;
		      bool doTrkRot     = 0, //0: use LS;   1: use track rotation
		      bool doConstrainFit   = 0,  //1: use constrain method
		      int useRef            = 1, // # 0 none 1: MC, 2: MB
		      bool plotpars         = 1, //1: plot parameters;   0: plot CMS label
		      const char* choseSampleCase = "regit", //
		      const char* outFigsDir      = "pdfOutput/",// figs dir outfile location
		      TString outDatsDir          = "txtOutput",// dats dir outfile location
		      const char* outFilePrefix   = "ratios", // yields, ratios
		      bool narrowMass             = false,
		      float vProbPick = 0.01
			 )
{
  // gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gROOT->Macro("cm/logon.C+");

  // input file
  TString finput;
  switch (choseSample) 
    {
   
  // when no ref is specified, then do free fit
  if(useRef==0)// no reference
    {
      fixSigma1 = 0;
      fixFSR    = 0;
    }
  //
  double muonEtaCut_min = muonEtaMin;
  double muonEtaCut_max = muonEtaMax; 
  // kinematic cuts

  cout << "Muon Eta Min is: "<< muonEtaMin << endl;
 case 1: // pp @ 7TeV - dimuon0-v1
      finput   = "../dimuonTree_upsiMiniTree_pp7tev_dimu0v1_Runa_trigBit1_allTriggers0_pt4.root";
      mass_l = 8.5;
      mass_h = 11.5;
      binw=0.05;
      break;
    case 2://pp @ 7TeV - non dimuon0-v1 (GG)
      finput   = "../dimuonTree_upsiMiniTree_pp7tev_glbglb_Runa_trigBit1_allTriggers0_pt4.root";
      mass_l = 8.5;
      mass_h = 11.5;
      binw=0.05;
      break;
    case 3://PbPb @ 2.76TeV regit
      //finput   = "../dimuonTree_upsiMiniTree_aa276tev_regitreco_glbglb_Runa_trigBit1_allTriggers0_pt4.root"; // cent 0-40 "cm"
      // finput = "../dimuonTree_upsiMiniTree_AA2p76tev_ptmu3_july09_Run2011-2011_trigBit1_allTriggers0.root";// cent0-40 "td"
       finput = "../dimuonTree_upsiMiniTree_AA276tevC0100_regit_ptmu4_Run210498-211631_trigBit1_allTriggers0.root"; // cent0-40 "nf"
      break;
    case 4://PbPb @ 2.76TeV pp reco trk-trk
      // finput   = "../dimuonTree_upsiMiniTree_aa276tev_50100ppreco__Runa_trigBit1_allTriggers0_pt4.root";
      // break;
      finput   = "../dimuonTree_upsiMiniTree_aa276tevC50100_ppofficial_trktrk_ptmu4_Runa_trigBit1_allTriggers0.root";
      break;
    case 5://PbPb @ 2.76TeV pp reco glbglb
      finput   = "../dimuonTree_upsiMiniTree_aa276tev_50100ppreco_glbglb_Runa_trigBit1_allTriggers0_pt4.root";
      break;
    case 6://PbPb @ 2.76TeV Zhen's tree glbglb
      //finput   = "../dimuonTree_upsiMiniTree_AA276tevC0100_regit_ptmu4_Run210498-211631_trigBit1_allTriggers0.root";
      finput   = "../dimuonTree_HI2011_fulldataset_trkRot.root";
      break;
    case 7://pp @ 2.76TeV
      finput   = "../dimuonTree_upsiMiniTree_pp276tev_5p41_Run211739-211831_trigBit1_allTriggers0_pt4.root";
      break;
    default:
      cout<<"You don't know what you are doing! Pick one of the available datasets in the choseSampleCases[] array"<<endl;
      break;
   }
  double upsYCut_min    = dimuYMin;
  double upsYCut_max    = dimuYMax;
 
  int centrality_max = centralityMax; 
  int centrality_min = centralityMin; 

  TString cut_ap(Form("(%d<=Centrality && Centrality<%d) && (%.2f<muPlusEta && muPlusEta < %.2f) && (%.2f<muMinusEta && muMinusEta < %.2f) && (%.2f<upsRapidity && upsRapidity<%.2f) && (vProb > %.2f) ",centrality_min,centrality_max,muonEtaMin,muonEtaMax,muonEtaMin,muonEtaMax,upsYCut_min,upsYCut_max,vProbPick));
 
  int centMin = centrality_min;
  int centMax = centrality_max;
  if(choseSample==3 || choseSample==6) 
    {
      centMin = (int)(centrality_min*2.5);
      centMax = (int)(centrality_max*2.5);
    }
  // output setting files:
  TString figsDir(Form("%s",outFigsDir)); //output fig location
  
  // if param are on, this gets filled in the naming
  TString paramOn_("");
  
  TString figName_(Form("%s_%s_cent%d-%d_bkgModel%d_muonEta%.2f%.2f_dimuY%.2f%.2f_trkRot%d_constrain%d_fsr%d_sigma%d_ref%d",
			  outFilePrefix,choseSampleCase,centMin,centMax,
			  bkgdModel, muonEtaMin,muonEtaMax,upsYCut_min,upsYCut_max,doTrkRot,doConstrainFit,fixFSR,fixSigma1,useRef));

  // output file names
  if(narrowMass)
    {
      mass_l =  8.5;
      mass_h = 11.5;
      TString figName_(Form("%s_%s_cent%d%d_bkgModel%d_muonEta%.2f%.2f_dimuY%.2f%.2f_trkRot%d_constrain%d_%d_%d_ref%d_mass8p511p5",
			    outFilePrefix,choseSampleCase,centMin,centMax,bkgdModel,
			    muonEtaMin,muonEtaMax,upsYCut_min,upsYCut_max,doTrkRot,doConstrainFit,fixFSR,fixSigma1,useRef)); // output file names 
      
      //dimupT%.1f%.1f_
      //upsPtCutMin,upsPtCutMax
      //double upsPtCutMin   = 0, //for ups_pT bins
      //double upsPtCutMax   = 150, //for ups_pT bins
      //latex1.DrawLatex(0.15,1.-0.05*5.5,Form("p_{T}^{#Upsilon} > %1.f",upsPt));
    } 
  figName_.ReplaceAll("-","M");
  figName_.ReplaceAll(".","");

  cout<<"Fitting: y["<< upsYCut_min <<","<<upsYCut_max<<"] and centrality ["<<centrality_min<<","<<centrality_max<<"]!!!!"<<endl;
  cout << "oniafitter processing"
       << "\n\tInput:  \t" << finput
       << "\n\tOutput: \t" << figName_
       << endl;
  // -----------

  TFile f(finput,"read");

  TTree* theTree       = (TTree*)gROOT->FindObject("UpsilonTree"); // OS --- all mass
  TTree* allsignTree   = (TTree*)gROOT->FindObject("UpsilonTree_allsign");//all sign and all mass
  TTree* trkRotTree=0;
  if (doTrkRot) trkRotTree = (TTree*)gROOT->FindObject("UpsilonTree_trkRot");
  
  RooRealVar* mass       = new RooRealVar("invariantMass","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}");
  RooRealVar* upsPt      = new RooRealVar("upsPt","p_{T}(#Upsilon)",0,60,"GeV");
  RooRealVar* upsEta     = new RooRealVar("upsEta",  "upsEta"  ,-10,10);
  RooRealVar* upsRapidity= new RooRealVar("upsRapidity",  "upsRapidity",-2.4, 2.4);
  RooRealVar* vProb      = new RooRealVar("vProb",  "vProb"  ,0.01,1.00);
  RooRealVar* QQsign     = new RooRealVar("QQsign",  "QQsign"  ,-1,5);
  RooRealVar* Centrality = new RooRealVar("Centrality","Centrality",0,100);
   
  RooRealVar* muPlusPt   = new RooRealVar("muPlusPt","muPlusPt",muonpTcut,50);
  RooRealVar* muMinusPt  = new RooRealVar("muMinusPt","muMinusPt",muonpTcut,50);
  RooRealVar* muPlusEta  = new RooRealVar("muPlusEta","muPlusEta",  -2.4,2.4);
  RooRealVar* muMinusEta = new RooRealVar("muMinusEta","muMinusEta",-2.4,2.4);
  RooRealVar* runNumber  = new RooRealVar("runNb",  "runNb",210000, 222000);

  // *************************************************** importing
  //##### import unlike-sign data set
  RooDataSet *data0 = new RooDataSet("data0","data0",theTree,
				     RooArgSet(*mass,*upsRapidity,*muPlusPt,*muMinusPt,*muPlusEta,*muMinusEta,*Centrality,*vProb));

  RooDataSet *data0_ap =  ( RooDataSet*)data0->reduce(Cut(cut_ap));
  data0_ap->SetName("data0_ap");
  data0_ap->Print();

  data0_ap->Print();

  RooDataSet *data = ( RooDataSet*)data0_ap;
  data->SetName("data");
  data->Print();
 
  // ###### import like-sign data set
  RooDataSet *likesignData0    = new RooDataSet("likesignData0","likesignData0",allsignTree,
						RooArgSet(*mass,*upsRapidity,*muPlusPt,*muMinusPt,*muPlusEta,*muMinusEta,*Centrality,*QQsign,*vProb));
  RooDataSet *likesignData0_ap = ( RooDataSet*)likesignData0->reduce(Cut(cut_ap+" && QQsign != 0"));
  likesignData0_ap->SetName("likesignData0_ap");
  likesignData0_ap->Print();

  RooDataSet *likesignData = ( RooDataSet*)likesignData0_ap;
  likesignData->SetName("likesignData");
  likesignData->Print();
  
  //import track-rotation data set
  RooDataSet *TrkRotData0, *TrkRotData;
  if (doTrkRot) 
    {
      TrkRotData0 = new RooDataSet("TrkRotData_ap","TrkRotData_ap",trkRotTree,
				   RooArgSet(*mass,*upsRapidity,*muPlusPt,*muMinusPt,*muPlusEta,*muMinusEta,*Centrality,*QQsign,*vProb));

      RooDataSet *TrkRotData0_ap = ( RooDataSet*)TrkRotData0->reduce(Cut(cut_ap));
      TrkRotData0_ap->SetName("TrkRotData0_ap");
      TrkRotData0_ap->Print();
    
      TrkRotData = ( RooDataSet*) TrkRotData0_ap;
      TrkRotData->SetName("TrkRotData");
      TrkRotData->Print();
    }
  
  // *************************************************** signal PDF
  const double M1S = 9.46;   //upsilon 1S pgd mass value
  const double M2S = 10.023;  //upsilon 2S pgd mass value
  const double M3S = 10.355;  //upsilon 3S pgd mass value
  
  RooRealVar  *mean = new RooRealVar("mass1S","#Upsilon mean",M1S,M1S-0.1,M1S+0.1);
  RooConstVar *rat2 = new RooConstVar("rat2", "rat2", M2S/M1S);
  RooConstVar *rat3 = new RooConstVar("rat3", "rat3", M3S/M1S);

  // scale mean and resolution by mass ratio
  RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",RooArgList(*mean));
  RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0*@1", RooArgList(*mean,*rat2));
  RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0*@1", RooArgList(*mean,*rat3));

 //detector resolution ?? where is this coming from?
  RooRealVar    *sigma1  = new RooRealVar("sigma1","#sigma_{1S}",0.0623,0.01,0.9); //  MC 5tev 1S pol2 
  RooFormulaVar *sigma1S = new RooFormulaVar("sigma1S","@0"   ,RooArgList(*sigma1));
  RooFormulaVar *sigma2S = new RooFormulaVar("sigma2S","@0*@1",RooArgList(*sigma1,*rat2));
  RooFormulaVar *sigma3S = new RooFormulaVar("sigma3S","@0*@1",RooArgList(*sigma1,*rat3));
  
  /// to describe final state radiation tail on the left of the peaks
  RooRealVar *alpha  = new RooRealVar("alpha","tail shift",1.78,0.4,20);    // MC 5tev 1S pol2 
  RooRealVar *npow   = new RooRealVar("npow","power order",1.52,1,100);    // MC 5tev 1S pol2 

  // ratios fit paramters:1 (pp@7.tev), 2 (pp@2.76TeV), 3(PbPb 2.76TeV)
  //no. all values except last column come from AN2011-455-v9
  // col.1 = MonteCarlo
  // col.2 = PbPb data nominal fit
  // col.3 = MC from 2013 vineet
  double alpha_mc[7]   = {0, 1.67 , 0.98 , 1.626};//
  double npow_mc[7]    = {0, 2.3  , 2.3  , 2.3};// 
  double sigma1_mc[7]  = {0, 0.090, 0.078,0.06506};// 

  alpha->setVal(alpha_mc[useRef]);
  npow->setVal(npow_mc[useRef]);
  switch (fixFSR) // 0: free;  1: both fixed 2: alpha fixed 3: npow fixed
    {      
    case 0:// all free
      alpha->setConstant(false);
      npow->setConstant(false);
      break;
    case 1:// fix both to MB or MC values
      alpha->setConstant(true);
      npow ->setConstant(true);
      break;
    case 2: // fix alpha
      alpha->setConstant(true);
      npow ->setConstant(false);
      break;
    case 3:// npow fix
      alpha->setConstant(false);
      npow ->setConstant(true);
      break;
    default:
      cout<<"Donno this choice! Pick somehting for FSR parameters that I know"<<endl;
      break;
    }

  switch (fixSigma1)
    {
    case 0: // free
      sigma1->setConstant(false);
      break;
    case 1: // fix to MB or MC values; 
      sigma1->setVal(sigma1_mc[useRef]);
      sigma1->setConstant(true);
      break;
    default:
      cout<<"Donno this choice! Pick somehting for sigma1! "<<endl;
      break;
    }


  // relative fraction of the two Gaussians components for each CB
  RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.3,0.,1.);
  sigmaFraction->setVal(0);
  sigmaFraction->setConstant(kTRUE);
  
  /// Upsilon 1S
  RooCBShape  *cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
 					   *mass,*mean1S,*sigma1,*alpha,*npow);
  RooCBShape  *cb1S_2    = new RooCBShape ("cb1S_2", "FSR cb 1s",
 					   *mass,*mean1S,*sigma1S,*alpha,*npow);
  RooAddPdf      *sig1S  = new RooAddPdf  ("sig1S","1S mass pdf",
 					   RooArgList(*cb1S_1,*cb1S_2),*sigmaFraction);
  // sig1S is then jsut cb1S2: c*pdf_1+(1-c)*pdf_2
  
  /// Upsilon 2S
  RooCBShape  *cb2S_1    = new RooCBShape ("cb2S_1", "FSR cb 2s", 
 					   *mass,*mean2S,*sigma1,*alpha,*npow); 
  RooCBShape  *cb2S_2    = new RooCBShape ("cb2S_2", "FSR cb 2s", 
 					   *mass,*mean2S,*sigma2S,*alpha,*npow); 
  RooAddPdf      *sig2S  = new RooAddPdf  ("sig2S","2S mass pdf",
 					   RooArgList(*cb2S_1,*cb2S_2),*sigmaFraction);
  
  /// Upsilon 3S
  RooCBShape  *cb3S_1    = new RooCBShape ("cb3S_1", "FSR cb 3s", 
 					   *mass,*mean3S,*sigma1,*alpha,*npow); 
  RooCBShape  *cb3S_2    = new RooCBShape ("cb3S_2", "FSR cb 3s", 
 					   *mass,*mean3S,*sigma3S,*alpha,*npow); 
  RooAddPdf      *sig3S  = new RooAddPdf  ("sig3S","3S mass pdf",
 					   RooArgList(*cb3S_1,*cb3S_2),*sigmaFraction); // = cb3S1*sigmaFrac + cb3S2*(1-sigmaFrac)
  
  // *************************************************** free param in the fit
  int nt = 100000;
  RooRealVar *nsig1f   = new RooRealVar("N_{#Upsilon(1S)}","nsig1S",nt*0.25,0,10*nt);
  switch (choseFitParams)
    {
    case 0://use the YIELDs of 2S and 3S as free parameters
      RooRealVar *nsig2f  = new RooRealVar("N2S","nsig2S",   nt*0.25,-1*nt,10*nt);
      RooRealVar *nsig3f  = new RooRealVar("N3S","nsig3S",   nt*0.25,-1*nt,10*nt);
      break;
    case 1:  //use the RATIOs of 2S and 3S as free parameters
      RooRealVar *f2Svs1S   = new RooRealVar("R_{#frac{2S}{1S}}","f2Svs1S",0.26,-0.1,1.0);
      RooRealVar *f3Svs1S   = new RooRealVar("R_{#frac{3S}{1S}}","f3Svs1S",0.13,-0.1,1.0);
      RooFormulaVar *nsig2f = new RooFormulaVar("N2S","@0*@1", RooArgList(*nsig1f,*f2Svs1S));
      RooFormulaVar *nsig3f = new RooFormulaVar("N3S","@0*@1", RooArgList(*nsig1f,*f3Svs1S));
      f2Svs1S->setConstant(kFALSE);
      f3Svs1S->setConstant(kFALSE);
      break;
    case 2:// do (2s+3s)/1s
      RooRealVar *f2Svs1S   = new RooRealVar("R_{#frac{2S}{1S}}","f2Svs1S",0.26,-0.1,1.0);
      RooRealVar *f23vs1S   = new RooRealVar("R_{#frac{2S+3S}{1S}}","f23vs1S",0.45,-0.1,1);
      RooFormulaVar *nsig2f = new RooFormulaVar("N2S","@0*@1", RooArgList(*nsig1f,*f2Svs1S));
      RooFormulaVar *nsig3f = new RooFormulaVar("N3S","@0*@2-@0*@1", 
 						RooArgList(*nsig1f,*f2Svs1S,*f23vs1S));
      break;
    case 3://do 2s/1s, 3s/1s, 3s/2s
      RooRealVar *f2Svs1S   = new RooRealVar("R_{#frac{2S}{1S}}","f2Svs1S",0.26,-0.1,1.0);
      RooRealVar *f3Svs1S   = new RooRealVar("R_{#frac{3S}{1S}}","f3Svs1S",0.13,-0.1,1.0);
      RooFormulaVar *nsig2f = new RooFormulaVar("N2S","@0*@1", RooArgList(*nsig1f,*f2Svs1S));
      //   RooFormulaVar *nsig3f = new RooFormulaVar("N3S","@0*@1", RooArgList(*nsig1f,*f3Svs1S));
      RooRealVar *f3Svs2S = new RooRealVar("R_{#frac{3S}{2S}}","f3Svs2S",0.5,-1,1);
      RooFormulaVar *nsig3f= new RooFormulaVar("N32S","@0/@1",RooArgList(*f3Svs1S,*f2Svs1S));
      f2Svs1S->setConstant(kFALSE);
      f3Svs1S->setConstant(kFALSE);
      break;
    default:
      cout<<"Make a pick from choseFitParams!!!"<<endl;
      break;

    }  

  // bkg Chebychev
  RooRealVar *nbkgd   = new RooRealVar("nBkgd","nbkgd",nt*0.75,-100,10*nt);
  RooRealVar *bkg_a1  = new RooRealVar("a1_bkg", "bkg_{a1}", 0, -2, 2);
  RooRealVar *bkg_a2  = new RooRealVar("a2_Bkg", "bkg_{a2}", 0, -2, 2);
  RooRealVar *bkg_a3  = new RooRealVar("a3_Bkg", "bkg_{a3}", 0, -0.5, 2);

  RooAbsPdf  *bkgPdf  = new RooChebychev("bkgPdf","bkgPdf",
 					 *mass, RooArgList(*bkg_a1,*bkg_a2,*bkg_a3));

  //  likesign
  RooRealVar *nLikesignbkgd = new RooRealVar("NLikesignBkg","nlikesignbkgd",nt*0.75,0,10*nt);
  if (doTrkRot) 
    {
      nLikesignbkgd->setVal(TrkRotData->sumEntries());
      nLikesignbkgd->setError(sqrt(TrkRotData->sumEntries()));
    }
  else 
    {
      nLikesignbkgd->setVal(likesignData->sumEntries());
      nLikesignbkgd->setError(sqrt(likesignData->sumEntries()));
    }
  nbkgd->setVal(data0->sumEntries());
  nbkgd->setError(sqrt(data0->sumEntries()));
  if (doConstrainFit) 
    {
      RooGaussian* nLikesignbkgd_constr = new RooGaussian("nLikesignbkgd_constr","nLikesignbkgd_constr",
 							  *nLikesignbkgd,
 							  RooConst(nLikesignbkgd->getVal()),    //mean
 							  RooConst(nLikesignbkgd->getError())); //sigma
    }
  else nLikesignbkgd->setConstant(kTRUE);
  
  RooFormulaVar *nResidualbkgd = new RooFormulaVar("NResidualBkg","@0-@1",
 						   RooArgList(*nbkgd,*nLikesignbkgd));
 

 // *************************************************** bkgModel
  RooRealVar turnOn("turnOn","turnOn", 6., 0., 13.);
  RooRealVar width("width","width", 4.1,2.3, 20.);// MB 2.63
  RooRealVar decay("decay","decay", 3.5, 1, 100.);// MB: 3.39

  width.setConstant(false);
  decay.setConstant(false);
  turnOn.setConstant(false);

  RooGaussian* turnOn_constr;
  RooGaussian* width_constr;
  RooGaussian* decay_constr;

  //thisPdf: form of the bkg pdf
  //pdf_combinedbkgd; // total bkg pdf. usually form*normalization
  switch (bkgdModel) 
    {
    case 1 :  //(erf*exp+pol2 ) to fit the SS, then fix the shape and fit OS, in case of constrain option
      bkg_a3->setConstant(true);
      RooGenericPdf *thisPdf = new  RooGenericPdf("thisPdf","thisPdf",
 						  "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
 						  RooArgList(*mass,turnOn,width,decay));
      RooFitResult* fit_1st  = thisPdf->fitTo(*likesignData,Save(),NumCPU(4)) ; // likesign data
      if (doTrkRot) fit_1st  = thisPdf->fitTo(*TrkRotData,Save(),NumCPU(4)) ;
        
      if (doConstrainFit) 
 	{ // allow parameters to vary within cenral value from above fit + their sigma
 	  turnOn_constr = new RooGaussian("turnOn_constr","turnOn_constr",
 					   turnOn,
 					   RooConst(turnOn.getVal()),
 					   RooConst(turnOn.getError()));
 	  width_constr   = new RooGaussian("width_constr","width_constr",
 					   width,					   RooConst(width.getVal()),
 					   RooConst(width.getError()));
 	  decay_constr    = new RooGaussian("decay_constr","decay_constr",
 					   decay,
 					   RooConst(decay.getVal()),
 					   RooConst(decay.getError()));
 	}
      else 
 	{
 	  turnOn.setConstant(kTRUE);
 	  width.setConstant(kTRUE);
 	  decay.setConstant(kTRUE);
 	}
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
 					  RooArgList(*bkgPdf,*thisPdf),
 					  RooArgList(*nResidualbkgd,*nLikesignbkgd));
      break;
    case 2 : //us eRooKeysPdf to smooth the SS, then fit OS with pol+keys
      bkg_a3->setConstant(true);
      RooKeysPdf *thisPdf        = new RooKeysPdf("thisPdf","thisPdf",*mass,*likesignData,
 						  RooKeysPdf::MirrorBoth, 1.4);
      if (doTrkRot) thisPdf     = new RooKeysPdf("thisPdf","thisPdf",*mass,*TrkRotData,
 						  RooKeysPdf::MirrorBoth, 1.4);
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
 						      RooArgList(*bkgPdf,*thisPdf),
 						      RooArgList(*nResidualbkgd,*nLikesignbkgd));
      break;
    case 3 : //use error function to fit the OS directly
      bkg_a3->setConstant(true);
      RooGenericPdf *thisPdf            = new  RooGenericPdf("thisPdf","thisPdf",
 							     "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
 							     RooArgList(*mass,turnOn,width,decay));
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf("pdf_combinedbkgd","total combined background pdf",
 						     RooArgList(*thisPdf),
 						     RooArgList(*nbkgd));
      break;
      
    case 4 : //use pol 2 to fit the OS directly
      bkg_a3->setConstant(true);
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
 						      RooArgList(*bkgPdf),
 						      RooArgList(*nbkgd));
      break;
    case 5 : //use ( error function + polynomial 2) to fit the OS directly
      bkg_a3->setConstant(true);
      RooGenericPdf *thisPdf     = new  RooGenericPdf("thisPdf","thisPdf",
 						      "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
 						      RooArgList(*mass,turnOn,width,decay));
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
 						      RooArgList(*bkgPdf,*thisPdf),
 						      RooArgList(*nResidualbkgd,*nLikesignbkgd));
      break;
    case 6: // pol 3 to fit OS dirrectly
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
 						      RooArgList(*bkgPdf),
 						      RooArgList(*nbkgd));
      break;
    default :
      cout<<"Donno what you are talking about! Pick anothe fit option!"<<endl;
      break;
    }
  
  //###### the nominal fit with default pdf 
  RooFitResult* fit_2nd;// fit results
  RooAbsPdf  *pdf; // nominal PDF
  RooAbsPdf  *pdf_unconstr;
  if (doConstrainFit) 
    {
      pdf_unconstr   = new RooAddPdf ("pdf_unconstr","total signal+background pdf",
 				      RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd),
 				      RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
      RooProdPdf *pdf            = new RooProdPdf ("pdf","total constr pdf",
 						   RooArgSet(*pdf_unconstr,*turnOn_constr,*width_constr,*decay_constr,*nLikesignbkgd_constr));
      fit_2nd      = pdf->fitTo(*data,Constrained(),Save(kTRUE),Extended(kTRUE),Minos(doMinos),NumCPU(6));
    }
  else 
    {
      RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total signal+background pdf",
 						   RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd),
 						   RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
      fit_2nd       = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    }

  // *************************************************** plotting
  TCanvas c; c.cd();
  int nbins = ceil((mass_h-mass_l)/binw); 
  RooPlot* frame = mass->frame(Bins(nbins),Range(mass_l,mass_h));
 
  data->plotOn(frame,Name("theData"),MarkerSize(0.8));
  pdf->plotOn(frame,Name("thePdf")); // signal + bkg pdf
  // pdf->plotOn(frame,Name("thePdf"),VisualizeError(*fit_2nd,1),FillColor(kOrange)); /// change this, to see errrors
  if (plotBkg) {
    if (doTrkRot) TrkRotData->plotOn(frame,Name("TrkRotData"),MarkerSize(0.8),MarkerColor(kMagenta),MarkerStyle(22));
    else likesignData->plotOn(frame,Name("theLikeSignData"),MarkerSize(0.8),MarkerColor(kRed),MarkerStyle(24));
    
  }
  RooArgSet * pars = pdf->getParameters(data);

  //draw the fit lines and save plots
  data->plotOn(frame,Name("theData"),MarkerSize(0.8));
  pdf->plotOn(frame,Components("thisPdf"),Name("theBkg"),LineStyle(5),LineColor(kGreen));// polX bkg
  pdf->plotOn(frame,Components("pdf_combinedbkgd"),LineStyle(kDashed));// total bkg, blue
  if (plotBkg) 
    {
      if (doTrkRot && bkgdModel!=4 && bkgdModel!=6) pdf->plotOn(frame,Components("thisPdf"),Name("theLikeSign"),LineStyle(9),LineColor(kMagenta));// trk rot
      else
 	{
 	  if(bkgdModel!=4 && bkgdModel!=6)
 	    {
 	      pdf->plotOn(frame,Components("thisPdf"),Name("theLikeSign"),LineStyle(9),LineColor(kRed)); // like sign component
 	      RooArgSet * pars = thisPdf->getParameters(likesignData);
 	    }
 	}
    }
  // need this re-plotting, so the pulls pick the right fit
  pdf->plotOn(frame,Name("thePdf")); // signal + bkg pdf
  //----------------------------------------------------------------------  
   frame->SetTitle( "" );
  frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->Draw();
  // cout<<"Chi2 value is: "<< frame->chiSquare("pdf","data",fit_2nd->floatParsFinal().getSize());

  //plot parameters
  TLatex latex1;
  latex1.SetNDC();
  // latex1.DrawLatex(0.55,1.-0.05*3,"CMS Preliminary");
  // latex1.DrawLatex(0.55,1.-0.05*4,Form("%s",choseSampleLegend[choseSample]));
  // 
  
  // latex1.SetTextSize(0.035);
  // if(choseSample!=1 && choseSample!=2) latex1.DrawLatex(0.2,1.-0.05*3,Form("Cent. %d-%d%%",centrality_min,centrality_max));
  // // latex1.DrawLatex(0.2,1.-0.05*4,Form("%.2f< y <%.2f",upsYCut_min,upsYCut_max)); 
  // latex1.DrawLatex(0.2,1.-0.05*4,Form("|y_{CM}| < 1")); 
  // latex1.DrawLatex(0.2,1.-0.05*5.5,Form("p_{T}^{#mu} > %.0f GeV/c",muonpTcut));
  
  // c.SaveAs(figsDir+figName_+paramOn_+".png");
  // c.SaveAs(figsDir+figName_+paramOn_+".pdf");

  // -------------------------------------- plot pulls
  TCanvas cm("cm","cm");
  cm.cd();
  TPad *pPad1 = new TPad("pPad1","pPad1",0.05,0.35,0.95,0.97);
  pPad1->SetBottomMargin(0.00);
  pPad1->Draw();
  pPad1->cd();
  if(plotpars) 
       pdf->paramOn(frame,Layout(0.15,0.6,0.4),Layout(0.6,0.935,0.97));
  frame->Draw();
 
  latex1.SetTextSize(0.035);
  latex1.DrawLatex(0.15,1.-0.05*1.5,Form("%s",choseSampleLegend[choseSample]));
  latex1.DrawLatex(0.15,1.-0.05*2.5,Form("%s",choseSampleLumi[choseSample])); 
  if(choseSample!=1 && choseSample!=2 &&choseSample!=7) latex1.DrawLatex(0.15,1.-0.05*5.5,Form("Cent. %d-%d%%",centMin,centMax));
  latex1.DrawLatex(0.15,1.-0.05*3.5,Form("%2.f < y < %.2f",dimuYMin,dimuYMax)); 
  latex1.DrawLatex(0.15,1.-0.05*4.5,Form("p_{T}^{#mu} > %.1f GeV/c",muonpTcut));

  cm.cd(0);
  TPad *pPad2 = new TPad("pPad2","pPad2",0.05,0.05,0.95,0.35);
  pPad2->SetTopMargin(0.0);
  pPad2->Draw();
  pPad2->cd();
  // **************** create pulls; change the chi2 calculation also
  double chi2FromRoo = frame->chiSquare(fit_2nd->floatParsFinal().getSize());
  cout<<"!!!!!!!! chi2 from simple pull= "<<frame->chiSquare()<<"\t chi2 from RooFit= "<<chi2FromRoo <<endl;
  RooHist *phPullm = frame->pullHist(0,0,true); // this calcualtes the pulls taking the integral of the fit in each bin, instead of the value in the middle of the bid
  phPullm->SetName("phPullm");
  double *ypull     = phPullm->GetY();

  TH1 *phData      = data->createHistogram("invariantMass",nbins);
  double Chi2       = 0;
  int nFullBinsPull = 0;
  for (int i=0; i < nbins; i++) 
    {
      if (phData->GetBinContent(i) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i],2);
    }

  // for writing on canvas
  int nFitParam     = fit_2nd->floatParsFinal().getSize();
  int Dof           = nFullBinsPull - nFitParam;
  double UnNormChi2 = Chi2;
  Chi2             /= (nFullBinsPull - nFitParam);

  cout<<"!!!!! nFullBinsPull="<<nFullBinsPull<<"\tnFitParam="<<nFitParam<<endl;
  // draw pulls
  pPad2->cd();
  double mFrameMax = 0;
  RooPlot* prpFramePull = mass->frame(Title("Pull"),Bins(nbins),Range(mass_l,mass_h));
  prpFramePull->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  prpFramePull->GetXaxis()->CenterTitle(kTRUE);
  prpFramePull->GetYaxis()->SetTitleOffset(1.3);
  prpFramePull->GetYaxis()->SetTitle("Pull");
  prpFramePull->addPlotable(phPullm,"PX");
 
  if (prpFramePull->GetMinimum()*-1 > prpFramePull->GetMaximum()) mFrameMax = prpFramePull->GetMinimum()*-1;
  else mFrameMax = prpFramePull->GetMaximum();
  prpFramePull->SetMaximum(mFrameMax); 
  prpFramePull->SetMinimum(-1*mFrameMax); 
  prpFramePull->Draw();

  latex1.SetTextSize(0.085);
  double myChi2 = chi2FromRoo*Dof;
  latex1.DrawLatex(0.7,1.-0.05*3.5,Form("#chi^{2}/ndf = %2.1f/%d",myChi2,Dof));
   
  //cm.SaveAs(figsDir+figName_+paramOn_+"_pulls.png");
  cm.SaveAs(figsDir+figName_+paramOn_+"_pulls.pdf");

   //-------------------------------
  // print the final pdf parameters!
  
  // Print fit results 
  cout << endl << "figure name: "<< figName_ << endl;
  cout << "the nominal fit with the default pdf " << endl ;
  cout<<"p-value = "<< TMath::Prob(UnNormChi2,Dof)<<endl;

  // ------ calculate the single yields
  // if (choseSample<3)
  //   { RooRealVar* n1s_fitresult    = (RooRealVar*)fit_2nd->floatParsFinal().find("N_{#Upsilon(1S)}");
  //     RooRealVar* n2on1s_fitresult = (RooRealVar*)fit_2nd->floatParsFinal().find("R_{#frac{2S}{1S}}");
  //     RooRealVar* n3on1s_fitresult = (RooRealVar*)fit_2nd->floatParsFinal().find("N3onN1S");
      
  //     double myn2s = computeSingle(*n1s_fitresult,*n2on1s_fitresult);
  //     double myn3s = computeSingle(*n1s_fitresult,*n3on1s_fitresult);
      
  //     double myn2s1s_corr = fit_2nd->correlation(*n1s_fitresult,*n2on1s_fitresult);
  //     double myn3s1s_corr = fit_2nd->correlation(*n1s_fitresult,*n3on1s_fitresult);
      
  //     double myn2s_err = computeSingleError(*n1s_fitresult,*n2on1s_fitresult,myn2s1s_corr); 
  //     double myn3s_err = computeSingleError(*n1s_fitresult,*n3on1s_fitresult,myn3s1s_corr);
  //   } else if (choseSample>=3)
  //   {RooRealVar* n1s_fitresult    = (RooRealVar*)fit_2nd->floatParsFinal().find("N_{#Upsilon(1S)}");
  //     RooRealVar* n2on1s_fitresult = (RooRealVar*)fit_2nd->floatParsFinal().find("R_{#frac{2S}{1S}}");
  //     RooRealVar* n23on1s_fitresult = (RooRealVar*)fit_2nd->floatParsFinal().find("R_{#frac{2S+3S}{1S}}");

  //     double myn2s = computeSingle(*n1s_fitresult,*n2on1s_fitresult);
  //     double myn23s = computeSingle(*n1s_fitresult,*n23on1s_fitresult);
  //     //   double myn3s = myn23s - myn2s;
      
  //     double myn2s1s_corr = fit_2nd->correlation(*n1s_fitresult,*n2on1s_fitresult);
  //     double myn23s1s_corr = fit_2nd->correlation(*n1s_fitresult,*n23on1s_fitresult);
  //     //  double myn3s1s_corr = fit_2nd->correlation(*n1s_fitresult,(*n23on1s_fitresult - *n2on1s_fitresult));

  //     double myn2s_err = computeSingleError(*n1s_fitresult,*n2on1s_fitresult,myn2s1s_corr); 
  //     double myn23s_err = computeSingleError(*n1s_fitresult,*n23on1s_fitresult,myn23s1s_corr);
  //     //  double myn3s_err = computeSingleError(*n1s_fitresult,(*n23on1s_fitresult - *n2on1s_fitresult),myn3s1s_corr);

  //   }
  //-------- done
 
  // fit results 
  float baseNll = fit_2nd->minNll();
  float estimatedDistance2Minimum = fit_2nd->edm();
  fit_2nd->Print();

  // write out the fitting params
  string outParameters = outDatsDir+"/"+figName_+".txt";
  cout<<"Output file: " << outParameters<<endl;
  ofstream outfileFitResults;
  outfileFitResults.open(outParameters.c_str(), ios_base::out);

  fit_2nd->printMultiline(cout,1) << endl;
  
  
  // cout << "N2S= "<< myn2s << "\t err= " << myn2s_err <<endl;
  // if (choseSample>=3)							
  //   {  cout << "N23S= "<< myn23s << "\t err= " << myn23s_err <<endl;
  //   } else {
  //   cout << "N3S= "<< myn3s << "\t err= " << myn3s_err <<endl;
  // }
  // outfileFitResults <<endl;
  switch(choseFitParams) {
  case 0:
    outfileFitResults<<figName_<<" "<<nsig1f->getVal()<<" "<<nsig1f->getError()<<" "<<nsig2f->getVal()<<" "<<nsig2f->getError()<<" "<<nsig3f->getVal()<<" "<<nsig3f->getError()<<" "<<npow->getVal()<<" "<<npow->getError()<<" "<<alpha->getVal()<<" "<<alpha->getError()<<" "<<sigma1->getVal()<<" "<<sigma1->getError()<<" "<<2*nFitParam+2*baseNll<<" "<<fit_2nd->edm()<<" "<<UnNormChi2<<" "<<UnNormChi2/Dof<<" "<<TMath::Prob(UnNormChi2,Dof)<<" "<<Dof<<" "<<nFitParam<<" "<<baseNll<< endl;
    break;
  case 1 :
 
  outfileFitResults<<figName_<<" "<<nsig1f->getVal()<<" "<<nsig1f->getError()<<" "<<f2Svs1S->getVal()<<" "<<f2Svs1S->getError()<<" "<<f3Svs1S->getVal()<<" "<<f3Svs1S->getError()<<" "<<npow->getVal()<<" "<<npow->getError()<<" "<<alpha->getVal()<<" "<<alpha->getError()<<" "<<sigma1->getVal()<<" "<<sigma1->getError()<<" "<<2*nFitParam+2*baseNll<<" "<<fit_2nd->edm()<<" "<<UnNormChi2<<" "<<UnNormChi2/Dof<<" "<<TMath::Prob(UnNormChi2,Dof)<<" "<<Dof<<" "<<nFitParam<<" "<<baseNll<< endl;
  break;
  case 2 :
     outfileFitResults<<figName_<<" "<<nsig1f->getVal()<<" "<<nsig1f->getError()<<" "<<f2Svs1S->getVal()<<" "<<f2Svs1S->getError()<<" "<<f23vs1S->getVal()<<" "<<f23vs1S->getError()<<" "<<npow->getVal()<<" "<<npow->getError()<<" "<<alpha->getVal()<<" "<<alpha->getError()<<" "<<sigma1->getVal()<<" "<<sigma1->getError()<<" "<<2*nFitParam+2*baseNll<<" "<<fit_2nd->edm()<<" "<<UnNormChi2<<" "<<UnNormChi2/Dof<<" "<<TMath::Prob(UnNormChi2,Dof)<<" "<<Dof<<" "<<nFitParam<<" "<<baseNll<< endl;
     break;
  case 3 :
    outfileFitResults<<figName_<<" "<<nsig1f->getVal()<<" "<<nsig1f->getError()<<" "<<f2Svs1S->getVal()<<" "<<f2Svs1S->getError()<<" "<<f3Svs1S->getVal()<<" "<<f3Svs1S->getError()<<" "<<f3Svs2S->getVal()<<" "<<f3Svs2S->getError()<<" "<<alpha->getVal()<<" "<<alpha->getError()<<" "<<sigma1->getVal()<<" "<<sigma1->getError()<<" "<<2*nFitParam+2*baseNll<<" "<<fit_2nd->edm()<<" "<<UnNormChi2<<" "<<UnNormChi2/Dof<<" "<<TMath::Prob(UnNormChi2,Dof)<<" "<<Dof<<" "<<nFitParam<<" "<<baseNll<< endl;
     break;
  default : break;
    
   
  }


 // << "FreeParam= "<< nFitParam <<" "<<baseNll<< " " <<estimatedDistance2Minimum<< " " << 2*nFitParam+2*baseNll<< " " << UnNormChi2 << " " << Dof<<" "<< TMath::Prob(UnNormChi2,Dof) << " " << endl;
 //  outfileFitResults <<endl;

 //  outfileFitResults<< "alpha= "<<alpha->printValue(outfileFitResults)  << "\t npow= "<<npow->printValue(outfileFitResults)<< "\t sigma1= "<<sigma1->printValue(outfileFitResults) << "\t turnOn= "<<turnOn.printValue(outfileFitResults) <<endl;

 //  outfileFitResults << "rooFitChSquare= "<< chi2FromRoo*Dof <<endl;
  outfileFitResults.close();
  // ConfidencInterval(0.683, f3Svs1S, data, pdf);
  
}


//_______________________________________________________________________
double computeSingle(RooRealVar& x, RooRealVar& y) 
{
  // pass the 1S(x) and xS/1S(y) ratios and calcualte the xS yield
   return x.getVal() * y.getVal();
}


double computeSingleError(RooRealVar& x, RooRealVar& y, double correlation = 0.) 
{
  // pass the 1S(x) and xS/1S(y) ratios and calcualte the xS yield error

  double err2 = (x.getError()*x.getError())/(x.getVal()*x.getVal()) 
    + (y.getError()*y.getError())/(y.getVal()*y.getVal()) 
    + 2.*(x.getError()*y.getError())/(x.getVal()*y.getVal())*correlation;
  
  return fabs(computeSingle(x,y))*sqrt(err2);
}


//_______________________________________________________________________


//calculate the confidence interval with RooStats
pair<double, double> ConfidencInterval(float CI, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf)  {  
	ProfileLikelihoodCalculator pl(*data,*pdf,*fnll);
	pl.SetConfidenceLevel(CI); 
	int ci = 100*CI;
	LikelihoodInterval* interval = pl.GetInterval();
	LikelihoodIntervalPlot plot(interval);
	TCanvas c4; c4.cd(); 
	plot.SetRange(0.,3.,-.05,0.5);
	plot.Draw();
	TLatex latexCI;
	latexCI.SetNDC();
	latexCI.SetTextSize(0.035);
	latexCI.DrawLatex(0.5,1.-0.05*2,Form("R_{#frac{3S}{1S}} %d % C.I.",ci));
	latexCI.DrawLatex(0.5,1.-0.05*3,Form("Upper limit: %f",interval->UpperLimit(*fnll)));
	TString intrvlName = fnll->GetTitle();
		// print out the iterval on the Parameter of Interest
	cout <<endl<< CI <<"\% interval on " <<fnll->GetName()<<" is : ["<<
		interval->LowerLimit(*fnll) << ", "<<
		interval->UpperLimit(*fnll) << "] "<<endl;
	pair<double, double> CnfdncIntrvl;
	CnfdncIntrvl.first  = interval->LowerLimit(*fnll);
	CnfdncIntrvl.second = interval->UpperLimit(*fnll);
	c4.SaveAs("pdfOutput/1806/HIMBProfileTest.pdf");

	return CnfdncIntrvl;
}

