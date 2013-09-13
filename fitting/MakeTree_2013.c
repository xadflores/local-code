#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>

#include <TH1.h>
#include <TH2D.h>

#include <TBranch.h>
#include <TCanvas.h>
#include "TClonesArray.h"
#include <TDirectory.h>
#include <TFile.h>
#include "TH1F.h"
#include <TLatex.h>
#include <TLegend.h>
#include "TLorentzVector.h"
#include <TMath.h>
#include "TRandom.h"
#include <TStyle.h>
#include <TSystem.h>
#include "TTree.h"
#include "TString.h"

// miscellaneous  
#include <fstream>
#include <map>
#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#endif



// this macro assumes all cuts are applied in the tree, 
/*
Legend: Sign pairs:  +- sign=0; ++ sign=1, -- sign=2
Legend: muon pairs: Glb-Glb id=0, Glb-Trk = 1, Trk-Trk = 2

n --> 2^n -1
0  theTriggerNames.push_back("NoTrigger");
1  theTriggerNames.push_back("HLT_PAL1DoubleMuOpen_v1");
2  theTriggerNames.push_back("HLT_PAL1DoubleMu0_HighQ_v1");
3  theTriggerNames.push_back("HLT_PAL2DoubleMu3_v1");
4  theTriggerNames.push_back("HLT_PAMu3_v1");
5  theTriggerNames.push_back("HLT_PAMu7_v1");
6  theTriggerNames.push_back("HLT_PAMu12_v1");

 */


void MakeTree_2013(const char* inputOniaTree = "/afs/cern.ch/user/n/nfilipov/Public/Jpsi_Histos_PbPb_RegIT_glbglb.root",
		   //  const char* runNumber = "211739-211831",
		   const char* runNumber = "2011-2011",
		   //const char* runNumber = "a",
		   float newVtxProbCut   = 0.01,// default in the input tree is already 0.01
		   float ptcut           = 4., // single muon cut pt
		   int nTriggerBit       = 1, 
		   // const char* dataSource= "upsiMiniTree_pp276tev_5p41_ptmu4",
		   const char* dataSource= "upsiMiniTree_AA2p76tev_ptmu4_july09",
		   // const char* dataSource= "upsiMiniTree_pp7tev_dimu0v1_ptmu4",
		   //   const char* dataSource= "upsiMiniTree_aa276tevC50100_ppofficial_trktrk_ptmu4",
		   bool isAArereco       = true,// only for newly processed AA; old tree is with triger v1,v2 included
		   bool bAllTriggers     = false,
		   bool addExtraCentrality = true,
		   bool addExtraTrack      = true
		   )
{
  gROOT->Macro("setTDRStyle_modified.C+");
 
  float mass_min = 7.0; 
  float mass_max = 14.0;
  int   nBins    = 70;
  
  TFile *f = TFile::Open(Form("%s",inputOniaTree));
 
  TTree *t = (TTree*)f->Get("myTree");
  Long64_t nentries = t->GetEntries();
  cout<<"Got the tree!"<<endl;

  const int NMAX=100;
  UInt_t eventNb;
  UInt_t runNb;
  Int_t Centrality;
  Int_t HLTriggers;
  Int_t Reco_QQ_matchNtrkVtx;
  Int_t Reco_QQ_size;
  Int_t Reco_QQ_trig[NMAX];
  Int_t Reco_QQ_type[NMAX];
  Int_t Reco_QQ_sign[NMAX];
  Int_t Reco_QQ_NtrkDeltaR03[NMAX];
  Int_t Reco_QQ_NtrkDeltaR04[NMAX];
  Int_t Reco_QQ_NtrkDeltaR05[NMAX];

  Float_t Reco_QQ_VtxProb[NMAX];
  TClonesArray *Reco_QQ_4mom = 0;
  TClonesArray *Reco_QQ_mupl_4mom = 0;
  TClonesArray *Reco_QQ_mumi_4mom = 0;

  Int_t Reco_mu_size;
  Int_t Reco_mu_type[NMAX];
  Int_t Reco_mu_charge[NMAX];
  TClonesArray *Reco_mu_4mom = 0;
  int nPlusMu;
  int nMinusMu;
  float muPt[NMAX];
  float muEta[NMAX];
  float muPhi[NMAX];
  //Int_t nReco_QQ;
  Float_t invariantMass;
  Int_t QQtrig;
  Int_t QQsign;
  Int_t QQTrkDr03;
  Int_t QQTrkDr04;
  Int_t QQTrkDr05;
  float weight;
  float weight2;
  Float_t upsPt;
  Float_t upsEta;
  Float_t upsPhi;
  Float_t upsRapidity;	
  Float_t vProb;
  Float_t muPlusPt;
  Float_t muMinusPt;
  Float_t muPlusEta;
  Float_t muMinusEta;
  Float_t muPlusPhi;
  Float_t muMinusPhi;
  Float_t RmuPlusPhi = 0;
  Float_t RmuMinusPhi = 0;
  Int_t upsMatchNtrkVtx;
  // centrliaty extra stuff
  Int_t Npix, NpixelTracks, Ntracks;
  Float_t zVtx;
  Float_t SumET_HF, SumET_HFplus, SumET_HFminus, SumET_HFplusEta4, SumET_HFminusEta4, SumET_EB, SumET_ET, SumET_EE, SumET_EEplus, SumET_EEminus, SumET_ZDC, SumET_ZDCplus, SumET_ZDCminus;

  // track extra stuff
  const int NTRKMAX=100000;
  int Reco_trk_size=0;
  TClonesArray *Reco_trk_4mom = 0;
  float trkPt[NTRKMAX];
  float trkEta[NTRKMAX];
  float trkPhi[NTRKMAX];
  //---------------------------
  TBranch *b_Centrality; //!
  TBranch *b_eventNb; //!
  TBranch *b_runNb; //!
  TBranch *b_HLTriggers; //!
  TBranch *b_Reco_QQ_size; //!
  TBranch *b_Reco_QQ_trig; //!
  TBranch *b_Reco_QQ_type; //!
  TBranch *b_Reco_QQ_sign; //!
  TBranch *b_Reco_QQ_VtxProb; //!
  TBranch *b_Reco_QQ_4mom; //!

  TBranch *b_Reco_QQ_mupl_4mom; //!
  TBranch *b_Reco_QQ_mumi_4mom; //!

  TBranch *b_Reco_mu_size; //!
  TBranch *b_Reco_mu_type; //!
  TBranch *b_Reco_mu_charge; //!
  TBranch *b_Reco_mu_4mom; //!

  TBranch *b_Reco_QQ_NtrkDeltaR03;//!
  TBranch *b_Reco_QQ_NtrkDeltaR04;//!
  TBranch *b_Reco_QQ_NtrkDeltaR05;//!

  TBranch *b_Reco_QQ_matchNtrkVtx;//!
  //------------------------------
  TBranch *b_zVtx; //!
  TBranch *b_Npix; //!
  TBranch *b_NpixelTracks; //!
  TBranch *b_Ntracks; //!
  TBranch *b_SumET_HF; //!
  TBranch *b_SumET_HFplus; //!
  TBranch *b_SumET_HFminus; //!
  TBranch *b_SumET_HFplusEta4; //!
  TBranch *b_SumET_HFminusEta4; //!
  TBranch *b_SumET_ET; //!
  TBranch *b_SumET_EE; //!
  TBranch *b_SumET_EB; //!
  TBranch *b_SumET_EEplus; //!
  TBranch *b_SumET_EEminus; //!
  TBranch *b_SumET_ZDC; //!
  TBranch *b_SumET_ZDCplus; //!
  TBranch *b_SumET_ZDCminus; //!

  //------------------------------
  TBranch *b_Reco_trk_size; //! 
  TBranch *b_Reco_trk_4mom;//!


  // ###### read input TTree
  t->SetBranchAddress("Centrality",    &Centrality,   &b_Centrality);
  t->SetBranchAddress("eventNb",       &eventNb,      &b_eventNb);
  t->SetBranchAddress("runNb",         &runNb,        &b_runNb);
  t->SetBranchAddress("HLTriggers",    &HLTriggers,   &b_HLTriggers);
  t->SetBranchAddress("Reco_QQ_size",  &Reco_QQ_size, &b_Reco_QQ_size);
  t->SetBranchAddress("Reco_QQ_trig",  &Reco_QQ_trig, &b_Reco_QQ_trig);
  t->SetBranchAddress("Reco_QQ_type",  &Reco_QQ_type, &b_Reco_QQ_type);
  t->SetBranchAddress("Reco_QQ_sign",  &Reco_QQ_sign, &b_Reco_QQ_sign);
  t->SetBranchAddress("Reco_QQ_VtxProb",   &Reco_QQ_VtxProb,   &b_Reco_QQ_VtxProb);
  t->SetBranchAddress("Reco_QQ_4mom",      &Reco_QQ_4mom,      &b_Reco_QQ_4mom);
  t->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  t->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  
  t->SetBranchAddress("Reco_mu_size",  &Reco_mu_size,   &b_Reco_mu_size);
  t->SetBranchAddress("Reco_mu_type",  &Reco_mu_type,   &b_Reco_mu_type);
  t->SetBranchAddress("Reco_mu_charge",&Reco_mu_charge, &b_Reco_mu_charge);
  t->SetBranchAddress("Reco_mu_4mom",  &Reco_mu_4mom,   &b_Reco_mu_4mom);

  t->SetBranchAddress("Reco_QQ_NtrkDeltaR03",  &Reco_QQ_NtrkDeltaR03, &b_Reco_QQ_NtrkDeltaR03);
  t->SetBranchAddress("Reco_QQ_NtrkDeltaR04",  &Reco_QQ_NtrkDeltaR04, &b_Reco_QQ_NtrkDeltaR04);
  t->SetBranchAddress("Reco_QQ_NtrkDeltaR05",  &Reco_QQ_NtrkDeltaR05, &b_Reco_QQ_NtrkDeltaR05);

  // extra centrality variables
  t->SetBranchAddress("zVtx",             &zVtx,              &b_zVtx);
  t->SetBranchAddress("Npix",             &Npix,              &b_Npix);
  t->SetBranchAddress("NpixelTracks",     &NpixelTracks,      &b_NpixelTracks);
  t->SetBranchAddress("Ntracks",          &Ntracks,           &b_Ntracks);
  t->SetBranchAddress("SumET_HF",         &SumET_HF,          &b_SumET_HF);
  t->SetBranchAddress("SumET_HFplus",     &SumET_HFplus,      &b_SumET_HFplus);
  t->SetBranchAddress("SumET_HFminus",    &SumET_HFminus,     &b_SumET_HFminus);
  t->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4,  &b_SumET_HFplusEta4);
  t->SetBranchAddress("SumET_HFminusEta4",&SumET_HFminusEta4, &b_SumET_HFminusEta4);
  t->SetBranchAddress("SumET_ET",         &SumET_ET,          &b_SumET_ET);
  t->SetBranchAddress("SumET_EE",         &SumET_EE,          &b_SumET_EE);
  t->SetBranchAddress("SumET_EB",         &SumET_EB,          &b_SumET_EB);
  t->SetBranchAddress("SumET_EEplus",     &SumET_EEplus,      &b_SumET_EEplus);
  t->SetBranchAddress("SumET_EEminus",    &SumET_EEminus,     &b_SumET_EEminus);
  t->SetBranchAddress("SumET_ZDC"    ,    &SumET_ZDC,         &b_SumET_ZDC);
  t->SetBranchAddress("SumET_ZDCplus",    &SumET_ZDCplus,     &b_SumET_ZDCplus);
  t->SetBranchAddress("SumET_ZDCminus",   &SumET_ZDCminus,    &b_SumET_ZDCminus);
  
  // extra track variable  
  t->SetBranchAddress("Reco_QQ_matchNtrkVtx",  &Reco_QQ_matchNtrkVtx,   &b_Reco_QQ_matchNtrkVtx);
  t->SetBranchAddress("Reco_trk_size",  &Reco_trk_size,   &b_Reco_trk_size);
  t->SetBranchAddress("Reco_trk_4mom",  &Reco_trk_4mom,   &b_Reco_trk_4mom);

  // #### define control histograms
  TH1F *h_QQ_mass   = new TH1F("h_QQ_mass","",nBins, mass_min,mass_max);  // all OS
  TH1F *h_QQ_mass_1 = new TH1F("h_QQ_mass_1","",nBins, mass_min,mass_max);// OS in acceptance
  TH1F *h_QQ_mass_2 = new TH1F("h_QQ_mass_2","",nBins, mass_min, mass_max);// SS 
  
  // ##### output file
  TFile *f1 = new TFile(Form("../dimuonTree_%s_Run%s_trigBit%d_allTriggers%d.root",dataSource,runNumber,nTriggerBit,bAllTriggers),"RECREATE");

  TTree *MuTree              = new TTree("MuTree","MuTree");
  TTree *UpsilonTree         = new TTree("UpsilonTree","UpsilonTree");
  TTree *UpsilonTree_allsign = new TTree("UpsilonTree_allsign","UpsilonTree_allsign");
  TTree *UpsilonTree_trkRot  = new TTree("UpsilonTree_trkRot","UpsilonTree_trkRot");

  MuTree->Branch("Reco_mu_size",   &Reco_mu_size,   "Reco_mu_size/I");
  MuTree->Branch("nPlusMu",        &nPlusMu,        "nPlusMu/I");
  MuTree->Branch("nMinusMu",       &nMinusMu,       "nMinusMu/I");
  MuTree->Branch("muPt",           &muPt,           "muPt[Reco_mu_size]/F");
  MuTree->Branch("muEta",          &muEta,          "muEta[Reco_mu_size]/F");
  MuTree->Branch("muPhi",          &muPhi,          "muPhi[Reco_mu_size]/F");
  MuTree->Branch("Reco_mu_type",   &Reco_mu_type,    "Reco_mu_type[Reco_mu_size]/I");
  MuTree->Branch("Reco_mu_charge", &Reco_mu_charge,  "Reco_mu_charge[Reco_mu_size]/I" );
  MuTree->Branch("eventNb",        &eventNb,        "eventNb/I");
  MuTree->Branch("runNb",          &runNb,          "runNb/I");

  //UpsilonTree->Branch("nReco_QQ", &nReco_QQ, "nReco_QQ/I");
  UpsilonTree->Branch("Centrality",    &Centrality,    "Centrality/I");
  UpsilonTree->Branch("HLTriggers",    &HLTriggers,    "HLTriggers/I");
  UpsilonTree->Branch("QQtrig",        &QQtrig,        "QQtrig/I");
  UpsilonTree->Branch("QQsign",        &QQsign,        "QQsign/I");
  UpsilonTree->Branch("QQTrkDr03",     &QQTrkDr03,     "QQTrkDr03/I");
  UpsilonTree->Branch("QQTrkDr04",     &QQTrkDr04,     "QQTrkDr04/I");
  UpsilonTree->Branch("QQTrkDr05",     &QQTrkDr05,     "QQTrkDr05/I");
  UpsilonTree->Branch("vProb",         &vProb,     "vProb/F");
  UpsilonTree->Branch("eventNb",       &eventNb,     "eventNb/I");
  UpsilonTree->Branch("runNb",         &runNb,     "runNb/I");
  UpsilonTree->Branch("invariantMass", &invariantMass, "invariantMass/F");
  UpsilonTree->Branch("upsPt",         &upsPt, "upsPt/F");
  UpsilonTree->Branch("upsEta",        &upsEta, "upsEta/F");
  UpsilonTree->Branch("upsPhi",        &upsPhi, "upsPhi/F");
  UpsilonTree->Branch("upsRapidity",   &upsRapidity, "upsRapidity/F");
  UpsilonTree->Branch("muPlusPt",      &muPlusPt, "muPlusPt/F");
  UpsilonTree->Branch("muMinusPt",     &muMinusPt, "muMinusPt/F");
  UpsilonTree->Branch("muPlusEta",     &muPlusEta, "muPlusEta/F");
  UpsilonTree->Branch("muMinusEta",    &muMinusEta, "muMinusEta/F");
  UpsilonTree->Branch("muPlusPhi",     &muPlusPhi, "muPlusPhi/F");
  UpsilonTree->Branch("muMinusPhi",    &muMinusPhi, "muMinusPhi/F");
  if(addExtraCentrality)
    {
      UpsilonTree->Branch("Npix",&Npix,"Npix/I");
      UpsilonTree->Branch("NpixelTracks",&NpixelTracks,"NpixelTracks/I");
      UpsilonTree->Branch("Ntracks", &Ntracks, "Ntracks/I");
      UpsilonTree->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
      UpsilonTree->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
      UpsilonTree->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
      UpsilonTree->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
      UpsilonTree->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
      UpsilonTree->Branch("SumET_ET",&SumET_ET,"SumET_ET/F");
      UpsilonTree->Branch("SumET_EE",&SumET_EE,"SumET_EE/F");
      UpsilonTree->Branch("SumET_EB",&SumET_EB,"SumET_EB/F");
      UpsilonTree->Branch("SumET_EEplus",&SumET_EEplus,"SumET_EEplus/F");
      UpsilonTree->Branch("SumET_EEminus",&SumET_EEminus,"SumET_EEminus/F");
      UpsilonTree->Branch("SumET_ZDC",&SumET_ZDC,"SumET_ZDC/F");
      UpsilonTree->Branch("SumET_ZDCplus",&SumET_ZDCplus,"SumET_ZDCplus/F");
      UpsilonTree->Branch("SumET_ZDCminus",&SumET_ZDCminus,"SumET_ZDCminus/F");
    }
  if(addExtraTrack)
    {
      UpsilonTree->Branch("Reco_trk_size",   &Reco_trk_size,   "Reco_trk_size/I");
      UpsilonTree->Branch("upsMatchNtrkVtx", &upsMatchNtrkVtx, "upsMatchNtrkVtx/I");
      UpsilonTree->Branch("trkPt", &trkPt,"trkPt[Reco_trk_size]/F");
      UpsilonTree->Branch("trkEta",&trkEta,"trkEta[Reco_trk_size]/F");
      UpsilonTree->Branch("trkPhi",&trkPhi,"trkPhi[Reco_trk_size]/F");
    }

  UpsilonTree_allsign->Branch("Centrality", &Centrality,    "Centrality/I");
  UpsilonTree_allsign->Branch("HLTriggers", &HLTriggers,    "HLTriggers/I");
  UpsilonTree_allsign->Branch("QQtrig",     &QQtrig,    "QQtrig/I");
  UpsilonTree_allsign->Branch("QQsign",     &QQsign,    "QQsign/I");
  UpsilonTree_allsign->Branch("QQTrkDr03",     &QQTrkDr03,     "QQTrkDr03/I");
  UpsilonTree_allsign->Branch("QQTrkDr04",     &QQTrkDr04,     "QQTrkDr04/I");
  UpsilonTree_allsign->Branch("QQTrkDr05",     &QQTrkDr05,     "QQTrkDr05/I");
  UpsilonTree_allsign->Branch("weight",     &weight,    "weight/F");
  UpsilonTree_allsign->Branch("weight2",    &weight2,    "weight2/F");
  UpsilonTree_allsign->Branch("vProb",      &vProb,     "vProb/F");
  UpsilonTree_allsign->Branch("eventNb",    &eventNb,     "eventNb/I");
  UpsilonTree_allsign->Branch("runNb",      &runNb,     "runNb/I");
  UpsilonTree_allsign->Branch("invariantMass", &invariantMass, "invariantMass/F");
  UpsilonTree_allsign->Branch("upsPt",      &upsPt, "upsPt/F");
  UpsilonTree_allsign->Branch("upsEta",     &upsEta, "upsEta/F");
  UpsilonTree_allsign->Branch("upsPhi",     &upsPhi, "upsPhi/F");
  UpsilonTree_allsign->Branch("upsRapidity",&upsRapidity, "upsRapidity/F");
  UpsilonTree_allsign->Branch("muPlusPt",   &muPlusPt, "muPlusPt/F");
  UpsilonTree_allsign->Branch("muMinusPt",  &muMinusPt, "muMinusPt/F");
  UpsilonTree_allsign->Branch("muPlusEta",  &muPlusEta, "muPlusEta/F");
  UpsilonTree_allsign->Branch("muMinusEta", &muMinusEta, "muMinusEta/F");
  UpsilonTree_allsign->Branch("muPlusPhi",  &muPlusPhi, "muPlusPhi/F");
  UpsilonTree_allsign->Branch("muMinusPhi", &muMinusPhi, "muMinusPhi/F");
  
   if(addExtraCentrality)
    {
      UpsilonTree_allsign->Branch("Npix",&Npix,"Npix/I");
      UpsilonTree_allsign->Branch("NpixelTracks",&NpixelTracks,"NpixelTracks/I");
      UpsilonTree_allsign->Branch("Ntracks", &Ntracks, "Ntracks/I");
      UpsilonTree_allsign->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
      UpsilonTree_allsign->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
      UpsilonTree_allsign->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
      UpsilonTree_allsign->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
      UpsilonTree_allsign->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
      UpsilonTree_allsign->Branch("SumET_ET",&SumET_ET,"SumET_ET/F");
      UpsilonTree_allsign->Branch("SumET_EE",&SumET_EE,"SumET_EE/F");
      UpsilonTree_allsign->Branch("SumET_EB",&SumET_EB,"SumET_EB/F");
      UpsilonTree_allsign->Branch("SumET_EEplus",&SumET_EEplus,"SumET_EEplus/F");
      UpsilonTree_allsign->Branch("SumET_EEminus",&SumET_EEminus,"SumET_EEminus/F");
      UpsilonTree_allsign->Branch("SumET_ZDC",&SumET_ZDC,"SumET_ZDC/F");
      UpsilonTree_allsign->Branch("SumET_ZDCplus",&SumET_ZDCplus,"SumET_ZDCplus/F");
      UpsilonTree_allsign->Branch("SumET_ZDCminus",&SumET_ZDCminus,"SumET_ZDCminus/F");
    }

   if(addExtraTrack)
     {
       UpsilonTree_allsign->Branch("upsMatchNtrkVtx", &upsMatchNtrkVtx, "upsMatchNtrkVtx/I");
       UpsilonTree_allsign->Branch("Reco_trk_size",   &Reco_trk_size,   "Reco_trk_size/I");
       UpsilonTree_allsign->Branch("trkPt",           &trkPt,           "trkPt[Reco_trk_size]/F");
       UpsilonTree_allsign->Branch("trkEta",          &trkEta,          "trkEta[Reco_trk_size]/F");
       UpsilonTree_allsign->Branch("trkPhi",          &trkPhi,          "trkPhi[Reco_trk_size]/F");
    }

  UpsilonTree_trkRot->Branch("Centrality", &Centrality,    "Centrality/I");
  UpsilonTree_trkRot->Branch("HLTriggers", &HLTriggers,    "HLTriggers/I");
  UpsilonTree_trkRot->Branch("QQtrig",     &QQtrig,    "QQtrig/I");
  UpsilonTree_trkRot->Branch("QQsign",     &QQsign,    "QQsign/I");
  UpsilonTree_trkRot->Branch("QQTrkDr03",     &QQTrkDr03,     "QQTrkDr03/I");
  UpsilonTree_trkRot->Branch("QQTrkDr04",     &QQTrkDr04,     "QQTrkDr04/I");
  UpsilonTree_trkRot->Branch("QQTrkDr05",     &QQTrkDr05,     "QQTrkDr05/I");

  UpsilonTree_trkRot->Branch("weight",     &weight,    "weight/F");
  UpsilonTree_trkRot->Branch("weight2",    &weight2,    "weight2/F");
  UpsilonTree_trkRot->Branch("vProb",      &vProb,     "vProb/F");
  UpsilonTree_trkRot->Branch("eventNb",    &eventNb,     "eventNb/I");
  UpsilonTree_trkRot->Branch("runNb",      &runNb,     "runNb/I");
  UpsilonTree_trkRot->Branch("invariantMass", &invariantMass, "invariantMass/F");
  UpsilonTree_trkRot->Branch("upsPt",      &upsPt, "upsPt/F");
  UpsilonTree_trkRot->Branch("upsEta",     &upsEta, "upsEta/F");
  UpsilonTree_trkRot->Branch("upsPhi",     &upsPhi, "upsPhi/F");
  UpsilonTree_trkRot->Branch("upsRapidity",&upsRapidity, "upsRapidity/F");
  UpsilonTree_trkRot->Branch("muPlusPt",   &muPlusPt, "muPlusPt/F");
  UpsilonTree_trkRot->Branch("muMinusPt",  &muMinusPt, "muMinusPt/F");
  UpsilonTree_trkRot->Branch("muPlusEta",  &muPlusEta, "muPlusEta/F");
  UpsilonTree_trkRot->Branch("muMinusEta", &muMinusEta, "muMinusEta/F");
  UpsilonTree_trkRot->Branch("muPlusPhi",  &muPlusPhi, "muPlusPhi/F");
  UpsilonTree_trkRot->Branch("muMinusPhi", &muMinusPhi, "muMinusPhi/F");
  if(addExtraCentrality)
    {
      UpsilonTree_trkRot->Branch("Npix",&Npix,"Npix/I");
      UpsilonTree_trkRot->Branch("NpixelTracks",&NpixelTracks,"NpixelTracks/I");
      UpsilonTree_trkRot->Branch("Ntracks", &Ntracks, "Ntracks/I");
      UpsilonTree_trkRot->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
      UpsilonTree_trkRot->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
      UpsilonTree_trkRot->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
      UpsilonTree_trkRot->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
      UpsilonTree_trkRot->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
      UpsilonTree_trkRot->Branch("SumET_ET",&SumET_ET,"SumET_ET/F");
      UpsilonTree_trkRot->Branch("SumET_EE",&SumET_EE,"SumET_EE/F");
      UpsilonTree_trkRot->Branch("SumET_EB",&SumET_EB,"SumET_EB/F");
      UpsilonTree_trkRot->Branch("SumET_EEplus",&SumET_EEplus,"SumET_EEplus/F");
      UpsilonTree_trkRot->Branch("SumET_EEminus",&SumET_EEminus,"SumET_EEminus/F");
      UpsilonTree_trkRot->Branch("SumET_ZDC",&SumET_ZDC,"SumET_ZDC/F");
      UpsilonTree_trkRot->Branch("SumET_ZDCplus",&SumET_ZDCplus,"SumET_ZDCplus/F");
      UpsilonTree_trkRot->Branch("SumET_ZDCminus",&SumET_ZDCminus,"SumET_ZDCminus/F");
    }

  if(addExtraTrack)
    {
      UpsilonTree_trkRot->Branch("upsMatchNtrkVtx", &upsMatchNtrkVtx, "upsMatchNtrkVtx/I");
      UpsilonTree_trkRot->Branch("Reco_trk_size",   &Reco_trk_size,   "Reco_trk_size/I");
      UpsilonTree_trkRot->Branch("trkPt",           &trkPt,           "trkPt[Reco_trk_size]/F");
      UpsilonTree_trkRot->Branch("trkEta",          &trkEta,          "trkEta[Reco_trk_size]/F");
      UpsilonTree_trkRot->Branch("trkPhi",          &trkPhi,          "trkPhi[Reco_trk_size]/F");
    }

  //____________________________________ loop over all events
  for (int i=0; i<nentries; i++) 
    {
      t->GetEntry(i);
      
      if (i%100000==0) cout<<i<<endl;

      // ##### single muon tree
      // countng + and - single muons
      nPlusMu=0;
      nMinusMu=0;
      for(int iMu = 0; iMu < Reco_mu_size; iMu++)
	{
	  if (Reco_mu_charge[iMu] == 1) nPlusMu++;
	  else nMinusMu++;
	  TLorentzVector *Reco_mu = (TLorentzVector *) Reco_mu_4mom->At(iMu);
	  muPt[iMu]=Reco_mu->Pt();
	  muEta[iMu]=Reco_mu->Eta();
	  muPhi[iMu]=Reco_mu->Phi();
	}
      MuTree->Fill();

      // ntrk loop
       for(int iTrk = 0; iTrk < Reco_trk_size; iTrk++)
	{
	  TLorentzVector *Reco_trk = (TLorentzVector *) Reco_trk_4mom->At(iTrk);
	  trkPt[iTrk]=Reco_trk->Pt();
	  trkEta[iTrk]=Reco_trk->Eta();
	  trkPhi[iTrk]=Reco_trk->Phi();
	}

      //----------------------------------------------------
       upsMatchNtrkVtx = Reco_QQ_matchNtrkVtx; 
      for(int iQQ = 0; iQQ < Reco_QQ_size; iQQ++)
	{
	  vProb  = Reco_QQ_VtxProb[iQQ];
	  QQtrig = Reco_QQ_trig[iQQ];
	  QQsign = Reco_QQ_sign[iQQ];
	  
	  QQTrkDr03 = Reco_QQ_NtrkDeltaR03[iQQ];
	  QQTrkDr04 = Reco_QQ_NtrkDeltaR04[iQQ];
	  QQTrkDr05 = Reco_QQ_NtrkDeltaR05[iQQ];

	  // have to red a bti about the weighting Zhen calculates ...
	  if (Reco_QQ_sign[iQQ] == 0) // opposite sign
	    { 
	      weight  = 1;
	      weight2 = 1;
	    }
	  else // same sign
	    {
	      weight  = -1;
	      float likesign_comb        = (float)nPlusMu*(nPlusMu-1.0)/2.0+(float)nMinusMu*(nMinusMu-1.0)/2.0; // number of combinatorial pairs C(nplus)^2, C(nminus)^2
	    float unlikesign_bkgd_comb = (float)nPlusMu*nMinusMu - (nPlusMu>nMinusMu?nMinusMu:nPlusMu);
	    weight2 = -1.0 * unlikesign_bkgd_comb/likesign_comb;
	    }
	  
	  // dimuon variables
	  TLorentzVector *Reco_QQ      = (TLorentzVector *) Reco_QQ_4mom->At(iQQ);
	  TLorentzVector *Reco_QQ_mupl = (TLorentzVector *) Reco_QQ_mupl_4mom->At(iQQ);
	  TLorentzVector *Reco_QQ_mumi = (TLorentzVector *) Reco_QQ_mumi_4mom->At(iQQ);
	  invariantMass = Reco_QQ->M();
	  upsPt         = Reco_QQ->Pt();
	  upsEta        = Reco_QQ->Eta();
	  upsPhi        = Reco_QQ->Phi();
	  upsRapidity   = Reco_QQ->Rapidity();
	  
	  // single muon variables
	  muMinusPt     = Reco_QQ_mumi->Pt();
	  muMinusEta    = Reco_QQ_mumi->Eta();
	  muMinusPhi    = Reco_QQ_mumi->Phi();
	  muPlusPt      = Reco_QQ_mupl->Pt();
	  muPlusEta     = Reco_QQ_mupl->Eta();
	  muPlusPhi     = Reco_QQ_mupl->Phi();
	  	  
      // apply extra selection

	  bool bProcess = false;
	 
	  if(!bAllTriggers) bProcess = ((HLTriggers&nTriggerBit)==nTriggerBit && (Reco_QQ_trig[iQQ]&nTriggerBit)==nTriggerBit && 
					vProb>newVtxProbCut && muMinusPt>ptcut && muPlusPt>ptcut);
	  // for PbPB sample, the v1 and v2 same trigger is in nNtriggerBit=1 and nTriggerBit=2 respectivelly
	  /* if(isAArereco) bProcess = (( ( (HLTriggers&nTriggerBit)==nTriggerBit && (Reco_QQ_trig[iQQ]&nTriggerBit)==nTriggerBit ) ||  */
	  /* 			       ( (HLTriggers&(nTriggerBit+1))==(nTriggerBit+1) && (Reco_QQ_trig[iQQ]&(nTriggerBit+1))==(nTriggerBit+1) ) ) */
	  /* 			     &&  vProb>newVtxProbCut && muMinusPt>ptcut && muPlusPt>ptcut); */
	  if(isAArereco) bProcess = ( ( (Reco_QQ_trig[iQQ]&nTriggerBit)==nTriggerBit ) 
				     &&  vProb>newVtxProbCut && muMinusPt>ptcut && muPlusPt>ptcut);
	  if (bProcess)
	    {
	      UpsilonTree_allsign->Fill();// all sign and all mass
	      if (Reco_QQ_sign[iQQ]==0) // opposite sign
		{
		  UpsilonTree->Fill();// OS and all mass
		  if (Reco_QQ->M()>mass_min && Reco_QQ->M()<mass_max) 
		    {
		      h_QQ_mass->Fill(Reco_QQ->M());// all upsilons in 7->14
		    }
		  if (Reco_QQ_mupl->Pt()>=ptcut && Reco_QQ_mumi->Pt()>=ptcut && Reco_QQ->M()>mass_min && Reco_QQ->M()<mass_max) 
		    {
		      h_QQ_mass_1->Fill(Reco_QQ->M()); // all OS upsilons in 7->14 and pt_mu>4GeV/c
		    }
		}//opposite sign
	      else // same sign in the acceptance
		if (Reco_QQ_mupl->Pt()>=ptcut && Reco_QQ_mumi->Pt()>=ptcut && Reco_QQ->M()>mass_min && Reco_QQ->M()<mass_max) 
		  {
		    h_QQ_mass_2->Fill(Reco_QQ->M());// all SS upsilons in 7->14 and pt_mu>4GeV/c
		  }
	      

	  //--------------------------------------------------------
	  // %%%%%%%%%% track rotation: redefine some of the variables, the others remain the same
	  Double_t ran = gRandom->Rndm();
	  
	  if(ran < 0.5 ) RmuPlusPhi = muPlusPhi + TMath::Pi();
	  else  RmuMinusPhi = muMinusPhi + TMath::Pi();
	  
	  TLorentzVector mu1;
	  mu1.SetPtEtaPhiM( muPlusPt, muPlusEta, RmuPlusPhi, 0.105);
	  TLorentzVector mu2;
	  mu2.SetPtEtaPhiM( muMinusPt, muMinusEta, RmuMinusPhi, 0.105);
	  
	  TLorentzVector dimuon;
	  dimuon = mu1 + mu2;
	  
	  invariantMass = dimuon.M();
	  upsPt         = dimuon.Pt();
	  upsEta        = dimuon.Eta();
	  upsRapidity   = dimuon.Rapidity();

	  UpsilonTree_trkRot->Fill();
	  
	  /*
	    for (int k=0;k<10;k++){
	    //track rotation
	    //if(ran < 0.5 ) RmuPlusPhi = muPlusPhi + TMath::Pi();
	    //else  RmuMinusPhi = muMinusPhi + TMath::Pi();
	    Double_t ran = gRandom->Rndm();
	    
	    if (k<5) {
	    muMinusPhi = muMinusPhi + TMath::Pi()*ran;
	    }   
	    else {
	    muPlusPhi = muPlusPhi + TMath::Pi()*ran;
	    }   
	    
	    TLorentzVector mu1;
	    mu1.SetPtEtaPhiM( muPlusPt, muPlusEta, muPlusPhi, 0.105);
	    TLorentzVector mu2;
	    mu2.SetPtEtaPhiM( muMinusPt, muMinusEta, muMinusPhi, 0.105);
	    
	    TLorentzVector dimuon;
	    dimuon = mu1 + mu2;
	    
	    invariantMass=dimuon.M();
	    upsPt=dimuon.Pt();
	    upsEta=dimuon.Eta();
	    upsRapidity=dimuon.Rapidity();
	    
	    UpsilonTree_trkRot->Fill();
	    } 
	  */  
	}// if bProcess
      }// for each QQ pair		
    }// for each event in the tree
  
  //  __________________________________________________________________
  // ###### plotting 
  TH1 * phMAxis = new TH1D("phMAxis",";m_{#mu#mu} [GeV/c^{2}];Events/(0.1 GeV/c^{2})",1,mass_min,mass_max);
  phMAxis->SetDirectory(0);
  phMAxis->SetMinimum(1);
  phMAxis->SetMaximum(4e4);

  TCanvas *c1 = new TCanvas("c1","c1");
  phMAxis->Draw();
  
  //h_QQ_mass->GetYaxis()->SetRangeUser(0,180);
  h_QQ_mass->SetMarkerColor(kRed);
  h_QQ_mass->SetMarkerStyle(22);
  h_QQ_mass->Draw("PEsame");
  
  h_QQ_mass_1->SetMarkerColor(kBlue);
  h_QQ_mass_1->SetMarkerStyle(20);
  h_QQ_mass_1->GetXaxis()->CenterTitle(kTRUE);
  h_QQ_mass_1->Draw("PEsame");
  
  //h_QQ_mass_2->GetYaxis()->SetRangeUser(0,400);
  h_QQ_mass_2->SetMarkerColor(kRed);
  h_QQ_mass_2->SetMarkerStyle(24);
  h_QQ_mass_2->GetXaxis()->CenterTitle(kTRUE);
  h_QQ_mass_2->Draw("PEsame");
  
  TLegend *legend = new TLegend(.4,.7,.8,.9);
  legend->SetTextSize(0.025);
  legend->AddEntry(h_QQ_mass,  Form("OS && M [%.1f,%.1f]GeV: %.0f",mass_min, mass_max,h_QQ_mass->Integral(1,nBins)),"P");
  legend->AddEntry(h_QQ_mass_1,Form("OS && M [%.1f,%.1f]GeV && p_{T}>%.1fGeV/c: %.0f",mass_min, mass_max,ptcut,h_QQ_mass_1->Integral(1,nBins)),"P");
  legend->AddEntry(h_QQ_mass_2,Form("SS && M [%.1f,%.1f]GeV && p_{T}>%.1fGeV/c:%.0f",mass_min, mass_max,ptcut,h_QQ_mass_2->Integral(1,nBins)),"P");

  legend->Draw();
  
  h_QQ_mass->Write();
  h_QQ_mass_1->Write();
  h_QQ_mass_2->Write();
  c1->SaveAs(Form("dimuonDistribution_%s_Run%s_trigBit%d_allTriggers%d.pdf",dataSource,runNumber,nTriggerBit,bAllTriggers));
  //  c1->SaveAs(Form("dimuonDistribution_%s_Run%s_trigBit%d_allTriggers%d.gif",dataSource,runNumber,nTriggerBit,bAllTriggers));
  c1->Write();

  f1->Write();
  
}


