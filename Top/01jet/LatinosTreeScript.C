#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>


// Constants
//------------------------------------------------------------------------------
const Bool_t runAtOviedo = true;

const Double_t ZMASS = 91.1876;  // GeV

const UInt_t numberMetCuts = 5;

const Double_t MetCut[numberMetCuts] = {20, 25, 30, 45, 1000};  // GeV

const UInt_t nLevels = 10;

const Int_t HwidthMVA_NoBin = 17 ;// sizeof(HwidthMVA_Bins);
const Double_t HwidthMVA_Bins[HwidthMVA_NoBin+1] = {-1.00, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00};

enum {
  TriggerLevel,
  MetCutLevel,
  LowMinvLevel,
  ZVetoLevel,
  pMetCutLevel,
  JetVetoLevel,
  DeltaPhiJetLevel,
  SoftMuVetoLevel,
  ExtraLeptonLevel,
  TopTaggingLevel
};

const TString sLevel[nLevels] = {
  "Trigger",
  "MetCut",
  "LowMinv",
  "ZVeto",
  "pMetCut",
  "JetVeto",
  "DeltaPhiJet",
  "SoftMuVeto",
  "ExtraLepton",
  "TopTagging"
};


// Member functions
//------------------------------------------------------------------------------
void    FillHistograms(UInt_t level,
		UInt_t check);

// Data members
//------------------------------------------------------------------------------
TH1F*   hWeff           [nLevels];
TH1F*   hW              [nLevels];
TH1F*   hPtLepton1      [nLevels];
TH1F*   hPtLepton2      [nLevels];
TH1F*   hBDT1           [nLevels];
TH1F*   hBDT2           [nLevels];
TH1F*   hHwidthMVA      [nLevels];
TH1F*   hPtDiLepton     [nLevels];
TH1F*   hMinv           [nLevels];
TH1F*   hMt             [nLevels];
TH1F*   hNJets30        [nLevels];
TH1F*   hpfMet          [nLevels];
TH1F*   hppfMet         [nLevels];
TH1F*   hchMet          [nLevels];
TH1F*   hpchMet         [nLevels];
TH1F*   hpminMet        [nLevels];
TH1F*   hDeltaRLeptons  [nLevels];
TH1F*   hDeltaPhiLeptons[nLevels];
TH1F*   hDPhiPtllJet    [nLevels];

TH1F*   hWeff_NM1           [nLevels];
TH1F*   hW_NM1              [nLevels];
TH1F*   hPtLepton1_NM1      [nLevels];
TH1F*   hPtLepton2_NM1      [nLevels];
TH1F*   hBDT1_NM1           [nLevels];
TH1F*   hBDT2_NM1           [nLevels];
TH1F*   hHwidthMVA_NM1      [nLevels];
TH1F*   hPtDiLepton_NM1     [nLevels];
TH1F*   hMinv_NM1           [nLevels];
TH1F*   hMt_NM1             [nLevels];
TH1F*   hNJets30_NM1        [nLevels];
TH1F*   hpfMet_NM1          [nLevels];
TH1F*   hppfMet_NM1         [nLevels];
TH1F*   hchMet_NM1          [nLevels];
TH1F*   hpchMet_NM1         [nLevels];
TH1F*   hpminMet_NM1        [nLevels];
TH1F*   hDeltaRLeptons_NM1  [nLevels];
TH1F*   hDeltaPhiLeptons_NM1[nLevels];
TH1F*   hDPhiPtllJet_NM1    [nLevels];

Float_t efficiencyW;
Float_t totalW;


// Tree variables
//------------------------------------------------------------------------------
Float_t baseW;
Float_t channel;
Float_t chmet;
Float_t dataset;
Float_t dphill;
Float_t dphilljet;
Float_t dphilljetjet;
Float_t drll;
Float_t effW;
Float_t fakeW;
Float_t eta1;
Float_t eta2;
Float_t detajj;
Float_t jeteta1;
Float_t jeteta2;
Float_t jeteta3;
Float_t jetpt1;
Float_t jetpt2;
Float_t jetpt3;
Float_t jettche1;
Float_t jettche2;
Float_t mctruth;
Float_t mll;
Float_t mjj;
Float_t mpmet;
Float_t mth;
Float_t nbjet;
Float_t nbjettche;
Float_t nextra;
Float_t njet;
Float_t nvtx;
Float_t pchmet;
Float_t pfmet;
Float_t ppfmet;
Float_t pt1;
Float_t pt2;
Float_t bdt1;
Float_t bdt2;
Float_t puW;
Float_t ch1;
Float_t ch2;
Float_t ptll;
Float_t softtche;
Float_t hardtche;
Float_t trigger;
Float_t triggW;
Int_t   bveto;
Int_t   bveto_ip;
Int_t   bveto_mu;
Int_t   bveto_nj30;
Int_t   dphiveto;
Int_t   sameflav;
Int_t   zveto;
UInt_t  event;
UInt_t  lumi;
UInt_t  run;
Float_t jetmva1;
Float_t HwidthMVAbkg;


//------------------------------------------------------------------------------
// LatinosTreeScript
//------------------------------------------------------------------------------
void LatinosTreeScript(Float_t luminosity,
                       Int_t   jetChannel,
		       TString flavorChannel,
		       TString theSample,
		       Bool_t  verbose)
{
  TH1::SetDefaultSumw2();
  
  TString path = Form("rootfiles/%djet/%s/", jetChannel, flavorChannel.Data());
  
  gSystem->mkdir(path, kTRUE);
  
  TFile* output = new TFile(path + theSample + ".root", "recreate");
  
  // Histograms
  //----------------------------------------------------------------------------
  for (UInt_t i=0; i<nLevels; i++) {
    hWeff           [i] = new TH1F("hWeff"            + sLevel[i], "",   3, 0,   3);
    hW              [i] = new TH1F("hW"               + sLevel[i], "",   3, 0,   3);
    hPtLepton1      [i] = new TH1F("hPtLepton1"       + sLevel[i], "", 200, 0, 200);
    hPtLepton2      [i] = new TH1F("hPtLepton2"       + sLevel[i], "", 200, 0, 200);
    hBDT1           [i] = new TH1F("hBDT1"            + sLevel[i], "", 100, 0,   1);
    hBDT2           [i] = new TH1F("hBDT2"            + sLevel[i], "", 100, 0,   1);
    hHwidthMVA      [i] = new TH1F("hHwidthMVA"       + sLevel[i], "", HwidthMVA_NoBin,HwidthMVA_Bins);
    hPtDiLepton     [i] = new TH1F("hPtDiLepton"      + sLevel[i], "", 200, 0, 200);
    hMinv           [i] = new TH1F("hMinv"            + sLevel[i], "", 200, 0, 200);
    hMt             [i] = new TH1F("hMt"              + sLevel[i], "", 250, 0, 250);
    hNJets30        [i] = new TH1F("hNJetsPF30"       + sLevel[i], "",  10, 0,  10);
    hpfMet          [i] = new TH1F("hpfMet"           + sLevel[i], "", 150, 0, 150);
    hppfMet         [i] = new TH1F("hppfMet"          + sLevel[i], "", 150, 0, 150);
    hchMet          [i] = new TH1F("hchMet"           + sLevel[i], "", 150, 0, 150);
    hpchMet         [i] = new TH1F("hpchMet"          + sLevel[i], "", 150, 0, 150);
    hpminMet        [i] = new TH1F("hpminMet"         + sLevel[i], "", 150, 0, 150);
    hDeltaRLeptons  [i] = new TH1F("hDeltaRLeptons"   + sLevel[i], "",  50, 0,   5);
    hDeltaPhiLeptons[i] = new TH1F("hDeltaPhiLeptons" + sLevel[i], "",  32, 0, 3.2);
    hDPhiPtllJet    [i] = new TH1F("hDPhiPtllJet"     + sLevel[i], "",  32, 0, 3.2);
    
    hWeff_NM1           [i] = new TH1F("hWeff"            + sLevel[i] + "_NM1", "",   3, 0,   3);
    hW_NM1              [i] = new TH1F("hW"               + sLevel[i] + "_NM1", "",   3, 0,   3);
    hPtLepton1_NM1      [i] = new TH1F("hPtLepton1"       + sLevel[i] + "_NM1", "", 200, 0, 200);
    hPtLepton2_NM1      [i] = new TH1F("hPtLepton2"       + sLevel[i] + "_NM1", "", 200, 0, 200);
    hBDT1_NM1           [i] = new TH1F("hBDT1"            + sLevel[i] + "_NM1", "", 100, 0,   1);
    hBDT2_NM1           [i] = new TH1F("hBDT2"            + sLevel[i] + "_NM1", "", 100, 0,   1);
    hHwidthMVA_NM1      [i] = new TH1F("hHwidthMVA"       + sLevel[i], "", HwidthMVA_NoBin,HwidthMVA_Bins);
    hPtDiLepton_NM1     [i] = new TH1F("hPtDiLepton"      + sLevel[i] + "_NM1", "", 200, 0, 200);
    hMinv_NM1           [i] = new TH1F("hMinv"            + sLevel[i] + "_NM1", "", 200, 0, 200);
    hMt_NM1             [i] = new TH1F("hMt"              + sLevel[i] + "_NM1", "", 250, 0, 250);
    hNJets30_NM1        [i] = new TH1F("hNJetsPF30"       + sLevel[i] + "_NM1", "",  10, 0,  10);
    hpfMet_NM1          [i] = new TH1F("hpfMet"           + sLevel[i] + "_NM1", "", 150, 0, 150);
    hppfMet_NM1         [i] = new TH1F("hppfMet"          + sLevel[i] + "_NM1", "", 150, 0, 150);
    hchMet_NM1          [i] = new TH1F("hchMet"           + sLevel[i] + "_NM1", "", 150, 0, 150);
    hpchMet_NM1         [i] = new TH1F("hpchMet"          + sLevel[i] + "_NM1", "", 150, 0, 150);
    hpminMet_NM1        [i] = new TH1F("hpminMet"         + sLevel[i] + "_NM1", "", 150, 0, 150);
    hDeltaRLeptons_NM1  [i] = new TH1F("hDeltaRLeptons"   + sLevel[i] + "_NM1", "",  50, 0,   5);
    hDeltaPhiLeptons_NM1[i] = new TH1F("hDeltaPhiLeptons" + sLevel[i] + "_NM1", "",  32, 0, 3.2);
    hDPhiPtllJet_NM1    [i] = new TH1F("hDPhiPtllJet"     + sLevel[i] + "_NM1", "",  32, 0, 3.2);
  }


  // Data-driven methods: Z+jets
  //----------------------------------------------------------------------------
  TH1F* hNinZevents     [numberMetCuts];
  TH1F* hNoutZevents    [numberMetCuts];
  TH1F* hNinLooseZevents[numberMetCuts];
  TH1F* hMassInZevents  [numberMetCuts];
  TH1F* hMassOutZevents [numberMetCuts];

  for (size_t nC=0; nC<numberMetCuts; nC++) {
    hNinZevents     [nC] = new TH1F(Form("hNinZevents%.2f",      MetCut[nC]), "",   3, 0,   3);
    hNoutZevents    [nC] = new TH1F(Form("hNoutZevents%.2f",     MetCut[nC]), "",   3, 0,   3);
    hNinLooseZevents[nC] = new TH1F(Form("hNinLooseZevents%.2f", MetCut[nC]), "",   3, 0,   3);
    hMassInZevents  [nC] = new TH1F(Form("hMassInZevents%.2f",   MetCut[nC]), "", 200, 0, 200);
    hMassOutZevents [nC] = new TH1F(Form("hMassOutZevents%.2f",  MetCut[nC]), "", 200, 0, 200);
  }


  // Data-driven methods: Top
  //----------------------------------------------------------------------------
  TH1F* hTopTaggedEvents            = new TH1F("hTopTaggedEvents",            "", 3, 0, 3);
  TH1F* hNTopControlRegion          = new TH1F("hNTopControlRegion",          "", 3, 0, 3);
  TH1F* hNTopTaggedTopControlRegion = new TH1F("hNTopTaggedTopControlRegion", "", 3, 0, 3);
  
  TH1F* hbTagDisTopTaggedEvents            = new TH1F("hbTagDisTopTaggedEvents",            "", 300, -10, 20);
  TH1F* hbTagDisNTopControlRegion          = new TH1F("hbTagDisNTopControlRegion",          "", 300, -10, 20);
  TH1F* hbTagDisNTopTaggedTopControlRegion = new TH1F("hbTagDisNTopTaggedTopControlRegion", "", 300, -10, 20);
  
  //TopBkgr estimation control plots
  TH1F* hbTagDisTopTaggedEvents_channel  = new TH1F("hbTagDisTopTaggedEvents_channel", "", 4,  0,  4);
  TH1F* hbTagDisTopTaggedEvents_dphill   = new TH1F("hbTagDisTopTaggedEvents_dphill",  "", 18, 0,  180);
  TH1F* hbTagDisTopTaggedEvents_jetpt1   = new TH1F("hbTagDisTopTaggedEvents_jetpt1",  "", 40, 0,  300);
  TH1F* hbTagDisTopTaggedEvents_jetmva1  = new TH1F("hbTagDisTopTaggedEvents_jetmva1", "", 10, -1, 1);
  TH1F* hbTagDisTopTaggedEvents_softtche = new TH1F("hbTagDisTopTaggedEvents_softtche","", 25, -5, 20);
  TH1F* hbTagDisTopTaggedEvents_hardtche = new TH1F("hbTagDisTopTaggedEvents_hardtche","", 25, -5, 20);
  TH1F* hbTagDisTopTaggedEvents_mll      = new TH1F("hbTagDisTopTaggedEvents_mll",     "", 20, 0,  200);
  TH1F* hbTagDisTopTaggedEvents_mth      = new TH1F("hbTagDisTopTaggedEvents_mth",     "", 30, 0,  300);
  TH1F* hbTagDisTopTaggedEvents_bdt1     = new TH1F("hbTagDisTopTaggedEvents_bdt1",    "", 100, 0,  1);
  TH1F* hbTagDisTopTaggedEvents_bdt2     = new TH1F("hbTagDisTopTaggedEvents_bdt2",    "", 100, 0,  1);
  TH1F* hbTagDisTopTaggedEvents_HwidthMVAbkg     = new TH1F("hbTagDisTopTaggedEvents_HwidthMVAbkg",    "", HwidthMVA_NoBin,HwidthMVA_Bins);
  
  TH1F* hbTagDisNTopControlRegion_channel  = new TH1F("hbTagDisNTopControlRegion_channel", "", 4,  0,  4);
  TH1F* hbTagDisNTopControlRegion_dphill   = new TH1F("hbTagDisNTopControlRegion_dphill",  "", 18, 0,  180);
  TH1F* hbTagDisNTopControlRegion_jetpt1   = new TH1F("hbTagDisNTopControlRegion_jetpt1",  "", 40, 0,  300);
  TH1F* hbTagDisNTopControlRegion_jetmva1  = new TH1F("hbTagDisNTopControlRegion_jetmva1", "", 10, -1, 1);
  TH1F* hbTagDisNTopControlRegion_softtche = new TH1F("hbTagDisNTopControlRegion_softtche","", 25, -5, 20);
  TH1F* hbTagDisNTopControlRegion_hardtche = new TH1F("hbTagDisNTopControlRegion_hardtche","", 25, -5, 20);
  TH1F* hbTagDisNTopControlRegion_mll      = new TH1F("hbTagDisNTopControlRegion_mll",     "", 20, 0,  200);
  TH1F* hbTagDisNTopControlRegion_mth      = new TH1F("hbTagDisNTopControlRegion_mth",     "", 30, 0,  300);
  TH1F* hbTagDisNTopControlRegion_bdt1     = new TH1F("hbTagDisNTopControlRegion_bdt1",    "", 100, 0,  1);
  TH1F* hbTagDisNTopControlRegion_bdt2     = new TH1F("hbTagDisNTopControlRegion_bdt2",    "", 100, 0,  1);
  TH1F* hbTagDisNTopControlRegion_HwidthMVAbkg     = new TH1F("hbTagDisNTopControlRegion_HwidthMVAbkg",    "", HwidthMVA_NoBin,HwidthMVA_Bins);
  
  TH1F* hbTagDisNTopTaggedTopControlRegion_channel  = new TH1F("hbTagDisNTopTaggedTopControlRegion_channel", "", 4,  0,  4);
  TH1F* hbTagDisNTopTaggedTopControlRegion_dphill   = new TH1F("hbTagDisNTopTaggedTopControlRegion_dphill",  "", 18, 0,  180);
  TH1F* hbTagDisNTopTaggedTopControlRegion_jetpt1   = new TH1F("hbTagDisNTopTaggedTopControlRegion_jetpt1",  "", 40, 0,  300);
  TH1F* hbTagDisNTopTaggedTopControlRegion_jetmva1  = new TH1F("hbTagDisNTopTaggedTopControlRegion_jetmva1", "", 10, -1, 1);
  TH1F* hbTagDisNTopTaggedTopControlRegion_softtche = new TH1F("hbTagDisNTopTaggedTopControlRegion_softtche","", 25, -5, 20);
  TH1F* hbTagDisNTopTaggedTopControlRegion_hardtche = new TH1F("hbTagDisNTopTaggedTopControlRegion_hardtche","", 25, -5, 20);
  TH1F* hbTagDisNTopTaggedTopControlRegion_mll      = new TH1F("hbTagDisNTopTaggedTopControlRegion_mll",     "", 20, 0,  200);
  TH1F* hbTagDisNTopTaggedTopControlRegion_mth      = new TH1F("hbTagDisNTopTaggedTopControlRegion_mth",     "", 30, 0,  300);
  TH1F* hbTagDisNTopTaggedTopControlRegion_bdt1     = new TH1F("hbTagDisNTopTaggedTopControlRegion_bdt1",    "", 100, 0,  1);
  TH1F* hbTagDisNTopTaggedTopControlRegion_bdt2     = new TH1F("hbTagDisNTopTaggedTopControlRegion_bdt2",    "", 100, 0,  1);
  TH1F* hbTagDisNTopTaggedTopControlRegion_HwidthMVAbkg     = new TH1F("hbTagDisNTopTaggedTopControlRegion_HwidthMVAbkg",    "", HwidthMVA_NoBin,HwidthMVA_Bins);
  
  //----------------------------------------------------------------------------
  // Input files
  //----------------------------------------------------------------------------
  TString filesPath;
  
  if (runAtOviedo) filesPath = "/afs/cern.ch/work/x/xjanssen/public/LatinoTrees/R53X_S1_V08_S2_V09_S3_V13/MoriondeffWPuWtriggW/MC_TightTight_DABCABC/";
  else             filesPath = "/afs/cern.ch/work/x/xjanssen/public/LatinoTrees/R53X_S1_V08_S2_V09_S3_V13/MoriondeffWPuWtriggW/MC_TightTight_DABCABC/";

  TChain* tree = new TChain("latino", "latino");

  if (theSample == "DataRun2012_Total") {
    tree->Add("/afs/cern.ch/work/j/jfernan2/public/LatinoTrees/R53X_S1_V08_S2_V09_S3_V13/Data_TightTight/4L/latino_RunA_892pbinv.root");
    tree->Add("/afs/cern.ch/work/j/jfernan2/public/LatinoTrees/R53X_S1_V08_S2_V09_S3_V13/Data_TightTight/4L/latino_RunB_4404pbinv.root");
    tree->Add("/afs/cern.ch/work/j/jfernan2/public/LatinoTrees/R53X_S1_V08_S2_V09_S3_V13/Data_TightTight/4L/latino_RunC_7032pbinv.root");
    tree->Add("/afs/cern.ch/work/j/jfernan2/public/LatinoTrees/R53X_S1_V08_S2_V09_S3_V13/Data_TightTight/4L/latino_RunD_7274pbinv.root");
  }
  else if (theSample == "WJetsFakes_Total") {
    tree->Add("/afs/cern.ch/work/c/calderon/public/LatinoTrees/Moriond13/LooseLoose/latino_LooseLoose_19.5fb.root");
  }
  else if (theSample == "ggWWto2L") {
    tree->Add(filesPath + "latino_001_GluGluToWWTo4L.root");
  }
  else if (theSample == "WWTo2L2Nu") {
    tree->Add(filesPath + "latino_000_WWJets2LMad.root");
  }
  else if (theSample == "WZ") {
    tree->Add(filesPath + "latino_074_WZJetsMad.root");
    tree->Add(filesPath + "latino_078_WZTo2L2QMad.root");
  }
  else if (theSample == "ZZ") {
    tree->Add(filesPath + "latino_075_ZZJetsMad.root");
    tree->Add(filesPath + "latino_079_ZZTo2L2QMad.root"); 
  }
  else if (theSample == "TTbar") {
    tree->Add(filesPath + "latino_019_TTTo2L2Nu2B.root");
  }
  else if (theSample == "TW") {
    tree->Add(filesPath + "latino_011_TtWFullDR.root");
    tree->Add(filesPath + "latino_012_TbartWFullDR.root");
  }
  else if (theSample == "DY") {
    tree->Add(filesPath + "latino_036_DY10toLLMad.root");
    tree->Add(filesPath + "latino_037_DY50toLLMad.root");
  }
  else if (theSample == "DYtautau") {
    tree->Add("/afs/cern.ch/work/x/xjanssen/public/LatinoTrees/R53X_S1_V08_S2_V09_S3_V13/MoriondeffWPuWtriggW/TauTau_btagvar/latino_DYtt_19.5fb.root");
  }
  else if (theSample == "WgammaNoStar") {
    tree->Add(filesPath + "latino_085_WgammaToLNuG.root");
    tree->Add(filesPath + "latino_008_GamGamWW.root");
  }
  else if (theSample == "WgammaStar") {
    tree->Add(filesPath + "latino_082_WGstarToElNuMad.root");
    tree->Add(filesPath + "latino_083_WGstarToMuNuMad.root");
    tree->Add(filesPath + "latino_084_WGstarToTauNuMad.root");
  }
  else if (theSample == "HWW125") { 
    tree->Add(filesPath + "latino_1125_ggToH125toWWTo2LAndTau2Nu.root");
    tree->Add(filesPath + "latino_2125_vbfToH125toWWTo2LAndTau2Nu.root");
  }
  else if (theSample == "Zgamma") { 
    tree->Add(filesPath +"latino_086_ZgammaToLLuG.root");
  }
  else {
    return;
  }

  // Declaration of leaf types
  //----------------------------------------------------------------------------
  tree->SetBranchAddress("baseW",        &baseW);
  tree->SetBranchAddress("channel",      &channel);
  tree->SetBranchAddress("chmet",        &chmet);
  tree->SetBranchAddress("dataset",      &dataset);
  tree->SetBranchAddress("dphill",       &dphill);
  tree->SetBranchAddress("dphilljet",    &dphilljet);
  tree->SetBranchAddress("dphilljetjet", &dphilljetjet);
  tree->SetBranchAddress("drll",         &drll);
  tree->SetBranchAddress("effW",         &effW);
  tree->SetBranchAddress("eta1",         &eta1);
  tree->SetBranchAddress("eta2",         &eta2);
  tree->SetBranchAddress("detajj",       &detajj);
  tree->SetBranchAddress("jeteta1",      &jeteta1);
  tree->SetBranchAddress("jeteta2",      &jeteta2);
  tree->SetBranchAddress("jeteta3",      &jeteta3);
  tree->SetBranchAddress("jetpt1",       &jetpt1);
  tree->SetBranchAddress("jetpt2",       &jetpt2);
  tree->SetBranchAddress("jetpt3",       &jetpt3);
  tree->SetBranchAddress("jettche1",     &jettche1);
  tree->SetBranchAddress("jettche2",     &jettche2);
  tree->SetBranchAddress("mctruth",      &mctruth);
  tree->SetBranchAddress("mjj",          &mjj);
  tree->SetBranchAddress("mll",          &mll);
  tree->SetBranchAddress("mpmet",        &mpmet); 
  tree->SetBranchAddress("mth",          &mth);
  tree->SetBranchAddress("nbjet",        &nbjet);
  tree->SetBranchAddress("nbjettche",    &nbjettche);
  tree->SetBranchAddress("nextra",       &nextra);
  tree->SetBranchAddress("njet",         &njet);
  tree->SetBranchAddress("nvtx",         &nvtx);
  tree->SetBranchAddress("pchmet",       &pchmet);
  tree->SetBranchAddress("pfmet",        &pfmet);
  tree->SetBranchAddress("ppfmet",       &ppfmet);
  tree->SetBranchAddress("pt1",          &pt1);
  tree->SetBranchAddress("pt2",          &pt2);
  tree->SetBranchAddress("bdt1",         &bdt1);
  tree->SetBranchAddress("bdt2",         &bdt2);
  tree->SetBranchAddress("ch1",          &ch1);
  tree->SetBranchAddress("ch2",          &ch2);
  tree->SetBranchAddress("ptll",         &ptll);
  tree->SetBranchAddress("hardtche",     &hardtche);
  tree->SetBranchAddress("softtche",     &softtche);
  tree->SetBranchAddress("trigger",      &trigger);
  tree->SetBranchAddress("triggW",       &triggW);
  tree->SetBranchAddress("bveto",        &bveto);
  tree->SetBranchAddress("bveto_ip",     &bveto_ip);
  tree->SetBranchAddress("bveto_mu",     &bveto_mu);
  tree->SetBranchAddress("bveto_nj30",   &bveto_nj30);
  tree->SetBranchAddress("dphiveto",     &dphiveto);
  tree->SetBranchAddress("sameflav",     &sameflav);
  tree->SetBranchAddress("zveto",        &zveto);
  tree->SetBranchAddress("event",        &event);
  tree->SetBranchAddress("lumi",         &lumi);
  tree->SetBranchAddress("run",          &run);
  tree->SetBranchAddress("jetmva1",      &jetmva1);
  tree->SetBranchAddress("HwidthMVAbkg", &HwidthMVAbkg);
  
  if (theSample.Contains("WJetsFakes"))
    tree->SetBranchAddress("fakeW", &fakeW);

  if (!theSample.Contains("WJetsFakes") && !theSample.Contains("Data"))
    tree->SetBranchAddress("puW", &puW);

  // Set the channel
  //----------------------------------------------------------------------------
  Float_t SelectedChannel = -999;
  
  if      (flavorChannel == "MuMu") SelectedChannel =  0;
  else if (flavorChannel == "EE"  ) SelectedChannel =  1;
  else if (flavorChannel == "EMu" ) SelectedChannel =  2;
  else if (flavorChannel == "MuE" ) SelectedChannel =  3;
  else if (flavorChannel == "All" ) SelectedChannel = -1;
  
  //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------
  for (int ievent=0; ievent<tree->GetEntries(); ievent++) {
    tree->GetEntry(ievent);

    efficiencyW = effW * triggW;
    totalW      = -999;

    if (theSample.Contains("Data"))
    {
      totalW = 1.0;
    }
    else if (theSample.Contains("WJetsFakes"))
    {
      totalW = fakeW;
    }
    else
    {
      efficiencyW = puW * effW * triggW;

      totalW = (1 + 0.6 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * luminosity;
    }

    // Help variables
    //--------------------------------------------------------------------------
    Int_t dphiv = ((njet <= 1 && dphiveto) || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));

    Float_t metvar = (njet <= 1) ? mpmet : pfmet;

    Float_t jetbin = njet;

    if (njet == 3) jetbin = 2;

    Int_t vbfsel = ((jetpt3 <= 30 || !(jetpt3 > 30 && ((jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))));
   
    //Hwidth Selection
    //--------------------------------------------------------------------------
    //https://github.com/latinos/HWWAnalysis/blob/master/ShapeAnalysis/python/hwwinfo.py#L952
    Int_t HwidthSelection = ( (njet==0 && ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>40) )) ||
	                      ((njet==1) && ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>40) )) ||
			      ((njet>=2) && ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && ( (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0))))) &&
				abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5 && detajj>2.5 && mjj>500 && ( (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20) ) )));

    //HwidthmthmllSelection selection
    //https://github.com/latinos/HWWAnalysis/blob/master/ShapeAnalysis/python/hwwinfo.py#L965
//    Int_t HwidthSelection = ( (njet==0  && ( (mll<70 && pt2<35 && ptll>45) || (mll>=70 && pt2>20 && pt1>50) )) ||
//                            ((njet==1) && ( (mll<70 && pt2<35 && ptll>45) || (mll>=70 && pt2>20 && pt1>50) )) ||
//                            ((njet>=2) && ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && ( (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0))))) &&
//                              abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5 && detajj>2.5 && mjj>500 && ( (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20) ) )));

    // The selection begins here
    //--------------------------------------------------------------------------
    Bool_t accept_WGstar = (chmet < (0.75*pt1+100) && chmet < (0.75*jetpt1+100));

    if ((dataset == 36 || dataset == 37) && mctruth == 2)                    continue;
    if (dataset == 82 && !accept_WGstar)                                     continue;
    if (dataset == 86 && (flavorChannel == "MuMu" || flavorChannel == "EE")) continue;
    if ((SelectedChannel != -1) && (channel != SelectedChannel))             continue;
    if (run == 201191)                                                       continue;
    if (!trigger)                                                            continue;
    if (pt1 <= 20)                                                           continue;
    if (pt2 <= 10)                                                           continue;
    if (ch1*ch2 > 0)                                                         continue;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Data-driven methods
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (pfmet > 20 && mpmet > 20 && mll > 12 && nextra == 0 && (dphiv || (!sameflav && HwidthSelection))) {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //
      //Z+jets
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (dphiv && jetbin == jetChannel && bveto_mu && (bveto_ip && (nbjettche == 0 || njet > 3))) {
			
	// Loop over the metvar bins
	//----------------------------------------------------------------------
	for (size_t mc=0; mc<numberMetCuts; mc ++) {
	  if (metvar > MetCut[mc] && fabs(mll - ZMASS) < 7.5) {
	    hNinLooseZevents[mc]->Fill(1,totalW);
	  }
	  
	  if (metvar > MetCut[mc] && metvar < MetCut[mc+1]) {   
	    if (fabs(mll - ZMASS) < 7.5) {
	      hNinZevents[mc]   ->Fill(  1, totalW);
	      
	      hMassInZevents[mc]->Fill(mll, totalW);
	    }
	    else if (fabs(mll - ZMASS) > 15) {
	      hNoutZevents[mc]   ->Fill(  1, totalW);
	      hMassOutZevents[mc]->Fill(mll, totalW);
	    }
	  }
	}
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //
      // Top
      //
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (zveto && metvar > (20 + 25*sameflav)) {
	
	// btag_eff denominator
	if ((jetChannel == 0 && njet == 1 && nbjet == 1) ||
	    (jetChannel == 1 && njet == 2 && jettche2 > 2.1) ||
	    (jetChannel == 2)) {
	  hNTopControlRegion                     -> Fill(1,                        totalW);
	  hbTagDisNTopControlRegion              -> Fill(jettche2,                 totalW);
	  hbTagDisNTopControlRegion_dphill       -> Fill(dphill*TMath::RadToDeg(), totalW);
	  hbTagDisNTopControlRegion_channel      -> Fill(channel,                  totalW);
	  hbTagDisNTopControlRegion_jetpt1       -> Fill(jetpt1,                   totalW);
	  hbTagDisNTopControlRegion_jetmva1      -> Fill(jetmva1,                  totalW);
	  hbTagDisNTopControlRegion_softtche     -> Fill(softtche,                 totalW);
	  hbTagDisNTopControlRegion_hardtche     -> Fill(hardtche,                 totalW);
	  hbTagDisNTopControlRegion_mll          -> Fill(mll,                      totalW);
	  hbTagDisNTopControlRegion_mth          -> Fill(mth,                      totalW);
	  hbTagDisNTopControlRegion_bdt1         -> Fill(bdt1,                     totalW);
	  hbTagDisNTopControlRegion_bdt2         -> Fill(bdt2,                     totalW);
	  hbTagDisNTopControlRegion_HwidthMVAbkg -> Fill(HwidthMVAbkg,             totalW);

	  // btag_eff numerator
	  if ((jetChannel == 0 && !bveto_nj30) ||
	      (jetChannel == 1 && jettche1 > 2.1) ||
	      (jetChannel == 2)) {
	    
	    hNTopTaggedTopControlRegion                     -> Fill(1,                        totalW);
	    hbTagDisNTopTaggedTopControlRegion              -> Fill(jettche2,                 totalW);
	    hbTagDisNTopTaggedTopControlRegion_dphill       -> Fill(dphill*TMath::RadToDeg(), totalW);
	    hbTagDisNTopTaggedTopControlRegion_channel      -> Fill(channel,                  totalW);
	    hbTagDisNTopTaggedTopControlRegion_jetpt1       -> Fill(jetpt1,                   totalW);
	    hbTagDisNTopTaggedTopControlRegion_jetmva1      -> Fill(jetmva1,                  totalW);
	    hbTagDisNTopTaggedTopControlRegion_softtche     -> Fill(softtche,                 totalW);
	    hbTagDisNTopTaggedTopControlRegion_hardtche     -> Fill(hardtche,                 totalW);
	    hbTagDisNTopTaggedTopControlRegion_mll          -> Fill(mll,                      totalW);
	    hbTagDisNTopTaggedTopControlRegion_mth          -> Fill(mth,                      totalW);
	    hbTagDisNTopTaggedTopControlRegion_bdt1         -> Fill(bdt1,                     totalW);
	    hbTagDisNTopTaggedTopControlRegion_bdt2         -> Fill(bdt2,                     totalW);
	    hbTagDisNTopTaggedTopControlRegion_HwidthMVAbkg -> Fill(HwidthMVAbkg,             totalW);
	  }
	}
      }
      
      // Top-tagged events for ttbar estimation
      //----------------------------------------------------------------------~~
      if (zveto && metvar > (20 + 25*sameflav)) {
	
	if ((jetChannel == 0 && njet == 0 && !bveto) ||
	    (jetChannel == 1 && njet == 1 && bveto && jettche1 > 2.1) ||
	    (jetChannel == 2)) {
	  
	  hTopTaggedEvents                     -> Fill(1,                        totalW);
	  hbTagDisTopTaggedEvents              -> Fill(jettche2,                 totalW);
	  hbTagDisTopTaggedEvents_dphill       -> Fill(dphill*TMath::RadToDeg(), totalW);
	  hbTagDisTopTaggedEvents_channel      -> Fill(channel,                  totalW);
	  hbTagDisTopTaggedEvents_jetpt1       -> Fill(jetpt1,                   totalW);
	  hbTagDisTopTaggedEvents_jetmva1      -> Fill(jetmva1,                  totalW);
	  hbTagDisTopTaggedEvents_softtche     -> Fill(softtche,                 totalW);
	  hbTagDisTopTaggedEvents_hardtche     -> Fill(hardtche,                 totalW);
	  hbTagDisTopTaggedEvents_mll          -> Fill(mll,                      totalW);
	  hbTagDisTopTaggedEvents_mth          -> Fill(mth,                      totalW);
	  hbTagDisTopTaggedEvents_bdt1         -> Fill(bdt1,                     totalW);
	  hbTagDisTopTaggedEvents_bdt2         -> Fill(bdt2,                     totalW);
	  hbTagDisTopTaggedEvents_HwidthMVAbkg -> Fill(HwidthMVAbkg,             totalW);
	}
      }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Main analysis
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (UInt_t ilevel=0; ilevel<nLevels; ilevel++) {
      
      FillHistograms(TriggerLevel, ilevel);

      if ((pfmet > 20 && mpmet > 20) || (ilevel == MetCutLevel)) {

	FillHistograms(MetCutLevel, ilevel);

	if ((mll > 12) || (ilevel == LowMinvLevel)) {

	  FillHistograms(LowMinvLevel, ilevel);

	  if (zveto || (ilevel == ZVetoLevel)) {

	    FillHistograms(ZVetoLevel, ilevel);

	    if ((metvar > (20 + 25*sameflav)) || (ilevel == pMetCutLevel)) {

	      FillHistograms(pMetCutLevel, ilevel);

	      if ((jetbin == jetChannel) || (ilevel == JetVetoLevel)) {

		FillHistograms(JetVetoLevel, ilevel);

		if (dphiv || (!sameflav && HwidthSelection) || (ilevel == DeltaPhiJetLevel)) {

		  FillHistograms(DeltaPhiJetLevel, ilevel);

		  if (bveto_mu || (ilevel == SoftMuVetoLevel)) {

		    FillHistograms(SoftMuVetoLevel, ilevel);

		    if ((nextra == 0) || (ilevel == ExtraLeptonLevel)) {

		      FillHistograms(ExtraLeptonLevel, ilevel);

		      Bool_t top_veto = (bveto_ip && (nbjettche == 0 || njet > 3) && (njet <= 1 || vbfsel));

		      if (top_veto || (ilevel == TopTaggingLevel)) {

			FillHistograms(TopTaggingLevel, ilevel);

			hWeff_NM1           [ilevel]->Fill(1,            efficiencyW);
			hW_NM1              [ilevel]->Fill(1,            totalW);
			hPtLepton1_NM1      [ilevel]->Fill(pt1,          totalW);
			hPtLepton2_NM1      [ilevel]->Fill(pt2,          totalW);
			hBDT1_NM1           [ilevel]->Fill(bdt1,         totalW);
			hBDT2_NM1           [ilevel]->Fill(bdt2,         totalW);
			hPtDiLepton_NM1     [ilevel]->Fill(ptll,         totalW);
			hMinv_NM1           [ilevel]->Fill(mll,          totalW);
			hMt_NM1             [ilevel]->Fill(mth,          totalW);
			hNJets30_NM1        [ilevel]->Fill(njet,         totalW);
			hpfMet_NM1          [ilevel]->Fill(pfmet,        totalW);
			hppfMet_NM1         [ilevel]->Fill(ppfmet,       totalW);
			hchMet_NM1          [ilevel]->Fill(chmet,        totalW);
			hpchMet_NM1         [ilevel]->Fill(pchmet,       totalW);
			hpminMet_NM1        [ilevel]->Fill(mpmet,        totalW);
			hDeltaRLeptons_NM1  [ilevel]->Fill(drll,         totalW);
			hDeltaPhiLeptons_NM1[ilevel]->Fill(dphill,       totalW);
			hDPhiPtllJet_NM1    [ilevel]->Fill(dphilljet,    totalW);
			hHwidthMVA_NM1      [ilevel]->Fill(HwidthMVAbkg, totalW);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }


  // Print
  //----------------------------------------------------------------------------
  if (verbose) {
    
    printf("\n Expected number of RAW events for %s\n", theSample.Data());
    printf(" -------------------+-----------\n");
    
    for (UInt_t i=0; i<nLevels; i++)
      printf(" %18s | %.0f\n", sLevel[i].Data(), hW[i]->GetEntries());
  
    printf("\n");

    if (!theSample.Contains("Data")) {

      printf("\n Normalized to %.3f 1/fb\n", luminosity);
      printf(" -------------------+-----------\n");

      for (UInt_t i=0; i<nLevels; i++)
	printf(" %18s | %.0f\n", sLevel[i].Data(), hW[i]->GetSumOfWeights());

      printf("\n");
    }
  }
  
  
  // Save the histograms
  //----------------------------------------------------------------------------
  output->cd();
  output->Write("", TObject::kOverwrite);
  output->Close();
}


//------------------------------------------------------------------------------
// FillHistograms
//------------------------------------------------------------------------------
void FillHistograms(UInt_t level, UInt_t check)
{
  if (check != 0) return;
  
  hWeff           [level]->Fill(1,            efficiencyW);
  hW              [level]->Fill(1,            totalW);
  hPtLepton1      [level]->Fill(pt1,          totalW);
  hPtLepton2      [level]->Fill(pt2,          totalW);
  hBDT1           [level]->Fill(bdt1,         totalW);
  hBDT2           [level]->Fill(bdt2,         totalW);
  hPtDiLepton     [level]->Fill(ptll,         totalW);
  hMinv           [level]->Fill(mll,          totalW);
  hMt             [level]->Fill(mth,          totalW);
  hNJets30        [level]->Fill(njet,         totalW);
  hpfMet          [level]->Fill(pfmet,        totalW);
  hppfMet         [level]->Fill(ppfmet,       totalW);
  hchMet          [level]->Fill(chmet,        totalW);
  hpchMet         [level]->Fill(pchmet,       totalW);
  hpminMet        [level]->Fill(mpmet,        totalW);
  hDeltaRLeptons  [level]->Fill(drll,         totalW);
  hDeltaPhiLeptons[level]->Fill(dphill,       totalW);
  hDPhiPtllJet    [level]->Fill(dphilljet,    totalW);
  hHwidthMVA      [level]->Fill(HwidthMVAbkg, totalW);
}
