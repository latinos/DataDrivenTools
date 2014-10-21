#include "TFile.h"
#include "TH1F.h"


const Double_t twScaleFactor = 1.0;   // A value of ~1.03 is computed in the 1-jet bin for HWW
const Double_t wjSyst        = 0.36;
const Double_t twSyst        = 0.17;  // 17% systematic for tt-tW cross-section difference


Double_t ratioValue(Double_t num,
		    Double_t den);

Double_t ratioError(Double_t num,
		    Double_t den,
		    Double_t numErr,
		    Double_t denErr);


//------------------------------------------------------------------------------
// Top
//------------------------------------------------------------------------------
void Top(Double_t &yield,
	 Double_t &statError,
	 Double_t &systError,
	 Double_t &topScaleFactor,
	 Double_t  dyScaleFactor,
	 Int_t     njet,
	 TString   channel,
	 TString   directory,
	 Bool_t    useDataDriven,
	 Int_t     printLevel)
{
  std::cout<<"Input parameters: "<<yield<<", "<<statError<<", "<<systError<<", "<<topScaleFactor<<", "<<dyScaleFactor<<", "<<njet<<", "<<channel<<", "<<directory<<", "<<useDataDriven<<", "<<printLevel<<std::endl;
  TString path = Form("%s/%djet/%s/", directory.Data(), njet, channel.Data());

  TFile* inputTT   = new TFile(path + "TTbar.root");  
  TFile* inputTW   = new TFile(path + "TW.root");  
  TFile* inputWW   = new TFile(path + "WW.root");  
  TFile* inputWZ   = new TFile(path + "WZ.root");  
  TFile* inputWG   = new TFile(path + "Wgamma.root");
  TFile* inputH125 = new TFile(path + "HWW125.root");
  TFile* inputZZ   = new TFile(path + "ZZ.root");  
  TFile* inputZjnt = new TFile(path + "DY.root");
  TFile* inputWj   = new TFile(path + "WJetsFakes_Total.root");
  TFile* inputData = new TFile(path + "DataRun2012_Total.root");  

  
  //----------------------------------------------------------------------------
  //
  // From LatinosTreeScript.C
  //
  // Float_t pMetCut        = 20 + 25*sameflav;
  // Bool_t  commonCutsTrue = (ptll > 45 && nextra == 0 && (dphilljet < 165.*TMath::Pi()/180. || !sameflav) && mll > 12);
  // Bool_t  topTagged_1jet = (commonCutsTrue && pfmet > 20 && zveto == 1 && mpmet > pMetCut && njet == 1 && nbjet == 1);
  // Bool_t  topTagged_0jet = (commonCutsTrue && pfmet > 20 && zveto == 1 && mpmet > pMetCut && njet  < 1 && !bveto);
  // Bool_t  analysisCuts   = (commonCutsTrue && pfmet > 20 && zveto == 1 && mpmet > pMetCut && njet == 0 && bveto_mu == 1 && bveto_ip == 1);
  // Bool_t  jetVetoCut     = (mll > 12       && pfmet > 20 && zveto == 1 && mpmet > pMetCut && njet == 0);
  //
  // From step3_cff.py
  //
  // nbjet is the number of jets with pt > 30 and jetBProbabilityBJetTags > 1.05
  //
  // bveto_nj30 is true if there are
  // 0 jets b-tagged (trackCountingHighEffBJetTags > 2.1) with 10 < pt <= 30 GeV
  // and
  // 0 soft muons with pt > 3 GeV and deltaR(muon, jets) >= 0.3, for jets with pt > 30 GeV
  //
  // bveto is true if there are 
  // 0 jets b-tagged (trackCountingHighEffBJetTags > 2.1) with 10 < pt < = 30 GeV
  // and
  // 0 soft muons with pt > 30 GeV
  //
  //----------------------------------------------------------------------------
  //
  //                               eff_2b                                  = 1b-tagging (TCHE or soft muon) efficiency in control (b-enriched) sample
  //                           1 - eff_2b                                  = 1b-veto efficiency in control (b-enriched) sample
  //                          (1 - eff_2b) * (1 - eff_2b)                  = 2b-veto efficiency in control (b-enriched) sample
  // eff_0b_tt  =         1 - (1 - eff_2b) * (1 - eff_2b)                  = efficiency for a      tt         event to pass the veto efficiency
  // eff_0b_top = f_tt * [1 - (1 - eff_2b) * (1 - eff_2b)] + f_tW * eff_2b = efficiency for a top (tt and tW) event to pass the veto efficiency
  //
  // N_top_0jet = (N_top_tagged_data - N_top_tagged_mc_bkg) * (1 - eff_0b_top) / eff_0b_top
  // 
  //----------------------------------------------------------------------------


  TH1F* hNTopTaggedEventsControlRegionTT = (TH1F*) inputTT->Get("hNTopTaggedTopControlRegion");  // if (topTagged_1jet && bveto_nj30 == 0)
  TH1F* hNTopControlRegionTT             = (TH1F*) inputTT->Get("hNTopControlRegion");           // if (topTagged_1jet)
  TH1F* hTopTaggedEventsTT               = (TH1F*) inputTT->Get("hTopTaggedEvents");             // if (topTagged_0jet)
  TH1F* hTopTaggingTT                    = (TH1F*) inputTT->Get("hWTopTagging");                 // if (analysisCuts)
  TH1F* hJetVetoTT                       = (TH1F*) inputTT->Get("hWJetVeto");                    // if (jetVetoCut)

  TH1F* hNTopTaggedEventsControlRegionTW = (TH1F*) inputTW->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionTW             = (TH1F*) inputTW->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsTW               = (TH1F*) inputTW->Get("hTopTaggedEvents");
  TH1F* hTopTaggingTW                    = (TH1F*) inputTW->Get("hWTopTagging");
  TH1F* hJetVetoTW                       = (TH1F*) inputTW->Get("hWJetVeto");

  TH1F* hNTopTaggedEventsControlRegionWW = (TH1F*) inputWW->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionWW             = (TH1F*) inputWW->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsWW               = (TH1F*) inputWW->Get("hTopTaggedEvents");

  TH1F* hNTopTaggedEventsControlRegionWZ = (TH1F*) inputWZ->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionWZ             = (TH1F*) inputWZ->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsWZ               = (TH1F*) inputWZ->Get("hTopTaggedEvents");

  TH1F* hNTopTaggedEventsControlRegionWG = (TH1F*) inputWG->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionWG             = (TH1F*) inputWG->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsWG               = (TH1F*) inputWG->Get("hTopTaggedEvents");

  TH1F* hNTopTaggedEventsControlRegionH125 = (TH1F*) inputH125->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionH125             = (TH1F*) inputH125->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsH125               = (TH1F*) inputH125->Get("hTopTaggedEvents");

  TH1F* hNTopTaggedEventsControlRegionZZ = (TH1F*) inputZZ->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionZZ             = (TH1F*) inputZZ->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsZZ               = (TH1F*) inputZZ->Get("hTopTaggedEvents");

  TH1F* hNTopTaggedEventsControlRegionWj = (TH1F*) inputWj->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionWj             = (TH1F*) inputWj->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsWj               = (TH1F*) inputWj->Get("hTopTaggedEvents");

  TH1F* hNTopTaggedEventsControlRegionZjnt = (TH1F*) inputZjnt->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionZjnt		   = (TH1F*) inputZjnt->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsZjnt               = (TH1F*) inputZjnt->Get("hTopTaggedEvents");

  TH1F* hNTopTaggedEventsControlRegionDATA = (TH1F*) inputData->Get("hNTopTaggedTopControlRegion");
  TH1F* hNTopControlRegionDATA             = (TH1F*) inputData->Get("hNTopControlRegion");
  TH1F* hTopTaggedEventsData               = (TH1F*) inputData->Get("hTopTaggedEvents");


  // Yields in (topTagged_1jet && bveto_nj30 == 0)
  //----------------------------------------------------------------------------
  Double_t NnumTT   = hNTopTaggedEventsControlRegionTT  ->GetBinContent(2);
  Double_t NnumTW   = hNTopTaggedEventsControlRegionTW  ->GetBinContent(2) * twScaleFactor;
  Double_t NnumWW   = hNTopTaggedEventsControlRegionWW  ->GetBinContent(2);
  Double_t NnumWZ   = hNTopTaggedEventsControlRegionWZ  ->GetBinContent(2);
  Double_t NnumWG   = hNTopTaggedEventsControlRegionWG  ->GetBinContent(2);
  Double_t NnumH125 = hNTopTaggedEventsControlRegionH125->GetBinContent(2);
  Double_t NnumZZ   = hNTopTaggedEventsControlRegionZZ  ->GetBinContent(2);
  Double_t NnumWj   = hNTopTaggedEventsControlRegionWj  ->GetBinContent(2);
  Double_t NnumZjnt = hNTopTaggedEventsControlRegionZjnt->GetBinContent(2) * dyScaleFactor;
  Double_t NnumData = hNTopTaggedEventsControlRegionDATA->GetBinContent(2);
  
  Double_t NnumBkg       = NnumTW + NnumWW + NnumWZ + NnumWG + NnumH125 + NnumZZ + NnumWj + NnumZjnt;
  Double_t NnumDataFinal = NnumData - NnumBkg;

  Double_t errNnumTT   = hNTopTaggedEventsControlRegionTT  ->GetBinError(2);
  Double_t errNnumTW   = hNTopTaggedEventsControlRegionTW  ->GetBinError(2) * twScaleFactor;
  Double_t errNnumWW   = hNTopTaggedEventsControlRegionWW  ->GetBinError(2);
  Double_t errNnumWZ   = hNTopTaggedEventsControlRegionWZ  ->GetBinError(2);
  Double_t errNnumWG   = hNTopTaggedEventsControlRegionWG  ->GetBinError(2);
  Double_t errNnumH125 = hNTopTaggedEventsControlRegionH125->GetBinError(2);
  Double_t errNnumZZ   = hNTopTaggedEventsControlRegionZZ  ->GetBinError(2);
  Double_t errNnumWj   = hNTopTaggedEventsControlRegionWj  ->GetBinError(2);
  Double_t errNnumZjnt = hNTopTaggedEventsControlRegionZjnt->GetBinError(2) * dyScaleFactor;
  Double_t errNnumData = hNTopTaggedEventsControlRegionDATA->GetBinError(2);

  Double_t errNnumBkg = 0;
  errNnumBkg += errNnumTW   * errNnumTW;
  errNnumBkg += errNnumWW   * errNnumWW;
  errNnumBkg += errNnumWZ   * errNnumWZ;
  errNnumBkg += errNnumWG   * errNnumWG;
  errNnumBkg += errNnumH125 * errNnumH125;
  errNnumBkg += errNnumZZ   * errNnumZZ;
  errNnumBkg += errNnumWj   * errNnumWj;
  errNnumBkg += errNnumZjnt * errNnumZjnt;
  errNnumBkg  = sqrt(errNnumBkg);

  Double_t errNnumDataFinal = sqrt(errNnumData*errNnumData + errNnumBkg*errNnumBkg);


  // Yields in (topTagged_1jet)
  //----------------------------------------------------------------------------
  Double_t NdenTT   = hNTopControlRegionTT  ->GetBinContent(2);
  Double_t NdenTW   = hNTopControlRegionTW  ->GetBinContent(2) * twScaleFactor;
  Double_t NdenWW   = hNTopControlRegionWW  ->GetBinContent(2);
  Double_t NdenWZ   = hNTopControlRegionWZ  ->GetBinContent(2);
  Double_t NdenWG   = hNTopControlRegionWG  ->GetBinContent(2);
  Double_t NdenH125 = hNTopControlRegionH125->GetBinContent(2);
  Double_t NdenZZ   = hNTopControlRegionZZ  ->GetBinContent(2);
  Double_t NdenWj   = hNTopControlRegionWj  ->GetBinContent(2);
  Double_t NdenZjnt = hNTopControlRegionZjnt->GetBinContent(2) * dyScaleFactor;
  Double_t NdenData = hNTopControlRegionDATA->GetBinContent(2);

  Double_t NdenBkg       = NdenTW + NdenWW + NdenWZ + NdenWG + NdenH125 + NdenZZ + NdenWj + NdenZjnt;
  Double_t NdenDataFinal = NdenData - NdenBkg;

  Double_t errNdenTT   = hNTopControlRegionTT  ->GetBinError(2);
  Double_t errNdenTW   = hNTopControlRegionTW  ->GetBinError(2) * twScaleFactor;
  Double_t errNdenWW   = hNTopControlRegionWW  ->GetBinError(2);
  Double_t errNdenWZ   = hNTopControlRegionWZ  ->GetBinError(2);
  Double_t errNdenWG   = hNTopControlRegionWG  ->GetBinError(2);
  Double_t errNdenH125 = hNTopControlRegionH125->GetBinError(2);
  Double_t errNdenZZ   = hNTopControlRegionZZ  ->GetBinError(2);
  Double_t errNdenWj   = hNTopControlRegionWj  ->GetBinError(2);
  Double_t errNdenZjnt = hNTopControlRegionZjnt->GetBinError(2) * dyScaleFactor;
  Double_t errNdenData = hNTopControlRegionDATA->GetBinError(2);

  Double_t errNdenBkg = 0;
  errNdenBkg += errNdenTW   * errNdenTW;
  errNdenBkg += errNdenWW   * errNdenWW;
  errNdenBkg += errNdenWZ   * errNdenWZ;
  errNdenBkg += errNdenWG   * errNdenWG;
  errNdenBkg += errNdenH125 * errNdenH125;
  errNdenBkg += errNdenZZ   * errNdenZZ;
  errNdenBkg += errNdenWj   * errNdenWj;
  errNdenBkg += errNdenZjnt * errNdenZjnt;
  errNdenBkg  = sqrt(errNdenBkg);

  Double_t errNdenDataFinal = sqrt(errNdenData*errNdenData + errNdenBkg*errNdenBkg);


  // tt and tW yields at jet veto (njet == 0) level
  //----------------------------------------------------------------------------
  Double_t NJetVetott         = hJetVetoTT->GetBinContent(2);
  Double_t NJetVetotW         = hJetVetoTW->GetBinContent(2);
  Double_t NJetVetoTopTotalMC = NJetVetott + NJetVetotW;
  Double_t fTT                = NJetVetott / NJetVetoTopTotalMC;


  // Compute efficiencies
  //----------------------------------------------------------------------------
  Double_t eff2bTT                 = NnumTT / NdenTT;
  Double_t eff2bTW                 = NnumTW / NdenTW;
  Double_t eff2bDataFinal          = NnumDataFinal / NdenDataFinal;
  Double_t eff0b                   = 1 - (1 - eff2bTT)        * (1 - eff2bTT);
  Double_t eff0bDataFinal          = 1 - (1 - eff2bDataFinal) * (1 - eff2bDataFinal);
  Double_t eff0b_fraction          = fTT * eff0b          + (1 - fTT) * eff2bTT;
  Double_t eff0bDataFinal_fraction = fTT * eff0bDataFinal + (1 - fTT) * eff2bDataFinal;

  Double_t erreff2bTT                 = ratioError(NnumTT, NdenTT, errNnumTT, errNdenTT);
  Double_t erreff2bTW                 = ratioError(NnumTW, NdenTW, errNnumTW, errNdenTW);
  Double_t erreff2bDataFinal          = ratioError(NnumDataFinal, NdenDataFinal, errNnumDataFinal, errNdenDataFinal);
  Double_t erreff0b                   = 2 * (1 - eff2bTT)        * erreff2bTT;
  Double_t erreff0bDataFinal          = 2 * (1 - eff2bDataFinal) * erreff2bDataFinal;
  Double_t erreff0b_fraction          = erreff0b;
  Double_t erreff0bDataFinal_fraction = erreff0bDataFinal;


  if (printLevel > 1) {
    printf("\n");
    printf(" N^{control}_{2b} (numerators)\n ------------------------------------\n");
    printf("       tt                = %9.4f +- %8.4f\n", NnumTT,        errNnumTT);
    printf("       tW                = %9.4f +- %8.4f\n", NnumTW,        errNnumTW);
    printf("       WW                = %9.4f +- %8.4f\n", NnumWW,        errNnumWW);
    printf("       WZ                = %9.4f +- %8.4f\n", NnumWZ,        errNnumWZ);
    printf("       Wg                = %9.4f +- %8.4f\n", NnumWG,        errNnumWG);
    printf("       H125              = %9.4f +- %8.4f\n", NnumH125,      errNnumH125);
    printf("       ZZ                = %9.4f +- %8.4f\n", NnumZZ,        errNnumZZ);
    printf("       Wj                = %9.4f +- %8.4f\n", NnumWj,        errNnumWj);
    printf("       Zjnt              = %9.4f +- %8.4f\n", NnumZjnt,      errNnumZjnt);
    printf("       bkg (includes tW) = %9.4f +- %8.4f\n", NnumBkg,       errNnumBkg);
    printf("       data              = %9.0f +- %8.0f\n", NnumData,      errNnumData);
    printf("       data - bkg        = %9.4f +- %8.4f\n", NnumDataFinal, errNnumDataFinal);

    printf("\n\n");
    printf(" N^{control} (denominators)\n ------------------------------------\n");
    printf("       tt                = %9.4f +- %8.4f\n", NdenTT,        errNdenTT);
    printf("       tW                = %9.4f +- %8.4f\n", NdenTW,        errNdenTW);
    printf("       WW                = %9.4f +- %8.4f\n", NdenWW,        errNdenWW);
    printf("       WZ                = %9.4f +- %8.4f\n", NdenWZ,        errNdenWZ);
    printf("       Wg                = %9.4f +- %8.4f\n", NdenWG,        errNdenWG);
    printf("       H125              = %9.4f +- %8.4f\n", NdenH125,      errNdenH125);
    printf("       ZZ                = %9.4f +- %8.4f\n", NdenZZ,        errNdenZZ);
    printf("       Wj                = %9.4f +- %8.4f\n", NdenWj,        errNdenWj);
    printf("       Zjnt              = %9.4f +- %8.4f\n", NdenZjnt,      errNdenZjnt);
    printf("       bkg (includes tW) = %9.4f +- %8.4f\n", NdenBkg,       errNdenBkg);
    printf("       data              = %9.0f +- %8.0f\n", NdenData,      errNdenData);
    printf("       data - bkg        = %9.4f +- %8.4f\n", NdenDataFinal, errNdenDataFinal);

    printf("\n\n");
    printf(" efficiencies\n ------------------------------------\n");
    printf("       eff^{soft}_{2b}_tW   = (%5.2f +- %4.2f)%s\n", 1e2*eff2bTW,                 1e2*erreff2bTW,                 "%");
    printf("       eff^{soft}_{2b}_tt   = (%5.2f +- %4.2f)%s\n", 1e2*eff2bTT,                 1e2*erreff2bTT,                 "%");
    printf("       eff^{soft}_{2b}_data = (%5.2f +- %4.2f)%s\n", 1e2*eff2bDataFinal,          1e2*erreff2bDataFinal,          "%");
    printf("       eff_{0b}_tt          = (%5.2f +- %4.2f)%s\n", 1e2*eff0b,                   1e2*erreff0b,                   "%");
    printf("       eff_{0b}_data        = (%5.2f +- %4.2f)%s\n", 1e2*eff0bDataFinal,          1e2*erreff0bDataFinal,          "%");
    printf(" >>>>> eff_{0b}_top,MC      = (%5.2f +- %4.2f)%s\n", 1e2*eff0b_fraction,          1e2*erreff0b_fraction,          "%");
    printf(" >>>>> eff_{0b}_top,data    = (%5.2f +- %4.2f)%s\n", 1e2*eff0bDataFinal_fraction, 1e2*erreff0bDataFinal_fraction, "%");
  }


  // tt and tW yields at analysis level
  //----------------------------------------------------------------------------
  Double_t Ntt         = hTopTaggingTT->GetBinContent(2);
  Double_t NtW         = hTopTaggingTW->GetBinContent(2);
  Double_t NTopTotalMC = Ntt + NtW;

  Double_t errNtt         = hTopTaggingTT->GetBinError(2);
  Double_t errNtW         = hTopTaggingTW->GetBinError(2);
  Double_t errNTopTotalMC = sqrt(errNtt*errNtt + errNtW*errNtW);


  // Yields in (topTagged_0jet)
  //----------------------------------------------------------------------------
  Double_t NTaggedTT    = hTopTaggedEventsTT  ->GetBinContent(2);
  Double_t NTaggedTW    = hTopTaggedEventsTW  ->GetBinContent(2);
  Double_t NTaggedWW    = hTopTaggedEventsWW  ->GetBinContent(2);
  Double_t NTaggedWZ    = hTopTaggedEventsWZ  ->GetBinContent(2);
  Double_t NTaggedWG    = hTopTaggedEventsWG  ->GetBinContent(2);
  Double_t NTaggedH125  = hTopTaggedEventsH125->GetBinContent(2);
  Double_t NTaggedZZ    = hTopTaggedEventsZZ  ->GetBinContent(2);
  Double_t NTaggedWj    = hTopTaggedEventsWj  ->GetBinContent(2);
  Double_t NTaggedZjnt  = hTopTaggedEventsZjnt->GetBinContent(2) * dyScaleFactor;
  Double_t NTaggedDATA  = hTopTaggedEventsData->GetBinContent(2);

  Double_t errNTaggedTT   = hTopTaggedEventsTT  ->GetBinError(2);
  Double_t errNTaggedTW   = hTopTaggedEventsTW  ->GetBinError(2);
  Double_t errNTaggedWW   = hTopTaggedEventsWW  ->GetBinError(2);
  Double_t errNTaggedWZ   = hTopTaggedEventsWZ  ->GetBinError(2);
  Double_t errNTaggedWG   = hTopTaggedEventsWG  ->GetBinError(2);
  Double_t errNTaggedH125 = hTopTaggedEventsH125->GetBinError(2);
  Double_t errNTaggedZZ   = hTopTaggedEventsZZ  ->GetBinError(2);
  Double_t errNTaggedWj   = hTopTaggedEventsWj  ->GetBinError(2);
  Double_t errNTaggedZjnt = hTopTaggedEventsZjnt->GetBinError(2) * dyScaleFactor;
  Double_t errNTaggedDATA = hTopTaggedEventsData->GetBinError(2);
				   
  Double_t Nbkg = NTaggedWW + NTaggedWZ + NTaggedWG + NTaggedH125 + NTaggedZZ + NTaggedWj + NTaggedZjnt;

  Double_t errNbkg = 0;
  errNbkg += errNTaggedWW     * errNTaggedWW;
  errNbkg += errNTaggedWZ     * errNTaggedWZ;
  errNbkg += errNTaggedWG     * errNTaggedWG;
  errNbkg += errNTaggedH125   * errNTaggedH125;
  errNbkg += errNTaggedZZ     * errNTaggedZZ;
  errNbkg += errNTaggedWj     * errNTaggedWj;
  errNbkg += wjSyst*NTaggedWj * wjSyst*NTaggedWj;
  errNbkg += errNTaggedZjnt   * errNTaggedZjnt;
  errNbkg = sqrt(errNbkg);

  
  // Estimate the top yield in the 0-jet bin
  //----------------------------------------------------------------------------
  Double_t NTopEstimatedDataFinal_fraction = (NTaggedDATA - Nbkg) * (1 - eff0bDataFinal_fraction) / eff0bDataFinal_fraction;

  Double_t err1DataFinal_f = errNTaggedDATA * (1 - eff0bDataFinal_fraction) / eff0bDataFinal_fraction;
  Double_t err2DataFinal_f = errNbkg        * (1 - eff0bDataFinal_fraction) / eff0bDataFinal_fraction; 
  Double_t err3DataFinal_f = (NTaggedDATA - Nbkg) * erreff0bDataFinal_fraction / eff0bDataFinal_fraction;
  Double_t err4DataFinal_f = NTopEstimatedDataFinal_fraction * erreff0bDataFinal_fraction / eff0bDataFinal_fraction;

  Double_t errNTopEstimatedDataFinal_fraction = 0;
  errNTopEstimatedDataFinal_fraction += err1DataFinal_f * err1DataFinal_f;
  errNTopEstimatedDataFinal_fraction += err2DataFinal_f * err2DataFinal_f;
  errNTopEstimatedDataFinal_fraction += err3DataFinal_f * err3DataFinal_f;
  errNTopEstimatedDataFinal_fraction += err4DataFinal_f * err4DataFinal_f;
  errNTopEstimatedDataFinal_fraction = sqrt(errNTopEstimatedDataFinal_fraction);

  Double_t SFDataFinal_f = NTopEstimatedDataFinal_fraction / NTopTotalMC;


  // Closure test on MC
  //----------------------------------------------------------------------------
  Double_t NTaggedTopMC    = NTaggedTT + NTaggedTW;
  Double_t errNTaggedTopMC = sqrt(errNTaggedTT*errNTaggedTT + errNTaggedTW*errNTaggedTW);

  Double_t NTopEstimatedClosureMC = NTaggedTopMC * (1 - eff0b_fraction) / eff0b_fraction;

  Double_t err1MC = errNTaggedTopMC * (1 - eff0b_fraction) / eff0b_fraction;
  Double_t err3MC = NTaggedTopMC * erreff0b / eff0b_fraction; 
  Double_t err4MC = NTaggedTopMC * (1 - eff0b_fraction) * erreff0b / (eff0b_fraction * eff0b_fraction);

  Double_t errNTopEstimatedClosureMC = sqrt(err1MC*err1MC + err3MC*err3MC + err4MC*err4MC);


  if (printLevel > 1) {
    printf("\n\n");
    printf(" yields in (topTagged_0jet)\n ------------------------------------\n");
    printf("       tt   = %9.4f +- %8.4f\n", NTaggedTT,   errNTaggedTW);
    printf("       tW   = %9.4f +- %8.4f\n", NTaggedTW,   errNTaggedTW);
    printf("       WW   = %9.4f +- %8.4f\n", NTaggedWW,   errNTaggedWW);
    printf("       WZ   = %9.4f +- %8.4f\n", NTaggedWZ,   errNTaggedWZ);
    printf("       Wg   = %9.4f +- %8.4f\n", NTaggedWG,   errNTaggedWG);
    printf("       H125 = %9.4f +- %8.4f\n", NTaggedH125, errNTaggedH125);
    printf("       ZZ   = %9.4f +- %8.4f\n", NTaggedZZ,   errNTaggedZZ);
    printf("       Wj   = %9.4f +- %8.4f\n", NTaggedWj,   errNTaggedWj);
    printf("       Zjnt = %9.4f +- %8.4f\n", NTaggedZjnt, errNTaggedZjnt);
    printf("       bkg  = %9.4f +- %8.4f\n", Nbkg,        errNbkg);
    printf("       data = %9.0f +- %8.0f\n", NTaggedDATA, errNTaggedDATA);

    printf("\n\n");
    printf(" top yield estimation\n ------------------------------------\n");
    printf("       N_top_data                 = %9.4f +- %8.4f\n", NTopEstimatedDataFinal_fraction, errNTopEstimatedDataFinal_fraction);
    printf("       N_top_mc^{cut based}       = %9.4f +- %8.4f\n", NTopTotalMC, errNTopTotalMC);
    printf("       N_top_mc^{data driven}     = %9.4f +- %8.4f\n", NTopEstimatedClosureMC, errNTopEstimatedClosureMC);
    printf(" >>>>> SF (data / mc^{cut based}) = %9.4f\n",          SFDataFinal_f);
  }


  // Final numbers
  //----------------------------------------------------------------------------
  Double_t xtw           = eff2bTW;
  Double_t finalfTT      = fTT + (1-fTT)*xtw;
  Double_t errorfinalfTT = twSyst*fTT*(1 - xtw);


  if (printLevel > 1) {
    printf("\n\n");
    printf(" tt and tW numbers before applying top veto cuts (jet veto level)\n");
    printf(" ------------------------------------\n");
    printf("          x =  %5.2f\n",   xtw);
    printf("       N_tt =  %5.2f\n",   NJetVetott);
    printf("       N_tW =  %5.2f\n",   NJetVetotW);
    printf("       f_tt =  %5.2f%s\n", 1e2*fTT,       "%");
    printf(" final f_tt = (%5.2f +- %4.2f)%s\n", 1e2*finalfTT, 1e2*errorfinalfTT, "%");

    printf("\n\n");
    printf(" final values\n ------------------------------------\n");
    printf("       eff_{0b}_top,data = (%5.2f +- %4.2f (stat.))%s", 1e2*eff0bDataFinal_fraction, 1e2*erreff0bDataFinal_fraction, "%");
    printf("                 (without tW correction)\n");
  }


  eff0bDataFinal_fraction = finalfTT*eff0bDataFinal + (1-finalfTT)*eff2bDataFinal;

  Double_t stat_error = erreff2bDataFinal * (1 + finalfTT*(1 - 2*eff2bDataFinal));
  Double_t syst_error = errorfinalfTT * (eff0bDataFinal - eff2bDataFinal);

  erreff0bDataFinal_fraction = sqrt(stat_error*stat_error + syst_error*syst_error);

  
  if (printLevel > 1) {
    printf("       eff_{0b}_top,data = (%5.2f +- %4.2f (stat.) +- %4.2f (syst.))%s",
	   1e2*eff0bDataFinal_fraction, 1e2*stat_error, 1e2*syst_error, "%");
    printf(" (with    tW correction)\n");
    printf("       eff_{0b}_top,data = (%5.2f +- %4.2f (total))%s\n",
	   1e2*eff0bDataFinal_fraction, 1e2*erreff0bDataFinal_fraction, "%");
  }


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  if (njet == 1) {

    eff0bDataFinal_fraction    = eff2bDataFinal;
    erreff0bDataFinal_fraction = erreff2bDataFinal;
  }
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  // Final top yield estimation
  //----------------------------------------------------------------------------
  NTopEstimatedDataFinal_fraction = (NTaggedDATA - Nbkg) * (1 - eff0bDataFinal_fraction) / eff0bDataFinal_fraction;

  err1DataFinal_f = errNTaggedDATA * (1 - eff0bDataFinal_fraction) / eff0bDataFinal_fraction;
  err2DataFinal_f =        errNbkg * (1 - eff0bDataFinal_fraction) / eff0bDataFinal_fraction;

  errNTopEstimatedDataFinal_fraction = sqrt(err1DataFinal_f*err1DataFinal_f + err2DataFinal_f*err2DataFinal_f);  // Final statistical error

  Double_t errtotal_syst = (NTaggedDATA - Nbkg) * erreff0bDataFinal_fraction / (eff0bDataFinal_fraction*eff0bDataFinal_fraction);  // Final systematic error


  // Save the result
  //----------------------------------------------------------------------------
  yield          = (useDataDriven) ? NTopEstimatedDataFinal_fraction : NTopTotalMC;
  statError      = errNTopEstimatedDataFinal_fraction;
  systError      = errtotal_syst;
  topScaleFactor = yield / NTopTotalMC;


  Double_t topScaleFactorError = ratioError(yield,
					    NTopTotalMC,
					    sqrt(statError*statError + systError*systError),
					    errNTopTotalMC);


  // Use this number for the note 14-192
  //----------------------------------------------------------------------------
  if (printLevel > 0) {
    printf("\n [%s] Top values for the note (AN-14-192)\n", channel.Data());
    printf(" -------------------------------------------------\n");
    printf(" N^{control}(bkg. sub.)               &   %4.0f\n",                   NdenDataFinal);
    printf(" N^{control}_{tag}(bkg. sub.)         &   %4.0f\n",                   NnumDataFinal);
    if (njet == 0) {
      printf(" eff_{tag}              (%s)           & %6.2f $\\pm$ %4.2f\n", " ", eff2bDataFinal,          erreff2bDataFinal);
      printf(" f_{tt}                 (%s)           &  %6.0f $\\pm$ %4.0f\n", "%", 1e2*finalfTT,            1e2*errorfinalfTT);
      printf(" eff_{2b}               (%s)           & %6.2f $\\pm$ %4.2f\n", " ", eff0bDataFinal,          erreff0bDataFinal);
    }
    printf(" eff_{top-tag}^{0-jet} (%s)           & %6.2f $\\pm$ %4.2f\n", "%", eff0bDataFinal_fraction, erreff0bDataFinal_fraction);
    printf(" top-tagged events in data            &   %4.0f\n",                   NTaggedDATA);
    printf(" background events in control region  &   %6.2f $\\pm$ %4.2f\n",      Nbkg,                    errNbkg);
    printf(" estimated top events in simulation   &   %6.2f $\\pm$ %4.2f\n",      NTopTotalMC,             errNTopTotalMC);
    printf(" data-driven top background estimate  &   %6.2f $\\pm$ %4.2f\n",      yield,                   statError);
    printf(" data/MC                              & %6.2f $\\pm$ %4.2f\n\n",    topScaleFactor,          topScaleFactorError);
    printf(" -------------------------------------------------\n");
  }

  // For the note
  //----------------------------------------------------------------------------
  if (printLevel > 0) {
    printf("\n [%s] Top values for the note\n", channel.Data());
    printf(" -------------------------------------------------\n");
    printf(" N^{control,data}                    &   %4.0f\n",                      NdenData);
    printf(" N^{control,data}_{tag}              &   %4.0f\n",                      NnumData);
    printf(" eff_{top-tag}^{data}  (%s)           & %6.2f $\\pm$ %4.2f\n", "%", 1e2*eff2bDataFinal,          1e2*erreff2bDataFinal);
    printf(" eff_{top-tag}^{tW}    (%s)           & %6.2f $\\pm$ %4.2f\n", "%", 1e2*eff2bTW,                 1e2*erreff2bTW);
    printf(" f_{tt}                (%s)           & %6.2f $\\pm$ %4.2f\n", "%", 1e2*finalfTT,                1e2*errorfinalfTT);
    printf(" eff_{top-tag}^{0-jet} (%s)           & %6.2f $\\pm$ %4.2f\n", "%", 1e2*eff0bDataFinal_fraction, 1e2*erreff0bDataFinal_fraction);
    printf(" data top-tagged events              &   %4.0f\n",                      NTaggedDATA);
    printf(" background events in control region & %6.2f $\\pm$ %4.2f\n",           Nbkg,                        errNbkg);
    printf(" estimated top events in simulation  & %6.2f $\\pm$ %4.2f\n",           NTopTotalMC,                 errNTopTotalMC);
    printf(" data-driven top estimate            & %6.2f $\\pm$ %4.2f\n",           yield,                       statError);
    printf(" data/MC                             & %6.2f $\\pm$ %4.2f\n\n",         topScaleFactor,              topScaleFactorError);
    printf("\n [%s] Top relative systematic uncertainties\n", channel.Data());
    printf(" -------------------------------------------------\n");
    printf(" top normalisation = %.0f (stat.) $\\bigoplus$ %.0f (syst.)\n\n",
	   1e2*statError/yield, 1e2*systError/yield);
    printf(" ---------------DataDriven numbers----------------------------------\n");
    printf(" NControl (data top-tagged events)            &   %4.0f\n",                     NTaggedDATA);
    printf(" FinalEstimation (data-driven top estimate)   & %6.1f $\\pm$ %4.1f\n",          yield,               statError);
    printf(" alpha (FinalEstimation/NControl)             & %6.19f $\\pm$ %4.19f\n",        yield/NTaggedDATA,   sqrt(statError*statError/(yield*yield)-1./NTaggedDATA)*yield/NTaggedDATA);
  }
}


//------------------------------------------------------------------------------
// ratioValue
//------------------------------------------------------------------------------
Double_t ratioValue(Double_t num,
		    Double_t den)
{
  if (den == 0) return -999;

  Double_t rValue = num / den;
  
  return rValue;
}


//------------------------------------------------------------------------------
// ratioError
//------------------------------------------------------------------------------
Double_t ratioError(Double_t num,
		    Double_t den,
		    Double_t numErr,
		    Double_t denErr)
{
  if (num == 0) return -999;

  Double_t rValue = ratioValue(num, den);

  Double_t rError = (numErr/num)*(numErr/num) + (denErr/den)*(denErr/den);

  rError = rValue * sqrt(rError);

  return rError;
}
