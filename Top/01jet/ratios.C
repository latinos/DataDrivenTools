//------------------------------------------------------------------------------
//
// Vector boson pair production at the LHC, arXiv:1105.0020
//
// All cross sections are in pb, unless stated otherwise.
//
//------------------------------------------------------------------------------
const Double_t xs7tev            = 47.04;
const Double_t xs7tev_upperError = 4.3 * xs7tev / 1e2;
const Double_t xs7tev_lowerError = 3.2 * xs7tev / 1e2;

const Double_t xs8tev            = 57.25;
const Double_t xs8tev_upperError = 4.1 * xs8tev / 1e2;
const Double_t xs8tev_lowerError = 2.8 * xs8tev / 1e2;


//------------------------------------------------------------------------------
//
// ratios
//
//------------------------------------------------------------------------------
const UInt_t nchannels = 3;


void ratios()
{
  gStyle->SetEndErrorSize(5);

  Double_t xsValue  [nchannels] = {   51.9,  52.4,  69.9};  // previous ATLAS 54.4
  Double_t statError[nchannels] = {    2.0,   2.0,   2.8};  // previous ATLAS  4.0
  Double_t systError[nchannels] = {    3.9,   4.5,   5.6};  // previous ATLAS  3.9
  Double_t lumiError[nchannels] = {    2.0,   1.2,   3.1};  // previous ATLAS  2.0
  Double_t lumiUsed [nchannels] = {    4.6,  4.92,  3.54};  // previous ATLAS  1.02
  Int_t    energy   [nchannels] = {      7,     7,     8};
  TString  label    [nchannels] = {"ATLAS", "CMS", "CMS"};


  // Loop
  //----------------------------------------------------------------------------
  TGraphErrors* gExp  = new TGraphErrors(nchannels);
  TGraphErrors* gTheo = new TGraphErrors(nchannels);

  for (UInt_t i=0; i<nchannels; i++) {

    Double_t nloValue = (i <= 1) ? xs7tev : xs8tev;
    Double_t nloError = (i <= 1) ? xs7tev_upperError : xs8tev_upperError;

    Double_t f = xsValue[i] / nloValue;

    Double_t expError = 0;
    expError += statError[i]*statError[i];
    expError += systError[i]*systError[i];
    expError += lumiError[i]*lumiError[i];

    expError = sqrt(expError);
    expError = f * expError / xsValue[i];

    Double_t theoError  = f * nloError / nloValue;
    Double_t totalError = sqrt(expError*expError + theoError*theoError);

    gExp ->SetPoint(i, f, i+1);
    gTheo->SetPoint(i, f, i+1);

    gExp ->SetPointError(i, expError,   0.0);
    gTheo->SetPointError(i, totalError, 0.0);
  }


  // Cosmetics
  //----------------------------------------------------------------------------
  gExp->SetLineWidth  (2);
  gExp->SetMarkerSize (1.3);
  gExp->SetMarkerStyle(kFullCircle);

  gTheo->SetLineColor  (kBlue);
  gTheo->SetLineWidth  (2);
  gTheo->SetMarkerSize (1.3);
  gTheo->SetMarkerStyle(kFullCircle);


  // Draw
  //----------------------------------------------------------------------------
  TCanvas* canvas = new TCanvas();

  canvas->SetLeftMargin(canvas->GetRightMargin());

  Double_t xmin = 0.70;
  Double_t xmax = 1.98;
  Double_t ymin = 0.50;
  Double_t ymax = 3.50;
  
  TH2F* dummy = new TH2F("dummy", "", 100, xmin, xmax, 100, ymin, ymax);

  dummy->Draw();


  // Vertical line at 1
  //----------------------------------------------------------------------------
  TLine* line = new TLine(1.0, ymin, 1.0, ymax);

  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);

  line->Draw("same");


  // Ratios
  //----------------------------------------------------------------------------
  gTheo->Draw("p||,same");
  gExp ->Draw("p,same");


  // Labels
  //----------------------------------------------------------------------------
  for (UInt_t i=0; i<nchannels; i++) {

    Double_t x = gExp->GetX()[i];
    Double_t y = gExp->GetY()[i];

    DrawTLatex(xmin+0.03, y+0.15, 0.035, 12, Form("%s",                label[i].Data()));
    DrawTLatex(xmin+0.03, y,      0.035, 12, Form("#sqrt{s} = %d TeV", energy[i]));
    DrawTLatex(xmin+0.03, y-0.15, 0.035, 12, Form("L = %.2f fb^{-1}",  lumiUsed[i]));

    Double_t expError   = gExp->GetErrorX(i);
    Double_t totalError = gTheo->GetErrorX(i);
    Double_t theoError  = sqrt(totalError*totalError - expError*expError);

    DrawTLatex(xmax-0.03, y, 0.035, 32, Form("%.2f #pm %.2f exp #pm %.2f theo",
					     x, expError, theoError));
  }

  dummy->GetXaxis()->CenterTitle();
  dummy->GetXaxis()->SetTitle("#sigma_{WW}^{exp} / #sigma_{WW}^{theo}");
  dummy->GetYaxis()->SetTitle("");


  // Remove y-axis labels
  //----------------------------------------------------------------------------
  TAxis* yaxis = dummy->GetYaxis();
  
  for (UInt_t j=1; j<yaxis->GetNbins(); j++) yaxis->SetBinLabel(j, "");


  // Save
  //----------------------------------------------------------------------------
  canvas->Update();
  canvas->GetFrame()->DrawClone();
  canvas->RedrawAxis();

  canvas->SaveAs("ratios.png");
}


//------------------------------------------------------------------------------
// DrawTLatex
//------------------------------------------------------------------------------
void DrawTLatex(Double_t    x,
		Double_t    y,
		Double_t    tsize,
		Short_t     align,
		const char* text)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetTextAlign(align);
  tl->SetTextFont (   42);
  tl->SetTextSize (tsize);

  tl->Draw("same");
}
