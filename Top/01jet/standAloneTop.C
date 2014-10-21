#include "DY.C"
#include "Top.C"


TString directory = "rootfiles";


//------------------------------------------------------------------------------
//
// standAloneTop
//
//------------------------------------------------------------------------------
void standAloneTop(Int_t  njet          = 0,
		   Bool_t useDataDriven = true,
		   Int_t  printLevel    = 2)
{
  Double_t dyYield;
  Double_t dyStatError;
  Double_t dySystError;
  Double_t dyScaleFactor;

  DY(dyYield,
     dyStatError,
     dySystError,
     dyScaleFactor,
     njet,
     "SF",
     directory,
     useDataDriven,
     printLevel);

  Double_t topYield;
  Double_t topStatError;
  Double_t topSystError;
  Double_t topScaleFactor;

  std::cout<<"dyScaleFactor: "<<dyScaleFactor<<std::endl;
  Top(topYield,
      topStatError,
      topSystError,
      topScaleFactor,
      dyScaleFactor,
      njet,
      "All",
      directory,
      useDataDriven,
      printLevel);
}
