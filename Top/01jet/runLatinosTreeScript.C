void runLatinosTreeScript(Double_t luminosity    = 19.468,
			  Int_t    jetChannel    = 0,
			  TString  flavorChannel = "EMu",
			  TString  theSample     = "WJetsFakes_Total",
			  Bool_t   verbose       = true)
{
  gInterpreter->LoadMacro("LatinosTreeScript.C+");

  LatinosTreeScript(luminosity, jetChannel, flavorChannel, theSample, verbose);
}
