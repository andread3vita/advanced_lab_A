#include "./../include/GrUtil.h"
#include "./../include/StateFile.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TPaveStats.h"

#include "TStyle.h"
#include <fstream>
#include <string>

double calibration(int n, double m, double q)
{
  /*
  arguments explanation: n, numer that we need to convert in ns
                         m, slope of the calibration function
                         q, intercept of the calibration function

  This function calibrates the value provided by arietta in arbitraty units into nanoseconds
  */
  return q + m * n;
}

TH1D *FitData(std::string filename)
{
  /*
  arguments explanation: filename, the path to the file that contains all the data

  This function takes all the data from arietta, calibrates them in ns, than fill an instogram.
  From this istogram a fit of an exponential + constant (random backgound) is made.
  The results of the fit are then saved into the result.txt statefile, together with a plot in .pdf
  */

  StateFile *calibration_parameters = new StateFile("./../cal_arietta/results.txt");

  TCanvas *can = new TCanvas("c", "c", 900, 600);
  can->SetGrid();

  // loading data
  TH1D *h = new TH1D("h", "h", 75, 0, 8000);
  h->GetXaxis()->SetTitle("time [#mus]");
  h->GetYaxis()->SetTitle(GrUtil::HLabel(h, "#mus").c_str());
  h->SetTitle("Fit muon lifetime");
  std::ifstream input(filename.c_str());
  float         value;
  while (input >> value)
  {
    h->Fill(calibration(value, std::stoi(calibration_parameters->ValueOf("m")), std::stoi(calibration_parameters->ValueOf("q"))));
  }

  // fit
  TF1 *f = new TF1("f", "[0]+[1]*exp(-x/[2])", 300, 8000);
  f->SetParNames("constant", "normalization", "decay rate [ns]");
  f->SetParameter(2, 2200);
  f->SetParLimits(0, 0, 200);
  h->Fit(f, "R");

  // graphics
  gStyle->SetOptStat("rm");
  gStyle->SetOptFit(1);
  TPaveStats *pv = (TPaveStats *)h->GetListOfFunctions()->FindObject("stats");
  pv->SetX1NDC(0.45);
  pv->SetX2NDC(0.9);
  pv->SetY1NDC(0.55);
  pv->SetY2NDC(0.9);
  GrUtil::SetHTextSize(h);
  GrUtil::SetCountDigits(h);
  can->SaveAs("./lifetime_fit.pdf");

  // state file output
  StateFile *results = new StateFile("./results.txt");
  results->UpdateValueOf("life_time", std::to_string(f->GetParameter(2)));
  results->UpdateValueOf("life_time_error", std::to_string(f->GetParError(2)));

  return h;
}