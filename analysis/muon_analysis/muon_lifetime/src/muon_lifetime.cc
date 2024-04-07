#include "./../../../../include/GrUtil.h"
#include "./../../../../include/StateFile.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TPaveStats.h"

#include "TStyle.h"
#include <fstream>
#include <string>

TH1D *h; //  histogram of events
TF1  *f; // funtion to fit

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

Double_t chi_sqr(Double_t *val, Double_t *par)
{
  /*
  This function is a user defined function for a TF2. Basically it computes manually the value of the chi^2
  of the fit of an histogram. The function fitted is an exponential in the form f(t)=c+x*exp(-t/y) where:
  - x is the normalization constant
  - t is the time
  - y the mean life time
  - c a constant value as background

  this function given the histogram h compute the chi^2 as sum over all fitted bins of ( binvalue(t)-f(t) )^2 / f(t).
  Where f(t) plays the role of the expected value and binvalue(t) is the hight of the bin with center t measured during
  the experiment.
  The function is a 2Dim function of (x,y) = (normalization,life time).
  */
  //  Check if the histogram to compute the chi² is not empty, otherwise exit
  if (h == nullptr)
  {
    std::cerr << "ERROR --> Trying to compute the chi² over an empy histogram, abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  Float_t x = val[0];
  Float_t y = val[1];

  Double_t sum = 0;

  //  loop over all the bin in   // double xMin = 2000;
  // double xMax = 16000;the histogram, but only in the range of the fit
  // double xMin = 2000;
  // double xMax = 16000;
  double xMin, xMax;
  f->GetRange(xMin, xMax);
  int bin_min = h->FindBin(xMin);
  int bin_max = h->FindBin(xMax) - 1;
  for (int i = bin_min; i <= bin_max; i++)
  {
    double expected = par[0] + x * exp(-1 * h->GetBinCenter(i) / y);
    sum += pow(h->GetBinContent(i) - expected, 2) / expected;
  }
  return sum;
}

TF1 *FitData(std::string filename, int xmin_fit = 2000, int xmax_fit = 16000, int total_bins = 75)
{
  /*
  arguments explanation: filename, the path to the file that contains all the data

  This function takes all the data from arietta, calibrates them in ns, than fill an instogram.
  From this istogram a fit of an exponential + constant (random backgound) is made.
  The results of the fit are then saved into the result.txt statefile, together with a plot in .pdf
  */

  int range_xmin = 0;
  int range_xmax = 16000;

  StateFile *calibration_parameters = new StateFile("./../../detector/arietta/results.txt");

  TCanvas *can = new TCanvas("c", "c", 900, 600);
  can->SetGrid();

  // loading data
  h = new TH1D("h", "h", total_bins, range_xmin, range_xmax);
  h->GetXaxis()->SetTitle("time [ns]");
  h->GetYaxis()->SetTitle(GrUtil::HLabel(h, "ns").c_str());
  h->SetTitle("Fit muon lifetime");
  std::ifstream input(filename.c_str());
  float         value;
  while (input >> value)
  {
    h->Fill(calibration(value, std::stod(calibration_parameters->ValueOf("m")), std::stod(calibration_parameters->ValueOf("q"))));
  }

  // fit
  f = new TF1("f", "[0]+[1]*exp(-x/[2])", xmin_fit, xmax_fit);
  f->SetParNames("constant", "normalization", "lifetime [ns]");
  f->SetParameter(2, 2200);
  f->SetParLimits(0, 0, 200);
  h->Fit(f, "RMQ");
  h->Draw("E1");

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
  StateFile *results = new StateFile("./results/results_single_fit.txt");

  results->UpdateValueOf("filename", filename);
  results->UpdateValueOf("xmin_fit", std::to_string(xmin_fit));
  results->UpdateValueOf("xmax_fit", std::to_string(xmax_fit));
  results->UpdateValueOf("total_bins", std::to_string(total_bins));
  results->UpdateValueOf("range_xmin", std::to_string(range_xmin));
  results->UpdateValueOf("range_xmax", std::to_string(range_xmax));

  results->UpdateValueOf("life_time", std::to_string(f->GetParameter(2)));
  results->UpdateValueOf("life_time_error", std::to_string(f->GetParError(2)));
  results->UpdateValueOf("normalization", std::to_string(f->GetParameter(1)));
  results->UpdateValueOf("normalization_error", std::to_string(f->GetParError(1)));
  results->UpdateValueOf("constant", std::to_string(f->GetParameter(0)));
  results->UpdateValueOf("constant_error", std::to_string(f->GetParError(0)));
  results->UpdateValueOf("total_number_events", std::to_string(h->Integral(0, total_bins)));
  results->UpdateValueOf("number_fitted_events", std::to_string(h->Integral(h->FindBin(xmin_fit), total_bins)));

  // ############### DRAW THE CHI² ###########################
  TCanvas *c_chi = new TCanvas("c_chi", "c_chi", 1000, 650);
  c_chi->SetLeftMargin(0.13);
  c_chi->SetRightMargin(0.13);
  c_chi->SetBottomMargin(0.11);
  c_chi->cd();
  c_chi->SetGrid();
  c_chi->SetLogz();
  TF2 *f_chi = new TF2("f_chi", chi_sqr, f->GetParameter(1) - 40 * f->GetParError(1), f->GetParameter(1) + 40 * f->GetParError(1), f->GetParameter(2) - 40 * f->GetParError(2), f->GetParameter(2) + 40 * f->GetParError(2), 1);
  f_chi->SetParameter(0, f->GetParameter(0));
  f_chi->SetNpx(100);
  f_chi->SetNpy(100);

  f_chi->GetXaxis()->SetTitle("normalization parameter");
  f_chi->GetYaxis()->SetTitle("#tau [ns]");
  f_chi->GetZaxis()->SetTitle("#chi^{2}");
  f_chi->SetTitle("Correlation plot of normalization parameter and #tau");

  f_chi->GetXaxis()->SetLabelSize(0.047);
  f_chi->GetXaxis()->SetLabelOffset(0.015);
  f_chi->GetXaxis()->SetTitleSize(0.05);
  f_chi->GetXaxis()->SetTitleOffset(1.01);

  f_chi->GetYaxis()->SetLabelSize(0.047);
  f_chi->GetYaxis()->SetTitleSize(0.05);
  f_chi->GetYaxis()->SetTitleOffset(1.3);

  f_chi->GetZaxis()->SetLabelSize(0.047);
  f_chi->GetZaxis()->SetTitleSize(0.05);
  f_chi->GetZaxis()->SetNdivisions(5, 10, 20, false);
  f_chi->GetZaxis()->SetTitleOffset(0.5);

  f_chi->Draw("CONT4Z");
  c_chi->SaveAs("./correlation_chi_sqr.png"); // I can't export a pdf since from there it can be seed that the plot is made point
                                              // by point. If I use png the resolution is not enought
                                              // and the overall plot looks smoother.

  return f;
}

void Fit_various_intervalas(int min_left, int max_left, int step, int max_right)
{
  std::ofstream out("output_intervalas.txt");
  for (int i = min_left; i <= max_left; i = i + step)
  {
    TF1 *res_tmp = FitData("./data/total.dat", i, max_right, 75);
    out << i << "-" << max_right << ":\t" << res_tmp->GetParameter(2) << "\t" << res_tmp->GetParError(2) << "\t" << res_tmp->GetChisquare() / res_tmp->GetNDF() << std::endl;
  }
  return;
}
