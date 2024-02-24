#include "./../../../../include/AnUtil.h"
#include "./../../../../include/FolderManager.h"
#include "./../../../../include/StateFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"

#include <TH1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// ======= HELP FUNCTIONS DECLARATION =========
TH1D *LoadData(std::string path_to_data_folder, int nbins, int xmin, int xmax);
TH1D *GetResiduals(const TH1D *h, const TF1 *f, double xmin, double xmax);

// ======= ANALYSIS FUNCTIONS ==========

// define functions outside the main method
TF1     *f_no_B;
Double_t f_fit(Double_t *x, Double_t *par)
{
  return f_no_B->Eval(x[0]) + f_no_B->GetParameter(1) * par[0] * TMath::Cos(par[1] * x[0] + par[2]);
}

void G2_estimation()
{

  // get all vlues from analysis without magnetic field
  StateFile input_parameters_file("./../muon_lifetime/results/results_single_fit.txt");

  int    range_xmin           = std::stoi(input_parameters_file.ValueOf("range_xmin"));
  int    range_xmax           = std::stoi(input_parameters_file.ValueOf("range_xmax"));
  int    xmin_fit             = std::stoi(input_parameters_file.ValueOf("xmin_fit"));
  int    xmax_fit             = std::stoi(input_parameters_file.ValueOf("xmax_fit"));
  int    total_bins           = std::stoi(input_parameters_file.ValueOf("total_bins"));
  int    total_number_events  = std::stoi(input_parameters_file.ValueOf("total_number_events"));
  int    number_fitted_events = std::stoi(input_parameters_file.ValueOf("number_fitted_events"));
  double bin_width            = (range_xmax - range_xmin) * 1.0 / total_bins;

  double constant      = std::stod(input_parameters_file.ValueOf("constant"));
  double normalization = std::stod(input_parameters_file.ValueOf("normalization"));
  double life_time     = std::stod(input_parameters_file.ValueOf("life_time"));

  // ================== PLOT WITHOUT MAGNETIC FIELD ==================================
  TCanvas *c_no_cos = new TCanvas("c_no_cos", "c_no_cos", 1800, 680);
  c_no_cos->Divide(2, 1);
  c_no_cos->cd(1);
  gPad->SetGrid();

  TH1D *h = LoadData("./data", total_bins, range_xmin, range_xmax);
  h->Draw();

  f_no_B = new TF1("f_no_B", "[0]+[1]*exp(-x/[2])", xmin_fit, xmax_fit);

  f_no_B->SetParName(0, "constant");
  f_no_B->SetParName(1, "normalization");
  f_no_B->SetParName(2, "life time [ns]");

  // we have to adjust the normalization of the function since the parameters in input
  // are estimated using another dataset.

  int n_events_in_range = h->Integral(h->FindBin(xmin_fit), h->FindBin(xmax_fit));

  f_no_B->FixParameter(0, constant * 1.0 * n_events_in_range / number_fitted_events);
  f_no_B->FixParameter(1, normalization * 1.0 * n_events_in_range / number_fitted_events);
  f_no_B->FixParameter(2, life_time);

  f_no_B->Draw("SAME");

  c_no_cos->cd(2);
  gPad->SetGrid();
  TH1D *h_res = GetResiduals(h, f_no_B, xmin_fit, xmax_fit);
  h_res->Draw();

  //================== PLOT WITH MAGNETIC FILED ==========================
  TCanvas *c_cos = new TCanvas("c_cos", "c_cos", 1800, 680);
  c_cos->Divide(2, 1);
  c_cos->cd(1);
  gPad->SetGrid();

  TH1D *h_cos = new TH1D("h_cos", "h_cos", total_bins, range_xmin, range_xmax);
  h_cos->Add(h);
  h_cos->Draw();

  TF1 *f_cos = new TF1("f_cos", f_fit, xmin_fit, xmax_fit, 3);
  f_cos->SetParName(0, "A");
  f_cos->SetParName(1, "#omega");
  f_cos->SetParName(2, "#phi");
  f_cos->SetParLimits(0, 0, 5);
  f_cos->SetParLimits(1, 0, 20);
  f_cos->SetParLimits(2, 0, 7);

  h_cos->Fit(f_cos, "MR");

  gStyle->SetOptStat("rm");
  gStyle->SetOptFit(1);

  c_cos->cd(2);
  gPad->SetGrid();
  TH1D *h_cos_res = GetResiduals(h_cos, f_cos, xmin_fit, xmax_fit);
  h_cos_res->Draw();

  return;
}

// ======= HELP FUNCTION DEFINITIONS =========
TH1D *LoadData(std::string path_to_data_folder, int nbins, int xmin, int xmax)
{
  // get calibration values from arietta
  StateFile *arietta_cal = new StateFile("../../detector/arietta/results.txt");
  double     m_arietta   = std::stod(arietta_cal->ValueOf("m"));
  double     q_arietta   = std::stod(arietta_cal->ValueOf("q"));

  // declare and fill an histogram with all the data fromt the data folder
  TH1D *data_hist = new TH1D("data_hist", "data_hist", nbins, xmin, xmax);

  FolderManager *data_folder = new FolderManager("./data");
  // Fill histogram
  std::vector<std::string> files = data_folder->GetAllObjectsPath();
  for (int i = 0; i < files.size(); i++)
  {
    AnUtil::ProgressBar((i + 1) * 1.0 / files.size(), "loading " + files[i].substr(files[i].find("dati"), 30), 1, 1, false);
    // std::cout << i * 1.0 / files.size() << std::endl;
    std::ifstream input_file(files[i].c_str());
    int           tmp;
    while (input_file >> tmp)
    {
      data_hist->Fill(m_arietta * tmp + q_arietta);
    }
    input_file.close();
  }

  return data_hist;
}

TH1D *GetResiduals(const TH1D *h, const TF1 *f, double xmin, double xmax)
{
  double      copy_xmin = h->GetXaxis()->GetXmin();
  double      copy_xmax = h->GetXaxis()->GetXmax();
  int         nbins     = h->GetNbinsX();
  std::string name(h->GetName());
  name = name + "_residuals";

  TH1D *h_res = new TH1D(name.c_str(), name.c_str(), nbins, copy_xmin, copy_xmax);
  // Calculate residuals
  for (int i = 1; i <= h_res->GetNbinsX(); ++i)
  {
    double x = h_res->GetBinCenter(i);
    if (x < xmin || x > xmax)
    {
      continue;
    }
    double dataValue = h->GetBinContent(i);
    double fitValue  = f->Eval(x);
    double residual  = dataValue - fitValue;
    h_res->SetBinContent(i, residual);
  }

  return h_res;
}
