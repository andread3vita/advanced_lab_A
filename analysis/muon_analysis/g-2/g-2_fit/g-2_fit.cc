#include "./../../../../include/AnUtil.h"
#include "./../../../../include/FolderManager.h"
#include "./../../../../include/GrUtil.h"
#include "./../../../../include/StateFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TSystem.h"

#include <RtypesCore.h>
#include <TH1.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

/*
 *
 */
// ======= HELP FUNCTIONS AND OBJECTS =========
TH1D    *LoadData(std::string path_to_data_folder, int nbins, int xmin, int xmax, int dataset_initial, int dataset_final);
TH1D    *GetResiduals(const TH1D *h, const TF1 *f, double xmin, double xmax);
TF1     *f_no_B;
Double_t f_fit(Double_t *x, Double_t *par)
{
  double N   = f_no_B->GetParameter(1);
  double tau = f_no_B->GetParameter(2);
  return f_no_B->Eval(x[0]) + N * par[0] * TMath::Exp(-x[0] / tau) * TMath::Cos(par[1] * x[0] + par[2]);
}
// ======= ANALYSIS FUNCTIONS ==========

// physical quantities
double e              = 1.602176634e-19;               // Coulomb
double c              = 299792458;                     // m/s
double muon_mass      = (105.6583755e6) * e / (c * c); // kg
double magnetic_filed = 5.4e-3;                        // tesla

std::string G2_estimation(int bin_number_B, int dataset_initial = 0, int dataset_final = 0)
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
  double bin_width_no_B       = (range_xmax - range_xmin) * 1.0 / total_bins;

  double constant      = std::stod(input_parameters_file.ValueOf("constant"));
  double normalization = std::stod(input_parameters_file.ValueOf("normalization"));
  double life_time     = std::stod(input_parameters_file.ValueOf("life_time"));

  // ================== PLOT WITHOUT MAGNETIC FIELD ==================================
  TCanvas *c_no_cos = new TCanvas("c_no_cos", "c_no_cos", 1800, 680);
  c_no_cos->Divide(2, 1);
  c_no_cos->cd(1);
  gPad->SetGrid();
  gPad->SetLogy();

  int    bin_B       = bin_number_B;
  double bin_width_B = (range_xmax - range_xmin) * 1.0 / bin_B;
  TH1D  *h           = LoadData("./data", bin_B, range_xmin, range_xmax, dataset_initial, dataset_final);
  h->SetStats(kFALSE);
  h->SetTitle("N(t)=c+N_{0}e^{-t/#tau} over dataset with B");
  h->GetXaxis()->SetTitle("time [ns]");
  h->GetYaxis()->SetTitle(GrUtil::HLabel(h, "ns").c_str());
  GrUtil::SetHTextSize(h);
  h->GetXaxis()->SetTitleOffset(0.96);
  h->GetYaxis()->SetTitleOffset(1.06);
  h->GetXaxis()->SetMaxDigits(3);
  h->GetXaxis()->SetNdivisions(8, 4, 0, kFALSE);
  h->Draw();

  f_no_B = new TF1("f_no_B", "[0]+[1]*exp(-x/[2])", xmin_fit, xmax_fit);

  f_no_B->SetParName(0, "constant");
  f_no_B->SetParName(1, "normalization");
  f_no_B->SetParName(2, "life time [ns]");

  // we have to adjust the normalization of the function since the parameters in input
  // are estimated using another dataset.

  int n_events_in_range = h->Integral(h->FindBin(xmin_fit), h->FindBin(xmax_fit));

  double k = 1.0 * n_events_in_range * bin_width_B / (number_fitted_events * bin_width_no_B);

  f_no_B->FixParameter(0, constant * 1.0 * k);
  f_no_B->FixParameter(1, normalization * 1.0 * k);
  f_no_B->FixParameter(2, life_time);

  f_no_B->Draw("SAME");

  auto h_legend = new TLegend(0.25, 0.7, 0.9, 0.9);
  h_legend->AddEntry(h, "data with B", "l");
  h_legend->AddEntry(f_no_B, "function obtained without B", "l");
  h_legend->Draw();

  c_no_cos->cd(2);
  gPad->SetGrid();
  TH1D *h_res = GetResiduals(h, f_no_B, xmin_fit, xmax_fit);
  h_res->SetStats(kFALSE);
  h_res->GetXaxis()->SetTitle("time [ns]");
  h_res->SetTitle("Residuals of the plot on the left");
  h_res->GetYaxis()->SetTitle(GrUtil::HLabel(h_res, "ns").c_str());
  h_res->GetXaxis()->SetTitleOffset(0.96);
  h_res->GetYaxis()->SetTitleOffset(1.06);
  h_res->GetXaxis()->SetMaxDigits(3);
  h_res->GetXaxis()->SetNdivisions(8, 4, 0, kFALSE);
  GrUtil::SetHTextSize(h_res);

  // get maximum and minimum to get a good plotting range
  double res_max = h_res->GetMaximum();
  double res_min = h_res->GetMinimum();
  h_res->GetYaxis()->SetRangeUser(res_min * 1.2, res_max * 1.2);

  h_res->Draw();

  c_no_cos->SaveAs("./g-2_fit/g-2_plot_olf_function.pdf");

  //================== PLOT WITH MAGNETIC FILED ==========================
  TCanvas *c_cos = new TCanvas("c_cos", "c_cos", 1800, 680);
  c_cos->Divide(2, 1);
  c_cos->cd(1);
  gPad->SetGrid();
  gPad->SetLogy();

  TH1D *h_cos = new TH1D("h_cos", "Fit data with B and N(t)=c+N_{0}e^{-t/#tau} [1+Acos(#omega_{a} t + #phi)]", bin_B, range_xmin, range_xmax);
  h_cos->Add(h);
  h_cos->SetStats(kTRUE);
  h_cos->GetXaxis()->SetTitle("time [ns]");
  h_cos->GetYaxis()->SetTitle(GrUtil::HLabel(h_cos, "ns").c_str());
  GrUtil::SetHTextSize(h_cos);
  h_cos->GetXaxis()->SetTitleOffset(0.96);
  h_cos->GetYaxis()->SetTitleOffset(1.06);
  h_cos->GetXaxis()->SetMaxDigits(3);
  h_cos->GetXaxis()->SetNdivisions(8, 4, 0, kFALSE);
  h_cos->Draw();

  TF1 *f_cos = new TF1("f_cos", f_fit, xmin_fit, xmax_fit, 3);
  f_cos->SetParName(0, "A");
  f_cos->SetParName(1, "#omega");
  f_cos->SetParName(2, "#phi");
  f_cos->SetParLimits(0, 0, 1);
  f_cos->SetParLimits(1, 0, 0.01);
  f_cos->SetParLimits(2, 0, 2 * TMath::Pi());
  f_cos->SetParameter(1, 4.7 * 1e-3);
  // f_cos->SetParameter(2, 0);

  // TF1 *f_cos = new TF1("f_cos", "[0]+[1]*exp(-x/[5])*(1+[2]*cos([3]*x+[4]))", xmin_fit, xmax_fit);
  // f_cos->SetParName(0, "c");
  // f_cos->SetParName(1, "N");
  // f_cos->SetParName(2, "A");
  // f_cos->SetParName(3, "#omega");
  // f_cos->SetParName(4, "#phi");
  // f_cos->SetParName(5, "#tau [ns]");
  //
  // f_cos->SetParLimits(2, 0, 1);
  // f_cos->SetParLimits(4, -TMath::Pi(), TMath::Pi());
  // f_cos->SetParLimits(3, 0, 0.00001);
  // f_cos->SetParameter(3, 1e-7);
  // f_cos->FixParameter(5, life_time);

  h_cos->Fit(f_cos, "MRQ");

  gStyle->SetOptStat("rm");
  gStyle->SetOptFit(1);

  TPaveStats *pv = (TPaveStats *)h_cos->GetListOfFunctions()->FindObject("stats");
  pv->SetX1NDC(0.45);
  pv->SetX2NDC(0.9);
  pv->SetY1NDC(0.55);
  pv->SetY2NDC(0.9);
  pv->Draw("SAME");

  c_cos->cd(2);
  gPad->SetGrid();
  TH1D *h_cos_res = GetResiduals(h_cos, f_cos, xmin_fit, xmax_fit);
  h_cos_res->SetTitle("Residuals of the plot on the left");
  h_cos_res->SetStats(kFALSE);
  h_cos_res->GetXaxis()->SetTitle("time [ns]");
  h_cos_res->GetYaxis()->SetTitle(GrUtil::HLabel(h_cos_res, "ns").c_str());
  GrUtil::SetHTextSize(h_cos_res);
  h_cos_res->GetXaxis()->SetTitleOffset(0.96);
  h_cos_res->GetYaxis()->SetTitleOffset(1.06);
  h_cos_res->GetXaxis()->SetMaxDigits(3);
  h_cos_res->GetXaxis()->SetNdivisions(8, 4, 0, kFALSE);
  h_cos_res->GetYaxis()->SetRangeUser(res_min * 1.2, res_max * 1.2);
  h_cos_res->GetYaxis()->SetRangeUser(res_min * 1.2, res_max * 1.2);
  h_cos_res->Draw();

  c_cos->SaveAs("./g-2_fit/g-2_plot_fit.pdf");
  // ====== summary on obtained value =========
  double omega_s       = f_cos->GetParameter(1) * 1e9; // 1/s
  double omega_s_error = f_cos->GetParError(1) * 1e9;  // 1/s
  double g_muon        = 2.0 * muon_mass * omega_s / (magnetic_filed * e);
  std::cout << " \n==== ANALYSIS RESULTS ====== " << std::endl;
  std::cout << "omega=" << omega_s << " s" << std::endl;
  std::cout << "muon mass=" << muon_mass << " kg" << std::endl;
  std::cout << "magnetic field=" << magnetic_filed << " T" << std::endl;
  std::cout << "electric charge=" << e << " C" << std::endl;
  std::cout << "giromagnetic factor g: " << g_muon << std::endl;

  // string outout
  std::ostringstream output_string;
  output_string << std::setw(8) << "omega = " << std::setw(12) << omega_s << " +- " << std::setw(12) << omega_s_error << "\t\t";
  output_string << std::setw(10) << "phi = " << std::setw(13) << f_cos->GetParameter(2) << " +- " << std::setw(10) << f_cos->GetParError(2) << "\t\t";
  output_string << std::setw(10) << "A = " << std::setw(10) << f_cos->GetParameter(0) << " +- " << std::setw(10) << f_cos->GetParError(0) << "\t\t";
  output_string << std::setw(6) << "chi² = " << std::setw(10) << f_cos->GetChisquare() << " / " << std::setw(4) << f_cos->GetNDF() << "\n";

  return output_string.str();
}

void G2_subset_estimation(int day_window, int bin_number_B)
{
  FolderManager           *data_folder = new FolderManager("./data");
  std::vector<std::string> files       = data_folder->GetAllObjectsPath();
  int                      first_day   = 0;
  int                      last_day    = files.size() - day_window;

  std::vector<std::string> results;

  for (int i = first_day; i < last_day; i++)
  {
    AnUtil::ProgressBar(i * 1.0 / last_day, "analysing " + files[i].substr(files[i].find("dati"), 30) + " --> " + files[i + day_window].substr(files[i + day_window].find("dati"), 30), 1, 1, false);
    gSystem->RedirectOutput("/dev/null");
    std::string s_tmp = G2_estimation(bin_number_B, i, i + day_window);
    gSystem->RedirectOutput(0);
    results.push_back(s_tmp);
  }
  AnUtil::ProgressBar(last_day * 1.0 / last_day, "analysing " + files[last_day - 1].substr(files[last_day - 1].find("dati"), 30) + " --> " + files[last_day - 1 + day_window].substr(files[last_day - 1 + day_window].find("dati"), 30), 1, 1, false);

  std::ofstream output("./g-2_fit/restults_fit.txt");
  for (auto s : results)
  {
    output << s;
  }
  return;
}

// ======= HELP FUNCTION DEFINITIONS =========
TH1D *LoadData(std::string path_to_data_folder, int nbins, int xmin, int xmax, int dataset_initial, int dataset_final)
{
  // get calibration values from arietta
  StateFile *arietta_cal = new StateFile("../../detector/arietta/results.txt");
  double     m_arietta   = std::stod(arietta_cal->ValueOf("m"));
  double     q_arietta   = std::stod(arietta_cal->ValueOf("q"));

  // declare and fill an histogram with all the data fromt the data folder
  TH1D *data_hist = new TH1D("data_hist", "data_hist", nbins, xmin, xmax);

  FolderManager *data_folder = new FolderManager("./data");
  // Fill histogram
  std::vector<std::string> files = data_folder->GetListOfObjectsPath(".dat");

  if (dataset_final == 0)
  {
    dataset_final = files.size();
  }
  for (int i = dataset_initial; i < dataset_final; i++)
  {
    AnUtil::ProgressBar((i + 1) * 1.0 / (dataset_final - dataset_initial), "loading " + files[i].substr(files[i].find("dati"), 30), 1, 1, false);
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
