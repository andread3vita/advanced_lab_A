#include "./../../../../include/AnUtil.h"
#include "./../../../../include/StateFile.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TVectorD.h"
#include <TGraph.h>
#include <TVectorDfwd.h>
#include <string>

void calibration()
{
  StateFile *results = new StateFile("./results.txt");
  gStyle->SetOptFit();
  // import data
  std::string file_data = "./data/calibration_data.txt";
  TVectorD    x         = AnUtil::LoadVector(file_data.c_str(), 2);
  TVectorD    y         = AnUtil::LoadVector(file_data.c_str(), 1);
  TVectorD    x_err     = AnUtil::LoadVector(file_data.c_str(), 3);
  TVectorD    y_err(x.GetNoElements());

  TCanvas *c = new TCanvas("c", "c", 1000, 600);
  c->SetGrid();

  TGraphErrors *g = new TGraphErrors(x, y, x_err, y_err);
  TF1          *f = new TF1("f", "[0]+x*[1]", 0, 70);
  f->SetParameters(7, 14);

  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlue);
  g->SetMarkerSize(1.5);
  g->GetXaxis()->SetTitle("A.u. arietta");
  g->GetYaxis()->SetTitle("time [ns]");
  g->SetTitle("");

  g->GetXaxis()->SetLabelSize(0.047);
  g->GetXaxis()->SetTitleSize(0.047);
  g->GetYaxis()->SetLabelSize(0.047);
  g->GetYaxis()->SetTitleSize(0.047);
  g->GetYaxis()->SetTitleOffset(1.03);

  g->Fit(f, "R");
  g->Draw("AP");

  c->Modified();
  c->Update();

  TPaveStats *pv = (TPaveStats *)g->GetListOfFunctions()->FindObject("stats");
  pv->SetX1NDC(0.1);
  pv->SetX2NDC(0.4);
  pv->SetY1NDC(0.7);
  pv->SetY2NDC(0.9);

  results->UpdateValueOf("q", std::to_string(f->GetParameter(0)));
  results->UpdateValueOf("q_err", std::to_string(f->GetParError(0)));
  results->UpdateValueOf("m", std::to_string(f->GetParameter(1)));
  results->UpdateValueOf("m_err", std::to_string(f->GetParError(1)));

  c->SaveAs("arietta_calibration.pdf");

  return;
}
