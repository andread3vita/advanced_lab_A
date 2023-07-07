#include <TGraph.h>
#include <TVectorDfwd.h>
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TVectorD.h"
#include "./../../include/AnUtil.h"
#include <string>

void calibration(){
  // import data
  std::string file_data = "./data/calibration_data.txt";
  TVectorD x = AnUtil::LoadVector(file_data.c_str(), 2);
  TVectorD y = AnUtil::LoadVector(file_data.c_str(), 1);
  TVectorD x_err = AnUtil::LoadVector(file_data.c_str(), 3);
  TVectorD y_err(x.GetNoElements()); 

  TGraphErrors *g = new TGraphErrors(x,y,x_err,y_err);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1);
  g->GetXaxis()->SetTitle("A.u. arietta");
  g->GetYaxis()->SetTitle("time [ns]");
  g->SetTitle("arietta time calibration");
  g->Fit("pol1");

  g->Draw("AP");

  return;
}
