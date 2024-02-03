#include "./../../../include/FolderManager.h"
#include "./../../../include/GrUtil.h"
#include "./../../../include/StateFile.h"
#include "TH1D.h"
#include <TCanvas.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TPaveStats.h"

void SpectrumPlot(int nbins = 75, double xmin = 0, double xmax = 18000)
{
  StateFile *arietta_cal = new StateFile("../../detector/arietta/results.txt");
  double     m_arietta   = std::stod(arietta_cal->ValueOf("m"));
  double     q_arietta   = std::stod(arietta_cal->ValueOf("q"));

  FolderManager *data_folder = new FolderManager("./data");

  // create the histogram
  TH1D *spectrum = new TH1D("spectrum", "spectrum", nbins, xmin, xmax);

  // Fill histogram
  std::vector<std::string> files = data_folder->GetAllObjectsPath();
  std::cout << files.size() << std::endl;
  for (std::string file : files)
  {
    std::cout << file << std::endl;
    std::ifstream input_file(file.c_str());
    int           tmp;
    while (input_file >> tmp)
    {
      spectrum->Fill(m_arietta * tmp + q_arietta);
    }
    input_file.close();
  }

  // graphs stuff
  TCanvas *c_spectrum = new TCanvas("c", "c", 900, 650);
  c_spectrum->cd();
  c_spectrum->SetGrid();
  spectrum->GetYaxis()->SetTitle(GrUtil::HLabel(spectrum, "ns", "counts").c_str());
  spectrum->GetXaxis()->SetTitle("time [ns]");
  spectrum->SetTitle("muon decay in B");
  spectrum->SetLineColor(kBlue);
  spectrum->SetLineWidth(2);
  GrUtil::SetCountDigits(spectrum);
  GrUtil::SetHTextSize(spectrum);
  spectrum->Draw();

  return;
}
