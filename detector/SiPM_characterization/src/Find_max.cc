// Copyright 2023 Nicolò Salimbeni
#include <RtypesCore.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TSpectrum.h>

#include <ostream>
#include <string>
#include <vector>

#include "../include/AnUtil.h"

void Find_peaks(const std::string file_name, const std::string histogram_name,
                int n_expected_peaks, int fit_width) {
  /*
    This macro loads an histogram (called like the input histogram_name)
    from a file and then fit its peaks with gaussians.
    To make it effective we have to insert how many peaks we
    expect (n_expected_peaks).
    Then to fit them we have to choose the range in which each
    gaussian must be fitted. The parameter fit_width refers to
    the numebr of bins below and above the detected peak we'll
    fit
  */

  // Open the TFile
  TFile* file = TFile::Open(file_name.c_str());

  // Read the histogram from the file
  TH1F* histogram = dynamic_cast<TH1F*>(file->Get("histogram"));

  // Create a TSpectrum object to find peaks
  TSpectrum* spectrum = new TSpectrum(n_expected_peaks);

  // Find the peaks in the histogram using TSpectrum
  int n_peaks = spectrum->Search(histogram, 2, "", 0.01);

  // Get the positions of the peaks
  double* peak_positions = spectrum->GetPositionX();

  // We don't now the order of the peaks in the array
  // I order them from the lower to the upper
  AnUtil::Reorder(peak_positions, n_peaks);

  // Fit gaussians
  TF1*     gaussians[n_peaks];  // 4 gaussians declaration
  Double_t bin_width = histogram->GetXaxis()->GetBinWidth(1);

  for (int i = 0; i < n_peaks; i++) {
    // take mean, and range values to fit
    Double_t mean      = peak_positions[i];
    Double_t range_inf = mean - fit_width * bin_width;  // 10 bin blow the mean
    Double_t range_max = mean + fit_width * bin_width;  // 10 bin above the mean

    gaussians[i] = new TF1(Form("gaussian%d", i), "gaus", range_inf, range_max);
    gaussians[i]->SetParameter(1, mean);
    histogram->Fit(gaussians[i], "R+");
  }

  // I save all the values of the gain in a vector
  std::vector<Double_t> gains;
  for (int i = 0; i < n_peaks - 1; i++) {
    Double_t gain =
        gaussians[i + 1]->GetParameter(1) - gaussians[i]->GetParameter(1);
    // remember: 1 is the parameter corresponding to the mean

    gains.push_back(gain);
  }

  // print the mean gain
  Double_t mean_gain = AnUtil::Mean(gains, gains.size());
  Double_t rms_gain  = AnUtil::Rms(gains, gains.size());
  std::cout << "Gain mean is:\t" << mean_gain << std::endl;
  std::cout << "Gain rms is:\t" << rms_gain << std::endl;
  return;
}
