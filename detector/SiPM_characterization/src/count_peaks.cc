#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TSpectrum.h"
int Smooth(int n = 3) {
  TFile *f = new TFile("./../test_peaks.root");
  TH1F  *h = reinterpret_cast<TH1F *>(f->Get("h"));

  const Int_t nbins = h->GetNbinsX();
  Double_t    xmin  = h->GetXaxis()->GetXmin();
  Double_t    xmax  = h->GetXaxis()->GetXmax();
  Double_t    source[nbins];
  gROOT->ForceStyle();

  double b_width = h->GetBinWidth(1);
  double min     = h->GetMinimum();
  for (int i = 0; i < nbins; i++) {
    h->SetBinContent(i, h->GetBinContent(i) - min);
  }

  for (int i = 0; i < nbins; i++) source[i] = h->GetBinContent(i);
  h->Draw("L");

  TSpectrum *s = new TSpectrum();

  TH1F *smooth = new TH1F("smooth", "smooth 7 iterations", nbins, xmin, xmax);
  smooth->SetLineColor(kRed);

  s->SmoothMarkov(source, nbins, n);  // 3, 7, 10
  for (int i = 0; i < nbins; i++) smooth->SetBinContent(i, source[i]);
  smooth->Draw("L SAME");

  TSpectrum *sp     = new TSpectrum();
  int        counts = sp->Search(smooth, 2, "", 0.01);
  return counts;
}
