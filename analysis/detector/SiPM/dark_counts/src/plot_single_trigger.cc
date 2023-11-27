#include "./../../../../include/Event.h"
#include "./../../../../include/InfoAcq.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include <iostream>

float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue)
{
  // function defined for GetEvent, I don't know what it does
  uint16_t inputRanges[12] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000};

  return (raw * inputRanges[rangeIndex]) * 1. / maxADCValue;
}

TH1F *SmoothHistogram(TH1F *h, int n)
{
  /*
  This function smooth an istogram h reducing the noise statistical
  fluctuations.
  It uses the markov algorithm in the TSpectrum class. The parameter n is
  the number of iterations for the algorithm, the bigger it is the smoother
  the histogram will be. If n is too large the output histogram does not
  rapresent the inout histogram anymore, so be careful.
  */
  const Int_t nbins = h->GetNbinsX();
  Double_t    xmin  = h->GetXaxis()->GetXmin();
  Double_t    xmax  = h->GetXaxis()->GetXmax();
  Double_t    source[nbins];

  double b_width = h->GetBinWidth(1);
  double min     = h->GetMinimum();
  for (int i = 0; i < nbins; i++)
  {
    h->SetBinContent(i, h->GetBinContent(i) - min);
  }
  for (int i = 0; i < nbins; i++)
    source[i] = h->GetBinContent(i);

  TSpectrum *s = new TSpectrum();

  TH1F *smooth = new TH1F("smooth", "smooth 4 iterations", nbins, xmin, xmax);
  smooth->SetLineColor(kRed);

  s->SmoothMarkov(source, nbins, n); // 3, 7, 10
  for (int i = 0; i < nbins; i++)
    smooth->SetBinContent(i, source[i]);

  return smooth;
}

TH1F *PlotSingleTrigger(TFile *input_file, int i, int channel)
{
  // diable all couts defined in the functions
  InfoAcq::chSettings       chSet1;
  InfoAcq::chSettings       chSet2;
  InfoAcq::samplingSettings sampSet;

  unsigned long long ID;
  int                samplesStored;
  long long          triggerInstant;
  short              timeUnit;
  short             *sample0;
  short             *sample1;
  short             *sample;

  unsigned long long waveformInBlock;
  unsigned long long elapsedTime;
  unsigned long long waveformStored;

  // leggo i trees
  TTree *treeCh   = reinterpret_cast<TTree *>(input_file->Get("Channels"));
  TTree *treeSamp = reinterpret_cast<TTree *>(input_file->Get("SampSets"));
  TTree *treeEvt  = reinterpret_cast<TTree *>(input_file->Get("Event"));
  TTree *treeRTI  = reinterpret_cast<TTree *>(input_file->Get("RTI"));

  // prelevo i branch con le info e li associo alle struct
  treeCh->SetBranchAddress("Ch1", &chSet1.enabled);
  treeCh->SetBranchAddress("Ch2", &chSet2.enabled);
  treeSamp->SetBranchAddress("Settings", &sampSet.max_adc_value);

  // leggo le entries
  // dichiaro l'oggetto InfoAcq e lo riempio
  InfoAcq *info = new InfoAcq();
  std::cout << "Riempio l'oggetto INFO\n";
  treeCh->GetEntry(0);
  treeSamp->GetEntry(0);
  info->FillSettings(&chSet1, &chSet2, &sampSet);
  info->PrintInfo();
  InfoAcq::chSettings chSet;
  chSet = channel == chA ? chSet1 : chSet2;

  // imposto i branches per gli eventi
  sample0 = new short[sampSet.samplesStoredPerEvent];
  sample1 = new short[sampSet.samplesStoredPerEvent];
  sample  = channel == chA ? sample0 : sample1;
  //
  treeEvt->SetBranchAddress("ID", &ID);
  treeEvt->SetBranchAddress("nSamp", &samplesStored);
  treeEvt->SetBranchAddress("Instant", &triggerInstant);
  treeEvt->SetBranchAddress("TimeUnit", &timeUnit);
  treeEvt->SetBranchAddress("WaveformsA", &sample0[0]);
  treeEvt->SetBranchAddress("WaveformsB", &sample1[0]);
  //
  treeRTI->SetBranchAddress("WaveformInBlock", &waveformInBlock);
  treeRTI->SetBranchAddress("ElapsedTime", &elapsedTime);
  treeRTI->SetBranchAddress("WaveformStored", &waveformStored);

  int daVedere = i;

  // leggo e disegno un evento
  // dichiaro l'oggetto Event e lo riempio
  Event *evt = new Event(sampSet.samplesStoredPerEvent);
  std::cout << "Riempio l'oggetto EVENT\n";

  treeRTI->GetEntry(0);
  treeRTI->GetEntry(daVedere / waveformInBlock);
  evt->FillRTI(waveformInBlock, elapsedTime, waveformStored);
  evt->PrintRTI();

  treeEvt->GetEntry(daVedere);
  evt->FillEvent(ID, triggerInstant, timeUnit, sample, channel);

  TH1F *signal = new TH1F("Event Plot", "Event Plot", sampSet.samplesStoredPerEvent, -sampSet.preTrig * sampSet.timeIntervalNanoseconds, (sampSet.samplesStoredPerEvent - sampSet.preTrig) * sampSet.timeIntervalNanoseconds);

  for (int jj = 0; jj < sampSet.samplesStoredPerEvent; jj++)
    signal->SetBinContent(jj, adc_to_mv(sample[jj], chSet.range, sampSet.max_adc_value));

  // Grafici
  signal->SetXTitle("time [ns]");
  signal->SetYTitle("Amplitude [mV]");
  signal->SetTitle("Single trigger event");
  signal->Scale(-1, "nosw2");
  signal->GetXaxis()->SetLabelSize(0.05);
  signal->GetXaxis()->SetTitleSize(0.05);
  signal->GetYaxis()->SetLabelSize(0.05);
  signal->GetYaxis()->SetTitleSize(0.05);
  signal->GetYaxis()->SetTitleOffset(1);
  signal->SetStats(false);
  signal->SetLineWidth(2);

  TCanvas *c_not_smooth = new TCanvas("c", "c", 1600, 1200);
  c_not_smooth->cd();
  signal->Draw();

  // now print the smooth histogram
  TH1F *h_smooth = SmoothHistogram(signal, 4);
  h_smooth->SetXTitle("time [ns]");
  h_smooth->SetYTitle("Amplitude [mV]");
  h_smooth->SetTitle("Single trigger event after smooth algorithm");
  h_smooth->GetXaxis()->SetLabelSize(0.05);
  h_smooth->GetXaxis()->SetTitleSize(0.05);
  h_smooth->GetYaxis()->SetLabelSize(0.05);
  h_smooth->GetYaxis()->SetTitleSize(0.05);
  h_smooth->GetYaxis()->SetTitleOffset(1);
  h_smooth->SetStats(false);
  h_smooth->SetLineWidth(2);

  TCanvas *c_smooth = new TCanvas("c_smoot", "c_smooth", 1600, 1200);
  c_smooth->cd();
  h_smooth->Draw();

  return signal;
}
