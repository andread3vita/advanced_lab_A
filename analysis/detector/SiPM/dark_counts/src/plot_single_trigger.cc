#include "./../../../../../include/Event.h"
#include "./../../../../../include/InfoAcq.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include <iostream>

int   supress_stdout();
void  resume_stdout(int fd);
float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue);
TH1F *SmoothHistogram(TH1F *h, int n);
int   GetNPeaksManual(TH1F *h, Double_t tr, int n_search);
int   GetNPeaks(TH1F *h, Double_t tr, Double_t sigma, Option_t *opt);

TH1F *PlotSingleTrigger(std::string path, int i, int channel, double threshold, int n_markov = 3)
{
  TFile *input_file = TFile::Open(path.c_str());

  std::cout.setstate(std::ios_base::failbit);
  int suppres_print = supress_stdout();

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

  std::cout.clear();
  resume_stdout(suppres_print);
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

  TCanvas *c = new TCanvas("c", "c", 1000, 800);
  c->Divide(1, 2);
  c->cd(1);
  signal->Draw();

  // now print the smooth histogram
  TH1F *h_smooth = SmoothHistogram(signal, n_markov);
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

  c->cd(2);
  h_smooth->Draw();

  // count numbers of peaks manually
  int n_manual = GetNPeaksManual(h_smooth, threshold, 5);
  std::cout << "detected peaks manually: " << n_manual << std::endl;
  int n_root = GetNPeaks(h_smooth, 0.01, 3, "");
  std::cout << "detected peaks bt root: " << n_root << std::endl;

  return signal;
}
