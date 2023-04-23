// Copyright 2023 nicolò salimbeni andrea de vita
#include <TSpectrum.h>
#include <TVectorDfwd.h>
#include <fcntl.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "./../../../include/Event.h"
#include "./../../../include/InfoAcq.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TList.h"
#include "TObject.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TVectorD.h"

/*
 In this file there are all the required functions to obtain the dark
 counts of left and right SiPM. Our .root files cointain in chA the left
 and in chB the right SiPM counts.
 At first we need to load the files, then we have to count the dark counts
 and the we fit the values as a function of the applied voltage.
 What we expect is a linear behaviour.

 function list:
 -list_files: return a vector with files inside a directory
 -GetEvent: return an event histogram
 -GetNEvents: return a int with the number of events
 -SmoothHistogram: return a smooth histogram without noise fluctuations
  using the markov algorithm
 -GetNPeaks: count the number of peaks in a histogram

 -supress_stdout: suppres printf outputs
 -resum_stdout: resume printf outputs
 */

/////////////// FUNCTION DECLARATIONS /////////////////
std::vector<std::string> list_files(const char *dirname = "C:/root/folder/",
                                    const char *ext     = ".root");
TH1F                    *GetEvent(TFile *input_file, int i, int channel = chB);
int                      GetNEvents(TFile *input_file);
TH1F                    *SmoothHistogram(TH1F *h, int n);
int GetNPeaks(TH1F *h, Double_t sigma = 2, Option_t *opt = "",
              Double_t tr = 0.01);

int  supress_stdout();
void resume_stdout(int fd);

/////////////// MAIN FUNCTION /////////////////////////

void DarkCountsAnalysis(const char *dirname = "./data/counts_A/",
                        const char *ext     = ".root") {
  // find the files in the proper directory
  std::vector<std::string> files_names = list_files(dirname, ext);
  std::cout << "Loaded files:" << std::endl;
  for (auto file_name : files_names) {
    std::cout << file_name << std::endl;
  }
  std::cout << std::endl;

  // define TVectorD for counts and voltages
  const int dim_array = files_names.size();
  TVectorD  counts_A(dim_array);
  TVectorD  counts_B(dim_array);
  TVectorD  voltages(dim_array);
  TVectorD  voltages_error(dim_array);

  // apply uniform distribution for voltage error
  for (int i = 0; i < dim_array; i++) {
    voltages_error[i] = 0.2 / std::sqrt(12);
  }

  // define TVectorD for the algorithm error
  TVectorD counts_A_up(dim_array);
  TVectorD counts_A_down(dim_array);
  TVectorD counts_B_up(dim_array);
  TVectorD counts_B_down(dim_array);

  // define vectors for error evaluation
  TVectorD err_A(dim_array);
  TVectorD err_B(dim_array);

  TVectorD sqr_A_gaus(dim_array);
  TVectorD sqr_B_gaus(dim_array);
  TVectorD err_A_algorithm(dim_array);
  TVectorD err_B_algorithm(dim_array);

  for (int i = 0; i < files_names.size(); i++) {
    std::string voltage_s = files_names[i];
    voltage_s.erase(0, 8);
    voltage_s.erase(3, 5);
    int    voltage      = std::stoi(voltage_s);
    double voltage_volt = voltage * 1.0 / 10;
    voltages[i]         = voltage_volt;
  }

  // loop over the files
  for (int j = 0; j < files_names.size(); j++) {
    // loop over the events inside this file
    TFile *root_file     = TFile::Open(files_names[j].c_str());
    int    number_events = GetNEvents(root_file);
    for (int i = 0; i < 5; i++) {
      TH1F *h_chA = GetEvent(root_file, i, chA);
      TH1F *h_chB = GetEvent(root_file, i, chB);

      // find the time windows in each event in micro seconds
      double x_min_A       = h_chA->GetXaxis()->GetXmin();  // in ns
      double x_max_A       = h_chA->GetXaxis()->GetXmax();  // in ns
      double time_window_A = (x_max_A - x_min_A) / 1000;    // in micro sec

      double x_min_B       = h_chB->GetXaxis()->GetXmin();  // in ns
      double x_max_B       = h_chB->GetXaxis()->GetXmax();  // in ns
      double time_window_B = (x_max_B - x_min_B) / 1000;    // in micro sec

      // we look for max instead of min, it's easier
      h_chA->Scale(-1, "nosw2");
      h_chB->Scale(-1, "nosw2");

      // smooth histograms
      // 3 iterations for our counts, 2 4 for algorithm error evaluation
      TH1F *h_chA_smooth      = SmoothHistogram(h_chA, 3);
      TH1F *h_chA_smooth_up   = SmoothHistogram(h_chA, 4);
      TH1F *h_chA_smooth_down = SmoothHistogram(h_chA, 2);

      TH1F *h_chB_smooth      = SmoothHistogram(h_chB, 3);
      TH1F *h_chB_smooth_up   = SmoothHistogram(h_chA, 4);
      TH1F *h_chB_smooth_down = SmoothHistogram(h_chA, 2);

      // count the dark events and add the value in the vectors
      int n_chA      = GetNPeaks(h_chA_smooth);
      int n_chA_up   = GetNPeaks(h_chA_smooth_up);
      int n_chA_down = GetNPeaks(h_chA_smooth_down);

      int n_chB      = GetNPeaks(h_chB_smooth);
      int n_chB_up   = GetNPeaks(h_chB_smooth_up);
      int n_chB_down = GetNPeaks(h_chB_smooth_down);

      counts_A[j] += (n_chA / time_window_A) / number_events;
      counts_B[j] += (n_chB / time_window_B) / number_events;

      sqr_A_gaus[j] += std::pow((n_chA / time_window_A), 2);
      sqr_B_gaus[j] += std::pow((n_chB / time_window_B), 2);

      counts_A_up[j] += (n_chA_up / time_window_A) / number_events;
      counts_B_up[j] += (n_chB_up / time_window_B) / number_events;

      counts_A_down[j] += (n_chA_down / time_window_A) / number_events;
      counts_B_down[j] += (n_chB_down / time_window_B) / number_events;
    }
    // now the position j of vectors is completed
    // we have to obtain the proper error for out counts
    err_A_algorithm[j] =
        std::abs(counts_A_up[j] - counts_A_down[j]) / std::sqrt(12);
    err_B_algorithm[j] =
        std::abs(counts_B_up[j] - counts_B_down[j]) / std::sqrt(12);

    sqr_A_gaus[j] = std::sqrt(
        std::abs(sqr_A_gaus[j] / number_events - pow(counts_A[j], 2)));
    sqr_B_gaus[j] = std::sqrt(
        std::abs(sqr_B_gaus[j] / number_events - pow(counts_B[j], 2)));

    err_A[j] =
        std::sqrt(std::pow(sqr_A_gaus[j], 2) + std::pow(err_A_algorithm[j], 2));
    err_B[j] =
        std::sqrt(std::pow(sqr_B_gaus[j], 2) + std::pow(err_B_algorithm[j], 2));
  }

  // print the graps
  TGraphErrors *g_A =
      new TGraphErrors(voltages, counts_A, voltages_error, err_A);

  TGraphErrors *g_B =
      new TGraphErrors(voltages, counts_B, voltages_error, err_B);

  TCanvas *c = new TCanvas("c", "c", 1600, 680);
  c->Divide(2, 1);
  c->cd(1);
  g_A->Draw("AP");

  c->cd(2);
  g_B->Draw("AP");
  return;
}

//////////// MINOR FUNCTION DEFINITIONS ///////////////////

std::vector<std::string> list_files(const char *dirname, const char *ext) {
  /*
    This function take as input the path to a folder and a file extention
    like .root
    Then it returns a std::vector containing std::string with the name of all
    the files inside that folder with the desired extention
  */

  TSystemDirectory         dir(dirname, dirname);
  TList                   *files = dir.GetListOfFiles();
  std::vector<std::string> v;
  if (files) {
    TSystemFile *file;
    TString      fname;
    TIter        next(files);
    for (TObject *obj : *files) {
      file  = reinterpret_cast<TSystemFile *>(obj);
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
        v.push_back(fname.Data());
      }
    }
  }
  return v;
}

int GetNEvents(TFile *input_file) {
  TTree *treeEvt = reinterpret_cast<TTree *>(input_file->Get("Event"));
  return treeEvt->GetEntries();
}

float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue) {
  // function defined for GetEvent, I don't know what it does
  uint16_t inputRanges[12] = {10,   20,   50,   100,   200,   500,
                              1000, 2000, 5000, 10000, 20000, 50000};

  return (raw * inputRanges[rangeIndex]) * 1. / maxADCValue;
}

TH1F *GetEvent(TFile *input_file, int i, int channel) {
  // diable all couts defined in the functions
  std::cout.setstate(std::ios_base::failbit);
  int suppres_print = supress_stdout();

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

  // controllo a schermo la struttura del file
  input_file->Print();

  // leggo i trees
  TTree *treeCh   = reinterpret_cast<TTree *>(input_file->Get("Channels"));
  TTree *treeSamp = reinterpret_cast<TTree *>(input_file->Get("SampSets"));
  TTree *treeEvt  = reinterpret_cast<TTree *>(input_file->Get("Event"));
  TTree *treeRTI  = reinterpret_cast<TTree *>(input_file->Get("RTI"));
  // TFile->Get() restituisce un oggetto generico che va
  // convertito esplicitamente anteponendo (TTree*)

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

  TH1F *signal =
      new TH1F("Event Plot", "Event Plot", sampSet.samplesStoredPerEvent,
               -sampSet.preTrig * sampSet.timeIntervalNanoseconds,
               (sampSet.samplesStoredPerEvent - sampSet.preTrig) *
                   sampSet.timeIntervalNanoseconds);

  for (int jj = 0; jj < sampSet.samplesStoredPerEvent; jj++)
    signal->SetBinContent(
        jj, adc_to_mv(sample[jj], chSet.range, sampSet.max_adc_value));

  // Grafici
  TCanvas *c0 = new TCanvas("c0");

  c0->cd(3);
  signal->SetXTitle("Instant (ns)");
  signal->SetYTitle("Amplitude (mV)");
  signal->Draw();

  // reactivate all the couts
  std::cout.clear();
  resume_stdout(suppres_print);
  return signal;
}

int supress_stdout() {
  fflush(stdout);

  int ret    = dup(1);
  int nullfd = open("/dev/null", O_WRONLY);
  // check nullfd for error omitted
  dup2(nullfd, 1);
  close(nullfd);

  return ret;
}

void resume_stdout(int fd) {
  fflush(stdout);
  dup2(fd, 1);
  close(fd);
}

TH1F *SmoothHistogram(TH1F *h, int n) {
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
  for (int i = 0; i < nbins; i++) {
    h->SetBinContent(i, h->GetBinContent(i) - min);
  }
  for (int i = 0; i < nbins; i++) source[i] = h->GetBinContent(i);

  TSpectrum *s = new TSpectrum();

  TH1F *smooth = new TH1F("smooth", "smooth 4 iterations", nbins, xmin, xmax);
  smooth->SetLineColor(kRed);

  s->SmoothMarkov(source, nbins, n);  // 3, 7, 10
  for (int i = 0; i < nbins; i++) smooth->SetBinContent(i, source[i]);

  return smooth;
}

int GetNPeaks(TH1F *h, Double_t sigma, Option_t *opt, Double_t tr) {
  TSpectrum *sp     = new TSpectrum();
  int        counts = sp->Search(h, sigma, opt, tr);
  return counts;
}
