// Copyright 2023 nicolò salimbeni andrea de vita
#include <fcntl.h>

#include <iostream>
#include <string>
#include <vector>

#include "./../../../include/Event.h"
#include "./../../../include/InfoAcq.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TList.h"
#include "TObject.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

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
 -GetDarkCounts: return an int with the number of dark counts

 -supress_stdout: suppres printf outputs
 -resum_stdout: resume printf outputs
 */

/////////////// FUNCTION DECLARATIONS /////////////////
std::vector<std::string> list_files(const char *dirname = "C:/root/folder/",
                                    const char *ext     = ".root");
TH1F                    *GetEvent(TFile *input_file, int i, int channel = chB);
int                      GetNEvents(TFile *input_file);
int                      GetDarkCounts(TH1F *h);

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
  // loop over the file
  for (auto file_name : files_names) {
    // loop over the events inside this file
    TFile *root_file     = TFile::Open(file_name.c_str());
    int    number_events = GetNEvents(root_file);
    for (int i = 0; i < number_events; i++) {
      TH1F *h = GetEvent(root_file, i);
      h->Scale(-1, "nosw2");
      //
      //
      // TODO(nicolo): find number of peaks in each file
      //
      //
    }
  }
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

int GetDarkCounts(TH1F *h) {
  int counts = 0;

  return counts;
}
