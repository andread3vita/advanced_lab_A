// Copyright 2023 Andrea De Vita - Nicolò Salimbeni

// landauAnalysis.cc
// Created on: Dec 23, 2023

// cpp libraries
#include <filesystem>
#include <fstream>
#include <ostream>
#include <stdint.h>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// root libraries
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSpectrum.h"
#include "TSystem.h"
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

// custom libraries
#include "../../../../include/AnUtil.h"
#include "../../../../include/Event.h"
#include "../../../../include/InfoAcq.h"
#include "../../../../src/InfoAcq.cc"

/*
 - landauSummaryFile:

 - adc_to_mv:

 - ReadTree:

 - landauFit:

*/

//////////// MAIN FUNCTION DECLARATIONS ///////////////////
void landauSummaryFile(const char *folderPath);

//////////// MINOR FUNCTION DECLARATIONS ///////////////////
float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue);

TH1F *ReadTree(const char *fileName, bool negative, int channel = chB, Long64_t nEvtMax = -1, int numBins = 100);

void landauFit(const char *fileName, int numBins);

//////////// FUNCTION DEFINITIONS ///////////////////
void landauSummaryFile(const char *folderPath)
{
    // Convert char pointers to strings for easier manipulation
    std::string path = folderPath;

    // Vector to store file names
    std::vector<std::string> fileNames;

    // Open a file to store fit results
    std::ofstream resultFile("./../results/landauResult.txt");
    resultFile << "SIPM\tConstant\tNormalization\tmean\tsigma\tErrConst\tErrNorm\tErrMean\tErrSigma\tQuantileProb"
               << std::endl;

    // Iterate through files in the specified directory
    for (const auto &entry : fs::directory_iterator(path))
    {
        // Exclude hidden files (those starting with '.')
        if (entry.path().filename().string()[0] != '.')
        {
            fileNames.push_back(entry.path().filename().string());
        }
    }

    // Sort file names for consistent processing order
    std::sort(fileNames.begin(), fileNames.end());

    std::vector<int> numBins = {120, 110, 120, 85, 60};

    // Process each file in the directory
    for (int i = 0; i < fileNames.size(); i++)
    {
        // Construct the full path to the current file
        std::string fullpath = path + "/" + fileNames[i];
        const char *file_path = fullpath.c_str();

        // Read data from a tree and create a histogram
        TH1F *histo = ReadTree(file_path, "chB", TRUE, -1, numBins[i]);

        // Create a flipped Landau function to fit the histogram
        TF1 *flippedLandauFit = new TF1("flippedLandauFit", "[3]+landau(-x)", -600, 0);

        // Set parameter names for the fit function
        flippedLandauFit->SetParNames("normalization", "mean", "sigma", "constant");

        // Set initial parameter values for the fit
        flippedLandauFit->SetParameter(0, 1000); // Amplitude
        flippedLandauFit->SetParameter(1, 100);  // Mean
        flippedLandauFit->SetParameter(2, 2);    // Sigma

        // Set limits for the constant parameter
        flippedLandauFit->SetParLimits(3, 0, 10000);

        // Fit the histogram with the flipped Landau function using the RM option
        histo->Fit("flippedLandauFit", "RMNQ");

        // Fill the output file
        resultFile << fileNames[i].substr(7, 4) << "\t";
        resultFile << flippedLandauFit->GetParameter(3) << "\t";
        resultFile << flippedLandauFit->GetParameter(0) << "\t";
        resultFile << flippedLandauFit->GetParameter(1) << "\t";
        resultFile << flippedLandauFit->GetParameter(2) << "\t";
        resultFile << flippedLandauFit->GetParError(3) << "\t";
        resultFile << flippedLandauFit->GetParError(0) << "\t";
        resultFile << flippedLandauFit->GetParError(1) << "\t";
        resultFile << flippedLandauFit->GetParError(2) << "\t";
        resultFile << flippedLandauFit->GetProb() << std::endl;
    }
}

float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue)
{
    uint16_t inputRanges[12] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000};

    return (raw * inputRanges[rangeIndex]) * 1. / maxADCValue;
}

TH1F *ReadTree(const char *fileName, bool negative, int channel = chB, Long64_t nEvtMax = -1, int numBins = 100)
{
    // struct declaration
    InfoAcq::chSettings chSet1;
    InfoAcq::chSettings chSet2;
    InfoAcq::samplingSettings sampSet;

    // event variables declaration
    unsigned long long ID;
    int samplesStored;
    long long triggerInstant;
    short timeUnit;
    short *sample0;
    short *sample1;

    // Open the file in read mode
    TFile *input_file = new TFile(fileName, "READ");

    // read the trees
    TTree *treeCh = (TTree *)input_file->Get("Channels");
    TTree *treeSamp = (TTree *)input_file->Get("SampSets");
    TTree *treeEvt = (TTree *)input_file->Get("Event");
    // TFile->Get() restituisce un oggetto generico che va
    // convertito esplicitamente anteponendo (TTree*)

    // prelevo i branch con le info e li associo alle struct
    treeCh->SetBranchAddress("Ch1", &chSet1.enabled);
    treeCh->SetBranchAddress("Ch2", &chSet2.enabled);
    treeSamp->SetBranchAddress("Settings", &sampSet.max_adc_value);

    // leggo le entries
    // dichiaro l'oggetto InfoAcq e lo riempio
    InfoAcq *info = new InfoAcq();
    treeCh->GetEntry(0);
    treeSamp->GetEntry(0);
    info->FillSettings(&chSet1, &chSet2, &sampSet);

    InfoAcq::chSettings chSet;
    chSet = channel == chA ? chSet1 : chSet2;

    if (chSet.enabled == false)
    {
        cout << " channel " << channel << " not enabled, stopping" << endl;
    }
    // imposto i branches per gli eventi
    sample0 = new short[sampSet.samplesStoredPerEvent];
    sample1 = new short[sampSet.samplesStoredPerEvent];
    treeEvt->SetBranchAddress("ID", &ID);
    treeEvt->SetBranchAddress("nSamp", &samplesStored);
    treeEvt->SetBranchAddress("Instant", &triggerInstant);
    treeEvt->SetBranchAddress("TimeUnit", &timeUnit);
    treeEvt->SetBranchAddress("WaveformsA", &sample0[0]);
    treeEvt->SetBranchAddress("WaveformsB", &sample1[0]);

    Long64_t nEvt = treeEvt->GetEntries();
    if (nEvtMax >= 0 && nEvtMax < nEvt)
    {
        nEvt = nEvtMax;
    }
    float maximum = 0.0;
    float minimum = 0.0;

    // energy spectrum

    // float xmin= negative?
    // adc_to_mv(sampSet.max_adc_value,chSet.range,-1*sampSet.max_adc_value) : 0 ;
    // float xmax = negative? 0 :
    // adc_to_mv(sampSet.max_adc_value,chSet.range,sampSet.max_adc_value) ;
    float xmin = -700;
    float xmax = 50;
    TH1F *spectrumMaximum = new TH1F("hMax", "Signal peaks distribution", numBins, xmin, xmax);

    for (Long64_t index = 0; index < nEvt; index++)
    {
        treeEvt->GetEntry(index);
        for (int ii = 0; ii < sampSet.samplesStoredPerEvent; ii++)
        {
            short sample = channel == chA ? sample0[ii] : sample1[ii];
            float value = adc_to_mv(sample, chSet.range, sampSet.max_adc_value);
            if (value > maximum)
                maximum = value;
            if (value < minimum)
                minimum = value;
        }

        // std::cout << minimum - maximum << std::endl;
        spectrumMaximum->Fill(negative ? minimum - maximum : maximum - minimum);
        maximum = 0.0;
        minimum = 0.0;
    }
    return (spectrumMaximum);
}

void landauFit(const char *fileName, int numBins)
{
    // Read data from a tree and create a histogram
    TH1F *histo = ReadTree(fileName, "chB", TRUE, -1, numBins);

    // Create a flipped Landau function to fit the histogram
    TF1 *flippedLandauFit = new TF1("flippedLandauFit", "[3]+landau(-x)", -600, 0);

    // Set parameter names for the fit function
    flippedLandauFit->SetParNames("normalization", "mean", "sigma", "constant");

    // Set initial parameter values for the fit
    flippedLandauFit->SetParameter(0, 1000); // Amplitude
    flippedLandauFit->SetParameter(1, 60);   // Mean
    flippedLandauFit->SetParameter(2, 20);   // Sigma

    // Set limits for the constant parameter
    flippedLandauFit->SetParLimits(3, 0, 10000);

    // Fit the histogram with the flipped Landau function using the RM option
    histo->Fit("flippedLandauFit", "RM");
}
