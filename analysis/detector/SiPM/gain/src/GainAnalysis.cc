// Copyright 2023 Andrea De Vita - Nicolò Salimbeni

// GainAnalysis.cc
// Created on: Nov 27, 2023

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
#include <TTree.h>

// custom libraries
#include "./../../../../../include/AnUtil.h"
#include "./../../../../../include/Event.h"
#include "./../../../../../include/InfoAcq.h"
#include "./../../../../../src/InfoAcq.cc"

/*
 In this file there are all the required functions to obtain the gain of left and right SiPM.
 Our .root files cointain in chA informations about all peaks that have been collected
 during the data acquisition interval. At first we need to load the files,
 then we have to fit mean value and rms error of the gain and the we fit the values as a function of
 the applied voltage. What we expect is a linear behaviour.

 function list:
 - adc_to_mv: it converts digital values taken by picoScope to mV;

 - ReadTree: it read .root files and it returns a pointer of the histogram of the peaks;

 - fitGain: it read all the .root files within a specific folder and returns a .txt file with
            the results of gain analysis. At each bias voltage gain and its rms error are computed;

 - gainVoltage: it reads .txt files and fit the linear behaviour of gain vs bias voltage;
 */

//////////// FUNCTION DECLARATIONS ///////////////////
float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue);
TH1F *ReadTree(const char *fileName, bool negative, int channel = chB, Long64_t nEvtMax = -1);
void fitGain(const char *dirpath, const char *txt_file, int n_expected_peaks = 3, int fit_width = 7);

//////////// MAIN FUNCTION DEFINITION ///////////////////
void fitGain(const char *dirpath, const char *txt_file, int n_expected_peaks = 3, int fit_width = 7)
{
    /*
        INPUT:
            - dirpath: path of the folder with all the .root files
            - txt_file: name of the output file
            - n_expected_peaks: number of expected peaks in the TSpectrum analysis
            - fit_width: expected distance between close peaks
        OUTPUT:
            - none
    */

    // Set ROOT to batch mode (non-interactive)
    gROOT->SetBatch(kTRUE);

    // Convert char pointers to strings for easier manipulation
    std::string path = dirpath;
    std::string out_file = txt_file;

    // Vector to store file names
    std::vector<std::string> fileNames;

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

    // Open the output file for writing
    std::ofstream file(out_file, std::ios::out);

    // Write the header to the output file
    file << "V[dV]\tGain[mV]\tSigmaV[dV]\tSigmaGain[mv]" << std::endl;

    // Process each file in the directory
    for (int i = 0; i < fileNames.size(); i++)
    {
        // Construct the full path to the current file
        std::string fullpath = path + "/" + fileNames[i];
        const char *file_path = fullpath.c_str();

        // Read the tree from the current file
        TH1F *Hist = ReadTree(file_path, true, chA, -1);

        // Create a TSpectrum object to find peaks
        TSpectrum *spectrum = new TSpectrum(n_expected_peaks);

        // Find the peaks in the histogram using TSpectrum
        int n_peaks = spectrum->Search(Hist, 2, "", 0.01);

        // Get the positions of the peaks
        double *peak_positions = spectrum->GetPositionX();

        // Order the peaks from lower to higher positions
        AnUtil::Reorder(peak_positions, n_peaks);

        // Fit Gaussians to the peaks
        TF1 *gaussians[n_peaks];
        Double_t bin_width = Hist->GetXaxis()->GetBinWidth(1);

        for (int i = 0; i < n_peaks; i++)
        {
            Double_t mean = peak_positions[i];
            Double_t range_inf = mean - fit_width * bin_width;
            Double_t range_max = mean + fit_width * bin_width;

            gaussians[i] = new TF1(Form("gaussian%d", i), "gaus", range_inf, range_max);
            gaussians[i]->SetParameter(1, mean);
            Hist->Fit(gaussians[i], "R+QN");
        }

        // Save the gain values in a vector
        std::vector<Double_t> gains;
        for (int i = 0; i < n_peaks - 1; i++)
        {
            Double_t gain = gaussians[i + 1]->GetParameter(1) - gaussians[i]->GetParameter(1);
            gains.push_back(gain);
        }

        // Print the individual gain values
        for (int j = 0; j < gains.size(); j++)
        {
            std::cout << gains[j] << std::endl;
        }

        // Calculate and print the mean and RMS of the gain values
        Double_t mean_gain = AnUtil::Mean(gains, gains.size());
        Double_t rms_gain = AnUtil::Rms(gains, gains.size());
        std::cout << "Gain mean is:\t" << mean_gain << std::endl;
        std::cout << "Gain rms is:\t" << rms_gain << std::endl;

        // Write the results to the output file
        file << fileNames[i].substr(10, 3) << "\t" << mean_gain << "\t"
             << "3"
             << "\t" << rms_gain << std::endl;

        // Clean up dynamically allocated memory
        delete spectrum;
        delete Hist;
        for (int j = 0; j < n_peaks; j++)
        {
            delete gaussians[j];
        }
    }

    // Close the output file
    file.close();
}

//////////// MINOR FUNCTION DEFINITIONS ///////////////////
float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue)
{
    uint16_t inputRanges[12] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000};

    return (raw * inputRanges[rangeIndex]) * 1. / maxADCValue;
}

TH1F *ReadTree(const char *fileName, bool negative, int channel = chB, Long64_t nEvtMax = -1)
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
    cout << "Riempio l'oggetto INFO\n";
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
    float xmin = -100;
    float xmax = 0;
    TH1F *spectrumMaximum = new TH1F("hMax", "Maxima Distribution Spectrum", 200, xmin, xmax);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //	TH1F* spectrumMaximum = new TH1F( "hMax", "Maximum Spectrum", 256,
    // xmin,xmax );
    cerr << " Range   : " << xmax - xmin << " mV " << endl;
    cerr << " #Events : " << nEvt << endl;
    cerr << " #Samples: " << sampSet.samplesStoredPerEvent << endl;
    cerr << " Timestamp: " << sampSet.timeIntervalNanoseconds << " ns" << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

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
