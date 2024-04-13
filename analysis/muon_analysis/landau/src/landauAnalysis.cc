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
    resultFile << "SIPM\tConstant\tNormalization\tmode\tsigma\tErrConst\tErrNorm\tErrMode\tErrSigma\tQuantileProb"
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
        flippedLandauFit->SetParNames("normalization", "MPV", "sigma", "constant");

        // Set initial parameter values for the fit
        flippedLandauFit->SetParameter(0, 1300); // Amplitude
        flippedLandauFit->SetParameter(1, 100);  // Mean
        flippedLandauFit->SetParameter(2, 40);   // Sigma

        // Set limits for the constant parameter
        flippedLandauFit->SetParLimits(3, 0, 20);

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

struct Vavilov_Func
{
    Vavilov_Func()
    {
    }

    double operator()(const double *x, const double *p)
    {
        double kappa = p[0];
        double beta2 = p[1];
        return p[4] * (pdf.Pdf(-(x[0] - p[2]) / p[3], kappa, beta2));
    }

    // double mode()(const double *p)
    // {
    //     double kappa = p[0];
    //     double beta2 = p[1];
    //     return pdf.Mode(kappa,beta2);
    // }

    ROOT::Math::VavilovAccurate pdf;
};

void landauFit(const char *fileName, int numBins)
{
    TCanvas *can = new TCanvas();
    // Read data from a tree and create a histogram
    TH1F *histo = ReadTree(fileName, "chB", TRUE, -1, numBins);

    histo->SetTitle("");
    // Create a flipped Landau function to fit the histogram

    Vavilov_Func *func = new Vavilov_Func();
    TF1 *flippedLandauFit = new TF1("flippedLandauFit", func, -600, 5, 5, "Vavilov_Func");

    flippedLandauFit->SetParameters(0.3, 0.05, -70, 20, histo->GetEntries());
    flippedLandauFit->SetParLimits(1, 0, 1);
    flippedLandauFit->SetParLimits(0, 0.001, 10);

    // Fit the histogram with the flipped Landau function using the RM option
    histo->Fit("flippedLandauFit", "RML");
    TFitResultPtr r = histo->Fit(flippedLandauFit, "SRML");
    TMatrixDSym cov = r->GetCovarianceMatrix();

    // Add x and y labels
    histo->GetXaxis()->SetTitle("Signal peak [mV]");
    histo->GetYaxis()->SetTitle(Form("Counts / %.0d mV", (750) / numBins));

    ROOT::Math::VavilovAccurate f1;
    double_t lambda_max = f1.Mode(flippedLandauFit->GetParameter(0), flippedLandauFit->GetParameter(1));

    Double_t mostProbableValue = flippedLandauFit->GetMaximumX();
    Double_t errorMode =
        sqrt(pow(lambda_max * flippedLandauFit->GetParError(3), 2) + pow(flippedLandauFit->GetParError(2), 2));

    Double_t euler = 0.5772156649;
    Double_t meanValue_temp =
        euler - 1 - std::log(flippedLandauFit->GetParameter(0)) - flippedLandauFit->GetParameter(1);

    Double_t meanValue = -meanValue_temp * flippedLandauFit->GetParameter(3) + flippedLandauFit->GetParameter(2);

    // error propagation
    Double_t err_tmp =
        sqrt((cov(0, 1) + cov(1, 0)) / flippedLandauFit->GetParameter(0) +
             (cov(0, 0)) / (flippedLandauFit->GetParameter(0) * flippedLandauFit->GetParameter(0)) + (cov(1, 1)));

    Double_t errorMean =
        sqrt(pow(meanValue_temp * flippedLandauFit->GetParError(3), 2) + pow(flippedLandauFit->GetParError(2), 2) +
             pow(err_tmp * flippedLandauFit->GetParameter(3), 2));

    // Create a TLine for mean value
    TLine *meanLine = new TLine(meanValue, 0, meanValue, flippedLandauFit->Eval(meanValue));
    meanLine->SetLineColor(kBlue);
    meanLine->SetLineWidth(2);
    meanLine->SetLineStyle(2);

    // Create a TLine for most probable value
    TLine *mostProbableLine =
        new TLine(mostProbableValue, 0, mostProbableValue, flippedLandauFit->Eval(mostProbableValue));
    mostProbableLine->SetLineColor(kRed);
    mostProbableLine->SetLineWidth(2);
    mostProbableLine->SetLineStyle(2);

    // Create a legend
    TLegend *legend = new TLegend(0.15, 0.7, 0.4, 0.85);
    legend->AddEntry(histo, "Data", "l");
    legend->AddEntry(flippedLandauFit, "Landau distribution (fit)", "l");
    legend->AddEntry(meanLine, Form("Mean: %.3f #pm %.3f", meanValue, errorMean), "l");
    legend->AddEntry(mostProbableLine, Form("Mode: %.3f #pm %.3f", mostProbableValue, errorMode), "l");

    gStyle->SetOptStat(0);

    // // Draw histogram with fit function
    histo->Draw("E1");
    flippedLandauFit->Draw("same");

    // Draw vertical lines
    meanLine->Draw("same");
    mostProbableLine->Draw("same");

    // Draw legend
    legend->Draw("same");

    string outputName = fileName;
    outputName.erase(outputName.length() - 5);
    size_t pos = outputName.find("landau");
    outputName.erase(0, pos);

    outputName += ".pdf";
    outputName = "../figures/" + outputName;
    can->SaveAs(outputName.c_str());

    delete flippedLandauFit;
    delete can;
}

void efficiency_calc(const char *fileName, int numBins, double threshold)
{
    // Read data from a tree and create a histogram
    TH1F *histo = ReadTree(fileName, "chB", TRUE, -1, numBins);

    double norm = histo->Integral(0, histo->GetNbinsX());
    double integral = histo->Integral(0, histo->FindFixBin(threshold));

    std::cout << "Filename: " << fileName << "\nNorm: " << norm << "\nthreshold: " << threshold << " mV"
              << "\nEfficiency: " << integral / norm * 100 << " %" << std::endl;
}