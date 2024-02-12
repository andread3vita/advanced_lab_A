// Copyright 2024 Andrea De Vita - Nicolò Salimbeni

// noiseCheck.cc
// Created on: Feb 3, 2024

// cpp libraries
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using std::string;

std::random_device rd;
std::mt19937 gen(rd()); // Mersenne Twister 19937 generator

// root libraries
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"

#include "./../../../../include/AnUtil.h"
#include "./../../../../include/StateFile.h"
#include "./../../../../src/StateFile.cc"

/*


*/

//////////// MAIN FUNCTION DECLARATIONS ///////////////////
void windowsCheck(const char *datafile, const int window_size, const double alpha = 0.05, bool verbose = TRUE,
                  const int N = 1000000);

//////////// MINOR FUNCTION DECLARATIONS ///////////////////
double upperConfidenceLimit(std::vector<double> samples, double alpha);
double CL_probability(const int N, const int window_size, const double alpha);

//////////// FUNCTION DEFINITIONS ///////////////////

void windowsCheck(const char *datafile, const int window_size, const double alpha = 0.05, bool verbose = TRUE,
                  const int N = 1000000)

{

    // ANSI escape codes for text colors
    const char *RED_TEXT = "\033[1;31m";
    const char *RESET_COLOR = "\033[0m";

    // Open data file
    std::ifstream infile(datafile);

    // Open Arietta calibration parameters file
    StateFile *calibration_parameters = new StateFile("../../../detector/arietta/results.txt");

    double a = std::stod(calibration_parameters->ValueOf("m"));
    double b = std::stod(calibration_parameters->ValueOf("q"));

    std::vector<double> converted_values; // Vector to store converted values

    // Read rows from the infile, convert, and store in the vector
    double x;
    while (infile >> x)
    {
        double converted_value = a * x + b;
        converted_values.push_back(converted_value);
    }

    std::vector<double> probabilities;
    probabilities.reserve(converted_values.size());

    double delta = a + b;
    double lifetime = 2196.981;
    for (double x : converted_values)
    {
        double x1 = x - delta;
        double x2 = x + delta;
        double prob = exp(-x1 / lifetime) - exp(-x2 / lifetime);
        probabilities.push_back(prob);
    }

    double CL_prob = CL_probability(N, window_size, alpha);

    std::cout << "Confidence level probability value:" << CL_prob << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    sleep(5);

    // Print sliding window of 5 elements
    int suspect = 0;
    for (size_t i = 0; i <= probabilities.size() - window_size; ++i)
    {

        double prob = 0;
        for (size_t j = i; j < i + window_size; ++j)
        {
            prob += probabilities[j];
        }

        if (prob < CL_prob)
        {
            if (verbose)
            {
                std::cout << "Window:" << i << "\t"
                          << "probWin:" << round(prob * 1000) / 1000 << std::endl;
            }
        }
        else
        {
            suspect += 1;
            if (verbose)
            {
                std::cout << RED_TEXT << "Window:" << i << "\t"
                          << "probWin:" << round(prob * 1000) / 1000 << RESET_COLOR << std::endl;
            }
        }
    }

    std::cout << RED_TEXT << "\n\nSuspecious windows:" << suspect << RESET_COLOR << std::endl;

    delete calibration_parameters; // Free memory for calibration_parameters
}

double upperConfidenceLimit(std::vector<double> samples, double alpha)
{
    // Compute mean and standard deviation of the samples
    double sum = 0.0;
    for (double sample : samples)
    {
        sum += sample;
    }
    double mean = sum / samples.size();

    double sum_squared_diff = 0.0;
    for (double sample : samples)
    {
        double diff = sample - mean;
        sum_squared_diff += diff * diff;
    }
    double variance = sum_squared_diff / (samples.size() - 1);
    double std_dev = sqrt(variance);

    double Z_score = TMath::NormQuantile(1 - alpha);

    double upper_limit = Z_score * std_dev + mean;

    return upper_limit;
}

double CL_probability(const int N, const int window_size, const double alpha)
{
    // Open Arietta calibration parameters file
    StateFile *calibration_parameters = new StateFile("../../../detector/arietta/results.txt");
    double lifetime = 2196.981;

    std::exponential_distribution<double> distribution(1 / lifetime);

    // Generate N exponential distribution samples
    std::vector<double> samples;
    for (int i = 0; i < N; ++i)
    {
        double sample = distribution(gen);
        samples.push_back(sample);
    }

    double a = std::stod(calibration_parameters->ValueOf("m"));
    double b = std::stod(calibration_parameters->ValueOf("q"));

    std::vector<double> probabilities;
    probabilities.reserve(samples.size());

    double delta = a + b;
    for (double x : samples)
    {
        double x1 = x - delta;
        double x2 = x + delta;
        double prob = exp(-x1 / lifetime) - exp(-x2 / lifetime);
        probabilities.push_back(prob);
    }

    std::vector<double> window_probabilities;
    // Print sliding window of 5 elements
    for (size_t i = 0; i <= probabilities.size() - window_size; ++i)
    {
        double prob = 0;
        for (size_t j = i; j < i + window_size; ++j)
        {
            prob += probabilities[j];
        }
        window_probabilities.push_back(prob);
    }

    double upper_limit = upperConfidenceLimit(window_probabilities, alpha);

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);

    // Create a histogram to draw
    TH1D *histogram = new TH1D("histogram", "Histogram of Window Probabilities", 70, 0, 1);

    // Fill the histogram with the values from the vector
    for (size_t i = 0; i < window_probabilities.size(); ++i)
    {
        histogram->Fill(window_probabilities[i]); // i+0.5 to center the bin at i
    }

    // Draw the histogram on the canvas
    histogram->Draw();

    // Draw a vertical line corresponding to the upper limit
    TLine *line = new TLine(upper_limit, 0, upper_limit, histogram->GetMaximum());
    line->SetLineColor(kRed); // Set line color to red
    line->Draw("same");       // Draw line on the same canvas

    return upper_limit;
}