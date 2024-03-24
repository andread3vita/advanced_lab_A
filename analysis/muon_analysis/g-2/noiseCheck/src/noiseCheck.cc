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
#include <utility>

#include "./../../../../../include/AnUtil.h"
#include "./../../../../../include/StateFile.h"
#include "./../../../../../src/StateFile.cc"

/*
 - windowsCheck:
 - upperConfidenceLimit:
 - CL_probability:

*/

//////////// MAIN FUNCTION DECLARATIONS ///////////////////
void windowsCheck(const char *datafile, const int window_size, const double alpha = 0.01, bool verbose = TRUE,
                  const int N = 1000000);

//////////// MINOR FUNCTION DECLARATIONS ///////////////////
std::pair<double, double> ConfidenceLimit(std::vector<double> samples, double alpha);
std::pair<double, double> MC_probability(const int N, const int window_size, const double alpha);

std::vector<double> removeIndices(const std::vector<double> &vec, const std::vector<int> &indices);

//////////// FUNCTION DEFINITIONS ///////////////////
void windowsCheck(const char *datafile, const int window_size, const double alpha = 0.01, bool verbose = TRUE,
                  const int N = 1000000)
{

    // ANSI escape codes for text colors
    const char *RED_TEXT = "\033[1;31m";
    const char *BLUE_TEXT = "\033[1;34m";
    const char *RESET_COLOR = "\033[0m";

    // Open data file
    std::ifstream infile(datafile);

    // Open Arietta calibration parameters file
    StateFile *calibration_parameters = new StateFile("../../../../detector/arietta/results.txt");

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

    std::pair<double, double> CL_prob = MC_probability(N, window_size, alpha);

    std::cout << "Confidence level sum of deltaTs value (left):" << CL_prob.first << " ns" << std::endl;
    std::cout << "Confidence level sum of deltaTs value (right):" << CL_prob.second << " ns" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    sleep(10);

    // Print sliding window of 5 elements
    int suspect_left = 0;
    int suspect_right = 0;
    std::vector<int> suspecious_index_left;
    std::vector<int> suspecious_index_right;

    std::vector<double> window_probabilities;
    for (size_t i = 0; i <= converted_values.size() - window_size; ++i)
    {

        double prob = 0;
        std::vector<int> wind_indices;
        for (size_t j = i; j < i + window_size; ++j)
        {
            wind_indices.push_back(j);
            prob += converted_values[j];
        }

        window_probabilities.push_back(prob);

        if (prob < CL_prob.first)
        {
            suspect_left += 1;
            suspecious_index_left.insert(suspecious_index_left.end(), wind_indices.begin(), wind_indices.end());

            if (verbose)
            {
                std::cout << RED_TEXT << "Window:" << i << "\t"
                          << "sum of deltaTs:" << prob << " ns" << RESET_COLOR << std::endl;
            }
        }
        else if (prob > CL_prob.second)
        {
            suspect_right += 1;
            suspecious_index_right.insert(suspecious_index_right.end(), wind_indices.begin(), wind_indices.end());

            if (verbose)
            {
                std::cout << BLUE_TEXT << "Window:" << i << "\t"
                          << "sum of deltaTs:" << prob << " ns" << RESET_COLOR << std::endl;
            }
        }
        else
        {
            if (verbose)
            {
                std::cout << "Window:" << i << "\t"
                          << "sum of deltaTs:" << prob << " ns" << std::endl;
            }
        }
    }

    int expected_outliers = (converted_values.size() - window_size) * (alpha / 2);

    std::cout << "\nExpected outliers windows :" << expected_outliers << std::endl;
    std::cout << RED_TEXT << "\n\nSuspecious windows (left):" << suspect_left << RESET_COLOR << std::endl;
    std::cout << BLUE_TEXT << "Suspecious windows (right):" << suspect_right << RESET_COLOR << std::endl;

    double start = 0;
    double end = 80e3;
    int bins = 100;
    TH1D *dataWin = new TH1D("", "", bins, start, end);

    double bin_size = (end-start)/bins;
    for (double el : window_probabilities)
    {
        dataWin->Fill(el);
    }

    TCanvas *c_win = new TCanvas("c_win", "c_win", 800, 600);
    dataWin->SetStats(0);
    //c_win->SetGrid();

    // Draw the histogram on the canvas
    dataWin->Draw();

    // Draw a vertical line corresponding to the upper limit
    TLine *line = new TLine(CL_prob.first, 0, CL_prob.first, dataWin->GetMaximum());
    TLine *line2 = new TLine(CL_prob.second, 0, CL_prob.second, dataWin->GetMaximum());
    line->SetLineColor(kRed); // Set line color to red
    line->Draw("same");       // Draw line on the same canvas

    line2->SetLineColor(kRed); // Set line color to red
    line2->Draw("same");       // Draw line on the same canvas

    // Add x and y labels
    dataWin->GetXaxis()->SetTitle("pesudo-probability [ns]");
    dataWin->GetYaxis()->SetTitle(Form("counts/(%.0f ns)",bin_size));

    int Cl = (1. - alpha) * 100;

    // Create a legend
    TLegend *legend = new TLegend(0.70, 0.8, 0.87, 0.87); // Adjust the position of the legend as needed
    legend->AddEntry(line, Form("%.0d%% CL interval", Cl), "l");
    legend->AddEntry("", Form("Window Size: %.0d", window_size), "");
    legend->SetFillColor(0); // Set the legend fill color to transparent
    legend->Draw("same");

    // Write the x value over the vertical lines
    TText *text1 = new TText(CL_prob.first, dataWin->GetMaximum()+10, Form("%.2f", CL_prob.first));
    TText *text2 = new TText(CL_prob.second, dataWin->GetMaximum()+10, Form("%.2f", CL_prob.second));
    text1->SetTextAlign(22); // Center alignment
    text2->SetTextAlign(22); // Center alignment
    text1->SetTextColor(kRed); // Set text color to red
    text2->SetTextColor(kRed); // Set text color to red
    text1->SetTextSize(0.03);  // Set text size smaller
    text2->SetTextSize(0.03);  // Set text size smaller
    text1->Draw("same");
    text2->Draw("same");

    // Update the canvas
    c_win->Update();
    c_win->SaveAs("../figures/probWinDistribution.pdf");
    // std::cout << expected_outliers << std::endl;
    // std::vector<int> indices;

    // if (suspect_left > expected_outliers)
    // {
    //     indices.insert(indices.end(), suspecious_index_left.begin(), suspecious_index_left.end());
    // }
    // if (suspect_right > expected_outliers)
    // {
    //     indices.insert(indices.end(), suspecious_index_right.begin(), suspecious_index_right.end());
    // }

    // std::vector<double> result = removeIndices(converted_values, indices);

    // TCanvas *can = new TCanvas("canvas", "Canvas", 800, 600);
    // TH1D *hist = new TH1D("histogram", "Histogram of Window Probabilities", 100, 0, 14000);

    // for (double val : result)
    // {
    //     hist->Fill(val);
    // }

    // hist->Draw();

    // std::ofstream outputFile("totalPrime.txt"); // Open the file for writing
    // if (outputFile.is_open())
    // {
    //     for (double num : result)
    //     {
    //         outputFile << num << std::endl; // Write each element to the file
    //     }
    //     outputFile.close(); // Close the file
    //     std::cout << "Output file filled successfully." << std::endl;
    // }

    delete calibration_parameters; // Free memory for calibration_parameters
}

std::pair<double, double> ConfidenceLimit(std::vector<double> samples, double alpha)
{
    std::vector<double> data = samples;
    sort(data.begin(), data.end());

    int n = data.size();
    int left_idx = static_cast<int>((alpha) / 2 * n);
    int right_idx = static_cast<int>((1 - alpha / 2) * n);

    return std::make_pair(data[left_idx], data[right_idx]);
}

std::pair<double, double> MC_probability(const int N, const int window_size, const double alpha)
{
    double lifetime = 2196.981;
    std::exponential_distribution<double> distribution(1 / lifetime);

    // Generate N exponential distribution samples
    std::vector<double> samples;
    for (int i = 0; i < N; ++i)
    {
        double sample = distribution(gen);
        samples.push_back(sample);
    }

    // prob = (2*sigma/tau)*exp(-t/tau)
    // win prob = prod(prob) = (2*sigma/tau)^N * exp(-sum(t)/tau)
    // log(win prob) = N*log(2*sigma/tau) - sum(t)/tau
    // check = sum(t)

    std::vector<double> window_probabilities;
    for (size_t i = 0; i <= samples.size() - window_size; ++i)
    {
        double prob = 0;
        for (size_t j = i; j < i + window_size; ++j)
        {
            prob += samples[j];
        }
        window_probabilities.push_back(prob);
    }

    std::pair<double, double> limit = ConfidenceLimit(window_probabilities, alpha);

    // Create a canvas and a histogram
    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);
    TH1D *histogram = new TH1D("histogram", "Histogram of Window Probabilities", 100, 0, 140000);

    // Fill the histogram with the values from the vector
    for (size_t i = 0; i < window_probabilities.size(); ++i)
    {
        histogram->Fill(window_probabilities[i]); // i+0.5 to center the bin at i
    }

    // Draw the histogram on the canvas
    histogram->Draw("E");

    // Draw a vertical line corresponding to the upper limit
    TLine *line = new TLine(limit.first, 0, limit.first, histogram->GetMaximum());
    TLine *line2 = new TLine(limit.second, 0, limit.second, histogram->GetMaximum());
    line->SetLineColor(kRed); // Set line color to red
    line->Draw("same");       // Draw line on the same canvas

    line2->SetLineColor(kRed); // Set line color to red
    line2->Draw("same");       // Draw line on the same canvas

    return limit;
}

std::vector<double> removeIndices(const std::vector<double> &vec, const std::vector<int> &indices)
{
    std::vector<double> result;
    result.reserve(vec.size() - indices.size());

    // Creare un vettore con indici da rimuovere
    std::vector<bool> toRemove(vec.size(), false);
    for (int index : indices)
    {
        if (index >= 0 && index < vec.size())
        {
            toRemove[index] = true;
        }
    }

    // Copiare gli elementi non rimossi nel risultato
    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (!toRemove[i])
        {
            result.push_back(vec[i]);
        }
    }

    return result;
}