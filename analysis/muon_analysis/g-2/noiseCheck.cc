// Copyright 2024 Andrea De Vita - Nicolò Salimbeni

// noiseCheck.cc
// Created on: Feb 3, 2024

// cpp libraries
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::string;

// root libraries
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"

#include "./../../../include/AnUtil.h"
#include "./../../../include/StateFile.h"
#include "./../../../src/StateFile.cc"

/*


*/

//////////// MAIN FUNCTION DECLARATIONS ///////////////////

//////////// MINOR FUNCTION DECLARATIONS ///////////////////

//////////// FUNCTION DEFINITIONS ///////////////////

void test(const char *datafile, const int window_size, const double interval)
{
    // Open data file
    std::ifstream infile(datafile);

    // Open Arietta calibration parameters file
    StateFile *calibration_parameters = new StateFile("../../detector/arietta/results.txt");

    double a = std::stod(calibration_parameters->ValueOf("m"));
    double b = std::stod(calibration_parameters->ValueOf("q"));

    std::vector<double> converted_values; // Vector to store converted values

    double lifetime = 2196.981;
    double prob = 1 - exp(-interval / lifetime);
    double expected_points = prob * window_size;

    // Read rows from the infile, convert, and store in the vector
    double x;
    while (infile >> x)
    {
        double converted_value = a * x + b;
        converted_values.push_back(converted_value);
    }

    // Print sliding window of 5 elements
    for (size_t i = 0; i <= converted_values.size() - window_size; ++i)
    {
        int count = 0;
        for (size_t j = i; j < i + window_size; ++j)
        {
            if (converted_values[j] < interval)
            {
                count += 1;
            }
        }
        if (count > int(expected_points))
        {
            std::cout << "Window " << i + 1 << ":" << count << std::endl;
        }
    }

    std::cout << "\nExpected points: " << int(expected_points) << std::endl;
    delete calibration_parameters; // Free memory for calibration_parameters
}