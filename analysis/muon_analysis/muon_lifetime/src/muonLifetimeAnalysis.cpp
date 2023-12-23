// Copyright 2023 Andrea De Vita - Nicolò Salimbeni

// muonLifetimeAnalysis.cpp
// Created on: Nov 27, 2023

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

#include "./../../../../include/AnUtil.h"
#include "./../../../../include/StateFile.h"
#include "./../../../../src/StateFile.cc"

/*

 - summaryPlot: it generates a comprehensive plot encompassing the fitted tau values, facilitating a side-by-side
                comparison with the theoretical value. Additionally, the weighted average fit value is calculated,
                incorporating its corresponding sigma;

 - residualsAnalysis: it fits the absolute value of the residuals w.r.t. lifetime between 10 ns to 2500 ns,
                      which is the region of negative muon lifetime regime;

 - biExponentialFit: it performs a fitting operation on the data histogram within the time range of 2000 to 15000 ns,
                     employing a negative exponential model. Subsequently, the obtained fit results serve as fixed
                     parameters for a second fitting process conducted within the 200 to 2000 ns range. This second fit
                     involves a combination of two exponentials, where one set of parameters is constrained based on the
                     outcomes of the previous fit, while the other set remains free.
                     The objective of this analysis is to investigate the initial segment of the histogram, aiming to
                     establish evidence for an additional contribution to the exponential decay, attributed to the decay
                     rate of negative muons absorbed in muonic states;

 - createFitFile: it creates a histogram from a data file, perform fits in different time ranges and
                  save results to a file.txt;

*/

//////////// MAIN FUNCTION DECLARATIONS ///////////////////
void summaryPlot(const char *datafile);
void residualsAnalysis(const char *datafile);
void biExponentialFit(const char *datafile);

//////////// MINOR FUNCTION DECLARATIONS ///////////////////
void createFitFile(const char *datafile);

//////////// FUNCTION DEFINITIONS ///////////////////
void summaryPlot(const char *resultFile)
{
    /////////////////
    // IMPORT DATA //
    /////////////////

    // Open the file
    std::ifstream infile(resultFile);

    // Define vectors for each column
    std::vector<string> range;
    std::vector<double> val, pos, err, prob;

    // Read the file line by line and fill the vectors
    double x, y, z, p;
    string r, columnNames;
    getline(infile, columnNames); // Read and discard column names
    while (infile >> r >> x >> y >> z >> p)
    {

        range.push_back(r);
        val.push_back(x);
        pos.push_back(y);
        err.push_back(z);
        prob.push_back(p);
    }

    // Close the file
    infile.close();

    // Print all values
    for (size_t i = 0; i < val.size(); ++i)
    {
        std::cout << "range: " << range[i] << "\t"
                  << "tau: " << val[i] << "\t"
                  << "sigma: " << err[i] << std::endl;
    }

    // Convert std::vectors to arrays
    double *x_arr = const_cast<double *>(val.data());
    double *y_arr = const_cast<double *>(pos.data());
    double *err_x_arr = const_cast<double *>(err.data());

    ////////////////
    // STATISTICS //
    ////////////////

    // Compute mean and sigma
    double mean = 0;
    double norm = 0;

    for (int i = 0; i < val.size(); ++i)
    {
        mean += x_arr[i] / pow(err_x_arr[i], 2);
        norm += pow(1.0 / err_x_arr[i], 2);
    }

    mean = mean / norm;
    double sigma = pow(norm, -0.5);

    //////////////
    // GRAPHICS //
    //////////////

    // Create TGraphErrors
    auto *c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->SetGrid();

    TGraphErrors *graph = new TGraphErrors(val.size(), x_arr, y_arr, err_x_arr);
    graph->SetMarkerStyle(kFullSquare);

    // Create a vertical line with a coloured region
    TLine *verticalLine = new TLine(mean, graph->GetYaxis()->GetBinLowEdge(graph->GetYaxis()->GetFirst()), mean,
                                    graph->GetYaxis()->GetBinUpEdge(graph->GetYaxis()->GetLast()));
    verticalLine->SetLineColor(42);
    verticalLine->SetLineWidth(2);

    // Create a vertical line with a coloured region
    double theory_val = 2196.9811;
    double err_theory = 0.0022;
    TLine *theory = new TLine(theory_val, graph->GetYaxis()->GetBinLowEdge(graph->GetYaxis()->GetFirst()), theory_val,
                              graph->GetYaxis()->GetBinUpEdge(graph->GetYaxis()->GetLast()));
    theory->SetLineColor(kRed);
    theory->SetLineWidth(2);
    theory->SetLineStyle(9);

    // Create the coloured region
    Double_t x1[5] = {mean - sigma, mean + sigma, mean + sigma, mean - sigma, mean - sigma};
    Double_t y1[5] = {graph->GetYaxis()->GetBinLowEdge(graph->GetYaxis()->GetFirst()),
                      graph->GetYaxis()->GetBinLowEdge(graph->GetYaxis()->GetFirst()),
                      graph->GetYaxis()->GetBinUpEdge(graph->GetYaxis()->GetLast()),
                      graph->GetYaxis()->GetBinUpEdge(graph->GetYaxis()->GetLast()),
                      graph->GetYaxis()->GetBinLowEdge(graph->GetYaxis()->GetFirst())};

    auto excl1 = new TGraph(5, x1, y1);
    excl1->SetFillColorAlpha(41, 0.5);
    excl1->SetFillStyle(1001);

    // Create multigraph
    TMultiGraph *multigraph = new TMultiGraph();
    multigraph->SetTitle("Fits performed at different intervals");

    multigraph->Add(excl1, "F");
    multigraph->Add(graph, "PE");

    multigraph->GetXaxis()->SetTitle("#tau (ns)");

    // Hide the tick labels on the TAxis
    multigraph->GetYaxis()->SetLabelSize(0);
    multigraph->GetYaxis()->SetTickLength(0);

    multigraph->Draw("A*");

    verticalLine->Draw("same");
    theory->Draw("same");

    // Draw text near each point
    for (int i = 0; i < val.size(); ++i)
    {
        // Adjust the position and content of the text as needed
        double textX = x_arr[i];
        double textY = y_arr[i] + 0.1; // Offset from the point

        TText *text = new TText(textX, textY, range[i].c_str());
        text->SetTextSize(0.035);
        text->Draw();
    }

    // Draw the Legend
    TLegend leg(.13, .67, .43, .87, "Fit results");
    leg.SetFillColor(0);
    leg.AddEntry(graph, "Fitted points");
    leg.AddEntry(theory, Form("Theoretical value:   %.3f #pm %.3f [ns]", theory_val, err_theory), "l");
    leg.AddEntry(verticalLine, Form("weighted average:   %.3f #pm %.3f [ns]", mean, sigma), "l");
    leg.AddEntry(excl1, "1#sigma region", "f");

    leg.DrawClone("Same");
}

void residualsAnalysis(const char *datafile)
{
    // Open data file
    std::ifstream infile(datafile);

    // Open Arietta calibration parameters file
    StateFile *calibration_parameters = new StateFile("../../../detector/arietta/results.txt");

    // Create a new canvas
    auto *c1 = new TCanvas();

    // Create a histogram for muon lifetime
    TH1D *histo = new TH1D("hist", "Muon lifetime", 200, 0, 16000);
    double x;

    double a = std::stod(calibration_parameters->ValueOf("m"));
    double b = std::stod(calibration_parameters->ValueOf("q"));

    // Fill the histogram with data from the file
    while (infile >> x)
    {
        histo->Fill(a * x + b);
    }
    infile.close();

    // Define a fit function
    TF1 *f = new TF1("f", "[0]+[1]*exp(-x/[2])", 2000, 15000);
    f->SetParNames("constant", "normalization", "decay rate [ns]");
    f->SetParLimits(0, 0, 200);
    f->SetParameter(2, 2200);

    // Perform the fit on the histogram
    histo->Fit(f, "RMN");

    // Create a graph for residuals
    TGraph *graph = new TGraph(30);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerSize(0.5);

    // Fill the graph with residuals
    for (int i = 1; i <= 30; ++i)
    {
        double xValue = histo->GetBinCenter(i);
        double yValueData = histo->GetBinContent(i);

        // Calculate the fitted value and residual
        double yValueFit = f->Eval(xValue);
        double residual = yValueData - yValueFit;

        // Add point to the graph
        graph->AddPoint(xValue, TMath::Abs(residual));
    }

    // Set titles for the graph axes
    graph->SetTitle("Residuals Analysis");
    graph->GetXaxis()->SetTitle("Lifetime [ns]");
    graph->GetYaxis()->SetTitle("Residuals");

    // Draw the graph
    graph->Draw("APL");

    // Define another fit function for a subrange
    TF1 *fitfunc = new TF1("fitfunc", "[0]+[1]*exp(-x/[2])", 20, 2500);
    fitfunc->SetParNames("constant", "normalization", "decay rate [ns]");
    fitfunc->SetParameter(2, 10000);
    fitfunc->SetParLimits(0, 0, 1000);

    // Perform the fit on the graph with the new fit function
    graph->Fit(fitfunc, "R+");

    // Draw default fit parameters legend
    gStyle->SetOptFit(1);

    delete calibration_parameters;
}

// Function for bi-exponential fit
void biExponentialFit(const char *datafile)
{
    // Open Arietta calibration parameters file
    StateFile *calibration_parameters = new StateFile("../../../detector/arietta/results.txt");

    // Open data file
    std::ifstream infile(datafile);

    // Initialize ROOT canvas
    auto *c1 = new TCanvas();

    // Create a histogram for muon lifetime
    TH1D *histo = new TH1D("hist", "Muon lifetime", 125, 0, 16000);

    // Read calibration parameters
    double x;
    double a = std::stod(calibration_parameters->ValueOf("m"));
    double b = std::stod(calibration_parameters->ValueOf("q"));

    // Fill the histogram using the calibrated values
    while (infile >> x)
    {
        histo->Fill(a * x + b);
    }
    infile.close();

    // Define and fit a single exponential function in the range [2000, 15000]
    TF1 *f = new TF1("f", "[0]+[1]*exp(-x/[2])", 2000, 15000);
    f->SetParNames("constant", "normalization", "decay rate [ns]");
    f->SetParLimits(0, 0, 200);
    f->SetParameter(2, 2200);
    histo->Fit(f, "RMN");

    // Define a combined function with fixed and free parameters for the second fit
    TF1 *combF = new TF1("combF", "[0]+[1]*exp(-x/[2])+[3]*exp(-x/[4])", 200, 2000);
    combF->FixParameter(2, f->GetParameter(2));
    combF->SetParameter(3, f->GetParameter(1));

    // Set limits and initial values for the parameters
    combF->SetParLimits(0, 0, 1000);
    combF->SetParLimits(1, 0, 10000);
    combF->SetParameter(4, 1000);

    // Fit the histogram with the combined function
    histo->Fit(combF, "RM");
}

// ///////////////// TRASH /////////////////
void createFitFile(const char *datafile)
{

    // Open the data file for reading
    std::ifstream infile(datafile);

    // Create a histogram for muon lifetime
    TH1D *histo = new TH1D("hist", "Muon lifetime", 75, 0, 16000);
    double x;

    // Open Arietta calibration parameters file
    StateFile *calibration_parameters = new StateFile("../../../detector/arietta/results.txt");

    // Retrieve calibration parameters
    double a = std::stod(calibration_parameters->ValueOf("m"));
    double b = std::stod(calibration_parameters->ValueOf("q"));

    // Fill the histogram with calibrated data
    while (infile >> x)
    {
        histo->Fill(a * x + b);
    }
    infile.close();

    // Open a file to store fit results
    std::ofstream resultFile("./../results/fitResult.txt");

    resultFile << "Range\tTau\tFit_number\tSigma\tQuantileProb" << std::endl;
    // Loop through different fit configurations
    for (int i = 1; i < 10; ++i)
    {
        // Create a fit function
        TF1 *f = new TF1("f", "[0]+[1]*exp(-x/[2])", 1000, 15000);
        f->SetParNames("constant", "normalization", "decay rate [ns]");

        // Set fit range based on iteration
        if (i < 5)
        {
            f->SetRange(15000 - 3700 * i, 15000);
        }
        else if (i == 5)
        {
            f->SetRange(200, 13000);
        }
        else if (i == 6)
        {
            f->SetRange(200, 2000);
        }
        else if (i == 7)
        {
            f->SetRange(5000, 15000);
        }
        else if (i == 8)
        {
            f->SetRange(3000, 9000);
        }
        else
        {
            f->SetRange(9000, 15000);
        }

        // Set parameter limits and initial values
        f->SetParLimits(0, 0, 200);
        f->SetParameter(2, 2200);

        // Fit the function to the histogram
        histo->Fit(f, "RMNQ");

        // Print fit results to the result file
        Double_t xmin, xmax;
        f->GetRange(xmin, xmax);

        resultFile << Form("%.0f-%.0f", xmin, xmax) << "\t";
        resultFile << f->GetParameter(2) << "\t";
        resultFile << i << "\t";
        resultFile << f->GetParError(2) << "\t";
        resultFile << f->GetProb() << std::endl;

        // Clean up memory by deleting the fit function object
        delete f;
    }
}
