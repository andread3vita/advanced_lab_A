#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include <TColor.h>
#include <fstream>
#include <iostream>
#include <string>

void plot_v_gain(std::string file_name = "./results/txt_files/A_SX.txt")
{
    // Open the data file
    std::ifstream file(file_name.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file." << std::endl;
        return;
    }

    // Create TGraphErrors to hold the data
    TGraphErrors *graph = new TGraphErrors();
    graph->SetTitle(";V_{bias} [dV];Gain [mV]");

    // Read data from the file and fill the graph
    double x, y, ex, ey;
    int pointIndex = -1;
    while (file >> x >> y >> ex >> ey)
    {
        graph->SetPoint(pointIndex, x, y);
        graph->SetPointError(pointIndex, ex, ey);
        pointIndex++;
    }

    // Close the file
    file.close();

    // Create a canvas and draw the graph
    TCanvas *canvas = new TCanvas("canvas", "Graph", 800, 600);
    canvas->SetGrid();
    graph->SetMarkerSize(1.4);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);
    graph->GetXaxis()->SetLimits(275, 340);
    graph->GetXaxis()->SetRangeUser(275, 340);
    graph->GetYaxis()->SetRangeUser(0, 18);
    graph->GetXaxis()->SetLabelSize(0.045);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);

    graph->Draw("AP"); // "AP" draws the axis and the points

    // Perform linear fit and draw fit function
    TF1 *linearFit = new TF1("linearFit", "[0]*x + [1]");
    linearFit->SetParNames("a", "b");

    // Perform the fit
    graph->Fit(linearFit, "");

    // Get the fit parameters
    double a = linearFit->GetParameter(0);
    double b = linearFit->GetParameter(1);
    double vb = -b / a;

    // Draw the fit function on the canvas
    linearFit->SetLineColor(kRed);
    linearFit->Draw("same");

    // Draw a vertical line at vb
    TLine *line = new TLine(vb, 0, vb, 18);
    line->SetLineStyle(2);     // Imposta uno stile di linea tratteggiato
    line->SetLineWidth(2);     // Imposta la larghezza della linea
    line->SetLineColor(kBlue); // Imposta il colore della linea
    line->Draw("same");
    gStyle->SetOptFit(1);

    TPaveStats *pv = (TPaveStats *)graph->GetListOfFunctions()->FindObject("stats");
    pv->SetX1NDC(0.5);
    pv->SetX2NDC(0.9);
    pv->SetY1NDC(0.1);
    pv->SetY2NDC(0.3);

    TPaveText *text = new TPaveText(vb + 1, 16, vb + 20, 17);
    text->SetFillColor(0);  // Imposta il colore di riempimento a bianco
    text->SetLineColor(0);  // Imposta il colore del bordo a trasparente
    text->SetTextAlign(22); // Allinea il testo al centro orizzontale
    text->AddText("Breakdown Voltage");
    text->SetTextColor(kBlue);
    text->Draw("same");

    // Show the canvas
    canvas->Draw();
    canvas->SaveAs("./results/Fit_gain_A_SX.pdf");
}
