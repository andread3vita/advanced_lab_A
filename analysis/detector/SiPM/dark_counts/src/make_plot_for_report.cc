#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <iostream>
#include <string>

void make_plot_for_report()
{
    TFile *file_B = TFile::Open("./figures/dark_counts_B.root");
    TFile *file_C = TFile::Open("./figures/dark_counts_C.root");
    TGraphErrors *g_B = (TGraphErrors *)file_B->Get("g_B;2");
    TGraphErrors *g_C = (TGraphErrors *)file_C->Get("g_B");

    TCanvas *c = new TCanvas("c", "c", 1500, 600);
    c->Divide(2);

    // make B plot on the left
    c->cd(1);
    g_B->SetMarkerStyle(20);
    g_B->SetMarkerSize(1.1);
    g_B->SetTitle("");
    g_B->SetMarkerColor(kBlue);
    g_B->GetXaxis()->SetTitle("applied voltage [V]");
    g_B->GetYaxis()->SetTitle("dark counts / #mus");
    gPad->SetGrid();

    g_B->GetXaxis()->SetTitleSize(0.05);
    g_B->GetXaxis()->SetTitleOffset(0.97);
    g_B->GetXaxis()->SetLabelSize(0.05);
    g_B->GetYaxis()->SetTitleSize(0.05);
    g_B->GetYaxis()->SetTitleOffset(0.99);
    g_B->GetYaxis()->SetLabelSize(0.05);
    g_B->GetYaxis()->SetRangeUser(0.1, 2.8);

    g_B->Draw("AP");

    c->cd(2);
    g_C->SetMarkerStyle(20);
    g_C->SetMarkerSize(1.1);
    g_C->SetTitle("");
    g_C->SetMarkerColor(kBlue);
    g_C->GetXaxis()->SetTitle("applied voltage [V]");
    g_C->GetYaxis()->SetTitle("dark counts / #mus");
    gPad->SetGrid();

    g_C->GetXaxis()->SetTitleSize(0.05);
    g_C->GetXaxis()->SetTitleOffset(0.97);
    g_C->GetXaxis()->SetLabelSize(0.05);
    g_C->GetYaxis()->SetTitleSize(0.05);
    g_C->GetYaxis()->SetTitleOffset(0.99);
    g_C->GetYaxis()->SetRangeUser(0.1, 2.8);
    g_C->GetYaxis()->SetLabelSize(0.05);

    g_C->Draw("AP");

    c->SaveAs("./figures/dark_counts_bdx_cdx.pdf");

    return;
}
