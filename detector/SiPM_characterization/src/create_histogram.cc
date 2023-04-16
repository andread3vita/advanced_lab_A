// Copyright 2023 chatGPT + Nicolò Salimbeni
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>

void create_histogram() {
  // Crea un generatore di numeri casuali con una distribuzione gaussiana
  TRandom random;

  // Crea un istogramma con un binning di 1000 bin da 0 a 5
  TH1F* histogram =
      new TH1F("histogram", "4 gaussiani con media 1, 2, 3 e 4", 1000, 0, 5);

  // Aggiungi i dati casuali alle gaussiane
  for (int i = 0; i < 10000; i++) {
    histogram->Fill(random.Gaus(1, 0.1));
    histogram->Fill(random.Gaus(2, 0.1));
    histogram->Fill(random.Gaus(3, 0.1));
    histogram->Fill(random.Gaus(4, 0.1));
  }

  // Crea un file TFile per salvare l'istogramma
  TFile* file = new TFile("histogram.root", "RECREATE");

  // Scrivi l'istogramma nel file
  histogram->Write();

  // Chiudi il file
  file->Close();

  // Dealloca la memoria dell'istogramma e del file
  delete histogram;
  delete file;
}
