#include "../include/GrUtil.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

GrUtil::GrUtil()
{
}

GrUtil::~GrUtil()
{
}

std::string GrUtil::HLabel(const TH1D *h, std::string unita_x, std::string title)
{
  std::string titolo;

  Double_t xmin, xmax, xbin;
  xmin = h->GetXaxis()->GetXmin();
  xmax = h->GetXaxis()->GetXmax();
  xbin = h->GetXaxis()->GetNbins();

  Double_t    bin_width       = (xmax - xmin) / xbin;
  std::string temp            = std::to_string(bin_width);
  int         count           = 0;
  bool        first           = false;
  int         count_interno_x = 0;
  for (auto c : temp)
  {
    count_interno_x++;
    if (c != '0' && c != '.' && !first)
    {
      count = count_interno_x;
      first = true;
    }
  }
  temp.resize(count + 1);

  titolo = title + "/" + temp + " " + unita_x;
  return titolo;
}

std::string GrUtil::HLabel(const TH2D *h, std::string unita_x, std::string unita_y, std::string title)
{
  std::string titolo;

  Double_t xmin, xmax, xbin;
  xmin = h->GetXaxis()->GetXmin();
  xmax = h->GetXaxis()->GetXmax();
  xbin = h->GetXaxis()->GetNbins();

  Double_t ymin, ymax, ybin;
  ymin = h->GetYaxis()->GetXmin();
  ymax = h->GetYaxis()->GetXmax();
  ybin = h->GetYaxis()->GetNbins();

  Double_t bin_width_x = (xmax - xmin) / xbin;
  Double_t bin_width_y = (ymax - ymin) / ybin;

  std::string temp_x          = std::to_string(bin_width_x);
  int         count_x         = 0;
  bool        first_x         = false;
  int         count_interno_x = 0;
  for (auto c : temp_x)
  {
    count_interno_x++;
    if (c != '0' && c != '.' && !first_x)
    {
      count_x = count_interno_x;
      first_x = true;
    }
  }
  temp_x.resize(count_x + 1);

  std::string temp_y          = std::to_string(bin_width_y);
  int         count_y         = 0;
  bool        first_y         = false;
  int         count_interno_y = 0;
  for (auto c : temp_y)
  {
    count_interno_y++;
    if (c != '0' && c != '.' && !first_y)
    {
      count_y = count_interno_y;
      first_y = true;
    }
  }
  temp_y.resize(count_y + 1);

  titolo = title + "/(" + temp_x + " " + unita_x + "#times" + temp_y + " " + unita_y + ")";
  return titolo;
}

void GrUtil::SetHTextSize(TH1 *h)
{
  h->GetXaxis()->SetLabelSize(0.047);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.85);

  h->GetYaxis()->SetLabelSize(0.047);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.99);
}

void GrUtil::SetHTextSize(TH2D *h)
{
  h->GetXaxis()->SetLabelSize(0.047);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.85);

  h->GetYaxis()->SetLabelSize(0.047);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.03);

  h->GetZaxis()->SetLabelSize(0.047);
  h->GetZaxis()->SetTitleSize(0.05);
  h->GetZaxis()->SetTitleOffset(0.85);
}

void GrUtil::SetCountDigits(TH1D *h)
{
  if (h->GetMaximum() > 1000 && h->GetMaximum() < 10000) // evita che vada fuori dal canvas
  {
    h->GetYaxis()->SetMaxDigits(3);
  }
  else if (h->GetMaximum() > 10000 && h->GetMaximum() < 100000)
  {
    h->GetYaxis()->SetMaxDigits(4);
  }
  else if (h->GetMaximum() > 100000 && h->GetMaximum() < 1000000)
  {
    h->GetYaxis()->SetMaxDigits(2);
  }
}