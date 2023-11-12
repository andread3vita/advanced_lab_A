// in questa classe ci saranno solamente funzioni static
// quando si va nel particolare per un'analisi specifica va creata una classe derivata

#ifndef GrUtil_h
#define GrUtil_h

#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

class GrUtil
{
  public:
  GrUtil();
  ~GrUtil();

  static std::string HLabel(const TH1D *h, std::string unita_x, std::string title = "counts");
  static std::string HLabel(const TH2D *h, std::string unita_x, std::string unita_y, std::string title = "counts");
  static void        SetHTextSize(TH1 *h);
  static void        SetHTextSize(TH2D *h);
  static void        SetCountDigits(TH1D *h);
  static void        SetCountDigits(TH2D *h);

  private:
};

#endif
