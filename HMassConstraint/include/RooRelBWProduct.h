#ifndef ROORELBWPRODUCT_H
#define ROORELBWPRODUCT_H

#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpin.h>
#else
#include "RooSpin.h"
#endif


class RooRelBWProduct : public RooSpin {
public:


  RooRelBWProduct(){};
  RooRelBWProduct(
    const char* name, const char* title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );
  RooRelBWProduct(const RooRelBWProduct& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooRelBWProduct(*this, newname); }
  inline virtual ~RooRelBWProduct(){}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const { return 0; }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const { return 1.; }

protected:

  Double_t evaluate() const;

};

#endif
