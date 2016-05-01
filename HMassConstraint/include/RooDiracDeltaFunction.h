#ifndef ROODIRACDELTAFUNCTION_H
#define ROODIRACDELTAFUNCTION_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsCategory.h"
#include "Riostream.h" 
#include <cmath>
#include <vector>
#include "TMath.h"


class RooDiracDeltaFunction : public RooAbsPdf {
public:

  RooDiracDeltaFunction(){};
  RooDiracDeltaFunction(
    const char* name, const char* title,
    RooAbsReal& var_,
    RooAbsReal& ref_
    );
  RooDiracDeltaFunction(const RooDiracDeltaFunction& other, const char* name=0);
  inline virtual ~RooDiracDeltaFunction(){}

  virtual TObject* clone(const char* newname) const { return new RooDiracDeltaFunction(*this, newname); }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

  RooRealProxy var;
  RooRealProxy ref;

  Double_t evaluate() const;
  Bool_t checkFundamentalType(const RooRealProxy& proxy) const;

};

#endif
