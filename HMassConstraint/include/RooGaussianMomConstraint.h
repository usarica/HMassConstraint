#ifndef ROOGAUSSIANMOMCONSTRAINT_H
#define ROOGAUSSIANMOMCONSTRAINT_H

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

#ifndef rmg_debug
#define rmg_debug 0
#endif


class RooGaussianMomConstraint : public RooAbsPdf {
public:

  enum CoordinateSystem{
    kXYZ,
    kRhoLambdaPhi // == pT, Lambda, Phi
  };
  enum CoordinatePrimes{
    prime_var1=2,
    prime_var2=3,
    prime_var3=5,
    prime_mean1=7,
    prime_mean2=11,
    prime_mean3=13
  };

  RooGaussianMomConstraint(){};
  RooGaussianMomConstraint(
    const char* name, const char* title,
    const RooArgList& variables_,
    const RooArgList& means_,
    const RooArgList& matrixElement_,
    RooGaussianMomConstraint::CoordinateSystem coordinates_=RooGaussianMomConstraint::kXYZ,
    Int_t fixCode_=1
    );
  RooGaussianMomConstraint(const RooGaussianMomConstraint& other, const char* name=0);
  inline virtual ~RooGaussianMomConstraint(){}

  virtual TObject* clone(const char* newname) const { return new RooGaussianMomConstraint(*this, newname); }

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
  virtual Double_t evaluate() const;

  virtual void fixVariable(Int_t code);

protected:

  RooListProxy variables;
  RooListProxy means;
  RooListProxy matrixElement;
  RooGaussianMomConstraint::CoordinateSystem coordinates;
  Int_t fixCode;

  virtual void setProxyList(const RooArgList& args, RooListProxy& target, Int_t checkDim=-1);

  virtual Double_t computeGaussian(const Int_t code) const;

  virtual Double_t computeCaseXYZ(const Int_t code) const;
  virtual Double_t computeCaseRhoLambdaPhi(const Int_t code) const;

};

#endif
