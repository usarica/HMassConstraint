#include "../include/RooDiracDeltaFunction.h"

using namespace std;


RooDiracDeltaFunction::RooDiracDeltaFunction(
  const char* name, const char* title,
  RooAbsReal& var_,
  RooAbsReal& ref_
  ) : RooAbsPdf(name, title),

  var("var", "var", this, (RooAbsReal&)var_),
  ref("ref", "ref", this, (RooAbsReal&)ref_)
{}


RooDiracDeltaFunction::RooDiracDeltaFunction(const RooDiracDeltaFunction& other, const char* name) :
RooAbsPdf(other, name),

var("var", this, other.var),
ref("ref", this, other.ref)
{}

Bool_t RooDiracDeltaFunction::checkFundamentalType(const RooRealProxy& proxy) const{
  RooAbsArg* arg = proxy.absArg();
  return (dynamic_cast<RooRealVar*>(arg)!=0);
}

Int_t RooDiracDeltaFunction::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  Int_t code=1;
  if (checkFundamentalType(var)){ if (matchArgs(allVars, analVars, var)) code *= 2; }
  if (checkFundamentalType(ref)){ if (matchArgs(allVars, analVars, ref)) code *= 3; }
  if (code==6 && var.arg().hasRange(rangeName) && !ref.arg().hasRange(rangeName)) code=2;
  else if (code==6 && !var.arg().hasRange(rangeName) && ref.arg().hasRange(rangeName)) code=3;
  else if (code==6 && !var.arg().hasRange(rangeName) && !ref.arg().hasRange(rangeName)) code=1;
  if (code==1) code=0;
  //cout << "RooDiracDeltaFunction::getAnalyticalIntegral code = " << code << endl;
  return code;
}

Double_t RooDiracDeltaFunction::evaluate() const{
  const Double_t epsilon = 0;
  Double_t value = 1.;
  if (var!=ref) value = epsilon;
  //cout << "RooDiracDeltaFunction::evaluate = " << value << endl;
  return value;
}
Double_t RooDiracDeltaFunction::analyticalIntegral(Int_t code, const char* rangeName) const{
  const Double_t epsilon = 0;
  Double_t value = epsilon;
  if (code%2==0 && code%3!=0){
    if (var.min(rangeName)<ref && var.max(rangeName)>ref) value = 1.;
    else if (var.min(rangeName)<ref && var.max(rangeName)>=ref) value = 0.5;
    else if (var.min(rangeName)<=ref && var.max(rangeName)>ref) value = 0.5;
  }
  else if (code%3==0 && code%2!=0){
    if (ref.min(rangeName)<var && ref.max(rangeName)>var) value = 1.;
    else if (ref.min(rangeName)<var && ref.max(rangeName)>=var) value = 0.5;
    else if (ref.min(rangeName)<=var && ref.max(rangeName)>var) value = 0.5;
  }
  else if (var.min(rangeName)==ref.min(rangeName) && var.max(rangeName)==ref.max(rangeName)){
    Double_t minVal = max(var.min(rangeName), ref.min(rangeName));
    Double_t maxVal = min(var.max(rangeName), ref.max(rangeName));
    value = maxVal-minVal;
  }
  //cout << "RooDiracDeltaFunction::analyticalIntegral = " << value << " (code=" << code << ")" << endl;
  return value;
}


