#include "../include/RooDiracDeltaFunction.h"

using namespace std;


RooDiracDeltaFunction::RooDiracDeltaFunction(
  const char* name, const char* title,
  const RooAbsReal& var_,
  const RooAbsReal& ref_
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
  if (checkFundamentalType(var)) return 2;
  else return 0; 
}

Double_t RooDiracDeltaFunction::evaluate() const{
  const Double_t epsilon = 0;
  if (var!=ref) return epsilon;
  else return 1.;
}
Double_t RooDiracDeltaFunction::analyticalIntegral(Int_t code, const char* rangeName) const{
  const Double_t epsilon = 0;
  if (code!=2) return epsilon;
  if (var.min(rangeName)<ref && var.max(rangeName)>ref) return 1.;
  else if (var.min(rangeName)<ref && var.max(rangeName)>=ref) return 0.5;
  else if (var.min(rangeName)<=ref && var.max(rangeName)>ref) return 0.5;
  return epsilon;
}


