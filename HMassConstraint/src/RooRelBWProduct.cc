#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooRelBWProduct.h>
#else
#include "../include/RooRelBWProduct.h"
#endif


RooRelBWProduct::RooRelBWProduct(
  const char *name, const char *title,
  modelMeasurables _measurables,
  modelParameters _parameters,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2
  ) : RooSpin(
  name, title,
  _measurables,
  _parameters,
  _Vdecay1, _Vdecay2
  )
{}


RooRelBWProduct::RooRelBWProduct(
  const RooRelBWProduct& other, const char* name
  ) : RooSpin(other, name)
{}


Double_t RooRelBWProduct::evaluate() const{
  Double_t mV;
  getMVGamV(&mV);
  bool isZZ = (mV >= 90.);
  Double_t epsilon=1e-15;
  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) m1_=0;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;
  if (isZZ/* && Vdecay1==Vdecay2*/) {
    if ((m1_+m2_) > m12 || (fabs(m2_-mV)<fabs(m1_-mV) && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) || (m2_ <= 0. && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) || (m1_ <= 0. && Vdecay1!=RooSpin::kVdecayType_GammaOnshell)) return epsilon;
  }
  else if ((m1_+m2_) > m12 || ((m2_ <= 0. || m1_ <= 0.) && Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)) return epsilon;

  Double_t betaValSq = (1.-(pow(m1_-m2_, 2)/pow(m12, 2)))*(1.-(pow(m1_+m2_, 2)/pow(m12, 2)));
  if (betaValSq<0) return epsilon;
  Double_t betaVal = sqrt(betaValSq);

  Double_t term1Coeff = 1;
  Double_t term2Coeff = 1;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) term1Coeff = 2.*m1_;
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) term2Coeff = 2.*m2_;

  Double_t propV1Re=1, propV2Re=1;
  Double_t propV1Im=0, propV2Im=0;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) calculatePropagator(propV1Re, propV1Im, m1_, 1);
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) calculatePropagator(propV2Re, propV2Im, m2_, 1);

  //cout << "m1, m2, m12 = " << m1_ << " " << m2_ << " " << m12 << endl;

  Double_t valueRe = propV1Re*propV2Re-propV1Im*propV2Im;
  Double_t valueIm = propV1Re*propV2Im+propV1Im*propV2Re;
  Double_t value = term1Coeff*term2Coeff*betaVal*sqrt(pow(valueRe, 2)+pow(valueIm, 2));

  if (!(value==value)) cout << "Evaluate NaN=" << value << endl;
  if (value<=0){
    cout << "Evaluated value<=0: " << value << endl;
    value=epsilon;
  }
  return value;
}
