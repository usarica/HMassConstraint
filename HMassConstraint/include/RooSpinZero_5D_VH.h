#ifndef ROOSPINZERO_5D_VH
#define ROOSPINZERO_5D_VH

#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpinZero.h>
#else
#include "RooSpinZero.h"
#endif


class RooSpinZero_5D_VH : public RooSpinZero {
public:

  RooSpinZero_5D_VH(){}
  RooSpinZero_5D_VH(
    const char *name, const char *title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    modelCouplings _couplings,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );

  RooSpinZero_5D_VH(const RooSpinZero_5D_VH& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSpinZero_5D_VH(*this, newname); }
  inline virtual ~RooSpinZero_5D_VH(){}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

  Double_t evaluate() const;

  void evaluatePolarizationTerms(Double_t& A00term, Double_t& Appterm, Double_t& Ammterm, Double_t& A00ppterm, Double_t& A00mmterm, Double_t& Appmmterm, const Int_t code, bool isGammaV1=false, bool isGammaV2=false) const;

};
 
#endif
