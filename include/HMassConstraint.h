#ifndef HMassConstraint_h
#define HMassConstraint_h

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cmath>
#include <map>
#include "TLorentzVector.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include <HMassConstraint/interface/ScalarPdfFactory_ggH.h>
#include <HMassConstraint/interface/TensorPdfFactory_XVV.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>


class HMassConstraint {
public:

  // Member functions
  HMassConstraint(
    const Double_t sqrts_, // 7, 8, 13, 14 etc. in units of TeV
    // See RooSpin.h for what enumerators are available:
    RooSpin::VdecayType Vdecay1_,
    RooSpin::VdecayType Vdecay2_,
    const Int_t X_spin_=0, // 0=spin-0, 2=spin-2
    // Do not touch the arguments below at the moment
    const Int_t intCodeStart_=RooSpin::prime_h1*RooSpin::prime_h2*RooSpin::prime_Phi*RooSpin::prime_hs*RooSpin::prime_Phi1
    );
  virtual ~HMassConstraint();

  void reset(); // { if(fit_res!=0) delete fit_res; }

  void addFermionsWithFSR(const std::vector<map<reco::Candidate*, pat::PFParticle*>>& FermionWithFSR); // To set the Lepton and photon arrays in a map form. Pass null-pointer if the photon does not exist.

  RooFitResult* fit(); // Does the fit and returns a pointer to the fitresult.

  // Data members
  RooSpin* spinPDF;
  SpinPdfFactory* pdfFactory;

private:

  // Data members
  const Double_t sqrts;
  RooSpin::VdecayType Vdecay1;
  RooSpin::VdecayType Vdecay2;
  const Int_t X_spin;
  const Int_t intCodeStart;

  RooRealVar* pT_ferm[2][2];
  RooRealVar* lambda_ferm[2][2]; // lambda = Pi/2 - theta
  RooRealVar* phi_ferm[2][2];
  RooRealVar* massbar_ferm[2][2];
  RooRealVar* pTbar_ferm[2][2];
  RooRealVar* lambdabar_ferm[2][2];
  RooRealVar* phibar_ferm[2][2];

  RooRealVar* pT_fsr[2][2];
  RooRealVar* lambda_fsr[2][2];
  RooRealVar* phi_fsr[2][2];
  RooRealVar* pTbar_fsr[2][2];
  RooRealVar* lambdabar_fsr[2][2];
  RooRealVar* phibar_fsr[2][2];

  RooFormulaVar* px_ferm[2][2];
  RooFormulaVar* py_ferm[2][2];
  RooFormulaVar* pz_ferm[2][2];
  RooFormulaVar* E_ferm[2][2];
  RooFormulaVar* px_fsr[2][2];
  RooFormulaVar* py_fsr[2][2];
  RooFormulaVar* pz_fsr[2][2];
  RooFormulaVar* E_fsr[2][2];
  RooFormulaVar* px_Hdaughter[2][2];
  RooFormulaVar* py_Hdaughter[2][2];
  RooFormulaVar* pz_Hdaughter[2][2];
  RooFormulaVar* E_Hdaughter[2][2];
  RooFormulaVar* m_Hdaughter[2][2]; // Use this for the extra phasespace multipliers

  RooFormulaVar* pTdiff_ferm[2][2];
  RooFormulaVar* lambdadiff_ferm[2][2];
  RooFormulaVar* phidiff_ferm[2][2];
  RooFormulaVar* pTdiff_fsr[2][2];
  RooFormulaVar* lambdadiff_fsr[2][2];
  RooFormulaVar* phidiff_fsr[2][2];

  RooFormulaVar* beta_Vdaughter[2];

  RooFormulaVar* m[3];
  RooRealVar* h1;
  RooRealVar* h2;
  RooRealVar* hs;
  RooRealVar* Phi;
  RooRealVar* Phi1;
  RooRealVar* Y;
  RooSpin::modelMeasurables measurables;

  ScalarPdfFactory_ggH* hvvFactory;
  TensorPdfFactory_XVV* xvvFactory;

  std::vector<map<reco::Candidate*, TMatrixDSym>> fermCovMatrix; // Maximum 4 x (3x3)
  std::vector<map<pat::PFParticle*, TMatrixDSym>> gamCovMatrix; // Maximum 4 x (3x3)

  TMatrixDSym fitCovMatrix;
  RooGaussian* momCovConstraints;
  RooFitResult* fitRes;


  // Member functions
  virtual void constructVariables();
  virtual void destroyVariables();

  virtual void constructPdfFactory();
  virtual void destroyPdfFactory();

};

#endif


