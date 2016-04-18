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
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include <HMassConstraint/interface/ScalarPdfFactory_ggH.h>
#include <HMassConstraint/interface/TensorPdfFactory_XVV.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>


class HMassConstraint {
public:

  // Member functions
  HMassConstraint(Double_t sqrts_, RooSpin::VdecayType Vdecay1_, RooSpin::VdecayType Vdecay2_); // See RooSpin.h for what enumerators are available.
  virtual ~HMassConstraint();

  void reset(); // { if(fit_res!=0) delete fit_res; }

  void addFermionsWithFSR(const std::vector<map<reco::Candidate*, pat::PFParticle*>>& FermionWithFSR); // To set the Lepton and photon arrays in a map form. Pass null-pointer if the photon does not exist.

  RooFitResult* fit(); // Does the fit and returns a pointer to the fitresult.

  // Data members
  RooSpin* pdf;
  SpinPdfFactory* pdfFactory;

private:

  RooRealVar* mass_ferm[2][2];
  RooRealVar* pT_ferm[2][2];
  RooRealVar* eta_ferm[2][2];
  RooRealVar* phi_ferm[2][2];
  RooRealVar* pTbar_ferm[2][2];
  RooRealVar* etabar_ferm[2][2];
  RooRealVar* phibar_ferm[2][2];

  RooRealVar* pT_fsr[2][2];
  RooRealVar* eta_fsr[2][2];
  RooRealVar* phi_fsr[2][2];
  RooRealVar* pTbar_fsr[2][2];
  RooRealVar* etabar_fsr[2][2];
  RooRealVar* phibar_fsr[2][2];

  RooFormulaVar* px_ferm[2][2];
  RooFormulaVar* py_ferm[2][2];
  RooFormulaVar* pz_ferm[2][2];
  RooFormulaVar* E_ferm[2][2];
  RooFormulaVar* px_fsr[2][2];
  RooFormulaVar* py_fsr[2][2];
  RooFormulaVar* pz_fsr[2][2];
  RooFormulaVar* E_fsr[2][2];

  RooFormulaVar* pTdiff_ferm[2][2];
  RooFormulaVar* etadiff_ferm[2][2];
  RooFormulaVar* phidiff_ferm[2][2];
  RooFormulaVar* pTdiff_fsr[2][2];
  RooFormulaVar* etadiff_fsr[2][2];
  RooFormulaVar* phidiff_fsr[2][2];

  RooFormulaVar* m1;
  RooFormulaVar* m2;
  RooFormulaVar* m12;
  RooRealVar* h1;
  RooRealVar* h2;
  RooRealVar* hs;
  RooRealVar* Phi;
  RooRealVar* Phi1;
  RooRealVar* Y;

  RooSpin::VdecayType Vdecay1;
  RooSpin::VdecayType Vdecay2;
  RooSpin::modelMeasurables measurables;
  const Double_t sqrts;
  const Int_t intCodeStart;

  std::vector<map<reco::Candidate*, TMatrixDSym>> fermCovMatrix; // Maximum 4 x (3x3)
  std::vector<map<pat::PFParticle*, TMatrixDSym>> gamCovMatrix; // Maximum 4 x (3x3)

  TMatrixDSym fitCovMatrix;
  RooGaussian* momCovConstraints;
  RooFitResult* fitRes;


  virtual void constructVariables();
  virtual void destroyVariables();

  virtual void constructPdfFactory();
  virtual void destroyPdfFactory();
};

#endif


