#ifndef HMassConstraint_H
#define HMassConstraint_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cmath>
#include <pair>
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
#include <RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h>


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

  void setJECUserFloatString(TString jecString_="jec_unc"){ jecString=jecString_; }
  void setPtEtaCuts(
    Double_t pTcut_muon_=0.,
    Double_t etacut_muon_=-1.,
    Double_t pTcut_electron_=0.,
    Double_t etacut_electron_=-1.,
    Double_t pTcut_jet_=0.,
    Double_t etacut_jet_=-1.,
    Double_t pTcut_fsr_=0.,
    Double_t etacut_fsr_=-1.
    );
  void reset(); // { if(fit_res!=0) delete fit_res; }
  void addDaughters(const std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR); // To set the Lepton and photon arrays in a pair form. Pass null-pointer if the photon does not exist.
  void fit(); // Does the fit

  RooAbsPdf* getPDF(){ return PDF; }
  RooAbsPdf* getPDFFactory(){ return pdffactory; }


protected:

  enum{
    pdgTau=15,
    pdgMu=13,
    pdgEle=11,
    pdgTauNu=16,
    pdgMuNu=14,
    pdgEleNu=12,
    pdgUp=2,
    pdgDown=1,
    pdgCharm=4,
    pdgStrange=3,
    pdgTop=6,
    pdgBottom=5,
    pdgJet=0,
    pdgUnknown=-99
  }

  // Data members
  const Double_t sqrts;
  RooSpin::VdecayType Vdecay1;
  RooSpin::VdecayType Vdecay2;
  const Int_t X_spin;
  const Int_t intCodeStart;
  TString jecString;

  Double_t pTcut_muon;
  Double_t lambdacut_muon;
  Double_t pTcut_electron;
  Double_t lambdacut_electron;
  Double_t pTcut_jet;
  Double_t lambdacut_jet;

  RooRealVar* pT_ferm[2][2];
  RooRealVar* lambda_ferm[2][2]; // lambda = Pi/2 - theta
  RooRealVar* phi_ferm[2][2];

  Int_t pdgid_ferm[2][2];
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

  std::vector<pair<TMatrixDSym, TMatrixDSym>> inputCovMatrix; // Maximum 4x2 x (3x3)
  TMatrixDSym prefitInverseCovMatrix;
  TMatrixDSym fitCovMatrix;

  SpinPdfFactory* pdfFactory;
  ScalarPdfFactory_ggH* hvvFactory;
  TensorPdfFactory_XVV* xvvFactory;
  RooSpin* spinPDF;
  RooGaussian* constraintsPDF;
  RooProdPdf* PDF;

  RooFitResult* fitResult;


  // Member functions
  virtual void constructVariables();
  virtual void destroyVariables();

  virtual void constructPdfFactory();
  virtual void destroyPdfFactory();

  virtual void constructPdfConstraint();
  virtual void destroyPdfConstraint();

  virtual void constructCompoundPdf();
  virtual void destroyCompoundPdf();

  TMatrixDSym getCovarianceMatrix(const reco::Candidate* particle);
  TMatrixDSym getCovarianceMatrix(const reco::GsfElectron* particle);
  TMatrixDSym getCovarianceMatrix(const pat::Muon* particle);
  TMatrixDSym getCovarianceMatrix(const reco::PFCandidate* particle);
};

#endif


