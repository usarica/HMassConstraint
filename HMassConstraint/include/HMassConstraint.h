#ifndef HMassConstraint_H
#define HMassConstraint_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cmath>
#include <utility>
#include <algorithm>
#include "TLorentzVector.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"

#include <HMassConstraint/HMassConstraint/include/ScalarPdfFactory_ggH.h>
#include <HMassConstraint/HMassConstraint/include/TensorPdfFactory_HVV.h>
#include <HMassConstraint/HMassConstraint/include/RooGaussianMomConstraint.h>

#include <DataFormats/TrackReco/interface/TrackBase.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/JetReco/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include <RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h>


namespace HMCtoolkit{
  const Double_t pi_val = 3.14159265358979323846;// TMath::Pi();
  const Double_t piovertwo_val = 1.57079632679489661923;// pi_val/2.;
  const Int_t pdgTau=15;
  const Int_t pdgMu=13;
  const Int_t pdgEle=11;
  const Int_t pdgTauNu=16;
  const Int_t pdgMuNu=14;
  const Int_t pdgEleNu=12;
  const Int_t pdgUp=2;
  const Int_t pdgDown=1;
  const Int_t pdgCharm=4;
  const Int_t pdgStrange=3;
  const Int_t pdgTop=6;
  const Int_t pdgBottom=5;
  const Int_t pdgJet=0;
  const Int_t pdgUnknown=-99;

  template<typename varType> void deletePtr(varType* ptr){ if (ptr!=0) delete ptr; ptr=0; }
}


class HMassConstraint {
public:

  enum FitStrategy{
    FullCov_All_pTLambdaPhi,
    FullCov_All_pTLambda,
    FullCov_All_pTPhi,
    FullCov_All_LambdaPhi,
    FullCov_All_pT,
    FullCov_All_Lambda,
    FullCov_All_Phi,
    FullCov_NoFSR_pTLambdaPhi,
    FullCov_NoFSR_pTLambda,
    FullCov_NoFSR_pTPhi,
    FullCov_NoFSR_LambdaPhi,
    FullCov_NoFSR_pT,
    FullCov_NoFSR_Lambda,
    FullCov_NoFSR_Phi,
    FullCov_OnlyFSR_pTLambdaPhi,
    FullCov_OnlyFSR_pTLambda,
    FullCov_OnlyFSR_pTPhi,
    FullCov_OnlyFSR_LambdaPhi,
    FullCov_OnlyFSR_pT,
    FullCov_OnlyFSRR_Lambda,
    FullCov_OnlyFSR_Phi,
    CovDiagonals_All_pTLambdaPhi,
    CovDiagonals_All_pTLambda,
    CovDiagonals_All_pTPhi,
    CovDiagonals_All_LambdaPhi,
    CovDiagonals_All_pT,
    CovDiagonals_All_Lambda,
    CovDiagonals_All_Phi,
    CovDiagonals_NoFSR_pTLambdaPhi,
    CovDiagonals_NoFSR_pTLambda,
    CovDiagonals_NoFSR_pTPhi,
    CovDiagonals_NoFSR_LambdaPhi,
    CovDiagonals_NoFSR_pT,
    CovDiagonals_NoFSR_Lambda,
    CovDiagonals_NoFSR_Phi,
    CovDiagonals_OnlyFSR_pTLambdaPhi,
    CovDiagonals_OnlyFSR_pTLambda,
    CovDiagonals_OnlyFSR_pTPhi,
    CovDiagonals_OnlyFSR_LambdaPhi,
    CovDiagonals_OnlyFSR_pT,
    CovDiagonals_OnlyFSR_Lambda,
    CovDiagonals_OnlyFSR_Phi,

    nFitStrategies // This is actually not a strategy, just a number!
  };

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
    Double_t pTcut_muon_=5.,
    Double_t etacut_muon_=2.4,
    Double_t pTcut_electron_=7.,
    Double_t etacut_electron_=2.5,
    Double_t pTcut_jet_=30.,
    Double_t etacut_jet_=4.7,
    Double_t pTcut_fsr_=2.,
    Double_t etacut_fsr_=2.4
    );
  void setM1M2Cuts(
    Double_t m1cut_=40.,
    Double_t m2cut_=12.,
    Double_t mFFcut_=4.
    );
  //void reset(); // { if(fit_res!=0) delete fit_res; }

  // Make sure each strategy is implemented correctly. Affects the behavior of covariance matrix extractions in addDaughters.
  void setFitStrategy(HMassConstraint::FitStrategy fitStrategy_=HMassConstraint::FullCov_All_pT/*FullCov_All_pTLambdaPhi*/);
  HMassConstraint::FitStrategy getFitStrategy();
  // Add the fermion-FSR pairs, FSR-being per-lepton. fitRetry=true prevents clearing of the protected objects container and allows another fit through different strategies in case the initial fit fails.
  void addDaughters(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR, bool fitRetry=false); // To set the Lepton and photon arrays in a pair form. Pass null-pointer if the photon does not exist.
  // Do the fit, retry if unsuccessful.
  void fit();

  RooAbsPdf* getPDF();
  SpinPdfFactory* getPDFFactory();

  Double_t getMassError(Int_t iZ) const; // iZ==0 is m1, iZ==1 is m2, iZ==2 is m12.

protected:

  // Data members
  const Double_t sqrts;
  RooSpin::VdecayType Vdecay1;
  RooSpin::VdecayType Vdecay2;
  const Int_t X_spin;
  const Int_t intCodeStart;
  TString jecString;

  HMassConstraint::FitStrategy fitStrategy;
  HMassConstraint::FitStrategy fitStrategy_final;

  Double_t pTcut_muon;
  Double_t lambdacut_muon;
  Double_t pTcut_electron;
  Double_t lambdacut_electron;
  Double_t pTcut_jet;
  Double_t lambdacut_jet;
  Double_t pTcut_fsr;
  Double_t lambdacut_fsr;
  RooRealVar* m1cut;
  RooRealVar* m2cut;
  RooRealVar* mFFcut;

  RooRealVar* varZero;
  RooRealVar* varOne;

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

  RooRealVar* invcov_ferm[2][2][9];
  RooRealVar* invcov_fsr[2][2][9];

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

  RooFormulaVar* beta_Vdaughter[2];

  RooFormulaVar* massCuts;

  RooFormulaVar* m[3];
  RooFormulaVar* mAB[2];
  RooRealVar* h1;
  RooRealVar* h2;
  RooRealVar* hs;
  RooRealVar* Phi;
  RooRealVar* Phi1;
  RooRealVar* Y;
  RooSpin::modelMeasurables measurables;

  std::vector<pair<const reco::Candidate*, const pat::PFParticle*>> inputRaw_Fermion_FSR;

  SpinPdfFactory* pdfFactory;
  ScalarPdfFactory_ggH* hvvFactory;
  TensorPdfFactory_HVV* xvvFactory;
  RooSpin* spinPDF;
  RooGaussianMomConstraint* gausConstraintsPDF[2][2][2];
  RooGenericPdf* auxilliaryConstraintsPDF;
  RooProdPdf* constraintsPDF;
  RooProdPdf* PDF;

  RooFitResult* fitResult;
  TMatrixDSym fitCovMatrix;

  // Member functions
  virtual void constructVariables();
  virtual void destroyVariables();

  virtual void constructPdfFactory();
  virtual void destroyPdfFactory();

  virtual void constructConstraintPdfs();
  virtual void destroyConstraintPdfs();

  virtual void constructCompoundPdf();
  virtual void destroyCompoundPdf();

  void setWorkingFitStrategy(HMassConstraint::FitStrategy fitStrategy_);
  void testFitStrategy(Int_t& useFullCov, Int_t& FermFSRType, Int_t& fitpT, Int_t& fitlambda, Int_t& fitphi) const;
  void decrementStrategy(HMassConstraint::FitStrategy& strategy_);
  void incrementStrategy(HMassConstraint::FitStrategy& strategy_);

  // Overloaded functions to compute input block-diagonal C(pT, lambda, phi)
  void sortGetCovarianceMatrix(double (&momCov)[9], const reco::Candidate* particle);
  void getCovarianceMatrix(double (&momCov)[9], const reco::GsfElectron* particle);
  //void getCovarianceMatrix(double (&momCov)[9], const pat::Electron* particle);
  void getCovarianceMatrix(double (&momCov)[9], const pat::Muon* particle);
  void getCovarianceMatrix(double (&momCov)[9], const reco::PFCandidate* particle);
  void getCovarianceMatrix(double (&momCov)[9], const pat::Jet* particle);

  void invertOneDimensional(Int_t includeIndex, double (&momCov)[9]);
  void invertTwoDimensional(Int_t omitIndex, double (&momCov)[9]);
  void invertThreeDimensional(double (&momCov)[9]);
  void strategicInvertCovarianceMatrix(const Int_t useFullCov, Int_t fitpT, Int_t fitlambda, Int_t fitphi, double (&momCov)[9]);
  void setInverseCovarianceMatrix(Int_t iZ, Int_t iferm, Int_t fsrindex, Double_t momCov[9]);

  RooDataSet* getDataset() const;

  // Derivative functions to use in contracting final C(pT, lambda, phi) to sigma(m1), sigma(m2) or sigma(m12)
  Double_t d_Ek_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_pjk_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const;
  Double_t d_Ek_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_pjk_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const;
  Double_t d_Ek_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_pjk_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const;
};

#endif


