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
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"

#include <HMassConstraint/HMassConstraint/include/ScalarPdfFactory_ggH.h>
#include <HMassConstraint/HMassConstraint/include/TensorPdfFactory_HVV.h>
#include <HMassConstraint/HMassConstraint/include/RooRelBWProduct.h>
#include <HMassConstraint/HMassConstraint/include/RooGaussianMomConstraint.h>
#include <HMassConstraint/HMassConstraint/include/RooDiracDeltaFunction.h>

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


#ifndef hmc_debug
#define hmc_debug 1
#endif


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

  enum FitMomentumStrategy{
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

    nFitMomentumStrategies // This is actually not a strategy, just a number!
  };
  enum FitVVStrategy{
    Fit_All_V1V2,
    Fit_All_V1, // = Fit_WorstTwo_V1
    Fit_All_V2, // = Fit_WorstTwo_V2

    nFitVVStrategies // This is actually not a strategy, just a number!
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
    Double_t m1lowcut_=40.,
    Double_t m2lowcut_=12.,
    Double_t m1highcut_=120.,
    Double_t m2highcut_=120.,
    Double_t mFFOScut_=4.,
    Double_t mFFSScut_=4.
    );

  // Set PDF to fast or thorough
  void setFastPDF(bool useFastPDF_=false);

  // Make sure each strategy is implemented correctly. Affects the behavior of covariance matrix extractions in addDaughters.
  void setFitMomentumStrategy(HMassConstraint::FitMomentumStrategy fitMomStrategy_=HMassConstraint::FullCov_All_pTLambdaPhi/*CovDiagonals_All_pT*/);
  void setFitVVStrategy(HMassConstraint::FitVVStrategy fitVVStrategy_=HMassConstraint::Fit_All_V1V2);
  HMassConstraint::FitMomentumStrategy getFitMomentumStrategy();

  // Do the fit for the fermion-FSR pairs, FSR-being per-fermion.
  void fitTo(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR);

  RooAbsPdf* getPDF();
  SpinPdfFactory* getPDFFactory();

  Double_t getRefittedMassError(Int_t imass) const; // imass==0 is m1, imass==1 is m2, imass==2 is m12.
  Double_t getRefittedMass(Int_t imass) const;
  TLorentzVector getRefittedMomentum(Int_t iZ, Int_t iferm, Int_t fsrindex) const;


protected:

  // Data members
  const Double_t sqrts;
  RooSpin::VdecayType Vdecay1;
  RooSpin::VdecayType Vdecay2;
  const Int_t X_spin;
  const Int_t intCodeStart;
  TString jecString;
  bool useFastPDF;

  HMassConstraint::FitMomentumStrategy fitMomStrategy;
  HMassConstraint::FitMomentumStrategy fitMomStrategy_final;

  HMassConstraint::FitVVStrategy fitVVStrategy; // There is no "final" counterpart to this one. Decrementation or incrementation is not allowed due to physics interest/reasons!

  Double_t pTcut_muon;
  Double_t lambdacut_muon;
  Double_t pTcut_electron;
  Double_t lambdacut_electron;
  Double_t pTcut_jet;
  Double_t lambdacut_jet;
  Double_t pTcut_fsr;
  Double_t lambdacut_fsr;
  RooRealVar* m1lowcut;
  RooRealVar* m2lowcut;
  RooRealVar* m1highcut;
  RooRealVar* m2highcut;
  RooRealVar* mFFOScut;
  RooRealVar* mFFSScut;

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

  RooRealVar* pTobs_ferm[2][2];
  RooRealVar* lambdaobs_ferm[2][2];
  RooRealVar* phiobs_ferm[2][2];
  RooRealVar* pTobs_fsr[2][2];
  RooRealVar* lambdaobs_fsr[2][2];
  RooRealVar* phiobs_fsr[2][2];

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
  RooFormulaVar* mAB[4];
  RooRealVar* h1;
  RooRealVar* h2;
  RooRealVar* hs;
  RooRealVar* Phi;
  RooRealVar* Phi1;
  RooRealVar* Y;
  RooSpin::modelMeasurables measurables;

  std::vector<pair<const reco::Candidate*, const pat::PFParticle*>> inputRaw_Fermion_FSR;

  // Delta functions of bar-obs
  RooDiracDeltaFunction* pTDeltaFcn_ferm[2][2];
  RooDiracDeltaFunction* lambdaDeltaFcn_ferm[2][2];
  RooDiracDeltaFunction* phiDeltaFcn_ferm[2][2];
  RooDiracDeltaFunction* pTDeltaFcn_fsr[2][2];
  RooDiracDeltaFunction* lambdaDeltaFcn_fsr[2][2];
  RooDiracDeltaFunction* phiDeltaFcn_fsr[2][2];
  RooProdPdf* DiracDeltaPDF;
  // J^CP PDFs
  SpinPdfFactory* pdfFactory;
  ScalarPdfFactory_ggH* hvvFactory;
  TensorPdfFactory_HVV* xvvFactory;
  RooSpin* spinPDF;
  // Fast propagator PDF
  RooRelBWProduct* bwProdPDF; // For fast PDF
  // PDF of constraintson {pT, lambda, phi}_i of fermion i
  RooGaussianMomConstraint* gausConstraintsPDF[2][2][2];
  // Other Heaviside functions
  RooGenericPdf* auxilliaryConstraintsPDF;
  // Products of constraints
  RooProdPdf* constraintsPDF;
  // Products constrantsPDF*auxilliaryConstraintsPDF*spinPDF or constrantsPDF*auxilliaryConstraintsPDF*bwProdPDF
  RooProdPdf* PDF;
  RooProdPdf* fastPDF;

  RooFitResult* fitResult;
  TMatrixDSym fitCovMatrix;

  // Member functions
  virtual void constructVariables();
  virtual void destroyVariables();

  virtual void constructDeltaFunctions();
  virtual void destroyDeltaFunctions();

  virtual void constructPdfFactory();
  virtual void destroyPdfFactory();

  virtual void constructConstraintPdfs();
  virtual void destroyConstraintPdfs();

  virtual void constructCompoundPdf();
  virtual void destroyCompoundPdf();

  // Add the fermion-FSR pairs, FSR-being per-fermion. fitRetry=true prevents clearing of the protected objects container and allows another fit through different strategies in case the initial fit fails.
  void addDaughters(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR, bool fitRetry=false); // To set the Lepton and photon arrays in a pair form. Pass null-pointer if the photon does not exist.
  // Do the fit, retry if unsuccessful.
  RooArgSet getDataVariables() const;
  RooDataSet* getDataset() const;
  void fit();

  // Momentum strategy functions
  void setWorkingFitMomentumStrategy(HMassConstraint::FitMomentumStrategy fitMomStrategy_);
  void testFitMomentumStrategy(Int_t& useFullCov, Int_t& FermFSRType, Int_t& fitpT, Int_t& fitlambda, Int_t& fitphi) const;
  void decrementMomentumStrategy(HMassConstraint::FitMomentumStrategy& strategy_);
  void incrementMomentumStrategy(HMassConstraint::FitMomentumStrategy& strategy_);
  // VV strategy functions
  void testFitVVStrategy(Int_t& fitV1, Int_t& fitV2) const;


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

  bool standardOrderedFinalCovarianceMatrix(const RooArgList& pars); // Re-order the covariance matrix from the fit, expand as necessary
  Int_t fitParameterCorrespondance(RooRealVar* par);


  // Derivative functions to use in contracting final C(pT, lambda, phi) to sigma(m1), sigma(m2) or sigma(m12)
  Double_t d_Ek_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_pjk_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const;
  Double_t d_Ek_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_pjk_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const;
  Double_t d_Ek_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_pjk_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const;

  Double_t d_m123_d_pTk(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_m123_d_lambdak(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const;
  Double_t d_m123_d_phik(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const;

};

#endif


