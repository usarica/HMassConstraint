#include <HMassConstraint/HMassConstraint/include/HMassConstraint.h>
#include <cassert>
#include "RooFit.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "RooMinuit.h"
#include "RooNameReg.h"
#include "RooCmdConfig.h"
#include "RooGlobalFunc.h"
#include "RooAddition.h"
#include "RooNumIntConfig.h"
#include "RooProjectedPdf.h"
#include "RooInt.h"
#include "RooCustomizer.h"
#include "RooConstraintSum.h"
#include "RooCachedReal.h"
#include "RooXYChi2Var.h"
#include "RooChi2Var.h"
#include "RooRealIntegral.h"
#include "Math/Minimizer.h"


using namespace std;
using namespace reco;
using namespace pat;
using namespace HMCtoolkit;


HMassConstraint::HMassConstraint(
  const Double_t sqrts_,
  RooSpin::VdecayType Vdecay1_,
  RooSpin::VdecayType Vdecay2_,
  const Int_t X_spin_,
  const Int_t intCodeStart_
  ) :
  sqrts(sqrts_*1000.),
  Vdecay1(Vdecay1_),
  Vdecay2(Vdecay2_),
  X_spin(X_spin_),
  intCodeStart(intCodeStart_)
{
  if (X_spin!=0 && X_spin!=2){
    cerr << "HMassConstraint::HMassConstraint: Only spin-0 or -2 are supported!" << endl;
    assert(0);
  }

  setFastPDF(); // By default, this is set to false. If true, it will attempt to use the product of two BWs.

  setFitMomentumStrategy(); // Set default momentum strategy
  setFitVVStrategy(); // Set default VV strategy
  setJECUserFloatString(); // Set JEC uncertainty default string

  setPtEtaCuts(); // Set default cuts on pT, eta and phi of leptons, jets and FSR. Note the the cut targets are numbers!
  constructVariables();
  setM1M2Cuts(); // Set default minimum cuts on m1, m2, mA, mB. Note that the cut targets are RooRealVars, so constructVariables needs to be called first!
  constructDeltaFunctions();
  constructPdfFactory();
  constructConstraintPdfs();
  constructCompoundPdf();
}

HMassConstraint::~HMassConstraint(){
  destroyCompoundPdf();
  destroyConstraintPdfs();
  destroyPdfFactory();
  destroyDeltaFunctions();
  destroyVariables();
}

RooAbsPdf* HMassConstraint::getPDF(){ return PDF; }
SpinPdfFactory* HMassConstraint::getPDFFactory(){ return pdfFactory; }

void HMassConstraint::constructVariables(){
  fitResult=0;
  varZero = new RooRealVar("varZero", "", 0.);
  varOne = new RooRealVar("varOne", "", 1.);

  m1cut = new RooRealVar("m1cut", "", 0., 0., sqrts);
  m2cut = new RooRealVar("m2cut", "", 0., 0., sqrts);
  mFFcut = new RooRealVar("mFFcut", "", 0., 0., sqrts);

  RooArgList m12_args;
  for (int iZ=0; iZ<2; iZ++){
    RooArgList mHdaughter_args;
    for (int iferm=0; iferm<2; iferm++){
      pT_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambda_ferm[iZ][iferm] = new RooRealVar(Form("lambdaRefit_Z%iFermion%i", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phi_ferm[iZ][iferm] = new RooRealVar(Form("phiRefit_Z%iFermion%i", iZ+1, iferm+1), "", -pi_val, pi_val);
      massbar_ferm[iZ][iferm] = new RooRealVar(Form("mass_Z%iFermion%i", iZ+1, iferm+1), "", 0., -sqrts, sqrts);
      pTbar_ferm[iZ][iferm] = new RooRealVar(Form("pTInit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdabar_ferm[iZ][iferm] = new RooRealVar(Form("lambdaInit_Z%iFermion%i", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phibar_ferm[iZ][iferm] = new RooRealVar(Form("phiInit_Z%iFermion%i", iZ+1, iferm+1), "", -pi_val, pi_val);

      pT_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambda_fsr[iZ][iferm] = new RooRealVar(Form("lambdaRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phi_fsr[iZ][iferm] = new RooRealVar(Form("phiRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -pi_val, pi_val);
      pTbar_fsr[iZ][iferm] = new RooRealVar(Form("pTInit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdabar_fsr[iZ][iferm] = new RooRealVar(Form("lambdaInit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phibar_fsr[iZ][iferm] = new RooRealVar(Form("phiInit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -pi_val, pi_val);

      pTobs_ferm[iZ][iferm] = new RooRealVar(Form("pTObs_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdaobs_ferm[iZ][iferm] = new RooRealVar(Form("lambdaObs_Z%iFermion%i", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phiobs_ferm[iZ][iferm] = new RooRealVar(Form("phiObs_Z%iFermion%i", iZ+1, iferm+1), "", -pi_val, pi_val);
      pTobs_fsr[iZ][iferm] = new RooRealVar(Form("pTObs_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdaobs_fsr[iZ][iferm] = new RooRealVar(Form("lambdaObs_Z%iFermion%iFSR", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phiobs_fsr[iZ][iferm] = new RooRealVar(Form("phiObs_Z%iFermion%iFSR", iZ+1, iferm+1), "", -pi_val, pi_val);

      massbar_ferm[iZ][iferm]->removeMin();
      massbar_ferm[iZ][iferm]->removeMax();

      // Covariance elements
      for (int ix=0; ix<3; ix++){
        TString title_ix="";
        if (ix==0) title_ix="pT";
        else if (ix==1) title_ix="lambda";
        else title_ix="phi";
        for (int iy=0; iy<3; iy++){
          TString title_iy="";
          if (iy==0) title_iy="pT";
          else if (iy==1) title_iy="lambda";
          else title_iy="phi";
          Double_t minval=-1e9, maxval=1e9;
          if (ix==iy) minval=0;
          invcov_ferm[iZ][iferm][3*ix+iy] = new RooRealVar(Form("%s_vs_%s_Z%iFermion%i", title_ix.Data(), title_iy.Data(), iZ+1, iferm+1), "", 0., minval, maxval);
          if (ix!=iy) invcov_ferm[iZ][iferm][3*ix+iy]->removeMin();
          invcov_ferm[iZ][iferm][3*ix+iy]->removeMax();
          invcov_fsr[iZ][iferm][3*ix+iy] = new RooRealVar(Form("%s_vs_%s_Z%iFermion%iFSR", title_ix.Data(), title_iy.Data(), iZ+1, iferm+1), "", 0., minval, maxval);
          if (ix!=iy) invcov_fsr[iZ][iferm][3*ix+iy]->removeMin();
          invcov_fsr[iZ][iferm][3*ix+iy]->removeMax();
        }
      }

      E_ferm[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%i", iZ+1, iferm+1), "sqrt( pow(@0,2)*TMath::Sign(1.,pow(@0,2)) + pow(@1/cos(@2),2) )", RooArgList(*(massbar_ferm[iZ][iferm]), *(pT_ferm[iZ][iferm]), *(lambda_ferm[iZ][iferm])));
      px_ferm[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      py_ferm[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      pz_ferm[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*tan(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(lambda_ferm[iZ][iferm])));
      E_fsr[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0/cos(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(lambda_fsr[iZ][iferm])));
      px_fsr[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      py_fsr[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      pz_fsr[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*tan(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(lambda_fsr[iZ][iferm])));

      E_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(E_ferm[iZ][iferm]), *(E_fsr[iZ][iferm])));
      px_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(px_ferm[iZ][iferm]), *(px_fsr[iZ][iferm])));
      py_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(py_ferm[iZ][iferm]), *(py_fsr[iZ][iferm])));
      pz_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(pz_ferm[iZ][iferm]), *(pz_fsr[iZ][iferm])));
      // There is no Sign needed for m_Hdaughter since it is only used in beta_Vdaughter, which does not care about the sign.
      m_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("m_Z%iDau%iRefit", iZ+1, iferm+1), "sqrt( abs(pow(@0,2)-pow(@1,2)-pow(@2,2)-pow(@3,2)) )", RooArgList(*(E_Hdaughter[iZ][iferm]), *(px_Hdaughter[iZ][iferm]), *(py_Hdaughter[iZ][iferm]), *(pz_Hdaughter[iZ][iferm])));

      mHdaughter_args.add(*(E_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(px_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(py_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(pz_Hdaughter[iZ][iferm]));
    }

    // Add mHdaughter arguments into m12 arguments
    m12_args.add(mHdaughter_args);

    // Construct m1/m2
    m[iZ] = new RooFormulaVar(Form("m%iRefit", iZ+1), "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mHdaughter_args);

    // This beta should multiply the spinPDF bc. having massive fermions have additional scale factors. These factors are even more relevant when FSR is present!
    beta_Vdaughter[iZ] = new RooFormulaVar(Form("betaV%iRefit", iZ+1), "(@0>0. ? sqrt( TMath::Max(1e-15, ( 1.-pow((@1+@2)/@0,2) )*( 1.-pow((@1-@2)/@0,2) ) ) ) : 1.)", RooArgList(*(m[iZ]), *(m_Hdaughter[iZ][0]), *(m_Hdaughter[iZ][1])));
  }

  // Construct m12
  m[2] = new RooFormulaVar("m12Refit", "sqrt( pow(@0+@4+@8+@12,2)-pow(@1+@5+@9+@13,2)-pow(@2+@6+@10+@14,2)-pow(@3+@7+@11+@15,2) )", m12_args);

  // Construct mA, mB
  RooArgList mA_args;
  RooArgList mB_args;
  for (int iZ=0; iZ<2; iZ++){
    int iferm = 0;
    if (iZ==0){
      mA_args.add(*(E_Hdaughter[iZ][iferm]));
      mA_args.add(*(px_Hdaughter[iZ][iferm]));
      mA_args.add(*(py_Hdaughter[iZ][iferm]));
      mA_args.add(*(pz_Hdaughter[iZ][iferm]));
    }
    else{
      mB_args.add(*(E_Hdaughter[iZ][iferm]));
      mB_args.add(*(px_Hdaughter[iZ][iferm]));
      mB_args.add(*(py_Hdaughter[iZ][iferm]));
      mB_args.add(*(pz_Hdaughter[iZ][iferm]));
    }
    iferm = 1;
    if (iZ==0){
      mA_args.add(*(E_Hdaughter[1-iZ][iferm]));
      mA_args.add(*(px_Hdaughter[1-iZ][iferm]));
      mA_args.add(*(py_Hdaughter[1-iZ][iferm]));
      mA_args.add(*(pz_Hdaughter[1-iZ][iferm]));
    }
    else{
      mB_args.add(*(E_Hdaughter[1-iZ][iferm]));
      mB_args.add(*(px_Hdaughter[1-iZ][iferm]));
      mB_args.add(*(py_Hdaughter[1-iZ][iferm]));
      mB_args.add(*(pz_Hdaughter[1-iZ][iferm]));
    }
  }
  mAB[0] = new RooFormulaVar("mARefit", "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mA_args);
  mAB[1] = new RooFormulaVar("mBRefit", "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mB_args);

  massCuts = new RooFormulaVar("mABCutParameterization", "( (@0>@4 && @1>@4 && @2>@5 && @3>@6) ? 1. : 1.e-15)", RooArgList(*(mAB[0]), *(mAB[1]), *(m[0]), *(m[1]), *mFFcut, *m1cut, *m2cut));

  // Variables integrated over
  // Should be re-written in terms of pT, lambda and phi at some point...
  hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{1}", -1, 1);
  h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{2}", -1, 1);
  Phi = new RooRealVar("Genhelphi", "#Phi", -pi_val, pi_val);
  Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{1}", -pi_val, pi_val);
  Y = new RooRealVar("GenY", "Y", 0);

  // Initialize the meaurables: Set those always integrated over to 0
  measurables.m1 = m[0];
  measurables.m2 = m[1];
  measurables.m12 = m[2];
  if (intCodeStart%RooSpin::prime_h1 != 0) measurables.h1 = h1;
  else measurables.h1 = 0;
  if (intCodeStart%RooSpin::prime_h2 != 0) measurables.h2 = h2;
  else measurables.h2 = 0;
  if (intCodeStart%RooSpin::prime_Phi != 0) measurables.Phi = Phi;
  else measurables.Phi = 0;
  if (intCodeStart%RooSpin::prime_hs != 0) measurables.hs = hs;
  else measurables.hs = 0;
  if (intCodeStart%RooSpin::prime_Phi1 != 0) measurables.Phi1 = Phi1;
  else measurables.Phi1 = 0;
  measurables.Y = Y;
}
void HMassConstraint::constructDeltaFunctions(){
  RooArgList diracdelta_args;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      pTDeltaFcn_ferm[iZ][iferm] = new RooDiracDeltaFunction(Form("%s_DiracDeltaFcn", pTobs_ferm[iZ][iferm]->GetName()), Form("%s_DiracDeltaFcn", pTobs_ferm[iZ][iferm]->GetName()), *(pTobs_ferm[iZ][iferm]), *(pTbar_ferm[iZ][iferm]));
      lambdaDeltaFcn_ferm[iZ][iferm] = new RooDiracDeltaFunction(Form("%s_DiracDeltaFcn", lambdaobs_ferm[iZ][iferm]->GetName()), Form("%s_DiracDeltaFcn", lambdaobs_ferm[iZ][iferm]->GetName()), *(lambdaobs_ferm[iZ][iferm]), *(lambdabar_ferm[iZ][iferm]));
      phiDeltaFcn_ferm[iZ][iferm] = new RooDiracDeltaFunction(Form("%s_DiracDeltaFcn", phiobs_ferm[iZ][iferm]->GetName()), Form("%s_DiracDeltaFcn", phiobs_ferm[iZ][iferm]->GetName()), *(phiobs_ferm[iZ][iferm]), *(phibar_ferm[iZ][iferm]));
      pTDeltaFcn_fsr[iZ][iferm] = new RooDiracDeltaFunction(Form("%s_DiracDeltaFcn", pTobs_fsr[iZ][iferm]->GetName()), Form("%s_DiracDeltaFcn", pTobs_fsr[iZ][iferm]->GetName()), *(pTobs_fsr[iZ][iferm]), *(pTbar_fsr[iZ][iferm]));
      lambdaDeltaFcn_fsr[iZ][iferm] = new RooDiracDeltaFunction(Form("%s_DiracDeltaFcn", lambdaobs_fsr[iZ][iferm]->GetName()), Form("%s_DiracDeltaFcn", lambdaobs_fsr[iZ][iferm]->GetName()), *(lambdaobs_fsr[iZ][iferm]), *(lambdabar_fsr[iZ][iferm]));
      phiDeltaFcn_fsr[iZ][iferm] = new RooDiracDeltaFunction(Form("%s_DiracDeltaFcn", phiobs_fsr[iZ][iferm]->GetName()), Form("%s_DiracDeltaFcn", phiobs_fsr[iZ][iferm]->GetName()), *(phiobs_fsr[iZ][iferm]), *(phibar_fsr[iZ][iferm]));

      diracdelta_args.add(*(pTDeltaFcn_ferm[iZ][iferm]));
      diracdelta_args.add(*(lambdaDeltaFcn_ferm[iZ][iferm]));
      diracdelta_args.add(*(phiDeltaFcn_ferm[iZ][iferm]));
      diracdelta_args.add(*(pTDeltaFcn_fsr[iZ][iferm]));
      diracdelta_args.add(*(lambdaDeltaFcn_fsr[iZ][iferm]));
      diracdelta_args.add(*(phiDeltaFcn_fsr[iZ][iferm]));
    }
  }
  DiracDeltaPDF = new RooProdPdf("DiracDeltaPDF", "DiracDeltaPDF", diracdelta_args);
}
void HMassConstraint::constructPdfFactory(){
  hvvFactory=0;
  xvvFactory=0;
  if (X_spin==0){
    hvvFactory = new ScalarPdfFactory_ggH(measurables, false, Vdecay1, Vdecay2, true); // First false for acceptance, then true for always-on-shell H
    hvvFactory->makeParamsConst(false); // So that we can play with couplings
    spinPDF = hvvFactory->getPDF();
    pdfFactory = hvvFactory;
  }
  else if (X_spin==2){
    xvvFactory = new TensorPdfFactory_HVV(measurables, Vdecay1, Vdecay2, true); // true for always-on-shell X
    xvvFactory->makeParamsConst(false); // So that we can play with couplings
    spinPDF = xvvFactory->getPDF();
    pdfFactory = xvvFactory;
  }
  spinPDF->alwaysIntegrate(intCodeStart);

  bwProdPDF = new RooRelBWProduct(
    "RelBWProductPDF", "RelBWProductPDF",
    pdfFactory->measurables,
    pdfFactory->parameters,
    Vdecay1, Vdecay2
    );

#if hmc_debug==1
  for (int iZ=0; iZ<2; iZ++) simpleBWPDF[iZ] = new RooGenericPdf(Form("Z%i_simpleBWPDF", iZ+1), "2*@0/( (@0**2-@1**2)**2 + (@1*@2)**2 )", RooArgList(*(m[iZ]), *((pdfFactory->parameters).mZ), *((pdfFactory->parameters).gamZ)));
#endif
}
void HMassConstraint::constructConstraintPdfs(){
  RooArgList constraints;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      RooArgList vars_ferm, means_ferm, me_ferm;
      vars_ferm.add(*(pT_ferm[iZ][iferm]));
      vars_ferm.add(*(lambda_ferm[iZ][iferm]));
      vars_ferm.add(*(phi_ferm[iZ][iferm]));
      means_ferm.add(*(pTbar_ferm[iZ][iferm]));
      means_ferm.add(*(lambdabar_ferm[iZ][iferm]));
      means_ferm.add(*(phibar_ferm[iZ][iferm]));
      for (int im=0; im<9; im++) me_ferm.add(*(invcov_ferm[iZ][iferm][im]));
      gausConstraintsPDF[iZ][iferm][0] = new RooGaussianMomConstraint(Form("gausConstraintsPDF_Z%iFermion%i", iZ+1, iferm+1), Form("gausConstraintsPDF_Z%iFermion%i", iZ+1, iferm+1), vars_ferm, means_ferm, me_ferm, RooGaussianMomConstraint::kRhoLambdaPhi);
      constraints.add(*(gausConstraintsPDF[iZ][iferm][0]));

      RooArgList vars_fsr, means_fsr, me_fsr;
      vars_fsr.add(*(pT_fsr[iZ][iferm]));
      vars_fsr.add(*(lambda_fsr[iZ][iferm]));
      vars_fsr.add(*(phi_fsr[iZ][iferm]));
      means_fsr.add(*(pTbar_fsr[iZ][iferm]));
      means_fsr.add(*(lambdabar_fsr[iZ][iferm]));
      means_fsr.add(*(phibar_fsr[iZ][iferm]));
      for (int im=0; im<9; im++) me_fsr.add(*(invcov_fsr[iZ][iferm][im]));
      gausConstraintsPDF[iZ][iferm][1] = new RooGaussianMomConstraint(Form("gausConstraintsPDF_Z%iFermion%i", iZ+1, iferm+1), Form("gausConstraintsPDF_Z%iFermion%i", iZ+1, iferm+1), vars_fsr, means_fsr, me_fsr, RooGaussianMomConstraint::kRhoLambdaPhi);
      constraints.add(*(gausConstraintsPDF[iZ][iferm][1]));
    }
  }
  auxilliaryConstraintsPDF = new RooGenericPdf("auxilliaryConstraintsPDF", "@0*@1*@2", RooArgList(*(beta_Vdaughter[0]), *(beta_Vdaughter[1]), *massCuts)); // Will need to add m1, m2 cuts here as well!
  //constraints.add(*(auxilliaryConstraintsPDF));

  constraints.add(*(DiracDeltaPDF));

  //TEST
  constraintsPDF = new RooProdPdf("constraintsPDF", "constraintsPDF", constraints);
  //constraintsPDF = new RooProdPdf("constraintsPDF", "constraintsPDF", RooArgList(*varOne, *auxilliaryConstraintsPDF));
  //constraintsPDF = new RooProdPdf("constraintsPDF", "constraintsPDF", RooArgList(*varOne, *varOne));
}
void HMassConstraint::constructCompoundPdf(){
  RooArgList pdfList(*spinPDF, *constraintsPDF);
  PDF = new RooProdPdf("HMassConstraint_PDF", "HMassConstraint_PDF", pdfList);

//#if hmc_debug==1
//  RooArgList fastpdfList(*(simpleBWPDF[0]), *(simpleBWPDF[1]), *constraintsPDF);
//  fastPDF = new RooProdPdf("HMassConstraint_FastPDF", "HMassConstraint_FastPDF", fastpdfList);
//#else
  RooArgList fastpdfList(*bwProdPDF, *constraintsPDF);
  fastPDF = new RooProdPdf("HMassConstraint_FastPDF", "HMassConstraint_FastPDF", fastpdfList);
//#endif
}

void HMassConstraint::destroyVariables(){
  // Destroy the fit result, the ultimate culmination of all evil!
  deletePtr(fitResult);

  // Destroy in ~reverse order of creation
  deletePtr(h1);
  deletePtr(h2);
  deletePtr(hs);
  deletePtr(Phi);
  deletePtr(Phi1);
  deletePtr(Y);

  deletePtr(massCuts);
  for (int iZ=1; iZ>=0; iZ--) deletePtr(mAB[iZ]);

  deletePtr(m[2]);
  for (int iZ=1; iZ>=0; iZ--){
    deletePtr(beta_Vdaughter[iZ]);
    deletePtr(m[iZ]);

    for (int iferm=1; iferm>=0; iferm--){
      deletePtr(E_Hdaughter[iZ][iferm]);
      deletePtr(px_Hdaughter[iZ][iferm]);
      deletePtr(py_Hdaughter[iZ][iferm]);
      deletePtr(pz_Hdaughter[iZ][iferm]);
      deletePtr(m_Hdaughter[iZ][iferm]);

      deletePtr(E_fsr[iZ][iferm]);
      deletePtr(px_fsr[iZ][iferm]);
      deletePtr(py_fsr[iZ][iferm]);
      deletePtr(pz_fsr[iZ][iferm]);
      deletePtr(E_ferm[iZ][iferm]);
      deletePtr(px_ferm[iZ][iferm]);
      deletePtr(py_ferm[iZ][iferm]);
      deletePtr(pz_ferm[iZ][iferm]);

      for (int ix=0; ix<3; ix++){
        for (int iy=0; iy<3; iy++){
          deletePtr(invcov_ferm[iZ][iferm][3*ix+iy]);
          deletePtr(invcov_fsr[iZ][iferm][3*ix+iy]);
        }
      }
      
      deletePtr(pTobs_ferm[iZ][iferm]);
      deletePtr(lambdaobs_ferm[iZ][iferm]);
      deletePtr(phiobs_ferm[iZ][iferm]);
      deletePtr(pTobs_fsr[iZ][iferm]);
      deletePtr(lambdaobs_fsr[iZ][iferm]);
      deletePtr(phiobs_fsr[iZ][iferm]);

      deletePtr(pT_fsr[iZ][iferm]);
      deletePtr(lambda_fsr[iZ][iferm]);
      deletePtr(phi_fsr[iZ][iferm]);
      deletePtr(pTbar_fsr[iZ][iferm]);
      deletePtr(lambdabar_fsr[iZ][iferm]);
      deletePtr(phibar_fsr[iZ][iferm]);

      deletePtr(pT_ferm[iZ][iferm]);
      deletePtr(lambda_ferm[iZ][iferm]);
      deletePtr(phi_ferm[iZ][iferm]);
      deletePtr(massbar_ferm[iZ][iferm]);
      deletePtr(pTbar_ferm[iZ][iferm]);
      deletePtr(lambdabar_ferm[iZ][iferm]);
      deletePtr(phibar_ferm[iZ][iferm]);
    }
  }

  deletePtr(m1cut);
  deletePtr(m2cut);
  deletePtr(mFFcut);

  deletePtr(varOne);
  deletePtr(varZero);
}
void HMassConstraint::destroyDeltaFunctions(){
  deletePtr(DiracDeltaPDF);
  for (int iZ=1; iZ>=0; iZ--){
    for (int iferm=1; iferm>=0; iferm--){
      deletePtr(pTDeltaFcn_ferm[iZ][iferm]);
      deletePtr(lambdaDeltaFcn_ferm[iZ][iferm]);
      deletePtr(phiDeltaFcn_ferm[iZ][iferm]);
      deletePtr(pTDeltaFcn_fsr[iZ][iferm]);
      deletePtr(lambdaDeltaFcn_fsr[iZ][iferm]);
      deletePtr(phiDeltaFcn_fsr[iZ][iferm]);
    }
  }
}
void HMassConstraint::destroyPdfFactory(){
  // Delete test PDF if it exists.
#if hmc_debug==1
  for (int iZ=1; iZ>=0; iZ--) deletePtr(simpleBWPDF[iZ]);
#endif

  // Delete the bwProdPDF first since the measurables and parameters come from the pdfFactory
  deletePtr(bwProdPDF);

  // Only one of these is true; no need to delete pdfFactory since it is simply a mother-pointer to either of these.
  pdfFactory=0;
  deletePtr(hvvFactory);
  deletePtr(xvvFactory);
}
void HMassConstraint::destroyConstraintPdfs(){
  deletePtr(constraintsPDF);
  deletePtr(auxilliaryConstraintsPDF);
  for (int iZ=1; iZ>=0; iZ--){
    for (int iferm=1; iferm>=0; iferm--){
      for (int ifsr=1; ifsr>=0; ifsr--) deletePtr(gausConstraintsPDF[iZ][iferm][ifsr]);
    }
  }
}
void HMassConstraint::destroyCompoundPdf(){
  deletePtr(fastPDF);
  deletePtr(PDF);
}

void HMassConstraint::setFastPDF(bool useFastPDF_){ useFastPDF = useFastPDF_; }
void HMassConstraint::setPtEtaCuts(
  Double_t pTcut_muon_,
  Double_t etacut_muon_,
  Double_t pTcut_electron_,
  Double_t etacut_electron_,
  Double_t pTcut_jet_,
  Double_t etacut_jet_,
  Double_t pTcut_fsr_,
  Double_t etacut_fsr_
  ){
  pTcut_muon=pTcut_muon_;
  pTcut_electron=pTcut_electron_;
  pTcut_jet=pTcut_jet_;
  pTcut_fsr=pTcut_fsr_;

  if (etacut_muon_>0.) lambdacut_muon=piovertwo_val-2.*atan(exp(-etacut_muon_));
  else lambdacut_muon=piovertwo_val;
  if (etacut_electron_>0.) lambdacut_electron=piovertwo_val-2.*atan(exp(-etacut_electron_));
  else lambdacut_electron=piovertwo_val;
  if (etacut_jet_>0.) lambdacut_jet=piovertwo_val-2.*atan(exp(-etacut_jet_));
  else lambdacut_jet=piovertwo_val;
  if (etacut_fsr_>0.) lambdacut_fsr=piovertwo_val-2.*atan(exp(-etacut_fsr_));
  else lambdacut_fsr=piovertwo_val;
  // Actual setting of ranges is done per-event
}
void HMassConstraint::setM1M2Cuts(
  Double_t m1cut_,
  Double_t m2cut_,
  Double_t mFFcut_
  ){
  m1cut->setConstant(false);
  m2cut->setConstant(false);
  mFFcut->setConstant(false);
  m1cut->setVal(m1cut_);
  m2cut->setVal(m2cut_);
  mFFcut->setVal(mFFcut_);
  m1cut->setConstant(true);
  m2cut->setConstant(true);
  mFFcut->setConstant(true);
}


HMassConstraint::FitMomentumStrategy HMassConstraint::getFitMomentumStrategy(){ return fitMomStrategy_final; }
void HMassConstraint::setWorkingFitMomentumStrategy(HMassConstraint::FitMomentumStrategy fitMomStrategy_){ fitMomStrategy_final=fitMomStrategy_; }
void HMassConstraint::setFitMomentumStrategy(HMassConstraint::FitMomentumStrategy fitMomStrategy_){ fitMomStrategy=fitMomStrategy_; setWorkingFitMomentumStrategy(fitMomStrategy_); }
void HMassConstraint::testFitMomentumStrategy(Int_t& useFullCov, Int_t& FermFSRType, Int_t& fitpT, Int_t& fitlambda, Int_t& fitphi) const{

  if (
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pT ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_Lambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_Phi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pT ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_Lambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_Phi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi
    ) useFullCov=1;
  else useFullCov=0;

  if (
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pT ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pT ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pT ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pT ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pT
    ) fitpT=1;
  else fitpT=0;

  if (
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_Lambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_Lambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_Lambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Lambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Lambda
    ) fitlambda=1;
  else fitlambda=0;

  if (
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_Phi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_Phi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_Phi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Phi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Phi
    ) fitphi=1;
  else fitphi=0;

  if (
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_pT ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_Lambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_All_Phi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pT ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_Lambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_All_Phi
    ) FermFSRType=2;
  else if (
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
    fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pT ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Lambda ||
    fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Phi
    ) FermFSRType=1;
  else FermFSRType=0;

/*
  fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
  fitMomStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_All_pT ||
  fitMomStrategy_final==HMassConstraint::FullCov_All_Lambda ||
  fitMomStrategy_final==HMassConstraint::FullCov_All_Phi ||
  fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
  fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_pT ||
  fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_Lambda ||
  fitMomStrategy_final==HMassConstraint::FullCov_NoFSR_Phi ||
  fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
  fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
  fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
  fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
  fitMomStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_All_pT ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_All_Lambda ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_All_Phi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambda ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pT ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Lambda ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Phi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pT ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Lambda ||
  fitMomStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Phi
*/

}
void HMassConstraint::decrementMomentumStrategy(HMassConstraint::FitMomentumStrategy& strategy_){
  if (strategy_==HMassConstraint::FullCov_All_pTLambdaPhi) strategy_ = HMassConstraint::FullCov_All_pTLambda;
  else if (strategy_==HMassConstraint::FullCov_All_pTLambda) strategy_ = HMassConstraint::FullCov_All_pTPhi;
  else if (strategy_==HMassConstraint::FullCov_All_pTPhi) strategy_ = HMassConstraint::FullCov_All_LambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_All_LambdaPhi) strategy_ = HMassConstraint::FullCov_All_pT;
  else if (strategy_==HMassConstraint::FullCov_All_pT) strategy_ = HMassConstraint::FullCov_All_Lambda;
  else if (strategy_==HMassConstraint::FullCov_All_Lambda) strategy_ = HMassConstraint::FullCov_All_Phi;
  else if (strategy_==HMassConstraint::FullCov_All_Phi) strategy_ = HMassConstraint::FullCov_NoFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pTLambdaPhi) strategy_ = HMassConstraint::FullCov_NoFSR_pTLambda;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pTLambda) strategy_ = HMassConstraint::FullCov_NoFSR_pTPhi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pTPhi) strategy_ = HMassConstraint::FullCov_NoFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_LambdaPhi) strategy_ = HMassConstraint::FullCov_NoFSR_pT;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pT) strategy_ = HMassConstraint::FullCov_NoFSR_Lambda;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_Lambda) strategy_ = HMassConstraint::FullCov_NoFSR_Phi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_Phi) strategy_ = HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi) strategy_ = HMassConstraint::FullCov_OnlyFSR_pTLambda;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pTLambda) strategy_ = HMassConstraint::FullCov_OnlyFSR_pTPhi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pTPhi) strategy_ = HMassConstraint::FullCov_OnlyFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_LambdaPhi) strategy_ = HMassConstraint::FullCov_OnlyFSR_pT;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pT) strategy_ = HMassConstraint::FullCov_OnlyFSRR_Lambda;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSRR_Lambda) strategy_ = HMassConstraint::FullCov_OnlyFSR_Phi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_Phi) strategy_ = HMassConstraint::CovDiagonals_All_pTLambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pTLambdaPhi) strategy_ = HMassConstraint::CovDiagonals_All_pTLambda;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pTLambda) strategy_ = HMassConstraint::CovDiagonals_All_pTPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pTPhi) strategy_ = HMassConstraint::CovDiagonals_All_LambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_LambdaPhi) strategy_ = HMassConstraint::CovDiagonals_All_pT;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pT) strategy_ = HMassConstraint::CovDiagonals_All_Lambda;
  else if (strategy_==HMassConstraint::CovDiagonals_All_Lambda) strategy_ = HMassConstraint::CovDiagonals_All_Phi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_Phi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pTLambda;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pTLambda) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pTPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pTPhi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pT;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pT) strategy_ = HMassConstraint::CovDiagonals_NoFSR_Lambda;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_Lambda) strategy_ = HMassConstraint::CovDiagonals_NoFSR_Phi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_Phi) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pTLambda;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pTPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pT;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pT) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_Lambda;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_Lambda) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_Phi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_Phi) strategy_ = HMassConstraint::nFitMomentumStrategies;
  else strategy_ = HMassConstraint::nFitMomentumStrategies;
}
void HMassConstraint::incrementMomentumStrategy(HMassConstraint::FitMomentumStrategy& strategy_){
  if (strategy_==HMassConstraint::FullCov_All_pTLambdaPhi) strategy_ = HMassConstraint::nFitMomentumStrategies;
  else if (strategy_==HMassConstraint::FullCov_All_pTLambda) strategy_ = HMassConstraint::FullCov_All_pTLambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_All_pTPhi) strategy_ = HMassConstraint::FullCov_All_pTLambda;
  else if (strategy_==HMassConstraint::FullCov_All_LambdaPhi) strategy_ = HMassConstraint::FullCov_All_pTPhi;
  else if (strategy_==HMassConstraint::FullCov_All_pT) strategy_ = HMassConstraint::FullCov_All_LambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_All_Lambda) strategy_ = HMassConstraint::FullCov_All_pT;
  else if (strategy_==HMassConstraint::FullCov_All_Phi) strategy_ = HMassConstraint::FullCov_All_Lambda;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pTLambdaPhi) strategy_ = HMassConstraint::FullCov_All_Phi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pTLambda) strategy_ = HMassConstraint::FullCov_NoFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pTPhi) strategy_ = HMassConstraint::FullCov_NoFSR_pTLambda;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_LambdaPhi) strategy_ = HMassConstraint::FullCov_NoFSR_pTPhi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_pT) strategy_ = HMassConstraint::FullCov_NoFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_Lambda) strategy_ = HMassConstraint::FullCov_NoFSR_pT;
  else if (strategy_==HMassConstraint::FullCov_NoFSR_Phi) strategy_ = HMassConstraint::FullCov_NoFSR_Lambda;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi) strategy_ = HMassConstraint::FullCov_NoFSR_Phi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pTLambda) strategy_ = HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pTPhi) strategy_ = HMassConstraint::FullCov_OnlyFSR_pTLambda;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_LambdaPhi) strategy_ = HMassConstraint::FullCov_OnlyFSR_pTPhi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_pT) strategy_ = HMassConstraint::FullCov_OnlyFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSRR_Lambda) strategy_ = HMassConstraint::FullCov_OnlyFSR_pT;
  else if (strategy_==HMassConstraint::FullCov_OnlyFSR_Phi) strategy_ = HMassConstraint::FullCov_OnlyFSRR_Lambda;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pTLambdaPhi) strategy_ = HMassConstraint::FullCov_OnlyFSR_Phi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pTLambda) strategy_ = HMassConstraint::CovDiagonals_All_pTLambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pTPhi) strategy_ = HMassConstraint::CovDiagonals_All_pTLambda;
  else if (strategy_==HMassConstraint::CovDiagonals_All_LambdaPhi) strategy_ = HMassConstraint::CovDiagonals_All_pTPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_pT) strategy_ = HMassConstraint::CovDiagonals_All_LambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_All_Lambda) strategy_ = HMassConstraint::CovDiagonals_All_pT;
  else if (strategy_==HMassConstraint::CovDiagonals_All_Phi) strategy_ = HMassConstraint::CovDiagonals_All_Lambda;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi) strategy_ = HMassConstraint::CovDiagonals_All_Phi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pTLambda) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pTPhi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pTLambda;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pTPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_pT) strategy_ = HMassConstraint::CovDiagonals_NoFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_Lambda) strategy_ = HMassConstraint::CovDiagonals_NoFSR_pT;
  else if (strategy_==HMassConstraint::CovDiagonals_NoFSR_Phi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_Lambda;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi) strategy_ = HMassConstraint::CovDiagonals_NoFSR_Phi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pTLambda;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pTPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_pT) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_Lambda) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_pT;
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_Phi) strategy_ = HMassConstraint::CovDiagonals_OnlyFSR_Lambda;
  else strategy_ = HMassConstraint::nFitMomentumStrategies;
}

void HMassConstraint::setFitVVStrategy(HMassConstraint::FitVVStrategy fitVVStrategy_){ fitVVStrategy=fitVVStrategy_; }
void HMassConstraint::testFitVVStrategy(Int_t& fitV1, Int_t& fitV2) const{
  if (
    fitVVStrategy==HMassConstraint::Fit_All_V1V2 ||
    fitVVStrategy==HMassConstraint::Fit_All_V1
    ) fitV1=1;
  else fitV1=0;
  if (
    fitVVStrategy==HMassConstraint::Fit_All_V1V2 ||
    fitVVStrategy==HMassConstraint::Fit_All_V2
    ) fitV2=1;
  else fitV2=0;
}

void HMassConstraint::addDaughters(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR, bool fitRetry){ // Candidate supports jets as well! FSR is also a reco::Candidate daughter.
  // If the current trial is fresh, reset relevant variables.
  if (!fitRetry){
    setWorkingFitMomentumStrategy(fitMomStrategy);
    inputRaw_Fermion_FSR.clear();
  }

  // Check the fit strategy
  Int_t useFullCov, FermFSRType, fitpT, fitlambda, fitphi;
  testFitMomentumStrategy(useFullCov, FermFSRType, fitpT, fitlambda, fitphi);
  Int_t fitV1, fitV2;
  testFitVVStrategy(fitV1, fitV2);

  int ndaughters=0;
  // Initialize PDG id's and bar-momenta
  for (int ix=0; ix<2; ix++){ for (int iy=0; iy<2; iy++) pdgid_ferm[ix][iy]=pdgUnknown; }

  int iZ = 0;
  for (std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>::iterator dau=FermionWithFSR.begin(); dau<FermionWithFSR.end(); dau++){
    const reco::Candidate* fermion = (*dau).first;
    if (fermion==0){ cerr << "HMassConstraint::addDaughters : Daughter " << ndaughters << " pointer is 0!" << endl; break; }
    else{
      ndaughters++;
      if (ndaughters>4){ cerr << "HMassConstraint::addDaughters : Number of daughters (" << ndaughters << ") exceeds 4!" << endl; break; }
      else{
        if (!fitRetry) inputRaw_Fermion_FSR.push_back(*dau);

        // ndaughters=1..4
        int pdgid = fermion->pdgId();
        int iferm = 0;
        if (pdgid<0) iferm=1;
        if (pdgid_ferm[iZ][iferm]!=pdgUnknown) iZ=1-iZ; // If pdgid_ferm[iZ][iferm] is occupied, fill pdgid_ferm[1-iZ][iferm] instead.
        if (pdgid_ferm[iZ][iferm]!=pdgUnknown){// If iZ is changed to 1-iZ and pdgid_ferm[iZ][iferm] is still occupied, something is wrong.
          cerr << "The current particle with PDG id " << pdgid << " cannot be placed in any slot! Something went wrong. Please check the order of immediate Higgs daughters passed." << endl;
          cerr << "The allocated PDG id's are:" << endl;
          for (int jZ=0; jZ<2; jZ++){
            for (int jferm=0; jferm<2; jferm++) cerr << pdgid_ferm[jZ][jferm] << '\t';
            cerr << endl;
          }
          break;
        }

        // Set PDG id
        pdgid_ferm[iZ][iferm] = pdgid;
        if (!(abs(pdgid_ferm[iZ][iferm])==pdgEle || abs(pdgid_ferm[iZ][iferm])==pdgMu || abs(pdgid_ferm[iZ][iferm])==pdgTau)) pdgid_ferm[iZ][iferm]=pdgJet; // Needs to be more thorough if jets are passsed
#if hmc_debug==1
        cout << "HMassConstraint::addDaughters : Daughter with PDG id " << pdgid << " is assigned to iZ=" << iZ << " and iferm=" << iferm << " with effective id " << pdgid_ferm[iZ][iferm] << endl;
        cout << "HMassConstraint::addDaughters : Daughter (px, py, pz, E) = (" << fermion->px() << '\t' << fermion->py() << '\t' << fermion->pz() << '\t' << fermion->energy() << endl;
        cout << "HMassConstraint::addDaughters : Daughter (pT, theta, phi, mass) = (" << fermion->pt() << '\t' << fermion->theta() << " (lambda=" << (piovertwo_val-fermion->theta()) << ")" << '\t' << fermion->phi() << '\t' << fermion->mass() << endl;
#endif

        // Set bar-momenta ranges
        massbar_ferm[iZ][iferm]->setConstant(false);
        pTbar_ferm[iZ][iferm]->setConstant(false);
        lambdabar_ferm[iZ][iferm]->setConstant(false);
        phibar_ferm[iZ][iferm]->setConstant(false);
        massbar_ferm[iZ][iferm]->setRange(-sqrts, sqrts);
        pTbar_ferm[iZ][iferm]->setRange(0., sqrts);
        lambdabar_ferm[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val);
        phibar_ferm[iZ][iferm]->setRange(-pi_val, pi_val);
        massbar_ferm[iZ][iferm]->setVal(fermion->mass());
        pTbar_ferm[iZ][iferm]->setVal(fermion->pt());
        lambdabar_ferm[iZ][iferm]->setVal(piovertwo_val-fermion->theta());
        phibar_ferm[iZ][iferm]->setVal(fermion->phi());
        massbar_ferm[iZ][iferm]->setRange(fermion->mass(), fermion->mass());
        pTbar_ferm[iZ][iferm]->setRange(fermion->pt(), fermion->pt());
        lambdabar_ferm[iZ][iferm]->setRange(piovertwo_val-fermion->theta(), piovertwo_val-fermion->theta());
        phibar_ferm[iZ][iferm]->setRange(fermion->phi(), fermion->phi());
        massbar_ferm[iZ][iferm]->setConstant(true);
        pTbar_ferm[iZ][iferm]->setConstant(true);
        lambdabar_ferm[iZ][iferm]->setConstant(true);
        phibar_ferm[iZ][iferm]->setConstant(true);

        // Set observed momenta, should be the same as bar-momenta
        pTobs_ferm[iZ][iferm]->setConstant(false);
        lambdaobs_ferm[iZ][iferm]->setConstant(false);
        phiobs_ferm[iZ][iferm]->setConstant(false);
        pTobs_ferm[iZ][iferm]->setVal(fermion->pt());
        lambdaobs_ferm[iZ][iferm]->setVal(piovertwo_val-fermion->theta());
        phiobs_ferm[iZ][iferm]->setVal(fermion->phi());
        pTobs_ferm[iZ][iferm]->setConstant(true);
        lambdaobs_ferm[iZ][iferm]->setConstant(true);
        phiobs_ferm[iZ][iferm]->setConstant(true);

        // Set refit fermion momenta ranges
        bool fixAll = (iZ==0 && fitV1==0) || (iZ==1 && fitV2==0);
        pT_ferm[iZ][iferm]->setConstant(false);
        lambda_ferm[iZ][iferm]->setConstant(false);
        phi_ferm[iZ][iferm]->setConstant(false);
        if (abs(pdgid_ferm[iZ][iferm])==pdgEle){
          pT_ferm[iZ][iferm]->setRange(pTcut_electron, sqrts);
          lambda_ferm[iZ][iferm]->setRange(-lambdacut_electron, lambdacut_electron);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgMu){
          pT_ferm[iZ][iferm]->setRange(pTcut_muon, sqrts);
          lambda_ferm[iZ][iferm]->setRange(-lambdacut_muon, lambdacut_muon);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgJet){
          pT_ferm[iZ][iferm]->setRange(pTcut_jet, sqrts);
          lambda_ferm[iZ][iferm]->setRange(-lambdacut_jet, lambdacut_jet);
        }
        phi_ferm[iZ][iferm]->setRange(-pi_val, pi_val);
        // Initialize refit fermion momenta to the values of fermion bar-momenta
        pT_ferm[iZ][iferm]->setVal(pTbar_ferm[iZ][iferm]->getVal());
        lambda_ferm[iZ][iferm]->setVal(lambdabar_ferm[iZ][iferm]->getVal());
        phi_ferm[iZ][iferm]->setVal(phibar_ferm[iZ][iferm]->getVal());
        if (fixAll){
          pT_ferm[iZ][iferm]->setConstant(true);
          lambda_ferm[iZ][iferm]->setConstant(true);
          phi_ferm[iZ][iferm]->setConstant(true);
        }
        else{
          if (fitpT==0 || FermFSRType==1) pT_ferm[iZ][iferm]->setConstant(true);
          if (fitlambda==0 || FermFSRType==1) lambda_ferm[iZ][iferm]->setConstant(true);
          if (fitphi==0 || FermFSRType==1) phi_ferm[iZ][iferm]->setConstant(true);
        }

        // Get fermion covariance matrices in terms of pT, lambda and phi
        Double_t coefMat_ferm[9] = { 0 };
        if (FermFSRType!=1 && !fixAll){
          sortGetCovarianceMatrix(coefMat_ferm, fermion);

#if hmc_debug==1
          cout << "HMassConstraint::addDaughters : Daughter " << iZ << " / " << iferm << " input covariance matrix:" << endl;
          for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++) cout << coefMat_ferm[3*ix+iy] << '\t'; cout << endl; }
#endif
          if (coefMat_ferm[3*0+0]==0. && coefMat_ferm[3*0+1]==0. && coefMat_ferm[3*0+2]==0.) pT_ferm[iZ][iferm]->setConstant(true);
          if (coefMat_ferm[3*1+1]==0. && coefMat_ferm[3*1+1]==0. && coefMat_ferm[3*1+2]==0.) lambda_ferm[iZ][iferm]->setConstant(true);
          if (coefMat_ferm[3*2+2]==0. && coefMat_ferm[3*1+2]==0. && coefMat_ferm[3*2+2]==0.) phi_ferm[iZ][iferm]->setConstant(true);
          if (!pT_ferm[iZ][iferm]->isConstant()) pT_ferm[iZ][iferm]->setRange(max(pT_ferm[iZ][iferm]->getMin(), pT_ferm[iZ][iferm]->getVal()-5.*sqrt(coefMat_ferm[3*0+0])), min(pT_ferm[iZ][iferm]->getMax(), pT_ferm[iZ][iferm]->getVal()+5.*sqrt(coefMat_ferm[3*0+0])));
          if (!lambda_ferm[iZ][iferm]->isConstant()) lambda_ferm[iZ][iferm]->setRange(max(lambda_ferm[iZ][iferm]->getMin(), lambda_ferm[iZ][iferm]->getVal()-5.*sqrt(coefMat_ferm[3*1+1])), min(lambda_ferm[iZ][iferm]->getMax(), lambda_ferm[iZ][iferm]->getVal()+5.*sqrt(coefMat_ferm[3*1+1])));
          if (!phi_ferm[iZ][iferm]->isConstant()) phi_ferm[iZ][iferm]->setRange(max(phi_ferm[iZ][iferm]->getMin(), phi_ferm[iZ][iferm]->getVal()-5.*sqrt(coefMat_ferm[3*2+2])), min(phi_ferm[iZ][iferm]->getMax(), phi_ferm[iZ][iferm]->getVal()+5.*sqrt(coefMat_ferm[3*2+2])));

          strategicInvertCovarianceMatrix(useFullCov, fitpT, fitlambda, fitphi, coefMat_ferm);

          // Fix the gaussian pdf variables as necessary
          Int_t gaussianCode=1;
          if (pT_ferm[iZ][iferm]->isConstant()) gaussianCode *= RooGaussianMomConstraint::prime_var1;
          if (lambda_ferm[iZ][iferm]->isConstant()) gaussianCode *= RooGaussianMomConstraint::prime_var2;
          if (phi_ferm[iZ][iferm]->isConstant()) gaussianCode *= RooGaussianMomConstraint::prime_var3;
          if (gausConstraintsPDF[iZ][iferm][0]!=0) gausConstraintsPDF[iZ][iferm][0]->fixVariable(gaussianCode);
        }
        else{
          if (gausConstraintsPDF[iZ][iferm][0]!=0) gausConstraintsPDF[iZ][iferm][0]->fixVariable(RooGaussianMomConstraint::prime_var1*RooGaussianMomConstraint::prime_var2*RooGaussianMomConstraint::prime_var3);
        }
        setInverseCovarianceMatrix(iZ, iferm, 0, coefMat_ferm);

        // Do FSR here
        const pat::PFParticle* gamma = (*dau).second;
        Double_t coefMat_fsr[9] ={ 0 };
        if (gamma!=0){
#if hmc_debug==1
          cout << "HMassConstraint::addDaughters : An FSR is assigned to iZ=" << iZ << " and iferm=" << iferm << endl;
#endif
          // Set FSR bar-momenta ranges
          pTbar_fsr[iZ][iferm]->setConstant(false);
          lambdabar_fsr[iZ][iferm]->setConstant(false);
          phibar_fsr[iZ][iferm]->setConstant(false);
          pTbar_fsr[iZ][iferm]->setRange(0., sqrts);
          lambdabar_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val);
          phibar_fsr[iZ][iferm]->setRange(-pi_val, pi_val);
          pTbar_fsr[iZ][iferm]->setVal(gamma->pt());
          lambdabar_fsr[iZ][iferm]->setVal(piovertwo_val-gamma->theta());
          phibar_fsr[iZ][iferm]->setVal(gamma->phi());
          pTbar_fsr[iZ][iferm]->setRange(gamma->pt(), gamma->pt());
          lambdabar_fsr[iZ][iferm]->setRange(piovertwo_val-gamma->theta(), piovertwo_val-gamma->theta());
          phibar_fsr[iZ][iferm]->setRange(gamma->phi(), gamma->phi());
          pTbar_fsr[iZ][iferm]->setConstant(true);
          lambdabar_fsr[iZ][iferm]->setConstant(true);
          phibar_fsr[iZ][iferm]->setConstant(true);

          // Set observed momenta, should be the same as bar-momenta
          pTobs_fsr[iZ][iferm]->setConstant(false);
          lambdaobs_fsr[iZ][iferm]->setConstant(false);
          phiobs_fsr[iZ][iferm]->setConstant(false);
          pTobs_fsr[iZ][iferm]->setVal(gamma->pt());
          lambdaobs_fsr[iZ][iferm]->setVal(piovertwo_val-gamma->theta());
          phiobs_fsr[iZ][iferm]->setVal(gamma->phi());
          pTobs_fsr[iZ][iferm]->setConstant(true);
          lambdaobs_fsr[iZ][iferm]->setConstant(true);
          phiobs_fsr[iZ][iferm]->setConstant(true);

          // Set fsr ranges within the cuts and initialize
          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., sqrts); pT_fsr[iZ][iferm]->setVal(pTbar_fsr[iZ][iferm]->getVal());
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-lambdacut_fsr, lambdacut_fsr); lambda_fsr[iZ][iferm]->setVal(lambdabar_fsr[iZ][iferm]->getVal());
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); phi_fsr[iZ][iferm]->setVal(phibar_fsr[iZ][iferm]->getVal());
          if (fixAll){
            pT_fsr[iZ][iferm]->setConstant(true);
            lambda_fsr[iZ][iferm]->setConstant(true);
            phi_fsr[iZ][iferm]->setConstant(true);
          }
          else{
            if (fitpT==0 || FermFSRType==0) pT_fsr[iZ][iferm]->setConstant(true);
            if (fitlambda==0 || FermFSRType==0) lambda_fsr[iZ][iferm]->setConstant(true);
            if (fitphi==0 || FermFSRType==0) phi_fsr[iZ][iferm]->setConstant(true);
          }

          // Get fsr covariance matrices in terms of pT, lambda and phi
          if (FermFSRType!=0 && !fixAll){
            sortGetCovarianceMatrix(coefMat_fsr, fermion);

            if (coefMat_fsr[3*0+0]==0. && coefMat_fsr[3*0+1]==0. && coefMat_fsr[3*0+2]==0.) pT_fsr[iZ][iferm]->setConstant(true);
            if (coefMat_fsr[3*1+1]==0. && coefMat_fsr[3*1+1]==0. && coefMat_fsr[3*1+2]==0.) lambda_fsr[iZ][iferm]->setConstant(true);
            if (coefMat_fsr[3*2+2]==0. && coefMat_fsr[3*1+2]==0. && coefMat_fsr[3*2+2]==0.) phi_fsr[iZ][iferm]->setConstant(true);
            if (!pT_fsr[iZ][iferm]->isConstant()) pT_fsr[iZ][iferm]->setRange(max(pT_fsr[iZ][iferm]->getMin(), pT_fsr[iZ][iferm]->getVal()-5.*sqrt(coefMat_fsr[3*0+0])), max(pT_fsr[iZ][iferm]->getMax(), pT_fsr[iZ][iferm]->getVal()+5.*sqrt(coefMat_fsr[3*0+0])));
            if (!lambda_fsr[iZ][iferm]->isConstant()) lambda_fsr[iZ][iferm]->setRange(max(lambda_fsr[iZ][iferm]->getMin(), lambda_fsr[iZ][iferm]->getVal()-5.*sqrt(coefMat_fsr[3*1+1])), max(lambda_fsr[iZ][iferm]->getMax(), lambda_fsr[iZ][iferm]->getVal()+5.*sqrt(coefMat_fsr[3*1+1])));
            if (!phi_fsr[iZ][iferm]->isConstant()) phi_fsr[iZ][iferm]->setRange(max(phi_fsr[iZ][iferm]->getMin(), phi_fsr[iZ][iferm]->getVal()-5.*sqrt(coefMat_fsr[3*2+2])), max(phi_fsr[iZ][iferm]->getMax(), phi_fsr[iZ][iferm]->getVal()+5.*sqrt(coefMat_fsr[3*2+2])));

            strategicInvertCovarianceMatrix(useFullCov, fitpT, fitlambda, fitphi, coefMat_fsr);

            // Fix the gaussian pdf variables as necessary
            Int_t gaussianCode=1;
            if (pT_fsr[iZ][iferm]->isConstant()) gaussianCode *= RooGaussianMomConstraint::prime_var1;
            if (lambda_fsr[iZ][iferm]->isConstant()) gaussianCode *= RooGaussianMomConstraint::prime_var2;
            if (phi_fsr[iZ][iferm]->isConstant()) gaussianCode *= RooGaussianMomConstraint::prime_var3;
            if (gausConstraintsPDF[iZ][iferm][1]!=0) gausConstraintsPDF[iZ][iferm][1]->fixVariable(gaussianCode);
          }
          else{
            if (gausConstraintsPDF[iZ][iferm][1]!=0) gausConstraintsPDF[iZ][iferm][1]->fixVariable(RooGaussianMomConstraint::prime_var1*RooGaussianMomConstraint::prime_var2*RooGaussianMomConstraint::prime_var3);
          }
          setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
        }
        else{
          pTbar_fsr[iZ][iferm]->setConstant(false); pTbar_fsr[iZ][iferm]->setRange(0., sqrts); pTbar_fsr[iZ][iferm]->setVal(0.); pTbar_fsr[iZ][iferm]->setRange(0., 0.); pTbar_fsr[iZ][iferm]->setConstant(true);
          lambdabar_fsr[iZ][iferm]->setConstant(false); lambdabar_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambdabar_fsr[iZ][iferm]->setVal(0.); lambdabar_fsr[iZ][iferm]->setRange(0., 0.); lambdabar_fsr[iZ][iferm]->setConstant(true);
          phibar_fsr[iZ][iferm]->setConstant(false); phibar_fsr[iZ][iferm]->setRange(-pi_val, pi_val); phibar_fsr[iZ][iferm]->setVal(0.); phibar_fsr[iZ][iferm]->setRange(0., 0.); phibar_fsr[iZ][iferm]->setConstant(true);

          pTobs_fsr[iZ][iferm]->setConstant(false);
          lambdaobs_fsr[iZ][iferm]->setConstant(false);
          phiobs_fsr[iZ][iferm]->setConstant(false);
          pTobs_fsr[iZ][iferm]->setVal(0.);
          lambdaobs_fsr[iZ][iferm]->setVal(0.);
          phiobs_fsr[iZ][iferm]->setVal(0.);
          pTobs_fsr[iZ][iferm]->setConstant(true);
          lambdaobs_fsr[iZ][iferm]->setConstant(true);
          phiobs_fsr[iZ][iferm]->setConstant(true);

          // setRange below resets range from previous iteration
          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., sqrts); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setConstant(true);
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setConstant(true);
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-pi_val, pi_val); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setConstant(true);

          if (gausConstraintsPDF[iZ][iferm][1]!=0) gausConstraintsPDF[iZ][iferm][1]->fixVariable(RooGaussianMomConstraint::prime_var1*RooGaussianMomConstraint::prime_var2*RooGaussianMomConstraint::prime_var3);

          setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
        }

      }
    }
  }

  // Set those non-existing particles
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      if (pdgid_ferm[iZ][iferm]!=pdgUnknown) continue;
#if hmc_debug==1
      cout << "HMassConstraint::addDaughters : Particle at iZ=" << iZ << " and iferm=" << iferm << " is not assigned. Fixing it to 0." << endl;
#endif

      massbar_ferm[iZ][iferm]->setConstant(false); massbar_ferm[iZ][iferm]->setRange(0., sqrts); massbar_ferm[iZ][iferm]->setVal(0.); massbar_ferm[iZ][iferm]->setRange(0., 0.); massbar_ferm[iZ][iferm]->setConstant(true);
      pTbar_ferm[iZ][iferm]->setConstant(false); pTbar_ferm[iZ][iferm]->setRange(0., sqrts); pTbar_ferm[iZ][iferm]->setVal(0.); pTbar_ferm[iZ][iferm]->setRange(0., 0.); pTbar_ferm[iZ][iferm]->setConstant(true);
      lambdabar_ferm[iZ][iferm]->setConstant(false); lambdabar_ferm[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambdabar_ferm[iZ][iferm]->setVal(0.); lambdabar_ferm[iZ][iferm]->setRange(0., 0.); lambdabar_ferm[iZ][iferm]->setConstant(true);
      phibar_ferm[iZ][iferm]->setConstant(false); phibar_ferm[iZ][iferm]->setRange(-pi_val, pi_val); phibar_ferm[iZ][iferm]->setVal(0.); phibar_ferm[iZ][iferm]->setRange(0., 0.); phibar_ferm[iZ][iferm]->setConstant(true);
      pTbar_fsr[iZ][iferm]->setConstant(false); pTbar_fsr[iZ][iferm]->setRange(0., sqrts); pTbar_fsr[iZ][iferm]->setVal(0.); pTbar_fsr[iZ][iferm]->setRange(0., 0.); pTbar_fsr[iZ][iferm]->setConstant(true);
      lambdabar_fsr[iZ][iferm]->setConstant(false); lambdabar_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambdabar_fsr[iZ][iferm]->setVal(0.); lambdabar_fsr[iZ][iferm]->setRange(0., 0.); lambdabar_fsr[iZ][iferm]->setConstant(true);
      phibar_fsr[iZ][iferm]->setConstant(false); phibar_fsr[iZ][iferm]->setRange(-pi_val, pi_val); phibar_fsr[iZ][iferm]->setVal(0.); phibar_fsr[iZ][iferm]->setRange(0., 0.); phibar_fsr[iZ][iferm]->setConstant(true);

      // Set observed momenta, should be the same as bar-momenta
      pTobs_ferm[iZ][iferm]->setConstant(false);
      lambdaobs_ferm[iZ][iferm]->setConstant(false);
      phiobs_ferm[iZ][iferm]->setConstant(false);
      pTobs_ferm[iZ][iferm]->setVal(0.);
      lambdaobs_ferm[iZ][iferm]->setVal(0.);
      phiobs_ferm[iZ][iferm]->setVal(0.);
      pTobs_ferm[iZ][iferm]->setConstant(true);
      lambdaobs_ferm[iZ][iferm]->setConstant(true);
      phiobs_ferm[iZ][iferm]->setConstant(true);
      pTobs_fsr[iZ][iferm]->setConstant(false);
      lambdaobs_fsr[iZ][iferm]->setConstant(false);
      phiobs_fsr[iZ][iferm]->setConstant(false);
      pTobs_fsr[iZ][iferm]->setVal(0.);
      lambdaobs_fsr[iZ][iferm]->setVal(0.);
      phiobs_fsr[iZ][iferm]->setVal(0.);
      pTobs_fsr[iZ][iferm]->setConstant(true);
      lambdaobs_fsr[iZ][iferm]->setConstant(true);
      phiobs_fsr[iZ][iferm]->setConstant(true);

      // setRange below resets range from previous iteration
      pT_ferm[iZ][iferm]->setConstant(false); pT_ferm[iZ][iferm]->setRange(0., sqrts); pT_ferm[iZ][iferm]->setVal(0.); pT_ferm[iZ][iferm]->setRange(0., 0.); pT_ferm[iZ][iferm]->setConstant(true);
      lambda_ferm[iZ][iferm]->setConstant(false); lambda_ferm[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambda_ferm[iZ][iferm]->setVal(0.); lambda_ferm[iZ][iferm]->setRange(0., 0.); lambda_ferm[iZ][iferm]->setConstant(true);
      phi_ferm[iZ][iferm]->setConstant(false); phi_ferm[iZ][iferm]->setRange(-pi_val, pi_val); phi_ferm[iZ][iferm]->setVal(0.); phi_ferm[iZ][iferm]->setRange(0., 0.); phi_ferm[iZ][iferm]->setConstant(true);
      pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., sqrts); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setConstant(true);
      lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setConstant(true);
      phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-pi_val, pi_val); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setConstant(true);

      for (int ifsr=0; ifsr<2; ifsr++){ if (gausConstraintsPDF[iZ][iferm][ifsr]!=0) gausConstraintsPDF[iZ][iferm][ifsr]->fixVariable(RooGaussianMomConstraint::prime_var1*RooGaussianMomConstraint::prime_var2*RooGaussianMomConstraint::prime_var3); }

      Double_t coefMat_ferm[9] ={ 0 };
      Double_t coefMat_fsr[9] ={ 0 };
      setInverseCovarianceMatrix(iZ, iferm, 0, coefMat_ferm);
      setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
    }
  }

#if hmc_debug==1
  cout << "=== SUMMARY OF FINAL STATES ===" << endl;
  cout << "| Fermions |" << endl;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      cout << "Z" << iZ << " daughter " << iferm;
      cout << " ";
      cout << "(px,py,pz,E) = ( ";
      cout << px_ferm[iZ][iferm]->getVal() << " " << py_ferm[iZ][iferm]->getVal() << " " << pz_ferm[iZ][iferm]->getVal() << " " << E_ferm[iZ][iferm]->getVal() << " )";
      cout << " = ";
      cout << "(pT,lambda,phi,m) = ( ";
      cout << pT_ferm[iZ][iferm]->getVal() << " " << lambda_ferm[iZ][iferm]->getVal() << " " << phi_ferm[iZ][iferm]->getVal() << " " << massbar_ferm[iZ][iferm]->getVal() << " )";
      cout << endl;
      cout << "Z" << iZ << " daughter " << iferm << " FSR";
      cout << " ";
      cout << "(px,py,pz,E) = ( ";
      cout << px_fsr[iZ][iferm]->getVal() << " " << py_fsr[iZ][iferm]->getVal() << " " << pz_fsr[iZ][iferm]->getVal() << " " << E_fsr[iZ][iferm]->getVal() << " )";
      cout << " = ";
      cout << "(pT,lambda,phi,m) = ( ";
      cout << pT_fsr[iZ][iferm]->getVal() << " " << lambda_fsr[iZ][iferm]->getVal() << " " << phi_fsr[iZ][iferm]->getVal() << " " << 0 << " )";
      cout << endl;
      cout << "Z" << iZ << " daughter " << iferm << " sum";
      cout << " ";
      cout << "(px,py,pz,E,m) = ( ";
      cout << px_Hdaughter[iZ][iferm]->getVal() << " " << py_Hdaughter[iZ][iferm]->getVal() << " " << pz_Hdaughter[iZ][iferm]->getVal() << " " << E_Hdaughter[iZ][iferm]->getVal() << " " << m_Hdaughter[iZ][iferm]->getVal() << " )";
      cout << endl;
    }
  }
  for (int iZ=0; iZ<2; iZ++) cout << "beta(" << iZ+1 << ") = " << beta_Vdaughter[iZ]->getVal() << endl;
  for (int iZ=0; iZ<2; iZ++) cout << "mAB(" << iZ+1 << ") = " << mAB[iZ]->getVal() << endl;
  for (int iZ=0; iZ<2; iZ++) cout << "m" << iZ+1 << " = " << m[iZ]->getVal() << endl;
  cout << "m12 = " << m[2]->getVal() << endl;
#endif
}

void HMassConstraint::fitTo(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR){
  addDaughters(FermionWithFSR);
  fit();
}
RooArgSet HMassConstraint::getDataVariables() const{
  RooArgSet data_args;
  if (intCodeStart%RooSpin::prime_h1 != 0) data_args.add(*(h1));
  if (intCodeStart%RooSpin::prime_h2 != 0) data_args.add(*(h2));
  if (intCodeStart%RooSpin::prime_hs != 0) data_args.add(*(hs));
  if (intCodeStart%RooSpin::prime_Phi != 0) data_args.add(*(Phi));
  if (intCodeStart%RooSpin::prime_Phi1 != 0) data_args.add(*(Phi1));
  //data_args.add(*(Y));

  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      if (!pT_ferm[iZ][iferm]->isConstant()) data_args.add(*(pTobs_ferm[iZ][iferm]));
      if (!lambda_ferm[iZ][iferm]->isConstant()) data_args.add(*(lambdaobs_ferm[iZ][iferm]));
      if (!phi_ferm[iZ][iferm]->isConstant()) data_args.add(*(phiobs_ferm[iZ][iferm]));
      if (!pT_fsr[iZ][iferm]->isConstant()) data_args.add(*(pTobs_fsr[iZ][iferm]));
      if (!lambda_fsr[iZ][iferm]->isConstant()) data_args.add(*(lambdaobs_fsr[iZ][iferm]));
      if (!phi_fsr[iZ][iferm]->isConstant()) data_args.add(*(phiobs_fsr[iZ][iferm]));
    }
  }
  return data_args;
}
RooDataSet* HMassConstraint::getDataset() const{
  RooArgSet data_args = getDataVariables();
  RooDataSet* data = 0;
  if (data_args.getSize()<=15 && data_args.getSize()>0){
    data = new RooDataSet("fittedHdaughters", "", data_args);
    data->add(data_args);
#if hmc_debug==1
    cout << "HMassConstraint::getDataset: Number of arguments: " << data_args.getSize() << endl;
    cout << "HMassConstraint::getDataset: Data:" << endl;
    data->Print("v");
    cout << endl;
#endif
  }
  return data;
}
void HMassConstraint::fit(){
#if hmc_debug==1
  //if (!(pdgid_ferm[0][0]==13 || pdgid_ferm[0][1]==13 || pdgid_ferm[1][0]==13 || pdgid_ferm[1][1]==13)) return;
  cout << "Begin HMassConstraint::fit()" << endl;
#endif
  // Get the data to fit
  RooDataSet* data = getDataset();

  // Get the PDF to use
  RooAbsPdf* activePDF;
  if (useFastPDF){
    activePDF = fastPDF;
#if hmc_debug==1
    cout << "HMassConstraint::fit: FastPDF option is active. The PDF being used is " << activePDF->GetName() << "." << endl;
    activePDF->Print("v"); cout << endl;
#endif
  }
  else{
    activePDF = PDF;
#if hmc_debug==1
    cout << "HMassConstraint::fit: FastPDF option is not active. The PDF being used is " << activePDF->GetName() << "." << endl;
    activePDF->Print("v"); cout << endl;
#endif
  }

  // Conditional variables (i.e. m12)
  RooArgSet conditionals;
  conditionals.add(*(m[2]));

  // Delete the fitResult in case it exists, just to avoid unwanted memory leaks.s
  deletePtr(fitResult);

  // Hold the factory parameters tight!
  if (hvvFactory!=0){ hvvFactory->makeParamsConst(true); hvvFactory->makeCouplingsConst(true); }
  else if (xvvFactory!=0){ xvvFactory->makeParamsConst(true); xvvFactory->makeCouplingsConst(true); }

  /******************************************** BEGIN FIT **************************************************/
  const Int_t minimizerSuccess=0;
  Int_t fitStatus=-999;
  // Set the fit commands
#if hmc_debug==1
  cout << "HMassConstraint::fit: Attempting first fit." << endl;
#endif
  RooLinkedList cmdList;
  RooCmdArg condObsArg = RooFit::ConditionalObservables(conditionals); cmdList.Add((TObject*)&condObsArg);
  RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*)&saveArg);
  RooCmdArg hesseArg = RooFit::Hesse(true); cmdList.Add((TObject*)&hesseArg);
  //RooCmdArg minimizerArg = RooFit::Minimizer("Minuit", "migrad"); cmdList.Add((TObject*)&minimizerArg);
  //RooCmdArg minimizerStrategyArg = RooFit::Strategy(0); cmdList.Add((TObject*)&minimizerStrategyArg);
  // Misc. options
#if hmc_debug==1
  RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*)&timerArg);
  RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*)&printlevelArg);
#else
  RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*)&printlevelArg);
#endif
  // Try with the default strategy
  if (data!=0){
    fitResult = activePDF->fitTo(*data, cmdList);
    fitStatus = fitResult->status();
  }
#if hmc_debug==1
  cout << "HMassConstraint::fit: First fit attempted." << endl;
  if (fitResult!=0){
    cout << "Fit status is " << fitStatus << endl;
    cout << "Fit properties:" << endl;
    fitResult->Print("v");
  }
#endif
  // If the default strategy fails, decrement it until there is no strategy.
  while (fitStatus!=minimizerSuccess){
    decrementMomentumStrategy(fitMomStrategy_final);
    if (fitMomStrategy_final==HMassConstraint::nFitMomentumStrategies) break;
    cerr << "HMassConstraint::fit: Fit did not converge or was not valid. Status changed from " << fitMomStrategy << " to " << fitMomStrategy_final << " to retry." << endl;

    vector<pair<const reco::Candidate*, const pat::PFParticle*>> pseudoinput = inputRaw_Fermion_FSR;
    addDaughters(pseudoinput, true); deletePtr(data); data = getDataset();
    deletePtr(fitResult);
    if (data!=0){
      fitResult = activePDF->fitTo(*data, cmdList);
      fitStatus = fitResult->status();
    }
    else fitStatus=-999;
  }
  // If decrementing the strategy fails, increment it instead until there is no strategy.
  if (fitStatus!=minimizerSuccess) setWorkingFitMomentumStrategy(fitMomStrategy);
  while (fitStatus!=minimizerSuccess){
    incrementMomentumStrategy(fitMomStrategy_final);
    if (fitMomStrategy_final==HMassConstraint::nFitMomentumStrategies) break;
    cerr << "HMassConstraint::fit: Fit did not converge. Status changed from " << fitMomStrategy << " to " << fitMomStrategy_final << " to retry." << endl;

    vector<pair<const reco::Candidate*, const pat::PFParticle*>> pseudoinput = inputRaw_Fermion_FSR;
    addDaughters(pseudoinput, true); deletePtr(data); data = getDataset();
    deletePtr(fitResult);
    if (data!=0){
      fitResult = activePDF->fitTo(*data, cmdList);
      fitStatus = fitResult->status();
    }
    else fitStatus=-999;
  }
  /********************************************  END FIT  **************************************************/

  Int_t nDimCovMat=0;
  Int_t nFinalPars=0;
  bool fitFinalSuccess = (fitResult!=0 && fitStatus==minimizerSuccess);

  if (fitFinalSuccess){
    TMatrixDSym mat_tmp = fitResult->covarianceMatrix();;
    nDimCovMat = mat_tmp.GetNcols();
#if hmc_debug==1
    cout << "Number of columns in the unprocessed covariance matrix is " << nDimCovMat << ". The covariance matrix is" << endl;
    for (int ix=0; ix<nDimCovMat; ix++){
      for (int iy=0; iy<nDimCovMat; iy++) cout << mat_tmp[ix][iy] << '\t';
      cout << endl;
    }
    cout << endl;
#endif
    fitCovMatrix.ResizeTo(nDimCovMat, nDimCovMat);
    fitCovMatrix = mat_tmp;
    fitFinalSuccess = fitFinalSuccess && (nDimCovMat>0);
  }
  if (!fitFinalSuccess) cout << "Fit did not converge after all trials. Default parameters are to be used." << endl;

  if (fitFinalSuccess){
    const RooArgList pars = fitResult->floatParsFinal();
    nFinalPars = pars.getSize();
    for (int ip=0; ip<nFinalPars; ip++){
      const RooAbsArg* arg = pars.at(ip);
      if (dynamic_cast<const RooRealVar*>(arg)==0){ cerr << "Parameter " << ip << " (" << arg->GetName() << ") is not a RooRealVar!" << endl; assert(0); }
#if hmc_debug==1
      else{
        cout << "Parameter " << arg->GetName() << " = " << dynamic_cast<const RooRealVar*>(arg)->getVal() << " +- " << dynamic_cast<const RooRealVar*>(arg)->getError() << endl;
      }
#endif
    }

    if (nFinalPars!=nDimCovMat){ cerr << "Number of columns in the covariance matrix (" << nDimCovMat << ") does not match with the number of final floating variables (" << nFinalPars << ")!" << endl; assert(0); }
    else{
#if hmc_debug==1
      cout << "Number of columns in the covariance matrix is " << nDimCovMat << ". The covariance matrix is" << endl;
      for (int ix=0; ix<nDimCovMat; ix++){
        for (int iy=0; iy<nDimCovMat; iy++) cout << fitCovMatrix[ix][iy] << '\t';
        cout << endl;
      }
      cout << endl;
#endif
      // Re-order the covariance matrix here
      standardOrderedFinalCovarianceMatrix(pars);
    }
  }

  // Relax the factory parameters and return to normal
  if (hvvFactory!=0){ hvvFactory->makeParamsConst(false); hvvFactory->makeCouplingsConst(false); }
  else if (xvvFactory!=0){ xvvFactory->makeParamsConst(false); xvvFactory->makeCouplingsConst(false); }
  deletePtr(data);
}

void HMassConstraint::sortGetCovarianceMatrix(double (&momCov)[9], const reco::Candidate* particle){
  const reco::GsfElectron* electron = dynamic_cast<const reco::GsfElectron*>(particle);
  const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(particle);
  const pat::Muon* muon = dynamic_cast<const pat::Muon*>(particle);
  //const pat::Electron* pat_electron = dynamic_cast<const pat::Electron*>(particle);
  const pat::Jet* jet = dynamic_cast<const pat::Jet*>(particle);

  if (jet!=0) return getCovarianceMatrix(momCov, jet);
  else if (muon!=0) return getCovarianceMatrix(momCov, muon);
  else if (pfcand!=0) return getCovarianceMatrix(momCov, pfcand); // This is a general PFCandidate, which is mostly for use as a photon
  //else if (pat_electron!=0) return getCovarianceMatrix(momCov, pat_electron);
  else if (electron!=0) return getCovarianceMatrix(momCov, electron);
  else{
#if hmc_debug==1
    cout << "HMassConstraint::sortGetCovarianceMatrix: Could not determine the type of particle " << particle << ". Setting all covariance matrices to 0." << endl;
#endif
    for(int i=0;i<9;i++) momCov[i]=0.;
  }
}
void HMassConstraint::getCovarianceMatrix(double (&momCov)[9], const reco::GsfElectron* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

#if hmc_debug==1
  cout << "Begin HMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle) with argument " << particle << endl;
#endif

  double lambda = piovertwo_val - particle->theta();
  double energyerr;
  if (particle->ecalDriven()) energyerr = particle->p4Error(reco::GsfElectron::P4_COMBINATION);
  else{
    double ecalEnergy = particle->correctedEcalEnergy();
    double err2 = 0.;
    if (particle->isEB()){
      err2 += (5.24e-02*5.24e-02)/ecalEnergy;
      err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
      err2 += 1.00e-02*1.00e-02;
    }
    else if (particle->isEE()){
      err2 += (1.46e-01*1.46e-01)/ecalEnergy;
      err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
      err2 += 1.94e-03*1.94e-03;
    }
    energyerr = ecalEnergy * sqrt(err2);
  }
#if hmc_debug==1
  cout << "HMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle): Energy error before GsfTrack: " << energyerr << endl;
#endif

  double pterr = energyerr*cos(lambda);
  const GsfTrack* gsftrack = &(*(particle->gsfTrack()));
  if (gsftrack!=0){
#if hmc_debug==1
    cout << "HMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle): GsfTrack " << gsftrack << " found!" << endl;
#endif
    double pterr_uncorrected = gsftrack->ptModeError();
    if (pterr_uncorrected==0.) pterr_uncorrected = pterr;
    double correction = 1.;
    if (pterr_uncorrected!=0.) correction = pow(pterr/pterr_uncorrected, 2);

    double trackCov[GsfTrack::dimensionMode*GsfTrack::dimensionMode];
    for (int ix=0; ix<GsfTrack::dimensionMode; ix++){
      for (int iy=0; iy<GsfTrack::dimensionMode; iy++){
        if (iy>=ix){
          trackCov[GsfTrack::dimensionMode*ix+iy] = gsftrack->covarianceMode(ix, iy);
          if ((ix==TrackBase::i_qoverp || ix==TrackBase::i_lambda) && (iy==TrackBase::i_qoverp || iy==TrackBase::i_lambda)) trackCov[GsfTrack::dimensionMode*ix+iy] *= correction;
          else if ((ix==TrackBase::i_qoverp || ix==TrackBase::i_lambda) || (iy==TrackBase::i_qoverp || iy==TrackBase::i_lambda)) trackCov[GsfTrack::dimensionMode*ix+iy] *= sqrt(correction);
        }
        else trackCov[GsfTrack::dimensionMode*ix+iy] = trackCov[GsfTrack::dimensionMode*iy+ix];
      }
    }

    double q = particle->charge();
    double qoverp = q/particle->p();
    double d_pT_d_qoverp;
    if (q==0.){
      q=1.;
      qoverp = q/particle->p();
      d_pT_d_qoverp = 1e12; // (1000 TeV)**2
    }
    else d_pT_d_qoverp = -q*cos(lambda)/pow(qoverp, 2); // ==-p*pT/q
    const double d_pT_d_lambda = -q*sin(lambda)/qoverp; // == -pz
    const double d_pT_d_phi = 0;
    const double d_lambda_d_qoverp = 0.;
    const double d_lambda_d_lambda = 1;
    const double d_lambda_d_phi = 0;
    const double d_phi_d_qoverp = 0;
    const double d_phi_d_lambda = 0;
    const double d_phi_d_phi = 1;

    momCov[3*0+0] = pterr; // pT, pT, no need to re-calculate
    momCov[3*0+1] = // pT, lambda
      d_pT_d_qoverp*d_lambda_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_lambda_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_lambda_d_lambda + d_pT_d_lambda*d_lambda_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_lambda) * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*0+2] = // pT, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*1+0] = momCov[3*0+1];// lambda, pT
    momCov[3*1+1] = // lambda, lambda
      d_lambda_d_qoverp*d_lambda_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_lambda_d_lambda*d_lambda_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_lambda_d_phi*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      2.*d_lambda_d_qoverp*d_lambda_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.*d_lambda_d_qoverp*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.*d_lambda_d_lambda*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*1+2] = // lambda, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*2+0] = momCov[3*0+2];// phi, pT
    momCov[3*2+1] = momCov[3*1+2];// phi, lambda
    momCov[3*2+2] = // phi, phi
      d_phi_d_qoverp*d_phi_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_phi_d_lambda*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_phi_d_phi*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      2.* d_phi_d_qoverp*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.* d_phi_d_qoverp*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.* d_phi_d_lambda*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;
  }
  else{
#if hmc_debug==1
    cout << "HMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle): No GsfTrack present." << endl;
#endif
    momCov[3*0+0] = pow(pterr, 2);
    // Everything else remains 0.
  }
}
void HMassConstraint::getCovarianceMatrix(double (&momCov)[9], const pat::Muon* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

  const double lambda = piovertwo_val - particle->theta();

  double pterr_uncorrected = particle->muonBestTrack()->ptError();
  double pterr = pterr_uncorrected;
  if(particle->hasUserFloat("correctedPtError")) pterr = particle->userFloat("correctedPtError");
  double correction = 1.;
  if (pterr_uncorrected!=0.) correction = pow(pterr/pterr_uncorrected, 2);

  double trackCov[TrackBase::dimension*TrackBase::dimension];
  for (int ix=0; ix<TrackBase::dimension; ix++){
    for (int iy=0; iy<TrackBase::dimension; iy++){
      if (iy>=ix) {
        trackCov[TrackBase::dimension*ix+iy] = particle->muonBestTrack()->covariance(ix, iy);
        if ((ix==TrackBase::i_qoverp || ix==TrackBase::i_lambda) && (iy==TrackBase::i_qoverp || iy==TrackBase::i_lambda)) trackCov[TrackBase::dimension*ix+iy] *= correction;
        else if ((ix==TrackBase::i_qoverp || ix==TrackBase::i_lambda) || (iy==TrackBase::i_qoverp || iy==TrackBase::i_lambda)) trackCov[TrackBase::dimension*ix+iy] *= sqrt(correction);
      }
      else trackCov[TrackBase::dimension*ix+iy] = trackCov[TrackBase::dimension*iy+ix];
    }
  }

  double q = particle->charge();
  double qoverp = q/particle->p();
  double d_pT_d_qoverp;
  if (q==0.){
    q=1.;
    qoverp = q/particle->p();
    d_pT_d_qoverp = 1e12; // (1000 TeV)**2
  }
  else d_pT_d_qoverp = -q*cos(lambda)/pow(qoverp, 2); // ==-p*pT/q
  const double d_pT_d_lambda = -q*sin(lambda)/qoverp; // == -pz
  const double d_pT_d_phi = 0;
  const double d_lambda_d_qoverp = 0.;
  const double d_lambda_d_lambda = 1;
  const double d_lambda_d_phi = 0;
  const double d_phi_d_qoverp = 0;
  const double d_phi_d_lambda = 0;
  const double d_phi_d_phi = 1;

  momCov[3*0+0] = pterr; // pT, pT, no need to re-calculate
  momCov[3*0+1] = // pT, lambda
    d_pT_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_lambda_d_lambda + d_pT_d_lambda*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;
  momCov[3*0+2] = // pT, phi
    d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  momCov[3*1+0] = momCov[3*0+1];// lambda, pT
  momCov[3*1+1] = // lambda, lambda
    d_lambda_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_lambda_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_lambda_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    2.*d_lambda_d_qoverp*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    2.*d_lambda_d_qoverp*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    2.*d_lambda_d_lambda*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;
  momCov[3*1+2] = // lambda, phi
    d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  momCov[3*2+0] = momCov[3*0+2];// phi, pT
  momCov[3*2+1] = momCov[3*1+2];// phi, lambda
  momCov[3*2+2] = // phi, phi
    d_phi_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_phi_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_phi_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    2.* d_phi_d_qoverp*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    2.* d_phi_d_qoverp*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    2.* d_phi_d_lambda*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;
}
void HMassConstraint::getCovarianceMatrix(double (&momCov)[9], const reco::PFCandidate* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

  double lambda = piovertwo_val - particle->theta();

  double energyerr = PFEnergyResolution().getEnergyResolutionEm(particle->energy(), particle->eta());
  double pterr = energyerr*cos(lambda);
  const reco::Track* track = &(*(particle->bestTrack()));
  if (track!=0){
#if hmc_debug==1
    cout << "HMassConstraint::getCovarianceMatrix(const reco::PFCandidate* particle): Track " << track << " found!" << endl;
#endif
    double pterr_uncorrected = track->ptError();
    if (pterr_uncorrected==0.) pterr_uncorrected = pterr;
    double correction = 1.;
    if (pterr_uncorrected!=0.) correction = pow(pterr/pterr_uncorrected, 2);

    double trackCov[TrackBase::dimension*TrackBase::dimension];
    for (int ix=0; ix<TrackBase::dimension; ix++){
      for (int iy=0; iy<TrackBase::dimension; iy++){
        if (iy>=ix){
          trackCov[TrackBase::dimension*ix+iy] = track->covariance(ix, iy);
          if ((ix==TrackBase::i_qoverp || ix==TrackBase::i_lambda) && (iy==TrackBase::i_qoverp || iy==TrackBase::i_lambda)) trackCov[TrackBase::dimension*ix+iy] *= correction;
          else if ((ix==TrackBase::i_qoverp || ix==TrackBase::i_lambda) || (iy==TrackBase::i_qoverp || iy==TrackBase::i_lambda)) trackCov[TrackBase::dimension*ix+iy] *= sqrt(correction);
        }
        else trackCov[TrackBase::dimension*ix+iy] = trackCov[TrackBase::dimension*iy+ix];
      }
    }

    double q = particle->charge();
    double qoverp = q/particle->p();
    double d_pT_d_qoverp;
    if (q==0.){
      q=1.;
      qoverp = q/particle->p();
      d_pT_d_qoverp = 1e12; // (1000 TeV)**2
    }
    else d_pT_d_qoverp = -q*cos(lambda)/pow(qoverp, 2); // ==-p*pT/q
    const double d_pT_d_lambda = -q*sin(lambda)/qoverp; // == -pz
    const double d_pT_d_phi = 0;
    const double d_lambda_d_qoverp = 0.;
    const double d_lambda_d_lambda = 1;
    const double d_lambda_d_phi = 0;
    const double d_phi_d_qoverp = 0;
    const double d_phi_d_lambda = 0;
    const double d_phi_d_phi = 1;

    momCov[3*0+0] = pterr; // pT, pT, no need to re-calculate
    momCov[3*0+1] = // pT, lambda
      d_pT_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_lambda_d_lambda + d_pT_d_lambda*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*0+2] = // pT, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*1+0] = momCov[3*0+1];// lambda, pT
    momCov[3*1+1] = // lambda, lambda
      d_lambda_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_lambda_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_lambda_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      2.*d_lambda_d_qoverp*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.*d_lambda_d_qoverp*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.*d_lambda_d_lambda*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*1+2] = // lambda, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*2+0] = momCov[3*0+2];// phi, pT
    momCov[3*2+1] = momCov[3*1+2];// phi, lambda
    momCov[3*2+2] = // phi, phi
      d_phi_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_phi_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_phi_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      2.* d_phi_d_qoverp*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.* d_phi_d_qoverp*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.* d_phi_d_lambda*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;
  }
  else{
#if hmc_debug==1
    cout << "HMassConstraint::getCovarianceMatrix(const reco::PFCandidate* particle): No track present." << endl;
#endif
  momCov[3*0+0] = pow(pterr, 2);
  // Everything else is 0. I know this cannot be correct, but let's work with it for now.
  }
}
void HMassConstraint::getCovarianceMatrix(double (&momCov)[9], const pat::Jet* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

  double lambda = piovertwo_val - particle->theta();
  double energyerr = particle->userFloat(jecString.Data());

  double C_p_p = energyerr*energyerr;
  // Everything else is 0. I know this cannot be correct, but let's work with it for now.
  // Should loop over track references if the object is PFJet using reco::TrackRefVector PFJet::getTrackRefs() and sum inverse covariance matrices
  momCov[3*0+0] = C_p_p*pow(cos(lambda), 2);
}

void HMassConstraint::invertOneDimensional(Int_t includeIndex, double (&momCov)[9]){
  double momCov_tmp[9]={ 0 };
  if (momCov[3*includeIndex+includeIndex]!=0.) momCov_tmp[3*includeIndex+includeIndex] = 1./momCov[3*includeIndex+includeIndex];
  for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++) momCov[3*ix+iy]=momCov_tmp[3*ix+iy]; }
}
void HMassConstraint::invertTwoDimensional(Int_t omitIndex, double (&momCov)[9]){
  double momCov_tmp[9]={ 0 };
  double determinant = 0;
  if (omitIndex==0){
    determinant = (momCov[3*1+1]*momCov[3*2+2]) - (momCov[3*1+2]*momCov[3*2+1]);
    momCov_tmp[3*1+1] = momCov[3*2+2];
    momCov_tmp[3*2+2] = momCov[3*1+1];
    momCov_tmp[3*1+2] = -momCov[3*1+2];
    momCov_tmp[3*2+1] = -momCov[3*2+1];
  }
  else if (omitIndex==1){
    determinant = (momCov[3*0+0]*momCov[3*2+2]) - (momCov[3*0+2]*momCov[3*2+0]);
    momCov_tmp[3*0+0] = momCov[3*2+2];
    momCov_tmp[3*2+2] = momCov[3*0+0];
    momCov_tmp[3*0+2] = -momCov[3*0+2];
    momCov_tmp[3*2+0] = -momCov[3*2+0];
  }
  else if (omitIndex==2){
    determinant = (momCov[3*0+0]*momCov[3*1+1]) - (momCov[3*0+1]*momCov[3*1+0]);
    momCov_tmp[3*0+0] = momCov[3*1+1];
    momCov_tmp[3*1+1] = momCov[3*0+0];
    momCov_tmp[3*0+1] = -momCov[3*0+1];
    momCov_tmp[3*1+0] = -momCov[3*1+0];
  }
  for (int ix=0; ix<3; ix++){
    for (int iy=0; iy<3; iy++){
      if (determinant!=0) momCov[3*ix+iy]=momCov_tmp[3*ix+iy]/determinant;
      else momCov[3*ix+iy]=0;
    }
  }
}
void HMassConstraint::invertThreeDimensional(double (&momCov)[9]){
  double momCov_tmp[9]={ 0 };
  double determinant = 0.
    + (momCov[3*0+0]*momCov[3*1+1]*momCov[3*2+2])
    + (momCov[3*0+1]*momCov[3*1+2]*momCov[3*2+0])
    + (momCov[3*0+2]*momCov[3*1+0]*momCov[3*2+1])
    - (momCov[3*0+2]*momCov[3*1+1]*momCov[3*2+0])
    - (momCov[3*0+1]*momCov[3*1+0]*momCov[3*2+2])
    - (momCov[3*0+0]*momCov[3*1+2]*momCov[3*2+1]);
  momCov_tmp[3*0+0] = momCov[3*1+1]*momCov[3*2+2] - momCov[3*1+2]*momCov[3*2+1];
  momCov_tmp[3*1+1] = momCov[3*0+0]*momCov[3*2+2] - momCov[3*0+2]*momCov[3*2+0];
  momCov_tmp[3*2+2] = momCov[3*0+0]*momCov[3*1+1] - momCov[3*0+1]*momCov[3*1+0];
  momCov_tmp[3*0+1] = momCov[3*0+2]*momCov[3*2+1] - momCov[3*0+1]*momCov[3*2+2];
  momCov_tmp[3*1+0] = momCov[3*1+2]*momCov[3*2+0] - momCov[3*2+2]*momCov[3*1+0];
  momCov_tmp[3*1+2] = momCov[3*1+0]*momCov[3*0+2] - momCov[3*0+0]*momCov[3*1+2];
  momCov_tmp[3*2+1] = momCov[3*2+0]*momCov[3*0+1] - momCov[3*2+1]*momCov[3*0+0];
  momCov_tmp[3*2+0] = momCov[3*2+1]*momCov[3*1+0] - momCov[3*1+1]*momCov[3*2+0];
  momCov_tmp[3*0+2] = momCov[3*0+1]*momCov[3*1+2] - momCov[3*0+2]*momCov[3*1+1];
  for (int ix=0; ix<3; ix++){
    for (int iy=0; iy<3; iy++){
      if (determinant!=0) momCov[3*ix+iy]=momCov_tmp[3*ix+iy]/determinant;
      else momCov[3*ix+iy]=0;
    }
  }
}
void HMassConstraint::strategicInvertCovarianceMatrix(Int_t useFullCov, Int_t fitpT, Int_t fitlambda, Int_t fitphi, double (&momCov)[9]){
  // Make sure that te diagonal elements are set to 0 before further computation if useFullCov==0
  if (useFullCov==0){ for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++){ if(ix!=iy) momCov[3*ix+iy]=0; } } }
  // Fit for only one of the observable types
  if((fitpT+fitlambda+fitphi)==1){
    if(fitpT==1) invertOneDimensional(0, momCov);
    else if(fitlambda==1) invertOneDimensional(1, momCov);
    else if(fitphi==1) invertOneDimensional(2, momCov);
  }
  // Constrain only one variable to initial values
  else if((fitpT+fitlambda+fitphi)==2){
    if(fitpT==0) invertTwoDimensional(0, momCov);
    else if(fitlambda==0) invertTwoDimensional(1, momCov);
    else if(fitphi==0) invertTwoDimensional(2, momCov);
  }
  // Release all three variables
  else invertThreeDimensional(momCov);
}
void HMassConstraint::setInverseCovarianceMatrix(Int_t iZ, Int_t iferm, Int_t fsrindex, Double_t momCov[9]){
  if (fsrindex==0){ for (int i=0; i<9; i++){ invcov_ferm[iZ][iferm][i]->setConstant(false); invcov_ferm[iZ][iferm][i]->setVal(momCov[i]); invcov_ferm[iZ][iferm][i]->setConstant(true); } }
  else{ for (int i=0; i<9; i++){ invcov_fsr[iZ][iferm][i]->setConstant(false); invcov_fsr[iZ][iferm][i]->setVal(momCov[i]); invcov_fsr[iZ][iferm][i]->setConstant(true); } }
#if hmc_debug==1
  cout << "Inverse of the covariance matrix for Z" << iZ << " daughter " << iferm << " is:" << endl;
  if (fsrindex==0){
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++) cout <<  invcov_ferm[iZ][iferm][3*i+j]->getVal() << '\t';
      cout << endl;
    }
  }
  else{
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++) cout <<  invcov_fsr[iZ][iferm][3*i+j]->getVal() << '\t';
      cout << endl;
    }
  }
#endif
}

bool HMassConstraint::standardOrderedFinalCovarianceMatrix(const RooArgList& pars){
  const int nDims = 24;
  TMatrixDSym finalMatrix(nDims);
  for (int ix=0; ix<nDims; ix++){
    for (int iy=0; iy<nDims; iy++) finalMatrix[ix][iy] = 0;
  }

  Int_t nFinalPars = pars.getSize();
  Int_t* order = new Int_t[nFinalPars];
  for (int ip=0; ip<nFinalPars; ip++){
    RooRealVar* arg = dynamic_cast<RooRealVar*>(pars.at(ip));
    Int_t index = fitParameterCorrespondance(arg);
    if (index<0){ delete[] order; return false; }
    order[ip]=index;
  }
  for (int ix=0; ix<nFinalPars; ix++){
    for (int iy=0; iy<nFinalPars; iy++){
      finalMatrix[order[ix]][order[iy]] = fitCovMatrix[ix][iy];
    }
  }
  fitCovMatrix.ResizeTo(nDims, nDims);
  fitCovMatrix=finalMatrix;

#if hmc_debug==1
  cout << "HMassConstraint::standardOrderedFinalCovarianceMatrix: Final covariance matrix is:" << endl;
  for (int ix=0; ix<nDims; ix++){
    for (int iy=0; iy<nDims; iy++) cout << fitCovMatrix[ix][iy] << '\t';
    cout << endl;
  }
#endif
  return true;
}
Int_t HMassConstraint::fitParameterCorrespondance(RooRealVar* par){
  Int_t index=-1;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      if (TString(par->GetName())==TString(pT_ferm[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+0;
      else if (TString(par->GetName())==TString(lambda_ferm[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+1;
      else if (TString(par->GetName())==TString(phi_ferm[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+2;
      else if (TString(par->GetName())==TString(pT_fsr[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+3;
      else if (TString(par->GetName())==TString(lambda_fsr[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+4;
      else if (TString(par->GetName())==TString(phi_fsr[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+5;
      if (index>=0) break;
    }
  }
  if (index<0) cerr << "HMassConstraint::fitParameterCorrespondance: Parameter " << par->GetName() << " not found!" << endl;
  return index;
}


Double_t HMassConstraint::d_Ek_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  Double_t pt_over_E=0;
  Double_t lambda=0;
  Double_t mass=0;
  Double_t pT=0;
  Double_t E=0;
  if (fsrindex==0){
    E = E_ferm[kZ][kferm]->getVal();
    pT = pT_ferm[kZ][kferm]->getVal();
    lambda = lambda_ferm[kZ][kferm]->getVal();
    mass = massbar_ferm[kZ][kferm]->getVal();
  }
  else{
    E = E_fsr[kZ][kferm]->getVal();
    pT = pT_fsr[kZ][kferm]->getVal();
    lambda = lambda_fsr[kZ][kferm]->getVal();
  }

  if (E==0. && mass!=0.) pt_over_E=0;
  else if (E==0. && mass==0.) pt_over_E = cos(lambda);
  else pt_over_E = pT/E;

  return (pt_over_E/pow(cos(lambda), 2));
}
Double_t HMassConstraint::d_pjk_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const{
  Double_t lambda=0;
  Double_t phi=0;
  if (fsrindex==0){
    phi = phi_ferm[kZ][kferm]->getVal();
    lambda = lambda_ferm[kZ][kferm]->getVal();
  }
  else{
    phi = phi_fsr[kZ][kferm]->getVal();
    lambda = lambda_fsr[kZ][kferm]->getVal();
  }

  if (j==0) return cos(phi);
  else if (j==1) return sin(phi);
  else return tan(lambda);
}
Double_t HMassConstraint::d_Ek_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  Double_t pz=0;
  if (fsrindex==0) pz = pz_ferm[kZ][kferm]->getVal();
  else pz = pz_fsr[kZ][kferm]->getVal();

  return (pz*d_Ek_d_pTk(kZ, kferm, fsrindex));
}
Double_t HMassConstraint::d_pjk_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const{
  if (j!=2) return 0.;
  Double_t pT=0;
  Double_t lambda=0;
  if (fsrindex==0){
    pT = pT_ferm[kZ][kferm]->getVal();
    lambda = lambda_ferm[kZ][kferm]->getVal();
  }
  else{
    pT = pT_fsr[kZ][kferm]->getVal();
    lambda = lambda_fsr[kZ][kferm]->getVal();
  }
  return (pT/pow(cos(lambda), 2));
}
Double_t HMassConstraint::d_Ek_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex) const{ return 0.; }
Double_t HMassConstraint::d_pjk_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const{
  if (j>=2) return 0.;
  Double_t pT=0;
  Double_t phi=0;
  if (fsrindex==0){
    pT = pT_ferm[kZ][kferm]->getVal();
    phi = phi_ferm[kZ][kferm]->getVal();
  }
  else{
    pT = pT_fsr[kZ][kferm]->getVal();
    phi = phi_fsr[kZ][kferm]->getVal();
  }
  if (j==0) return (-pT*sin(phi));
  else return (pT*cos(phi));
}

Double_t HMassConstraint::d_m123_d_pTk(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  if((imass<2 && kZ!=imass) || imass<0 || imass>2) return 0;

  Double_t mass=m[imass]->getVal();
  Double_t pj[4]={ 0 };
  for(int iferm=0; iferm<2; iferm++){
    pj[0] += px_Hdaughter[kZ][iferm]->getVal();
    pj[1] += py_Hdaughter[kZ][iferm]->getVal();
    pj[2] += pz_Hdaughter[kZ][iferm]->getVal();
    pj[3] += E_Hdaughter[kZ][iferm]->getVal();
  }
  if(imass==2){
    for(int iferm=0; iferm<2; iferm++){
      pj[0] += px_Hdaughter[1-kZ][iferm]->getVal();
      pj[1] += py_Hdaughter[1-kZ][iferm]->getVal();
      pj[2] += pz_Hdaughter[1-kZ][iferm]->getVal();
      pj[3] += E_Hdaughter[1-kZ][iferm]->getVal();
    }
  }

  Double_t value = d_Ek_d_pTk(kZ, kferm, fsrindex);
  if (mass!=0.){
    value *= pj[3]/mass;
    for (int j=0; j<3; j++) value -= pj[j]/mass*d_pjk_d_pTk(kZ, kferm, fsrindex, j);
  }
  return value;
}
Double_t HMassConstraint::d_m123_d_lambdak(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  if((imass<2 && kZ!=imass) || imass<0 || imass>2) return 0;

  Double_t mass=m[imass]->getVal();;
  Double_t pj[4]={ 0 };
  for(int iferm=0; iferm<2; iferm++){
    pj[0] += px_Hdaughter[kZ][iferm]->getVal();
    pj[1] += py_Hdaughter[kZ][iferm]->getVal();
    pj[2] += pz_Hdaughter[kZ][iferm]->getVal();
    pj[3] += E_Hdaughter[kZ][iferm]->getVal();
  }
  if(imass==2){
    for(int iferm=0; iferm<2; iferm++){
      pj[0] += px_Hdaughter[1-kZ][iferm]->getVal();
      pj[1] += py_Hdaughter[1-kZ][iferm]->getVal();
      pj[2] += pz_Hdaughter[1-kZ][iferm]->getVal();
      pj[3] += E_Hdaughter[1-kZ][iferm]->getVal();
    }
  }

  Double_t value = d_Ek_d_lambdak(kZ, kferm, fsrindex);
  if (mass!=0.){
    value *= pj[3]/mass;
    for (int j=0; j<3; j++) value -= pj[j]/mass*d_pjk_d_lambdak(kZ, kferm, fsrindex, j);
  }
  return value;
}
Double_t HMassConstraint::d_m123_d_phik(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  if((imass<2 && kZ!=imass) || imass<0 || imass>2) return 0;

  Double_t mass=m[imass]->getVal();;
  Double_t pj[4]={ 0 };
  for(int iferm=0; iferm<2; iferm++){
    pj[0] += px_Hdaughter[kZ][iferm]->getVal();
    pj[1] += py_Hdaughter[kZ][iferm]->getVal();
    pj[2] += pz_Hdaughter[kZ][iferm]->getVal();
    pj[3] += E_Hdaughter[kZ][iferm]->getVal();
  }
  if(imass==2){
    for(int iferm=0; iferm<2; iferm++){
      pj[0] += px_Hdaughter[1-kZ][iferm]->getVal();
      pj[1] += py_Hdaughter[1-kZ][iferm]->getVal();
      pj[2] += pz_Hdaughter[1-kZ][iferm]->getVal();
      pj[3] += E_Hdaughter[1-kZ][iferm]->getVal();
    }
  }

  Double_t value = d_Ek_d_phik(kZ, kferm, fsrindex);
  if (mass!=0.){
    value *= pj[3]/mass;
    for (int j=0; j<3; j++) value -= pj[j]*d_pjk_d_phik(kZ, kferm, fsrindex, j);
  }
  return value;
}

Double_t HMassConstraint::getRefittedMassError(Int_t imass) const{ // imass==0 is m1, imass==1 is m2, imass==2 is m12.
  Double_t value = 0;
  const Int_t NColCovMat=24;
  if(fitCovMatrix.GetNcols()!=NColCovMat) return value; // This means the fit was not run.
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      for (int ifsr=0; ifsr<2; ifsr++){
        Int_t ipt = 6*(2*iZ+iferm)+3*ifsr+0;
        Int_t ilambda = 6*(2*iZ+iferm)+3*ifsr+1;
        Int_t iphi = 6*(2*iZ+iferm)+3*ifsr+2;
        for (int jZ=0; jZ<2; jZ++){
          for (int jferm=0; jferm<2; jferm++){
            for (int jfsr=0; jfsr<2; jfsr++){
              Int_t jpt = 6*(2*jZ+jferm)+3*jfsr+0;
              Int_t jlambda = 6*(2*jZ+jferm)+3*jfsr+1;
              Int_t jphi = 6*(2*jZ+jferm)+3*jfsr+2;

              value += (fitCovMatrix[ipt][jpt])*d_m123_d_pTk(imass, iZ, iferm, ifsr)*d_m123_d_pTk(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[ipt][jlambda])*d_m123_d_pTk(imass, iZ, iferm, ifsr)*d_m123_d_lambdak(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[ipt][jphi])*d_m123_d_pTk(imass, iZ, iferm, ifsr)*d_m123_d_phik(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[ilambda][jpt])*d_m123_d_lambdak(imass, iZ, iferm, ifsr)*d_m123_d_pTk(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[ilambda][jlambda])*d_m123_d_lambdak(imass, iZ, iferm, ifsr)*d_m123_d_lambdak(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[ilambda][jphi])*d_m123_d_lambdak(imass, iZ, iferm, ifsr)*d_m123_d_phik(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[iphi][jpt])*d_m123_d_phik(imass, iZ, iferm, ifsr)*d_m123_d_pTk(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[iphi][jlambda])*d_m123_d_phik(imass, iZ, iferm, ifsr)*d_m123_d_lambdak(imass, jZ, jferm, jfsr);
              value += (fitCovMatrix[iphi][jphi])*d_m123_d_phik(imass, iZ, iferm, ifsr)*d_m123_d_phik(imass, jZ, jferm, jfsr);
            }
          }
        }
      }
    }
  }
  if (value<0.) value=0;
  value = sqrt(value);
  return value;
}
Double_t HMassConstraint::getRefittedMass(Int_t imass) const{
  if(imass>2) return 0;
  return m[imass]->getVal();
}
TLorentzVector HMassConstraint::getRefittedMomentum(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  TLorentzVector result;
  if(fsrindex==0) result.SetXYZT(px_ferm[iZ][iferm]->getVal(),py_ferm[iZ][iferm]->getVal(),pz_ferm[iZ][iferm]->getVal(),E_ferm[iZ][iferm]->getVal());
  else result.SetXYZT(px_fsr[iZ][iferm]->getVal(),py_fsr[iZ][iferm]->getVal(),pz_fsr[iZ][iferm]->getVal(),E_fsr[iZ][iferm]->getVal());
  return result;
}

