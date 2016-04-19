#include "HMassConstraint/include/HMassConstraint.h"


using namespace std;
using namespace reco;
using namespace pat;


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

  constructVariables();
  constructPdfFactory();
}

~HMassConstraint::HMassConstraint(){
  destroyPdfFactory();
  destroyVariables();
}


void HMassConstraint::constructVariables(){
  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;

  RooArgList m12_args;
  for (int iZ=0; iZ<2; iZ++){
    RooArgList mHdaughter_args;
    for (int iferm=0; iferm<2; iferm++){
      pT_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambda_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phi_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", -pi_val, pi_val);
      massbar_ferm[iZ][iferm] = new RooRealVar(Form("mass_Z%iFermion%i", iZ+1, iferm+1), "", 0., -sqrts, sqrts);
      pTbar_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdabar_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phibar_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", -pi_val, pi_val);

      pT_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambda_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phi_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -pi_val, pi_val);
      pTbar_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdabar_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phibar_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -pi_val, pi_val);

      massbar_ferm[iZ][iferm]->removeMin();
      massbar_ferm[iZ][iferm]->removeMax();

      pT_ferm[iZ][iferm]->removeMax();
      pTbar_ferm[iZ][iferm]->removeMax();
      pT_fsr[iZ][iferm]->removeMax();
      pTbar_fsr[iZ][iferm]->removeMax();

      E_ferm[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%i", iZ+1, iferm+1), "sqrt( pow(@0,2)+pow(@1/cos(@2),2) )", RooArgList(*(massbar_ferm[iZ][iferm]), *(pT_ferm[iZ][iferm]), *(lambda_ferm[iZ][iferm])));
      px_ferm[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      py_ferm[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      pz_ferm[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*tan(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(lambda_ferm[iZ][iferm])));
      E_fsr[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0/cos(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(lambda_fsr[iZ][iferm])));
      px_fsr[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      py_fsr[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      pz_fsr[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*tan(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(lambda_fsr[iZ][iferm])));
      E_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(E_ferm[iZ][iferm]), *(E_fsr[iZ][iferm])));
      px_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(px_ferm[iZ][iferm]), *(px_fsr[iZ][iferm])));
      py_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(py_ferm[iZ][iferm]), *(py_fsr[iZ][iferm])));
      pz_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(pz_ferm[iZ][iferm]), *(pz_fsr[iZ][iferm])));
      m_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("m_Z%iFermion%iRefit", iZ+1, iferm+1), "sqrt( abs(pow(@0,2)-pow(@1,2)-pow(@2,2)-pow(@3,2)) )*TMath::Sign(1.,pow(@0,2)-pow(@1,2)-pow(@2,2)-pow(@3,2))", RooArgList(*(E_Hdaughter[iZ][iferm]), *(px_Hdaughter[iZ][iferm]), *(py_Hdaughter[iZ][iferm]), *(pz_Hdaughter[iZ][iferm])));

      mHdaughter_args.add(*(E_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(px_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(py_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(pz_Hdaughter[iZ][iferm]));
      m12_args.add(mHdaughter_args);

      pTdiff_ferm[iZ][iferm] = new RooFormulaVar(Form("pTDiff_Z%iFermion%i", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(pT_ferm[iZ][iferm]), *(pTbar_ferm[iZ][iferm])));
      lambdadiff_ferm[iZ][iferm] = new RooFormulaVar(Form("lambdaDiff_Z%iFermion%i", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(lambda_ferm[iZ][iferm]), *(lambdabar_ferm[iZ][iferm])));
      phidiff_ferm[iZ][iferm] = new RooFormulaVar(Form("phiDiff_Z%iFermion%i", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(phi_ferm[iZ][iferm]), *(phibar_ferm[iZ][iferm])));
      pTdiff_fsr[iZ][iferm] = new RooFormulaVar(Form("pTDiff_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(pT_fsr[iZ][iferm]), *(pTbar_fsr[iZ][iferm])));
      lambdadiff_fsr[iZ][iferm] = new RooFormulaVar(Form("lambdaDiff_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(lambda_fsr[iZ][iferm]), *(lambdabar_fsr[iZ][iferm])));
      phidiff_fsr[iZ][iferm] = new RooFormulaVar(Form("phiDiff_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(phi_fsr[iZ][iferm]), *(phibar_fsr[iZ][iferm])));
    }

    // Construct m1/m2
    m[iZ] = new RooFormulaVar(Form("m%iRefit", iZ+1), "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mHdaughter_args);

    // This beta should multiply the spinPDF bc. having massive fermions have additional scale factors. These factors are even more relevant when FSR is present!
    beta_Vdaughter[iZ] = new RooFormulaVar(Form("betaV%iRefit", iZ+1), "sqrt( TMath::Max(1e-15, ( 1.-pow((@1+@2)/@0,2) )*( 1.-pow((@1-@2)/@0,2) ) ) )", RooArgList(*(m[iZ]), *(m_Hdaughter[iZ][0]), *(m_Hdaughter[iZ][1])));
  }
  // Construct m12
  m12 = new RooFormulaVar("m12Refit", "sqrt( pow(@0+@4+@8+@12,2)-pow(@1+@5+@9+@13,2)-pow(@2+@6+@10+@14,2)-pow(@3+@7+@11+@15,2) )", m12_args);

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
    xvvFactory = new TensorPdfFactory_XVV(measurables, Vdecay1, Vdecay2, true); // true for always-on-shell X
    xvvFactory->makeParamsConst(false); // So that we can play with couplings
    spinPDF = xvvFactory->getPDF();
    pdfFactory = xvvFactory;
  }
}

void HMassConstraint::destroyVariables(){
  // Destroy in ~reverse order of creation
  delete h1;
  delete h2;
  delete hs;
  delete Phi;
  delete Phi1;
  delete Y;

  delete m[2];
  for (int iZ=1; iZ>=0; iZ--){
    delete beta_Vdaughter[iZ];
    delete m[iZ];

    for (int iferm=1; iferm>=0; iferm--){
      delete pTdiff_ferm[iZ][iferm];
      delete lambdadiff_ferm[iZ][iferm];
      delete phidiff_ferm[iZ][iferm];
      delete pTdiff_fsr[iZ][iferm];
      delete lambdadiff_fsr[iZ][iferm];
      delete phidiff_fsr[iZ][iferm];

      delete E_Hdaughter[iZ][iHdaughter];
      delete px_Hdaughter[iZ][iHdaughter];
      delete py_Hdaughter[iZ][iHdaughter];
      delete pz_Hdaughter[iZ][iHdaughter];
      delete m_Hdaughter[iZ][iferm];

      delete E_fsr[iZ][iferm];
      delete px_fsr[iZ][iferm];
      delete py_fsr[iZ][iferm];
      delete pz_fsr[iZ][iferm];
      delete E_ferm[iZ][iferm];
      delete px_ferm[iZ][iferm];
      delete py_ferm[iZ][iferm];
      delete pz_ferm[iZ][iferm];

      delete pT_fsr[iZ][iferm];
      delete lambda_fsr[iZ][iferm];
      delete phi_fsr[iZ][iferm];
      delete pTbar_fsr[iZ][iferm];
      delete lambdabar_fsr[iZ][iferm];
      delete phibar_fsr[iZ][iferm];

      delete pT_ferm[iZ][iferm];
      delete lambda_ferm[iZ][iferm];
      delete phi_ferm[iZ][iferm];
      delete massbar_ferm[iZ][iferm];
      delete pTbar_ferm[iZ][iferm];
      delete lambdabar_ferm[iZ][iferm];
      delete phibar_ferm[iZ][iferm];
    }
  }
}
void HMassConstraint::destroyPdfFactory(){
  // Only one of these is true; no need to delete pdfFactory since it is simply a mother-pointer to either of these.
  pdfFactory=0;
  if (hvvFactory!=0){ delete hvvFactory; hvvFactory=0; }
  if (xvvFactory!=0){ delete xvvFactory; xvvFactory=0; }
}


