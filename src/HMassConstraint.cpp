#include "HMassConstraint/include/HMassConstraint.h"


using namespace std;
using namespace reco;
using namespace pat;


HMassConstraint::HMassConstraint(
  sqrts(sqrts_),
  RooSpin::VdecayType Vdecay1_,
  RooSpin::VdecayType Vdecay2_
  ) :
  Vdecay1(Vdecay1_),
  Vdecay2(Vdecay2_),
  intCodeStart(RooSpin::prime_h1*RooSpin::prime_h2*RooSpin::prime_hs*RooSpin::prime_Phi*RooSpin::prime_Phi1)
{
  constructVariables();
}

void HMassConstraint::constructVariables(){
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      mass_ferm[iZ][iferm] = new RooRealVar(Form("mass_Z%iFermion%i", iZ+1, iferm+1), "", 0., -sqrts, sqrts);
      pT_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      eta_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      phi_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      pTbar_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      etabar_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      phibar_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);

      pT_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      eta_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      phi_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      pTbar_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      etabar_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      phibar_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);

      mass_ferm[iZ][iferm]->removeMin();
      mass_ferm[iZ][iferm]->removeMax();

      pT_ferm[iZ][iferm]->removeMax();
      eta_ferm[iZ][iferm]->removeMax();
      phi_ferm[iZ][iferm]->removeMax();
      pTbar_ferm[iZ][iferm]->removeMax();
      etabar_ferm[iZ][iferm]->removeMax();
      phibar_ferm[iZ][iferm]->removeMax();
      pT_fsr[iZ][iferm]->removeMax();
      eta_fsr[iZ][iferm]->removeMax();
      phi_fsr[iZ][iferm]->removeMax();
      pTbar_fsr[iZ][iferm]->removeMax();
      etabar_fsr[iZ][iferm]->removeMax();
      phibar_fsr[iZ][iferm]->removeMax();

      E_ferm[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%i", iZ+1, iferm+1), "sqrt( pow(@0,2)+pow(@1*cosh(@2),2) )", RooArgList(*(mass_ferm[iZ][iferm]), *(pT_ferm[iZ][iferm]), *(eta_ferm[iZ][iferm])));
      px_ferm[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      py_ferm[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      pz_ferm[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*sinh(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(eta_ferm[iZ][iferm])));
      E_fsr[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%iFSR", iZ+1, iferm+1), "sqrt( pow(@0*cosh(@1),2) )", RooArgList(*(pT_fsr[iZ][iferm]), *(eta_fsr[iZ][iferm])));
      px_fsr[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      py_fsr[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      pz_fsr[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*sinh(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(eta_fsr[iZ][iferm])));

      pTdiff_ferm[iZ][iferm] = new RooFormulaVar(Form("pTDiff_Z%iFermion%i", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(pT_ferm[iZ][iferm]), *(pTbar_ferm[iZ][iferm])));
      etadiff_ferm[iZ][iferm] = new RooFormulaVar(Form("etaDiff_Z%iFermion%i", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(eta_ferm[iZ][iferm]), *(etabar_ferm[iZ][iferm])));
      phidiff_ferm[iZ][iferm] = new RooFormulaVar(Form("phiDiff_Z%iFermion%i", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(phi_ferm[iZ][iferm]), *(phibar_ferm[iZ][iferm])));
      pTdiff_fsr[iZ][iferm] = new RooFormulaVar(Form("pTDiff_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(pT_fsr[iZ][iferm]), *(pTbar_fsr[iZ][iferm])));
      etadiff_fsr[iZ][iferm] = new RooFormulaVar(Form("etaDiff_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(eta_fsr[iZ][iferm]), *(etabar_fsr[iZ][iferm])));
      phidiff_fsr[iZ][iferm] = new RooFormulaVar(Form("phiDiff_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0-@1)", RooArgList(*(phi_fsr[iZ][iferm]), *(phibar_fsr[iZ][iferm])));
    }
  }

  // Construct m1
  RooArgList m1_args;
  for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(E_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(E_fsr[iZ][iferm])); } };
  for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(px_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(px_fsr[iZ][iferm])); } };
  for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(py_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(py_fsr[iZ][iferm])); } };
  for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(pz_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<1; iZ++){ for (int iferm=0; iferm<2; iferm++){ m1_args.add(*(pz_fsr[iZ][iferm])); } };
  m1 = new RooFormulaVar("m1Refit", "sqrt( pow(@0+@1+@2+@3,2)-pow(@4+@5+@6+@7,2)-pow(@8+@9+@10+@11,2)-pow(@12+@13+@14+@15,2) )", m1_args);

  // Construct m2
  RooArgList m2_args;
  for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(E_ferm[iZ][iferm])); } }; for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(E_fsr[iZ][iferm])); } };
  for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(px_ferm[iZ][iferm])); } }; for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(px_fsr[iZ][iferm])); } };
  for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(py_ferm[iZ][iferm])); } }; for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(py_fsr[iZ][iferm])); } };
  for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(pz_ferm[iZ][iferm])); } }; for (int iZ=1; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m2_args.add(*(pz_fsr[iZ][iferm])); } };
  m2 = new RooFormulaVar("m2Refit", "sqrt( pow(@0+@1+@2+@3,2)-pow(@4+@5+@6+@7,2)-pow(@8+@9+@10+@11,2)-pow(@12+@13+@14+@15,2) )", m2_args);

  // Construct m12
  RooArgList m12_args;
  for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(E_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(E_fsr[iZ][iferm])); } };
  for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(px_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(px_fsr[iZ][iferm])); } };
  for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(py_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(py_fsr[iZ][iferm])); } };
  for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(pz_ferm[iZ][iferm])); } }; for (int iZ=0; iZ<2; iZ++){ for (int iferm=0; iferm<2; iferm++){ m12_args.add(*(pz_fsr[iZ][iferm])); } };
  m12 = new RooFormulaVar("m12Refit", "sqrt( pow(@0+@1+@2+@3+@4+@5+@6+@7,2)-pow(@8+@9+@10+@11+@12+@13+@14+@15,2)-pow(@16+@17+@18+@19+@20+@21+@22+@23,2)-pow(@24+@25+@26+@27+@28+@29+@30+@31,2) )", m12_args);

  // Variables integrated over
  hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{1}", -1, 1);
  h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{2}", -1, 1);
  Phi = new RooRealVar("Genhelphi", "#Phi", -TMath::Pi(), TMath::Pi());
  Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{1}", -TMath::Pi(), TMath::Pi());
  Y = new RooRealVar("GenY", "Y", 0);

  // Initialize the meaurables
  measurables.m1 = m1;
  measurables.m2 = m2;
  measurables.m12 = m12;
  measurables.h1 = h1;
  measurables.h2 = h2;
  measurables.Phi = Phi;
  measurables.hs = hs;
  measurables.Phi1 = Phi1;
  measurables.Y = Y;
}
void HMassConstraint::constructPdfFactory(){
  pdfFactory = new ScalarPdfFactory_ggH(measurables, false, Vdecay1, Vdecay2);
  ((ScalarPdfFactory_ggH*)pdfFactory)->makeParamsConst(false);
}

void HMassConstraint::destroyVariables(){
  // Destroy in reverse order of creation
  delete h1;
  delete h2;
  delete hs;
  delete Phi;
  delete Phi1;
  delete Y;
  delete m1;
  delete m2;
  delete m12;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      delete pTdiff_ferm[iZ][iferm];
      delete etadiff_ferm[iZ][iferm];
      delete phidiff_ferm[iZ][iferm];
      delete pTdiff_fsr[iZ][iferm];
      delete etadiff_fsr[iZ][iferm];
      delete phidiff_fsr[iZ][iferm];

      delete px_ferm[iZ][iferm];
      delete py_ferm[iZ][iferm];
      delete pz_ferm[iZ][iferm];
      delete E_ferm[iZ][iferm];
      delete px_fsr[iZ][iferm];
      delete py_fsr[iZ][iferm];
      delete pz_fsr[iZ][iferm];
      delete E_fsr[iZ][iferm];

      delete pT_ferm[iZ][iferm];
      delete eta_ferm[iZ][iferm];
      delete phi_ferm[iZ][iferm];
      delete pTbar_ferm[iZ][iferm];
      delete etabar_ferm[iZ][iferm];
      delete phibar_ferm[iZ][iferm];

      delete pT_fsr[iZ][iferm];
      delete eta_fsr[iZ][iferm];
      delete phi_fsr[iZ][iferm];
      delete pTbar_fsr[iZ][iferm];
      delete etabar_fsr[iZ][iferm];
      delete phibar_fsr[iZ][iferm];
    }
  }
}
void HMassConstraint::destroyPdfFactory(){
  delete pdfFactory;
}

virtual ~HMassConstraint::HMassConstraint(){
  destroyPdfFactory();
  destroyVariables();
}


