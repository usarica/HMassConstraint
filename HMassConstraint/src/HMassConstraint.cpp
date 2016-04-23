#include "HMassConstraint/HMassConstraint/include/HMassConstraint.h"

#ifndef hmc_debug
#define hmc_debug 1
#endif

using namespace std;
using namespace reco;
using namespace pat;
using namespace HMCconstants;


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

  setFitStrategy();
  setJECUserFloatString();

  constructVariables();
  setM1M2Cuts();
  constructPdfFactory();
  constructConstraintPdfs();
  constructCompoundPdf();
}

HMassConstraint::~HMassConstraint(){
  destroyCompoundPdf();
  destroyConstraintPdfs();
  destroyPdfFactory();
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

      RooArgList sum_product_ferm_fsr_args;
      for (int ix=0; ix<3; ix++){
        for(int iy=0;iy<3;iy++){
          RooArgList product_ferm_args;
          RooArgList product_fsr_args;

          product_ferm_args.add(*(invcov_ferm[iZ][iferm][3*ix+iy]));
          if(ix==0) product_ferm_args.add(*(pTdiff_ferm[iZ][iferm]));
          else if(ix==1) product_ferm_args.add(*(lambdadiff_ferm[iZ][iferm]));
          else product_ferm_args.add(*(phidiff_ferm[iZ][iferm]));
          if(ix==iy) diffproducts_ferm[iZ][iferm][3*ix+iy] = new RooFormulaVar(Form("%s_times_%s", product_ferm_args.at(ix)->GetName(), product_ferm_args.at(iy)->GetName()), "(@0*@1*@1)", product_ferm_args);
          else{
            if(iy==0) product_ferm_args.add(*(pTdiff_ferm[iZ][iferm]));
            else if(iy==1) product_ferm_args.add(*(lambdadiff_ferm[iZ][iferm]));
            else product_ferm_args.add(*(phidiff_ferm[iZ][iferm]));
            diffproducts_ferm[iZ][iferm][3*ix+iy] = new RooFormulaVar(Form("%s_times_%s", product_ferm_args.at(ix)->GetName(), product_ferm_args.at(iy)->GetName()), "(@0*@1*@2)", product_ferm_args);
          }
          sum_product_ferm_fsr_args.add(product_ferm_args);

          product_fsr_args.add(*(invcov_fsr[iZ][iferm][3*ix+iy]));
          if (ix==0) product_fsr_args.add(*(pTdiff_fsr[iZ][iferm]));
          else if(ix==1) product_fsr_args.add(*(lambdadiff_fsr[iZ][iferm]));
          else product_fsr_args.add(*(phidiff_fsr[iZ][iferm]));
          if(ix==iy) diffproducts_fsr[iZ][iferm][3*ix+iy] = new RooFormulaVar(Form("%s_times_%s", product_fsr_args.at(ix)->GetName(), product_fsr_args.at(iy)->GetName()), "(@0*@1*@1)", product_fsr_args);
          else{
            if(iy==0) product_fsr_args.add(*(pTdiff_fsr[iZ][iferm]));
            else if(iy==1) product_fsr_args.add(*(lambdadiff_fsr[iZ][iferm]));
            else product_fsr_args.add(*(phidiff_fsr[iZ][iferm]));
            diffproducts_fsr[iZ][iferm][3*ix+iy] = new RooFormulaVar(Form("%s_times_%s", product_fsr_args.at(ix)->GetName(), product_fsr_args.at(iy)->GetName()), "(@0*@1*@2)", product_fsr_args);
          }
          sum_product_ferm_fsr_args.add(product_fsr_args);
        }
      }
      // Notice that the sum below is multiplied a a factor -1/2!
      sumdiffproducts_ferm_fsr[iZ][iferm] = new RooFormulaVar(Form("sumdiffproducts_Z%iFermion%i", iZ+1, iferm+1), "-(@0+@1+@2+@3+@4+@5+@6+@7+@8+@9+@10+@11+@12+@13+@14+@15+@16+@17)/2.", sum_product_ferm_fsr_args);
    }

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
    int iferm = 1-iZ;
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

  // Construct the Gaussian contraint sums over all Zs and daughters
  sumdiffproducts_ferm_fsr_combined = new RooFormulaVar("sumdiffproducts_combined", "(@0+@1+@2+@3)", RooArgList(*(sumdiffproducts_ferm_fsr[0][0]), *(sumdiffproducts_ferm_fsr[0][1]), *(sumdiffproducts_ferm_fsr[1][0]), *(sumdiffproducts_ferm_fsr[1][1])));

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
    xvvFactory = new TensorPdfFactory_HVV(measurables, Vdecay1, Vdecay2, true); // true for always-on-shell X
    xvvFactory->makeParamsConst(false); // So that we can play with couplings
    spinPDF = xvvFactory->getPDF();
    pdfFactory = xvvFactory;
  }
}
void HMassConstraint::constructConstraintPdfs(){
  gausConstraintsPDF = new RooExponential("gausConstraintsPDF", "gausConstraintsPDF", *sumdiffproducts_ferm_fsr_combined, *varOne);
  auxilliaryConstraintsPDF = new RooGenericPdf("auxilliaryConstraintsPDF", "@0*@1*@2", RooArgList(*(beta_Vdaughter[0]), *(beta_Vdaughter[1]), *massCuts)); // Will need to add m1, m2 cuts here as well!
  constraintsPDF = new RooProdPdf("constraintsPDF", "constraintsPDF", RooArgList(*gausConstraintsPDF, *auxilliaryConstraintsPDF));
}
void HMassConstraint::constructCompoundPdf(){
  RooArgList pdfList(*spinPDF, *constraintsPDF);
  PDF = new RooProdPdf("HMassConstraint_PDF", "HMassConstraint_PDF", pdfList);
}

void HMassConstraint::destroyVariables(){
  // Destroy the fit result, the ultimate culmination of all evil!
  if (fitResult!=0) delete fitResult;

  // Destroy sums of (inv. cov.)^i^j*dx_idx_j
  delete sumdiffproducts_ferm_fsr_combined;

  // Destroy in ~reverse order of creation
  delete h1;
  delete h2;
  delete hs;
  delete Phi;
  delete Phi1;
  delete Y;

  delete massCuts;
  for (int iZ=1; iZ>=0; iZ--) delete mAB[iZ];

  delete m[2];
  for (int iZ=1; iZ>=0; iZ--){
    delete beta_Vdaughter[iZ];
    delete m[iZ];

    for (int iferm=1; iferm>=0; iferm--){

      delete sumdiffproducts_ferm_fsr[iZ][iferm];
      for(int ix=0;ix<3;ix++){
        for(int iy=0;iy<3;iy++){
          delete diffproducts_ferm[iZ][iferm][3*ix+iy];
          delete diffproducts_fsr[iZ][iferm][3*ix+iy];
        }
      }

      delete pTdiff_ferm[iZ][iferm];
      delete lambdadiff_ferm[iZ][iferm];
      delete phidiff_ferm[iZ][iferm];
      delete pTdiff_fsr[iZ][iferm];
      delete lambdadiff_fsr[iZ][iferm];
      delete phidiff_fsr[iZ][iferm];

      delete E_Hdaughter[iZ][iferm];
      delete px_Hdaughter[iZ][iferm];
      delete py_Hdaughter[iZ][iferm];
      delete pz_Hdaughter[iZ][iferm];
      delete m_Hdaughter[iZ][iferm];

      delete E_fsr[iZ][iferm];
      delete px_fsr[iZ][iferm];
      delete py_fsr[iZ][iferm];
      delete pz_fsr[iZ][iferm];
      delete E_ferm[iZ][iferm];
      delete px_ferm[iZ][iferm];
      delete py_ferm[iZ][iferm];
      delete pz_ferm[iZ][iferm];

      for (int ix=0; ix<3; ix++){
        for (int iy=0; iy<3; iy++){
          delete invcov_ferm[iZ][iferm][3*ix+iy];
          delete invcov_fsr[iZ][iferm][3*ix+iy];
        }
      }

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

  delete m1cut;
  delete m2cut;
  delete mFFcut;

  delete varOne;
  delete varZero;
}
void HMassConstraint::destroyPdfFactory(){
  // Only one of these is true; no need to delete pdfFactory since it is simply a mother-pointer to either of these.
  pdfFactory=0;
  if (hvvFactory!=0){ delete hvvFactory; hvvFactory=0; }
  if (xvvFactory!=0){ delete xvvFactory; xvvFactory=0; }
}
void HMassConstraint::destroyConstraintPdfs(){
  delete constraintsPDF;
  delete auxilliaryConstraintsPDF;
  delete gausConstraintsPDF;
}
void HMassConstraint::destroyCompoundPdf(){
  if (PDF!=0) delete PDF;
}

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


HMassConstraint::FitStrategy HMassConstraint::getFitStrategy(){ return fitStrategy_final; }
void HMassConstraint::setWorkingFitStrategy(HMassConstraint::FitStrategy fitStrategy_){ fitStrategy_final=fitStrategy_; }
void HMassConstraint::setFitStrategy(HMassConstraint::FitStrategy fitStrategy_){ fitStrategy=fitStrategy_; setWorkingFitStrategy(fitStrategy_); }
void HMassConstraint::testFitStrategy(Int_t& useFullCov, Int_t& FermFSRType, Int_t& fitpT, Int_t& fitlambda, Int_t& fitphi) const{

  if (
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pT ||
    fitStrategy_final==HMassConstraint::FullCov_All_Lambda ||
    fitStrategy_final==HMassConstraint::FullCov_All_Phi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pT ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_Lambda ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_Phi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi
    ) useFullCov=1;
  else useFullCov=0;

  if (
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pT ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pT ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pT ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pT ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pT
    ) fitpT=1;
  else fitpT=0;

  if (
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_Lambda ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_Lambda ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_Lambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Lambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Lambda
    ) fitlambda=1;
  else fitlambda=0;

  if (
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_Phi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_NoFSR_Phi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_Phi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Phi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Phi
    ) fitphi=1;
  else fitphi=0;

  if (
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_All_pT ||
    fitStrategy_final==HMassConstraint::FullCov_All_Lambda ||
    fitStrategy_final==HMassConstraint::FullCov_All_Phi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_pT ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_Lambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_All_Phi
    ) FermFSRType=2;
  else if (
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
    fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pT ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Lambda ||
    fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Phi
    ) FermFSRType=1;
  else FermFSRType=0;

/*
  fitStrategy_final==HMassConstraint::FullCov_All_pTLambdaPhi ||
  fitStrategy_final==HMassConstraint::FullCov_All_pTLambda ||
  fitStrategy_final==HMassConstraint::FullCov_All_pTPhi ||
  fitStrategy_final==HMassConstraint::FullCov_All_LambdaPhi ||
  fitStrategy_final==HMassConstraint::FullCov_All_pT ||
  fitStrategy_final==HMassConstraint::FullCov_All_Lambda ||
  fitStrategy_final==HMassConstraint::FullCov_All_Phi ||
  fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
  fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTLambda ||
  fitStrategy_final==HMassConstraint::FullCov_NoFSR_pTPhi ||
  fitStrategy_final==HMassConstraint::FullCov_NoFSR_LambdaPhi ||
  fitStrategy_final==HMassConstraint::FullCov_NoFSR_pT ||
  fitStrategy_final==HMassConstraint::FullCov_NoFSR_Lambda ||
  fitStrategy_final==HMassConstraint::FullCov_NoFSR_Phi ||
  fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
  fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTLambda ||
  fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pTPhi ||
  fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
  fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_pT ||
  fitStrategy_final==HMassConstraint::FullCov_OnlyFSRR_Lambda ||
  fitStrategy_final==HMassConstraint::FullCov_OnlyFSR_Phi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambdaPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_All_pTLambda ||
  fitStrategy_final==HMassConstraint::CovDiagonals_All_pTPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_All_LambdaPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_All_pT ||
  fitStrategy_final==HMassConstraint::CovDiagonals_All_Lambda ||
  fitStrategy_final==HMassConstraint::CovDiagonals_All_Phi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTLambda ||
  fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pTPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_pT ||
  fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Lambda ||
  fitStrategy_final==HMassConstraint::CovDiagonals_NoFSR_Phi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
  fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
  fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_pT ||
  fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Lambda ||
  fitStrategy_final==HMassConstraint::CovDiagonals_OnlyFSR_Phi
*/

}
void HMassConstraint::decrementStrategy(HMassConstraint::FitStrategy& strategy_){
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
  else if (strategy_==HMassConstraint::CovDiagonals_OnlyFSR_Phi) strategy_ = HMassConstraint::nFitStrategies;
  else strategy_ = HMassConstraint::nFitStrategies;
}
void HMassConstraint::incrementStrategy(HMassConstraint::FitStrategy& strategy_){
  if (strategy_==HMassConstraint::FullCov_All_pTLambdaPhi) strategy_ = HMassConstraint::nFitStrategies;
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
  else strategy_ = HMassConstraint::nFitStrategies;
}

void HMassConstraint::addDaughters(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR, bool fitRetry){ // Candidate supports jets as well! FSR is also a reco::Candidate daughter.
  // If the current trial is fresh, reset relevant variables.
  if (!fitRetry){
    setWorkingFitStrategy(fitStrategy);
    inputRaw_Fermion_FSR.clear();
  }

  // Check the fit strategy
  Int_t useFullCov, FermFSRType, fitpT, fitlambda, fitphi;
  testFitStrategy(useFullCov, FermFSRType, fitpT, fitlambda, fitphi);

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
#endif

        // Set bar-momenta
        massbar_ferm[iZ][iferm]->setConstant(false);
        pTbar_ferm[iZ][iferm]->setConstant(false);
        lambdabar_ferm[iZ][iferm]->setConstant(false);
        phibar_ferm[iZ][iferm]->setConstant(false);
        massbar_ferm[iZ][iferm]->setVal(fermion->mass());
        pTbar_ferm[iZ][iferm]->setVal(fermion->pt());
        lambdabar_ferm[iZ][iferm]->setVal(piovertwo_val-fermion->theta());
        phibar_ferm[iZ][iferm]->setVal(fermion->phi());
        massbar_ferm[iZ][iferm]->setConstant(true);
        pTbar_ferm[iZ][iferm]->setConstant(true);
        lambdabar_ferm[iZ][iferm]->setConstant(true);
        phibar_ferm[iZ][iferm]->setConstant(true);

        // Set refit fermion momenta ranges
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

        // Get fermion covariance matrices in terms of pT, lambda and phi
        Double_t coefMat_ferm[9] = { 0 };
        if (FermFSRType!=1){
          sortGetCovarianceMatrix(coefMat_ferm, fermion);
          strategicInvertCovarianceMatrix(useFullCov, fitpT, fitlambda, fitphi, coefMat_ferm);
        }
        setInverseCovarianceMatrix(iZ, iferm, 0, coefMat_ferm);

        // Do FSR here
        const pat::PFParticle* gamma = (*dau).second;
        Double_t coefMat_fsr[9] ={ 0 };
        if (gamma!=0){
#if hmc_debug==1
          cout << "HMassConstraint::addDaughters : An FSR is assigned to iZ=" << iZ << " and iferm=" << iferm << endl;
#endif
          pTbar_fsr[iZ][iferm]->setConstant(false); pTbar_fsr[iZ][iferm]->setVal(gamma->pt()); pTbar_fsr[iZ][iferm]->setConstant(true);
          lambdabar_fsr[iZ][iferm]->setConstant(false); lambdabar_fsr[iZ][iferm]->setVal(piovertwo_val-gamma->theta()); lambdabar_fsr[iZ][iferm]->setConstant(true);
          phibar_fsr[iZ][iferm]->setConstant(false); phibar_fsr[iZ][iferm]->setVal(gamma->phi()); phibar_fsr[iZ][iferm]->setConstant(true);

          // Set fsr ranges within the cuts and initialize
          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., sqrts); pT_fsr[iZ][iferm]->setVal(pTbar_fsr[iZ][iferm]->getVal());
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-lambdacut_fsr, lambdacut_fsr); lambda_fsr[iZ][iferm]->setVal(lambdabar_fsr[iZ][iferm]->getVal());
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); phi_fsr[iZ][iferm]->setVal(phibar_fsr[iZ][iferm]->getVal());

          // Get fsr covariance matrices in terms of pT, lambda and phi
          if (FermFSRType!=0){
            sortGetCovarianceMatrix(coefMat_fsr, fermion);
            strategicInvertCovarianceMatrix(useFullCov, fitpT, fitlambda, fitphi, coefMat_fsr);
          }
          setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
        }
        else{
          pTbar_fsr[iZ][iferm]->setConstant(false); pTbar_fsr[iZ][iferm]->setVal(0.); pTbar_fsr[iZ][iferm]->setConstant(true);
          lambdabar_fsr[iZ][iferm]->setConstant(false); lambdabar_fsr[iZ][iferm]->setVal(0.); lambdabar_fsr[iZ][iferm]->setConstant(true);
          phibar_fsr[iZ][iferm]->setConstant(false); phibar_fsr[iZ][iferm]->setVal(0.); phibar_fsr[iZ][iferm]->setConstant(true);

          // setRange below resets range from previous iteration
          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setConstant(true);
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setConstant(true);
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setConstant(true);

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

      massbar_ferm[iZ][iferm]->setConstant(false); massbar_ferm[iZ][iferm]->setVal(0.); massbar_ferm[iZ][iferm]->setConstant(true);
      pTbar_ferm[iZ][iferm]->setConstant(false); pTbar_ferm[iZ][iferm]->setVal(0.); pTbar_ferm[iZ][iferm]->setConstant(true);
      lambdabar_ferm[iZ][iferm]->setConstant(false); lambdabar_ferm[iZ][iferm]->setVal(0.); lambdabar_ferm[iZ][iferm]->setConstant(true);
      phibar_ferm[iZ][iferm]->setConstant(false); phibar_ferm[iZ][iferm]->setVal(0.); phibar_ferm[iZ][iferm]->setConstant(true);
      pTbar_fsr[iZ][iferm]->setConstant(false); pTbar_fsr[iZ][iferm]->setVal(0.); pTbar_fsr[iZ][iferm]->setConstant(true);
      lambdabar_fsr[iZ][iferm]->setConstant(false); lambdabar_fsr[iZ][iferm]->setVal(0.); lambdabar_fsr[iZ][iferm]->setConstant(true);
      phibar_fsr[iZ][iferm]->setConstant(false); phibar_fsr[iZ][iferm]->setVal(0.); phibar_fsr[iZ][iferm]->setConstant(true);

      // setRange below resets range from previous iteration
      pT_ferm[iZ][iferm]->setConstant(false); pT_ferm[iZ][iferm]->setRange(0., 0.); pT_ferm[iZ][iferm]->setVal(0.); pT_ferm[iZ][iferm]->setConstant(true);
      lambda_ferm[iZ][iferm]->setConstant(false); lambda_ferm[iZ][iferm]->setRange(0., 0.); lambda_ferm[iZ][iferm]->setVal(0.); lambda_ferm[iZ][iferm]->setConstant(true);
      phi_ferm[iZ][iferm]->setConstant(false); phi_ferm[iZ][iferm]->setRange(0., 0.); phi_ferm[iZ][iferm]->setVal(0.); phi_ferm[iZ][iferm]->setConstant(true);
      pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setConstant(true);
      lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setConstant(true);
      phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setConstant(true);

      Double_t coefMat_ferm[9] ={ 0 };
      Double_t coefMat_fsr[9] ={ 0 };
      setInverseCovarianceMatrix(iZ, iferm, 0, coefMat_ferm);
      setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
    }
  }

}

RooDataSet* HMassConstraint::getDataset() const{
  RooArgSet data_args;
  if (intCodeStart%RooSpin::prime_h1 != 0) data_args.add(*(h1));
  if (intCodeStart%RooSpin::prime_h2 != 0) data_args.add(*(h2));
  if (intCodeStart%RooSpin::prime_hs != 0) data_args.add(*(hs));
  if (intCodeStart%RooSpin::prime_Phi != 0) data_args.add(*(Phi));
  if (intCodeStart%RooSpin::prime_Phi1 != 0) data_args.add(*(Phi1));
  //data_args.add(*(Y));

  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      data_args.add(*(pT_ferm[iZ][iferm]));
      data_args.add(*(lambda_ferm[iZ][iferm]));
      data_args.add(*(phi_ferm[iZ][iferm]));
    }
  }
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      data_args.add(*(pT_fsr[iZ][iferm]));
      data_args.add(*(lambda_fsr[iZ][iferm]));
      data_args.add(*(phi_fsr[iZ][iferm]));
    }
  }

  RooDataSet* data = new RooDataSet("fittedHdaughters", "", data_args);
  data->add(data_args);
  return data;
}
void HMassConstraint::fit(){
  RooDataSet* data = getDataset();

  RooArgSet conditionals;
  conditionals.add(*(m[2]));

  if (fitResult!=0) delete fitResult;

  // Hold the factory parameters tight!
  if (hvvFactory!=0){ hvvFactory->makeParamsConst(true); hvvFactory->makeCouplingsConst(true); }
  else if (xvvFactory!=0){ xvvFactory->makeParamsConst(true); xvvFactory->makeCouplingsConst(true); }

  /******************************************** BEGIN FIT **************************************************/
  const Int_t minimizerSuccess=2;
  // Try with the default strategy
  fitResult = PDF->fitTo(*data, RooFit::ConditionalObservables(conditionals), RooFit::Save(true), RooFit::PrintLevel(-1));
  // If the default strategy fails, decrement it until there is no strategy.
  while (!(fitResult->status()>=minimizerSuccess) && fitResult!=0){
    decrementStrategy(fitStrategy_final);
    if (fitStrategy_final==HMassConstraint::nFitStrategies) break;
    cout << "Fit did not converge. Status changed from " << fitStrategy << " to " << fitStrategy << "to retry." << endl;

    vector<pair<const reco::Candidate*, const pat::PFParticle*>> pseudoinput = inputRaw_Fermion_FSR;
    addDaughters(pseudoinput, true); if(data!=0) delete data; data = getDataset();
    delete fitResult;
    fitResult = PDF->fitTo(*data, RooFit::ConditionalObservables(conditionals), RooFit::Save(true), RooFit::PrintLevel(-1));
  }
  // If decrementing the strategy fails, increment it instead until there is no strategy.
  if (!(fitResult->status()>=minimizerSuccess) && fitResult!=0) setWorkingFitStrategy(fitStrategy);
  while (!(fitResult->status()>=minimizerSuccess) && fitResult!=0){
    incrementStrategy(fitStrategy_final);
    if (fitStrategy_final==HMassConstraint::nFitStrategies) break;
    cout << "Fit did not converge. Status changed from " << fitStrategy << " to " << fitStrategy << "to retry." << endl;

    vector<pair<const reco::Candidate*, const pat::PFParticle*>> pseudoinput = inputRaw_Fermion_FSR;
    addDaughters(pseudoinput, true); if (data!=0) delete data; data = getDataset();
    delete fitResult;
    fitResult = PDF->fitTo(*data, RooFit::ConditionalObservables(conditionals), RooFit::Save(true), RooFit::PrintLevel(-1));
  }
  /********************************************  END FIT  **************************************************/
  if (fitResult!=0 && fitResult->status()>=minimizerSuccess) fitCovMatrix = fitResult->covarianceMatrix();
  else{
    cout << "Fit did not converge after all trials. Default parameters are to be used." << endl;
  }

  // Relax the factory parameters.
  if (hvvFactory!=0){ hvvFactory->makeParamsConst(false); hvvFactory->makeCouplingsConst(false); }
  else if (xvvFactory!=0){ xvvFactory->makeParamsConst(false); xvvFactory->makeCouplingsConst(false); }
  if (data!=0) delete data;

  cout << "Number of columns in the covariance matrix is " << fitCovMatrix.GetNcols() << endl;
  // LEFT HERE
}


void HMassConstraint::sortGetCovarianceMatrix(double(&momCov)[9], const reco::Candidate* particle){
  const reco::GsfElectron* electron = dynamic_cast<const reco::GsfElectron*>(particle);
  const pat::Muon* muon = dynamic_cast<const pat::Muon*>(particle);
  const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(particle);
  const pat::Jet* jet = dynamic_cast<const pat::Jet*>(particle);

  if (jet!=0) return getCovarianceMatrix(momCov, jet);
  else if (muon!=0) return getCovarianceMatrix(momCov, muon);
  else if (pfcand!=0) return getCovarianceMatrix(momCov, pfcand); // This is a general PFCandidate, which is mostly for use as a photon
  else if (electron!=0) return getCovarianceMatrix(momCov, electron);
  else{ for(int i=0;i<9;i++) momCov[i]=0.; }
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
void HMassConstraint::invertThreeDimensional(double(&momCov)[9]){
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
  else if (E==0. && mass==0.) pt_over_E = sqrt(1./(1.+pow(tan(lambda), 2)));
  else pt_over_E = pT/E;

  return (pt_over_E*(1.+pow(tan(lambda), 2)));
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
  Double_t pt_over_E=0;
  Double_t lambda=0;
  Double_t mass=0;
  Double_t pT=0;
  Double_t pz=0;
  Double_t E=0;
  if (fsrindex==0){
    E = E_ferm[kZ][kferm]->getVal();
    pz = pz_ferm[kZ][kferm]->getVal();
    pT = pT_ferm[kZ][kferm]->getVal();
    lambda = lambda_ferm[kZ][kferm]->getVal();
    mass = massbar_ferm[kZ][kferm]->getVal();
  }
  else{
    E = E_fsr[kZ][kferm]->getVal();
    pz = pz_fsr[kZ][kferm]->getVal();
    pT = pT_fsr[kZ][kferm]->getVal();
    lambda = lambda_fsr[kZ][kferm]->getVal();
  }

  if (E==0. && mass!=0.) pt_over_E=0;
  else if (E==0. && mass==0.) pt_over_E = sqrt(1./(1.+pow(tan(lambda), 2)));
  else pt_over_E = pT/E;

  return (pt_over_E*pz/pow(cos(lambda), 2));
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


