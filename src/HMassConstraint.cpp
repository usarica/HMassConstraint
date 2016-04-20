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

  setJECUserFloatString();

  constructVariables();
  constructPdfFactory();
  constructPdfConstraint();
  constructCompoundPdf();
}

~HMassConstraint::HMassConstraint(){
  destroyCompoundPdf();
  destroyPdfConstraint();
  destroyPdfFactory();
  destroyVariables();
}


void HMassConstraint::constructVariables(){
  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;

  varZero = new RooRealVar("varZero", "", 0.);
  varOne = new RooRealVar("varOne", "", 1.);

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
          invcov_ferm[iZ][iferm][3*ix+iy] = new RooRealVar(Form("%s_vs_%s_Z%iFermion%i", title_ix.Data(), title_iy.Data(), iZ+1, iferm+1), "", -1e9, 1e9);
          invcov_ferm[iZ][iferm][3*ix+iy]->removeMin();
          invcov_ferm[iZ][iferm][3*ix+iy]->removeMax();
          invcov_fsr[iZ][iferm][3*ix+iy] = new RooRealVar(Form("%s_vs_%s_Z%iFermion%iFSR", title_ix.Data(), title_iy.Data(), iZ+1, iferm+1), "", -1e9, 1e9);
          invcov_fsr[iZ][iferm][3*ix+iy]->removeMin();
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
      sumdiffproducts_ferm_fsr[iZ][iferm] = new RooFormulaVar(Form("sumdiffproducts_Z%iFermion%i", iZ+1, iferm+1), "(@0+@1+@2+@3+@4+@5+@6+@7+@8+@9+@10+@11+@12+@13+@14+@15+@16+@17)/2.", sum_product_ferm_fsr_args);
    }

    // Construct m1/m2
    m[iZ] = new RooFormulaVar(Form("m%iRefit", iZ+1), "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mHdaughter_args);

    // This beta should multiply the spinPDF bc. having massive fermions have additional scale factors. These factors are even more relevant when FSR is present!
    beta_Vdaughter[iZ] = new RooFormulaVar(Form("betaV%iRefit", iZ+1), "(@0>0. ? sqrt( TMath::Max(1e-15, ( 1.-pow((@1+@2)/@0,2) )*( 1.-pow((@1-@2)/@0,2) ) ) ) : 1.)", RooArgList(*(m[iZ]), *(m_Hdaughter[iZ][0]), *(m_Hdaughter[iZ][1])));
  }
  // Construct m12
  m[2] = new RooFormulaVar("m12Refit", "sqrt( pow(@0+@4+@8+@12,2)-pow(@1+@5+@9+@13,2)-pow(@2+@6+@10+@14,2)-pow(@3+@7+@11+@15,2) )", m12_args);

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
    xvvFactory = new TensorPdfFactory_XVV(measurables, Vdecay1, Vdecay2, true); // true for always-on-shell X
    xvvFactory->makeParamsConst(false); // So that we can play with couplings
    spinPDF = xvvFactory->getPDF();
    pdfFactory = xvvFactory;
  }
}
void HMassConstraint::constructConstraintPdfs(){
  gausConstraintsPDF = new RooExponential("gausConstraintsPDF", "gausConstraintsPDF", sumdiffproducts_ferm_fsr_combined, varOne);
  auxilliaryConstraintsPDF = new RooGenericPdf("auxilliaryConstraintsPDF", "@0*@1", RooArgList(*(beta_Vdaughter[0]), *(beta_Vdaughter[1]))); // Will need to add m1, m2 cuts here as well!
  constraintsPDF = new RooProdPdf("constraintsPDF", "constraintsPDF", RooArgList(*gausConstraintsPDF, *auxilliaryConstraintsPDF));
}
void HMassConstraint::constructCompoundPdf(){
  RooArgList pdfList(*spinPDF, *constraintsPDF);
  PDF = RooProdPdf("HMassConstraint_PDF", "HMassConstraint_PDF", pdfList);
}

void HMassConstraint::destroyVariables(){
  // Destroy sums of inv. cov.*dxidxj
  delete sumdiffproducts_ferm_fsr_combined;

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

  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;

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

void HMassConstraint::addDaughters(const std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR){ // Candidate supports jets as well! FSR is also a reco::Candidate daughter.
  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;
  const Double_t empty_matrix_coefficients[3*3] = { 0. };

  int ndaughters=0;
  // Initialize PDG id's and bar-momenta
  for (int ix=0; ix<2; ix++){ for (int iy=0; iy<2; iy++) pdgid_ferm[ix][iy]=pdgUnknown; }

  for (std::vector::iterator<pair<reco::Candidate*, pat::PFParticle*>> dau=FermionWithFSR.begin(); dau<FermionWithFSR.end(); dau++){
    const reco::Candidate* fermion = dau.first();
    if (fermion==0){ cerr << "HMassConstraint::addDaughters : Daughter " << ndaughters << " pointer is 0!" << endl; break; }
    else{
      ndaughters++;
      if (ndaughter>4){ cerr << "HMassConstraint::addDaughters : Number of daughters (" << ndaughters << ") exceeds 4!" << endl; break; }
      else{
        // ndaughters=1..4
        int iZ = (ndaughters>2 ? 1 : 0);
        int iferm = 1-(ndaughters%2); // 0, 1
        int pdgid = fermion->pdgId();
        if ((iferm==0 && pdgid<0) || (iferm==1 && pdgid>0)) iferm=1-iferm;

        // Set PDG id
        pdgid_ferm[iZ][iferm] = fermion->pdgId();
        if (!(abs(pdgid_ferm[iZ][iferm])==pdgEle || abs(pdgid_ferm[iZ][iferm])==pdgMu || abs(pdgid_ferm[iZ][iferm])==pdgTau)) pdgid_ferm[iZ][iferm]=pdgJet; // Needs to be more thorough if jets are passsed

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
          lambda__ferm[iZ][iferm]->setRange(-lambdacut_electron, lambdacut_electron);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgMuon){
          pT_ferm[iZ][iferm]->setRange(pTcut_muon, sqrts);
          lambda__ferm[iZ][iferm]->setRange(-lambdacut_muon, lambdacut_muon);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgJet){
          pT_ferm[iZ][iferm]->setRange(pTcut_jet, sqrts);
          lambda__ferm[iZ][iferm]->setRange(-lambdacut_jet, lambdacut_jet);
        }
        phi_ferm[iZ][iferm]->setRange(-pi_val, pi_val);
        // Initialize refit fermion momenta to the values of fermion bar-momenta
        pT_ferm[iZ][iferm]->setVal(pTbar_ferm[iZ][iferm]->getVal());
        lambda_ferm[iZ][iferm]->setVal(lambdabar_ferm[iZ][iferm]->getVal());
        phi_ferm[iZ][iferm]->setVal(phibar_ferm[iZ][iferm]->getVal());

        // Get fermion covariance matrices in terms of pT, lambda and phi
        invCovMatrix_ferm[iZ][iferm] = TMatrixDSym(empty_matrix_coefficients);
        TMatrixDSym covMatrix_ferm = getCovarianceMatrix(fermion);
        Double_t determinant = 0;
        invCovMatrix_ferm[iZ][iferm] = covMatrix_ferm.Inverse(&determinant);
        if(determinant==0. && covMatrix_ferm[0][0]*covMatrix_ferm[1][1]*covMatrix_ferm[2][2]!=0.){ // This means matrix inversion failed. Remove off-diagonal terms and retry.
          for(int ix=0;ix<3;ix++){ for(int iy=0;iy<3;iy++){ if(ix!=iy) covMatrix_ferm[ix][iy]=0.; } }
          invCovMatrix_ferm[iZ][iferm] = covMatrix_ferm.Inverse(&determinant);
        }
        else invCovMatrix_ferm[iZ][iferm] = TMatrixDSym(empty_matrix_coefficients);

        // Do FSR here
        const pat::PFParticle* gamma = dau.second();
        if (gamma!=0){
          pTbar_fsr[iZ][iferm]->setConstant(false);
          lambdabar_fsr[iZ][iferm]->setConstant(false);
          phibar_fsr[iZ][iferm]->setConstant(false);
          pTbar_fsr[iZ][iferm]->setVal(gamma->pt());
          lambdabar_fsr[iZ][iferm]->setVal(piovertwo_val-gamma->theta());
          phibar_fsr[iZ][iferm]->setVal(gamma->phi());
          pTbar_fsr[iZ][iferm]->setConstant(true);
          lambdabar_fsr[iZ][iferm]->setConstant(true);
          phibar_fsr[iZ][iferm]->setConstant(true);

          // Set fsr ranges within the cuts and initialize
          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., sqrts); pT_fsr[iZ][iferm]->setVal(pTbar_fsr[iZ][iferm]->getVal());
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-lambdacut_fsr, lambdacut_fsr); lambda_fsr[iZ][iferm]->setVal(lambdabar_fsr[iZ][iferm]->getVal());
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); phi_fsr[iZ][iferm]->setVal(phibar_fsr[iZ][iferm]->getVal());

          // Get fsr covariance matrices in terms of pT, lambda and phi
          TMatrixDSym covMatrix_fsr = getCovarianceMatrix(gamma);
          determinant = 0;
          invCovMatrix_fsr[iZ][iferm] = covMatrix_fsr.Inverse(&determinant);
          if(determinant==0. && covMatrix_fsr[0][0]*covMatrix_fsr[1][1]*covMatrix_fsr[2][2]!=0.){ // This means matrix inversion failed. Remove off-diagonal terms and retry.
            for(int ix=0;ix<3;ix++){ for(int iy=0;iy<3;iy++){ if(ix!=iy) covMatrix_fsr[ix][iy]=0.; } }
            invCovMatrix_fsr[iZ][iferm] = covMatrix_fsr.Inverse(&determinant);
          }
          else invCovMatrix_fsr[iZ][iferm] = TMatrixDSym(empty_matrix_coefficients);
        }
        else{
          pTbar_fsr[iZ][iferm]->setConstant(false);
          lambdabar_fsr[iZ][iferm]->setConstant(false);
          phibar_fsr[iZ][iferm]->setConstant(false);
          pTbar_fsr[iZ][iferm]->setVal(0.);
          lambdabar_fsr[iZ][iferm]->setVal(0.);
          phibar_fsr[iZ][iferm]->setVal(0.);
          pTbar_fsr[iZ][iferm]->setConstant(true);
          lambdabar_fsr[iZ][iferm]->setConstant(true);
          phibar_fsr[iZ][iferm]->setConstant(true);

          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setConstant(true);
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setConstant(true);
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setConstant(true);

          invCovMatrix_fsr[iZ][iferm] = TMatrixDSym(empty_matrix_coefficients);
        }

      }
    }
  }

  // Set those non-existing particles
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      if (pdgid_ferm[iZ][iferm]!=pdgUnknown) continue;

      massbar_ferm[iZ][iferm]->setConstant(false); massbar_ferm[iZ][iferm]->setVal(0.); massbar_ferm[iZ][iferm]->setConstant(true);
      pTbar_ferm[iZ][iferm]->setConstant(false); pTbar_ferm[iZ][iferm]->setVal(0.); pTbar_ferm[iZ][iferm]->setConstant(true);
      lambdabar_ferm[iZ][iferm]->setConstant(false); lambdabar_ferm[iZ][iferm]->setVal(0.); lambdabar_ferm[iZ][iferm]->setConstant(true);
      phibar_ferm[iZ][iferm]->setConstant(false); phibar_ferm[iZ][iferm]->setVal(0.); phibar_ferm[iZ][iferm]->setConstant(true);
      pTbar_fsr[iZ][iferm]->setConstant(false); pTbar_fsr[iZ][iferm]->setVal(0.); pTbar_fsr[iZ][iferm]->setConstant(true);
      lambdabar_fsr[iZ][iferm]->setConstant(false); lambdabar_fsr[iZ][iferm]->setVal(0.); lambdabar_fsr[iZ][iferm]->setConstant(true);
      phibar_fsr[iZ][iferm]->setConstant(false); phibar_fsr[iZ][iferm]->setVal(0.); phibar_fsr[iZ][iferm]->setConstant(true);

      pT_ferm[iZ][iferm]->setConstant(false); pT_ferm[iZ][iferm]->setRange(0., 0.); pT_ferm[iZ][iferm]->setVal(0.); pT_ferm[iZ][iferm]->setConstant(true);
      lambda_ferm[iZ][iferm]->setConstant(false); lambda_ferm[iZ][iferm]->setRange(0., 0.); lambda_ferm[iZ][iferm]->setVal(0.); lambda_ferm[iZ][iferm]->setConstant(true);
      phi_ferm[iZ][iferm]->setConstant(false); phi_ferm[iZ][iferm]->setRange(0., 0.); phi_ferm[iZ][iferm]->setVal(0.); phi_ferm[iZ][iferm]->setConstant(true);
      pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setConstant(true);
      lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setConstant(true);
      phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setConstant(true);

      invCovMatrix_ferm[iZ][iferm] = TMatrixDSym(empty_matrix_coefficients);
      invCovMatrix_fsr[iZ][iferm] = TMatrixDSym(empty_matrix_coefficients);
    }
  }

  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      // Propagate the inverse covariance matrix elements to the RooFormulaVars
      for (int ix=0; ix<3; ix++){
        for (int iy=0; iy<3; iy++){
          invcov_ferm[iZ][iferm][3*ix+iy]->setConstant(false); invcov_ferm[iZ][iferm][3*ix+iy]->setVal((invCovMatrix_ferm[iZ][iferm])[ix][iy]); invcov_ferm[iZ][iferm][3*ix+iy]->setConstant(true);
          invcov_fsr[iZ][iferm][3*ix+iy]->setConstant(false); invcov_fsr[iZ][iferm][3*ix+iy]->setVal((invCovMatrix_fsr[iZ][iferm])[ix][iy]); invcov_fsr[iZ][iferm][3*ix+iy]->setConstant(true);
        }
      }
    }
  }

}

void HMassConstraint::fit(){
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

  RooDataSet data("fittedHdaughters","",data_args);
  data.add(data_args);

  RooArgSet conditionals;
  conditionals.add(*(m[2]));
  fitResult = PDF->fitTo(data, ConditionalObservables(conditionals), Save(true), RooFit::PrintLevel(-1));
  fitCovMatrix = fitResult->covarianceMatrix();

  cout << "Number of columns in the covariance mtrix is " << fitCovMatrix->GetNCols() << endl;

}

TMatrixDSym HMassConstraint::getCovarianceMatrix(const reco::Candidate* particle){
  TMatrixDSym covMatrix_empty(3);

  const reco::GsfElectron* electron = dynamic_cast<const reco::GsfElectron*>(fermion);
  const pat::Muon* muon = dynamic_cast<const pat::Muon*>(fermion);
  const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(fermion);

  if (electron!=0) return getCovarianceMatrix(electron);
  else if (muon!=0) return getCovarianceMatrix(muon);
  else if (pfcand!=0) return getCovarianceMatrix(pfcand);
  else return covMatrix_empty;
}
TMatrixDSym HMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle){
  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;

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

  double C_pt_pt = pow(energyerr/(particle->pt()/particle->energy()*(1.+pow(tan(lambda), 2))), 2);
  double C_lambda_lambda = 0;
  if (fabs(lambda)>1e-5) C_lambda_lambda = pow(energyerr/(cos(lambda)*(1.+pow(tan(lambda), 2))), 2);
  else C_lambda_lambda = pow(energyerr/(cos(1e-5)*(1.+pow(tan(1e-5*sign(1., lambda)), 2))), 2);

  // Everything else is 0. I know this cannot be correct, but let's work with it for now.

  double momCov[9]={ 0 };
  momCov[3*0+0] = C_pt_pt;
  momCov[3*1+1] = C_lambda_lambda;
  TMatrixDSym cov(3, momCov);
  return cov;
}
TMatrixDSym HMassConstraint::getCovarianceMatrix(const pat::Muon* particle){
  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;

  double lambda = piovertwo_val - particle->theta();
  double pterr = particle->userFloat("correctedPtError");
  double pterr_uncorrected = particle->muonBestTrack()->ptError();
  double correction = pterr/pterr_uncorrected;

  double trackCov[TrackBase::dimension*TrackBase::dimension];
  for (int ix=0; ix<TrackBase::dimension; ix++){
    for (int iy=0; iy<TrackBase::dimension; iy++){
      trackCov[5*ix+iy] = particle->muonBestTrack()->covariance(ix, iy);
      if ((ix==TrackBase::i_qoverp || ix==TrackBase::i_lambda) && (iy==TrackBase::i_qoverp || iy==TrackBase::i_lambda)) trackCov[5*ix+iy] *= correction;
    }
  }

  double d_pT_d_qoverp = -particle->pt()*particle->p()/particle->charge();
  double d_pT_d_lambda = -particle->pz();
  const double d_pT_d_phi = 0;
  double d_lambda_d_qoverp;
  if (fabs(lambda)>1e-5) d_lambda_d_qoverp = -particle->p()/particle->charge()/tan(lambda);
  else d_lambda_d_qoverp = -particle->p()/particle->charge()/tan(1e-5*sign(1., lambda));
  const double d_lambda_d_lambda = 1;
  const double d_lambda_d_phi = 0;
  const double d_phi_d_qoverp = 0;
  const double d_phi_d_lambda = 0;
  const double d_phi_d_phi = 1;



  double momCov[9]={ 0 };
  momCov[3*0+0] = // pT, pT
    d_pT_d_qoverp*d_pT_d_qoverp * trackCov[5*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_pT_d_lambda * trackCov[5*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_pT_d_phi * trackCov[5*TrackBase::i_phi + TrackBase::i_phi] +
    2.*d_pT_d_qoverp*d_pT_d_lambda * trackCov[5*TrackBase::i_qoverp + TrackBase::i_lambda] +
    2.*d_pT_d_qoverp*d_pT_d_phi * trackCov[5*TrackBase::i_qoverp + TrackBase::i_phi] +
    2.*d_pT_d_lambda*d_pT_d_phi * trackCov[5*TrackBase::i_lambda + TrackBase::i_phi]
    ;
  momCov[3*0+1] = // pT, lambda
    d_pT_d_qoverp*d_lambda_d_qoverp * trackCov[5*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_lambda_d_lambda * trackCov[5*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_lambda_d_phi * trackCov[5*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_lambda_d_lambda + d_pT_d_lambda*d_lambda_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_lambda) * trackCov[5*TrackBase::i_lambda + TrackBase::i_phi]
    ;
  momCov[3*0+2] = // pT, phi
    d_pT_d_qoverp*d_phi_d_qoverp * trackCov[5*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_phi_d_lambda * trackCov[5*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_phi_d_phi * trackCov[5*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[5*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  momCov[3*1+0] = momCov[3*0+1];// lambda, pT
  momCov[3*1+1] = // lambda, lambda
    d_lambda_d_qoverp*d_lambda_d_qoverp * trackCov[5*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_lambda_d_lambda*d_lambda_d_lambda * trackCov[5*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_lambda_d_phi*d_lambda_d_phi * trackCov[5*TrackBase::i_phi + TrackBase::i_phi] +
    (d_lambda_d_qoverp*d_lambda_d_lambda + d_lambda_d_lambda*d_lambda_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_lambda_d_qoverp*d_lambda_d_phi + d_lambda_d_phi*d_lambda_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_lambda_d_lambda*d_lambda_d_phi + d_lambda_d_phi*d_lambda_d_lambda) * trackCov[5*TrackBase::i_lambda + TrackBase::i_phi]
    ;
  momCov[3*1+2] = // lambda, phi
    d_pT_d_qoverp*d_phi_d_qoverp * trackCov[5*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_phi_d_lambda * trackCov[5*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_phi_d_phi * trackCov[5*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[5*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  momCov[3*2+0] = momCov[3*0+2];// phi, pT
  momCov[3*2+1] = momCov[3*1+2];// phi, lambda
  momCov[3*2+2] = // phi, phi
    d_phi_d_qoverp*d_phi_d_qoverp * trackCov[5*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_phi_d_lambda*d_phi_d_lambda * trackCov[5*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_phi_d_phi*d_phi_d_phi * trackCov[5*TrackBase::i_phi + TrackBase::i_phi] +
    (d_phi_d_qoverp*d_phi_d_lambda + d_phi_d_lambda*d_phi_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_phi_d_qoverp*d_phi_d_phi + d_phi_d_phi*d_phi_d_qoverp) * trackCov[5*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_phi_d_lambda*d_phi_d_phi + d_phi_d_phi*d_phi_d_lambda) * trackCov[5*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  TMatrixDSym cov(3, momCov);
  return cov;
s}
TMatrixDSym HMassConstraint::getCovarianceMatrix(const reco::PFCandidate* particle){
  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;

  double energyerr = PFEnergyResolution().getEnergyResolutionEm(particle->energy(), particle->eta());
  double lambda = piovertwo_val - particle->theta();

  double C_pt_pt = pow(energyerr/(particle->pt()/particle->energy()*(1.+pow(tan(lambda), 2))), 2);
  double C_lambda_lambda = 0;
  if (fabs(lambda)>1e-5) C_lambda_lambda = pow(energyerr/(cos(lambda)*(1.+pow(tan(lambda), 2))), 2);
  else C_lambda_lambda = pow(energyerr/(cos(1e-5)*(1.+pow(tan(1e-5*sign(1.,lambda)), 2))), 2);

  // Everything else is 0. I know this cannot be correct, but let's work with it for now.

  double momCov[9]={ 0 };
  momCov[3*0+0] = C_pt_pt;
  momCov[3*1+1] = C_lambda_lambda;
  TMatrixDSym cov(3, momCov);
  return cov;
}
TMatrixDSym HMassConstraint::getCovarianceMatrix(const pat::Jet* particle){
  const Double_t pi_val = 3.14159265358979323846;//TMath::Pi();
  const Double_t piovertwo_val = pi_val/2.;

  double energyerr = particle->userFloat(jecString.Data());
  double lambda = piovertwo_val - particle->theta();

  double C_pt_pt = pow(energyerr/(particle->pt()/particle->energy()*(1.+pow(tan(lambda), 2))), 2);
  double C_lambda_lambda = 0;
  if (fabs(lambda)>1e-5) C_lambda_lambda = pow(energyerr/(cos(lambda)*(1.+pow(tan(lambda), 2))), 2);
  else C_lambda_lambda = pow(energyerr/(cos(1e-5)*(1.+pow(tan(1e-5*sign(1., lambda)), 2))), 2);

  // Everything else is 0. I know this cannot be correct, but let's work with it for now.

  double momCov[9]={ 0 };
  momCov[3*0+0] = C_pt_pt;
  momCov[3*1+1] = C_lambda_lambda;
  TMatrixDSym cov(3, momCov);
  return cov;
}



