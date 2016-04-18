#ifndef HMassConstraint_h
#define HMassConstraint_h

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cmath>
#include <map>
#include "TString.h"
#include "TLorentzVector.h"
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


using namespace std;
using namespace reco;
using namespace pat;

class HMassConstraint {
public:

   // Member functions
   HMassConstraint(RooSpin::VdecayType Vdecay1_, RooSpin::VdecayType Vdecay2_); // See RooSpin.h for what enumerators are available.
   void reset(); // { if(fit_res!=0) delete fit_res; }
   void addLeptons(const vector<map<Candidate*, PFParticle*>>& FermionFSR); // To set the Lepton and photon arrays in a map form. Pass null-pointer if the photon does not exist.
   RooFitResult* fit(); // Does the fit and returns a pointer to the fitresult.

   // Data members
   RooSpin* pdf;
   SpinPdfFactory* pdfFactory;

private:

   RooSpin::modelMeasurables measurables;
   RooFitResult* fitRes;

   vector<map<Candidate*,TMatrixDSym>> fermCovMatrix; // Maximum 4 x (3x3)
   vector<map<PFParticle*,TMatrixDSym>> gamCovMatrix; // Maximum 4 x (3x3)

   TMatrixDSym fitCovMatrix;
   RooGaussian* momCovConstraints;
};

#endif


