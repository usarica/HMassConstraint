#include <HMassConstraint/HMassConstraint/include/RooGaussianMomConstraint.h>
#include <cassert> 

using namespace std;


RooGaussianMomConstraint::RooGaussianMomConstraint(
  const char* name, const char* title,
  const RooArgList& variables_,
  const RooArgList& means_,
  const RooArgList& matrixElement_,
  RooGaussianMomConstraint::CoordinateSystem coordinates_,
  Int_t fixCode_
  ) : RooAbsPdf(name, title),

  variables("matrixElement", "List of inverse covariance matrix elements", this),
  means("matrixElement", "List of inverse covariance matrix elements", this),
  matrixElement("matrixElement", "List of inverse covariance matrix elements", this),
  coordinates(coordinates_)
{
  setProxyList(variables_, variables, 3);
  setProxyList(means_, means, 3);
  setProxyList(matrixElement_, matrixElement, 3*3);
  fixVariable(fixCode_);
}


RooGaussianMomConstraint::RooGaussianMomConstraint(const RooGaussianMomConstraint& other, const char* name) : RooAbsPdf(other, name),
variables(other.variables),
means(other.means),
matrixElement(other.matrixElement),
coordinates(other.coordinates),
fixCode(other.fixCode)
{}

void RooGaussianMomConstraint::fixVariable(Int_t code){
  fixCode=1;
  if ((code%prime_var1)==0)fixCode*=prime_var1;
  if ((code%prime_var2)==0)fixCode*=prime_var2;
  if ((code%prime_var3)==0)fixCode*=prime_var3;
}


void RooGaussianMomConstraint::setProxyList(const RooArgList& args, RooListProxy& target, Int_t checkDim){
  TIterator* argIter = args.createIterator();
  RooAbsArg* arg;
  Int_t nargs=0;
  while ((arg = (RooAbsArg*)argIter->Next())) {
    if (arg!=0){
      if (!dynamic_cast<RooRealVar*>(arg)) {
        coutE(InputArguments) << "ERROR::RooGaussianMomConstraint(" << GetName() << ")::setProxyList: Proxy " << arg->GetName() << " is not of type RooRealVar." << endl;
        assert(0);
      }
      target.add(*arg);
      nargs++;
    }
  }
  delete argIter;
  if (checkDim>0 && checkDim!=nargs){
    coutE(InputArguments) << "ERROR::RooGaussianMomConstraint(" << GetName() << ")::setProxyList: Proxies (size=" << nargs << ") do not have the required size " << checkDim << "." << endl;
    assert(0);
  }
  if (nargs==0){
    coutE(InputArguments) << "ERROR::RooGaussianMomConstraint(" << GetName() << ")::setProxyList: Number of proxy arguments are 0!" << endl;
    assert(0);
  }
}


Double_t RooGaussianMomConstraint::computeGaussian(const Int_t code) const{ // The code can be either for var1, 2 or 3, but it can only be for one of them.
  Int_t intVar=-1;
  if ((code%prime_var1)==0)intVar=0;
  else if ((code%prime_var2)==0)intVar=1;
  else if ((code%prime_var3)==0)intVar=2;


  //dynamic_cast<const RooRealVar*>(variables.at(intVar))->getVal()
  Double_t value=0;
  if (intVar<0){ // Just do regular computation
    const Double_t epsilon = 1e-15;

    for (int ix=0; ix<3; ix++){
      if (
        ((fixCode%prime_var1)==0 && ix==0) ||
        ((fixCode%prime_var2)==0 && ix==1) ||
        ((fixCode%prime_var3)==0 && ix==2)
        ) continue;

      Double_t val_ix = dynamic_cast<const RooRealVar*>(variables.at(ix))->getVal();
      Double_t valbar_ix = dynamic_cast<const RooRealVar*>(means.at(ix))->getVal();
      Double_t diff_ix=val_ix-valbar_ix;

      for (int iy=0; iy<3; iy++){
        if (
          ((fixCode%prime_var1)==0 && iy==0) ||
          ((fixCode%prime_var2)==0 && iy==1) ||
          ((fixCode%prime_var3)==0 && iy==2)
          ) continue;
        Double_t invCovMat = dynamic_cast<const RooRealVar*>(matrixElement.at(3*ix+iy))->getVal();
        if (invCovMat!=0.){
          Double_t val_iy = dynamic_cast<const RooRealVar*>(variables.at(iy))->getVal();
          Double_t valbar_iy = dynamic_cast<const RooRealVar*>(means.at(iy))->getVal();
          Double_t diff_iy=val_iy-valbar_iy;
          value += invCovMat*diff_ix*diff_iy;
        }
      }

    }
    if (!(value==value) || value<epsilon) return epsilon;
  }
  else{
    const Double_t epsilon = 1e-10;

    Double_t xmin = dynamic_cast<const RooRealVar*>(variables.at(intVar))->getMin();
    Double_t xmax = dynamic_cast<const RooRealVar*>(variables.at(intVar))->getMax();

    Double_t y = dynamic_cast<const RooRealVar*>(variables.at((intVar+1)%3))->getVal();
    Double_t z = dynamic_cast<const RooRealVar*>(variables.at((intVar+2)%3))->getVal();
    Double_t xb = dynamic_cast<const RooRealVar*>(means.at(intVar))->getVal();
    Double_t yb = dynamic_cast<const RooRealVar*>(means.at((intVar+1)%3))->getVal();
    Double_t zb = dynamic_cast<const RooRealVar*>(means.at((intVar+2)%3))->getVal();

    Double_t dxmax = xmax-xb;
    Double_t dxmin = xmin-xb;
    Double_t dy = y-yb;
    Double_t dz = z-zb;

    Double_t invCovMat[3][3]={ { 0 } };

    for (int ix=0; ix<3; ix++){
      Int_t rx = ((intVar+ix)%3);
      if (
        ((fixCode%prime_var1)==0 && rx==0) ||
        ((fixCode%prime_var2)==0 && rx==1) ||
        ((fixCode%prime_var3)==0 && rx==2)
        ) continue;
      for (int iy=0; iy<3; iy++){
        Int_t ry = ((intVar+iy)%3);
        if (
          ((fixCode%prime_var1)==0 && ry==0) ||
          ((fixCode%prime_var2)==0 && ry==1) ||
          ((fixCode%prime_var3)==0 && ry==2)
          ) continue;
        invCovMat[rx][ry] = dynamic_cast<const RooRealVar*>(matrixElement.at(3*rx+ry))->getVal();
      }
    }

    Double_t offxdiff = ((invCovMat[0][1]+invCovMat[1][0])*dy + (invCovMat[0][2]+invCovMat[2][0])*dz)/2.;
    Double_t diffothers = invCovMat[1][1]*pow(dy, 2) + (invCovMat[1][2]+invCovMat[2][1])*dy*dz + invCovMat[2][2]*pow(dz, 2);
    const Double_t tolerance = epsilon*0.01;
    if (fabs(invCovMat[0][0])>tolerance){
      Double_t exponent = 0.5*(
        pow(offxdiff, 2)/invCovMat[0][0]
        - diffothers
        );

      value = 1./sqrt(invCovMat[0][0]) * sqrt(Pi()/2.) * exp(exponent) * (
        Erf(2.*(invCovMat[0][0]*dxmax + offxdiff)/(2.*sqrt(2.)*sqrt(invCovMat[0][0])))
        - Erf(2.*(invCovMat[0][0]*dxmin + offxdiff)/(2.*sqrt(2.)*sqrt(invCovMat[0][0])))
        );
    }
    else if (fabs(offxdiff)>tolerance){
      value =
        exp(
        -0.5*diffothers
        )*(
        exp(-offxdiff*dxmin) - exp(-offxdiff*dxmax)
        )/offxdiff;
    }
    else{
      value =
        exp(
        -0.5*diffothers
        )*(xmax-xmin);
    }

    if (!(value==value) || value<epsilon) return epsilon;
  }
  return value;
}


Double_t RooGaussianMomConstraint::evaluate() const{
  return computeGaussian(1);
}

Double_t RooGaussianMomConstraint::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
  return computeGaussian(code);
}

Int_t RooGaussianMomConstraint::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = 1;
  if ((fixCode%prime_var1)!=0 && matchArgs(allVars, analVars, RooArgSet(*(variables.at(0))))) code *= prime_var1;
  else if ((fixCode%prime_var2)!=0 && matchArgs(allVars, analVars, RooArgSet(*(variables.at(1))))) code *= prime_var2;
  else if ((fixCode%prime_var3)!=0 && matchArgs(allVars, analVars, RooArgSet(*(variables.at(2))))) code *= prime_var3;
  if (code==1) code=0;
  return code;
}
