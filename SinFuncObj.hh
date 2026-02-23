#ifndef SINFUNCOBJ_HH
#define SINFUNCOBJ_HH

#include "Math/IFunction.h"
#include <math.h>
#include "Math/ChebyshevApprox.h"
#include "gsl/gsl_chebyshev.h"
#include "TH1D.h"

class SinFuncObj: public ROOT::Math::IBaseFunctionOneDim {
public:
  ROOT::Math::IBaseFunctionOneDim* Clone() const { return new SinFuncObj(); }  
private:
  double DoEval( double x ) const {
    return sin(x);
  }
};


class TH1DFuncObj: public ROOT::Math::IBaseFunctionOneDim {
public:
  TH1DFuncObj() : hist(0) {}
  TH1DFuncObj( TH1D* h ) : hist(h) {}
  ROOT::Math::IBaseFunctionOneDim* Clone() const { return new TH1DFuncObj(); }  
private:
  double DoEval( double x ) const {
    return hist->Interpolate( x );
  }
  TH1D* hist;
};


#endif
