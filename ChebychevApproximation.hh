#ifndef CHEBYCHEVAPPROXIMATION_HH
#define CHEBYCHEVAPPROXIMATION_HH

#include <vector>
#include "TH1D.h"
#include <cmath>
#include <iostream>
#include <memory>

// Basic Chebychev polynomial
double chebychevPolynomial( double x, size_t nn ) {
  double result= 0.0;
  if( nn == 0 ) {
    result= 1.0;
  }
  else if( nn == 1 ) {
    result= x;
  }
  else {
    double tnminusone= 1.0;
    double tn= x;
    double tnplusone= -99.0;
    for( size_t i= 2; i < nn+1; i++ ) {
      tnplusone= 2.0*x*tn - tnminusone;
      tnminusone= tn;
      tn= tnplusone;
    }
    result= tnplusone;
  }
  return result;
}
// Transformations [-1,1] <-> [a,b]
double tTox( double t, double a, double b ) {
  return (b-a)/2.0*t + (a+b)/2.0;
}
double xTot( double x, double a, double b ) {
  return (x-(a+b)/2.0)/(b-a)*2.0;
}

// Helper classes to evaluate input and output objects
class InputEvaluator {
public:
  virtual double evaluate( double x ) const = 0;
};
// For TH1D objects
class InputEvaluatorTH1D : public InputEvaluator {
  TH1D* hist;
public:
  InputEvaluatorTH1D( TH1D* h ) : hist(h) {}
  virtual double evaluate( double x ) const {
    return hist->Interpolate( x );
  }
};
// For TH1F objects
class InputEvaluatorTH1F : public InputEvaluator {
  TH1F* hist;
public:
  InputEvaluatorTH1F( TH1F* h ) : hist(h) {}
  virtual double evaluate( double x ) const {
    return hist->Interpolate( x );
  }
};
// Unmodified coefficients
class OutputEvaluator {
public:
  virtual double evaluate( double x, double* pars, double a, double b, size_t n ) const {
    double value= 0.0;
    for( size_t j= 0; j < n; j++ ) {
      value+= pars[j]*chebychevPolynomial( xTot( x, a, b ), j );
    }
    return value;    
  }
};
// Normalised coefficients
class OutputEvaluatorN : public OutputEvaluator {
public:
  virtual double evaluate( double x, double* pars, double a, double b, size_t n ) const {
    double value= 1.0;
    for( size_t j= 1; j < n; j++ ) {
      value+= pars[j]*chebychevPolynomial( xTot( x, a, b ), j );
    }
    value*= pars[0];
    return value;    
  }
};

// Chebychev approximation
class ChebychevApproximation {
  double a;
  double b;
  size_t n;
  double* coeffs= nullptr;
  std::shared_ptr<InputEvaluator> ie;
  std::shared_ptr<OutputEvaluator> oe;
  void initialiseCoefficients( bool lnorm ) {
    coeffs= new double[n];
    createOutputEvaluator( lnorm );
    // Interpolation nodes on [-1,1] and [a,b]
    std::vector<double> tnodes;
    std::vector<double> xnodes;
    for( size_t k= 1; k < n+1; k++ ) {
      double tnode= cos( (double(k)-0.5)/double(n)*3.14159 );
      double xnode= tTox( tnode, a , b );
      tnodes.push_back( tnode );
      xnodes.push_back( xnode );
    }
    // Coefficients c0 to cn:
    for( size_t j= 0; j < n; j++ ) {
      double coeff= 0.0;
      for( size_t k= 0; k < n; k++ ) {
	double xnode= xnodes[k];
	double tnode= tnodes[k];
	coeff+= ie->evaluate( xnode )*chebychevPolynomial( tnode, j );
      }
      coeff*= 2.0/double(n);
      coeffs[j]= coeff;
    }
    coeffs[0]/= 2.0;
    if( lnorm ) normaliseCoefficients();
    // Check on nodes
    std::cout << "ChebychevApproximation::initialiseCoefficients: degree " << n
	      << " approximation on nodes"
	      << std::endl;
    for( size_t k= 0; k < n; k++ ) {
      double xnode= xnodes[k];
      std::cout << xnode << " " << ie->evaluate( xnode ) << " " << evaluate( xnode )
		<< std::endl;
    }
    return;
  }
  void normaliseCoefficients() {
    double norm= coeffs[0];
    coeffs[0]= norm;
    for( size_t i= 1; i < n; i++ ) { coeffs[i]/= norm; }
  }
  void createOutputEvaluator( bool lnorm ) {
    if( lnorm ) {
      oe= std::make_shared<OutputEvaluatorN>();
    }
    else {
      oe= std::make_shared<OutputEvaluator>();
    }
  }
public:
  ChebychevApproximation( TH1D* h, double aa, double bb, size_t nn, bool lnorm=false ) :
    a(aa), b(bb), n(nn) {
    ie= std::make_shared<InputEvaluatorTH1D>( h );
    initialiseCoefficients( lnorm );
  }
  ChebychevApproximation( TH1F* h, double aa, double bb, size_t nn, bool lnorm=false ) :
    a(aa), b(bb), n(nn) {
    ie= std::make_shared<InputEvaluatorTH1F>( h );
    initialiseCoefficients( lnorm );
  }
  ChebychevApproximation( const ChebychevApproximation & rhs ) {
    a= rhs.a;
    b= rhs.b;
    n= rhs.n;
    coeffs= new double[n];
    for( size_t i = 0; i < n; i++ ) { coeffs[i]= rhs.coeffs[i]; }
    ie= rhs.ie;
    oe= rhs.oe;
  }
  virtual ~ChebychevApproximation() {
    delete [] coeffs;
  }
  double evaluate( double x ) const {
    return this->operator()( &x, coeffs );
  }
  virtual double operator()( double* x, double* pars ) const {
    double xx= x[0];
    return oe->evaluate( xx, pars, a, b, n );
  }
  double* getCoefficients() const {
    return coeffs;
  }
};

#endif
