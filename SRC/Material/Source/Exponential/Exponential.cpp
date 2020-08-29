// Exponential.cpp: implementation of the Exponential class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./Exponential.h"

namespace Material
{
namespace Source
{
namespace Exponential
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

  Data::Data() : _x(0), _y(0), _z(0), _A(0), _b(0), _n(0) {}

  Data::Data(double *x_, 
	     double *y_,
	     double *z_,
	     double *A_,
	     double *b_,
	     int n_) : 
    _x(x_), _y(y_), _z(z_), _A(A_), _b(b_), _n(n_) {}

Data::~Data(){}

Data::Data(const Data &rhs) {}

Data &
Data::operator = (const Data &rhs)
{
  if (this == &rhs) return *this;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::LocalState *
LocalState::Clone() const 
{
  return new LocalState(*this);
}

LocalState::LocalState(Data *Prop_) : Prop(Prop_) {}

LocalState::LocalState(const LocalState &rhs) : Prop(rhs.Prop) {}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
  return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &qp) const
{
  int m = qp.size();
  double w = 0.0;
  Data *D = LS->Prop;
  double * A = D->_A;
  double * x = D->_x;
  double * y = D->_y;
  double * z = D->_z;
  double * b = D->_b;
  double dx, dy, dz;
  double factor = 2.0;
  for ( int i = 1; i < m; ++i ) factor *= 2.0;

  for ( int i = 0; i < D->_n; ++i ) {
    double psi = 0.0;
    dx = qp[0] - x[i];
    dy = 0.;
    dz = 0.;

    if ( m == 1 ) {
      psi = -b[i]*dx*dx;
    }
    else if ( m == 2 ) {
      dy = qp[1] - y[i];
      psi = -b[i]*(dx*dx + dy*dy);
    }
    else {
      dy = qp[1] - y[i];
      dz = qp[2] - z[i];
      psi = -b[i]*(dx*dx + dy*dy + dz*dz);
    }

    w += A[i] * b[i] * exp(psi) * (4.0 * psi + factor);
  }

  return w;
}


}

}

}
