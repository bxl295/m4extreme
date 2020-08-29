// HalfSpace.cpp: implementation of the HalfSpace class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <map>
#include <algorithm>
#include <cassert>
#include "./HalfSpace.h"

namespace Potential
{
namespace Field
{
namespace HalfSpace
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data():Finite(false)  {}

Data::Data(const Set::VectorSpace::Vector &y_, 
	   const Set::VectorSpace::Vector &N_, 
	   const double &C_) : y(y_), N(N_), C(C_), Finite(false)
{
	assert(C >= 0.0); 
	assert(Norm(N) > 0);
	assert(y.size() == N.size());
	N /= Norm(N);
}

// for 2-dimensional problem
Data::Data(const Set::VectorSpace::Vector &y_, 
	   const Set::VectorSpace::Vector &yend_,
	   const Set::VectorSpace::Vector &N_, 
	   const double &C_) : y(y_), N(N_), C(C_), Finite(true)
{
	assert(C >= 0.0); 
	assert(Norm(N) > 0);
	assert(y.size() == N.size());
	N /= Norm(N);       
	y0_max = y_[0] > yend_[0] ? y_[0] : yend_[0];
	y0_min = y_[0] > yend_[0] ? yend_[0] : y_[0];
	y1_max = y_[1] > yend_[1] ? y_[1] : yend_[1];
	y1_min = y_[1] > yend_[1] ? yend_[1] : y_[1];
}

// for 3-dimensional problem (Note: only works for convex polygon in 3-dimensional Euclidean space)
Data::Data(const vector<Set::VectorSpace::Vector *> & V_, 	     
	   const Set::VectorSpace::Vector & N_, 
	   const double & C_): N(N_), C(C_), y(3), anchor(3), Finite(true)
{
  assert(C >= 0.0); 
  assert(N_.size() == 3 && Norm(N) > 0);

  const int numofVertices = V_.size();
  assert( numofVertices > 2 );

  // calculate the barycenter of the convex polygon
  for ( int i = 0; i < numofVertices; ++i ) {
    y += *V_[i];
  }

  y /= (double)numofVertices;

  // calculate the rays pointing from the barycenter to the vertices
  map<int, Set::VectorSpace::Vector> rays;
  map<int, Set::VectorSpace::Vector>::iterator pR;
  for ( int i = 0; i < numofVertices; ++i ) {
    rays.insert( make_pair(i, *V_[i] - y) );
  }

  anchor = rays[0];

  // sort the vertices in clockwise and store the wedges and angles
  double fval;
  int start_id = 0;
  int count = 0;
  int idx[numofVertices];
  idx[0] = start_id;
  const double * pN = N_.begin();
  double w[3];

  while ( count < numofVertices-1 ) {   

    map<double, int> tloc;

    Set::VectorSpace::Vector u(rays[start_id]);
    double norm_u = Norm(u);
    rays.erase(start_id);
    const double * pu = u.begin();

    for ( pR = rays.begin(); pR != rays.end(); pR++ ) {
      const Set::VectorSpace::Vector & v = pR->second;
      const double * pv = v.begin();

      // cross product of u and v
      w[0] = *(pu+1) * *(pv+2) - *(pu+2) * *(pv+1);
      w[1] = *(pu+2) * *(pv) - *(pu) * *(pv+2);
      w[2] = *(pu) * *(pv+1) - *(pu+1) * *(pv);

      // determine the cosine of the angle between u and v
      fval = u(v);
      fval /= norm_u;
      fval /= Norm(v);
      fval = fval >  1.0 ?  1.0 : fval;
      fval = fval < -1.0 ? -1.0 : fval;

      // based on the orientation to calculate the true angle
      if ( *pN * w[0] + *(pN+1) * w[1] + *(pN+2) * w[2] > 0. ) {
	tloc.insert( make_pair(6.28318530718 - acos(fval), pR->first) );    
      }
      else {
	tloc.insert( make_pair(acos(fval), pR->first) );
      }
    }
    
    theta.push_back( (tloc.begin())->first );
    start_id = (tloc.begin())->second;
    idx[++count] = start_id;
  }

  for ( int i = 1; i < numofVertices-1; ++i ) {
    theta[i] += theta[i-1];
  }

  for ( int i = 0; i < numofVertices; ++i ) {
    vertices.push_back(*V_[idx[i]]);
  }

}

Data::~Data(){}

Data::Data(const Data &rhs) : y(rhs.y), N(rhs.N), C(rhs.C) {}

Data &
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	y = rhs.y; N = rhs.N; C = rhs.C; 
	return *this;
}

void 
Data::Randomize()
{
	C = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Potential::Field::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Data *Prop_) : Prop(Prop_) {}

LocalState::LocalState(const LocalState &rhs) : Prop(rhs.Prop) {}

bool 
LocalState::IsOutofBound(double xn, const Set::VectorSpace::Vector & dx) const {
  if (Prop->Finite) {
    Set::VectorSpace::Vector xt = dx - xn*Prop->N;
    int dim = xt.size();
    
    if ( dim == 2 ) {
      xt += Prop->y;
      if ( xt[0] > Prop->y0_max ||
	   xt[0] < Prop->y0_min ||
	   xt[1] > Prop->y1_max ||
	   xt[1] < Prop->y1_min ) return true;
    }
    else if ( dim == 3 ) {
      double norm_xt = Norm(xt);
      if ( norm_xt == 0. ) return false;

      const Set::VectorSpace::Vector & u = Prop->anchor;
      const vector<double> & theta = Prop->theta;
      double w[3];
      const double * pu = u.begin();
      const double * pxt = xt.begin();

      // cross product of u and xt
      w[0] = *(pu+1) * *(pxt+2) - *(pu+2) * *(pxt+1);
      w[1] = *(pu+2) * *(pxt) - *(pu) * *(pxt+2);
      w[2] = *(pu) * *(pxt+1) - *(pu+1) * *(pxt);

      // find the angle between u and xt
      double fval = u(xt);
      fval /= Norm(u);
      fval /= norm_xt;
      fval = fval >  1.0 ?  1.0 : fval;
      fval = fval < -1.0 ? -1.0 : fval;

      vector<double>::const_iterator low;
      const double * pN = Prop->N.begin();

      if ( *pN * w[0] + *(pN+1) * w[1] + *(pN+2) * w[2] > 0. ) {
	low = std::lower_bound(theta.begin(), theta.end(), 6.28318530718 - acos(fval));
      }
      else {
	low = std::lower_bound(theta.begin(), theta.end(), acos(fval));
      }

      int start, end;
      if ( low != theta.end() ) {
	start = low - theta.begin();
	end = start + 1;
      }
      else {
	start = theta.size();
	end = 0;
      }

      Set::VectorSpace::Vector y_end   = Prop->vertices[end]   - Prop->y;
      Set::VectorSpace::Vector y_start = Prop->vertices[start] - Prop->y;
      Set::VectorSpace::Vector x_end   = y_end - xt;
      Set::VectorSpace::Vector x_start = y_start - xt;

      if( y_end(x_end) * y_start(x_start) < y_end(x_start) * y_start(x_end) ) return true;
    }
  }
  
  return false;
}
  
void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Potential::Field::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Vector dx = x - LS->Prop->y;
	double xn = LS->Prop->N(dx);
	if (xn >= 0.0 || LS->IsOutofBound(xn, dx)) return 0.0;
	else return 0.5*LS->Prop->C*xn*xn;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy() {}

Potential::Field::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Vector dx = x - LS->Prop->y;
	Set::VectorSpace::Vector f(x.size());
	double xn = LS->Prop->N(dx);
	if (xn >= 0.0 || LS->IsOutofBound(xn, dx) ) return f;
	else return LS->Prop->C*xn*LS->Prop->N;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy() {}

Potential::Field::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Vector dx = x - LS->Prop->y;
	double xn = LS->Prop->N(dx);
	Set::VectorSpace::Hom df(x.size());
	if (xn >= 0.0 || LS->IsOutofBound(xn, dx) ) return df;
	else return Dyadic(LS->Prop->C*LS->Prop->N,LS->Prop->N);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Potential::Field::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Vector dx = x - LS->Prop->y;
	double xn = LS->Prop->N(dx);
	Set::VectorSpace::Vector f(x.size());
	if (xn >= 0.0 || LS->IsOutofBound(xn, dx) ) return make_pair(0.0,f);
	else {
	  double fn = LS->Prop->C*xn;
	  return make_pair(0.5*fn*xn, fn*LS->Prop->N);
	}
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Potential::Field::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Vector dx = x - LS->Prop->y;
	double xn = LS->Prop->N(dx);
	Set::VectorSpace::Vector f(x.size());
	Set::VectorSpace::Hom df(x.size());
	if (xn >= 0.0 || LS->IsOutofBound(xn, dx) ) {
	  return make_pair(f,df);
	}
	else {
	  f = LS->Prop->C*xn*LS->Prop->N;
	  df = Dyadic(LS->Prop->C*LS->Prop->N,LS->Prop->N);
	  return make_pair(f,df);
	}

}

}

}

}
