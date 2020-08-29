// Hausdorff.cpp: implementation of the Hausdorff class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <map>
#include <algorithm>
#include <cassert>
#include "./Hausdorff.h"

namespace Potential
{
  namespace Field
  {
    namespace Hausdorff
    {
      //////////////////////////////////////////////////////////////////////
      // Class Data
      //////////////////////////////////////////////////////////////////////

      Data::Data():isDeformable(false)  {}
  
      Data::Data(vector<set<Set::Manifold::Point*> > * BC_,	  // connectivity of each face  
		 map<Set::Manifold::Point*, Set::VectorSpace::Vector> * yemb_, // coordinates of all the nodes
		 vector<Set::VectorSpace::Vector> * N_, // normal direction of each face
		 const double & C_, // sprint stiffness
		 bool  isDeformable_):
	BC(BC_), yemb(yemb_), N(N_), 
	C(C_), isDeformable(isDeformable_),
	xmin(N_->size()), xmax(N_->size()),
	ymin(N_->size()), ymax(N_->size()),
	isUpdated(false) {
	assert( C >= 0.0 );
      }
  
      Data::~Data(){}

      Data::Data(const Data &rhs) :
	y(rhs.y), N(rhs.N), C(rhs.C), yemb(rhs.yemb), 
	isDeformable(rhs.isDeformable), BC(rhs.BC),
	xmin(rhs.xmin), xmax(rhs.xmax), 
	ymin(rhs.ymin), ymax(rhs.ymax),
	isUpdated(rhs.isUpdated) {}

      Data &
      Data::operator = (const Data &rhs) {
	if (this == &rhs) return *this;
	y = rhs.y; N = rhs.N; C = rhs.C;
	yemb = rhs.yemb; BC = rhs.BC;
	isDeformable = rhs.isDeformable;
	xmin = rhs.xmin; xmax = rhs.xmax; 
	ymin = rhs.ymin; ymax = rhs.ymax;
	isUpdated = rhs.isUpdated;
	return *this;
      }

      void
      Data::Initialize() {
	assert(!yemb->empty());
	for (int i = 0; i < BC->size(); ++i ) {
	  y.push_back(&yemb->find(*(*BC)[i].begin())->second);
	}

	_calculateBounds();

	return;
      }

      void
      Data::ResetBounds() {
	if ( !isUpdated && isDeformable ) {
	  _calculateBounds();
	  isUpdated = true;
	}

	return;
      }

      void 
      Data::_calculateBounds() {
	assert(!yemb->empty());

	if ( (*N)[0].size() == 2 ) {
	  set<Set::Manifold::Point*>::const_iterator pE;
	  for ( int i = 0; i < BC->size(); ++i ) {
	    pE = (*BC)[i].begin();	    
	    const Set::VectorSpace::Vector & ya = yemb->find(*pE)->second;
	    pE++;
	    const Set::VectorSpace::Vector & yb = yemb->find(*pE)->second;
	    if ( ya[0] > yb[0] ) {
	      xmax[i] = ya[0];
	      xmin[i] = yb[0];
	    }
	    else {
	      xmax[i] = yb[0];
	      xmin[i] = ya[0];
	    }
	    
	    if ( ya[1] > yb[1] ) {
	      ymax[i] = ya[1];
	      ymin[i] = yb[1];
	    }
	    else {
	      ymax[i] = yb[1];
	      ymin[i] = ya[1];
	    }	   
	  }
	}
	else {
	  cerr << "_calculateBounds should not be called in dimensions other than 2" << endl;
	}
	
	return;
      }

      bool & 
      Data::IsDeformable() {
	return isDeformable;
      }

      bool &
      Data::IsUpdated() {
	return isUpdated;
      }      

      bool
      Data::IsDeformable() const {
	return isDeformable;
      }

      bool
      Data::IsUpdated() const {
	return isUpdated;
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

      LocalState::LocalState(Data *Prop_) : 
	Prop(Prop_)  {
      }

      LocalState::LocalState(const LocalState &rhs) : 
	Prop(rhs.Prop) {
      }

      double 
      LocalState::_det(const Set::VectorSpace::Vector & u, 
		       const Set::VectorSpace::Vector & v) const {
	return u[0]*v[1] - u[1]*v[0];
      }

      bool
      LocalState::FindNearestSurface(const Set::VectorSpace::Vector &x,
				     int & index, double & xn, 
				     Set::VectorSpace::Vector & dx) const {
	const vector<Set::VectorSpace::Vector> & Nloc = *Prop->N;
	const vector<Set::VectorSpace::Vector*> & yloc = Prop->y;

	bool hasNoContact = true;
	double min_xn = 1.0e16;
	for ( int i = 0; i < Nloc.size(); ++i ) {
	  Set::VectorSpace::Vector dxloc = x - *yloc[i];
	  double xnloc = dxloc(Nloc[i]);
	  if ( xnloc < 0.0 ) hasNoContact = false;
	  if ( fabs(xnloc) < min_xn ) {
	    min_xn = fabs(xnloc);
	    index = i;
	    dx = dxloc;
	    xn = xnloc;
	  }
	}

	return hasNoContact;
      }

      bool 
      LocalState::IsOutofBound(int index, double xn, const Set::VectorSpace::Vector & dx) const {

	Set::VectorSpace::Vector xt = dx - xn*(*Prop->N)[index];
	xt += *(Prop->y[index]);

	int dim = dx.size();
	if ( dim == 2 ) {	 
	  if ( xt[0] > Prop->xmax[index] ||
	       xt[0] < Prop->xmin[index] ||
	       xt[1] > Prop->ymax[index] ||
	       xt[1] < Prop->ymin[index] ) return true;	  
	}
	else if ( dim == 3 ) {
	  const set<Set::Manifold::Point*> & surface = (*Prop->BC)[index];
	  assert(surface.size() == 3);
	  set<Set::Manifold::Point*>::const_iterator pS = surface.begin();
	  const Set::VectorSpace::Vector & v0 = Prop->yemb->find(*pS)->second;

	  pS++;
	  Set::VectorSpace::Vector v1 = Prop->yemb->find(*pS)->second - v0;

	  pS++;
	  Set::VectorSpace::Vector v2 = Prop->yemb->find(*pS)->second - v0;

	  double v1v2 = _det(v1, v2);
	  double a = (_det(xt, v2) - _det(v0, v2))/v1v2;
	  double b = (_det(xt, v1) - _det(v0, v1))/v1v2;
	  if ( a > 0.0 && b > 0.0 && a+b < 1.0 ) {
	    return false;
	  }
	  else {
	    return true;
	  }
	}
	else {
	  return true;
	}
    
      }
  
      void 
      LocalState::operator ++ () {
	Prop->ResetBounds();
      }

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
	int index = 0;
	double xn;
	Set::VectorSpace::Vector dx(x.size());

	if ( LS->FindNearestSurface(x, index, xn, dx) ) {
	  return 0.0;
	}
	else {
	  if ( xn >= 0.0 ) {
	    return 0.0;
	  }
	  else {
	    if ( LS->IsOutofBound(index, xn, dx) ) {
	      return 0.0;
	    }
	    else {
	      return 0.5 * LS->Prop->C * xn * xn;
	    }
	  }
	}
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
	int index = 0;
	double xn;
	Set::VectorSpace::Vector dx(x.size());
	Set::VectorSpace::Vector f(x.size());	
	if ( LS->FindNearestSurface(x, index, xn, dx) ) {
	  return f;
	}
	else {
	  if ( xn >= 0.0 ) {
	    return f;
	  }
	  else {
	    if ( LS->IsOutofBound(index, xn, dx) ) {
	      return f;
	    }
	    else {
	      return LS->Prop->C * xn * (*LS->Prop->N)[index];
	    }
	  }
	}
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
	int index = 0;
	double xn;
	Set::VectorSpace::Vector dx(x.size());
	Set::VectorSpace::Hom df(x.size());

	if ( LS->FindNearestSurface(x, index, xn, dx) ) {
	  return df;
	}
	else {
	  if ( xn >= 0.0 ) {
	    return df;
	  }
	  else {
	    if ( LS->IsOutofBound(index, xn, dx) ) {
	      return df;
	    }
	    else {
	      const Set::VectorSpace::Vector & Nloc = (*LS->Prop->N)[index];
	      return Dyadic(LS->Prop->C*Nloc, Nloc);
	    }
	  }
	}
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
	int index = 0;
	double xn;
	Set::VectorSpace::Vector dx(x.size());
	Set::VectorSpace::Vector f(x.size());	
	if ( LS->FindNearestSurface(x, index, xn, dx) ) {
	  return make_pair(0.0, f);
	}
	else {
	  if ( xn >= 0.0 ) {
	    return make_pair(0.0, f);
	  }
	  else {
	    if ( LS->IsOutofBound(index, xn, dx) ) {
	      return make_pair(0.0, f);
	    }
	    else {
	      double fn = LS->Prop->C*xn;
	      return make_pair(0.5*fn*xn, fn*(*LS->Prop->N)[index]);
	    }
	  }
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
	int index = 0;
	double xn;
	Set::VectorSpace::Vector dx(x.size());
	Set::VectorSpace::Vector f(x.size());
	Set::VectorSpace::Hom df(x.size());
	
	if ( LS->FindNearestSurface(x, index, xn, dx) ) {
	  return make_pair(f, df);
	}
	else {
	  if ( xn >= 0.0 ) {
	    return make_pair(f, df);
	  }
	  else {
	    if ( LS->IsOutofBound(index, xn, dx) ) {
	      return make_pair(f, df);
	    }
	    else {
	      const Set::VectorSpace::Vector & Nloc = (*LS->Prop->N)[index];
	      f = LS->Prop->C*xn*Nloc;
	      df = Dyadic(LS->Prop->C*Nloc, Nloc);
	      return make_pair(f,df);
	    }
	  }
	}
      }

    }// namespace Hausdorff
    
  }//namespace Field
  
}//namespace Potential
