// PairPotential.cpp: implementation for the PairPotential class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "./PairPotential.h"

namespace Potential
{
namespace FinnisSinclair
{
namespace PairPotential
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}

Data::~Data() {}

Data::Data(const double * const rhs) 
	: c(rhs[0]), c0(rhs[1]), c1(rhs[2]), c2(rhs[3]), 
          fnd(rhs[4]), B(rhs[5]), alpha(rhs[6]) {}

Data::Data(
	const double &rhs_c, 
	const double &rhs_c0,
	const double &rhs_c1,
	const double &rhs_c2,
        const double &rhs_fnd,
        const double &rhs_B,
        const double &rhs_alpha) 
	: c(rhs_c), c0(rhs_c0), c1(rhs_c1), c2(rhs_c2), 
          fnd(rhs_fnd), B(rhs_B), alpha(rhs_alpha) {}

Data::Data(const Data &rhs) 
	: c(rhs.c), c0(rhs.c0), c1(rhs.c1), c2(rhs.c2), 
          fnd(rhs.fnd), B(rhs.B), alpha(rhs.alpha) {}

Data & Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this; 
	c = rhs.c; c0 = rhs.c0; c1 = rhs.c1; c2 = rhs.c2; 
        fnd = rhs.fnd; B = rhs.B; alpha = rhs.alpha; 
	return *this;
}

void Data::Randomize()
{
	c  = (double)rand()/(double)RAND_MAX;
	c0 = (double)rand()/(double)RAND_MAX;
	c1 = (double)rand()/(double)RAND_MAX;
	c2 = (double)rand()/(double)RAND_MAX;
        fnd = (double)rand()/(double)RAND_MAX;
        B = (double)rand()/(double)RAND_MAX;
        alpha = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState() {Prop = 0;}

LocalState::~LocalState() {}

Potential::Radial::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(const Data *rhs_Prop)
	: Prop(rhs_Prop) {}

LocalState::LocalState(const LocalState &) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	Prop = rhs.Prop; 
	return *this;
}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy() {LS=0;}

Energy<0>::~Energy(){}

Potential::Radial::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS){}

const Energy<0> & 
Energy<0>::operator = (const Energy<0> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<0>::operator () (const double &r) const
{
	double V = 0.0;
	if (r < LS->Prop->c) 
        {
		V = pow(r-LS->Prop->c,2)*
		(LS->Prop->c0 + LS->Prop->c1*r + LS->Prop->c2*r*r);
                if (r < LS->Prop->fnd)
                    V += LS->Prop->B*pow(LS->Prop->fnd-r,3)*
                         exp(-LS->Prop->alpha*r);
        }
	return 0.5*V;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {LS=0;}

Energy<1>::~Energy(){}

Potential::Radial::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS){}

const Energy<1> & 
Energy<1>::operator = (const Energy<1> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<1>::operator () (const double &r) const
{
	double dV = 0.0;
	if (r < LS->Prop->c)
        {    
		dV = 2.0*(r-LS->Prop->c)*
		(LS->Prop->c0 + LS->Prop->c1*r + LS->Prop->c2*r*r)
		+ pow(r-LS->Prop->c,2)*(LS->Prop->c1 + 2.0*LS->Prop->c2*r);
                if (r < LS->Prop->fnd)
                    dV -= 3.0*LS->Prop->B*pow(LS->Prop->fnd-r,2)*exp(-LS->Prop->alpha*r) +
                          LS->Prop->B*LS->Prop->alpha*pow(LS->Prop->fnd-r,3)*exp(-LS->Prop->alpha*r);
        }
	return 0.5*dV;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {LS=0;}

Energy<2>::~Energy(){}

Potential::Radial::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS){}

const Energy<2> & 
Energy<2>::operator = (const Energy<2> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<2>::operator () (const double &r) const
{
	double ddV = 0.0;
        if ( fabs(r-LS->Prop->c) < 1e-5) std::cout << "WARNING, Energy<2>: r->c" << std::endl;
	if (r < LS->Prop->c)
        {
		ddV = 2.0*(LS->Prop->c0 + LS->Prop->c1*r + LS->Prop->c2*r*r)
			+ 4.0*(r-LS->Prop->c)*(LS->Prop->c1 + 2.0*LS->Prop->c2*r)
			+ 2.0*LS->Prop->c2*pow(r-LS->Prop->c,2);
                if (r < LS->Prop->fnd)
                    ddV += LS->Prop->B*exp(-LS->Prop->alpha*r)*(LS->Prop->fnd-r) * 
                           ( 6.0 + 6.0*LS->Prop->alpha*(LS->Prop->fnd-r) +
                             pow(LS->Prop->alpha*(LS->Prop->fnd-r),2) );
        }
	return 0.5*ddV;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {LS=0;}

Jet<0>::~Jet() {}

Potential::Radial::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

Jet<0> & 
Jet<0>::operator = (const Jet<0> &rhs)
{
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
}

pair<double,double>
Jet<0>::operator () (const double &r) const
{
	double V = 0.0, DV = 0.0;
	if (r < LS->Prop->c) 
	{
		double A = pow(r-LS->Prop->c,2);
		double B = LS->Prop->c0 + LS->Prop->c1*r + LS->Prop->c2*r*r;
		double DA = 2.0*(r-LS->Prop->c);
		double DB = LS->Prop->c1 + 2.0*LS->Prop->c2*r;
		V = A*B; DV = DA*B + A*DB;
                if (r < LS->Prop->fnd) 
                {    
                   V += LS->Prop->B*pow(LS->Prop->fnd-r,3)*
                         exp(-LS->Prop->alpha*r);
                   DV -= 3.0*LS->Prop->B*pow(LS->Prop->fnd-r,2)*exp(-LS->Prop->alpha*r) +
                         LS->Prop->B*LS->Prop->alpha*pow(LS->Prop->fnd-r,3)*exp(-LS->Prop->alpha*r);
                }  
	}
	return make_pair(0.5*V,0.5*DV);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {LS=0;}

Jet<1>::~Jet() {}

Potential::Radial::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

Jet<1> & 
Jet<1>::operator = (const Jet<1> &rhs)
{
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
}

pair<double,double> 
Jet<1>::operator () (const double &r) const
{
	double DV = 0.0, DDV=0.0;
	if (r < LS->Prop->c) 
	{
		double A = pow(r-LS->Prop->c,2);
		double B = LS->Prop->c0 + LS->Prop->c1*r + LS->Prop->c2*r*r;
		double DA = 2.0*(r-LS->Prop->c);
		double DB = LS->Prop->c1 + 2.0*LS->Prop->c2*r;
		double DDB = 2.0*LS->Prop->c2*pow(r-LS->Prop->c,2);
		DV = DA*B + A*DB;
		DDV = 2.0*B + 2.0*DA*DB + DDB;
                if (r < LS->Prop->fnd)
                {    
                    DV -= 3.0*LS->Prop->B*pow(LS->Prop->fnd-r,2)*exp(-LS->Prop->alpha*r) +
                          LS->Prop->B*LS->Prop->alpha*pow(LS->Prop->fnd-r,3)*exp(-LS->Prop->alpha*r);
                    DDV += LS->Prop->B*exp(-LS->Prop->alpha*r)*(LS->Prop->fnd-r) * 
                           ( 6.0 + 6.0*LS->Prop->alpha*(LS->Prop->fnd-r) +
                             pow(LS->Prop->alpha*(LS->Prop->fnd-r),2) );
                }    
	}
	return make_pair(0.5*DV,0.5*DDV);
}

}

}

}
