// Embedding.cpp: implementation for the Embedding class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "./Embedding.h"

namespace Potential
{
namespace Johnson
{
namespace Embedding
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////


Data::Data(){}

Data::~Data(){}

Data::Data(const double * const rhs) 
 : Ec(rhs[0]), alpha(rhs[1]), beta(rhs[2]), gamma(rhs[3]), fe(rhs[4]), phie(rhs[5]) {}

Data::Data(const double &rhs_Ec,
           const double &rhs_alpha,
	   const double &rhs_beta,
	   const double &rhs_gamma,
	   const double &rhs_fe,
	   const double &rhs_phie) 
 : Ec(rhs_Ec), alpha(rhs_alpha), beta(rhs_beta), gamma(rhs_gamma), fe(rhs_fe), phie(rhs_phie) {}
	
Data::Data(const Data &rhs) 
 : Ec(rhs.Ec), alpha(rhs.alpha), beta(rhs.beta), gamma(rhs.gamma), fe(rhs.fe), phie(rhs.phie) {}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	Ec=rhs.Ec; alpha=rhs.alpha; beta=rhs.beta; gamma=rhs.gamma; fe=rhs.fe; phie=rhs.phie;
	return *this;
}

void 
Data::Randomize()
{
	Ec    = (double)rand()/(double)RAND_MAX;
	alpha = (double)rand()/(double)RAND_MAX;
	beta  = (double)rand()/(double)RAND_MAX;
	gamma = (double)rand()/(double)RAND_MAX;
	fe    = (double)rand()/(double)RAND_MAX;
	phie  = (double)rand()/(double)RAND_MAX;
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
Energy<0>::operator () (const double &ro) const
{
	double F = 0.0;
	double phi = 6.0*LS->Prop->phie;
	double roe = 12.0*LS->Prop->fe;
	double aob = LS->Prop->alpha/LS->Prop->beta;
	double gob = LS->Prop->gamma/LS->Prop->beta;
	F = -LS->Prop->Ec*( 1.0 - aob*log(ro/roe) )*pow(ro/roe,aob)-
	    phi*pow(ro/roe,gob);
	return F;
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
Energy<1>::operator () (const double &ro) const
{
        if (ro < 1.0e-16) std::cout << "Error in Energy<1>: ro->0" << std::endl;
        
        double dF = 0.0;
	double phi = 6.0*LS->Prop->phie;
	double roe = 12.0*LS->Prop->fe;
	double aob = LS->Prop->alpha/LS->Prop->beta;
        double gob = LS->Prop->gamma/LS->Prop->beta;
	dF = LS->Prop->Ec*aob/ro*pow(ro/roe,aob) -
             phi*gob/roe*pow(ro/roe,gob-1.) -
	     LS->Prop->Ec*aob/roe*pow(ro/roe,aob-1.)*( 1.-aob*log(ro/roe) );
	return dF;
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
Energy<2>::operator () (const double &ro) const
{
        if (ro < 1.0e-16) std::cout << "Error in Energy<2>: ro->0" << std::endl;
        
        double ddF = 0.0;
	double phi = 6.0*LS->Prop->phie;
	double roe = 12.0*LS->Prop->fe;
	double aob = LS->Prop->alpha/LS->Prop->beta;
	double gob = LS->Prop->gamma/LS->Prop->beta;
	ddF = -gob*phi/pow(roe,2.)*(gob-1.)*pow(ro/roe,gob-2.) -
               LS->Prop->Ec*
	    	( aob/pow(ro,2.)*pow(ro/roe,aob) -
	    	  2.*pow(aob,2.)/(ro*roe)*pow(ro/roe,aob-1.) +
	    	  aob/pow(roe,2.)*(aob-1.)*pow(ro/roe,aob-2.)*(1.-aob*log(ro/roe)) );
	return ddF;
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
Jet<0>::operator () (const double &ro) const
{
        if (ro < 1.0e-16) std::cout << "Error in Jet<0>: ro->0" << std::endl;
      
	double F = 0.0, dF = 0.0;
	double phi = 6.0*LS->Prop->phie;
	double roe = 12.0*LS->Prop->fe;
	double aob = LS->Prop->alpha/LS->Prop->beta;
	double gob = LS->Prop->gamma/LS->Prop->beta;
	F = -LS->Prop->Ec*( 1.0 - aob*log(ro/roe) )*pow(ro/roe,aob)-
	     phi*pow(ro/roe,gob);
	dF = LS->Prop->Ec*aob/ro*pow(ro/roe,aob) -
	     phi*gob/roe*pow(ro/roe,gob-1.) -
	     LS->Prop->Ec*aob/roe*pow(ro/roe,aob-1.)*( 1.-aob*log(ro/roe) );			 
        return make_pair(F,dF);
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
Jet<1>::operator () (const double &ro) const
{
	if (ro < 1.0e-16) std::cout << "Error in Jet<1>: ro->0" << std::endl;
             
	double dF = 0.0, ddF = 0.0;
	double phi = 6.0*LS->Prop->phie;
	double roe = 12.0*LS->Prop->fe;
	double aob = LS->Prop->alpha/LS->Prop->beta;
	double gob = LS->Prop->gamma/LS->Prop->beta;
	dF = LS->Prop->Ec*aob/ro*pow(ro/roe,aob) -
             phi*gob/roe*pow(ro/roe,gob-1.) -
             LS->Prop->Ec*aob/roe*pow(ro/roe,aob-1.)*( 1.-aob*log(ro/roe) );
	ddF = -gob*phi/pow(roe,2.)*(gob-1.)*pow(ro/roe,gob-2.) -
               LS->Prop->Ec*
	    	( aob/pow(ro,2.)*pow(ro/roe,aob) -
	    	  2.*pow(aob,2.)/(ro*roe)*pow(ro/roe,aob-1.) +
	    	  aob/pow(roe,2.)*(aob-1.)*pow(ro/roe,aob-2.)*(1.-aob*log(ro/roe)) );
        return make_pair(dF,ddF);
}

}

}

}
