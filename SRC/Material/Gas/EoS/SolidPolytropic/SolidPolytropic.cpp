// SolidPolytropic.cpp: implementation for the SolidPolytropic class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./SolidPolytropic.h"

namespace Material
{
namespace Gas
{
namespace EoS
{
namespace SolidPolytropic
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data(): C(1.0),Gamma(1.0) {}
	
Data::Data(const double &C_,const double &Gamma_,const double &Rho0_) 
: C(C_),Gamma(Gamma_),Rho0(Rho0_){}
	
Data::Data(const Data &rhs) : C(rhs.C),Gamma(rhs.Gamma),Rho0(rhs.Rho0){}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	C = rhs.C; 
        Gamma = rhs.Gamma;
	Rho0 = rhs.Rho0;
	return *this;
}

const double &
Data::GetC() const
{
	return C;
}

double & 
Data::GetC()
{
	return C;
}

const double &
Data::GetGamma() const
{
	return Gamma;
}

double & 
Data::GetGamma()
{
	return Gamma;
}

const double &
Data::GetRho0() const
{
	return Rho0;
}

double & 
Data::GetRho0()
{
	return Rho0;
}

void 
Data::Randomize()
{
	C = (double)rand()/(double)RAND_MAX;
	Gamma = (double)rand()/(double)RAND_MAX;
	Rho0 = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::Gas::EoS::LocalState *
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

Energy<0>::~Energy(){}

Material::Gas::EoS::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_){}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS){}

double 
Energy<0>::operator () (const double &J) const
{
	const double &C = LS->Prop->C;
	const double &Gamma = LS->Prop->Gamma;
	const double &Rho0 = LS->Prop->Rho0;
	double G = Gamma-1.0;	
        return C * pow(Rho0, Gamma) * ( pow(1.0/fabs(J), G) / G + fabs(J) - Gamma / G);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy(){}

Material::Gas::EoS::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_) : LS(LS_){}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS){}

double
Energy<1>::operator () (const double &J) const
{
	const double &C = LS->Prop->C;
	const double &Gamma = LS->Prop->Gamma;
	const double &Rho0 = LS->Prop->Rho0;
	
        double p = -C * pow(Rho0, Gamma) * ( pow(1.0/fabs(J), Gamma) - 1.0 );

	LS->_pressure = p;
        
	return p;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy(){}

Material::Gas::EoS::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_) : LS(LS_){}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS){}

double 
Energy<2>::operator () (const double &J) const
{
	const double &C = LS->Prop->C;
	const double &Gamma = LS->Prop->Gamma;
	const double &Rho0 = LS->Prop->Rho0;
	return (C*Gamma/Rho0)*pow((Rho0/J),Gamma+1.0);  
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Material::Gas::EoS::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_) : LS(LS_){}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS){}

pair<double,double>
Jet<0>::operator () (const double &J) const 
{
	const double &C = LS->Prop->C;
	const double &Gamma = LS->Prop->Gamma;
	const double &Rho0G = pow(LS->Prop->Rho0, Gamma);
	double G = Gamma-1.0;	
        return make_pair( C * Rho0G * ( pow(1.0/fabs(J), G) / G + fabs(J) - Gamma / G),
			  -C * Rho0G * (pow(1.0/fabs(J), Gamma) - 1.0) );
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Material::Gas::EoS::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_){}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS){}

pair<double,double>
Jet<1>::operator () (const double &J) const
{
	const double &C = LS->Prop->C;
	const double &Gamma = LS->Prop->Gamma;
	const double &Rho0G = pow(LS->Prop->Rho0, Gamma);
	double G = Gamma-1.0;	
        return make_pair( -C * Rho0G * ( pow(1.0/fabs(J), Gamma) - 1.0 ), 
			  C * Rho0G * Gamma* pow(1.0/fabs(J), Gamma+1.0) );
}

}

}

}

}
