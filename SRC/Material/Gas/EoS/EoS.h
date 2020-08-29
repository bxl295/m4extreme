// EoS.h: interface for the EoS class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_LAW_GAS_EOS__INCLUDED_)
#define MATERIAL_LAW_GAS_EOS__INCLUDED_

#pragma once

#include <utility>
#include "../../Material.h"

using namespace std;

namespace Material
{
namespace Gas
{
namespace EoS
{
class Data;
class LocalState;
template<unsigned int> class Energy;
template<unsigned int> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

class Data : public Material::Data
{
public: 

	Data(){}
	virtual ~Data(){}
};


//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState
{
public: 

	LocalState(){}
	virtual ~LocalState(){}
	virtual LocalState *Clone() const=0;
	virtual void operator ++ ()=0;
	virtual void Reset(double T) {}
        virtual double GetTemperature() const { return 0.0; }
        double GetPressure() const { return _pressure; }
protected:
        double _pressure;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy
{
public: 

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<p> *Clone() const=0;
	virtual double operator () (const double &) const=0;
	virtual double operator () (double, double) const {return 0.0;}
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet
{
public: 

	Jet(){}
	virtual ~Jet(){}
	virtual Jet<p> *Clone() const=0;
	virtual pair<double,double> operator () (const double &) const=0;
        virtual pair<double,double> operator () (double, double) const {assert(false);}
};

}

}

}

#endif // !defined(MATERIAL_LAW_GAS_EOS__INCLUDED_)
