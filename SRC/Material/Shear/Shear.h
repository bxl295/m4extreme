// Shear.h: interface for the Shear class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_SHEAR__INCLUDED_)
#define MATERIAL_SHEAR__INCLUDED_

#pragma once

#include <utility>
#include "../Material.h"
#include "../../Utils/Triplet/Triplet.h"

using namespace std;

namespace Material
{
namespace Shear
{
class LocalState;
template<unsigned int> class Modulus;
template<unsigned int> class Jet;
template<unsigned int> class JetJet;

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
};

//////////////////////////////////////////////////////////////////////
// Class Modulus<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Modulus
{
public: 

	Modulus(){}
	virtual ~Modulus(){}
	virtual Modulus<p> *Clone() const=0;
	virtual double operator () (const double &) const=0;
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
};

//////////////////////////////////////////////////////////////////////
// Class JetJet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class JetJet
{
public: 

	JetJet(){}
	virtual ~JetJet(){}
	virtual JetJet<p> *Clone() const=0;
	virtual triplet<double,double,double> operator () (const double &) const=0;
};

}

}

#endif // !defined(MATERIAL_SHEAR__INCLUDED_)
