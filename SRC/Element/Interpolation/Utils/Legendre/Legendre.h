// Legendre.h: interface for the Legendre class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_LEGENDRE_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_LEGENDRE_H__INCLUDED_

#pragma once

#include <map>
#include <vector>

using namespace std;

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class Legendre<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Legendre;

//////////////////////////////////////////////////////////////////////
// Class Legendre<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Legendre<0>
{
friend class Legendre<1>;

public: 

	Legendre();
	Legendre(const unsigned int &);
	virtual ~Legendre();
	vector<double> operator () (const double &) const;

private:

	unsigned int n;

private:

	Legendre(const Legendre<0> &);
	Legendre<0> & operator = (const Legendre<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Legendre<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Legendre<1>
{
public: 

	Legendre();
	Legendre(const unsigned int &);
	virtual ~Legendre();
	vector<double> operator () (const double &) const;

private:

	unsigned int n;

private:

	Legendre(const Legendre<1> &);
	Legendre<1> & operator = (const Legendre<1> &);
};

//////////////////////////////////////////////////////////////////////
// Class AssociatedLegendre<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class AssociatedLegendre;

//////////////////////////////////////////////////////////////////////
// Class AssociatedLegendre<0>
//////////////////////////////////////////////////////////////////////

template<> 
class AssociatedLegendre<0>
{
friend class AssociatedLegendre<1>;

public: 

	AssociatedLegendre();
	AssociatedLegendre(const unsigned int &);
	virtual ~AssociatedLegendre();
	vector<map<int,double> > operator () (const double &) const;

private:

	unsigned int n;

private:

	AssociatedLegendre(const AssociatedLegendre<0> &);
	AssociatedLegendre<0> & operator = (const AssociatedLegendre<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class AssociatedLegendre<1>
//////////////////////////////////////////////////////////////////////

template<> 
class AssociatedLegendre<1>
{
public: 

	AssociatedLegendre();
	AssociatedLegendre(const unsigned int &);
	virtual ~AssociatedLegendre();
	vector<map<int,double> > operator () (const double &) const;

private:

	unsigned int n;

private:

	AssociatedLegendre(const AssociatedLegendre<1> &);
	AssociatedLegendre<1> & operator = (const AssociatedLegendre<1> &);
};

}

}

}

//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

unsigned int DoubleFactorial(const unsigned int &);

unsigned int Factorial(const unsigned int &);

#endif // !defined(ELEMENT_INTERPOLATION_LEGENDRE_H__INCLUDED_

