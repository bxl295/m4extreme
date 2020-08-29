// Cholesky.h: interface for the Cholesky class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER_LINEAR_CHOLESKY__INCLUDED_)
#define SOLVER_LINEAR_CHOLESKY__INCLUDED_

#pragma once

#include <map>
#include "../Linear.h"
#include "../../../Set/Algebraic/VectorSpace/Vector/Vector.h"
#include "../../../Set/Algebraic/VectorSpace/Category/Category.h"
#include "../../../Utils/LinearAlgebra/Cholesky/Cholesky.h"

using namespace std;

namespace Solver
{
namespace Linear
{
template <class key>
class Cholesky : public System <key> 
{
public:

//	Constructors/destructors:

	Cholesky();
	Cholesky(const set <key> &);
	virtual ~Cholesky();

//	Accessors/mutators:

	double Get(key);
	double Get(key, key);
	void Set(key, double);
	void Set(key, key, double);

//	General methods:

	void SetToZero1();
	void SetToZero2();
	void Add(key, double);
	void Add(key, key, double);
	void Solve();
	double Norm();

private:

	unsigned int n;

	map <key, unsigned int> KeyMap;

	Set::VectorSpace::Hom a;
	Set::VectorSpace::Vector b;

	LinearAlgebra::Cholesky LS;

private:

	Cholesky(Cholesky &);
	void operator=(Cholesky &);
};
template <class key>
Cholesky<key>::Cholesky() : n(0) {}

template <class key>
Cholesky<key>::Cholesky(const set <key> &Keys)
: n(Keys.size()), a(n), b(n), LS(n)
{
	unsigned int i;
	typename set <key>::const_iterator iK;
	for (iK=Keys.begin(), i=0; iK!=Keys.end(); iK++, i++) 
		KeyMap[*iK] = i;
}

template <class key>
Cholesky<key>::~Cholesky() {}

template <class key>
double Cholesky<key>::Get(key K)
{
	return b[KeyMap[K]];
}

template <class key>
double Cholesky<key>::Get(key K1, key K2)
{
	return a[KeyMap[K2]][KeyMap[K1]];
}

template <class key>
void Cholesky<key>::Set(key K, double Input)
{
	b[KeyMap[K]] = Input;
}

template <class key>
void Cholesky<key>::Set(key K1, key K2, double Input)
{
	a[KeyMap[K2]][KeyMap[K1]] = Input;
}

template <class key>
void Cholesky<key>::SetToZero1()
{
	Null(b);
}

template <class key>
void Cholesky<key>::SetToZero2()
{
	Null(a);
}

template <class key>
void Cholesky<key>::Add(key K, double Input)
{
	b[KeyMap[K]] += Input;
}

template <class key>
double Cholesky<key>::Norm()
{
	unsigned int i; double bn2 = 0.0;
	for (i=0; i<n; i++) bn2 += b[i]*b[i];
	return sqrt(bn2);
}

template <class key>
void Cholesky<key>::Add(key K1, key K2, double Input)
{
	a[KeyMap[K2]][KeyMap[K1]] += Input;
}

template <class key>
void Cholesky<key>::Solve()
{
	LS.Solve(a,b);
}

}

}

#endif // !defined(SOLVER_LINEAR_CHOLESKY__INCLUDED_)
