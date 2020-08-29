// Axisymmetric.cpp: Implementation of the Axisymmetric class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Axisymmetric.h"

namespace Set
{
namespace Euclidean
{	
namespace Axisymmetric
{
//////////////////////////////////////////////////////////////////////
// Class Embedding<0>
//////////////////////////////////////////////////////////////////////

Embedding<0>::Embedding() {}

Embedding<0>::Embedding(const Embedding<0> &f) : y(f.y) {}

Embedding<0>::~Embedding() {}

Embedding<0> & 
Embedding<0>::operator = (const Embedding<0> &f)
{
	if (this == &f) return *this; 
	y = f.y; return *this;
}

Set::Manifold::Map * 
Embedding<0>::Clone()
{
	return new Embedding<0>;
}
	
Set::Manifold::TMap *
Embedding<0>::Diff()
{
	return new Embedding<1>;
}

void 
Embedding<0>::Randomize() {}

const Set::Euclidean::Cylindrical::Point &
Embedding<0>::operator () (const Set::Euclidean::Orthonormal::Point &x)
{
	assert(x.size() == 2);
	y[0] = x[0]; 
	y[1] = 0.0;
	y[2] = x[1]; 
	return y;
}

unsigned int 
Embedding<0>::size1() const
{
	return 2;

}
unsigned int 
Embedding<0>::size2() const
{
	return 3;
}

const Set::Manifold::Point & 
Embedding<0>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Embedding<1>
//////////////////////////////////////////////////////////////////////

Embedding<1>::Embedding() : A(3,2) {}

Embedding<1>::Embedding(const Embedding<0> &f) : A(3,2) {}

Embedding<1>::Embedding(const Embedding<1> &f) : A(f.A) {}

Embedding<1>::~Embedding(){}

Embedding<1> & 
Embedding<1>::operator = (const Embedding<1> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Embedding<1>::Clone()
{
	return new Embedding<1>;
}
	
Set::Manifold::TMap *
Embedding<1>::Diff()
{
	return new Embedding<2>;
}

void 
Embedding<1>::Randomize() {}

const Set::VectorSpace::Hom &
Embedding<1>::operator () (const Set::Euclidean::Orthonormal::Point &x)
{
	assert(x.size() == 2);
	Null(A);
	A[0][0] = 1.0; 
	A[1][2] = 1.0; 
	return A;
}

unsigned int 
Embedding<1>::size1() const
{
	return 2;

}
unsigned int 
Embedding<1>::size2() const
{
	return 6;
}

const Set::VectorSpace::Hom & 
Embedding<1>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Embedding<2>
//////////////////////////////////////////////////////////////////////

Embedding<2>::Embedding() : A(6,2) {}

Embedding<2>::Embedding(const Embedding<1> &f) : A(6,2) {}

Embedding<2>::Embedding(const Embedding<2> &f) : A(f.A) {}

Embedding<2>::~Embedding(){}

Embedding<2> & 
Embedding<2>::operator = (const Embedding<2> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Embedding<2>::Clone()
{
	return new Embedding<2>;
}

void 
Embedding<2>::Randomize() {}

const Set::VectorSpace::HomZero &
Embedding<2>::operator () (const Set::Euclidean::Orthonormal::Point &x)
{
	assert(x.size() == 2);
	return A;
}

unsigned int 
Embedding<2>::size1() const
{
	return 2;

}
unsigned int 
Embedding<2>::size2() const
{
	return 12;
}

const Set::VectorSpace::Hom & 
Embedding<2>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<0>
//////////////////////////////////////////////////////////////////////

Submersion<0>::Submersion() : x(2) {}

Submersion<0>::Submersion(const Submersion<0> &f) : x(f.x) {}

Submersion<0>::~Submersion(){}

Submersion<0> & 
Submersion<0>::operator = (const Submersion<0> &f)
{
	if (this == &f) return *this;
	x = f.x; return *this;
}

Set::Manifold::Map * 
Submersion<0>::Clone()
{
	return new Submersion<0>(*this);
}
	
Set::Manifold::TMap *
Submersion<0>::Diff()
{
	return new Submersion<1>;
}

void 
Submersion<0>::Randomize() {}

const Set::Euclidean::Orthonormal::Point &
Submersion<0>::operator () (const Set::Euclidean::Cylindrical::Point &y)
{
	x[0] = y[0]; 
	x[1] = y[2]; 
	return x;
}

unsigned int 
Submersion<0>::size1() const
{
	return 3;

}
unsigned int 
Submersion<0>::size2() const
{
	return 2;
}

const Set::Manifold::Point & 
Submersion<0>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<1>
//////////////////////////////////////////////////////////////////////

Submersion<1>::Submersion() : A(2,3) {}

Submersion<1>::Submersion(const Submersion<0> &f) : A(2,3) {}

Submersion<1>::Submersion(const Submersion<1> &f) : A(f.A) {}

Submersion<1>::~Submersion(){}

Submersion<1> & 
Submersion<1>::operator = (const Submersion<1> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Submersion<1>::Clone()
{
	return new Submersion<1>;
}
	
Set::Manifold::TMap *
Submersion<1>::Diff()
{
	return new Submersion<2>;
}

void 
Submersion<1>::Randomize() {}

const Set::VectorSpace::Hom &
Submersion<1>::operator () (const Set::Euclidean::Cylindrical::Point &y)
{
	A[0][0] = 1.0; 
	A[2][1] = 1.0; 
	return A;
}

unsigned int 
Submersion<1>::size1() const
{
	return 3;

}
unsigned int 
Submersion<1>::size2() const
{
	return 6;
}

const Set::VectorSpace::Hom & 
Submersion<1>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Cylindrical::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<2>
//////////////////////////////////////////////////////////////////////

Submersion<2>::Submersion() : A(6,3) {}

Submersion<2>::Submersion(const Submersion<1> &f) : A(6,3) {}

Submersion<2>::Submersion(const Submersion<2> &f) : A(f.A) {}

Submersion<2>::~Submersion(){}

Submersion<2> & 
Submersion<2>::operator = (const Submersion<2> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Submersion<2>::Clone()
{
	return new Submersion<2>(*this);
}

void 
Submersion<2>::Randomize() {}

const Set::VectorSpace::HomZero &
Submersion<2>::operator () (const Set::Euclidean::Cylindrical::Point &y)
{
	return A;
}

unsigned int 
Submersion<2>::size1() const
{
	return 3;

}
unsigned int 
Submersion<2>::size2() const
{
	return 18;
}

const Set::VectorSpace::Hom & 
Submersion<2>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Cylindrical::Point &)P);
}

}

}

}
