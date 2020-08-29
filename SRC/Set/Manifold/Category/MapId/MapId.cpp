// MapId.cpp: Implementation of the MapId class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "MapId.h"

namespace Set
{
namespace Manifold
{
namespace Category
{
//////////////////////////////////////////////////////////////////////
// MapId<0> class
//////////////////////////////////////////////////////////////////////

MapId<0>::MapId() : n(0){}

MapId<0>::MapId(const unsigned int &rhs_n) : n(rhs_n) {}

MapId<0>::~MapId(){}

MapId<0>::MapId(const MapId<0> &rhs) : n(rhs.n) {}

MapId<0> & 
MapId<0>::operator = (const MapId<0> &rhs)
{
	if (this == &rhs) return *this;
	n = rhs.n; return *this;
}

Set::Manifold::Map * 
MapId<0>::Clone()
{
	return new MapId<0>(*this);
}
	
Set::Manifold::TMap *
MapId<0>::Diff()
{
	return new MapId<1>(*this);
}

void 
MapId<0>::Randomize(){}

const Point & 
MapId<0>::operator () (const Point &x){return x;}

unsigned int 
MapId<0>::size1() const{return n;}

unsigned int 
MapId<0>::size2() const{return n;}

//////////////////////////////////////////////////////////////////////
// MapId<1> class
//////////////////////////////////////////////////////////////////////

MapId<1>::MapId() : n(0){}

MapId<1>::MapId(const unsigned int &rhs_n) : n(rhs_n), I(n) {}

MapId<1>::~MapId(){}

MapId<1>::MapId(const MapId<0> &rhs) : n(rhs.size1()), I(rhs.size1()) {}

MapId<1>::MapId(const MapId<1> &rhs) : n(rhs.n), I(n) {}

MapId<1> & 
MapId<1>::operator = (const MapId<1> &rhs)
{
	if (this == &rhs) return *this;
	n = rhs.n; return *this;
}

Set::Manifold::TMap * 
MapId<1>::Clone()
{
	return new MapId<1>(*this);
}
	
Set::Manifold::TMap *
MapId<1>::Diff()
{
	return new MapId<2>(*this);
}

void 
MapId<1>::Randomize(){}

const Set::VectorSpace::Hom & 
MapId<1>::operator () (const Point &){return I;}

unsigned int 
MapId<1>::size1() const{return n;}

unsigned int 
MapId<1>::size2() const{return n;}

//////////////////////////////////////////////////////////////////////
// MapId<2> class
//////////////////////////////////////////////////////////////////////

MapId<2>::MapId() : n(0){}

MapId<2>::MapId(const unsigned int &rhs_n) : n(rhs_n), O(n,n*n) {}

MapId<2>::~MapId(){}

MapId<2>::MapId(const MapId<1> &rhs) 
: n(rhs.size1()), O(rhs.size1(),rhs.size1()*rhs.size2()) {}

MapId<2>::MapId(const MapId<2> &rhs) : n(rhs.n), O(n,n*n) {}

MapId<2> & 
MapId<2>::operator = (const MapId<2> &rhs)
{
	if (this == &rhs) return *this;
	n = rhs.n; return *this;
}

Set::Manifold::TMap * 
MapId<2>::Clone()
{
	return new MapId<2>(*this);
}

void 
MapId<2>::Randomize(){}

const Set::VectorSpace::Hom & 
MapId<2>::operator () (const Point &){return O;}

unsigned int 
MapId<2>::size1() const{return n;}

unsigned int 
MapId<2>::size2() const{return n;}

}

}

}
