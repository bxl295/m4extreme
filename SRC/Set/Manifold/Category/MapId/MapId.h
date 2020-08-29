// MapId.h: Interface for the MapId class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_MANIFOLD_CATEGORY_MAPID_H__INCLUDED_)
#define SET_MANIFOLD_CATEGORY_MAPID_H__INCLUDED_

#pragma once

#include "../Category.h"

namespace Set
{
namespace Manifold
{
namespace Category
{
//////////////////////////////////////////////////////////////////////
// MapId class
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class MapId;

//////////////////////////////////////////////////////////////////////
// MapId<0> class
//////////////////////////////////////////////////////////////////////

template<> class MapId<0>  : public Map
{
public: 

	MapId();
	MapId(const unsigned int &);
	virtual ~MapId();
	MapId(const MapId<0> &);
	MapId<0> & operator = (const MapId<0> &);
	Set::Manifold::Map *Clone();
	Set::Manifold::TMap *Diff();
	void Randomize();
	const Point & operator () (const Point &);
	unsigned int size1() const;
	unsigned int size2() const;

private:

	unsigned int n;
};

//////////////////////////////////////////////////////////////////////
// MapId<1> class
//////////////////////////////////////////////////////////////////////

template<> class MapId<1>  : public TMap
{
public: 

	MapId();
	MapId(const unsigned int &);
	virtual ~MapId();
	MapId(const MapId<0> &);
	MapId(const MapId<1> &);
	MapId<1> & operator = (const MapId<1> &);
	Set::Manifold::TMap *Clone();
	Set::Manifold::TMap *Diff();
	void Randomize();
	const Set::VectorSpace::Hom & operator () (const Point &);
	unsigned int size1() const;
	unsigned int size2() const;

private:

	unsigned int n;
	Set::VectorSpace::HomId I;
};

//////////////////////////////////////////////////////////////////////
// MapId<2> class
//////////////////////////////////////////////////////////////////////

template<> class MapId<2>  : public TMap
{
public: 

	MapId();
	MapId(const unsigned int &);
	virtual ~MapId();
	MapId(const MapId<1> &);
	MapId(const MapId<2> &);
	MapId<2> & operator = (const MapId<2> &);
	Set::Manifold::TMap *Clone();
	void Randomize();
	const Set::VectorSpace::Hom & operator () (const Point &);
	unsigned int size1() const;
	unsigned int size2() const;

private:

	unsigned int n;
	Set::VectorSpace::HomZero O;

private:

	Set::Manifold::TMap *Diff(){return 0;}
};

}

}

}

#endif // !defined(SET_MANIFOLD_CATEGORY_MAPID_H__INCLUDED_)
