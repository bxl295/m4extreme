// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_EUCLIDEAN_CARTESIAN_FACTORY_H__INCLUDED_)
#define SET_EUCLIDEAN_CARTESIAN_FACTORY_H__INCLUDED_

#pragma once

#include <map>
#include "./Cartesian.h"
#include "../../Factory.h"

namespace Set
{
namespace Euclidean
{	
namespace Cartesian
{

//////////////////////////////////////////////////////////////////////
// Class FixedVector
//////////////////////////////////////////////////////////////////////

class FixedVector : public Set::Manifold::Factory
{
public: 

	FixedVector() : Emb(0), DEmb(0), DDEmb(0), Sub(0), DSub(0), DDSub(0) {}

	virtual ~FixedVector() 
	{
		if (DDSub != 0) delete DDSub;
		if (DSub != 0) delete DSub;
		if (Sub != 0) delete Sub;
		if (DDEmb != 0) delete DDEmb;
		if (DEmb != 0) delete DEmb;
		if (Emb != 0) delete Emb;
	} 

	void operator += (const Set::Euclidean::Orthonormal::Vector &u)
	{
		if (Emb != 0) Emb->RangeOrigin() += u;
		if (Sub != 0) Sub->DomainOrigin() += u;
	}

	void operator -= (const Set::Euclidean::Orthonormal::Vector &u)
	{
		if (Emb != 0) Emb->RangeOrigin() -= u;
		if (Sub != 0) Sub->DomainOrigin() -= u;
	}

	Set::Manifold::Map *GetEmb(const Set::Manifold::Point &Q)
	{
		if (Emb == 0)
		{
			Set::Euclidean::Cartesian::Point P(0); 
			Set::Euclidean::Orthonormal::Point O(Q.size());
			Set::VectorSpace::Hom A(Q.size(),(unsigned int)0);
			Emb = new Set::Euclidean::Cartesian::Embedding<0>(P,O,A);
		}
		return Emb;
	}

	Set::Manifold::TMap *GetDEmb(const Set::Manifold::Point &Q)
	{
		if (DEmb == 0)
		{
			Set::VectorSpace::Hom A(Q.size(),(unsigned int)0);
			DEmb = new Set::Euclidean::Cartesian::Embedding<1>(A);
		}
		return DEmb;
	}

	Set::Manifold::TMap *GetDDEmb(const Set::Manifold::Point &Q)
	{
		if (DDEmb == 0) DDEmb = new Set::Euclidean::Cartesian::Embedding<2>();
		return DDEmb;
	}

	Set::Manifold::Map *GetSub(const Set::Manifold::Point &Q)
	{
		if (Sub == 0)
		{
			Set::Euclidean::Cartesian::Point P(0); 
			Set::Euclidean::Orthonormal::Point O(Q.size());
			Set::VectorSpace::Hom A((unsigned int)0,Q.size());
			Sub = new Set::Euclidean::Cartesian::Submersion<0>(O,P,A);
		}
		return Sub;
	}

	Set::Manifold::TMap *GetDSub(const Set::Manifold::Point &Q)
	{
		if (DSub == 0)
		{
			Set::VectorSpace::Hom A((unsigned int)0,Q.size());
			DSub = new Set::Euclidean::Cartesian::Submersion<1>(A);
		}
		return DSub;
	}

	Set::Manifold::TMap *GetDDSub(const Set::Manifold::Point &Q)
	{
		if (DDSub == 0) DDSub = new Set::Euclidean::Cartesian::Submersion<2>();
		return DDSub;
	}

	void operator += (const Set::VectorSpace::Vector &u)
	{
		this->operator +=((const Set::Euclidean::Orthonormal::Vector &)u);
	}

	void operator -= (const Set::VectorSpace::Vector &u)
	{
		this->operator -=((const Set::Euclidean::Orthonormal::Vector &)u);
	}

private:

	Set::Euclidean::Cartesian::Embedding<0> *Emb;
	Set::Euclidean::Cartesian::Embedding<1> *DEmb;
	Set::Euclidean::Cartesian::Embedding<2> *DDEmb;
	Set::Euclidean::Cartesian::Submersion<0> *Sub;
	Set::Euclidean::Cartesian::Submersion<1> *DSub;
	Set::Euclidean::Cartesian::Submersion<2> *DDSub;

private:

	FixedVector(const FixedVector &);
	FixedVector & operator = (const FixedVector &);
};

//////////////////////////////////////////////////////////////////////
// Class FixedPoint
//////////////////////////////////////////////////////////////////////

class FixedPoint : public Set::Manifold::Factory
{
public: 

    typedef Set::Euclidean::Orthonormal::Point coor_type;
    typedef Set::Euclidean::Cartesian::Embedding<0> emb_type;
    typedef Set::Euclidean::Cartesian::Submersion<0> sub_type;

    FixedPoint() : DEmb(0), DDEmb(0), DSub(0), DDSub(0) {}
    virtual ~FixedPoint() 
    {
		map<const Set::Manifold::Point *, sub_type *>::iterator pSub = Sub.begin();
		while(pSub!=Sub.end()){delete pSub->second; pSub++;}
		if (DDSub != 0) delete DDSub;
		if (DSub != 0) delete DSub;
		map<const Set::Manifold::Point *, emb_type *>::iterator pEmb = Emb.begin();
		while(pEmb!=Emb.end()){delete pEmb->second; pEmb++;}
		if (DDEmb != 0) delete DDEmb;
		if (DEmb != 0) delete DEmb;
    } 

    Set::Manifold::Map *GetEmb(const Set::Manifold::Point &Q)
    {
		const Set::Manifold::Point *pQ = &Q;
		map<const Set::Manifold::Point *, emb_type *>::const_iterator pEmb;
		pEmb = Emb.find(pQ);
		if (pEmb==Emb.end())
		{
			Set::Euclidean::Cartesian::Point P(0); 
			Set::VectorSpace::Hom A(Q.size(),(unsigned int)0);
			emb_type *Embloc = 
			new emb_type(P,(coor_type &)Q,A);
			Emb.insert(make_pair(pQ, Embloc));
			return Embloc;
		}
		else
		{
			return pEmb->second;
		}
    }

    Set::Manifold::TMap *GetDEmb(const Set::Manifold::Point &Q)
    {
		if (DEmb == 0)
		{
			Set::VectorSpace::Hom A(Q.size(),(unsigned int)0);
			DEmb = new Set::Euclidean::Cartesian::Embedding<1>(A);
		}
		return DEmb;
    }

    Set::Manifold::TMap *GetDDEmb(const Set::Manifold::Point &Q)
    {
		if (DDEmb == 0) DDEmb = new Set::Euclidean::Cartesian::Embedding<2>();
		return DDEmb;
    }

    Set::Manifold::Map *GetSub(const Set::Manifold::Point &Q)
    {
		const Set::Manifold::Point *pQ = &Q;
		map<const Set::Manifold::Point *, sub_type *>::iterator pSub;
		pSub = Sub.find(pQ);
		if (pSub==Sub.end())
		{
			Set::Euclidean::Cartesian::Point P(0); 
			Set::VectorSpace::Hom A((unsigned int)0,Q.size());
			sub_type *Subloc = 
			new sub_type((coor_type &)Q,P,A);
			Sub.insert(make_pair(pQ, Subloc));
			return Subloc;
		}
		else
		{
			return pSub->second;
		}
    }

    Set::Manifold::TMap *GetDSub(const Set::Manifold::Point &Q)
    {
		if (DSub == 0)
		{
			Set::VectorSpace::Hom A((unsigned int)0,Q.size());
			DSub = new Set::Euclidean::Cartesian::Submersion<1>(A);
		}
		return DSub;
    }

    Set::Manifold::TMap *GetDDSub(const Set::Manifold::Point &Q)
    {
		if (DDSub == 0) DDSub = new Set::Euclidean::Cartesian::Submersion<2>();
		return DDSub;
    }

    void operator += (const Set::Euclidean::Orthonormal::Vector &u)
    {
		map<const Set::Manifold::Point *, emb_type *>::iterator pEmb;	
		map<const Set::Manifold::Point *, sub_type *>::iterator pSub;	
		for(pEmb = Emb.begin(); pEmb!=Emb.end(); pEmb++)
		{
			assert(u.size() == pEmb->first->size()); 
			pEmb->second->RangeOrigin() += u;
		}
		for(pSub = Sub.begin(); pSub!=Sub.end(); pSub++)
		{
			assert(u.size() == pSub->first->size()); 
			pSub->second->DomainOrigin() += u;
		}
    }

    void operator -= (const Set::Euclidean::Orthonormal::Vector &u)
    {
		map<const Set::Manifold::Point *, emb_type *>::iterator pEmb;	
		map<const Set::Manifold::Point *, sub_type *>::iterator pSub;	
		for(pEmb = Emb.begin(); pEmb!=Emb.end(); pEmb++)
		{
			assert(u.size() == pEmb->first->size()); 
			pEmb->second->RangeOrigin() -= u;
		}
		for(pSub = Sub.begin(); pSub!=Sub.end(); pSub++)
		{
			assert(u.size() == pSub->first->size()); 
			pSub->second->DomainOrigin() -= u;
		}
    }

	void operator += (const Set::VectorSpace::Vector &u)
	{
		this->operator +=((const Set::Euclidean::Orthonormal::Vector &)u);
	}

	void operator -= (const Set::VectorSpace::Vector &u)
	{
		this->operator -=((const Set::Euclidean::Orthonormal::Vector &)u);
	}

private:

	map<const Set::Manifold::Point *, Set::Euclidean::Cartesian::Embedding<0> *>  Emb;
	Set::Euclidean::Cartesian::Embedding<1> *DEmb;
	Set::Euclidean::Cartesian::Embedding<2> *DDEmb;
	map<const Set::Manifold::Point *, Set::Euclidean::Cartesian::Submersion<0> *> Sub;
	Set::Euclidean::Cartesian::Submersion<1> *DSub;
	Set::Euclidean::Cartesian::Submersion<2> *DDSub;

private:

	FixedPoint(const FixedPoint &);
	FixedPoint & operator = (const FixedPoint &);
};

//////////////////////////////////////////////////////////////////////
// Class FixedField
//////////////////////////////////////////////////////////////////////

class FixedField : public Set::Manifold::Factory
{
public: 

	FixedField() : n(0), Emb(0), DEmb(0), Sub(0), DSub(0) {}
	virtual ~FixedField() 
	{
		if (DSub != 0) delete DSub;
		if (Sub != 0) delete Sub;
		if (DEmb != 0) delete DEmb;
		if (Emb != 0) delete Emb;
	} 

	FixedField(const unsigned int &n_) : n(n_), 
		Q(Set::Euclidean::Orthonormal::Point(n)),
		Emb(0), DEmb(0), Sub(0), DSub(0) {}

	Set::Manifold::Map *GetEmb()
	{
		if (Emb == 0)
		{
			Set::Euclidean::Cartesian::Point P(0); 
			Set::VectorSpace::Hom A(n,(unsigned int)0);
			Emb = new Set::Euclidean::Cartesian::Embedding<0>(P,Q,A);
		}
		return Emb;
	}

	Set::Manifold::TMap *GetDEmb()
	{
		if (DEmb == 0)
		{
			Set::VectorSpace::Hom A(n,(unsigned int)0);
			DEmb = new Set::Euclidean::Cartesian::Embedding<1>(A);
		}
		return DEmb;
	}

    Set::Manifold::TMap *GetDDEmb()
    {
		if (DDEmb == 0) DDEmb = new Set::Euclidean::Cartesian::Embedding<2>();
		return DDEmb;
    }

	Set::Manifold::Map *GetSub()
	{
		if (Sub == 0)
		{
			Set::Euclidean::Cartesian::Point P(0); 
			Set::VectorSpace::Hom A((unsigned int)0,n);
			Sub = new Set::Euclidean::Cartesian::Submersion<0>(Q,P,A);
		}
		return Sub;
	}

	Set::Manifold::TMap *GetDSub()
	{
		if (DSub == 0)
		{
			Set::VectorSpace::Hom A((unsigned int)0,n);
			DSub = new Set::Euclidean::Cartesian::Submersion<1>(A);
		}
		return DSub;
	}

    Set::Manifold::TMap *GetDDSub()
    {
		if (DDSub == 0) DDSub = new Set::Euclidean::Cartesian::Submersion<2>();
		return DDSub;
    }

	void operator += (const Set::Euclidean::Orthonormal::Vector &u)
	{
		assert(u.size() == n); Q += u;
		if (Emb != 0) Emb->RangeOrigin() += u;
		if (Sub != 0) Sub->DomainOrigin() += u;
	}

	void operator -= (const Set::Euclidean::Orthonormal::Vector &u)
	{
		assert(u.size() == n); Q -= u;
		if (Emb != 0) Emb->RangeOrigin() -= u;
		if (Sub != 0) Sub->DomainOrigin() -= u;
	}

	Set::Manifold::Map *GetEmb(const Set::Manifold::Point &){return GetEmb();}
	Set::Manifold::TMap *GetDEmb(const Set::Manifold::Point &){return GetDEmb();}
	Set::Manifold::TMap *GetDDEmb(const Set::Manifold::Point &){return GetDDEmb();}
	Set::Manifold::Map *GetSub(const Set::Manifold::Point &){return GetSub();}
	Set::Manifold::TMap *GetDSub(const Set::Manifold::Point &){return GetDSub();}
	Set::Manifold::TMap *GetDDSub(const Set::Manifold::Point &){return GetDDSub();}

	void operator += (const Set::VectorSpace::Vector &u)
	{
		this->operator +=((const Set::Euclidean::Orthonormal::Vector &)u);
	}

	void operator -= (const Set::VectorSpace::Vector &u)
	{
		this->operator -=((const Set::Euclidean::Orthonormal::Vector &)u);
	}

private:

	unsigned int n;
	Set::Euclidean::Orthonormal::Point Q;

	Set::Euclidean::Cartesian::Embedding<0> *Emb;
	Set::Euclidean::Cartesian::Embedding<1> *DEmb;
	Set::Euclidean::Cartesian::Embedding<2> *DDEmb;
	Set::Euclidean::Cartesian::Submersion<0> *Sub;
	Set::Euclidean::Cartesian::Submersion<1> *DSub;
	Set::Euclidean::Cartesian::Submersion<2> *DDSub;

private:

	FixedField(const FixedField &);
	FixedField & operator = (const FixedField &);
};

}

}

}

#endif // !defined(SET_EUCLIDEAN_CARTESIAN_FACTORY_H__INCLUDED_
