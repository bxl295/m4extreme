// Conforming.cpp: implementation of the Conforming class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include <utility>
#include "./Conforming.h"
#include "../../Utils/Indexing/Indexing.h"

extern map<Set::Manifold::Point *, unsigned int> _gIndexSet;

namespace Element
{
namespace Conforming
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState() {}

LocalState::~LocalState() {}

Element::LocalState *
LocalState::Clone() const
{
	return new LocalState(*this);
}

LocalState::LocalState(const vector<Material::LocalState *> &MatLS_) 
	: MatLS(MatLS_) {}

LocalState::LocalState(
	const vector<Material::LocalState *> &MatLS_,
	const vector<double> &QW_) 
	: MatLS(MatLS_), QW(QW_) {}

LocalState::LocalState(
	const vector<Material::LocalState *> &MatLS_,
	const vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > &DN_,
	const vector<double> &QW_) : 
	MatLS(MatLS_), DN(DN_), QW(QW_) {}

LocalState::LocalState(const LocalState &rhs) : 
	MatLS(rhs.MatLS), DN(rhs.DN), QW(rhs.QW) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	MatLS = rhs.MatLS; DN = rhs.DN; QW = rhs.QW;
	return *this;
}

const vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > &
LocalState::GetDN() const
{
	return DN;
}

vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > &
LocalState::GetDN() 
{
	return DN;
}

const vector<double> &
LocalState::GetQW() const
{
	return QW;
}

vector<double> &
LocalState::GetQW()
{
	return QW;
}

vector<Set::VectorSpace::Vector>
LocalState::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	unsigned int q; 
	unsigned int m = x.begin()->second.size();
	unsigned int n = DN[0].begin()->second.size();
	Set::VectorSpace::Vector O(m*n);
	vector<Set::VectorSpace::Vector> Dx(QW.size(),O);
	for (q=0; q<QW.size(); q++)
	{
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
		for (pN=DN[q].begin(), px=x.begin(); pN!=DN[q].end(); pN++, px++)
			Dx[q] += Dyadic(pN->second,px->second);
	}
	return Dx;
}

void 
LocalState::operator ++ ()
{
	vector<Material::LocalState *>::iterator pMatLS;
	for (pMatLS=MatLS.begin(); pMatLS!=MatLS.end(); pMatLS++) ++(*(*pMatLS));
}

set<Set::Manifold::Point *>
LocalState::GetNodes() const
{
	set<Set::Manifold::Point *> Nodes;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pN;
	for (pN=DN[0].begin(); pN!=DN[0].end(); pN++)
		Nodes.insert(pN->first);
	return Nodes;
}

void
LocalState::GetNodes(set<Set::Manifold::Point *> & Nodes) const
{
	if ( !Nodes.empty() ) Nodes.clear();

	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pN;
	for (pN=DN[0].begin(); pN!=DN[0].end(); pN++)
		Nodes.insert(pN->first);
	return;
}

set<pair<Set::Manifold::Point *,Set::Manifold::Point *> >
LocalState::GetNodePairs() const
{
	set<pair<Set::Manifold::Point *,Set::Manifold::Point *> > NodePairs;
	set<Set::Manifold::Point *> Nodes = GetNodes();
	set<Set::Manifold::Point *>::iterator pM;
	set<Set::Manifold::Point *>::iterator pN;
	for (pM=Nodes.begin(); pM!=Nodes.end(); pM++)
		for (pN=Nodes.begin(); pN!=Nodes.end(); pN++)
			NodePairs.insert(make_pair(*pM,*pN));
	return NodePairs;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy() {}

Energy<0>::~Energy() {}

Element::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Element::LocalState *
Energy<0>::GetLocalState() const
{
	return LS;
}

Energy<0>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<0> *> &rhs_W) 
	: LS(rhs_LS), W(rhs_W) {}

Energy<0>::Energy(const Energy<0> &rhs)
: LS(rhs.LS), W(rhs.W) {}

Energy<0> & 
Energy<0>::operator = (const Energy<0> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; W = rhs.W;
	return *this;
}

double 
Energy<0>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	unsigned int q; double E = 0.0;
	vector<Set::VectorSpace::Vector> Dx = (*LS)(x);
	for (q=0; q<LS->QW.size(); q++) 
		E += LS->QW[q]*(*W[q])(Dx[q]);
	return E;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {}

Energy<1>::~Energy() {}

Element::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Element::LocalState *
Energy<1>::GetLocalState() const
{
	return LS;
}

Energy<1>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<1> *> &rhs_DW) 
	: LS(rhs_LS), DW(rhs_DW) {
  Element::Energy<1>::_ELS = LS;
}

Energy<1>::Energy(const Energy<1> &rhs) 
: LS(rhs.LS), DW(rhs.DW) {}

Energy<1> & 
Energy<1>::operator = (const Energy<1> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; DW = rhs.DW;
	return *this;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
Energy<1>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	unsigned int q;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> DE;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
	unsigned int m = x.begin()->second.size();
	unsigned int n = LS->DN[0].begin()->second.size();
	for (px=x.begin(); px!=x.end(); px++)
	{
		assert(m == px->second.size());
		DE.insert(range_type::value_type(px->first,
			Set::VectorSpace::Vector(px->second.size())));
	}
	vector<Set::VectorSpace::Vector> Dx = (*LS)(x);
	for (q=0; q<LS->QW.size(); q++)
	{
		Set::VectorSpace::Hom P(m,n); P = (*DW[q])(Dx[q]);
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDE;
		for (pN=LS->DN[q].begin(), pDE=DE.begin(); pN!=LS->DN[q].end(); pN++, pDE++)
			pDE->second += LS->QW[q]*P*pN->second;
	}
	return DE;
}

void	
Energy<1>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x,
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DE) const
{           
	if ( !DE.empty() ) DE.clear();
	unsigned int q;
	
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
	unsigned int m = x.begin()->second.size();
	unsigned int n = LS->DN[0].begin()->second.size();
	for (px=x.begin(); px!=x.end(); px++)
	{
		assert(m == px->second.size());
		DE.insert(range_type::value_type(px->first,
			Set::VectorSpace::Vector(px->second.size())));
	}
        
        
#ifdef _M4EXTREME_DEBUG_II
            std::cout << "\nCompute the internal force at Conforming element: [\t";
            for (px = x.begin(); px != x.end(); px++) {
                std::cout << _gIndexSet.find(px->first)->second << "\t";
            }
            std::cout << "]" << std::endl;
#endif
            
	vector<Set::VectorSpace::Vector> Dx = (*LS)(x);
	for (q=0; q<LS->QW.size(); q++)
	{
		Set::VectorSpace::Hom P(m,n); P = (*DW[q])(Dx[q]);
                
#ifdef _M4EXTREME_DEBUG_II
            std::cout << "\nquadrature point "<<q<<": F=[" << Dx[q] << "]\tP=[" << P << "]" << std::endl;
#endif
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDE;
		for (pN=LS->DN[q].begin(), pDE=DE.begin(); pN!=LS->DN[q].end(); pN++, pDE++)
			pDE->second += LS->QW[q]*P*pN->second;
	}

	return;
}


//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {}

Energy<2>::~Energy() {}

Element::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Element::LocalState *
Energy<2>::GetLocalState() const
{
	return LS;
}

Energy<2>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<2> *> &rhs_DDW) 
	: LS(rhs_LS), DDW(rhs_DDW) {}

Energy<2>::Energy(const Energy<2> &rhs)
: LS(rhs.LS), DDW(rhs.DDW) {}

Energy<2> & 
Energy<2>::operator = (const Energy<2> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; DDW = rhs.DDW;
	return *this;
}

map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom>
Energy<2>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	unsigned int i, j, k, l, q; 
	unsigned int m = x.begin()->second.size();
	unsigned int n = LS->DN[0].begin()->second.size();
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
	for (px=x.begin(); px!=x.end(); px++)
		assert(m == px->second.size());
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px1;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px2;
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> DDE;
	for (px1=x.begin(); px1!=x.end(); px1++)
		for (px2=x.begin(); px2!=x.end(); px2++)
			DDE.insert(range_type::value_type(
				make_pair(px1->first,px2->first),Set::VectorSpace::Hom(m)));
	vector<Set::VectorSpace::Vector> Dx = (*LS)(x);
	for (q=0; q<LS->QW.size(); q++)
	{
		Set::VectorSpace::Hom C = (*DDW[q])(Dx[q]);
		double **** CT = Indexing::New(C.begin(),m,n,m,n);
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN1;
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
			Set::VectorSpace::Hom>::iterator pDDE;
		for (pN1=LS->DN[q].begin(), pDDE=DDE.begin(); pN1!=LS->DN[q].end(); pN1++)
		{
			Set::VectorSpace::Vector B1 = pN1->second;
			map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN2;
			for (pN2=LS->DN[q].begin(); pN2!=LS->DN[q].end(); pN2++, pDDE++)
			{
				Set::VectorSpace::Vector B2 = pN2->second;
				Set::VectorSpace::Hom A(m,m);
				for (k=0; k<m; k++)
					for (i=0; i<m; i++)
						for (j=0; j<n; j++)
							for (l=0; l<n; l++)
								A[k][i] +=  CT[l][k][j][i]*B2[l]*B1[j];
				pDDE->second += LS->QW[q]*A;
			}
		}
		Indexing::Delete(CT,m,n,m,n);
	}
	return DDE;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {}

Jet<0>::~Jet() {} 

Element::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Element::LocalState *
Jet<0>::GetLocalState() const
{
	return LS;
}

Jet<0>::Jet(
	LocalState *rhs_LS, 
	const vector<Material::Jet<0> *> &rhs_J)
	: LS(rhs_LS), J(rhs_J) {}

Jet<0>::Jet(const Jet<0> &rhs)
	: LS(rhs.LS), J(rhs.J) {}

Jet<0> & 
Jet<0>::operator = (const Jet<0> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; J = rhs.J;
	return *this;
}

pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> > 
Jet<0>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const 
{
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> de_type;
	unsigned int q; double E = 0.0;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> DE;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
	unsigned int m = x.begin()->second.size();
	unsigned int n = LS->DN[0].begin()->second.size();
	for (px=x.begin(); px!=x.end(); px++)
	{
		assert(m == px->second.size());
		DE.insert(de_type::value_type(px->first,
			Set::VectorSpace::Vector(px->second.size())));
	}
	vector<Set::VectorSpace::Vector> Dx = (*LS)(x);
	for (q=0; q<LS->QW.size(); q++) 
	{
		pair<double,Set::VectorSpace::Vector> Jq = (*J[q])(Dx[q]);
		E += LS->QW[q]*Jq.first;
		Set::VectorSpace::Hom P(m,n); P = Jq.second;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDE;
		for (pN=LS->DN[q].begin(), pDE=DE.begin(); pN!=LS->DN[q].end(); pN++, pDE++)
			pDE->second += LS->QW[q]*P*pN->second;
	}
	return make_pair(E,DE);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {}

Jet<1>::~Jet() {}

Element::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Element::LocalState *
Jet<1>::GetLocalState() const
{
	return LS;
}

Jet<1>::Jet(
	LocalState *rhs_LS, 
	const vector<Material::Jet<1> *> &rhs_DJ)
	: LS(rhs_LS), DJ(rhs_DJ) {}

Jet<1>::Jet(const Jet<1> &rhs)
	: LS(rhs.LS), DJ(rhs.DJ) {}

Jet<1> & 
Jet<1>::operator = (const Jet<1> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; DJ = rhs.DJ;
	return *this;
}

pair<map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
	Set::VectorSpace::Hom> >
Jet<1>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> de_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> dde_type;
	unsigned int i, j, k, l, q; 
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> DE;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
	unsigned int m = x.begin()->second.size();
	unsigned int n = LS->DN[0].begin()->second.size();
	for (px=x.begin(); px!=x.end(); px++)
	{
		assert(m == px->second.size());
		DE.insert(de_type::value_type(px->first,
		Set::VectorSpace::Vector(px->second.size())));
	}
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> DDE;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px1;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px2;
	for (px1=x.begin(); px1!=x.end(); px1++)
		for (px2=x.begin(); px2!=x.end(); px2++)
			DDE.insert(dde_type::value_type(
			make_pair(px1->first,px2->first),
			Set::VectorSpace::Hom(m)));
	vector<Set::VectorSpace::Vector> Dx = (*LS)(x);
	for (q=0; q<LS->QW.size(); q++)
	{
		pair<Set::VectorSpace::Vector, Set::VectorSpace::Hom> 
			DJq = (*DJ[q])(Dx[q]);
		Set::VectorSpace::Hom P(m,n); P = DJq.first;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDE;
		for (pN=LS->DN[q].begin(), pDE=DE.begin(); pN!=LS->DN[q].end(); pN++, pDE++)
			pDE->second += LS->QW[q]*P*pN->second;
		Set::VectorSpace::Hom C = DJq.second;
		double **** CT = Indexing::New(C.begin(),m,n,m,n);
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN1;
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
			Set::VectorSpace::Hom>::iterator pDDE;
		for (pN1=LS->DN[q].begin(), pDDE=DDE.begin(); pN1!=LS->DN[q].end(); pN1++)
		{
			Set::Manifold::Point *e1 = pN1->first;
			Set::VectorSpace::Vector B1 = pN1->second;
			map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pN2;
			for (pN2=LS->DN[q].begin(); pN2!=LS->DN[q].end(); pN2++, pDDE++)
			{
				Set::Manifold::Point *e2 = pN2->first;
				Set::VectorSpace::Vector B2 = pN2->second;
				Set::VectorSpace::Hom A(m,m);
				for (k=0; k<m; k++)
					for (i=0; i<m; i++)
						for (j=0; j<n; j++)
							for (l=0; l<n; l++)
								A[k][i] +=  CT[l][k][j][i]*B2[l]*B1[j];
				pDDE->second += LS->QW[q]*A;
			}
		}
		Indexing::Delete(CT,m,n,m,n);
	}
	return make_pair(DE,DDE);
}

}

}
