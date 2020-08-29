// Conforming.cpp: implementation of the Conforming class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./LumpedMass.h"

namespace Element
{
namespace Conforming
{
LumpedMass::LumpedMass() {}

LumpedMass::~LumpedMass() {}

LumpedMass::LumpedMass(const double &rhs_rho) 
	: rho(rhs_rho) {}

LumpedMass::LumpedMass(
	const double &rhs_rho,
	const vector<map<Set::Manifold::Point *, double> > &rhs_N,
	const vector<double> &rhs_QW) :
	rho(rhs_rho), N(rhs_N), QW(rhs_QW) {}

LumpedMass::LumpedMass(const LumpedMass &rhs) :
	rho(rhs.rho), N(rhs.N), QW(rhs.QW) {}

LumpedMass & 
LumpedMass::operator = (const LumpedMass &rhs)
{
	if (this == &rhs) return *this; 
	rho = rhs.rho; N = rhs.N; QW = rhs.QW;
	return *this;
}

set<Set::Manifold::Point *>
LumpedMass::GetNodes() const
{
	set<Set::Manifold::Point *> Nodes;
	map<Set::Manifold::Point *, double>::const_iterator pN;
	for (pN=N[0].begin(); pN!=N[0].end(); pN++)
		Nodes.insert(pN->first);
	return Nodes;
}

map<Set::Manifold::Point *, double>
LumpedMass::GetMass() const
{
	map<Set::Manifold::Point *, double> Mass;
	map<Set::Manifold::Point *, double>::const_iterator pN;
	for (pN=N[0].begin(); pN!=N[0].end(); pN++)
		Mass.insert(make_pair(pN->first,0.0));
	map<Set::Manifold::Point *, double>::iterator pM;
	for (unsigned int a=0; a<QW.size(); a++)
		for (pN=N[a].begin(), pM=Mass.begin(); pN!=N[a].end(); pN++, pM++)
			pM->second += rho*QW[a]*pN->second;
	return Mass;
}

void
LumpedMass::GetMass(map<Set::Manifold::Point *, double> & Mass) const
{
	if (!Mass.empty()) Mass.clear();

	const map<Set::Manifold::Point *, double> & N0 = N[0];
	map<Set::Manifold::Point *, double>::const_iterator pN;	  

	for (pN=N0.begin(); pN!=N0.end(); pN++) {
	  Mass.insert(make_pair(pN->first, 0.0)); 
	}

	for ( int q = 0; q < QW.size(); ++q ) {
	  const map<Set::Manifold::Point *, double> & Nloc = N[q];
	  for (pN=Nloc.begin(); pN!=Nloc.end(); pN++) {
	    Mass.find(pN->first)->second += rho * QW[q] * pN->second;
	  }
	}

	return;
}

}

}
