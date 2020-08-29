// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#include "./LumpedMass.h"

namespace Element
{
namespace MaterialPoint
{
LumpedMass::LumpedMass() {}

LumpedMass::~LumpedMass() {}

  //
  // calculate the mass at the reference configuration 
  // because the mass of a material point remains constant during the deformation
  // however rho and volume are dynamically changing
  //
LumpedMass::LumpedMass(Data *D_, const double &rho_) : D(D_) {
  for ( int q = 0; q < D->QW.size(); ++q ) {
    dm_0.push_back(rho_ * (D->QW)[q]);
  }
}

LumpedMass::LumpedMass(LocalState * pLS, const double &rho_) : _pLS(pLS) {
  D = _pLS->GetData();
  for ( int q = 0; q < D->QW.size(); ++q ) {
    dm_0.push_back(rho_ * (D->QW)[q]);
  }
}

LumpedMass::LumpedMass(const LumpedMass &rhs) : D(rhs.D), dm_0(rhs.dm_0) {}

LumpedMass & 
LumpedMass::operator = (const LumpedMass &rhs)
{
	if (this == &rhs) return *this; 
	D = rhs.D; dm_0 = rhs.dm_0;
	return *this;
}

set<Set::Manifold::Point *>
LumpedMass::GetNodes() const
{
	set<Set::Manifold::Point *> Nodes;
	const map<Set::Manifold::Point *, double> & Nloc = (D->N)[0];
	map<Set::Manifold::Point *, double>::const_iterator pN;
	for (pN=Nloc.begin(); pN!=Nloc.end(); pN++)
		Nodes.insert(pN->first); 
	return Nodes;
}

map<Set::Manifold::Point *, double>
LumpedMass::GetMass() const
{
	map<Set::Manifold::Point *, double> Mass;

 	const double qw0 = dm_0[0];
	const map<Set::Manifold::Point *, double> & N0 = (D->N)[0];
	map<Set::Manifold::Point *, double>::const_iterator pN;	  
#if defined(_M4EXTREME_EIGEN_FRACTURE_I_)        
	if ( _pLS->isActivated() ) { 
#endif
	  for (pN=N0.begin(); pN!=N0.end(); pN++) {
	    Mass.insert(make_pair(pN->first, qw0 * pN->second)); 
	  }
	  
	  for ( int q = 1; q < dm_0.size(); ++q ) {
	    const double & qwloc = dm_0[q];
	    const map<Set::Manifold::Point *, double> & Nloc = (D->N)[q];
	    for (pN=Nloc.begin(); pN!=Nloc.end(); pN++) {
	      Mass.find(pN->first)->second += qwloc * pN->second;
	    }
	  }

#if defined(_M4EXTREME_EIGEN_FRACTURE_I_)  
	}
	else {
	  for (pN=N0.begin(); pN!=N0.end(); pN++) {
	    Mass.insert(make_pair(pN->first,0.0));
	  }
	}
#endif	
       
	return Mass;
}

void
LumpedMass::GetMass(map<Set::Manifold::Point *, double> & Mass) const
{
	if (!Mass.empty()) Mass.clear();

	const double qw0 = dm_0[0];
	const map<Set::Manifold::Point *, double> & N0 = (D->N)[0];
	map<Set::Manifold::Point *, double>::const_iterator pN;	  
#if defined(_M4EXTREME_EIGEN_FRACTURE_I_)        
	if ( _pLS->isActivated() ) { 
#endif
	  for (pN=N0.begin(); pN!=N0.end(); pN++) {
	    Mass.insert(make_pair(pN->first, qw0 * pN->second)); 
	  }
	  
	  for ( int q = 1; q < dm_0.size(); ++q ) {
	    const double & qwloc = dm_0[q];
	    const map<Set::Manifold::Point *, double> & Nloc = (D->N)[q];
	    for (pN=Nloc.begin(); pN!=Nloc.end(); pN++) {
	      Mass.find(pN->first)->second += qwloc * pN->second;
	    }
	  }

#if defined(_M4EXTREME_EIGEN_FRACTURE_I_)  
	}
	else {
	  for (pN=N0.begin(); pN!=N0.end(); pN++) {
	    Mass.insert(make_pair(pN->first,0.0));
	  }
	}
#endif	

	return;
}

}

}
