// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <utility>
#include <iomanip>
#include "./MaterialPoint.h"
#include "Utils/Indexing/Indexing.h"

namespace Element 
{
  namespace MaterialPoint 
  {

    //////////////////////////////////////////////////////////////////////
    // Class Data
    //////////////////////////////////////////////////////////////////////

    Data::Data()
      : _maxRange(0.0), _cutoff(0.0), _nbs(NULL), _carrier(NULL), _blacklist(NULL), 
	_updateNeighbor(false), _readflag(0), _adaptiveSearch(false), 
	_searchRange(1.5), _MAX_RESET(5), _volume(0.0), _centroid(0), 
        _isSharedInterface(false) {}

    Data::~Data() {}

    Data::Data(const double & Rho0_,
	       const double & Cutoff_,
	       const double & Ratio_,
	       const vector<double> & QW_,
	       const vector<vector_type> & QP_,
	       const map<dof_type *, vector_type> & x_,
	       bool updateNeighbor_,
	       Geometry::Search<dof_type*> ** nbs) 
      : _nbs(nbs), QW(QW_), QP(QP_), N(QW_.size()), DN(QW_.size()),
	_cutoff(Cutoff_), _MAX_RESET(50), _carrier(NULL), _blacklist(NULL), _readflag(0), 
	_adaptiveSearch(false),	_searchRange(1.5), _volume(0.0), _isSharedInterface(false),
	_updateNeighbor(updateNeighbor_), _centroid(x_.begin()->second.size()) {          

      int numofMPT = QW.size();
      assert(numofMPT > 0 );
      
      for ( int q = 0; q < numofMPT; ++q ) {
	_volume += QW[q];
	_centroid += QP[q];
	Mass.push_back(Rho0_ * QW[q]);
	_Jold.push_back(1.0);
	_DJ.push_back(1.0);
      }
      _centroid /= (double)numofMPT;

      double H0 = 0.0;
      map<dof_type*, vector_type>::const_iterator px = x_.begin();
      while ( px != x_.end() ) {
	vector_type dx = px->second - _centroid;
	rxaq.insert( make_pair(px->first, dx) );
	double d2 = dx(dx);
	if ( d2 > H0 ) H0 = d2;
	px++;
      }

      _maxRange = Ratio_ * sqrt(H0);
    }

    Data::Data(const Data &rhs)
      : QW(rhs.QW), QP(rhs.QP), _readflag(0), Mass(rhs.Mass),
	N(rhs.N), DN(rhs.DN), _Jold(rhs._Jold), _blacklist(rhs._blacklist),
	_DJ(rhs._DJ), _nbs(rhs._nbs), _cutoff(rhs._cutoff),
	_adaptiveSearch(rhs._adaptiveSearch), _volume(rhs._volume), 
	_updateNeighbor(rhs._updateNeighbor), _centroid(rhs._centroid),
	_searchRange(rhs._searchRange), _carrier(rhs._carrier) {}

    Data & Data::operator =(const Data &rhs) {

      if (this == &rhs) {
	return *this;
      }
            
      _nbs = rhs._nbs;
      _cutoff = rhs._cutoff;
      _readflag = 0;
      _carrier = rhs._carrier;
      _blacklist = rhs._blacklist;
      _volume = rhs._volume;
      _centroid = rhs._centroid;
      _adaptiveSearch = rhs._adaptiveSearch;
      _searchRange = rhs._searchRange;
      _updateNeighbor = rhs._updateNeighbor;
      _Jold = rhs._Jold;
      _DJ = rhs._DJ;

      Mass = rhs.Mass;
      QW = rhs.QW;      
      QP = rhs.QP;
      N = rhs.N;
      DN = rhs.DN;

      return *this;
    }

    void Data::SetSearch(Geometry::Search<dof_type*> ** nbs) {
      _nbs = nbs;
    }

    void Data::SetAdaptiveSearch(bool adaptive_search, double searchRange) {
      _adaptiveSearch = adaptive_search;
      _searchRange = searchRange;
    }
        
    void Data::SetCarrier(set<dof_type*> * carrier) {
      _carrier = carrier;
    }

    void Data::SetBlackList(set<dof_type*> * blacklist) {
      _blacklist = blacklist;
    }

    void Data::SetMaxReset( unsigned int maxreset ) {
      _MAX_RESET = maxreset;
    }

    void Data::Reset(const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) {

      unsigned int m = x.begin()->second.size();
      unsigned int n = DN[0].begin()->second.size();

      assert(m == n);
      assert(_carrier != NULL);

      int numofMPT = DN.size();
      set<Set::Manifold::Point*> oldNeighbor, newNeighbor, cutNeighbor;

      /**
       * recompute the incremental deformation gradient,
       * the weight and the location of the material point
       */
      for ( int q = 0; q < numofMPT; ++q ) {
	Set::VectorSpace::Hom Dx(m, n);
	const dshape_type & DNloc = DN[q];
	nodeset_type::const_iterator pDN, DNend = DNloc.end();
	for (pDN = DNloc.begin(); pDN != DNend; pDN++) {
	  nodeset_type::const_iterator px = x.find(pDN->first);
	  if ( px != x.end() ) {	 
	    DyadicSum(pDN->second, px->second, Dx);
	  } else {
	    throw(1);
	  }
	}
	
	double Jac = Jacobian(Dx);
	if ( Jac < 1.0e-3 || Jac > 1.0e3 || Jac != Jac ) {
#ifdef _M4EXTREME_DEBUG_II
	  cerr << "error at recalculating the deformation gradient in Element::MaterialPoint::Data::Reset function:" << endl
	       << "Dx = [" << Dx << "], Jacobian=" << Jac << endl; 
#endif
	  assert(false);
	}
	
	QW[q] *= Jac;
	//QW[q] *= _DJ[q];
	
#ifdef _M4EXTREME_DEBUG_II
	cout << "\nReset Element [";
#endif
	
	vector_type & xq = QP[q];
	vector_type xpold = xq;
	Null(xq);
	
	const shape_type & Nloc = N[q];
	shape_type::const_iterator pN;

	for (pN = Nloc.begin(); pN != Nloc.end(); pN++) {  
	  nodeset_type::const_iterator px = x.find(pN->first);	  
	  if ( px != x.end() ) {
	    oldNeighbor.insert(pN->first);
	    xq += pN->second * px->second;
	    if ( fabs(pN->second) > _cutoff
#if defined(_M4EXTREME_SHARED_INTERFACE_) // only include nodes belonging to the same body.
		                       // contact is realized by shared nodes on the interface
		 && _carrier->find(pN->first) != _carrier->end()
#elif defined(_M4EXTREME_MIXING_SEPERATION_) // include nodes beloning to the contact list.
                                          // for material points moving away from each other,
                                          // nodes from different body will be excluded                  
		 && _nodesInTension.find(pN->first) == _nodesInTension.end()
#endif
		 ) {
	      newNeighbor.insert(pN->first);
	    }
	    else {
	      cutNeighbor.insert(pN->first);
#ifdef _M4EXTREME_DEBUG_II
	      cout << "("<<px->second <<")\t"; ////
#endif
	    }
	  } else {
	    throw(1);
	  }
	  
	}
	
#ifdef _M4EXTREME_DEBUG_II
	cout << "]" << endl;
	cout << "new position for material point is [" << xq << "] with new weight "<< QW
	     << " and _maxRange=" << _maxRange << endl;
#endif
      }

      /**
       * update the volume and barycenter of the element
       */
      _volume = 0.0;
      Null(_centroid);
      for ( int q = 0; q < numofMPT; ++q ) {
	_volume += QW[q];
	_centroid += QP[q];
      }
      _centroid /= (double)numofMPT;

      /**
       * reset the neighborhood for all the material points
       */
      nodeset_type xloc;
      if ( _updateNeighbor ) {   
	assert(x.size() > m + 1 && _nbs != NULL);
                
	// Rset the support of all the material points in the element 
	// by using their barycenter as the center of the neighborhood
	double maxRange2 = _maxRange * _maxRange;
	double dxq, xq_distance;
	const double *xq_head = _centroid.begin();
	bool needReset = false;
	int numReset = 0;
	set<Set::Manifold::Point*>::const_iterator pNH, newNeighborEnd = newNeighbor.end();                

	do {
                    
	  needReset = false;
                    
#ifdef _M4EXTREME_DEBUG_II
	  cout << "new neighbor is: {"; ////
#endif                    
	  xloc.clear();
	  pNH = newNeighbor.begin();
                    
	  // the new neighbor must be a neighboring node of the previour neighbors
	  // note that _nbs only contains the attached nodes, if pNH is detached,
	  // this function call will throw out an exception and the node won't be
	  // part of the new neighborhood.
	  if (_adaptiveSearch) {
	    double factorial = 1.0;
	    for (int i = 2; i <= m; ++i) factorial *= (double) i;
	    double radius = _searchRange * pow(factorial * _volume, 1.0 / (double) m);
	    vector<Set::Manifold::Point*> nhloc;

	    while (pNH != newNeighborEnd) {
	      Set::Euclidean::Orthonormal::Point center(m, x.find(*pNH)->second.begin());
	      (**_nbs)(center, radius, nhloc);
                            
	      if (!nhloc.empty()) {
		// the new neighbor must belong to the support region of _centroid
		vector<Set::Manifold::Point*>::const_iterator pNHloc = nhloc.begin(), nhEnd = nhloc.end();
		while (pNHloc != nhEnd) {
		  if (cutNeighbor.find(*pNHloc) == cutNeighbor.end()
#if defined(_M4EXTREME_SHARED_INTERFACE_) // only include nodes belonging to the same body.
		                       // contact is realized by shared nodes on the interface
		      && _carrier->find(*pNHloc) != _carrier->end()
#elif defined(_M4EXTREME_MIXING_SEPERATION_) // include nodes beloning to the contact list.
                                          // for material points moving away from each other,
                                          // nodes from different body will be excluded                  
		      && _nodesInTension.find(*pNHloc) == _nodesInTension.end()
#endif      
		      ) {

		    if ( x.find(*pNHloc) != x.end() ) {
		      const Set::VectorSpace::Vector & px = x.find(*pNHloc)->second;		    
		    
		      xq_distance = 0.0;
		      const double *px_head = px.begin();					
		      for ( int k = 0; k < m; ++k ) {
			dxq = *(px_head+k) - *(xq_head+k);
			dxq *= dxq;
			xq_distance += dxq;
		      }

		      if ( xq_distance < maxRange2 ) xloc.insert(make_pair(*pNHloc, px));
		    }
		    else {
		       assert(false);
		    }
		  }

		  pNHloc++;
		}

		// include itself if pNH is an attached node
		xloc.insert(make_pair(*pNH, x.find(*pNH)->second));
			    
	      } else {
		cutNeighbor.insert(*pNH);
#ifdef _M4EXTREME_DEBUG_II
		cerr << "bad old neighbor removed." << endl;
#endif
	      }

	      pNH++;
	    }
	  } 
	  else {
	    while (pNH != newNeighborEnd) {
	      const set<Set::Manifold::Point*> * nhloc = (**_nbs)(*pNH);

	      if (nhloc != NULL) {
		// the new neighbor must belong to the support region of the element
		set<Set::Manifold::Point*>::const_iterator pNHloc = nhloc->begin(), nhEnd = nhloc->end();
		while (pNHloc != nhEnd) {
		  if (cutNeighbor.find(*pNHloc) == cutNeighbor.end()
#if defined(_M4EXTREME_SHARED_INTERFACE_) // only include nodes belonging to the same body.
		                       // contact is realized by shared nodes on the interface
		      && _carrier->find(*pNHloc) != _carrier->end()
#elif defined(_M4EXTREME_MIXING_SEPERATION_) // include nodes beloning to the contact list.
                                          // for material points moving away from each other,
                                          // nodes from different body will be excluded                  
		      && _nodesInTension.find(*pNHloc) == _nodesInTension.end()
#endif
		      ) {

		    if ( x.find(*pNHloc) != x.end() ) {
		      const Set::VectorSpace::Vector & px = x.find(*pNHloc)->second;
		      xq_distance = 0.0;
		      const double *px_head = px.begin();					
		      for ( int k = 0; k < m; ++k ) {
			dxq = *(px_head+k) - *(xq_head+k);
			dxq *= dxq;
			xq_distance += dxq;
		      }

		      if ( xq_distance < maxRange2 ) xloc.insert(make_pair(*pNHloc, px));
		    }
		    // else {
		    //   assert(false);
		    // }
		  }

		  pNHloc++;
		}

		// include itself if pNH is an attached node
		xloc.insert(make_pair(*pNH, x.find(*pNH)->second));
	      } else {
		cutNeighbor.insert(*pNH);
#ifdef _M4EXTREME_DEBUG_II
		cerr << "bad old neighbor removed." << endl;
#endif
	      }

	      pNH++;
	    }
	  }

#ifdef _M4EXTREME_DEBUG_II
	  cout << "initial lambda is [" << _Base->GetLambda() << "]" << endl;
#endif

	  if ( _blacklist != NULL ) { 
	    for (nodeset_type::iterator pxloc = xloc.begin(); pxloc != xloc.end(); ) {
	      if ( _blacklist->find(pxloc->first) != _blacklist->end() ) {
		xloc.erase(pxloc++);
	      }
	      else {
		++pxloc;
	      }
	    }
	  }

	  if (xloc.size() < m + 1) { // neighbor size too small
	    needReset = true;
	  }
	  else {
	    try {
	      this->operator ()(&xloc);                    
	    } catch (...) {                        
	      needReset = true;
	    }
	  }
                    
	  if (needReset) {
	    this->enlargeSupport();
	    _maxRange *= 1.2;
	    if (_adaptiveSearch) {
	      _searchRange *= 1.2;
	    }
	  }

	} while (needReset && numReset++ < _MAX_RESET);
                
	if ( numReset >= _MAX_RESET ) {
	  // cerr << "numReset >= " << _MAX_RESET 
	  //      <<": Failed to solve the new maxent shape functions@Element::MaterialPoint::Data::Reset()" << endl;
	  // int count = 1;
	  // for ( nodeset_type::iterator pxloc = xloc.begin(); 
	  // 	pxloc != xloc.end(); pxloc++, count++ ) {
	  //   cerr << "x(" << count << ")={[" << pxloc->second << "]};" << endl;
	  // }
	  set<Set::Manifold::Point*>::const_iterator pNHloc = oldNeighbor.begin();
	  if ( _blacklist != NULL ) { 
	    for ( ; pNHloc != oldNeighbor.end(); pNHloc++ ){
	      if ( _blacklist->find(*pNHloc) == _blacklist->end() ) {
		xloc.insert( make_pair(*pNHloc, x.find(*pNHloc)->second) );
	      }	  
	      else {
		continue;
	      }	      
	    }
	  }
	  else {
	    for ( ; pNHloc != oldNeighbor.end(); pNHloc++ ){
	      xloc.insert( make_pair(*pNHloc, x.find(*pNHloc)->second) );	      	    
	    }
	  }
	  try {
	    this->operator ()(&xloc);                    
	  } catch (...) {                        
	    throw(1);
	  }
	}
      }
      else {
	if ( newNeighbor.size() < m+1 ) {
	  cerr << "support is too small:" << newNeighbor.size() << endl;
	  newNeighbor = oldNeighbor;
	}

	set<Set::Manifold::Point*>::const_iterator pNHloc = newNeighbor.begin();
	if ( _blacklist != NULL ) {
	  for ( ; pNHloc != newNeighbor.end(); pNHloc++ ){
	    if( _blacklist->find(*pNHloc) == _blacklist->end() ) {
	      xloc.insert( make_pair(*pNHloc, x.find(*pNHloc)->second) );
	    }
	    else {
	      continue;
	    }	    
	  }
	}
	else {
	  for ( ; pNHloc != newNeighbor.end(); pNHloc++ ){
	    xloc.insert( make_pair(*pNHloc, x.find(*pNHloc)->second) );
	  }
	}

	try {
	  this->operator ()(&xloc);                    
	} catch (...) {                        
	  throw(1);
	}
      }

      rxaq.clear();
      for ( nodeset_type::iterator pxloc = xloc.begin(); 
	    pxloc != xloc.end(); pxloc++ ) {
	rxaq.insert( make_pair(pxloc->first, pxloc->second - _centroid) );
      }
      
      return;
    }

    void
    Data::operator ++() {
      for ( int q = 0; q < QW.size(); ++q ) {
	_Jold[q] *= _DJ[q];
      }
    }

    void 
    Data::_check_consistency(const map<dof_type *, vector_type> & x) {
      // zero and first-order consistency check
      unsigned int dim = QP[0].size();

      for ( int i = 0; i < QP.size(); ++i ) {
	double Nsum = 0.;
	vector_type NXsum(dim);
	vector_type DNsum(dim);
	Set::VectorSpace::Hom DNXsum(dim);
	for ( map<dof_type*, vector_type>::const_iterator px = x.begin();
	      px != x.end(); px++ ) {
	  dof_type * pxloc = px->first;
	  double Nloc = N[i].find(pxloc)->second;
	  Nsum += Nloc;
	  NXsum += Nloc * px->second;
	  
	  const vector_type & DNloc = DN[i].find(pxloc)->second;
	  DNsum += DNloc;
	  DyadicSum(DNloc, px->second, DNXsum);
	}
	
	double TOL = 1.0e-5;
	double value = Nsum;
	if ( fabs(value-1.0) > TOL ) {
	  cerr << "zeroth order consistency failed for the shape function sum(Na(x))=" << value << endl;
	  //assert(false);
	}
	
	value = Norm(NXsum - QP[i]);
	if ( value > TOL ) {
	  cerr << "first order consistency failed for the shape function ||sum(Na(x)Xa) - x||=" << value << endl;
	  //assert(false);
	}
	
	value = Norm(DNsum);
	if ( value > TOL ) {
	  cerr << "zeroth order consistency failed for the derivative of the shape function ||sum(dNa/dx)||=" << value << endl;
	  //assert(false);
	}
	
	value = Norm(DNXsum);
	if ( fabs(value-1.0) > TOL ){
	  cerr << "first order consistency failed for the derivative of the shape function ||sum(dNa/dx* Xa)||=" << value << endl;
	  //assert(false);
	}
      }
    }

    void
    Data::Remesh(const map<dof_type *, vector_type> & x) {
      if ( QW.size() > 1 ) return;

      const unsigned int q = 0;
      vector_type & xq = QP[q];
      const shape_type & Nloc = N[q];
      shape_type::const_iterator pN;
      nodeset_type xloc;
      double shape = 1.0 / (double)Nloc.size();
      for (pN = Nloc.begin(); pN != Nloc.end(); pN++) {  
	nodeset_type::const_iterator px = x.find(pN->first);	  
	if ( px != x.end() ) {
	  xloc.insert( make_pair(pN->first, px->second) );
	  xq += shape * px->second; 
	}
	else {
	  cerr << "Remeshing failed" << endl;
	  throw(1);
	}
      }

      _centroid = xq;
      this->operator()(&xloc);
      
      return;
    }

    //////////////////////////////////////////////////////////////////////
    // Class LMEData
    //////////////////////////////////////////////////////////////////////
    LMEData::LMEData():_Base(NULL) {}
        
    LMEData::LMEData(const double & Rho0_,
		     const double & Beta_,
		     const double & Cutoff_,
		     const double & Ratio_,
		     const vector<double> & QW_,
		     const vector<vector_type> & QP_,
		     const map<dof_type *, vector_type> & x_,
		     bool updateNeighbor_,
		     Geometry::Search<dof_type*> ** nbs) 
      : Data(Rho0_, Cutoff_, Ratio_, QW_, QP_, x_, updateNeighbor_, nbs) {

      _Base = new Element::Interpolation::MaxEnt::Base(Beta_, x_);

      Element::Interpolation::MaxEnt::Jet < 0 > J(_Base);
      
      for (int i = 0; i < QP.size(); ++i ) {
	try {
	  J(QP[i], N[i], DN[i]);
	}
	catch(const char* code) {
	  throw code;
	}
      }

      _check_consistency(x_);

      return;
    }
        
    LMEData::~LMEData() {
      if ( _Base != NULL ) delete _Base;
    }
        
    LMEData::LMEData(const LMEData &rhs)
      : _Base(rhs._Base){}

    LMEData & LMEData::operator =(const LMEData &rhs) {

      if (this == &rhs) {
	return *this;
      }
            
      _Base = rhs._Base;
      Data::operator =(rhs);
      return *this;
    }
        
    void LMEData::operator() (const nodeset_type * xloc) {
      _Base->SetNodes(xloc);
      Element::Interpolation::MaxEnt::Jet < 0 > J(_Base);

      for (int i = 0; i < QP.size(); ++i ) {
	try {
	  J(QP[i], N[i], DN[i]);
	}
	catch(const char* code) {
	  throw code;
	}
      }

      return;
    } 
        
    void LMEData::enlargeSupport() {
      double beta = _Base->GetBeta();
      _Base->SetBeta(0.7 * beta);
      return;
    }

    //////////////////////////////////////////////////////////////////////
    // Class MLSData
    //////////////////////////////////////////////////////////////////////
    MLSData::MLSData():_Base(NULL) {}
        
    MLSData::MLSData( const double & Rho0_,
		      const int & m_,
		      const double & rho_,
		      const double & Cutoff_,
		      const double & Ratio_,
		      const vector<double> & QW_,
		      const vector<vector_type> & QP_,
		      const map<dof_type *, vector_type> & x_,
		      bool updateNeighbor_,
		      Geometry::Search<dof_type*> ** nbs) 
      : Data(Rho0_, Cutoff_, Ratio_, QW_, QP_, x_, updateNeighbor_, nbs) {

      int dim = QP[0].size();
      _Base = new Element::Interpolation::MLS::Base(dim, m_, rho_, x_);

      Element::Interpolation::MLS::Jet < 0 > J(_Base);
      for (int i = 0; i < QP.size(); ++i ) {
	try {
	  J(QP[i], N[i], DN[i]);
	}
	catch(const char* code) {
	  throw code;
	}
      }

      _check_consistency(x_);

      return;
    }

        
    MLSData::~MLSData() {
      if ( _Base != NULL ) delete _Base;
    }
        
    MLSData::MLSData(const MLSData &rhs)
      : _Base(rhs._Base){}

    MLSData & MLSData::operator =(const MLSData &rhs) {

      if (this == &rhs) {
	return *this;
      }
            
      _Base = rhs._Base;
      Data::operator =(rhs);
      return *this;
    }
        
    void MLSData::operator() (const nodeset_type * xloc) {
      _Base->SetNodes(xloc);
      Element::Interpolation::MLS::Jet < 0 > J(_Base);
      for (int i = 0; i < QP.size(); ++i ) {
	try {
	  J(QP[i], N[i], DN[i]);
	}
	catch(const char* code) {
	  throw code;
	}
      }

      return;
    } 
        
    void MLSData::enlargeSupport() {	
      return;
    }

      
    //////////////////////////////////////////////////////////////////////
    // Class MixedData
    //////////////////////////////////////////////////////////////////////
      
    MixedData::MixedData():_m(1), _rho(0.), _type(MixedData::FEA), 
			    _MLSBase(NULL), _FEAnodes(NULL), _MLSnodes(NULL) {}
      
    MixedData::MixedData(const double & Rho0_,
			 const int & m_,
			 const double & rho_,			 
			 const double & Cutoff_,
			 const double & Ratio_,
			 const vector<double> & QW_,
			 const vector<vector_type> & QP_,
			 const map<dof_type *, vector_type> & x_,
			 vector<dof_type *> FirstNeighbor_,
			 const set<dof_type *> * FEAnodes_,
			 const set<dof_type *> * MLSnodes_,
			 bool updateNeighbor_,
			 Geometry::Search<dof_type*> ** nbs) 
      : _FEAnodes(FEAnodes_), _MLSnodes(MLSnodes_), 
	_FirstNeighbor(FirstNeighbor_), _MLSBase(NULL),
	_type(MixedData::FEA), _m(m_), _rho(rho_),
	Data(Rho0_, Cutoff_, Ratio_, QW_, QP_, x_, updateNeighbor_, nbs) {

      unsigned int count = 0;
      map<dof_type*, vector_type>::const_iterator px;
      for ( px = x_.begin(); px != x_.end(); px++ ) {
	if ( FEAnodes_->find(px->first) == FEAnodes_->end() ) {
	  _type = MIXED;
	  count++;
	}
      }

      if ( count == x_.size() ) {
	_type = MLS;
      }

      if ( _type != FEA ) {
	_MLSBase = new Element::Interpolation::MLS::Base(QP[0].size(), m_, rho_, x_);
      }

      compute(x_);
    }
        
    MixedData::~MixedData() {
       if ( _MLSBase != NULL ) delete _MLSBase;
    }
      
    MixedData::MixedData(const MixedData &rhs)
      :_m(rhs._m), _rho(rhs._rho), _type(rhs._type), 
       _FirstNeighbor(rhs._FirstNeighbor), _MLSBase(rhs._MLSBase),
       _FEAnodes(rhs._FEAnodes), _MLSnodes(rhs._MLSnodes){}
      
    MixedData & MixedData::operator =(const MixedData &rhs) {

      if (this == &rhs) {
	return *this;
      }
	
      _m = rhs._m;
      _rho = rhs._rho;
      _type = rhs._type;
      _FirstNeighbor = rhs._FirstNeighbor;
      _MLSBase = rhs._MLSBase;
      _FEAnodes = rhs._FEAnodes;
      _MLSnodes = rhs._MLSnodes;
	
      Data::operator =(rhs);
      return *this;
    }
        
    void MixedData::operator() (const nodeset_type * xloc) {
      unsigned int count = 0;
      nodeset_type::const_iterator px;
      for ( px = xloc->begin(); px != xloc->end(); px++ ) {
	if ( _FEAnodes->find(px->first) == _FEAnodes->end() ) {
	  _type = MIXED;
	  count++;
	}
      }

      if ( count == xloc->size() ) {
	_type = MLS;
      }

      if ( _type != FEA ) {
	if ( _MLSBase != NULL ) {
	  _MLSBase->SetNodes(xloc);
	}
	else {
	  _MLSBase = new Element::Interpolation::MLS::Base(QP[0].size(), _m, _rho, *xloc);
	}
      }

      compute(*xloc);
    }

    void MixedData::compute(const nodeset_type & x) {
      typedef Set::VectorSpace::Hom hom_type;
      map<dof_type *, vector_type>::const_iterator px;
	
      int numofMPT = QP.size();

      switch(_type) {
      default:
	cerr << "TYPE is not defined @MixedData" << endl;
	assert(false);
      case FEA:
	{
	  for ( int i = 0; i < numofMPT; ++i ) {
	      computeFE(x, QP[i], N[i], DN[i]);
	  }
	}
	break;
      case MLS:
	{
	  Element::Interpolation::MLS::Jet < 0 > JMLS(_MLSBase);
	  for ( int i = 0; i < numofMPT; ++i ) {
	    try {
	      JMLS(QP[i], N[i], DN[i]);
	    }
	    catch(...) {
	      throw(1);
	    }
	  }
	}
	break;
      case MIXED:
	{
	  int dim = QP[0].size();
	  Element::Interpolation::MLS::Jet < 0 > JMLS(_MLSBase);
	  for ( int i = 0; i < numofMPT; ++i ) {
	    const vector_type & xq = QP[i];
	    map<dof_type *, double> FEAN; 
	    map<dof_type *, vector_type> FEADN;
	    computeFE(x, xq, FEAN, FEADN);
	    
	    map<dof_type *, double> MLSN; 
	    map<dof_type *, vector_type> MLSDN;
	    map<dof_type *, vector_type> AIPW;
	    map<dof_type *, hom_type> AIxPW;
	    map<dof_type *, hom_type> AIPWx;
	    try {
	      JMLS(xq, MLSN, MLSDN, AIPW, AIxPW, AIPWx);
	    }
	    catch(const char* code) {
	      throw code;
	    }
	    
	    // assemble FEAN and MLSN
	    N[i].clear(); DN[i].clear();
	    
	    int k = _MLSBase->GetNumberOfTerms();
	    vector_type Pi(k);
	    map<dof_type *, vector_type> Pbs;	
	    for ( px = x.begin(); px != x.end(); px++ ) {
	      if ( _FEAnodes->find(px->first) != _FEAnodes->end() ) {
		monomials(px->second, Pi);
		Pbs.insert( make_pair(px->first, Pi) );
	      }
	    }
	    
	    vector_type NP(k);
	    hom_type DNP(dim, k);
	    for ( map<dof_type *, vector_type>::const_iterator pNode = Pbs.begin();
		  pNode != Pbs.end(); pNode++ ) {
	      NP += FEAN.find(pNode->first)->second * pNode->second;
	      DNP += Dyadic(pNode->second, FEADN.find(pNode->first)->second);
	    }
	    
	    //loop over all the nodes
	    for ( px = x.begin(); px != x.end(); px++ ) {
	      dof_type * pxloc = px->first;
	      if ( _MLSnodes->find(pxloc) != _MLSnodes->end() ) {
		double Nloc = MLSN.find(pxloc)->second;
		Nloc -= NP(AIPW.find(pxloc)->second);
		
		vector_type DNloc = MLSDN.find(pxloc)->second;
		DNloc -= DNP*AIPW.find(pxloc)->second;
		DNloc -= (AIxPW.find(pxloc)->second + AIPWx.find(pxloc)->second)*NP;
		
		if ( _FEAnodes->find(pxloc) != _FEAnodes->end() ) {
		  N[i].insert ( make_pair(pxloc, Nloc + FEAN.find(pxloc)->second) );
		  DN[i].insert ( make_pair(pxloc, DNloc + FEADN.find(pxloc)->second) );
		}
		else {
		  N[i].insert ( make_pair(pxloc, Nloc) );
		  DN[i].insert ( make_pair(pxloc, DNloc) );
		}
	      }
	      else {
		N[i].insert( make_pair(pxloc, FEAN.find(pxloc)->second) );
		DN[i].insert( make_pair(pxloc, FEADN.find(pxloc)->second) );
	      }
	    }
	  }	  
	}

	break;
      }

      _check_consistency(x);

      return;
    }
      
    void MixedData::computeFE(const nodeset_type & x, 
			      const vector_type &xq,
			      map<dof_type *, double> & FEAN,
			      map<dof_type *, vector_type> & FEADN ) const {
      int dim = xq.size();
      int npp = dim+1;
      assert(_FirstNeighbor.size() == npp );

      // regularize the data from the first neighbor for further calculations
      map<dof_type*, int> index;
      vector<double> xa, ya, za;
      for ( int i = 0; i< _FirstNeighbor.size(); ++i ) {
	dof_type * pE = _FirstNeighbor[i];
	const vector_type & xloc = x.find(pE)->second;
	
	if ( _FEAnodes->find(pE) != _FEAnodes->end() ) {
	  index.insert( make_pair(pE, i) );
	}

	xa.push_back(xloc[0]);
	if ( dim > 1 ) ya.push_back(xloc[1]);
	if ( dim > 2 ) za.push_back(xloc[2]);
      }

      nodeset_type::const_iterator px;
      map<dof_type*, int>::const_iterator pI;
      if ( index.empty() ) {
	for ( px = x.begin(); px != x.end(); px++ ) {
	  FEAN.insert( make_pair(px->first, 0.0) );
	  FEADN.insert( make_pair(px->first, vector_type(dim)) );
	}
      }
      else {
	switch(dim) {
	case 1:
	  break;
	case 2:
	  {
	    double Ae = xa[1]*ya[2] - xa[2]*ya[1] - xa[0]*ya[2] + xa[0]*ya[1] + ya[0] * xa[2] - ya[0] * xa[1];
	    double Nloc[npp];
	    Nloc[0] = (xa[1]*ya[2] - xa[2]*ya[1] + ya[1]*xq[0] - ya[2]*xq[0] + xa[2]*xq[1] - xa[1]*xq[1]) / Ae;
	    Nloc[1] = (xa[2]*ya[0] - xa[0]*ya[2] + ya[2]*xq[0] - ya[0]*xq[0] + xa[0]*xq[1] - xa[2]*xq[1]) / Ae;
	    Nloc[2] = (xa[0]*ya[1] - xa[1]*ya[0] + ya[0]*xq[0] - ya[1]*xq[0] + xa[1]*xq[1] - xa[0]*xq[1]) / Ae;
	    
	    vector<vector_type> DNloc(npp, vector_type(dim));
	    DNloc[0][0] = (ya[1] - ya[2]) / Ae; DNloc[0][1] = (xa[2] - xa[1]) / Ae;
	    DNloc[1][0] = (ya[2] - ya[0]) / Ae; DNloc[1][1] = (xa[0] - xa[2]) / Ae;
	    DNloc[2][0] = (ya[0] - ya[1]) / Ae; DNloc[2][1] = (xa[1] - xa[0]) / Ae;

	    for ( px = x.begin(); px != x.end(); px++ ) {
	      if ( _FEAnodes->find(px->first) != _FEAnodes->end() &&
		   (pI = index.find(px->first)) != index.end() ) {
		int k = pI->second;
		FEAN.insert( make_pair(px->first, Nloc[k]) );
		FEADN.insert( make_pair(px->first, DNloc[k]) );
	      }
	      else {
		FEAN.insert( make_pair(px->first, 0.0) );
		FEADN.insert( make_pair(px->first, vector_type(dim)) );
	      }
	    }
	  }
	  break;
	case 3:
	  break;
	default:
	  assert(false);
	}
      }
      return;
    }

    void MixedData::monomials(const vector_type & x, vector_type & P) const {
      int dim = _MLSBase->GetDim();
      int m   = _MLSBase->GetOrder();
      int k   = _MLSBase->GetNumberOfTerms();
      assert( P.size() == k && x.size() == dim);
	
      // 1D
      if (dim == 1) {
	double xloc = x[0];
	P[0] = 1;
	for (int i=1;i<k;i++) {
	  P[i] = P[i-1]*xloc;
	}
      }
	
      // 2D
      if (dim == 2) {
	double xloc = x[0];
	double yloc = x[1];	    
	double tab[m+1][m+1];
	// define 1st line and first column
	tab[0][0] = 1.;
	for (int i=1;i<m+1;i++) {
	  tab[0][i] = tab[0][i-1]*xloc;
	  tab[i][0] = tab[i-1][0]*yloc;
	}
			
	if (m>1){
	  for (int i=1;i<m+1;i++) {
	    for (int j=1;j<m+1;j++) {
	      tab[i][j] = tab[0][j]*tab[i][0];
	    }
	  }
	}
			
	int inc = 0;
	for (int i=0;i<m+1;i++) {
	  for (int j=0;j<i+1;j++) {
	    P[inc] = tab[j][i-j];
	    inc++;
	  }
	}
      }

      // 3D
      if (dim == 3) {			
	if (m==1){
	  P[0] = 1;
	  P[1] = x[0];
	  P[2] = x[1];
	  P[3] = x[2];
	}
	if (m == 2) {
	  double xloc = x[0];
	  double yloc = x[1];
	  double zloc = x[2];
	    
	  P[0] = 1;
	  P[1] = xloc;
	  P[2] = yloc;
	  P[3] = zloc;
	  P[4] = xloc*xloc;
	  P[5] = xloc*yloc;
	  P[6] = yloc*yloc;
	  P[7] = yloc*zloc;
	  P[8] = zloc*zloc;
	  P[9] = xloc*zloc;
	}
      }
    }
	  
    //////////////////////////////////////////////////////////////////////
    // Class LocalState
    //////////////////////////////////////////////////////////////////////

    LocalState::LocalState(): D(NULL) {
    }

    LocalState::~LocalState() {
    }

    Element::LocalState *
    LocalState::Clone() const {
      return new LocalState(*this);
    }

    LocalState::LocalState(Data *D_, const vector<Material::LocalState *> & MLS_) : 
      hourglass_modulus(0.0),
      POld(D_->_centroid.size()),
      friction_coefficient(0.0), Normal(D_->GetCentroid().size()), D(D_), MLS(MLS_){      
      for ( int q = 0; q < MLS.size(); ++q ) {
	SourceLS.push_back(0);
      }
    }
    
    LocalState::LocalState(Data *D_, const vector<Material::LocalState *> & MLS_,
			   const vector<Material::LocalState *> & SourceLS_) : 
      hourglass_modulus(0.0),
      POld(D_->_centroid.size()),
      friction_coefficient(0.0), Normal(D_->GetCentroid().size()), D(D_),
      MLS(MLS_), SourceLS(SourceLS_){
    }

    LocalState::LocalState(const LocalState &rhs) :
      hourglass_modulus(rhs.hourglass_modulus),
      POld(rhs.POld),
      friction_coefficient(rhs.friction_coefficient), Normal(rhs.Normal), D(rhs.D),
      MLS(rhs.MLS), SourceLS(rhs.SourceLS) {}
    
    LocalState &
    LocalState::operator =(const LocalState &rhs) {
      if (this == &rhs) return *this;
      hourglass_modulus = rhs.hourglass_modulus;
      POld = rhs.POld;
      friction_coefficient = rhs.friction_coefficient;
      Normal = rhs.Normal;
      D = rhs.D;
      MLS = rhs.MLS;
      SourceLS = rhs.SourceLS;
      return *this;
    }

    LocalState::range_type
    LocalState::operator () ( const LocalState::domain_type &x ) const {

      const vector<dshape_type> & DN = D->DN;
      unsigned int m = x.begin()->second.size();
      unsigned int n = DN[0].begin()->second.size();
      assert(m == n);

      range_type F;
      for ( int q = 0; q < DN.size(); ++q ) {
	Set::VectorSpace::Vector Dx(m * n);	            
	dshape_type::const_iterator pDN;
	for (pDN = DN[q].begin(); pDN != DN[q].end(); pDN++) {
	  Dx += Dyadic(pDN->second, x.find(pDN->first)->second);
	}

	F.push_back(Dx);
      }

      return F;
    }

    void LocalState::operator () ( const LocalState::domain_type &x,
				   LocalState::range_type & F ) const {
      const vector<dshape_type> & DN = D->DN;
      unsigned int m = x.begin()->second.size();
      unsigned int n = DN[0].begin()->second.size();
      assert(m == n);

      for ( int q = 0; q < DN.size(); ++q ) {
	Set::VectorSpace::Hom Dx(m, n);
	dshape_type::const_iterator pDN;
	for ( pDN = DN[q].begin(); pDN != DN[q].end(); pDN++ ) {
	  DyadicSum(pDN->second, x.find(pDN->first)->second, Dx);
	  //Dx += Dyadic(pDN->second, x.find(pDN->first)->second);
	}
	
	F.push_back(Dx);	

#ifdef _M4EXTREME_DEBUG_II
	std::cout<<"\nElement::MaterialPoint::LocalState computes the deformation gradient with"<<std::endl;
	for (pDN = D->DN.begin(); pDN != D->DN.end(); pDN++) {
	  std::cout<<"DN=["<<pDN->second<<"]"<<std::endl;
	}
#endif
      }
      
      //D->GetJnew() = Jacobian(F);
      
      return;
    }

    void LocalState::GetDisplacement( const LocalState::domain_type &x,
				      vector<Set::VectorSpace::Vector> & Phi ) const {
      const vector<shape_type> & N = D->N;
      unsigned int m = x.begin()->second.size();
      
      for ( int q = 0; q < N.size(); ++q ) {
	Set::VectorSpace::Vector xq(m);
	shape_type::const_iterator pN = N[q].begin(), Nend = N[q].end();
	for ( ; pN != Nend; pN++ ) {
	  xq +=  pN->second * x.find(pN->first)->second;
	}

	Phi.push_back(xq);
      }
      
      return;
    }    
    
    void
    LocalState::operator ++() {
      ++(*D);
      for ( int i = 0; i < MLS.size(); ++i ) {
	if ( MLS[i] != NULL ) MLS[i]->operator++();
	if ( SourceLS[i] != NULL ) SourceLS[i]->operator++();
      }
    }
    
    set<Set::Manifold::Point *>
    LocalState::GetNodes() const {
      set<Set::Manifold::Point *> Nodes;
      const dshape_type & DN = D->DN[0];
      dshape_type::const_iterator pDN;
      for (pDN = DN.begin(); pDN != DN.end(); pDN++) {
	Nodes.insert(pDN->first);
      }

      return Nodes;
    }

    set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >
    LocalState::GetNodePairs() const {
      set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NodePairs;
      set<Set::Manifold::Point *> Nodes;
      GetNodes(Nodes);

      set<Set::Manifold::Point *>::iterator pM;
      set<Set::Manifold::Point *>::iterator pN;
      for (pM = Nodes.begin(); pM != Nodes.end(); pM++)
	for (pN = Nodes.begin(); pN != Nodes.end(); pN++)
	  NodePairs.insert(make_pair(*pM, *pN));
            
      return NodePairs;
    }

    void LocalState::GetNodes(set<Set::Manifold::Point *> & Nodes) const {
      if (!Nodes.empty()) {
	Nodes.clear();
      }
      const dshape_type & DN = D->DN[0];
      dshape_type::const_iterator pDN;
      for (pDN = DN.begin(); pDN != DN.end(); pDN++) {
	Nodes.insert(pDN->first);
      }

      return;
    }

    void
    LocalState::GetNodePairs(set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > & NodePairs) const {
      if ( !NodePairs.empty() ) {
	NodePairs.clear();
      }

      set<Set::Manifold::Point *> Nodes;
      GetNodes(Nodes);

      set<Set::Manifold::Point *>::iterator pM;
      set<Set::Manifold::Point *>::iterator pN;
      for (pM = Nodes.begin(); pM != Nodes.end(); pM++)
	for (pN = Nodes.begin(); pN != Nodes.end(); pN++)
	  NodePairs.insert(make_pair(*pM, *pN));
            
      return;
    }

    void
    LocalState::Reset(
		      const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) {

      try {
	D->Reset(x);
      }
      catch(...) {
	_isActivated = false;
	return;

	//cerr << "\nFailure point @ Element::MaterialPoint::LocalState::Reset()"<<endl;
	//assert(false);
      }

    }

    //////////////////////////////////////////////////////////////////////
    // Class Energy<0>
    //////////////////////////////////////////////////////////////////////

    Energy < 0 > ::Energy() {
    }

    Energy < 0 > ::~Energy() {
    }

    Element::Energy < 0 > *
    Energy < 0 > ::Clone() const {
      return new Energy < 0 > (*this);
    }

    Element::LocalState *
    Energy < 0 > ::GetLocalState() const {
      return LS;
    }

    const vector<Material::Energy<0> *> &
    Energy< 0 >::GetW() const {
      return W;
    }

    Energy < 0 > ::Energy(LocalState *LS_,
			  const vector<Material::Energy < 0 > *> & W_) : LS(LS_), W(W_) {
      for ( int q = 0; q < W.size(); ++q ) {
	SourceW.push_back(0);
      }
    }
    
    Energy < 0 > ::Energy(LocalState *LS_,
			  const vector<Material::Energy < 0 > *> & W_,
			  const vector<Material::Energy < 0 > *> & SourceW_) :
      LS(LS_), W(W_), SourceW(SourceW_) {
    }
    
    Energy < 0 > ::Energy(const Energy < 0 > &rhs) :
      LS(rhs.LS), W(rhs.W), SourceW(rhs.SourceW) {
    }

    Energy < 0 > &
    Energy < 0 > ::operator =(const Energy < 0 > &rhs) {
      if (this == &rhs) return *this;
      LS = rhs.LS;
      W = rhs.W;
      SourceW = rhs.SourceW;
      return *this;
    }

    double
    Energy < 0 > ::operator () (
				const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const {

      vector<Set::VectorSpace::Vector> Dx;
      vector<Set::VectorSpace::Vector> Phi;
      
      LS->operator()(x, Dx);
      if ( SourceW[0] != NULL ) {
	LS->GetDisplacement(x, Phi);
      }

      unsigned int m = x.begin()->second.size();
      unsigned int n = (LS->D->DN)[0].begin()->second.size();

      double E = 0.0;
      const vector<double> & QW = LS->D->QW;
      const vector<double> & Mass = LS->D->Mass;
      
      for ( int q = 0; q < QW.size(); ++q ) {
	Set::VectorSpace::Hom F(m, n, Dx[q].begin());
	LS->D->_DJ[q] = Jacobian(F);

	if ( SourceW[q] != NULL ) {
	  E += QW[q] * W[q]->operator()(Dx[q]) + Mass[q] * SourceW[q]->operator()(Phi[q]);
	}
	else {
	  E += QW[q] * W[q]->operator()(Dx[q]);
	}
      }

      return E;

    }

    //////////////////////////////////////////////////////////////////////
    // Class Energy<1>
    //////////////////////////////////////////////////////////////////////

    Energy < 1 > ::Energy() {
    }

    Energy < 1 > ::~Energy() {
    }

    Element::Energy < 1 > *
    Energy < 1 > ::Clone() const {
      return new Energy < 1 > (*this);
    }

    inline
    Element::LocalState *
    Energy < 1 > ::GetLocalState() const {
      return LS;
    }

    const vector<Material::Energy<1> *> &
    Energy< 1 >::GetDW() const {
      return DW;
    }

    Energy < 1 > ::Energy(LocalState *LS_,
			  const vector<Material::Energy < 1 > *> & DW_) : LS(LS_), DW(DW_) {
      Element::Energy<1>::_ELS = LS;
      for ( int q = 0; q < DW.size(); ++q ) {
	SourceDW.push_back(0);
      }      
    }
    
    Energy < 1 > ::Energy(LocalState *LS_,
			  const vector<Material::Energy < 1 > *> & DW_,
			  const vector<Material::Energy < 1 > *> & SourceDW_) :
      LS(LS_), DW(DW_), SourceDW(SourceDW_) {
      Element::Energy<1>::_ELS = LS;
    }
    
    Energy < 1 > ::Energy(const Energy < 1 > &rhs)
      : LS(rhs.LS), DW(rhs.DW), SourceDW(rhs.SourceDW) {
    }

    Energy < 1 > &
    Energy < 1 > ::operator =(const Energy < 1 > &rhs) {
      if (this == &rhs) return *this;
      LS = rhs.LS;
      DW = rhs.DW;
      SourceDW = rhs.SourceDW;
      return *this;
    }

    map<Set::Manifold::Point *, Set::VectorSpace::Vector>
    Energy < 1 > ::operator () (
				const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const {
      map<Set::Manifold::Point *, Set::VectorSpace::Vector> DE;
      map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
      
      unsigned int m = x.begin()->second.size();
      unsigned int n = (LS->D->DN)[0].begin()->second.size();
      
      assert(m < 4 && n < 4);
           
      for (px = x.begin(); px != x.end(); px++) {
	assert(m == px->second.size());
	DE.insert(range_type::value_type(px->first,
					 Set::VectorSpace::Vector(px->second.size())));
      }
      
      vector<Set::VectorSpace::Vector> Dx;
      Dx = (*LS)(x);           
      
      const vector<double> & QW = LS->D->QW;
      const vector<double> & Mass = LS->D->Mass;
	    
      for ( int q = 0; q < QW.size(); ++q ) {
	Set::VectorSpace::Hom P(m, n);
	Set::VectorSpace::Hom F(m, n, Dx[q].begin());
	
	// check jocobian
	double J = Jacobian(F);
	J *= LS->D->_Jold[q];
	LS->D->_DJ[q] = J;
	if ( J < 1.0e-3 || J > 1.0e3 || J != J ) {
	  LS->_isActivated = false;                  
	}
	else {
	  try {
	    P = (*DW[q])(F);
	  }
	  catch(...) {
	    LS->_isActivated = false;
	    Null(P);
	  }	  
	}
	
#ifdef _M4EXTREME_DEBUG_III
	std::cout << "\nF=[" << F << "]\tP=[" << P << "]" << std::endl;
	
	//        double xq = _gxqs0.find(LS)->second;
	//        std::cout<<"\nt = "<<_gT<<" exact solution at ["<<xq<<"]-->["<<_gxsol(xq, _gT) <<"] && ["<<(LS->D->GetQP())[0] <<"] is:\nF="<<_gfsol(xq, _gT)
	//                <<"\tFdot="<<_gfdotsol(xq, _gT)
	//                <<"\tP="<<_gpsol(xq, _gT)<<std::endl;
#endif
	
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DN = (LS->D->DN)[q];
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pDN;
	
#if defined(_M4EXTREME_MIXING_SEPERATION_)
	// update the contact list
	if (LS->D->_carrier != NULL) {
	  if (P(F) > 0.0) {
	    for (pDN = DN.begin(); pDN != DN.end(); pDN++) {
	      if (LS->D->_carrier->find(pDN->first) == LS->D->_carrier->end()) {
		LS->D->_nodesInTension.insert(pDN->first);
	      }
	    }
	  } else {
	    if (!LS->D->_nodesInTension.empty()) {
	      LS->D->_nodesInTension.clear();
	    }
	  }
	}
	else {
	  assert("body can't be empty");
	}
#endif
        
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDE;
	for (pDE = DE.begin(); pDE != DE.end(); pDE++) {
	  if ( (pDN = DN.find(pDE->first)) != DN.end() ) {
	    pDE->second += QW[q] * P * pDN->second;
	  }
	  else {
	    cerr << "cannot find the derivative of the shape function defined at the node" \
	      "@Element::MaterialPoint::Energy<1>" << endl;
	    assert(false);
	  }
	}
	
      }
      
      return DE;
    }

    void Energy < 1 > ::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x,
				     map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DE) const {

      if ( __builtin_expect(!DE.empty(), 1) ) {
	DE.clear();
      }

      map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
      unsigned int m = x.begin()->second.size();
      unsigned int n = (LS->D->DN)[0].begin()->second.size();

      assert(m < 4 && n < 4 );

      // for stress
      vector<Set::VectorSpace::Vector> Dx;
      LS->operator()(x, Dx);

      // for body forces
      vector<Set::VectorSpace::Vector> Phi;
      if ( SourceDW[0] != NULL ) {
	LS->GetDisplacement(x, Phi);
      }
	
      const vector<double> & QW = LS->D->QW;
      const vector<double> & Mass = LS->D->Mass;
	    
      for ( int q = 0; q < QW.size(); ++q ) {
	Set::VectorSpace::Hom P(m, n);
	Set::VectorSpace::Hom F(m, n, Dx[q].begin());

	// check jocobian
	double J = Jacobian(F);
	J *= LS->D->_DJ[q];
	LS->D->_DJ[q] = J;
	if ( J < 0.001 || J > 1000. || J != J ) {
	  
	  LS->_isActivated = false;
	  
#ifdef _M4EXTREME_DEBUG_II
	  cerr << "Material Point got bad Jacobian value J=" << J << " @xq=[" << LS->D->QP[q] << "];" << endl
	       << "qw="  << QW[q] <<";" << endl
	       << "Search range: " << LS->D->GetMaxRange() << endl
	       << "F=[" << F << "];" << endl;		   
	  
	  int count=1;
	  for ( map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px = x.begin();
	  	px != x.end(); px++, count++) {
	    cerr << "x(" << count << ")={[" << px->second << "]};" << endl
	  	 << "N(" << count << ")="  << (LS->D->N)[q].find(px->first)->second << ";" << endl
	  	 << "DN(" << count << ")={[" << (LS->D->DN)[q].find(px->first)->second << "]};" << endl;
	  } 
#endif
	  
	  Null(P);
	  LS->_isActivated = false;
	}
	else {
	  try {
	    P = (*DW[q])(F);
	    LS->POld = P;
	  }
	  catch(...) {
	    P = LS->POld;
	    LS->_isActivated = false;  
	    
#ifdef _M4EXTREME_DEBUG_II
	    cout << "\nFailure point @ Element::MaterialPoint::Energy<1>::Operator():"<<endl
		 << "material point: ["<<(LS->D->QP)[q]<<"]"<<endl
		 << "search range for the first neighbors is " << LS->D->_maxRange << endl
		 << "volume: " << QW[q] << endl;
	    
	    for ( map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px = x.begin();
		  px != x.end(); px++ ) {
	      cerr << "[" << px->second << "]\t" << (LS->D->N)[q].find(px->first)->second 
		   << "\t" << (LS->D->DN)[q].find(px->first)->second << endl;
	    }
	    //exit(1);
#endif
	  }	    
	}
            
	if ( 0.1 < J && J < 10.0 ) {	  
	  Set::VectorSpace::Hom E = 0.5 * Adjoint(F) * F;	  
	  switch ( m ) {
	  case 1:
	    break;
	  case 2:
	    {
	      E[0][0] -= 0.5; E[1][1] -= 0.5;
	      double I1 = E[0][0] + E[1][1];
	      double I2 = E[0][0] * E[1][1] - E[0][1] * E[0][1];
	      double c = sqrt(I1*I1 - 4.0 * I2);
	      double e1 = I1; 
	      e1 += c; 
	      e1 *= 0.5;
	      double e2 = I1;
	      e2 -= c;
	      e2 *= 0.5;
	      double maxEps = e1 > e2 ? e1 : e2;
	      
	      if ( maxEps > 0. )
		LS->D->_maxRange *= (1.0 + 0.1 * maxEps);	      
	    }
	    break;
	  case 3:
	    {		
	      E[0][0] -= 0.5; E[1][1] -= 0.5; E[2][2] -= 0.5; // E=0.5*(C-I)		
	      
	      double I1 = E[0][0] + E[1][1] + E[2][2];
	      double I2 = E[0][0] * E[1][1] - E[0][1] * E[0][1] 
		+ E[0][0] * E[2][2] - E[0][2] * E[0][2]
		+ E[1][1] * E[2][2] - E[1][2] * E[1][2];
	      double I3 = Jacobian(E);
	      
	      double c = I1*I1 - 3*I2;
	      
	      if ( fabs(c) > 1.0e-16 ) {                
		
		double alpha = acos(0.5 * (2 * I1 * I1 * I1 - 9 * I1 * I2 + 27 * I3) / pow(c, 1.5)) / 3.0;
		
		double sqrtTerm = 2.0 * sqrt(c) / 3.0;
		double e1 = I1 / 3.0 + sqrtTerm * cos(alpha);
		double e2 = I1 / 3.0 + sqrtTerm * cos(alpha + 2.0943951);
		double e3 = I1 / 3.0 + sqrtTerm * cos(alpha + 4.18879);
		double maxEps = e1 > e2 ? e1 : e2;
		maxEps = maxEps > e3 ? maxEps : e3;
		if ( maxEps > 0. )
		  LS->D->_maxRange *= (1.0 + 0.1 * maxEps);
	      }
	    }
	    break;
	  default:
	    break;
	  }
	}
	
#ifdef _M4EXTREME_DEBUG_II
	std::cout << "\nF=[" << F << "]\tP=[" << P << "]" << std::endl;
#endif

	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DN = (LS->D->DN)[q];
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pDN;

	const map<Set::Manifold::Point *, double> & N = (LS->D->N)[q];
	map<Set::Manifold::Point *, double>::const_iterator pN;
	
#if defined(_M4EXTREME_MIXING_SEPERATION_)         
	// update the contact list
	if (LS->D->_carrier != NULL) {
	  if (P(F) > 0.0) {
	    for (pDN = DN.begin(); pDN != DN.end(); pDN++) {
	      if (LS->D->_carrier->find(pDN->first) == LS->D->_carrier->end()) {
		LS->D->_nodesInTension.insert(pDN->first);
	      }
	    }
	  } else {
	    if (!LS->D->_nodesInTension.empty()) {
	      LS->D->_nodesInTension.clear();
	    }
	  }
	}
#endif

	// total nodal forces
	Set::VectorSpace::Vector force(m);
	for (px = x.begin(); px != x.end(); px++) {
	  pDN = DN.find(px->first);
	  if ( pDN != DN.end() ) {
	    if ( SourceDW[q] != NULL ) {
	      pN = N.find(px->first);
	      DE.insert( range_type::value_type(px->first,
						P(QW[q], pDN->second, force)
						+ SourceDW[q]->operator()(Phi[q])*(Mass[q]*pN->second)) );
	    }
	    else {
	      DE.insert( range_type::value_type(px->first,
						P(QW[q], pDN->second, force)) );
	    }
	  }
	  else {
	    Null(force);
	    DE.insert( range_type::value_type(px->first, force) );
	  }
	}

        if (LS->hourglass_modulus > 0.0) {
            // input: predicted x_{a, k+1}
            // 0) Determine the new location of the material point x_{p,k+1} = x_{a,k+1}N_a(x_{p, k})
            // 1) Compute the true radial vectors at t_{k+1}: R_{a,k+1} = x_{a,k+1} - x_{p,k+1}
            // 2) Compute the linearly transformed radial vectors: Rt_{a,k+1} = F_{k->k+1}R_{a,k}
            //    where R_{a,k} is stored at time t_k, and
            //    F_{k->k+1} = x_{a,k+1} \tensor \dN_a(x_{p,k})/dx
            // 3) Calculate the error: Err_{a,k+1} = (R_{a,k+1} - Rt_{a,k+1})/||R_{a,k}||
            // 4) Compute the correction force: f_a = \epsilon * sum_p {Err_{a,k+1}N_a(x_p)}
            Set::VectorSpace::Vector xqnew(m);
            for (px = x.begin(); px != x.end(); px++) {
              xqnew += N.find(px->first)->second * px->second;
            }

            const map<Set::Manifold::Point*, Set::VectorSpace::Vector> & RX = LS->D->rxaq;
            map<Set::Manifold::Point*, Set::VectorSpace::Vector>::const_iterator pRX;
            for (px = x.begin(); px != x.end(); px++) {
              pRX = RX.find(px->first);
              if ( pRX != RX.end() ) {
                Set::VectorSpace::Vector error = px->second - xqnew;
                double rxold = Norm(pRX->second);
                error -= F(pRX->second);
                error *= LS->hourglass_modulus * N.find(px->first)->second / rxold;
                Set::VectorSpace::Vector & floc = DE.find(px->first)->second;
                floc -= error;		
              }
            }
        }
      }

      return;
    }


    //////////////////////////////////////////////////////////////////////
    // Class Energy<2>
    //////////////////////////////////////////////////////////////////////

    Energy < 2 > ::Energy() {
    }

    Energy < 2 > ::~Energy() {
    }

    Element::Energy < 2 > *
    Energy < 2 > ::Clone() const {
      return new Energy < 2 > (*this);
    }

    Element::LocalState *
    Energy < 2 > ::GetLocalState() const {
      return LS;
    }

    const vector<Material::Energy<2> *> &
    Energy< 2 >::GetDDW() const {
      return DDW;
    }

    Energy < 2 > ::Energy(LocalState *LS_,
			  const vector<Material::Energy < 2 > *> & DDW_) : LS(LS_), DDW(DDW_) {
      for ( int q = 0; q < DDW.size(); ++q ) {
	SourceDDW.push_back(0);
      }
    }
    
    Energy < 2 > ::Energy(LocalState *LS_,
			  const vector<Material::Energy < 2 > *> & DDW_,
			  const vector<Material::Energy < 2 > *> & SourceDDW_) :
      LS(LS_), DDW(DDW_), SourceDDW(SourceDDW_) {
    }
    
    Energy < 2 > ::Energy(const Energy < 2 > &rhs) :
      LS(rhs.LS), DDW(rhs.DDW), SourceDDW(rhs.SourceDDW) {
    }

    Energy < 2 > &
    Energy < 2 > ::operator =(const Energy < 2 > &rhs) {
      if (this == &rhs) return *this;
      LS = rhs.LS;
      DDW = rhs.DDW;
      SourceDDW = rhs.SourceDDW;
      return *this;
    }

    map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom>
    Energy < 2 > ::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const {
      assert(false);
    }

    //////////////////////////////////////////////////////////////////////
    // Class Jet<0>
    //////////////////////////////////////////////////////////////////////

    Jet < 0 > ::Jet() {
    }

    Jet < 0 > ::~Jet() {
    }

    Element::Jet < 0 > *
    Jet < 0 > ::Clone() const {
      return new Jet < 0 > (*this);
    }

    Element::LocalState *
    Jet < 0 > ::GetLocalState() const {
      return LS;
    }

    Jet < 0 > ::Jet(LocalState *LS_,
		    Material::Jet < 0 > *J_) : LS(LS_), J(J_) {
      SourceJ = 0;
    }
    
    Jet < 0 > ::Jet(LocalState *LS_,
		    Material::Jet < 0 > *J_,
		    Material::Jet < 0 > *SourceJ_) :
      LS(LS_), J(J_), SourceJ(SourceJ_) {
    }
    
    Jet < 0 > ::Jet(const Jet < 0 > &rhs) : LS(rhs.LS), J(rhs.J), SourceJ(rhs.SourceJ) {
    }

    Jet < 0 > &
    Jet < 0 > ::operator =(const Jet < 0 > &rhs) {
      if (this == &rhs) return *this;
      LS = rhs.LS;
      J = rhs.J;
      SourceJ = rhs.SourceJ;
      return *this;
    }

    pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> >
    Jet < 0 > ::operator () (
			     const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const {
      assert(false);
    }

    //////////////////////////////////////////////////////////////////////
    // Class Jet<1>
    //////////////////////////////////////////////////////////////////////

    Jet < 1 > ::Jet() {
    }

    Jet < 1 > ::~Jet() {
    }

    Element::Jet < 1 > *
    Jet < 1 > ::Clone() const {
      return new Jet < 1 > (*this);
    }

    Element::LocalState *
    Jet < 1 > ::GetLocalState() const {
      return LS;
    }

    Jet < 1 > ::Jet(LocalState *LS_,
		    Material::Jet < 1 > *DJ_) : LS(LS_), DJ(DJ_) {
      SourceDJ = 0;
    }
    
    Jet < 1 > ::Jet(LocalState *LS_,
		    Material::Jet < 1 > *DJ_,
		    Material::Jet < 1 > *SourceDJ_) :
      LS(LS_), DJ(DJ_), SourceDJ(SourceDJ_) {
    }
    
    Jet < 1 > ::Jet(const Jet < 1 > &rhs) : LS(rhs.LS), DJ(rhs.DJ), SourceDJ(rhs.SourceDJ) {
    }

    Jet < 1 > &
    Jet < 1 > ::operator =(const Jet < 1 > &rhs) {
      if (this == &rhs) return *this;
      LS = rhs.LS;
      DJ = rhs.DJ;
      SourceDJ = rhs.SourceDJ;
      return *this;
    }

    pair<map<Set::Manifold::Point *, Set::VectorSpace::Vector>,
	 map<pair<Set::Manifold::Point *, Set::Manifold::Point *>,
	     Set::VectorSpace::Hom> >
    Jet < 1 > ::operator () (
			     const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const {
      assert(false);
    }

  }

}
