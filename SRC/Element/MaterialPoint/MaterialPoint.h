// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_MATERIALPOINT_H__INCLUDED_)
#define ELEMENT_MATERIALPOINT_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include "Element/Element.h"
#include "Element/Interpolation/MaxEnt/MaxEnt.h"
#include "Element/Interpolation/MLS/MLS.h"
#include "Material/Material.h"
#include "Geometry/Search/SearchLib.h"

#define _M4EXTREME_MIXING_SEPERATION_

using namespace std;

namespace Element
{
  namespace MaterialPoint
  {
    class Data;
    class LocalState;
    template <unsigned int> class Energy;
    template <unsigned int> class Jet;

    //////////////////////////////////////////////////////////////////////
    // Class Data
    //////////////////////////////////////////////////////////////////////

    class Data
    {
      friend class LumpedMass;
      friend class LocalState;
      friend class Energy<0>;
      friend class Energy<1>;
      friend class Energy<2>;
      friend class Jet<0>;
      friend class Jet<1>;

    public: 

      typedef Set::Manifold::Point          dof_type;
      typedef Set::VectorSpace::Vector      vector_type;
      typedef map<dof_type *, vector_type>  domain_type; 
      typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> nodeset_type;       
      typedef nodeset_type dshape_type;
      typedef map<Set::Manifold::Point *, double> shape_type;

      Data();
      virtual ~Data();

      Data( const double &,
	    const double &,
	    const double &,
	    const vector<double> &,	      
	    const vector<vector_type> &,
	    const domain_type &,
	    bool = false,
	    Geometry::Search<dof_type*> ** = NULL );
        
      Data(const Data &);
      Data & operator = (const Data &);

      const vector<shape_type> & GetN() const { return N; }
      vector<shape_type> & GetN() { return N; }
      const vector<dshape_type> & GetDN() const { return DN; }
      vector<dshape_type> & GetDN() { return DN; }
      const vector<double> & GetQW() const { return QW; }
      vector<double> & GetQW() { return QW; }
      const double & GetJ(int q) const { return _Jold[q]; }
      double & GetJ(int q) { return _Jold[q]; }
      const vector<vector_type> & GetQP() const { return QP; }
      vector<vector_type> & GetQP() { return QP; }
      const double & GetVolume() const { return _volume; }
      double & GetVolume() { return _volume; }
      const vector_type & GetCentroid() const { return _centroid; }
      vector_type & GetCentroid() { return _centroid; }
      double GetMaxRange() const { return _maxRange; }
      double & GetMaxRange() { return _maxRange; }   
      double GetCutoff() const { return _cutoff; }

      void SetAdaptiveSearch(bool, double);        
      void SetCarrier ( set<dof_type*> *);
      set<dof_type*> * GetCarrier() const { return _carrier; }
      void SetBlackList ( set<dof_type*> *);
      set<dof_type*> * GetBlackList() const { return _blacklist; }
      set<dof_type*> & GetContactList() { return _nodesInTension; }
      void SetMaxReset ( unsigned int );
      void SetSearch ( Geometry::Search<dof_type*> **);
     
      void operator ++ ();         
      void Reset(const domain_type &);
      void Remesh(const domain_type &);
      virtual void operator() (const nodeset_type *) = 0;
      virtual void enlargeSupport() {}

    protected:
      void _check_consistency(const nodeset_type &);

    protected:
      map<dof_type *, vector_type> rxaq;

      vector<double> Mass;
      vector<double> QW;
      vector<vector_type> QP;
      vector<shape_type> N; 
      vector<dshape_type> DN;

    private:
      double _volume;
      vector_type _centroid;
      vector<double> _Jold;
      vector<double> _DJ;
      int _readflag;
      double _maxRange;
      double _cutoff;
      bool _adaptiveSearch;
      bool _updateNeighbor;
      double _searchRange;
      unsigned int _MAX_RESET;
        
      Geometry::Search<dof_type*> ** _nbs;
      set<dof_type*>   *  _carrier;        
      set<dof_type*>   *  _blacklist;        
      set<dof_type*> _nodesInTension;	
      bool _isSharedInterface;
    };

    //////////////////////////////////////////////////////////////////////
    // Class LMEData
    //////////////////////////////////////////////////////////////////////

    class LMEData : public Data 
    {
    public:
      LMEData();
      LMEData( const double &,
	       const double &, 
	       const double &,
	       const double &,
	       const vector<double> &,	      
	       const vector<vector_type>  &,
	       const domain_type &,
	       bool = true,
	       Geometry::Search<dof_type*> ** = NULL );
      virtual ~LMEData();
    	
      LMEData(const LMEData &);
      LMEData & operator = (const LMEData &);
    
      virtual void operator() (const nodeset_type *);
      virtual void enlargeSupport();

    private:
      Element::Interpolation::MaxEnt::Base * _Base;
    };

    //////////////////////////////////////////////////////////////////////
    // Class MLSData
    //////////////////////////////////////////////////////////////////////

    class MLSData : public Data 
    {
    public:
      MLSData();
      MLSData( const double &,
	       const int &,
	       const double &, 
	       const double &,
	       const double &,
	       const vector<double> &,	      
	       const vector<vector_type> &,
	       const domain_type &,
	       bool = true,
	       Geometry::Search<dof_type*> ** = NULL );
      virtual ~MLSData();
    	
      MLSData(const MLSData &);
      MLSData & operator = (const MLSData &);
    
      virtual void operator() (const nodeset_type *);
      virtual void enlargeSupport();

    private:
      Element::Interpolation::MLS::Base * _Base;
    };

    //////////////////////////////////////////////////////////////////////
    // Class MixedData
    //////////////////////////////////////////////////////////////////////

    class MixedData : public Data 
    {
    public:
      enum TYPE { FEA, MIXED, MLS };

    public:
      MixedData();
      MixedData(const double &,
		const int &,
		const double &,
		const double &,
		const double &,
		const vector<double> &,	      
		const vector<vector_type> &,
		const domain_type &,
		vector<dof_type*>,
		const set<dof_type*> *,
		const set<dof_type*> *,
		bool = true,
		Geometry::Search<dof_type*> ** = NULL );
      virtual ~MixedData();
    	
      MixedData(const MixedData &);
      MixedData & operator = (const MixedData &);
    
      virtual void operator() (const nodeset_type *);
    private:
      void monomials(const vector_type &, vector_type &) const;
      void compute(const nodeset_type &);
      void computeFE(const nodeset_type &, 
		     const vector_type &,
		     map<dof_type *, double> &,
		     map<dof_type *, vector_type> &) const;
    
    private:
      double _volume;
      unsigned int _m;
      double _rho;
      TYPE _type; 
      vector<dof_type*> _FirstNeighbor;
      const set<dof_type*> * _FEAnodes;
      const set<dof_type*> * _MLSnodes;
      Element::Interpolation::MLS::Base * _MLSBase;
    };

    //////////////////////////////////////////////////////////////////////
    // Class LocalState
    //////////////////////////////////////////////////////////////////////

    class LocalState : public Element::LocalState
      {
	friend class Energy<0>;
	friend class Energy<1>;
	friend class Energy<2>;
	friend class Jet<0>;
	friend class Jet<1>;

      public: 
	typedef Set::Manifold::Point          dof_type;
	typedef map<dof_type *, Set::VectorSpace::Vector> domain_type;
	typedef vector<Set::VectorSpace::Vector> range_type;  
	typedef domain_type dshape_type;
	typedef map<dof_type *, double> shape_type;

	using Element::LocalState::Reset;

	LocalState();
	virtual ~LocalState();
	Element::LocalState *Clone() const;
	LocalState(Data *, const vector<Material::LocalState *> &);
	LocalState(Data *, const vector<Material::LocalState *> &, const vector<Material::LocalState *> &);
	LocalState(const LocalState &);
	LocalState & operator = (const LocalState &);
	
	virtual range_type operator () (const domain_type &) const;
	virtual void operator ++ (); 
        set<dof_type *> GetNodes() const;
	set<pair<Set::Manifold::Point *, 
	  Set::Manifold::Point *> > GetNodePairs() const;
	virtual void Reset(const domain_type &);

        virtual void operator () (const domain_type &, range_type&) const;
	virtual void GetDisplacement(const domain_type &, vector<Set::VectorSpace::Vector> &) const;
        virtual void GetNodes(set<dof_type *> &) const;
	virtual void GetNodePairs(set<pair<dof_type *, dof_type *> > &) const; 

	const vector<shape_type> & GetN() const { return D->GetN(); }
	vector<shape_type> & GetN() { return D->GetN(); }
	const vector<dshape_type> & GetDN() const { return D->GetDN(); }
	vector<dshape_type> & GetDN() { return D->GetDN(); }
	const vector<double> & GetQW() const { return D->GetQW(); } 
        vector<double> & GetQW() { return D->GetQW(); }        
        const vector<vector_type> & GetQP() const { return D->GetQP(); }
        vector<vector_type> & GetQP() { return D->GetQP(); }

	const double & GetVolume() const { return D->GetVolume(); }
	const vector_type & GetCentroid() const { return D->GetCentroid(); }
	const vector_type & GetNormal() const { return Normal; }
	vector_type & GetNormal() { return Normal; }
	const double & GetFriction() const { return friction_coefficient; }
	double & GetFriction() { return friction_coefficient; }
        const Data * GetData() const{return D;}
        Data * GetData() {return D;}   
	const vector<Material::LocalState *> & GetMLS() const {return MLS;}

	double & GetHourglassModulus() { return hourglass_modulus; }
	const double & GetHourglassModulus() const { return hourglass_modulus; }
	
      protected:
	Set::VectorSpace::Hom POld;
	vector_type Normal;
	double friction_coefficient;
	Data *D;
	vector<Material::LocalState *> MLS;
	vector<Material::LocalState *> SourceLS;
	double hourglass_modulus;
      };

    //////////////////////////////////////////////////////////////////////
    // Class Energy<p>
    //////////////////////////////////////////////////////////////////////

    template<unsigned int p> class Energy;

    //////////////////////////////////////////////////////////////////////
    // Class Energy<0>
    //////////////////////////////////////////////////////////////////////

    template<>
      class Energy<0> : public Element::Energy<0>
      {
      public: 

	typedef Element::MaterialPoint::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<0> *Clone() const;
	Element::LocalState *GetLocalState() const;
	const vector<Material::Energy<0> *> & GetW() const;
	Energy(LocalState *, const vector<Material::Energy<0> *> &);
	Energy(LocalState *, const vector<Material::Energy<0> *> &, const vector<Material::Energy<0> *> &);
	Energy(const Energy<0> &);
	virtual Energy<0> & operator = (const Energy<0> &);
	virtual range_type operator () (const domain_type &) const;

      protected:

	LocalState *LS;
	vector<Material::Energy<0> *> W;
	vector<Material::Energy<0> *> SourceW;
      };

    //////////////////////////////////////////////////////////////////////
    // Class Energy<1>
    //////////////////////////////////////////////////////////////////////

    template<>
      class Energy<1> : public Element::Energy<1>
      {
      public: 

	typedef Element::MaterialPoint::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<1> *Clone() const;
	Element::LocalState *GetLocalState() const;
	const vector<Material::Energy<1> *> & GetDW() const;
	Energy(LocalState *, const vector<Material::Energy<1> *> &);
	Energy(LocalState *, const vector<Material::Energy<1> *> &, const vector<Material::Energy<1> *> &);
	Energy(const Energy<1> &);
	virtual Energy<1> & operator = (const Energy<1> &);
	virtual range_type operator () (const domain_type &) const;
        virtual void operator () (const domain_type &, range_type &) const;
        
      protected:

	LocalState *LS;
	vector<Material::Energy<1> *> DW;
	vector<Material::Energy<1> *> SourceDW;
      };

    //////////////////////////////////////////////////////////////////////
    // Class Energy<2>
    //////////////////////////////////////////////////////////////////////

    template<>
      class Energy<2> : public Element::Energy<2>
      {
      public: 

	typedef Element::MaterialPoint::Energy<3> tangent_type;
	typedef map<Set::Manifold::Point *, 
	  Set::VectorSpace::Vector> domain_type;
	typedef map<pair<Set::Manifold::Point *, 
	  Set::Manifold::Point *>, 
	  Set::VectorSpace::Hom> range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<2> *Clone() const;
	Element::LocalState *GetLocalState() const;
	const vector<Material::Energy<2> *> & GetDDW() const;
	Energy(LocalState *, const vector<Material::Energy<2> *> &);
	Energy(LocalState *, const vector<Material::Energy<2> *> &, const vector<Material::Energy<2> *> &);
	Energy(const Energy<2> &);
	virtual Energy<2> & operator = (const Energy<2> &);
	virtual range_type operator () (const domain_type &) const;

      protected:

	LocalState *LS;
	vector<Material::Energy<2> *> DDW;
	vector<Material::Energy<2> *> SourceDDW;
      };

    //////////////////////////////////////////////////////////////////////
    // Class Jet<p>
    //////////////////////////////////////////////////////////////////////

    template<unsigned int p> class Jet;

    //////////////////////////////////////////////////////////////////////
    // Class Jet<0>
    //////////////////////////////////////////////////////////////////////

    template<> 
      class Jet<0> : public Element::Jet<0>
      {
      public: 

	typedef Element::Jet<1> tangent_type;
	typedef map<Set::Manifold::Point *, 
	  Set::VectorSpace::Vector> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, 
	  Set::VectorSpace::Vector> > range_type;

	Jet();
	virtual ~Jet();
	Element::Jet<0> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Jet(LocalState *, Material::Jet<0> *);
	Jet(LocalState *, Material::Jet<0> *, Material::Jet<0> *);
	Jet(const Jet<0> &);
	virtual Jet<0> & operator = (const Jet<0> &);
	virtual range_type operator () (const domain_type &) const;

      protected:

	LocalState *LS;
	Material::Jet<0> *J;
	Material::Jet<0> *SourceJ;
      };

    //////////////////////////////////////////////////////////////////////
    // Class Jet<1>
    //////////////////////////////////////////////////////////////////////

    template <> 
      class Jet<1> : public Element::Jet<1>
      {
      public: 

	typedef Element::Jet<2> tangent_type;
	typedef map<Set::Manifold::Point *, 
	  Set::VectorSpace::Vector> domain_type;
	typedef pair<
	  map<Set::Manifold::Point *, 
	  Set::VectorSpace::Vector>, 
	  map<pair<Set::Manifold::Point *, 
	  Set::Manifold::Point *>, 
	  Set::VectorSpace::Hom> > range_type;

	Jet();
	virtual ~Jet();
	Element::Jet<1> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Jet(LocalState *, Material::Jet<1> *);
	Jet(LocalState *, Material::Jet<1> *, Material::Jet<1> *);
	Jet(const Jet<1> &);
	virtual Jet<1> & operator = (const Jet<1> &);
	virtual range_type operator () (const domain_type &) const;

      protected:

	LocalState *LS;
	Material::Jet<1> *DJ;
	Material::Jet<1> *SourceDJ;
      };

  }

}

#endif // !defined(ELEMENT_MATERIALPOINT_H__INCLUDED_
