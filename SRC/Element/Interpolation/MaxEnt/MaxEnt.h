// MaxEnt.h: interface for the MaxEnt class.
// Copyright (c) 2006 by Bo Li - All rights reserved
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_MAXENT_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_MAXENT_H__INCLUDED_

#define _M4EXTREME_MODIFIED_NEWTON_ 
#define _M4EXTREME_MAXENT_SYMMETRIC_CASE_ 
#define _M4EXTREME_VISCOUS_REGULARIZATION_

#pragma once

#include <math.h>
#include <map>
#include <set>
#include "Element/Interpolation/Interpolation.h"

using namespace std;

namespace Element
{
  namespace Interpolation
  {
    namespace MaxEnt
    {
      class Base;
      template<unsigned int p> class Shape;

      //////////////////////////////////////////////////////////////////////
      // Class Base
      //////////////////////////////////////////////////////////////////////

      class Base
      {
	friend class Shape<0>;
	friend class Shape<1>;

      public: 

        typedef Set::Manifold::Point          dof_type;
        typedef Set::VectorSpace::Vector      vector_type;
	typedef map<dof_type *, vector_type>  nodeset_type;       
        
	Base();
	Base(const double &, const map<dof_type *, vector_type> &);
	virtual ~Base();
        
	void operator() (const vector_type &);
        
        void  SetNodes(const nodeset_type *);
        const nodeset_type * GetNodes() const;
        
        void  SetBeta(const double &);
        const double & GetBeta() const;

        const Set::VectorSpace::Vector & GetLambda() const;
        const vector<double> & GetN() const;
        const Set::VectorSpace::Hom * GetJInverse() const ;

	void write(ostream & os) const {
	  os.write((char*)&_Beta, sizeof(double));
	  _Lambda->write(os);
	}

	void read(istream & is) {
	  is.read((char*)&_Beta, sizeof(double));

	  _Lambda = new Set::VectorSpace::Vector();
	  _Lambda->read(is);

	  _JInv = new Set::VectorSpace::Hom(_Lambda->size());
	}
        
      private: 
        Base(const Base &);
	Base & operator = (const Base &);
	
	void _computeExactLambda(const Set::VectorSpace::Vector &);
	void _computeLambda1D(const Set::VectorSpace::Vector &);
        void _computeLambda2D(const Set::VectorSpace::Vector &);
        void _computeLambda3D(const Set::VectorSpace::Vector &);
 
#if defined (_M4EXTREME_VISCOUS_REGULARIZATION_)
	void _get_psum(double** p, double * psum, int ndim) const;

	// 2D
	void _computeLambda2D_Regularized(const Set::VectorSpace::Vector &qp);
	bool _logZ(double & value,
		   double * Lambda,
		   double beta,
		   const vector<double> & dx,
		   const vector<double> & dy,
		   const vector<double> & dxdx) const;
	double _NMtry(double ** p, double * y, double * psum, int ihi, double fac,
		      double beta,
		      const vector<double> & dx,
		      const vector<double> & dy,
		      const vector<double> & dxdx) const;
	void _NelderMead(double & Z,
			 double beta,
			 double ftol,
			 const vector<double> & dx,
			 const vector<double> & dy,
			 const vector<double> & dxdx);

	// 3D
	void _computeLambda3D_Regularized(const Set::VectorSpace::Vector &qp);
	bool _logZ(double & value,
		   double * Lambda,
		   double beta,
		   const vector<double> & dx,
		   const vector<double> & dy,
		   const vector<double> & dz,
		   const vector<double> & dxdx) const;
	double _NMtry(double ** p, double * y, double * psum, int ihi, double fac,
		      double beta,
		      const vector<double> & dx,
		      const vector<double> & dy,
		      const vector<double> & dz,
		      const vector<double> & dxdx) const;
	void _NelderMead(double & Z,
			 double beta,
			 double ftol,
			 const vector<double> & dx,
			 const vector<double> & dy,
			 const vector<double> & dz,
			 const vector<double> & dxdx);
#endif       
        double _computeLineSearchf(double, 
				   const vector<double> &, 
				   const vector<double> &, 
				   const vector<double> & ) const;
        
        double _bisection(double, 
			  const vector<double> &, 
			  const vector<double> &, 
			  const vector<double> & ) const;

      private:
        double _Beta;
        const nodeset_type *       _x;
        Set::VectorSpace::Vector * _Lambda;
        vector<double> _NNodes;
        Set::VectorSpace::Hom *  _JInv;

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
	class _loop_list {
	public:
	  _loop_list():_pointer(0){
	    _value[0] = _value[1] = 0.0;
	  }
	  
	  double & head() {
	    return _value[(int)_pointer];
	  }

	  void operator++() {
	    _pointer = !_pointer;
	  }

	  void reset() {
	    _value[0] = _value[1] = 0;
	    _pointer = 0;
	  }

	private:
	  double _value[2];
	  bool _pointer;
	};
#endif
	
      };


      //////////////////////////////////////////////////////////////////////
      // Class Shape<0>
      //////////////////////////////////////////////////////////////////////

      template<> class Shape<0> : public Element::Interpolation::Shape<0>
	{
	  friend class Shape<1>;

	public: 

	  typedef Element::Interpolation::MaxEnt::Shape<1> tangent_type;

	  Shape();
	  Shape(Base *);
	  virtual ~Shape();
	  map<Set::Manifold::Point *, double> 
	    operator () (const Set::VectorSpace::Vector &) const;

	private:

	  Base *B;

	private:

	  Shape(const Shape<0> &);
	  Shape<0> & operator = (const Shape<0> &);
	};

      //////////////////////////////////////////////////////////////////////
      // Class Shape<1>
      //////////////////////////////////////////////////////////////////////

      template<> class Shape<1> : public Element::Interpolation::Shape<1>
	{
	public: 

	  typedef Element::Interpolation::MaxEnt::Shape<2> tangent_type;

	  Shape();
	  Shape(Base *);
	  Shape(const Shape<0> &);
	  virtual ~Shape();
	  map<Set::Manifold::Point *, Set::VectorSpace::Vector> 
	    operator () (const Set::VectorSpace::Vector &) const;

	private:

	  Base *B;

	private:

	  Shape(const Shape<1> &);
	  Shape<1> & operator = (const Shape<1> &);
	};

      /////////////////////////////////////////////////////////////////////
      // Class Jet<p>
      //////////////////////////////////////////////////////////////////////

      template<unsigned int p> class Jet;

      //////////////////////////////////////////////////////////////////////
      // Class Jet<0>
      //////////////////////////////////////////////////////////////////////

      template<>
	class Jet<0> {
      public:
	typedef map<Set::Manifold::Point *, double>                   shape_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> shape_derivative_type;

	typedef Set::VectorSpace::Vector                              domain_type;
	typedef pair<shape_type, shape_derivative_type>               range_type;

	Jet();
	Jet(Base *);
	virtual ~Jet();
	void operator () ( const domain_type &, shape_type &, shape_derivative_type &) const;

      private:
	Base * _BNode;

      private:
	Jet<0>(const Jet<0> &);
	Jet<0> & operator = (const Jet<0> &);
      };

    }

  }

}

#endif // !defined(ELEMENT_INTERPOLATION_MAXENT_H__INCLUDED_
