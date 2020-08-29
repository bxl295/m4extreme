// MLS.h: interface for the MLS class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_MLS_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_MLS_H__INCLUDED_

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
    namespace MLS
    {
    
      class Base;
      template<unsigned int p> class Shape;
    
      //////////////////////////////////////////////////////////////////////
      // Class Base
      //////////////////////////////////////////////////////////////////////

      class Base{
	


      public: 
	
	typedef Set::Manifold::Point          dof_type;
	typedef Set::VectorSpace::Vector      vector_type;
	typedef map<dof_type *, vector_type>  nodeset_type;       
        
      Base():x0(NULL) {}
      Base(int dim1, int m1, double beta, const nodeset_type & xloc)
	:dim(dim1), x0(&xloc), m0(m1), beta0(beta) {
	  if (dim == 1) k0 = m0+1;
	  if (dim == 2) k0 = (m0+1)*(m0+2)/2;
	  if (dim == 3) k0 = (m0+1)*(m0+2)*(m0+3)/6;
	}

	virtual ~Base(){}

	unsigned int GetDim() const {
	  return dim;
	}

	
	const map<dof_type *, vector_type> * GetNodes() const {
	  return x0;
	}

	int GetOrder() const {
	  return m0;
	}

	double GetRadius() const {
	  return beta0;
	}
		
	int GetNumberOfTerms() const {
	  return k0;
	}

	void SetNodes(const map<Set::Manifold::Point *, Set::VectorSpace::Vector> * xloc) {
	  x0 = xloc;
	}
    
      protected:
	
	unsigned int dim;
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> * x0;
	int m0;
	int k0;
	double beta0;
		
      private: 
	
	Base(const Base &);
	Base & operator = (const Base &);
	
      };

      // Class Shape<0>

      template<> class Shape<0> : public Element::Interpolation::Shape<0>
	{
	  friend class Shape<1>;

	public: 

	  typedef Element::Interpolation::MLS::Shape<1> tangent_type;

	  Shape();
	  Shape(const Base *);
	  virtual ~Shape();
	  map<Set::Manifold::Point *, double> 
	    operator () (const Set::VectorSpace::Vector &) const;
				
	protected:

	  const Base *B;

	private:

	  Shape(const Shape<0> &);
	  Shape<0> & operator = (const Shape<0> &);
	};

      // Class Shape<1>

      template<> class Shape<1> : public Element::Interpolation::Shape<1>
	{
	public: 

	  typedef Element::Interpolation::MLS::Shape<2> tangent_type;

	  Shape();
	  Shape(const Base *);
	  virtual ~Shape();
	  map<Set::Manifold::Point *, Set::VectorSpace::Vector> 
	    operator () (const Set::VectorSpace::Vector &) const;


	protected:

	  const Base *B;

	private:

	  Shape(const Shape<1> &);
	  Shape<1> & operator = (const Shape<1> &);
	};

      // Class Jet<0>

      template<unsigned int p> class Jet;

      template<> class Jet<0> 
	{
	
	public:
	
	  typedef map<Set::Manifold::Point *, double> shape_type;
	  typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> shape_derivative_type;
		
	  typedef Set::VectorSpace::Vector domain_type;
	  typedef pair<shape_type, shape_derivative_type> range_type;
	
	  Jet();
	  Jet(const Base *);
	  virtual ~Jet();
	  void operator () ( const domain_type &, shape_type &, shape_derivative_type &) const;
	  void operator () ( const domain_type &, shape_type &, shape_derivative_type &, 
			     map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
			     map<Set::Manifold::Point *, Set::VectorSpace::Hom> &,
			     map<Set::Manifold::Point *, Set::VectorSpace::Hom> &);

	private:

	  const Base * B;

	private:

	  Jet(const Jet<0> &);
	  Jet<0> & operator = (const Jet<0> &);
      
	};

    }

  }

}

#endif // !defined(ELEMENT_INTERPOLATION_MLS_H__INCLUDED_
