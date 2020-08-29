// Cholesky.cpp: Implementation of the Cholesky class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "Cholesky.h"

namespace LinearAlgebra  
{
  Cholesky::Cholesky() : n(0), p(0){}

  Cholesky::Cholesky(const unsigned int &n0) 
  : n(n0), p(new double [n])
  {
    for (unsigned int i=0; i<n; i++) p[i]=0.0;
  }

  Cholesky::~Cholesky()
  {
    delete [] p;
  }

  void 
  Cholesky::Decomposition(Set::VectorSpace::Hom &a)
  {
    assert (n == a.size1());
    assert (n == a.size2());

    unsigned int i, j, k; 
    double sum;

    for (i=0; i<n; i++)
      {
	for (j=i; j<n; j++)
	  {
	    sum = a[j][i];
	    for (k=0; k<i; k++) {
	      sum -= a[k][i]*a[k][j];
	    }

	    if (i == j)
	      {
		if (sum <= 0.0) {
		  throw(0);
		}
		p[i] = sqrt(sum);
	      }
	    else {
	      if ( fabs(p[i]) < 1.0e-12 ) {
		cerr << "bad value in Cholesky solver::Decomposition" << endl;
	      }

	      a[i][j] = sum/p[i];
	    }
		
	  }
      }
  }

  void 
  Cholesky::Substitution(const Set::VectorSpace::Hom &a, 
			 Set::VectorSpace::Vector &b)
  {
    assert (n == a.size1());
    assert (n == a.size2());
    assert (n == b.size());

    int i, k; int nn = (int) n;
    double sum;

    for (i=0; i<nn; i++)
      {
	sum = b[i];
	for (k=i-1; k>=0; k--) {
	  sum -= a[k][i]*b[k];
	}

	if ( fabs(p[i]) < 1.0e-12 ) {
	  cerr << "bad value in Cholesky solver::Substitution" << endl;
	}
	
	b[i] = sum/p[i];
      }

    for (i=nn-1; i>=0; i--)
      {
	sum = b[i];
	for (k=i+1; k<nn; k++) {
	  sum -= a[i][k]*b[k];
	}

	if ( fabs(p[i]) < 1.0e-12 ) {
	  cerr << "bad value in Cholesky solver::Substitution" << endl;
	}

	b[i] = sum/p[i];
      }
  }

  void 
  Cholesky::Substitution(const Set::VectorSpace::Hom &a, 
			 const Set::VectorSpace::Vector &b, 
			 Set::VectorSpace::Vector &x)
  {
    assert (n == a.size1());
    assert (n == a.size2());
    assert (n == b.size());
    assert (n == x.size());

    int i; 
    int nn = (int) n;
    for (i=0; i<nn; i++) {
      x[i] = b[i];
    }

    Substitution(a,x);
  }

  void 
  Cholesky::Solve(const Set::VectorSpace::Hom &a0, 
		  Set::VectorSpace::Vector &b)
  {
    assert (n == a0.size1());
    assert (n == a0.size2());

    Set::VectorSpace::Hom a(a0);

    try{
      Decomposition(a);
    }
    catch(int flag){
      throw(flag);
    }

    Substitution(a,b);
  }

  void 
  Cholesky::Solve(const Set::VectorSpace::Hom &a0, 
		  const Set::VectorSpace::Vector &b0, 
		  Set::VectorSpace::Vector &x)
  {
    assert (n == a0.size1());
    assert (n == a0.size2());
    assert (n == b0.size());
    assert (n == x.size());

    Set::VectorSpace::Hom a(a0);
    Set::VectorSpace::Vector b(b0);

    try{
      Decomposition(a);
    }
    catch(int flag){
      throw(flag);
    }

    Substitution(a,b,x);
  }

  void 
  Cholesky::Invert(const Set::VectorSpace::Hom &a0, 
		   Set::VectorSpace::Hom &ainv)
  {
    assert (n == a0.size1());
    assert (n == a0.size2());
    assert (n == ainv.size1());
    assert (n == ainv.size2());

    unsigned int i, j;
    Set::VectorSpace::Vector b(n);
    Set::VectorSpace::Hom a(a0);

    try{
      Decomposition(a);
    }
    catch(int flag){
      throw(flag);
    }

    for (j=0; j<n; j++)
      {
	for (i=0; i<n; i++) {
	  b[i] = 0.0;
	}
	b[j] = 1.0;
	Substitution(a,b,ainv[j]);
      }
  }

  double 
  Cholesky::Determinant(const Set::VectorSpace::Hom &a0)
  {
    assert (n == a0.size1());
    assert (n == a0.size2());

    unsigned int i; 
    double det;
    Set::VectorSpace::Hom a(a0);

    try{
      Decomposition(a);
    }
    catch(int flag){
      if (flag==0) {
	return 0.0;
      }
    }

    det = 1.0;
    for (i=0; i<n; i++) {
      det *= p[i]*p[i];
    }

    return det;
  }

}
