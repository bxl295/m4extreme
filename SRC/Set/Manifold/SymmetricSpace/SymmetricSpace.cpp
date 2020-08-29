// SymmetricSpace.cpp: Implementation of the SymmetricSpace class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "SymmetricSpace.h"

namespace Set
{
namespace SymmetricSpace
{
//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

Point::Point() {}

Point::Point(const unsigned int &n) : 
	Set::VectorSpace::Sym(n) {}

Point::Point(const unsigned int &n, double * const &u) : 
	Set::VectorSpace::Sym(n,u) {}

Set::Manifold::Point *Point::Clone() const
{
	return new Point(*this);
}

Point::~Point(){}

Point::Point(const Point &A) : 
	Set::VectorSpace::Sym(A) {}

Point::Point(const Set::VectorSpace::Sym &A) :
	Set::VectorSpace::Sym(A) {}

Point::Point(const Set::VectorSpace::Hom &A) :
	Set::VectorSpace::Sym(A) {}

Point & 
Point::operator = (const Point &A)
{
	if (this == &A) return *this; 
	Set::VectorSpace::Sym::operator = (A);
	return *this;
}

Point & 
Point::operator = (const Set::VectorSpace::Sym &A)
{
	if (this == &A) return *this; 
	Set::VectorSpace::Sym::operator = (A);
	return *this;
}

Point & 
Point::operator = (const Set::VectorSpace::Hom &A)
{
	Set::VectorSpace::Sym::operator = (A);
	return *this;
}

Set::Manifold::Point & 
Point::operator = (const Set::Manifold::Point &A)
{
	if (this == &A) return *this; 
	operator = ((Point &)A);
	return *this;
}

void 
Point::operator += (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Sym::operator += (A);
}

void 
Point::operator -= (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Sym::operator -= (A);
}

bool 
Point::operator != (const Set::Manifold::Point &P) const
{
	Point &Q = (Point &)P;
	return ((Set::VectorSpace::Sym &)(*this) != (Set::VectorSpace::Sym &)Q); 
}

bool 
Point::operator == (const Set::Manifold::Point &P) const
{
	Point &Q = (Point &)P;
	return ((Set::VectorSpace::Sym &)(*this) == (Set::VectorSpace::Sym &)Q); 
}

void 
Point::Randomize()
{
	Set::VectorSpace::Hom F(n1); F.Randomize();
	F *= 0.1; F += Set::VectorSpace::HomId(n1);
	Set::VectorSpace::Sym C(Adjoint(F)*F);
	*this = C; 
}

unsigned int 
Point::size() const
{
	return n;
}

void 
Point::print(ostream *os)
{
	Set::VectorSpace::Sym::print(os);
}

//////////////////////////////////////////////////////////////////////
// Class Vector
//////////////////////////////////////////////////////////////////////

Vector::Vector() {}

Vector::Vector(const unsigned int &n) : 
	Set::VectorSpace::Sym(n) {}

Vector::Vector(const unsigned int &n, double * const &u) : 
	Set::VectorSpace::Sym(n,u) {}

Set::Manifold::Point *Vector::Clone() const
{
	return new Vector(*this);
}

Vector::~Vector() {}

Vector::Vector(const Vector &A) : 
	Set::VectorSpace::Sym(A) {}

Vector::Vector(const Set::VectorSpace::Sym &A) : 
	Set::VectorSpace::Sym(A) {}

Vector & 
Vector::operator = (const Vector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Sym::operator = (A);
	return *this;
}

Vector & 
Vector::operator = (const Set::VectorSpace::Sym &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Sym::operator = (A);
	return *this;
}

Set::Manifold::Point & 
Vector::operator = (const Set::Manifold::Point &A)
{
	if (this == &A) return *this; 
	operator = ((Vector &)A);
	return *this;
}

void 
Vector::operator += (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Sym::operator += (A);
}

void 
Vector::operator -= (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Sym::operator -= (A);
}

void 
Vector::operator += (const Set::VectorSpace::Sym &A)
{
	Set::VectorSpace::Sym::operator += (A);
}

void 
Vector::operator -= (const Set::VectorSpace::Sym &A)
{
	Set::VectorSpace::Sym::operator -= (A);
}

bool 
Vector::operator != (const Set::Manifold::Point &P) const
{
	return (*this != (Vector &)P); 
}

bool 
Vector::operator == (const Set::Manifold::Point &P) const
{
	return (*this == (Vector &)P); 
}

void 
Vector::Randomize()
{
	Set::VectorSpace::Sym::Randomize();
}

unsigned int 
Vector::size() const
{
	return n;
}

void 
Vector::print(ostream *os)
{
	Set::VectorSpace::Sym::print(os);
}

//////////////////////////////////////////////////////////////////////
// Class VectorZero
//////////////////////////////////////////////////////////////////////

VectorZero::VectorZero() : 
	Set::VectorSpace::SymZero() {}

VectorZero::VectorZero(const unsigned int &n) : 
	Set::VectorSpace::SymZero(n) {}

VectorZero::VectorZero(const unsigned int &n, double * const &u) : 
	Set::VectorSpace::SymZero(n) {}

VectorZero::~VectorZero(){}

VectorZero::VectorZero(const VectorZero &A) : 
	Set::VectorSpace::SymZero(A) {}

VectorZero::VectorZero(const Set::VectorSpace::SymZero &A) : 
	Set::VectorSpace::SymZero(A) {}

VectorZero::VectorZero(const CovectorZero &A) : 
	Set::VectorSpace::SymZero(A) {}

VectorZero & 
VectorZero::operator = (const VectorZero &A)
{
	assert(n == A.n); 
	if (this == &A) return *this;
	return *this;
}

VectorZero & 
VectorZero::operator = (const Set::VectorSpace::SymZero &A)
{
	assert(n == A.size()); 
	return *this;
}

//////////////////////////////////////////////////////////////////////
// Class Covector
//////////////////////////////////////////////////////////////////////

Covector::Covector() {}

Covector::Covector(const unsigned int &n) : 
	Set::VectorSpace::Sym(n) {}

Covector::Covector(const unsigned int &n, double * const &u) : 
	Set::VectorSpace::Sym(n,u) {}

Covector::~Covector() {}

Covector::Covector(const Covector &A) : 
	Set::VectorSpace::Sym(A) {}

Covector::Covector(const Set::VectorSpace::Sym &A) : 
	Set::VectorSpace::Sym(A) {}

Covector::Covector(const Set::SymmetricSpace::Vector &A) : 
	Set::VectorSpace::Sym(A) {}

Covector & Covector::operator = (const Set::VectorSpace::Sym &A)
{
	Set::VectorSpace::Sym::operator = (A);
	return *this;
}

Covector & Covector::operator = (const Covector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Sym::operator = (A);
	return *this;
}

//////////////////////////////////////////////////////////////////////
// Class CovectorZero
//////////////////////////////////////////////////////////////////////

CovectorZero::CovectorZero() : 
	Set::VectorSpace::SymZero() {}

CovectorZero::CovectorZero(const unsigned int &n) : 
	Set::VectorSpace::SymZero(n) {}

CovectorZero::CovectorZero(const unsigned int &n, double * const &u) : 
	Set::VectorSpace::SymZero(n) {}

CovectorZero::~CovectorZero(){}

CovectorZero::CovectorZero(const Set::VectorSpace::SymZero &A) : 
	Set::VectorSpace::SymZero(A) {}

CovectorZero::CovectorZero(const CovectorZero &A) : 
	Set::VectorSpace::SymZero(A) {}

CovectorZero::CovectorZero(const Set::SymmetricSpace::VectorZero &A) : 
	Set::VectorSpace::SymZero(A) {}

CovectorZero & 
CovectorZero::operator = (const CovectorZero &A)
{
	assert(n == A.n); 
	if (this == &A) return *this;
	return *this;
}

CovectorZero & 
CovectorZero::operator = (const Set::VectorSpace::SymZero &A)
{
	assert(n == A.size()); 
	return *this;
}

//////////////////////////////////////////////////////////////////////
// Class ExpMap
//////////////////////////////////////////////////////////////////////

ExpMap<0>::ExpMap() {}

ExpMap<0>::~ExpMap() {}

Set::SymmetricSpace::Point
ExpMap<0>::operator () (const Set::SymmetricSpace::Vector &E)
{
	unsigned int n = E.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(E,Lambda,V);
	return this->operator () (Lambda,V);
}

Set::SymmetricSpace::Point
ExpMap<0>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	return V*Exp(Lambda)*Adjoint(V);
}

Set::SymmetricSpace::Point
ExpMap<0>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Vector E(A);
	return this->operator () (E);
}

//////////////////////////////////////////////////////////////////////
// Class ExpMap<1>
//////////////////////////////////////////////////////////////////////

ExpMap<1>::ExpMap() {}

ExpMap<1>::~ExpMap() {}

ExpMap<1>::ExpMap(const ExpMap<0> &) {} 

Set::VectorSpace::Hom
ExpMap<1>::operator () (const Set::SymmetricSpace::Vector &E)
{
	unsigned int d = E.size1();
	Set::VectorSpace::Diagonal Lambda(d);
	Set::VectorSpace::Hom V(d); 
	LinearAlgebra::EigenSym ES; ES(E,Lambda,V);
	return this->operator () (Lambda,V);
}

Set::VectorSpace::Hom
ExpMap<1>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	unsigned int d = V.size1();
	Set::VectorSpace::Diagonal LamExp = Exp(Lambda);

	unsigned int i,j,k,l,m,n;

	Set::Table f(d,d);
	for (k=0; k<d; k++) f[k][k] = LamExp[k];

	if (d > 1) 
	{
		for (j=0; j<d-1; j++) 
		{
			for (i=1; i<d; i++) 
			{
				if (Lambda[i] != Lambda[j]) 
				{
					f[j][i] = (LamExp[j] - LamExp[i])/(Lambda[j] - Lambda[i]);
					f[i][j] = f[j][i];
				}
				else 
				{
					f[j][i] = LamExp[i];
					f[i][j] = f[j][i];
				}
			}
		}
	}

	Set::VectorSpace::Hom DU(d*d); 
	double ****C = Indexing::New(DU.begin(),d,d,d,d);
	Set::VectorSpace::Hom W = Adjoint(V);

	for (i=0; i<d; i++) 
	{
		for (j=0; j<d; j++) 
		{
			for (k=0; k<d; k++) 
			{
				for (l=0; l<d; l++) 
				{
					for (m=0; m<d; m++) 
					{
						for (n=0; n<d; n++) 
						{
							C[l][k][j][i] += 
								f[n][m]*V[m][i]*W[j][n]*W[k][m]*V[n][l];
						}
					}
				}
			}
		}
	}

	Indexing::Delete(C,d,d,d,d);

	Set::VectorSpace::SymSub<1> dg(d);
	return dg(DU);
}

Set::VectorSpace::Hom
ExpMap<1>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Vector E(A);
	return this->operator () (E);
}

//////////////////////////////////////////////////////////////////////
// Class ExpMap<2>
//////////////////////////////////////////////////////////////////////

ExpMap<2>::ExpMap() {}

ExpMap<2>::~ExpMap() {}

ExpMap<2>::ExpMap(const ExpMap<1> &) {} 

Set::VectorSpace::Hom
ExpMap<2>::operator () (const Set::SymmetricSpace::Vector &E)
{
	unsigned int d = E.size1();
	Set::VectorSpace::Diagonal Lambda(d);
	Set::VectorSpace::Hom V(d); 
	LinearAlgebra::EigenSym ES; ES(E,Lambda,V);
	return this->operator () (Lambda,V);
}

Set::VectorSpace::Hom
ExpMap<2>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	unsigned int d = V.size1();
	Set::VectorSpace::Diagonal LamExp = Exp(Lambda);
//
	unsigned int i,j,k,l,m,n,ia,ib,ic;

	Set::Array F(d*d*d);
	double ***g = Indexing::New(F.begin(),d,d,d);

	for (i=0; i<d; i++) 
	{
		for (j=0; j<d; j++) 
		{
			for (k=0; k<d; k++) 
			{
				if ((Lambda[i] != Lambda[j]) && 
					(Lambda[i] != Lambda[k]) && 
					(Lambda[j] != Lambda[k])) 
				{
						g[k][j][i] = 
							(Lambda[j]*LamExp[i] - Lambda[k]*LamExp[i] - 
							 Lambda[i]*LamExp[j] + Lambda[k]*LamExp[j] + 
							 Lambda[i]*LamExp[k] - Lambda[j]*LamExp[k])/
							((Lambda[i] - Lambda[j])*(Lambda[i] - Lambda[k])
							*(Lambda[j] - Lambda[k]));
				}

				else if ((Lambda[i] == Lambda[j]) && 
					     (Lambda[i] != Lambda[k]) && 
					     (Lambda[j] != Lambda[k])) 
				{
						g[k][j][i] = 
							(-LamExp[i] + Lambda[i]*LamExp[i] 
							- Lambda[k]*LamExp[i] + LamExp[k])/
							((Lambda[i] - Lambda[k])*(Lambda[i] - Lambda[k]));
				}

				else if ((Lambda[i] != Lambda[j]) && 
					     (Lambda[i] == Lambda[k]) && 
					     (Lambda[j] != Lambda[k])) 
				{
						g[k][j][i] = 
							(-LamExp[i] + Lambda[i]*LamExp[i] 
							- Lambda[j]*LamExp[i] + LamExp[j])/
							((Lambda[i] - Lambda[j])*(Lambda[i] - Lambda[j]));
				}	

				else if ((Lambda[i] != Lambda[j]) && 
					     (Lambda[i] != Lambda[k]) && 
					     (Lambda[j] == Lambda[k])) 
				{
						g[k][j][i] = 
							(-LamExp[j] + Lambda[j]*LamExp[j] 
							- Lambda[i]*LamExp[j] + LamExp[i])/
							((Lambda[i] - Lambda[j])*(Lambda[i] - Lambda[j]));
				}

				else if ((Lambda[i] == Lambda[j]) && 
					     (Lambda[i] == Lambda[k]) &&
					     (Lambda[j] == Lambda[k])) 
				{
						g[k][j][i] = LamExp[i]/2.0;
				}
			}
		}
	}

	Set::VectorSpace::Hom DDU(d*d*d*d,d*d); 
	double ******C = Indexing::New(DDU.begin(),d,d,d,d,d,d);
	Set::VectorSpace::Hom W = Adjoint(V);

	for (i=0; i<d; i++) 
	{
		for (j=0; j<d; j++) 
		{
			for (k=0; k<d; k++) 
			{
				for (l=0; l<d; l++) 
				{
					for (m=0; m<d; m++) 
					{
						for (n=0; n<d; n++) 
						{
							for (ia=0; ia<d; ia++) 
							{
								for (ib=0; ib<d; ib++) 
								{
									for (ic=0; ic<d; ic++) 
									{
										C[n][m][l][k][j][i] += 
										  g[ic][ib][ia]*V[ia][m]*W[n][ib]*
										 (W[l][ia]*V[ib][j]*W[i][ic]*V[ic][k]
										+ W[i][ia]*V[ib][k]*W[l][ic]*V[ic][j]);
									}
								}
							}
						}
					}
				}
			}
		}
	}
//
	Indexing::Delete(C,d,d,d,d,d,d);
	Indexing::Delete(g,d,d,d);

	Set::VectorSpace::SymSub<2> ddg(d);
	return ddg(DDU);
}

Set::VectorSpace::Hom
ExpMap<2>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Vector E(A);
	return this->operator () (E);
}

//////////////////////////////////////////////////////////////////////
// Class ExpJet
//////////////////////////////////////////////////////////////////////

ExpJet<0>::ExpJet() {}

ExpJet<0>::~ExpJet() {}

pair<Set::SymmetricSpace::Point,Set::VectorSpace::Hom>
ExpJet<0>::operator () (const Set::SymmetricSpace::Vector &E)
{
	unsigned int n = E.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(E,Lambda,V);
	return this->operator () (Lambda,V);
}

pair<Set::SymmetricSpace::Point,Set::VectorSpace::Hom>
ExpJet<0>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	Set::SymmetricSpace::ExpMap<0> Exp; 
	Set::SymmetricSpace::ExpMap<1> DExp; 
	return make_pair(Exp(Lambda,V),DExp(Lambda,V));
}

pair<Set::SymmetricSpace::Point,Set::VectorSpace::Hom>
ExpJet<0>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Vector E(A);
	return this->operator () (E);
}

//////////////////////////////////////////////////////////////////////
// Class ExpJet
//////////////////////////////////////////////////////////////////////

ExpJet<1>::ExpJet() {}

ExpJet<1>::~ExpJet() {}

pair<Set::VectorSpace::Hom,Set::VectorSpace::Hom>
ExpJet<1>::operator () (const Set::SymmetricSpace::Vector &E)
{
	unsigned int n = E.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(E,Lambda,V);
	return this->operator () (Lambda,V);
}

pair<Set::VectorSpace::Hom,Set::VectorSpace::Hom>
ExpJet<1>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	Set::SymmetricSpace::ExpMap<1> DExp; 
	Set::SymmetricSpace::ExpMap<2> DDExp; 
	return make_pair(DExp(Lambda,V),DDExp(Lambda,V));
}

pair<Set::VectorSpace::Hom,Set::VectorSpace::Hom>
ExpJet<1>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Vector E(A);
	return this->operator () (E);
}

//////////////////////////////////////////////////////////////////////
// Class ExpJetJet
//////////////////////////////////////////////////////////////////////

ExpJetJet<0>::ExpJetJet() {}

ExpJetJet<0>::~ExpJetJet() {}

triplet<Set::SymmetricSpace::Point,
	 Set::VectorSpace::Hom,
	 Set::VectorSpace::Hom>
ExpJetJet<0>::operator () (const Set::SymmetricSpace::Vector &E)
{
	unsigned int n = E.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(E,Lambda,V);
	return this->operator () (Lambda,V);
}

triplet<Set::SymmetricSpace::Point,
	 Set::VectorSpace::Hom,
	 Set::VectorSpace::Hom>
ExpJetJet<0>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	Set::SymmetricSpace::ExpMap<0> Exp; 
	Set::SymmetricSpace::ExpMap<1> DExp; 
	Set::SymmetricSpace::ExpMap<2> DDExp; 
	return make_triplet(Exp(Lambda,V),DExp(Lambda,V),DDExp(Lambda,V));
}

triplet<Set::SymmetricSpace::Point,
	 Set::VectorSpace::Hom,
	 Set::VectorSpace::Hom>
ExpJetJet<0>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Vector E(A);
	return this->operator () (E);
}

//////////////////////////////////////////////////////////////////////
// Class LogMap<0>
//////////////////////////////////////////////////////////////////////

LogMap<0>::LogMap() {}

LogMap<0>::~LogMap() {}

Set::VectorSpace::Sym
LogMap<0>::operator () (const Set::SymmetricSpace::Point &U)
{
	unsigned int n = U.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(U,Lambda,V);
	return this->operator () (Lambda,V);
}

Set::VectorSpace::Sym
LogMap<0>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
  if ( Lambda.size() == 3 ) {
    double * Lh = Lambda.begin();
    double Lg0 = log(*Lh);
    double Lg1 = log(*(Lh+1));
    double Lg2 = log(*(Lh+2));

    double * Vh = V.begin();
    double v00 = *Vh,     v10 = *(Vh+1), v20 = *(Vh+2);
    double v01 = *(Vh+3), v11 = *(Vh+4), v21 = *(Vh+5);
    double v02 = *(Vh+6), v12 = *(Vh+7), v22 = *(Vh+8);

    double S[9];
    S[0] = Lg0 * v00 * v00 + Lg1 * v01 * v01 + Lg2 * v02 * v02;
    S[1] = S[3] = Lg0 * v00 * v10 + Lg1 * v11 * v01 + Lg2 * v12 * v02;
    S[2] = S[6] = Lg0 * v00 * v20 + Lg1 * v21 * v01 + Lg2 * v22 * v02;
    S[4] = Lg0 * v10 * v10 + Lg1 * v11 * v11 + Lg2 * v12 * v12;
    S[5] = S[7] = Lg0 * v10 * v20 + Lg1 * v11 * v21 + Lg2 * v12 * v22;
    S[8] = Lg0 * v20 * v20 + Lg1 * v21 * v21 + Lg2 * v22 * v22;
        
    return Set::VectorSpace::Hom(3, S); 
  }
  else{
    return V*Log(Lambda)*Adjoint(V);
  }
}

Set::VectorSpace::Sym
LogMap<0>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Point U(A);
	return this->operator () (U);
}

//////////////////////////////////////////////////////////////////////
// Class LogMap<1>
//////////////////////////////////////////////////////////////////////

LogMap<1>::LogMap() {}

LogMap<1>::~LogMap() {}

LogMap<1>::LogMap(const LogMap<0> &) {} 

Set::VectorSpace::Hom
LogMap<1>::operator () (const Set::SymmetricSpace::Point &U)
{
	unsigned int d = U.size1();
	Set::VectorSpace::Diagonal Lambda(d);
	Set::VectorSpace::Hom V(d); 
	LinearAlgebra::EigenSym ES; ES(U,Lambda,V);
	return this->operator () (Lambda,V);
}

 Set::VectorSpace::Hom
 LogMap<1>::operator () (
 	const Set::VectorSpace::Diagonal &Lambda,
 	const Set::VectorSpace::Hom &V)
{
  if ( V.size1() == 3 ) {
    return operator3D(Lambda, V);
  }
  else {
    return operatorGeneral(Lambda, V);
  }
}

Set::VectorSpace::Hom
LogMap<1>::operatorGeneral (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	unsigned int d = V.size1();
	Set::VectorSpace::Diagonal LamLog = Log(Lambda);

	unsigned int i,j,k,l,m,n;

	Set::Table f(d,d);
	for (k=0; k<d; k++) f[k][k] = 1.0/Lambda[k];

	if (d > 1) {
		for (j=0; j<d-1; j++) {
			for (i=1; i<d; i++) {
				if (Lambda[i] != Lambda[j]) {
					f[j][i] = (LamLog[j] - LamLog[i])/(Lambda[j] - Lambda[i]);
					f[i][j] = f[j][i];
				}
				else {
					f[j][i] = 1.0/Lambda[i];
					f[i][j] = f[j][i];
				}
			}
		}
	}

	Set::VectorSpace::Hom DE(d*d); 

	double ****C = Indexing::New(DE.begin(),d,d,d,d);

	Set::VectorSpace::Hom W = Adjoint(V);
	
	for (i=0; i<d; i++) {
		for (j=0; j<d; j++) {
			for (k=0; k<d; k++) {
				for (l=0; l<d; l++) {
					for (m=0; m<d; m++) {
						for (n=0; n<d; n++) {
							C[l][k][j][i] += 
							  f[n][m]*V[m][i]*W[j][n]*W[k][m]*V[n][l];
						}
					}
				}
			}
		}
	}

	Indexing::Delete(C,d,d,d,d);	

	Set::VectorSpace::SymSub<1> dg(d);
	return dg(DE);
}


Set::VectorSpace::Hom
LogMap<1>::operator3D (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	Set::VectorSpace::Diagonal LamLog = Log(Lambda);

	unsigned int i,j,k,l,m,n;

	Set::Table f(3,3);
	for (k=0; k<3; k++) f[k][k] = 1.0/Lambda[k];


	for (j=0; j<2; j++) {
	  for (i=1; i<3; i++) {
	    if (Lambda[i] != Lambda[j]) {
	      f[j][i] = (LamLog[j] - LamLog[i])/(Lambda[j] - Lambda[i]);
	      f[i][j] = f[j][i];
	    }
	    else {
	      f[j][i] = 1.0/Lambda[i];
	      f[i][j] = f[j][i];
	    }
	  }
	}
	
	double f00 = f[0][0], f01 = f[0][1], f02 = f[0][2];
	double f10 = f[1][0], f11 = f[1][1], f12 = f[1][2];
	double f20 = f[2][0], f21 = f[2][1], f22 = f[2][2];

	Set::VectorSpace::Hom W = Adjoint(V);
	
	double **Vd = Indexing::New(V.begin(), 3, 3);
	double **Wd = Indexing::New(W.begin(), 3, 3);

	double * Vd0 = Vd[0];
	double * Vd1 = Vd[1];
	double * Vd2 = Vd[2];

	double *Wdj, *Wdk;
	double Wdj0, Wdj1, Wdj2;
	double Wdk0, Wdk1, Wdk2;
	double Vd0i, Vd1i, Vd2i;
	double Ctemp = 0.0;

	double CArray[81];

	for (i=0; i<3; i++) {
	  Vd0i = Vd0[i];
	  Vd1i = Vd1[i];
	  Vd2i = Vd2[i];
	  for (j=0; j<3; j++) {
	    Wdj = Wd[j];
	    Wdj0 = Wdj[0];
	    Wdj1 = Wdj[1];
	    Wdj2 = Wdj[2];
	    for (k=0; k<3; k++) {
	      Wdk = Wd[k];
	      Wdk0 = Wdk[0];
	      Wdk1 = Wdk[1];
	      Wdk2 = Wdk[2];
	      for (l=0; l<3; l++) {
		Ctemp = 0.0;
		
		//m=0
		////n=0
		Ctemp += f00*Vd0i*Wdj0*Wdk0*Vd0[l];
		////n=1
		Ctemp += f10*Vd0i*Wdj1*Wdk0*Vd1[l];	
		////n=2
		Ctemp += f20*Vd0i*Wdj2*Wdk0*Vd2[l];		
		
		//m=1
		////n=0
		Ctemp += f01*Vd1i*Wdj0*Wdk1*Vd0[l];
		////n=1
		Ctemp += f11*Vd1i*Wdj1*Wdk1*Vd1[l];	
		////n=2
		Ctemp += f21*Vd1i*Wdj2*Wdk1*Vd2[l];		
		
		//m=2
		////n=0
		Ctemp += f02*Vd2i*Wdj0*Wdk2*Vd0[l];
		////n=1
		Ctemp += f12*Vd2i*Wdj1*Wdk2*Vd1[l];	
		////n=2
		Ctemp += f22*Vd2i*Wdj2*Wdk2*Vd2[l];		
		
		*(CArray+l*27+k*9+j*3+i) = Ctemp;  
	      }
	    }
	  }
	}

	Indexing::Delete(Wd, 3, 3);
	Indexing::Delete(Vd, 3, 3);

	Set::VectorSpace::SymSub<1> dg(3);
 	Set::VectorSpace::Hom B(6);
 	dg.operator3D(CArray, B);	
 	return B;
}


Set::VectorSpace::Hom
LogMap<1>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Point U(A);
	return this->operator () (U);
}

//////////////////////////////////////////////////////////////////////
// Class LogMap<2>
//////////////////////////////////////////////////////////////////////

LogMap<2>::LogMap() {}

LogMap<2>::~LogMap() {}

LogMap<2>::LogMap(const LogMap<1> &) {} 

Set::VectorSpace::Hom
LogMap<2>::operator () (const Set::SymmetricSpace::Point &U)
{
	unsigned int d = U.size1();
	Set::VectorSpace::Diagonal Lambda(d);
	Set::VectorSpace::Hom V(d); 
	LinearAlgebra::EigenSym ES; ES(U,Lambda,V);
	return this->operator () (Lambda,V);
}

 Set::VectorSpace::Hom
 LogMap<2>::operator () (
 	const Set::VectorSpace::Diagonal &Lambda,
 	const Set::VectorSpace::Hom &V)
{
  if ( V.size1() == 3 ) {
    return operator3D(Lambda, V);
  }
  else {
    return operatorGeneral(Lambda, V);
  }
}

Set::VectorSpace::Hom
LogMap<2>::operatorGeneral (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	unsigned int d = V.size1();
	unsigned int n1 = d*d;
//
	unsigned int i,j,k,l,m,n,ia,ib,ic;

	Set::VectorSpace::Diagonal LamLog = Log(Lambda);
	Set::VectorSpace::Hom W = Adjoint(V);

	Set::Table f(d,d);
	for (k=0; k<d; k++) f[k][k] = 1.0/Lambda[k];

	if (d > 1) {
		for (j=0; j<d-1; j++) {
			for (i=1; i<d; i++) {
				if (Lambda[i] != Lambda[j]) {
					f[j][i] = (LamLog[j] - LamLog[i])/(Lambda[j] - Lambda[i]);
					f[i][j] = f[j][i];
				}
				else {
					f[j][i] = 1.0/Lambda[i];
					f[i][j] = f[j][i];
				}
			}
		}
	}

	Set::VectorSpace::Hom DL(d*d); 
	double ****B = Indexing::New(DL.begin(),d,d,d,d);

	for (i=0; i<d; i++) {
		for (j=0; j<d; j++) {
			for (k=0; k<d; k++) {
				for (l=0; l<d; l++) {
					for (m=0; m<d; m++) {
						for (n=0; n<d; n++) {
							B[l][k][j][i] += f[n][m]*V[m][i]*W[j][n]*W[k][m]*V[n][l];
						}
					}
				}
			}
		}
	}

	Indexing::Delete(B,d,d,d,d);
//
	Set::Array F(d*d*d);
	double ***g = Indexing::New(F.begin(),d,d,d);
	double **dlogmap = Indexing::New(DL.begin(),n1,n1);

	for (i=0; i<d; i++) 
	{
		for (j=0; j<d; j++) 
		{
			for (k=0; k<d; k++) 
			{
				if ((LamLog[i] != LamLog[j]) && 
					(LamLog[i] != LamLog[k]) && 
					(LamLog[j] != LamLog[k])) 
				{
						g[k][j][i] = 
							(LamLog[j]*Lambda[i] - LamLog[k]*Lambda[i] - 
							 LamLog[i]*Lambda[j] + LamLog[k]*Lambda[j] + 
							 LamLog[i]*Lambda[k] - LamLog[j]*Lambda[k])/
							((LamLog[i] - LamLog[j])*(LamLog[i] - LamLog[k])
							*(LamLog[j] - LamLog[k]));
				}

				else if ((LamLog[i] == LamLog[j]) && 
					     (LamLog[i] != LamLog[k]) && 
					     (LamLog[j] != LamLog[k])) 
				{
						g[k][j][i] = 
							(-Lambda[i] + LamLog[i]*Lambda[i] 
							- LamLog[k]*Lambda[i] + Lambda[k])/
							((LamLog[i] - LamLog[k])*(LamLog[i] - LamLog[k]));
				}

				else if ((LamLog[i] != LamLog[j]) && 
					     (LamLog[i] == LamLog[k]) && 
					     (LamLog[j] != LamLog[k])) 
				{
						g[k][j][i] = 
							(-Lambda[i] + LamLog[i]*Lambda[i] 
							- LamLog[j]*Lambda[i] + Lambda[j])/
							((LamLog[i] - LamLog[j])*(LamLog[i] - LamLog[j]));
				}	

				else if ((LamLog[i] != LamLog[j]) && 
					     (LamLog[i] != LamLog[k]) && 
					     (LamLog[j] == LamLog[k])) 
				{
						g[k][j][i] = 
							(-Lambda[j] + LamLog[j]*Lambda[j] 
							- LamLog[i]*Lambda[j] + Lambda[i])/
							((LamLog[i] - LamLog[j])*(LamLog[i] - LamLog[j]));
				}

				else if ((LamLog[i] == LamLog[j]) && 
					     (LamLog[i] == LamLog[k]) &&
					     (LamLog[j] == LamLog[k])) 
				{
						g[k][j][i] = Lambda[i]/2.0;
				}
			}
		}
	}

	Set::VectorSpace::Hom DDE(d*d*d*d,d*d); 
	double ******C = Indexing::New(DDE.begin(),d,d,d,d,d,d);

	for (i=0; i<d; i++) 
	{
		for (j=0; j<d; j++) 
		{
			for (k=0; k<d; k++) 
			{
				for (l=0; l<d; l++) 
				{
					for (m=0; m<d; m++) 
					{
						for (n=0; n<d; n++) 
						{
							for (ia=0; ia<d; ia++) 
							{
								for (ib=0; ib<d; ib++) 
								{
									for (ic=0; ic<d; ic++) 
									{
									  C[n][m][l][k][j][i] +=
 										  g[ic][ib][ia]*V[ia][m]*W[n][ib]*
 										 (W[l][ia]*V[ib][j]*W[i][ic]*V[ic][k]
 										+ W[i][ia]*V[ib][k]*W[l][ic]*V[ic][j]);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	Indexing::Delete(C,d,d,d,d,d,d);
	Indexing::Delete(g,d,d,d);
//
	double ***ddexpmap = Indexing::New(DDE.begin(),n1,n1,n1);

	Set::VectorSpace::Hom C1(n1*n1,n1);
	double ***dum1 = Indexing::New(C1.begin(),n1,n1,n1);

	for (i=0; i<n1; i++) 
	{
		for (j=0; j<n1; j++) 
		{
			for (k=0; k<n1; k++)
			{
				for (l=0; l<n1; l++)
				{
					dum1[k][j][i] -= ddexpmap[l][j][i]*dlogmap[k][l];
				}
			}
		}
	}

	Set::VectorSpace::Hom C2(n1*n1,n1);
	double ***dum2 = Indexing::New(C2.begin(),n1,n1,n1);

	for (i=0; i<n1; i++) 
	{
		for (j=0; j<n1; j++) 
		{
			for (k=0; k<n1; k++)
			{
				for (l=0; l<n1; l++)
				{
					dum2[k][j][i] += dum1[k][l][i]*dlogmap[j][l];
				}
			}
		}
	}

	Set::VectorSpace::Hom C3(n1*n1,n1);
	double ***dum3 = Indexing::New(C3.begin(),n1,n1,n1);

	for (i=0; i<n1; i++) 
	{
		for (j=0; j<n1; j++) 
		{
			for (k=0; k<n1; k++)
			{
				for (l=0; l<n1; l++)
				{
					dum3[k][j][i] += dum2[k][j][l]*dlogmap[l][i];
				}
			}
		}
	}

	Indexing::Delete(dum3,n1,n1,n1);
	Indexing::Delete(dum2,n1,n1,n1);
	Indexing::Delete(dum1,n1,n1,n1);
	Indexing::Delete(ddexpmap,n1,n1,n1);
	Indexing::Delete(dlogmap,n1,n1);
//
	Set::VectorSpace::SymSub<2> ddg(d);
	return ddg(C3);
}

Set::VectorSpace::Hom
LogMap<2>::operator3D(
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	unsigned int i,j,k,l,m,n,ia,ib,ic;

	Set::VectorSpace::Diagonal LamLog = Log(Lambda);
	Set::VectorSpace::Hom W = Adjoint(V);

	Set::Table f(3,3);
	for (k=0; k<3; k++) f[k][k] = 1.0/Lambda[k];

	for (j=0; j<2; j++) {
	  for (i=1; i<3; i++) {
	    if (Lambda[i] != Lambda[j]) {
	      f[j][i] = (LamLog[j] - LamLog[i])/(Lambda[j] - Lambda[i]);
	      f[i][j] = f[j][i];
	    }
	    else {
	      f[j][i] = 1.0/Lambda[i];
	      f[i][j] = f[j][i];
	    }
	  }
	}
	

	Set::VectorSpace::Hom DL(9); 
	double ****B = Indexing::New(DL.begin(),3,3,3,3);

	double **Vd = Indexing::New(V.begin(), 3, 3);
	double **Wd = Indexing::New(W.begin(), 3, 3);

	double xtemp = 0.0;
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			for (k=0; k<3; k++) {
				for (l=0; l<3; l++) {
				        xtemp = 0.0;
					for (m=0; m<3; m++) {
						for (n=0; n<3; n++) {
							xtemp += 
							  f[n][m]*Vd[m][i]*Wd[j][n]*Wd[k][m]*Vd[n][l];
						}
					}
					B[l][k][j][i] = xtemp;
				}
			}
		}
	}

	Indexing::Delete(B,3,3,3,3);
//
	Set::Array F(27);
	double ***g = Indexing::New(F.begin(),3,3,3);
	double **dlogmap = Indexing::New(DL.begin(),9,9);

	for (i=0; i<3; i++) 
	{
		for (j=0; j<3; j++) 
		{
			for (k=0; k<3; k++) 
			{
				if ((LamLog[i] != LamLog[j]) && 
					(LamLog[i] != LamLog[k]) && 
					(LamLog[j] != LamLog[k])) 
				{
						g[k][j][i] = 
							(LamLog[j]*Lambda[i] - LamLog[k]*Lambda[i] - 
							 LamLog[i]*Lambda[j] + LamLog[k]*Lambda[j] + 
							 LamLog[i]*Lambda[k] - LamLog[j]*Lambda[k])/
							((LamLog[i] - LamLog[j])*(LamLog[i] - LamLog[k])
							*(LamLog[j] - LamLog[k]));
				}

				else if ((LamLog[i] == LamLog[j]) && 
					     (LamLog[i] != LamLog[k]) && 
					     (LamLog[j] != LamLog[k])) 
				{
						g[k][j][i] = 
							(-Lambda[i] + LamLog[i]*Lambda[i] 
							- LamLog[k]*Lambda[i] + Lambda[k])/
							((LamLog[i] - LamLog[k])*(LamLog[i] - LamLog[k]));
				}

				else if ((LamLog[i] != LamLog[j]) && 
					     (LamLog[i] == LamLog[k]) && 
					     (LamLog[j] != LamLog[k])) 
				{
						g[k][j][i] = 
							(-Lambda[i] + LamLog[i]*Lambda[i] 
							- LamLog[j]*Lambda[i] + Lambda[j])/
							((LamLog[i] - LamLog[j])*(LamLog[i] - LamLog[j]));
				}	

				else if ((LamLog[i] != LamLog[j]) && 
					     (LamLog[i] != LamLog[k]) && 
					     (LamLog[j] == LamLog[k])) 
				{
						g[k][j][i] = 
							(-Lambda[j] + LamLog[j]*Lambda[j] 
							- LamLog[i]*Lambda[j] + Lambda[i])/
							((LamLog[i] - LamLog[j])*(LamLog[i] - LamLog[j]));
				}

				else if ((LamLog[i] == LamLog[j]) && 
					     (LamLog[i] == LamLog[k]) &&
					     (LamLog[j] == LamLog[k])) 
				{
						g[k][j][i] = Lambda[i]/2.0;
				}
			}
		}
	}

	Set::VectorSpace::Hom DDE(81,9); 
	double ******C = Indexing::New(DDE.begin(),3,3,3,3,3,3);

	double *Wdn, *Wdl, *Wdi;
	double *Vd0 = Vd[0];
	double *Vd1 = Vd[1];
	double *Vd2 = Vd[2];

	double Wdi0, Wdi1, Wdi2, Wdl0, Wdl1, Wdl2;
	double Vd0j, Vd1j, Vd2j;
	
	for (i=0; i<3; i++) 
	  {
	    Wdi = Wd[i];
	    Wdi0 = Wdi[0];
	    Wdi1 = Wdi[1];
	    Wdi2 = Wdi[2];
	    for (j=0; j<3; j++) 
	      {
		Vd0j = Vd0[j];
		Vd1j = Vd1[j];
		Vd2j = Vd2[j];
		for (k=0; k<3; k++) 
		  {
		    for (l=0; l<3; l++) 
		      {
			Wdl = Wd[l];
			Wdl0 = Wdl[0];
			Wdl1 = Wdl[1];
			Wdl2 = Wdl[2];
			for (m=0; m<3; m++) 
			  {
			    for (n=0; n<3; n++) 
			      {
				Wdn = Wd[n];
				xtemp = 0.0;
				//ia=0
				   ////ib=0
				      //////ic=0
				      xtemp += g[0][0][0]*Vd0[m]*Wdn[0]
					* (Wdl0*Vd0j*Wdi0*Vd0[k]
					   + Wdi0*Vd0[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][0][0]*Vd0[m]*Wdn[0]
					* (Wdl0*Vd0j*Wdi1*Vd1[k]
					   + Wdi0*Vd0[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][0][0]*Vd0[m]*Wdn[0]
					* (Wdl0*Vd0j*Wdi2*Vd2[k]
					   + Wdi0*Vd0[k]*Wdl2*Vd2j);

				  ////ib=1
				      //////ic=0
				      xtemp += g[0][1][0]*Vd0[m]*Wdn[1]
					* (Wdl0*Vd1j*Wdi0*Vd0[k]
					   + Wdi0*Vd1[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][1][0]*Vd0[m]*Wdn[1]
					* (Wdl0*Vd1j*Wdi1*Vd1[k]
					   + Wdi0*Vd1[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][1][0]*Vd0[m]*Wdn[1]
					* (Wdl0*Vd1j*Wdi2*Vd2[k]
					   + Wdi0*Vd1[k]*Wdl2*Vd2j);

				  ////ib=2
				      //////ic=0
				      xtemp += g[0][2][0]*Vd0[m]*Wdn[2]
					* (Wdl0*Vd2j*Wdi0*Vd0[k]
					   + Wdi0*Vd2[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][2][0]*Vd0[m]*Wdn[2]
					* (Wdl0*Vd2j*Wdi1*Vd1[k]
					   + Wdi0*Vd2[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][2][0]*Vd0[m]*Wdn[2]
					* (Wdl0*Vd2j*Wdi2*Vd2[k]
					   + Wdi0*Vd2[k]*Wdl2*Vd2j);
				    
			      //ia=1
				  ////ib=0
				      //////ic=0
				      xtemp += g[0][0][1]*Vd1[m]*Wdn[0]
					* (Wdl1*Vd0j*Wdi0*Vd0[k]
					   + Wdi1*Vd0[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][0][1]*Vd1[m]*Wdn[0]
					* (Wdl1*Vd0j*Wdi1*Vd1[k]
					   + Wdi1*Vd0[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][0][1]*Vd1[m]*Wdn[0]
					* (Wdl1*Vd0j*Wdi2*Vd2[k]
					   + Wdi1*Vd0[k]*Wdl2*Vd2j);

				  ////ib=1
				      //////ic=0
				      xtemp += g[0][1][1]*Vd1[m]*Wdn[1]
					* (Wdl1*Vd1j*Wdi0*Vd0[k]
					   + Wdi1*Vd1[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][1][1]*Vd1[m]*Wdn[1]
					* (Wdl1*Vd1j*Wdi1*Vd1[k]
					   + Wdi1*Vd1[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][1][1]*Vd1[m]*Wdn[1]
					* (Wdl1*Vd1j*Wdi2*Vd2[k]
					   + Wdi1*Vd1[k]*Wdl2*Vd2j);

				  ////ib=2
				      //////ic=0
				      xtemp += g[0][2][1]*Vd1[m]*Wdn[2]
					* (Wdl1*Vd2j*Wdi0*Vd0[k]
					   + Wdi1*Vd2[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][2][1]*Vd1[m]*Wdn[2]
					* (Wdl1*Vd2j*Wdi1*Vd1[k]
					   + Wdi1*Vd2[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][2][1]*Vd1[m]*Wdn[2]
					* (Wdl1*Vd2j*Wdi2*Vd2[k]
					   + Wdi1*Vd2[k]*Wdl2*Vd2j);

				
			      //ia=2
				  ////ib=0
				      //////ic=0
				      xtemp += g[0][0][2]*Vd2[m]*Wdn[0]
					* (Wdl2*Vd0j*Wdi0*Vd0[k]
					   + Wdi2*Vd0[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][0][2]*Vd2[m]*Wdn[0]
					* (Wdl2*Vd0j*Wdi1*Vd1[k]
					   + Wdi2*Vd0[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][0][2]*Vd2[m]*Wdn[0]
					* (Wdl2*Vd0j*Wdi2*Vd2[k]
					   + Wdi2*Vd0[k]*Wdl2*Vd2j);

				  ////ib=1
				      //////ic=0
				      xtemp += g[0][1][2]*Vd2[m]*Wdn[1]
					* (Wdl2*Vd1j*Wdi0*Vd0[k]
					   + Wdi2*Vd1[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][1][2]*Vd2[m]*Wdn[1]
					* (Wdl2*Vd1j*Wdi1*Vd1[k]
					   + Wdi2*Vd1[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][1][2]*Vd2[m]*Wdn[1]
					* (Wdl2*Vd1j*Wdi2*Vd2[k]
					   + Wdi2*Vd1[k]*Wdl2*Vd2j);

				  ////ib=2
				      //////ic=0
				      xtemp += g[0][2][2]*Vd2[m]*Wdn[2]
					* (Wdl2*Vd2j*Wdi0*Vd0[k]
					   + Wdi2*Vd2[k]*Wdl0*Vd0j);
				      //////ic=1
				      xtemp += g[1][2][2]*Vd2[m]*Wdn[2]
					* (Wdl2*Vd2j*Wdi1*Vd1[k]
					   + Wdi2*Vd2[k]*Wdl1*Vd1j);
				      //////ic=2
				      xtemp += g[2][2][2]*Vd2[m]*Wdn[2]
					* (Wdl2*Vd2j*Wdi2*Vd2[k]
					   + Wdi2*Vd2[k]*Wdl2*Vd2j);
			    
				 C[n][m][l][k][j][i] = xtemp;
			      }
			  }
		      }
		  }
	      }
	  }

	Indexing::Delete(Wd, 3, 3);
	Indexing::Delete(Vd, 3, 3);

	Indexing::Delete(C,3,3,3,3,3,3);
	Indexing::Delete(g,3,3,3);
//
	double ***ddexpmap = Indexing::New(DDE.begin(),9,9,9);

	Set::VectorSpace::Hom C1(81,9);
	double ***dum1 = Indexing::New(C1.begin(),9,9,9);

	for (i=0; i<9; i++) 
	{
		for (j=0; j<9; j++) 
		{
			for (k=0; k<9; k++)
			{
				for (l=0; l<9; l++)
				{
					dum1[k][j][i] -= ddexpmap[l][j][i]*dlogmap[k][l];
				}
			}
		}
	}

	Set::VectorSpace::Hom C2(81,9);
	double ***dum2 = Indexing::New(C2.begin(),9,9,9);

	for (i=0; i<9; i++) 
	{
		for (j=0; j<9; j++) 
		{
			for (k=0; k<9; k++)
			{
				for (l=0; l<9; l++)
				{
					dum2[k][j][i] += dum1[k][l][i]*dlogmap[j][l];
				}
			}
		}
	}

	Set::VectorSpace::Hom C3(81,9);
	double ***dum3 = Indexing::New(C3.begin(),9,9,9);

	for (i=0; i<9; i++) 
	{
		for (j=0; j<9; j++) 
		{
			for (k=0; k<9; k++)
			{
				for (l=0; l<9; l++)
				{
					dum3[k][j][i] += dum2[k][j][l]*dlogmap[l][i];
				}
			}
		}
	}

	Indexing::Delete(dum3,9,9,9);
	Indexing::Delete(dum2,9,9,9);
	Indexing::Delete(dum1,9,9,9);
	Indexing::Delete(ddexpmap,9,9,9);
	Indexing::Delete(dlogmap,9,9);
//
	Set::VectorSpace::SymSub<2> ddg(3);
	return ddg(C3);
}

Set::VectorSpace::Hom
LogMap<2>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Point U(A);
	return this->operator () (U);
}

//////////////////////////////////////////////////////////////////////
// Class LogJet
//////////////////////////////////////////////////////////////////////

LogJet<0>::LogJet() {}

LogJet<0>::~LogJet() {}

pair<Set::VectorSpace::Sym,Set::VectorSpace::Hom>
LogJet<0>::operator () (const Set::SymmetricSpace::Point &U)
{
	unsigned int n = U.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(U,Lambda,V);
	return this->operator () (Lambda,V);
}

pair<Set::VectorSpace::Sym,Set::VectorSpace::Hom>
LogJet<0>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	Set::SymmetricSpace::LogMap<0> Log; 
	Set::SymmetricSpace::LogMap<1> DLog; 
	return make_pair(Log(Lambda,V),DLog(Lambda,V));
}

pair<Set::VectorSpace::Sym,Set::VectorSpace::Hom>
LogJet<0>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Point U(A);
	return this->operator () (U);
}

//////////////////////////////////////////////////////////////////////
// Class LogJet
//////////////////////////////////////////////////////////////////////

LogJet<1>::LogJet() {}

LogJet<1>::~LogJet() {}

pair<Set::VectorSpace::Hom,Set::VectorSpace::Hom>
LogJet<1>::operator () (const Set::SymmetricSpace::Point &U)
{
	unsigned int n = U.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(U,Lambda,V);
	return this->operator () (Lambda,V);
}

pair<Set::VectorSpace::Hom,Set::VectorSpace::Hom>
LogJet<1>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	Set::SymmetricSpace::LogMap<1> DLog; 
	Set::SymmetricSpace::LogMap<2> DDLog; 
	return make_pair(DLog(Lambda,V),DDLog(Lambda,V));
}

pair<Set::VectorSpace::Hom,Set::VectorSpace::Hom>
LogJet<1>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Point U(A);
	return this->operator () (U);
}

//////////////////////////////////////////////////////////////////////
// Class LogJetJet
//////////////////////////////////////////////////////////////////////

LogJetJet<0>::LogJetJet() {}

LogJetJet<0>::~LogJetJet() {}

triplet<Set::VectorSpace::Sym,
	 Set::VectorSpace::Hom,
	 Set::VectorSpace::Hom>
LogJetJet<0>::operator () (const Set::SymmetricSpace::Point &U)
{
	unsigned int n = U.size1();
	Set::VectorSpace::Diagonal Lambda(n);
	Set::VectorSpace::Hom V(n);
	LinearAlgebra::EigenSym ES; ES(U,Lambda,V);
	return this->operator () (Lambda,V);
}

triplet<Set::VectorSpace::Sym,
	 Set::VectorSpace::Hom,
	 Set::VectorSpace::Hom>
LogJetJet<0>::operator () (
	const Set::VectorSpace::Diagonal &Lambda,
	const Set::VectorSpace::Hom &V)
{
	Set::SymmetricSpace::LogMap<0> Log; 
	Set::SymmetricSpace::LogMap<1> DLog; 
	Set::SymmetricSpace::LogMap<2> DDLog; 
	return make_triplet(Log(Lambda,V),DLog(Lambda,V),DDLog(Lambda,V));
}

triplet<Set::VectorSpace::Sym,
	 Set::VectorSpace::Hom,
	 Set::VectorSpace::Hom>
LogJetJet<0>::operator () (const Set::VectorSpace::Vector &A)
{
	const Set::SymmetricSpace::Point U(A);
	return this->operator () (U);
}

}

}

Set::SymmetricSpace::Point
operator + (const Set::SymmetricSpace::Point &P,
			const Set::SymmetricSpace::Vector &A)
{
	Set::SymmetricSpace::Point Q = P; 
	Q += A; return Q;
}

Set::SymmetricSpace::Point
operator - (const Set::SymmetricSpace::Point &P,
			const Set::SymmetricSpace::Vector &A)
{
	Set::SymmetricSpace::Point Q = P; 
	Q -= A; return Q;
}

Set::SymmetricSpace::Vector 
operator - (const Set::SymmetricSpace::Point &Q, 
			const Set::SymmetricSpace::Point &P)
{
	Set::SymmetricSpace::Vector A(P); 
	Set::SymmetricSpace::Vector B(Q); 
	return B-A;
}

Set::SymmetricSpace::Point
operator + (const Set::SymmetricSpace::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::SymmetricSpace::Point Q = P; 
	Q += A; return Q;
}

Set::SymmetricSpace::Point
operator - (const Set::SymmetricSpace::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::SymmetricSpace::Point Q = P; 
	Q -= A; return Q;
}

void Random(Set::SymmetricSpace::Point &P){P.Randomize();}
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &A,
			const Set::VectorSpace::Vector &B)
{
	Set::SymmetricSpace::Vector C = A;
	C += B; return C;
}

Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &A,
			const Set::VectorSpace::Vector &B)
{
	Set::SymmetricSpace::Vector C = A;
	C -= B; return C;
}
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &A,
			const Set::VectorSpace::Sym &B)
{
	Set::SymmetricSpace::Vector C = A;
	C += B; return C;
}

Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &A,
			const Set::VectorSpace::Sym &B)
{
	Set::SymmetricSpace::Vector C = A;
	C -= B; return C;
}
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &A,
			const Set::VectorSpace::SymZero &B)
{
	return A;
}

Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &A,
			const Set::VectorSpace::SymZero &B)
{
	return A;
}
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &A,
			const Set::SymmetricSpace::Vector &B)
{
	Set::SymmetricSpace::Vector C = A;
	C += B; return C;
}

Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &A,
			const Set::SymmetricSpace::Vector &B)
{
	Set::SymmetricSpace::Vector C = A;
	C -= B; return C;
}

