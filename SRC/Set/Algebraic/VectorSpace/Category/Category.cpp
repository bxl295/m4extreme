// Category.cpp: implementation of the Category class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "Category.h" 
#include "../../../../Utils/LinearAlgebra/Crout/Crout.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "jama_eig.h"

using namespace std;

namespace Set
{
namespace VectorSpace
{
void Copy(const Set::VectorSpace::Hom &A,
		  TNT::Array2D<double> &B)
{
	unsigned int i1, n1 = A.size1(); int j1;
	unsigned int i2, n2 = A.size2(); int j2;
	for (i2=0, j2=0; i2<n2; i2++, j2++)
		for (i1=0, j1=0; i1<n1; i1++, j1++)
			B[j1][j2] = A[i2][i1];
}

//////////////////////////////////////////////////////////////////////
// Class Hom
//////////////////////////////////////////////////////////////////////

void 
Hom::print(ostream *os) const
{
	for (unsigned int i=0; i<n2; i++) *os << *L[i];
}


void Hom::operator += (const HomZero &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
}


void Hom::operator -= (const HomZero &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
}


double Hom::operator () (const HomZero &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
	return 0.0;
}


void Hom::operator += (const HomId &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
	unsigned int m1 = n1+1; 
	double *p;
	for (p=Vector::head; p<Vector::tail; p+=m1) *p += 1.0;
}


void Hom::operator -= (const HomId &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
	unsigned int m1 = n1+1; 
	double *p;
	for (p=Vector::head; p<Vector::tail; p+=m1) *p -= 1.0;
}


double Hom::operator () (const HomId &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
	unsigned int m1 = n1+1; 
	double a = 0.0; double *p;
	for (p=Vector::head; p<Vector::tail; p+=m1) a += *p;
	return a;
}


void Hom::write(ostream & os) const 
{
  os.write((char*)&n1, sizeof(int));
  os.write((char*)&n2, sizeof(int));

  Array::write(os);
}


void Hom::read(istream & is)
{
  is.read((char*)&n1, sizeof(int));
  is.read((char*)&n2, sizeof(int));  
  Array::read(is);

  L = new Vector * [n2];
  double *p = this->head;
  for (size_t i=0; i != n2; i++, p+=n1) L[i] = new Vector(n1,p);
}
	

//////////////////////////////////////////////////////////////////////
// Class HomZero
//////////////////////////////////////////////////////////////////////

HomZero::HomZero() : Hom(){}

HomZero::HomZero(const unsigned int &m1, 
				 const unsigned int &m2) : Hom(m1,m2){}

HomZero::HomZero(const unsigned int &m1, 
				 const unsigned int &m2, 
				 double * const &u) : Hom(m1,m2,u){}

HomZero::HomZero(const unsigned int &m) : Hom(m){}

HomZero::HomZero(const unsigned int &m, 
				 double * const &u) : Hom(m,u){}

HomZero::~HomZero(){}

HomZero::HomZero(const Hom &A) : Hom(A){}

HomZero::HomZero(const HomZero &A) : Hom(A){}

HomZero & 
HomZero::operator = (const HomZero &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
	return *this;
}

HomZero & 
HomZero::operator = (const VectorZero &A)
{
	assert(Vector::n == A.size()); 
	return *this;
}

void 
HomZero::operator += (const HomZero &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
}

void 
HomZero::operator -= (const HomZero &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
}

void 
HomZero::operator *= (const double &){}

void 
HomZero::operator /= (const double &){}

VectorZero
HomZero::operator () (const Vector &A) const
{
	assert (A.size() == n2); 
	Set::VectorSpace::VectorZero B(n1); return B;
}

double 
HomZero::operator () (const Vector &A, const Vector &B) const
{
	assert (A.size() == n2); 
	assert (B.size() == n1); 
	return 0.0;
}

double 
HomZero::operator () (const Hom &A)
{
	assert(n1 == A.size1()); 
	assert(n2 == A.size2()); 
	return 0.0;
}

double 
HomZero::operator () (const HomId &A)
{
	assert(n1 == A.size1()); 
	assert(n2 == A.size2()); 
	return 0.0;
}

//////////////////////////////////////////////////////////////////////
// Class HomId
//////////////////////////////////////////////////////////////////////

HomId::HomId() : Hom(){}

HomId::HomId(const unsigned int &m) : Hom(m)
{
	unsigned int m1=n1+1; double *p;
	for (p=this->head; p<this->tail; p+=m1) *p = 1.0;
}

HomId::HomId(const unsigned int &m, double * const &u) : Hom(m,u)
{
	unsigned int m1=n1+1; double *p;
	for (p=this->head; p<this->tail; p+=m1) *p = 1.0;
}

HomId::~HomId(){}

HomId::HomId(const HomId &A) : Hom(A)
{
	unsigned int m1=n1+1; double *p;
	for (p=this->head; p<this->tail; p+=m1) *p = 1.0;
}

HomId & 
HomId::operator = (const HomId &A)
{
	assert(n1 == A.n1); 
	assert(n2 == A.n2); 
	return *this;
}

void 
HomId::operator += (const HomZero &A)
{
	assert(n1 == A.size1()); 
	assert(n2 == A.size2()); 
}

void 
HomId::operator -= (const HomZero &A)
{
	assert(n1 == A.size1()); 
	assert(n2 == A.size2()); 
}

Vector 
HomId::operator () (const Vector &A) const
{
	assert (A.size() == n2); 
	return A;
}

double 
HomId::operator () (const Vector &A, const Vector &B) const
{
	return B(A);
}

double 
HomId::operator () (const Hom &A)
{
	assert(n1 == A.size1()); 
	assert(n2 == A.size2()); 
	return Trace(A);
}

double 
HomId::operator () (const HomZero &A)
{
	assert(n1 == A.size1()); 
	assert(n2 == A.size2()); 
	return 0.0;
}

}

}

//////////////////////////////////////////////////////////////////////
// Class Hom
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Hom &A, 
	    const Set::VectorSpace::Hom &B)
{
        unsigned int n = B.size1();
	assert(n == A.size2());
	unsigned int n1=A.size1(), n2=B.size2();
	Set::VectorSpace::Hom C(n1,n2);

	if ( n1 == n2 && n1 == 3 && n == 3 ) {
	  double * Ah = A.begin();
	  double * Bh = B.begin();
	  double * Ch = C.begin();
  
	  *Ch     = *Ah     * *Bh     + *(Ah+3) * *(Bh+1) + *(Ah+6) * *(Bh+2);
	  *(Ch+1) = *(Ah+1) * *Bh     + *(Ah+4) * *(Bh+1) + *(Ah+7) * *(Bh+2);   
	  *(Ch+2) = *(Ah+2) * *Bh     + *(Ah+5) * *(Bh+1) + *(Ah+8) * *(Bh+2);
	  *(Ch+3) = *Ah     * *(Bh+3) + *(Ah+3) * *(Bh+4) + *(Ah+6) * *(Bh+5);
	  *(Ch+4) = *(Ah+1) * *(Bh+3) + *(Ah+4) * *(Bh+4) + *(Ah+7) * *(Bh+5);   
	  *(Ch+5) = *(Ah+2) * *(Bh+3) + *(Ah+5) * *(Bh+4) + *(Ah+8) * *(Bh+5);
	  *(Ch+6) = *Ah     * *(Bh+6) + *(Ah+3) * *(Bh+7) + *(Ah+6) * *(Bh+8);
	  *(Ch+7) = *(Ah+1) * *(Bh+6) + *(Ah+4) * *(Bh+7) + *(Ah+7) * *(Bh+8);   
	  *(Ch+8) = *(Ah+2) * *(Bh+6) + *(Ah+5) * *(Bh+7) + *(Ah+8) * *(Bh+8);
	}  
	else {
	  unsigned int i1, i2, k;
	  for (i2=0; i2<n2; i2++) 
	    for (i1=0; i1<n1; i1++) 
	      for (k=0; k<n; k++) 
		C[i2][i1] += B[i2][k]*A[k][i1];
	}

	return C;
}

void Identity(Set::VectorSpace::Hom &A)
{
	assert (A.size1() == A.size2());
	unsigned int i, j, n=A.size1(); 
	for (j=0; j<n; j++) 
		for (i=0; i<n; i++)
			A[j][i] = 0.0; 
	for (i=0; i<n; i++) A[i][i] = 1.0;
}
/*
Set::VectorSpace::Hom
Inverse(const Set::VectorSpace::Hom &A)
{
	assert (A.size1() == A.size2());
	Set::VectorSpace::Hom B(A.size1());
	LinearAlgebra::Crout LA(A.size1());
	LA.Invert(A,B); return B;
}
*/
Set::VectorSpace::Hom
Inverse(const Set::VectorSpace::Hom &A)
{
	unsigned int n=A.size1();
	assert (n == A.size2());
	Set::VectorSpace::Hom B(n);
	if (n == 3)
	{
		double Det = 
		A[0][0]*A[1][1]*A[2][2] + A[1][0]*A[2][1]*A[0][2] + 
		A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0] - 
		A[0][1]*A[1][0]*A[2][2] - A[1][2]*A[2][1]*A[0][0];
		if (Det == 0.0) throw(0);
		B[0][0] = (A[1][1]*A[2][2] - A[2][1]*A[1][2])/Det;
		B[1][0] = (A[2][0]*A[1][2] - A[1][0]*A[2][2])/Det;
		B[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0])/Det;
		B[0][1] = (A[2][1]*A[0][2] - A[0][1]*A[2][2])/Det;
		B[1][1] = (A[0][0]*A[2][2] - A[2][0]*A[0][2])/Det;
		B[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1])/Det;
		B[0][2] = (A[0][1]*A[1][2] - A[1][1]*A[0][2])/Det;
		B[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2])/Det;
		B[2][2] = (A[0][0]*A[1][1] - A[1][0]*A[0][1])/Det;
		return B;
	}
	if (n == 2)
	{
		double Det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
		if (Det == 0.0) throw(0);
		B[0][0] =  A[1][1]/Det;
		B[0][1] = -A[0][1]/Det;
		B[1][0] = -A[1][0]/Det;
		B[1][1] =  A[0][0]/Det;
		return B;
	}
	if (n == 1)
	{
		if (A[0][0] == 0.0) throw(0);
		B[0][0] =  1.0/A[0][0];
		return B;
	}
	if (n == 0) return B;
	LinearAlgebra::Crout LA(n);
	try{LA.Invert(A,B);} catch(int code){throw(code);} 
	return B;
}

double Norm(const Set::VectorSpace::Hom &A)
{
	unsigned int n1 = A.size1();
	unsigned int n2 = A.size2();
	if ((n1 == 0) || (n2 == 0)) return 0.0;
	if (n2 < n1)
	{
		unsigned int i; int ii; 
		double rho, rhomax=0.0;
		Set::VectorSpace::Hom B=Adjoint(A)*A;
		TNT::Array2D<double> C(n2,n2); Copy(B,C);
		JAMA::Eigenvalue<double> Eigen(C);
		TNT::Array2D<double> D(n2,n2); Eigen.getD(D);
		for (i=0; i<n2; i++)
		{
			ii = (int)i; rho=D[ii][ii];
			if (rho > rhomax) rhomax = rho;
		}
		return sqrt(rhomax);
	}
	else      
	{
		unsigned int i; int ii;
		double rho, rhomax=0.0;
		Set::VectorSpace::Hom B=A*Adjoint(A);
		TNT::Array2D<double> C(n1,n1); Copy(B,C);
		JAMA::Eigenvalue<double> Eigen(C);
		TNT::Array2D<double> D(n1,n1); Eigen.getD(D);
		for (i=0; i<n1; i++)
		{
			ii = (int)i; rho=D[ii][ii];
			if (rho > rhomax) rhomax = rho;
		}
		return sqrt(rhomax);
	}
}

double 
VectorNorm(const Set::VectorSpace::Hom &A)
{
	double N2=0.0; double *p;
	for (p=A.begin(); p!=A.end(); p++) N2 += (*p)*(*p);
	return sqrt(N2);
}

// double
// Jacobian(const Set::VectorSpace::Hom &A)
// {
// 	unsigned int n=A.size1();
// 	assert (n == A.size2());
// 	if (n == 3) return 
// 		A[0][0]*A[1][1]*A[2][2] + A[1][0]*A[2][1]*A[0][2] + 
// 		A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0] - 
// 		A[0][1]*A[1][0]*A[2][2] - A[1][2]*A[2][1]*A[0][0];
// 	if (n == 2) return A[0][0]*A[1][1] - A[0][1]*A[1][0];
// 	if (n == 1) return A[0][0];
// 	if (n == 0) return 1.0;
// 	LinearAlgebra::Crout LA(n);
// 	return LA.Determinant(A);
// }

double
Jacobian(const Set::VectorSpace::Hom &A)
{
	unsigned int n=A.size1();
	assert (n == A.size2());
	const double *p = A.begin();

	switch(n) {
	default:
	  {
	    LinearAlgebra::Crout LA(n);
	    return LA.Determinant(A);
	  }
	case 3:
	  return *p * *(p+4) * *(p+8) +  *(p+1) * *(p+5) * *(p+6) + 
	         *(p+3) * *(p+7) * *(p+2) -  *(p+6) * *(p+4) * *(p+2) - 
		 *(p+3) * *(p+1) * *(p+8) -  *(p+7) * *(p+5) * *p;
	case 2: 	  
	  return *p * *(p+3) -  *(p+1) * *(p+2);	  
	case 1:
	  return *p;
	case 0:
	  return 1.0;
	}
}

//////////////////////////////////////////////////////////////////////
// Class HomZero
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::HomZero 
operator + (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert(A.size1() == B.size1()); 
	assert(A.size2() == B.size2()); 
	return A;
}

Set::VectorSpace::Hom 
operator + (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert(A.size1() == B.size1()); 
	assert(A.size2() == B.size2()); 
	return A;
}

Set::VectorSpace::Hom 
operator + (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::Hom &B)
{
	assert(A.size1() == B.size1()); 
	assert(A.size2() == B.size2()); 
	return B;
}

Set::VectorSpace::HomZero 
operator - (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert(A.size1() == B.size1()); 
	assert(A.size2() == B.size2()); 
	return A;
}

Set::VectorSpace::Hom 
operator - (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert(A.size1() == B.size1()); 
	assert(A.size2() == B.size2()); 
	return A;
}

Set::VectorSpace::Hom 
operator - (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::Hom &B)
{
	assert(A.size1() == B.size1()); 
	assert(A.size2() == B.size2()); 
	return -B;
}

Set::VectorSpace::HomZero 
operator - (const Set::VectorSpace::HomZero &A)
{
	return A;
}

Set::VectorSpace::HomZero 
operator * (const Set::VectorSpace::HomZero &A,
			const double &p)
{
	return A;
}

Set::VectorSpace::HomZero 
operator * (const double &p,
			const Set::VectorSpace::HomZero &A)
{
	return A;
}

Set::VectorSpace::HomZero 
operator / (const Set::VectorSpace::HomZero &A, 
			const double &p)
{
	return A;
}

Set::VectorSpace::VectorZero
operator * (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::Vector &B)
{
	assert (B.size() == A.size2()); 
	Set::VectorSpace::VectorZero O(A.size1()); 
	return O;
}

Set::VectorSpace::HomZero
operator * (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert(B.size1() == A.size2());
	Set::VectorSpace::HomZero O(A.size1(),B.size2());
	return O;
}

Set::VectorSpace::HomZero
operator * (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::Hom &B)
{
	assert(B.size1() == A.size2());
	Set::VectorSpace::HomZero O(A.size1(),B.size2());
	return O;
}

Set::VectorSpace::HomZero
operator * (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert(B.size1() == A.size2());
	Set::VectorSpace::HomZero O(A.size1(),B.size2());
	return O;
}

double 
Norm(const Set::VectorSpace::HomZero &)
{
	return 0.0;
}

double 
VectorNorm(const Set::VectorSpace::HomZero &)
{
	return 0.0;
}

Set::VectorSpace::HomZero
Adjoint(const Set::VectorSpace::HomZero &A)
{
	Set::VectorSpace::HomZero B(A.size2(),A.size1());
	return B;
}

ostream & 
operator<<(ostream &os, const Set::VectorSpace::HomZero &A)
{
	A.print(&os); return os;
}

double 
Trace(const Set::VectorSpace::HomZero &)
{
	return 0.0;
}

//////////////////////////////////////////////////////////////////////
// Class HomId
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::HomId &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	Set::VectorSpace::Hom C(A.size1(),A.size2());
	unsigned int m1=A.size1()+1; double *p;
	for (p=C.begin(); p<C.end(); p+=m1) *p = 2.0;
	return C;
}

Set::VectorSpace::HomZero 
operator - (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::HomId &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	Set::VectorSpace::HomZero C(A.size1(),A.size2()); 
	return C;
}

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::HomId &A, 
			const double &b)
{
	Set::VectorSpace::Hom C(A.size1(),A.size2());
	unsigned int m1=A.size1()+1; double *p;
	for (p=C.begin(); p<C.end(); p+=m1) *p = b;
	return C;
}

Set::VectorSpace::Hom
operator * (const double &a, 
			const Set::VectorSpace::HomId &B)
{
	Set::VectorSpace::Hom C(B.size1(),B.size2());
	unsigned int m1=B.size1()+1; double *p;
	for (p=C.begin(); p<C.end(); p+=m1) *p = a;
	return C;
}

Set::VectorSpace::Hom
operator / (const Set::VectorSpace::HomId &A, 
			const double &b)
{
	Set::VectorSpace::Hom C(A.size1(),A.size2());
	unsigned int m1=A.size1()+1; 
	double c = 1.0/b; double *p;
	for (p=C.begin(); p<C.end(); p+=m1) *p = c;
	return C;
}

Set::VectorSpace::HomId 
operator + (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	return A;
}

Set::VectorSpace::HomId 
operator + (const Set::VectorSpace::HomZero &A, 
			const Set::VectorSpace::HomId &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	return B;
}

Set::VectorSpace::HomId 
operator - (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::HomZero &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	return A;
}

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::Hom &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	Set::VectorSpace::Hom C=B; C+=A; return C;
}

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::HomId &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	Set::VectorSpace::Hom C=A; C+=B; return C;
}

Set::VectorSpace::Hom
operator - (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::Hom &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	Set::VectorSpace::Hom C=B; C-=A; return -C;
}

Set::VectorSpace::Hom
operator - (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::HomId &B)
{
	assert (A.size1() == B.size1());
	assert (A.size2() == B.size2());
	Set::VectorSpace::Hom C=A; C-=B; return C;
}

Set::VectorSpace::Vector
operator * (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::Vector &B)
{
	assert (B.size() == A.size2()); 
	return B;
}

Set::VectorSpace::HomId
operator * (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::HomId &B)
{
	assert(B.size1() == A.size2());
	return A;
}

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::HomId &A, 
			const Set::VectorSpace::Hom &B)
{
	assert(B.size1() == A.size2());
	return B;
}

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::HomId &B)
{
	assert(B.size1() == A.size2());
	return A;
}

double 
Norm(const Set::VectorSpace::HomId &A)
{
	return 1.0;
}

double 
VectorNorm(const Set::VectorSpace::HomId &A)
{
	return sqrt((double)A.size1());
}

double
Jacobian(const Set::VectorSpace::HomId &A)
{
	return 1.0;
}

double 
Trace(const Set::VectorSpace::HomId &A)
{
	return ((double)A.size1());
}
