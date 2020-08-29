// ChainComplex.cpp: implementation of the ChainComplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "ChainComplex.h"

namespace Geometry
{
//////////////////////////////////////////////////////////////////////
// Class Cell
//////////////////////////////////////////////////////////////////////

Cell::Cell() : n(0) {}

Cell::Cell(const int &n_) 
: n(n_), bo(n_-1), co(n_+1) {}

Cell::~Cell() {}

Cell::Cell(const Cell &A) 
: n(A.n), bo(A.bo), co(A.co){}

Cell & Cell::operator = (const Cell &A)
{
	assert (n == A.n);
	if (this == &A) return *this;
	n = A.n; bo = A.bo; co = A.co;
	return *this;
}

const Chain &
Cell::Boundary() const
{
	return bo;
}

Chain & 
Cell::Boundary()
{
	return bo;
}

const Cochain &
Cell::Coboundary() const
{
	return co;
}

Cochain & 
Cell::Coboundary()
{
	return co;
}

void 
Cell::print(ostream *os)
{
	*os << (size_t) this;
}

const int &
Cell::dim() const
{
	return n;
}

void 
Cell::Bo2Co()
{
	Chain::iterator pbo;
	for (pbo=bo.begin(); pbo!=bo.end(); pbo++) 
	{
		Cell * const f=pbo->first;
		f->Coboundary()[this]=pbo->second;
	}
}

void 
Cell::Co2Bo()
{
	Cochain::iterator pco;
	for (pco=co.begin(); pco!=co.end(); pco++) 
	{
		Cell * const f=pco->first;
		f->Boundary()[this]=pco->second;
	}
}

//////////////////////////////////////////////////////////////////////
// Class Chain
//////////////////////////////////////////////////////////////////////

Chain::Chain() : n(0) {}

Chain::Chain(const int &n_) : n(n_) {}

Chain::~Chain(){}

Chain::Chain(const Chain &A) 
: n(A.n), map<Cell *,int>(A){}

Chain & 
Chain::operator = (const Chain &A)
{
	assert (n == A.n);
	if (this == &A) return *this;
	n = A.n; map<Cell *,int>::operator = (A);
	return *this;
}

void 
Chain::operator += (const Chain &A)
{
	assert (n == A.n);
	Chain::const_iterator p;
	for (p=A.begin(); p!=A.end(); p++) 
		(*this)[p->first] += p->second;
}

void 
Chain::operator -= (const Chain &A)
{
	assert (n == A.n);
	Chain::const_iterator p;
	for (p=A.begin(); p!=A.end(); p++) 
		(*this)[p->first] -= p->second;
}

void 
Chain::operator *= (const int &a)
{
	Chain::iterator p;
	for (p=this->begin(); p!=this->end(); p++) 
		p->second *= a;
}

int 
Chain::operator () (const Cochain &A)
{
	assert (n == A.dim());
	int product=0; Cochain B=A;
	Chain::iterator p;
	for (p=this->begin(); p!=this->end(); p++)
		product += B[p->first]*p->second;
	return product;
}

Chain 
Chain::Boundary() const
{
	Chain B(n-1);
	Chain::const_iterator p;
	for (p=this->begin(); p!=this->end(); p++) 
	{
		Chain C(n-1);
		C = p->first->Boundary(); 
		C *= p->second; B += C;
	}
	return B;
}

void 
Chain::print(ostream *os)
{
	Chain::iterator p;
	for (p=this->begin(); p!=this->end(); p++) 
		*os << (size_t) p->first << ", " << p->second << endl;
}

const int &
Chain::dim() const
{
	return n;
}

//////////////////////////////////////////////////////////////////////
// Class Cochain
//////////////////////////////////////////////////////////////////////

Cochain::Cochain() : n(0) {}

Cochain::Cochain(const int &n_) : n(n_) {}

Cochain::~Cochain() {}

Cochain::Cochain(const Cochain &A) 
: n(A.n), map<Cell *,int>(A){}

Cochain & 
Cochain::operator = (const Cochain &A)
{
	assert (n == A.n);
	if (this == &A) return *this;
	n = A.n; map<Cell *,int>::operator = (A);
	return *this;
}

void 
Cochain::operator += (const Cochain &A)
{
	assert (n == A.n);
	Cochain::const_iterator p;
	for (p=A.begin(); p!=A.end(); p++) 
		(*this)[p->first] += p->second;
}
void 
Cochain::operator -= (const Cochain &A)
{
	assert (n == A.n);
	Cochain::const_iterator p;
	for (p=A.begin(); p!=A.end(); p++) 
		(*this)[p->first] -= p->second;
}
void 
Cochain::operator *= (const int &a)
{
	Cochain::iterator p;
	for (p=this->begin(); p!=this->end(); p++) 
		p->second *= a;
}

int 
Cochain::operator () (const Chain &A)
{
	assert (n == A.dim());
	int product=0; Chain B=A;
	Cochain::iterator p;
	for (p=this->begin(); p!=this->end(); p++)
		product += B[p->first]*p->second;
	return product;
}

Cochain 
Cochain::Coboundary() const
{
	Cochain B(n+1);
	Cochain::const_iterator p; 
	for (p=this->begin(); p!=this->end(); p++) 
	{
		Cochain C(n+1);
		C = p->first->Coboundary(); 
		C *= p->second; B += C;
	}
	return B;
}

void 
Cochain::print(ostream *os)
{
	Cochain::iterator p; 
	for (p=this->begin(); p!=this->end(); p++) 
		*os << (size_t) p->first << ", " << p->second << endl;
}

const int &
Cochain::dim() const
{
	return n;
}

}

//////////////////////////////////////////////////////////////////////
// Class Cell
//////////////////////////////////////////////////////////////////////

Geometry::Chain 
Boundary(const Geometry::Cell &A)
{
	return A.Boundary();
}

Geometry::Cochain 
Coboundary(const Geometry::Cell &A)
{
	return A.Coboundary();
}


ostream & 
operator<<(ostream &os, Geometry::Cell &A)
{
	A.print(&os); return os;
}

//////////////////////////////////////////////////////////////////////
// Class Chain
//////////////////////////////////////////////////////////////////////


Geometry::Chain 
operator - (const Geometry::Chain &A)
{
	Geometry::Chain B=A; 
	Geometry::Chain::iterator p; 
	for (p=B.begin(); p!=B.end(); p++) 
		p->second = - p->second;
	return B;
}


Geometry::Chain 
operator + (const Geometry::Chain &A, 
			const Geometry::Chain &B)
{
	Geometry::Chain C=A; C += B; return C;
}


Geometry::Chain 
operator - (const Geometry::Chain &A, 
			const Geometry::Chain &B)
{
	Geometry::Chain C=A; C -= B; return C;
}


Geometry::Chain 
operator * (const Geometry::Chain &A, 
			const int &a)
{
	Geometry::Chain C=A; C *= a; return C;
}


Geometry::Chain 
operator * (const int &a, 
			const Geometry::Chain &A)
{
	Geometry::Chain C=A; C *= a; return C;
}


Geometry::Chain 
Boundary(Geometry::Chain &A)
{
	Geometry::Chain B=A.Boundary(); return B;
}


void 
Random(Geometry::Chain &A)
{
	Geometry::Chain::iterator p; 
	for (p=A.begin(); p!=A.end(); p++) 
		p->second = 2*rand() - RAND_MAX;
}


ostream & 
operator<<(ostream &os, Geometry::Chain &A)
{
	A.print(&os); return os;
}

//////////////////////////////////////////////////////////////////////
// Class Cochain
//////////////////////////////////////////////////////////////////////


Geometry::Cochain 
operator - (const Geometry::Cochain &A)
{
	Geometry::Cochain B=A; 
	Geometry::Cochain::iterator p; 
	for (p=B.begin(); p!=B.end(); p++) 
		p->second = - p->second;
	return B;
}


Geometry::Cochain 
operator + (const Geometry::Cochain &A, 
			const Geometry::Cochain &B)
{
	Geometry::Cochain C=A; C += B; return C;
}


Geometry::Cochain 
operator - (const Geometry::Cochain &A, 
			const Geometry::Cochain &B)
{
	Geometry::Cochain C=A; C -= B; return C;
}


Geometry::Cochain 
operator * (const Geometry::Cochain &A, 
			const int &a)
{
	Geometry::Cochain C=A; C *= a; return C;
}


Geometry::Cochain 
operator * (const int &a, 
			const Geometry::Cochain &A)
{
	Geometry::Cochain C=A; C *= a; return C;
}


Geometry::Cochain 
Coboundary(Geometry::Cochain &A)
{
	Geometry::Cochain B=A.Coboundary(); return B;
}


void 
Random(Geometry::Cochain &A)
{
	Geometry::Cochain::iterator p;
	for (p=A.begin(); p!=A.end(); p++) 
		p->second = 2*rand() - RAND_MAX;
}


ostream & 
operator<<(ostream &os, Geometry::Cochain &A)
{
	A.print(&os); return os;
}
