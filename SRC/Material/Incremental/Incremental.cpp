// Incremental.cpp: implementation of the Incremental class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./Incremental.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Material
{
namespace Incremental
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Material::LocalState *LS_) :
	LS(LS_), Fold(Set::VectorSpace::HomId(3)),
        F(Set::VectorSpace::HomId(3)), Sigma(Set::VectorSpace::Hom(3)) {}

LocalState::LocalState(const LocalState &rhs) : 
	Sigma(rhs.Sigma),
	LS(rhs.LS), Fold(rhs.Fold), F(rhs.F){}

void 
LocalState::operator ++ () 
{
    Fold = F;
    ++(*LS);
}

void 
LocalState::Reset(double T_) {
    LS->Reset(T_);
}

void
LocalState::write(ostream & os) const
{
    Fold.write(os);
    F.write(os);

    Sigma.write(os);

    LS->write(os);
}

void
LocalState::read(istream & is)
{
    Fold.read(is);
    F.read(is);

    Sigma.read(is);

    LS->read(is);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_, Material::Energy<0> *W_) : LS(LS_), W(W_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS), W(rhs.W) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
        
        LS->F = Finc*Fold;       

        double J = Jacobian(LS->F);
	if ( J != J ) throw(0);

#ifdef _M4EXTREME_DEBUG_II
        if ( J < 1.0e-3 || J > 1.0e3 || J != J ) {
            cerr << "Material::Incremental::Energy<0> received nonsense input J = " << J << endl;
#ifdef _M4EXTREME_DEBUG_III
            cerr << "Fold=[" << Fold
                 << "]\tFinc=[" << Finc
                 << "]\tF=[" << LS->F
                 << "]\tJ=" << J << std::endl;
            throw(0);
#endif
        }
#endif

	try {
	  return (*W)(LS->F);
	}
        catch(...) {
#ifdef _M4EXTREME_DEBUG_II
            cerr << "Material::Incremental::Energy<0> received exceptions from constitutive model" << endl;
#ifdef _M4EXTREME_DEBUG_III
            cerr << "Fold=["<<Fold<<"]\tFinc=["<<Finc<<"]\tF=["<<LS->F<<"]"<<std::endl;
#endif
	    throw(0);
#endif
        }        	  
}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy, double T) const
{
    	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
      
        return (*W)(Finc*Fold, T);        
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy() {}

Material::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_, Material::Energy<1> *DW_) : LS(LS_), DW(DW_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS), DW(rhs.DW) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
	LS->F = Finc*Fold;

        double J = Jacobian(LS->F);

        if ( J < 1.0e-4 || J > 1.0e4 || J != J ) {
            cerr << "Material::Incremental::Energy<1> received nonsense input Jacobian=" << J << endl;
#ifdef _M4EXTREME_DEBUG_III
            cerr << "Fold=[" << Fold
                 << "]\tFinc=[" << Finc
                 << "]\tF=[" << LS->F
                 << std::endl;            
#endif
            throw(0);
        }
        
#ifdef _M4EXTREME_DEBUG_II
        std::cout<<"\nFold=["<<Fold<<"]\tFinc=["<<Finc<<"]\tF=["<<LS->F<<"]"<<std::endl;
#endif
	
        try {
	  Set::VectorSpace::Hom P(3);
	  P = (*DW)(LS->F);

#ifdef _M4EXTREME_DEBUG_II       
	  std::cout<<"\nP=[" << P << "]"<<std::endl;
#endif        

	  LS->Sigma = P * Adjoint(LS->F);
	  LS->Sigma /= J;

	  return P * Adjoint(Fold);
        }
        catch(...) {
            cerr << "Material::Incremental::Energy<1> received exceptions from constitutive model" << endl;
#ifdef _M4EXTREME_DEBUG_II
            cerr << "Fold=["<<Fold<<"]\tFinc=["<<Finc<<"]\tF=["<<LS->F<<"]"<<std::endl;
#endif
	    throw(0);
        }        

}

double
Energy<1>::operator () (const Set::VectorSpace::Vector &Dy, double T) const
{
	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
            
        return (*DW)(Finc*Fold, T);        
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy() {}

Material::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_, Material::Energy<2> *DDW_) : 
	LS(LS_), DDW(DDW_) {}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS), DDW(rhs.DDW) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
	Set::VectorSpace::Hom F = Finc*Fold;
        
        LS->F = F;
        
	Set::VectorSpace::Hom DDWDFincDFinc(9);
	{
		unsigned int i, j, k, l, m, n;
		Set::VectorSpace::Hom DDWDFDF = (*DDW)(F);
		double ****C = Indexing::New(DDWDFincDFinc.begin(),3,3,3,3);
		double ****D = Indexing::New(DDWDFDF.begin(),3,3,3,3);
		for (l=0; l<3; l++)
			for (k=0; k<3; k++)
				for (j=0; j<3; j++)
					for (i=0; i<3; i++)
						for (m=0; m<3; m++)
							for (n=0; n<3; n++)
								C[l][k][j][i] += D[n][k][m][i]*Fold[m][j]*Fold[n][l];
		Indexing::Delete(C,3,3,3,3);
		Indexing::Delete(D,3,3,3,3);
	}

	return DDWDFincDFinc;
}

double
Energy<2>::operator () (const Set::VectorSpace::Vector &Dy, double T) const
{
	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
        
        return (*DDW)(Finc*Fold, T);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Material::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_, Material::Jet<0> *J_) : LS(LS_), J(J_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS), J(rhs.J) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
	Set::VectorSpace::Hom F = Finc*Fold;
	pair<double,Set::VectorSpace::Vector> K = (*J)(F);
	Set::VectorSpace::Hom P(3); P = K.second;
	return make_pair(K.first,P*Adjoint(Fold));
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Material::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_, Material::Jet<1> *DJ_) : LS(LS_), DJ(DJ_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS), DJ(rhs.DJ) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9); 
	Set::VectorSpace::Hom Finc(3,Dy.begin());
	const Set::VectorSpace::Hom &Fold = LS->Fold;
	Set::VectorSpace::Hom F = Finc*Fold;
	pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> DK = (*DJ)(F);
	Set::VectorSpace::Hom P(3); P = DK.first;
	Set::VectorSpace::Hom DWDF = P*Adjoint(Fold);
	Set::VectorSpace::Hom DDWDFincDFinc(9);
	{
		unsigned int i, j, k, l, m, n;
		Set::VectorSpace::Hom DDWDFDF = DK.second;
		double ****C = Indexing::New(DDWDFincDFinc.begin(),3,3,3,3);
		double ****D = Indexing::New(DDWDFDF.begin(),3,3,3,3);
		for (l=0; l<3; l++)
			for (k=0; k<3; k++)
				for (j=0; j<3; j++)
					for (i=0; i<3; i++)
						for (m=0; m<3; m++)
							for (n=0; n<3; n++)
								C[l][k][j][i] += D[n][k][m][i]*Fold[m][j]*Fold[n][l];
		Indexing::Delete(C,3,3,3,3);
		Indexing::Delete(D,3,3,3,3);
	}
	return make_pair(DWDF,DDWDFincDFinc);
}

}

}
