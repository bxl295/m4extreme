// NeoHookean.cpp: implementation of the NeoHookean class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./NeoHookean.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Material
{
  namespace NeoHookean
  {
    //////////////////////////////////////////////////////////////////////
    // Class Data
    //////////////////////////////////////////////////////////////////////

    Data::Data() : Lambda(0.0), Mu(0.0) {}

    Data::Data(const double &rhs_Lambda, const double &rhs_Mu) 
      : Lambda(rhs_Lambda), Mu(rhs_Mu) {}

    Data::~Data(){}

    Data::Data(const Data &rhs) 
      : Lambda(rhs.Lambda), Mu(rhs.Mu) {}

    Data &
    Data::operator = (const Data &rhs)
    {
      if (this == &rhs) return *this;
      Lambda = rhs.Lambda; Mu = rhs.Mu;
      return *this;
    }

    const double &
    Data::GetLambda() const
    {
      return Lambda;
    }

    double & 
    Data::GetLambda()
    {
      return Lambda;
    }

    const double &
    Data::GetMu() const
    {
      return Mu;
    }

    double & 
    Data::GetMu()
    {
      return Mu;
    }

    void 
    Data::Randomize()
    {
      Lambda = (double)rand()/(double)RAND_MAX;
      Mu     = (double)rand()/(double)RAND_MAX;
    }

    //////////////////////////////////////////////////////////////////////
    // Class LocalState
    //////////////////////////////////////////////////////////////////////

    LocalState::~LocalState() {}

    Material::LocalState *
    LocalState::Clone() const 
    {
      return new LocalState(*this);
    }

    LocalState::LocalState(Data *Properties_)
      : Properties(Properties_) {}

    LocalState::LocalState(const LocalState &rhs) : 
      Properties(rhs.Properties) {}

    void 
    LocalState::operator ++ () {}

    //////////////////////////////////////////////////////////////////////
    // Class Energy<0>
    //////////////////////////////////////////////////////////////////////

    Energy<0>::~Energy() {}

    Material::Energy<0> *
    Energy<0>::Clone() const
    {
      return new Energy<0>(*this);
    }

    Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

    Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

    double 
    Energy<0>::operator () (const Set::VectorSpace::Vector &Dy) const
    {
      assert (Dy.size() == 9);
      const double &Lambda = LS->Properties->Lambda;
      const double &Mu = LS->Properties->Mu;
      Set::VectorSpace::Hom F(3,3,Dy.begin());
      Set::VectorSpace::Hom C = Adjoint(F)*F;
      double J = Jacobian(F);
      if ( J > 1.0e-4 ) {
	double Theta = log(J);
	return 0.5*Lambda*pow(Theta,2) - Mu*Theta + 0.5*Mu*(Trace(C) - 3.0);
      }
      else {
	return 0.5*Mu*(Trace(C) - 3.0); 
      }
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

    Energy<1>::Energy(LocalState *LS_) : LS(LS_) {}

    Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS) {}

    Set::VectorSpace::Vector 
    Energy<1>::operator () (const Set::VectorSpace::Vector &Dy) const
    {
      assert (Dy.size() == 9);
      const double &Lambda = LS->Properties->Lambda;
      const double &Mu = LS->Properties->Mu;
      Set::VectorSpace::Hom F(3,3,Dy.begin());
      double J = Jacobian(F);

      if ( J > 1.0e-4 ) {
	double Theta = log(J);
	Set::VectorSpace::Hom Finv = Inverse(F);
	return (Lambda*Theta - Mu)*Adjoint(Finv) + Mu*F;
      }
      else {
	return Mu*F;
      }
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

    Energy<2>::Energy(LocalState *LS_) : LS(LS_) {}

    Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS) {}

    Set::VectorSpace::Hom 
    Energy<2>::operator () (const Set::VectorSpace::Vector &Dy) const
    {
      assert (Dy.size() == 9);
      const double &Lambda = LS->Properties->Lambda;
      const double &Mu = LS->Properties->Mu;
      Set::VectorSpace::Hom F(3,3,Dy.begin());
      double J = Jacobian(F);

      Set::VectorSpace::Hom DDW(9);
      double **** C = Indexing::New(DDW.begin(),3,3,3,3);

      if ( J > 1.0e-4 ) {
	double Theta = log(J);
	Set::VectorSpace::Hom Finv = Inverse(F);

	unsigned int i, j, k, l;
	double factor, tangent;
	  
	for (i=0; i<3; i++)
	  for (j=0; j<3; j++) 
	    {
	      factor = Lambda*Finv[i][j];
	      for (k=0; k<3; k++)
		for (l=0; l<3; l++) 
		  C[l][k][j][i] = factor*Finv[k][l];
	    }
	  
	tangent = Mu - Lambda*Theta;
	  
	for (k=0; k<3; k++) 
	  for (j=0; j<3; j++)
	    {
	      factor = tangent*Finv[k][j];
	      for (i=0; i<3; i++) 
		for (l=0; l<3; l++)
		  C[l][k][j][i] += factor*Finv[i][l];
	    }
	  
	for (k=0; k<3; k++)
	  for (i=0; i<3; i++) 
	    C[k][i][k][i] += Mu;
      }
      else {
	for (int k=0; k<3; k++)
	  for (int i=0; i<3; i++) 
	    C[k][i][k][i] = Mu;    
      }

      Indexing::Delete(C,3,3,3,3);
      return DDW;
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

    Jet<0>::Jet(LocalState *LS_) : LS(LS_) {}

    Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

    pair<double,Set::VectorSpace::Vector>
    Jet<0>::operator () (const Set::VectorSpace::Vector &Dy) const
    {
      assert (Dy.size() == 9);
      const double &Lambda = LS->Properties->Lambda;
      const double &Mu = LS->Properties->Mu;
      Set::VectorSpace::Hom F(3,3,Dy.begin());
      Set::VectorSpace::Hom C = Adjoint(F)*F;

      double J = Jacobian(F);

      if ( J > 1.0e-4 ) {
	double Theta = log(J);
	Set::VectorSpace::Hom Finv = Inverse(F);
	double W = 0.5*Lambda*pow(Theta,2) - Mu*Theta + 0.5*Mu*(Trace(C) - 3.0);
      
	Set::VectorSpace::Vector DW = (Lambda*Theta - Mu)*Adjoint(Finv) + Mu*F;
	return make_pair(W,DW);
      }
      else {
	return make_pair(0.5*Mu*(Trace(C) - 3.0), Mu*F);
      }
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

    Jet<1>::Jet(LocalState *LS_) : LS(LS_) {}

    Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

    pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
    Jet<1>::operator () (const Set::VectorSpace::Vector &Dy) const
    {
      assert (Dy.size() == 9);
      const double &Lambda = LS->Properties->Lambda;
      const double &Mu = LS->Properties->Mu;
      Set::VectorSpace::Hom F(3,3,Dy.begin());
      double J = Jacobian(F);

      Set::VectorSpace::Hom DDW(9);
      double **** C = Indexing::New(DDW.begin(),3,3,3,3);

      if ( J > 1.0e-4 ) {
	double Theta = log(J);
	Set::VectorSpace::Hom Finv = Inverse(F);
	Set::VectorSpace::Vector DW = (Lambda*Theta - Mu)*Adjoint(Finv) + Mu*F;
	
	unsigned int i, j, k, l;
	double factor, tangent;
	
	for (i=0; i<3; i++)
	  for (j=0; j<3; j++) 
	    {
	      factor = Lambda*Finv[i][j];
	      for (k=0; k<3; k++)
		for (l=0; l<3; l++) 
		  C[l][k][j][i] = factor*Finv[k][l];
	    }
	
	tangent = Mu - Lambda*Theta;
	
	for (k=0; k<3; k++) 
	  for (j=0; j<3; j++)
	    {
	      factor = tangent*Finv[k][j];
	      for (i=0; i<3; i++) 
		for (l=0; l<3; l++)
		  C[l][k][j][i] += factor*Finv[i][l];
	    }
	
	for (k=0; k<3; k++)
	  for (i=0; i<3; i++) 
	    C[k][i][k][i] += Mu;

	Indexing::Delete(C,3,3,3,3);
	
	return make_pair(DW,DDW);
      }
      else {
	for (int k=0; k<3; k++)
	  for (int i=0; i<3; i++) 
	    C[k][i][k][i] = Mu; 

	Indexing::Delete(C,3,3,3,3);
	
	return make_pair(Mu*F, DDW);   
      }
    }

  }

}
