// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#include "./TwoStages.h"

namespace Material
{
namespace TwoStages
{

// Class LocalState

LocalState::LocalState() : LScomp(0), LSten(0){}

LocalState::~LocalState() {}

LocalState::LocalState(
		       Material::LocalState *LScomp_,
		       Material::LocalState *LSten_) : 
  LScomp(LScomp_), LSten(LSten_) {}

Material::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(const LocalState &rhs) :
  LScomp(rhs.LScomp), LSten(rhs.LSten) {}

void 
LocalState::operator ++ () {
  ++(*LScomp);
  ++(*LSten);
}

// Class Energy<0>

Energy<0>::Energy() {}

Energy<0>::~Energy(){}

Material::Energy<0> *
Energy<0>::Clone() const
{
      return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState * LS_,
		  Material::Energy<0> * Wcomp_, 
		  Material::Energy<0> * Wten_) : 
  LS(LS_), Wcomp(Wcomp_), Wten(Wten_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS), Wcomp(rhs.Wcomp), Wten(rhs.Wten) {}

double
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
  Set::VectorSpace::Hom F(3,3,Dy.begin());
  double J = Jacobian(F);
  if ( J > 1.0 )
  {
    return (*Wten)(Dy);
  }
  else
  {
    return (*Wcomp)(Dy);
  }
}
  

// Class Energy<1>

Energy<1>::Energy() {}

Energy<1>::~Energy(){}

Material::Energy<1> *
Energy<1>::Clone() const
{
        return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState * LS_,
		  Material::Energy<1> * DWcomp_, 
		  Material::Energy<1> * DWten_) : 
  LS(LS_), DWcomp(DWcomp_), DWten(DWten_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS), DWcomp(rhs.DWcomp), DWten(rhs.DWten) {}

Set::VectorSpace::Vector
Energy<1>::operator ()  (const Set::VectorSpace::Vector &Dy) const
{
  Set::VectorSpace::Hom F(3,3,Dy.begin());
  double J = Jacobian(F);
  if ( J > 1.0 )
  { 
    return (*DWten)(Dy);
  }
  else
  {
    return (*DWcomp)(Dy);
  }
}



// Class Energy<2>

Energy<2>::Energy() {}

Energy<2>::~Energy(){}

Material::Energy<2> *
Energy<2>::Clone() const
{
        return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState * LS_,		   
		  Material::Energy<2> * DDWcomp_, 
		  Material::Energy<2> * DDWten_) : 
  LS(LS_), DDWcomp(DDWcomp_), DDWten(DDWten_) {}

Energy<2>::Energy(const Energy<2> &rhs) : 
  LS(rhs.LS), DDWcomp(rhs.DDWcomp), DDWten(rhs.DDWten) {}

Set::VectorSpace::Hom 
Energy<2>::operator ()  (const Set::VectorSpace::Vector &Dy) const
{
  Set::VectorSpace::Hom F(3,3,Dy.begin());
  double J = Jacobian(F);
  if ( J > 1.0 )
  { 
    return (*DDWten)(Dy);
  }
  else
  {
    return (*DDWcomp)(Dy);
  }
}


// Class Jet<0>

Jet<0>::Jet() {}

Jet<0>::~Jet() {}

Material::Jet<0> *
Jet<0>::Clone() const
{
        return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState * LS_,
	    Material::Jet<0> * Jcomp_, 
	    Material::Jet<0> * Jten_) : 
  LS(LS_), Jcomp(Jcomp_), Jten(Jten_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS), Jcomp(rhs.Jcomp), Jten(rhs.Jten) {}

pair<double, Set::VectorSpace::Vector>
Jet<0>::operator ()  (const Set::VectorSpace::Vector &Dy) const
{
  Set::VectorSpace::Hom F(3,3,Dy.begin());
  double J = Jacobian(F);
  if ( J > 1.0 )
  {
    return (*Jten)(Dy);
  }
  else
  {
    return (*Jcomp)(Dy);
  }
}



// Class Jet<1>

Jet<1>::Jet() {}

Jet<1>::~Jet() {}

Material::Jet<1> *
Jet<1>::Clone() const
{
         return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState * LS_,
	    Material::Jet<1> * DJcomp_, 
	    Material::Jet<1> * DJten_) : 
  LS(LS_), DJcomp(DJcomp_), DJten(DJten_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS), DJcomp(rhs.DJcomp), DJten(rhs.DJten) {}

pair<Set::VectorSpace::Vector, Set::VectorSpace::Hom> 
Jet<1>::operator ()  (const Set::VectorSpace::Vector &Dy) const
{
  Set::VectorSpace::Hom F(3,3,Dy.begin());
  double J = Jacobian(F);
  if ( J > 1.0 )
  {
    return (*DJten)(Dy);
  }
  else
  {
    return (*DJcomp)(Dy); 
  }
}


}

}


