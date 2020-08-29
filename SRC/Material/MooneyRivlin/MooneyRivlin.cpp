// MooneyRivlin.cpp: implementation of the MooneyRivlin class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./MooneyRivlin.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Material
{
namespace MooneyRivlin
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() : K(0.0), Mu1(0.0), Mu2(0.0) {}

Data::Data(const double &rhs_K, const double &rhs_Mu1, const double &rhs_Mu2) 
	: K(rhs_K), Mu1(rhs_Mu1), Mu2(rhs_Mu2) {}

Data::~Data(){}

Data::Data(const Data &rhs) 
	: K(rhs.K), Mu1(rhs.Mu1), Mu2(rhs.Mu2) {}

Data &
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	K = rhs.K; Mu1 = rhs.Mu1; Mu2 = rhs.Mu2;
	return *this;
}

const double &
Data::GetK() const
{
	return K;
}

double & 
Data::GetK()
{
	return K;
}

const double &
Data::GetMu1() const
{
	return Mu1;
}

double & 
Data::GetMu1()
{
	return Mu1;
}

const double &
Data::GetMu2() const
{
	return Mu2;
}

double & 
Data::GetMu2()
{
	return Mu2;
}

void 
Data::Randomize()
{
	K = (double)rand()/(double)RAND_MAX;
	Mu1 = (double)rand()/(double)RAND_MAX;
        Mu2 = (double)rand()/(double)RAND_MAX;
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
	const double &K = LS->Properties->K;
	const double &Mu1 = LS->Properties->Mu1;
        const double &Mu2 = LS->Properties->Mu2;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	Set::VectorSpace::Hom Fdev = F/pow(J,1.0/3.0);
	Set::VectorSpace::Hom Cdev = Adjoint(Fdev)*Fdev;
        Set::VectorSpace::Hom CdevCdev = Cdev*Cdev;
        double I1 = Trace(Cdev);
        double I2 = 0.5*I1*I1 - 0.5*Trace(CdevCdev);
	return 0.5*Mu1*(I1 - 3.0) + 0.5*Mu2*(I2 - 3.0) + 0.5*K*pow(J-1,2);
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
	const double &K = LS->Properties->K;
	const double &Mu1 = LS->Properties->Mu1;
        const double &Mu2 = LS->Properties->Mu2;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
        Set::VectorSpace::Hom Fdev = F/pow(J,1.0/3.0);
	Set::VectorSpace::Hom Cdev = Adjoint(Fdev)*Fdev;
        Set::VectorSpace::Hom CdevCdev = Cdev*Cdev;
        double I1 = Trace(Cdev);
        double I2 = 0.5*I1*I1 - 0.5*Trace(CdevCdev);
	double fact1 = Mu1/pow(J,2.0/3.0);
        double fact2 = Mu2/pow(J,4.0/3.0);
	Set::VectorSpace::Hom C = Adjoint(F)*F;
	Set::VectorSpace::Hom Finv = Inverse(F);
	return fact1*(F - (Trace(C)/3.0)*Adjoint(Finv)) - 2.0/3.0*Mu2*I2*Adjoint(Finv) + fact2*Trace(C)*F - fact2*F*Adjoint(F)*F + K*(J-1)*J*Adjoint(Finv);
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
        const double &K = LS->Properties->K;
	const double &Mu1 = LS->Properties->Mu1;
        const double &Mu2 = LS->Properties->Mu2;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
        Set::VectorSpace::Hom Fdev = F/pow(J,1.0/3.0);
	Set::VectorSpace::Hom Cdev = Adjoint(Fdev)*Fdev;
        Set::VectorSpace::Hom CdevCdev = Cdev*Cdev;
        double I1 = Trace(Cdev);
        double I2 = 0.5*I1*I1 - 0.5*Trace(CdevCdev);
	Set::VectorSpace::Hom Finv = Inverse(F);
	Set::VectorSpace::Hom C = Adjoint(F)*F;
        Set::VectorSpace::Hom BF = F*C;
        Set::VectorSpace::Vector Finvt(9);
        Finvt = Adjoint(Finv);
	
      	unsigned int i, j, k, l, q;
	double fact = Mu1/pow(J,2.0/3.0);
        double fact_ = Mu2/pow(J,4.0/3.0);
	double t = Trace(C)/3.0;

        Set::VectorSpace::Hom DDW(9);
        DDW = K*(2*J-1)*J*Dyadic(Finvt,Finvt);
        DDW += 8.0/9.0*Mu2*I2*Dyadic(Finvt,Finvt);
        DDW += fact_*2*Dyadic(Dy,Dy);
        DDW -= 4.0*fact_*t*(Dyadic(Dy,Finvt)+Dyadic(Finvt,Dy)); 
        DDW += 4.0/3.0*fact_*(Dyadic(BF,Finvt)+Dyadic(Finvt,BF));
	double **** D = Indexing::New(DDW.begin(),3,3,3,3);      
        
        for (j=0; j<3; j++){
           for (k=0; k<3; k++){
               double fact2 = K*(J-1)*J*Finv[k][j];
               for (i=0; i<3; i++){
                    for (l=0; l<3; l++){ 
                        D[l][k][j][i] -= fact2*Finv[i][l];
                    }
                }
            }
        }
        
        
	{
		double fact2 = fact*t, fact3;

		for (k=0; k<3; k++) 
			for (j=0; j<3; j++)
			{
				fact3 = fact2*Finv[k][j];
				for (i=0; i<3; i++) 
					for (l=0; l<3; l++)
						D[l][k][j][i] += fact3*Finv[i][l];
			}
	}

	{
		double fact2 = fact*(2.0/3.0)*t, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] += fact3*Finv[i][j];
			}
	}

	{
		double fact2 = (2.0/3.0)*fact, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*F[j][i];
			}

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*F[l][k];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*Finv[i][j];
			}
	}

	for (k=0; k<3; k++)
		for (i=0; i<3; i++) 
			D[k][i][k][i] += fact + fact_*3.0*t;
        
       for (j=0; j<3; j++){
           for (k=0; k<3; k++){
               double fact2 = 2.0/3.0*Mu2*I2*Finv[k][j];
               for (i=0; i<3; i++){
                    for (l=0; l<3; l++){ 
                        D[l][k][j][i] += fact2*Finv[i][l];
                    }
                }
            }
        }

        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
		for (l=0; l<3; l++)
                    for (q=0; q<3; q++)
			D[l][i][j][i] -= fact_*F[j][q]*F[l][q];
        
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
		for (k=0; k<3; k++)
                    for (l=0; l<3; l++)
			D[l][k][j][i] -= fact_*F[j][k]*F[l][i];
        
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
		for (k=0; k<3; k++)
                    for (q=0; q<3; q++)
			D[j][k][j][i] -= fact_*F[q][k]*F[q][i];
        
        
	Indexing::Delete(D,3,3,3,3);
                
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
	const double &K = LS->Properties->K;
	const double &Mu1 = LS->Properties->Mu1;
        const double &Mu2 = LS->Properties->Mu2;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
        Set::VectorSpace::Hom Fdev = F/pow(J,1.0/3.0);
	Set::VectorSpace::Hom Cdev = Adjoint(Fdev)*Fdev;
        Set::VectorSpace::Hom CdevCdev = Cdev*Cdev;
        double I1 = Trace(Cdev);
        double I2 = 0.5*I1*I1 - 0.5*Trace(CdevCdev);
	double fact1 = Mu1/pow(J,2.0/3.0);
        double fact2 = Mu2/pow(J,4.0/3.0);
	Set::VectorSpace::Hom C = Adjoint(F)*F;
	Set::VectorSpace::Hom Finv = Inverse(F);
        double W = 0.5*Mu1*(I1 - 3.0) + 0.5*Mu2*(I2 - 3.0) + 0.5*K*pow(J-1,2);
        Set::VectorSpace::Vector DW(9);
        DW = fact1*(F - (Trace(C)/3.0)*Adjoint(Finv)) - 2.0/3.0*Mu2*I2*Adjoint(Finv) + fact2*Trace(C)*F - fact2*F*Adjoint(F)*F + K*(J-1)*J*Adjoint(Finv);
        return make_pair(W,DW);
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
        const double &K = LS->Properties->K;
	const double &Mu1 = LS->Properties->Mu1;
        const double &Mu2 = LS->Properties->Mu2;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
        Set::VectorSpace::Hom Fdev = F/pow(J,1.0/3.0);
	Set::VectorSpace::Hom Cdev = Adjoint(Fdev)*Fdev;
        Set::VectorSpace::Hom CdevCdev = Cdev*Cdev;
        double I1 = Trace(Cdev);
        double I2 = 0.5*I1*I1 - 0.5*Trace(CdevCdev);
	Set::VectorSpace::Hom Finv = Inverse(F);
	Set::VectorSpace::Hom C = Adjoint(F)*F;
        Set::VectorSpace::Hom BF = F*C;
        Set::VectorSpace::Vector Finvt(9);
        Finvt = Adjoint(Finv);
	
      	unsigned int i, j, k, l, q;
	double fact = Mu1/pow(J,2.0/3.0);
        double fact_ = Mu2/pow(J,4.0/3.0);
	double t = Trace(C)/3.0;

        Set::VectorSpace::Vector DW(9);
        DW = fact*(F - (Trace(C)/3.0)*Adjoint(Finv)) - 2.0/3.0*Mu2*I2*Adjoint(Finv) + fact_*Trace(C)*F - fact_*F*Adjoint(F)*F + K*(J-1)*J*Adjoint(Finv);

        Set::VectorSpace::Hom DDW(9);
        DDW = K*(2*J-1)*J*Dyadic(Finvt,Finvt);
        DDW += 8.0/9.0*Mu2*I2*Dyadic(Finvt,Finvt);
        DDW += fact_*2*Dyadic(Dy,Dy);
        DDW -= 4.0*fact_*t*(Dyadic(Dy,Finvt)+Dyadic(Finvt,Dy)); 
        DDW += 4.0/3.0*fact_*(Dyadic(BF,Finvt)+Dyadic(Finvt,BF));
	double **** D = Indexing::New(DDW.begin(),3,3,3,3);      
        
        for (j=0; j<3; j++){
           for (k=0; k<3; k++){
               double fact2 = K*(J-1)*J*Finv[k][j];
               for (i=0; i<3; i++){
                    for (l=0; l<3; l++){ 
                        D[l][k][j][i] -= fact2*Finv[i][l];
                    }
                }
            }
        }
        
        
	{
		double fact2 = fact*t, fact3;

		for (k=0; k<3; k++) 
			for (j=0; j<3; j++)
			{
				fact3 = fact2*Finv[k][j];
				for (i=0; i<3; i++) 
					for (l=0; l<3; l++)
						D[l][k][j][i] += fact3*Finv[i][l];
			}
	}

	{
		double fact2 = fact*(2.0/3.0)*t, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] += fact3*Finv[i][j];
			}
	}

	{
		double fact2 = (2.0/3.0)*fact, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*F[j][i];
			}

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*F[l][k];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*Finv[i][j];
			}
	}

	for (k=0; k<3; k++)
		for (i=0; i<3; i++) 
			D[k][i][k][i] += fact + fact_*3.0*t;
        
       for (j=0; j<3; j++){
           for (k=0; k<3; k++){
               double fact2 = 2.0/3.0*Mu2*I2*Finv[k][j];
               for (i=0; i<3; i++){
                    for (l=0; l<3; l++){ 
                        D[l][k][j][i] += fact2*Finv[i][l];
                    }
                }
            }
        }

        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
		for (l=0; l<3; l++)
                    for (q=0; q<3; q++)
			D[l][i][j][i] -= fact_*F[j][q]*F[l][q];
        
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
		for (k=0; k<3; k++)
                    for (l=0; l<3; l++)
			D[l][k][j][i] -= fact_*F[j][k]*F[l][i];
        
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
		for (k=0; k<3; k++)
                    for (q=0; q<3; q++)
			D[j][k][j][i] -= fact_*F[q][k]*F[q][i];
        
        
	Indexing::Delete(D,3,3,3,3);

	return make_pair(DW,DDW);
}

}

}
