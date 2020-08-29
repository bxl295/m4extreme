// Legendre.cpp: implementation of the Legendre class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "./Legendre.h"

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class Legendre<0>
//////////////////////////////////////////////////////////////////////

Legendre<0>::Legendre() {}

Legendre<0>::~Legendre() {}

Legendre<0>::Legendre(const unsigned &n_) : n(n_) {}

vector<double> 
Legendre<0>::operator () (const double &x) const
{
//	return P[l], l=0,...,n

	vector<double> P;

	P.push_back(1.0);
	if (n == 0) return P;

	P.push_back(x);
	if (n == 1) return P;

	unsigned int l;
	for (l=2; l<=n; l++)
	{
		unsigned int l1 = l-1;
		unsigned int l2 = l-2;
		double r = (double)(l);
		double fact1 = (2.0*r-1.0)/r;
		double fact2 = (r-1.0)/r;
		P.push_back(fact1*x*P[l1]-fact2*P[l2]);
	}

	return P;
}

//////////////////////////////////////////////////////////////////////
// Class Legendre<1>
//////////////////////////////////////////////////////////////////////

Legendre<1>::Legendre() {}

Legendre<1>::~Legendre() {}

Legendre<1>::Legendre(const unsigned &n_) : n(n_) {}

vector<double> 
Legendre<1>::operator () (const double &x) const
{
 // return P'[l], l=0,...,n

	vector<double> DP;

	DP.push_back(0.0);
	if (n == 0) return DP;

	DP.push_back(1.0);
	if (n == 1) return DP;

	vector<double> P;
	P.push_back(1.0);
	P.push_back(x);

	unsigned int l;
	for (l=2; l<=n; l++)
	{
		unsigned int l1 = l-1;
		unsigned int l2 = l-2;
		double r = (double)(l);
		double fact1 = (2.0*r-1.0)/r;
		double fact2 = (r-1.0)/r;
		P.push_back(fact1*x*P[l1]-fact2*P[l2]);
		DP.push_back(fact1*P[l1] + fact1*x*DP[l1]-fact2*DP[l2]);
	}

	return DP;
}

//////////////////////////////////////////////////////////////////////
// Class AssociatedLegendre<0>
//////////////////////////////////////////////////////////////////////

AssociatedLegendre<0>::AssociatedLegendre() {}

AssociatedLegendre<0>::~AssociatedLegendre() {}

AssociatedLegendre<0>::AssociatedLegendre(const unsigned &n_) : n(n_) {}

vector<map<int,double> >
AssociatedLegendre<0>::operator () (const double &x) const
{
//	return P[l][m], l=0,...,n, m=-l,...,l

	vector<map<int,double> > P;

        //P[0][0]
        map<int,double> P0map;        
        P0map.insert(pair<int,double>(0,1));
        P.push_back(P0map);
        if (n == 0) return P;
        
        unsigned int l, m;
        //P[l][l]
        for(l=1; l<=n; l++){
            int fact1 = 1;
            if(l % 2) fact1 = -1;
            map<int,double> PLmap;
            double Pll = (double)fact1*DoubleFactorial(2*l-1)*pow(1-x*x,0.5*l);
            PLmap.insert(pair<int,double>(l,Pll));
            P.push_back(PLmap);
        }
        //P[l+1][l]
        for(l=0; l<n; l++){
            P[l+1].insert(pair<int,double>(l,x*(2*l+1)*P[l][l]));
        }
        //P[l][m]
        for(m=0; m<=n-2 && n>1; m++){
            for(l=m+2; l<=n; l++){
                unsigned int l1 = l-1;
                unsigned int l2 = l-2;
                double fact1 = (2.0*l-1.0)/((double)(l-m));
                double fact2 = (l+m-1.0)/((double)(l-m));
                double Plm = fact1*x*P[l1][m]-fact2*P[l2][m];
                P[l].insert(pair<int,double>(m,Plm));
            }
        }
        //P[l][-m]
        for(l=1; l<=n; l++){
            for(m=1; m<=l; m++){
                int fact1 = 1;
                if(m % 2) fact1 = -1;
                double Plmn = (double)fact1*Factorial(l-m)/Factorial(l+m)*P[l][m];
                P[l].insert(pair<int,double>(-m,Plmn));
            }
        }
        
	return P;
}

//////////////////////////////////////////////////////////////////////
// Class AssociatedLegendre<1>
//////////////////////////////////////////////////////////////////////

AssociatedLegendre<1>::AssociatedLegendre() {}

AssociatedLegendre<1>::~AssociatedLegendre() {}

AssociatedLegendre<1>::AssociatedLegendre(const unsigned &n_) : n(n_) {}

vector<map<int,double> >
AssociatedLegendre<1>::operator () (const double &x) const
{
//	return P'[l][m], l=0,...,n, m=-l,...,l

	vector<map<int,double> > DP;

        Element::Interpolation::Utils::AssociatedLegendre<0> LP(n);
        vector<map<int,double> > P = LP(x);

        //DP[0][0]
        map<int,double> DP0map;        
        DP0map.insert(pair<int,double>(0,0));
        DP.push_back(DP0map);
        if (n == 0) return DP;
        
        unsigned int l, m;
        //DP[l][l]
        for(l=1; l<=n; l++){
            int fact1 = -1;
            if(l % 2) fact1 = 1;
            map<int,double> DPLmap;
            double DPll = (double)fact1*DoubleFactorial(2*l-1)*l*x*pow(1-x*x,0.5*l-1);
            DPLmap.insert(pair<int,double>(l,DPll));
            DP.push_back(DPLmap);
        }
        //DP[l+1][l]
        for(l=0; l<n; l++){
            DP[l+1].insert(pair<int,double>(l,(2*l+1)*P[l][l]+(2*l+1)*x*DP[l][l]));
        }
        //DP[l][m]
        for(m=0; m<=n-2 && n>1; m++){
            for(l=m+2; l<=n; l++){
                unsigned int l1 = l-1;
                unsigned int l2 = l-2;
                double fact1 = (2.0*l-1.0)/((double)(l-m));
                double fact2 = (l+m-1.0)/((double)(l-m));
                double DPlm = fact1*P[l1][m]+fact1*x*DP[l1][m]-fact2*DP[l2][m];
                DP[l].insert(pair<int,double>(m,DPlm));
            }
        }
        //DP[l][-m]
        for(l=1; l<=n; l++){
            for(m=1; m<=l; m++){
                int fact1 = 1;
                if(m % 2) fact1 = -1;
                double DPlmn = (double)fact1*Factorial(l-m)/Factorial(l+m)*DP[l][m];
                DP[l].insert(pair<int,double>(-m,DPlmn));
            }
        }
       
	return DP;
}

}

}

}



//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

unsigned int DoubleFactorial(const unsigned int &n)
{
    if(n == 0 || n ==1) return 1;
    unsigned int m=n*DoubleFactorial(n-2);
    return m;
}

unsigned int Factorial(const unsigned int &n)
{
    if(n == 0 || n ==1) return 1;
    unsigned int m=n*Factorial(n-1);
    return m;
}

