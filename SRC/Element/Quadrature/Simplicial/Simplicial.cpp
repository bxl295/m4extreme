// Simplicial.cpp: implementation of the Rule class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Simplicial.h"
#include "../Gaussian/Gaussian.h"

namespace Element
{
namespace Quadrature
{
namespace Simplicial
{
//////////////////////////////////////////////////////////////////////
// Class Rule
//////////////////////////////////////////////////////////////////////

Rule::Rule() {}

Rule::~Rule() {}

Rule::Rule(
	const unsigned int &rhs_n, 
	const unsigned int &rhs_N)
	: n(rhs_n), N(rhs_N) {}

vector<Set::VectorSpace::Vector> 
Rule::GetQ() const
{
	assert((n == 1) || (n == 2) || (n == 3));

	if (n == 1)
	{
		Element::Quadrature::Gaussian::Rule G(1,N);
		vector<Set::VectorSpace::Vector> Q = G.GetQ();
		for (unsigned int q=0; q<Q.size(); q++)
			Q[q][0] = 0.5*(1.0 + Q[q][0]);
		return Q;
	}

	if (n == 2)
	{
		switch (N)
		{
			case 0: return TriQ1();
			case 1: return TriQ1();
			case 2: return TriQ2();
			case 3: return TriQ3();
			case 4: return TriQ4();
			case 5: return TriQ5();
			case 6: return TriQ6();
			case 7: return TriQ7();
			case 8: return TriQ8();
			case 9: return TriQ9();
			default: throw(0);
		}
	}

	if (n == 3)
	{
		switch (N)
		{
			case 0: return TetQ1();
			case 1: return TetQ1();
			case 2: return TetQ2();
			case 3: return TetQ3();
			case 4: return TetQ4();
			case 5: return TetQ5();
			case 6: return TetQ6();
			default: throw(0);
		}
	}

	throw(0);
}

vector<double> 
Rule::GetW() const
{
	assert((n == 1) || (n == 2) || (n == 3));

	if (n == 1)
	{
		Element::Quadrature::Gaussian::Rule G(1,N);
		vector<double> W = G.GetW();
		for (unsigned int q=0; q<W.size(); q++) W[q] *= 0.5;
		return W;
	}

	if (n == 2)
	{
		switch (N)
		{
			case 0: return TriW1();
			case 1: return TriW1();
			case 2: return TriW2();
			case 3: return TriW3();
			case 4: return TriW4();
			case 5: return TriW5();
			case 6: return TriW6();
			case 7: return TriW7();
			case 8: return TriW8();
			case 9: return TriW9();
			default: throw(0);
		}
	}

	if (n == 3)
	{
		switch (N)
		{
			case 0: return TetW1();
			case 1: return TetW1();
			case 2: return TetW2();
			case 3: return TetW3();
			case 4: return TetW4();
			case 5: return TetW5();
			case 6: return TetW6();
			default: throw(0);
		}
	}

	throw(0);
}

//////////////////////////////////////////////////////////////////////
// Triangle
//////////////////////////////////////////////////////////////////////

const vector<Set::VectorSpace::Vector>
Rule::TriQ1() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(1,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(1,T0);

	const double OneThird=1.e0/3.e0;

	T[0][0] = OneThird;
	T[0][1] = OneThird;
	T[0][2] = OneThird;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

  // second set of 3-point formula
// const vector<Set::VectorSpace::Vector>
// Rule::TriQ2() const
// {
// 	Set::VectorSpace::Vector Q0(2);
// 	vector<Set::VectorSpace::Vector> Q(3,Q0);
// 	Set::VectorSpace::Vector T0(3);
// 	vector<Set::VectorSpace::Vector> T(3,T0);

// 	const double OneSixth=1.e0/6.e0;

// 	T[0][0] = 0.0e0;
// 	T[0][1] = 0.5e0;
// 	T[0][2] = 0.5e0;

// 	T[1][0] = 0.5e0;
// 	T[1][1] = 0.0e0;
// 	T[1][2] = 0.5e0;

// 	T[2][0] = 0.5e0;
// 	T[2][1] = 0.5e0;
// 	T[2][2] = 0.0e0;

// 	for (unsigned int k=0; k<Q.size(); k++) 
// 	{
// 		Q[k][0] = T[k][0]; 
// 		Q[k][1] = T[k][1]; 
// 	}

// 	return Q;
// }

const vector<Set::VectorSpace::Vector>
Rule::TriQ2() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(3,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(3,T0);

	const double OneSixth=1.e0/6.e0;
	const double TwoThird=2.0e0/3.0e0;

	T[0][0] = OneSixth;
	T[0][1] = OneSixth;
	T[0][2] = TwoThird;

	T[1][0] = TwoThird;
	T[1][1] = OneSixth;
	T[1][2] = OneSixth;

	T[2][0] = OneSixth;
	T[2][1] = TwoThird;
	T[2][2] = OneSixth;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TriQ3() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(6,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(6,T0);

	const double OneSixth=1.e0/6.e0;
	const double TwoThirds=2.e0/3.e0;
	const double OneSixtieth=1.e0/60.e0;

	T[0][0] = 0.0e0;
	T[0][1] = 0.5e0;
	T[0][2] = 0.5e0;

	T[1][0] = 0.5e0;
	T[1][1] = 0.0e0;
	T[1][2] = 0.5e0;

	T[2][0] = 0.5e0;
	T[2][1] = 0.5e0;
	T[2][2] = 0.0e0;

	T[3][0] = TwoThirds;
	T[3][1] = OneSixth;
	T[3][2] = OneSixth;

	T[4][0] = OneSixth;
	T[4][1] = TwoThirds;
	T[4][2] = OneSixth;

	T[5][0] = OneSixth;
	T[5][1] = OneSixth;
	T[5][2] = TwoThirds;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TriQ4() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(6,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(6,T0);

	T[0][1] = 0.0915762135097707434595714634022015e0;
	T[0][2] = 0.0915762135097707434595714634022015e0;
	T[0][0] = 1.e0 - T[0][1] - T[0][2];

	T[1][0] = 0.0915762135097707434595714634022015e0;
	T[1][2] = 0.0915762135097707434595714634022015e0;
	T[1][1] = 1.e0 - T[1][0] - T[1][2];

	T[2][0] = 0.0915762135097707434595714634022015e0;
	T[2][1] = 0.0915762135097707434595714634022015e0;
	T[2][2] = 1.e0 - T[2][0] - T[2][1];

	T[3][1] = 0.445948490915964886318329253883051e0;
	T[3][2] = 0.445948490915964886318329253883051e0;
	T[3][0] = 1.e0 - T[3][1] - T[3][2];

	T[4][0] = 0.445948490915964886318329253883051e0;
	T[4][2] = 0.445948490915964886318329253883051e0;
	T[4][1] = 1.e0 - T[4][0] - T[4][2];

	T[5][0] = 0.445948490915964886318329253883051e0;
	T[5][1] = 0.445948490915964886318329253883051e0;
	T[5][2] = 1.e0 - T[5][0] - T[5][1];

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TriQ5() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(7,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(7,T0);

	const double OneThird=1.e0/3.e0;

	T[0][0] = OneThird;
	T[0][1] = OneThird;
	T[0][2] = OneThird;

	T[1][1] = 0.101286507323456338800987361915123e0;
	T[1][2] = 0.101286507323456338800987361915123e0;
	T[1][0] = 1.e0 - T[1][1] - T[1][2];

	T[2][0] = 0.101286507323456338800987361915123e0;
	T[2][2] = 0.101286507323456338800987361915123e0;
	T[2][1] = 1.e0 - T[2][0] - T[2][2];

	T[3][0] = 0.101286507323456338800987361915123e0;
	T[3][1] = 0.101286507323456338800987361915123e0;
	T[3][2] = 1.e0 - T[3][0] - T[3][1];

	T[4][1] = 0.470142064105115089770441209513447e0;
	T[4][2] = 0.470142064105115089770441209513447e0;
	T[4][0] = 1.e0 - T[4][1] - T[4][2];

	T[5][0] = 0.470142064105115089770441209513447e0;
	T[5][2] = 0.470142064105115089770441209513447e0;
	T[5][1] = 1.e0 - T[5][0] - T[5][2];

	T[6][0] = 0.470142064105115089770441209513447e0;
	T[6][1] = 0.470142064105115089770441209513447e0;
	T[6][2] = 1.e0 - T[6][0] - T[6][1];

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TriQ6() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(12,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(12,T0);

    T[0][1] = 0.0630890144915022283403316028708191e0;
    T[0][2] = 0.0630890144915022283403316028708191e0;
    T[0][0] = 1.e0 - T[0][1] - T[0][2];

    T[1][0] = 0.0630890144915022283403316028708191e0;
    T[1][2] = 0.0630890144915022283403316028708191e0;
    T[1][1] = 1.e0 - T[1][0] - T[1][2];

    T[2][0] = 0.0630890144915022283403316028708191e0;
    T[2][1] = 0.0630890144915022283403316028708191e0;
    T[2][2] = 1.e0 - T[2][0] - T[2][1];

    T[3][1] = 0.249286745170910421291638553107019e0;
    T[3][2] = 0.249286745170910421291638553107019e0;
    T[3][0] = 1.e0 - T[3][1] - T[3][2];

    T[4][0] = 0.249286745170910421291638553107019e0;
    T[4][2] = 0.249286745170910421291638553107019e0;
    T[4][1] = 1.e0 - T[4][0] - T[4][2];

    T[5][0] = 0.249286745170910421291638553107019e0;
    T[5][1] = 0.249286745170910421291638553107019e0;
    T[5][2] = 1.e0 - T[5][0] - T[5][1];

    T[6][1] = 0.0531450498448169473532496716313981e0;
    T[6][2] = 0.310352451033784405416607733956552e0;
    T[6][0] = 1.e0 - T[6][1] - T[6][2];

    T[7][0] = 0.0531450498448169473532496716313981e0;
    T[7][2] = 0.310352451033784405416607733956552e0;
    T[7][1] = 1.e0 - T[7][0] - T[7][2];

    T[8][0] = 0.0531450498448169473532496716313981e0;
    T[8][1] = 0.310352451033784405416607733956552e0;
    T[8][2] = 1.e0 - T[8][0] - T[8][1];

    T[9][1] = 0.310352451033784405416607733956552e0;
    T[9][2] = 0.0531450498448169473532496716313981e0;
    T[9][0] = 1.e0 - T[9][1] - T[9][2];

    T[10][0] = 0.310352451033784405416607733956552e0;
    T[10][2] = 0.0531450498448169473532496716313981e0;
    T[10][1] = 1.e0 - T[10][0] - T[10][2];

    T[11][0] = 0.310352451033784405416607733956552e0;
    T[11][1] = 0.0531450498448169473532496716313981e0;
    T[11][2] = 1.e0 - T[11][0] - T[11][1];

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TriQ7() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(13,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(13,T0);

    const double OneThird=1.e0/3.e0;

    T[0][0] = OneThird;
    T[0][1] = OneThird;
    T[0][2] = OneThird;

    T[1][1] = 0.260345966079039826926242469139236e0;
    T[1][2] = 0.260345966079039826926242469139236e0;
    T[1][0] = 1.e0 - T[1][1] - T[1][2];

    T[2][0] = 0.260345966079039826926242469139236e0;
    T[2][2] = 0.260345966079039826926242469139236e0;
    T[2][1] = 1.e0 - T[2][0] - T[2][2];

    T[3][0] = 0.260345966079039826926242469139236e0;
    T[3][1] = 0.260345966079039826926242469139236e0;
    T[3][2] = 1.e0 - T[3][0] - T[3][1];

    T[4][1] = 0.0651301029022158115380259063119754e0;
    T[4][2] = 0.0651301029022158115380259063119754e0;
    T[4][0] = 1.e0 - T[4][1] - T[4][2];

    T[5][0] = 0.0651301029022158115380259063119754e0;
    T[5][2] = 0.0651301029022158115380259063119754e0;
    T[5][1] = 1.e0 - T[5][0] - T[5][2];

    T[6][0] = 0.0651301029022158115380259063119754e0;
    T[6][1] = 0.0651301029022158115380259063119754e0;
    T[6][2] = 1.e0 - T[6][0] - T[6][1];

    T[7][1] = 0.0486903154253164117930215585284131e0;
    T[7][2] = 0.312865496004873861406644476768401e0;
    T[7][0] = 1.e0 - T[7][1] - T[7][2];

    T[8][0] = 0.0486903154253164117930215585284131e0;
    T[8][2] = 0.312865496004873861406644476768401e0;
    T[8][1] = 1.e0 - T[8][0] - T[8][2];

    T[9][0] = 0.0486903154253164117930215585284131e0;
    T[9][1] = 0.312865496004873861406644476768401e0;
    T[9][2] = 1.e0 - T[9][0] - T[9][1];

    T[10][1] = 0.312865496004873861406644476768401e0;
    T[10][2] = 0.0486903154253164117930215585284131e0;
    T[10][0] = 1.e0 - T[10][1] - T[10][2];

    T[11][0] = 0.312865496004873861406644476768401e0;
    T[11][2] = 0.0486903154253164117930215585284131e0;
    T[11][1] = 1.e0 - T[11][0] - T[11][2];

    T[12][0] = 0.312865496004873861406644476768401e0;
    T[12][1] = 0.0486903154253164117930215585284131e0;
    T[12][2] = 1.e0 - T[12][0] - T[12][1];

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TriQ8() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(16,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(16,T0);

    const double OneThird=1.e0/3.e0;

    T[0][0] = OneThird;
    T[0][1] = OneThird;
    T[0][2] = OneThird;

    T[1][1] = 0.476665439382152376067071748878465e0;
    T[1][2] = 0.476665439382152376067071748878465e0;
    T[1][0] = 1.e0 - T[1][1] - T[1][2];

    T[2][0] = 0.476665439382152376067071748878465e0;
    T[2][2] = 0.476665439382152376067071748878465e0;
    T[2][1] = 1.e0 - T[2][0] - T[2][2];

    T[3][0] = 0.476665439382152376067071748878465e0;
    T[3][1] = 0.476665439382152376067071748878465e0;
    T[3][2] = 1.e0 - T[3][0] - T[3][1];

    T[4][1] = 0.0337718440544803428071248668994412e0;
    T[4][2] = 0.0337718440544803428071248668994412e0;
    T[4][0] = 1.e0 - T[4][1] - T[4][2];

    T[5][0] = 0.0337718440544803428071248668994412e0;
    T[5][2] = 0.0337718440544803428071248668994412e0;
    T[5][1] = 1.e0 - T[5][0] - T[5][2];

    T[6][0] = 0.0337718440544803428071248668994412e0;
    T[6][1] = 0.0337718440544803428071248668994412e0;
    T[6][2] = 1.e0 - T[6][0] - T[6][1];

    T[7][1] = 0.270347889165403620102711192274472e0;
    T[7][2] = 0.270347889165403620102711192274472e0;
    T[7][0] = 1.e0 - T[7][1] - T[7][2];

    T[8][0] = 0.270347889165403620102711192274472e0;
    T[8][2] = 0.270347889165403620102711192274472e0;
    T[8][1] = 1.e0 - T[8][0] - T[8][2];

    T[9][0] = 0.270347889165403620102711192274472e0;
    T[9][1] = 0.270347889165403620102711192274472e0;
    T[9][2] = 1.e0 - T[9][0] - T[9][1];

    T[10][1] = 0.745829490767251371671683009263721e0;
    T[10][2] = 0.0514643354866615258021682284492792e0;
    T[10][0] = 1.e0 - T[10][1] - T[10][2];

    T[11][0] = 0.745829490767251371671683009263721e0;
    T[11][2] = 0.0514643354866615258021682284492792e0;
    T[11][1] = 1.e0 - T[11][0] - T[11][2];

    T[12][0] = 0.745829490767251371671683009263721e0;
    T[12][1] = 0.0514643354866615258021682284492792e0;
    T[12][2] = 1.e0 - T[12][0] - T[12][1];

    T[13][1] = 0.0514643354866615258021682284492792e0;
    T[13][2] = 0.745829490767251371671683009263721e0;
    T[13][0] = 1.e0 - T[13][1] - T[13][2];

    T[14][0] = 0.0514643354866615258021682284492792e0;
    T[14][2] = 0.745829490767251371671683009263721e0;
    T[14][1] = 1.e0 - T[14][0] - T[14][2];

    T[15][0] = 0.0514643354866615258021682284492792e0;
    T[15][1] = 0.745829490767251371671683009263721e0;
    T[15][2] = 1.e0 - T[15][0] - T[15][1];

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TriQ9() const
{
	Set::VectorSpace::Vector Q0(2);
	vector<Set::VectorSpace::Vector> Q(19,Q0);
	Set::VectorSpace::Vector T0(3);
	vector<Set::VectorSpace::Vector> T(19,T0);

    const double OneThird=1.e0/3.e0;

    T[0][0] = OneThird;
    T[0][1] = OneThird;
    T[0][2] = OneThird;

    T[1][1] = 0.489682519198737627783706924836192e0;
    T[1][2] = 0.489682519198737627783706924836192e0;
    T[1][0] = 1.e0 - T[1][1] - T[1][2];

    T[2][0] = 0.489682519198737627783706924836192e0;
    T[2][2] = 0.489682519198737627783706924836192e0;
    T[2][1] = 1.e0 - T[2][0] - T[2][2];

    T[3][0] = 0.489682519198737627783706924836192e0;
    T[3][1] = 0.489682519198737627783706924836192e0;
    T[3][2] = 1.e0 - T[3][0] - T[3][1];

    T[4][1] = 0.437089591492936637269930364435354e0;
    T[4][2] = 0.437089591492936637269930364435354e0;
    T[4][0] = 1.e0 - T[4][1] - T[4][2];

    T[5][0] = 0.437089591492936637269930364435354e0;
    T[5][2] = 0.437089591492936637269930364435354e0;
    T[5][1] = 1.e0 - T[5][0] - T[5][2];

    T[6][0] = 0.437089591492936637269930364435354e0;
    T[6][1] = 0.437089591492936637269930364435354e0;
    T[6][2] = 1.e0 - T[6][0] - T[6][1];

    T[7][1] = 0.188203535619032730240961280467335e0;
    T[7][2] = 0.188203535619032730240961280467335e0;
    T[7][0] = 1.e0 - T[7][1] - T[7][2];

    T[8][0] = 0.188203535619032730240961280467335e0;
    T[8][2] = 0.188203535619032730240961280467335e0;
    T[8][1] = 1.e0 - T[8][0] - T[8][2];

    T[9][0] = 0.188203535619032730240961280467335e0;
    T[9][1] = 0.188203535619032730240961280467335e0;
    T[9][2] = 1.e0 - T[9][0] - T[9][1];

    T[10][1] = 0.0447295133944527098651065899662763e0;
    T[10][2] = 0.0447295133944527098651065899662763e0;
    T[10][0] = 1.e0 - T[10][1] - T[10][2];

    T[11][0] = 0.0447295133944527098651065899662763e0;
    T[11][2] = 0.0447295133944527098651065899662763e0;
    T[11][1] = 1.e0 - T[11][0] - T[11][2];

    T[12][0] = 0.0447295133944527098651065899662763e0;
    T[12][1] = 0.0447295133944527098651065899662763e0;
    T[12][2] = 1.e0 - T[12][0] - T[12][1];

    T[13][1] = 0.741198598784498020690079873523423e0;
    T[13][2] = 0.0368384120547362836348175987833851e0;
    T[13][0] = 1.e0 - T[13][1] - T[13][2];

    T[14][0] = 0.741198598784498020690079873523423e0;
    T[14][2] = 0.0368384120547362836348175987833851e0;
    T[14][1] = 1.e0 - T[14][0] - T[14][2];

    T[15][0] = 0.741198598784498020690079873523423e0;
    T[15][1] = 0.0368384120547362836348175987833851e0;
    T[15][2] = 1.e0 - T[15][0] - T[15][1];

    T[16][1] = 0.0368384120547362836348175987833851e0;
    T[16][2] = 0.741198598784498020690079873523423e0;
    T[16][0] = 1.e0 - T[16][1] - T[16][2];

    T[17][0] = 0.0368384120547362836348175987833851e0;
    T[17][2] = 0.741198598784498020690079873523423e0;
    T[17][1] = 1.e0 - T[17][0] - T[17][2];

    T[18][0] = 0.0368384120547362836348175987833851e0;
    T[18][1] = 0.741198598784498020690079873523423e0;
    T[18][2] = 1.e0 - T[18][0] - T[18][1];

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
	}

	return Q;
}

const vector<double>
Rule::TriW1() const
{
	vector<double> W(1);

	W[0]    = 0.5e0;

	return W;
}

const vector<double>
Rule::TriW2() const
{
	vector<double> W(3);

	const double OneSixth=1.e0/6.e0;

	W[0]    = OneSixth;
	W[1]    = OneSixth;
	W[2]    = OneSixth;

	return W;
}

const vector<double>
Rule::TriW3() const
{
	vector<double> W(6);

	const double OneSixth=1.e0/6.e0;
	const double TwoThirds=2.e0/3.e0;
	const double OneSixtieth=1.e0/60.e0;

	W[0]    = OneSixtieth;
	W[1]    = OneSixtieth;
	W[2]    = OneSixtieth;
	W[3]    = 0.15e0;
	W[4]    = 0.15e0;
	W[5]    = 0.15e0;

	return W;
}

const vector<double>
Rule::TriW4() const
{
	vector<double> W(6);

	W[0]    = 0.0549758718276609338191631624501052e0;
	W[1]    = 0.0549758718276609338191631624501052e0;
	W[2]    = 0.0549758718276609338191631624501052e0;
	W[3]    = 0.111690794839005732847503504216561e0;
	W[4]    = 0.111690794839005732847503504216561e0;
	W[5]    = 0.111690794839005732847503504216561e0;

	return W;
}

const vector<double>
Rule::TriW5() const
{
	vector<double> W(7);

	const double OneThird=1.e0/3.e0;

	W[0]    = 0.1125e0;
	W[1]    = 0.0629695902724135762978419727500906e0;
	W[2]    = 0.0629695902724135762978419727500906e0;
	W[3]    = 0.0629695902724135762978419727500906e0;
	W[4]    = 0.0661970763942530903688246939165759e0;
	W[5]    = 0.0661970763942530903688246939165759e0;
	W[6]    = 0.0661970763942530903688246939165759e0;

	return W;
}

const vector<double>
Rule::TriW6() const
{
	vector<double> W(12);

    W[0]    = 0.0254224531851034084604684045534344e0;
    W[1]    = 0.0254224531851034084604684045534344e0;
    W[2]    = 0.0254224531851034084604684045534344e0;
    W[3]    = 0.0583931378631896830126448056927897e0;
    W[4]    = 0.0583931378631896830126448056927897e0;
    W[5]    = 0.0583931378631896830126448056927897e0;
    W[6]    = 0.0414255378091867875967767282102212e0;
    W[7]    = 0.0414255378091867875967767282102212e0;
    W[8]    = 0.0414255378091867875967767282102212e0;
    W[9]    = 0.0414255378091867875967767282102212e0;
    W[10]   = 0.0414255378091867875967767282102212e0;
    W[11]   = 0.0414255378091867875967767282102212e0;

	return W;
}

const vector<double>
Rule::TriW7() const
{
	vector<double> W(13);

    const double OneThird=1.e0/3.e0;

    W[0]    = -0.0747850222338408753148556277365486e0;
    W[1]    = 0.0878076287166039058767597055781635e0;
    W[2]    = 0.0878076287166039058767597055781635e0;
    W[3]    = 0.0878076287166039058767597055781635e0;
    W[4]    = 0.0266736178044192456349936444495335e0;
    W[5]    = 0.0266736178044192456349936444495335e0;
    W[6]    = 0.0266736178044192456349936444495335e0;
    W[7]    = 0.0385568804451285701299325962755762e0;
    W[8]    = 0.0385568804451285701299325962755762e0;
    W[9]    = 0.0385568804451285701299325962755762e0;
    W[10]   = 0.0385568804451285701299325962755762e0;
    W[11]   = 0.0385568804451285701299325962755762e0;
    W[12]   = 0.0385568804451285701299325962755762e0;

	return W;
}

const vector<double>
Rule::TriW8() const
{
	vector<double> W(16);

    const double OneThird=1.e0/3.e0;

    W[0]    = -0.141709192555694106139608305269904e0;
    W[1]    = 0.0349534809663262894329607391832926e0;
    W[2]    = 0.0349534809663262894329607391832926e0;
    W[3]    = 0.0349534809663262894329607391832926e0;
    W[4]    = 0.854545633580044472243722969322347e-2;
    W[5]    = 0.854545633580044472243722969322347e-2;;
    W[6]    = 0.854545633580044472243722969322347e-2;;
    W[7]    = 0.109414941165223740315377978302527e0;
    W[8]    = 0.109414941165223740315377978302527e0;
    W[9]    = 0.109414941165223740315377978302527e0;
    W[10]   = 0.0304945928589404471212134106221289e0;
    W[11]   = 0.0304945928589404471212134106221289e0;
    W[12]   = 0.0304945928589404471212134106221289e0;
    W[13]   = 0.0304945928589404471212134106221289e0;
    W[14]   = 0.0304945928589404471212134106221289e0;
    W[15]   = 0.0304945928589404471212134106221289e0;

	return W;
}

const vector<double>
Rule::TriW9() const
{
	vector<double> W(19);

    const double OneThird=1.e0/3.e0;

    W[0]    = 0.0485678981413994169096209912536443e0;
    W[1]    = 0.0156673501135695352684274156436046e0;
    W[2]    = 0.0156673501135695352684274156436046e0;
    W[3]    = 0.0156673501135695352684274156436046e0;
    W[4]    = 0.0389137705023871396583696781497019e0;
    W[5]    = 0.0389137705023871396583696781497019e0;
    W[6]    = 0.0389137705023871396583696781497019e0;
    W[7]    = 0.0398238694636051265164458871320226e0;
    W[8]    = 0.0398238694636051265164458871320226e0;
    W[9]    = 0.0398238694636051265164458871320226e0;
    W[10]   = 0.0127888378293490156308393992794999e0;
    W[11]   = 0.0127888378293490156308393992794999e0;
    W[12]   = 0.0127888378293490156308393992794999e0;
    W[13]   = 0.0216417696886446886446886446886446e0;
    W[14]   = 0.0216417696886446886446886446886446e0;
    W[15]   = 0.0216417696886446886446886446886446e0;
    W[16]   = 0.0216417696886446886446886446886446e0;
    W[17]   = 0.0216417696886446886446886446886446e0;
    W[18]   = 0.0216417696886446886446886446886446e0;

	return W;
}

//////////////////////////////////////////////////////////////////////
// Tetrahedron
//////////////////////////////////////////////////////////////////////

const vector<Set::VectorSpace::Vector>
Rule::TetQ1() const
{
	Set::VectorSpace::Vector Q0(3);
	vector<Set::VectorSpace::Vector> Q(1,Q0);
	Set::VectorSpace::Vector T0(4);
	vector<Set::VectorSpace::Vector> T(1,T0);

    T[0][0]  = 0.25e0;
    T[0][1]  = 0.25e0;
    T[0][2]  = 0.25e0;
    T[0][3]  = 0.25e0;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
		Q[k][2] = T[k][2]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TetQ2() const
{
	Set::VectorSpace::Vector Q0(3);
	vector<Set::VectorSpace::Vector> Q(4,Q0);
	Set::VectorSpace::Vector T0(4);
	vector<Set::VectorSpace::Vector> T(4,T0);

    T[0][0] = 0.5854101966249685e0;
    T[0][1] = 0.1381966011250150e0;
    T[0][2] = 0.1381966011250150e0;
    T[0][3] = 0.1381966011250150e0;

    T[1][0] = 0.1381966011250150e0;
    T[1][1] = 0.5854101966249685e0;
    T[1][2] = 0.1381966011250150e0;
    T[1][3] = 0.1381966011250150e0;

    T[2][0] = 0.1381966011250150e0;
    T[2][1] = 0.1381966011250150e0;
    T[2][2] = 0.5854101966249685e0;
    T[2][3] = 0.1381966011250150e0;

    T[3][0] = 0.1381966011250150e0;
    T[3][1] = 0.1381966011250150e0;
    T[3][2] = 0.1381966011250150e0;
    T[3][3] = 0.5854101966249685e0;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
		Q[k][2] = T[k][2]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TetQ3() const
{
	Set::VectorSpace::Vector Q0(3);
	vector<Set::VectorSpace::Vector> Q(5,Q0);
	Set::VectorSpace::Vector T0(4);
	vector<Set::VectorSpace::Vector> T(5,T0);

    const double OneSixth=1.0/6.0;

    T[0][0] = 0.25e0;
    T[0][1] = 0.25e0;
    T[0][2] = 0.25e0;
    T[0][3] = 0.25e0;

    T[1][0] = 0.5e0;
    T[1][1] = OneSixth;
    T[1][2] = OneSixth;
    T[1][3] = OneSixth;

    T[2][0] = OneSixth;
    T[2][1] = 0.5e0;
    T[2][2] = OneSixth;
    T[2][3] = OneSixth;

    T[3][0] = OneSixth;
    T[3][1] = OneSixth;
    T[3][2] = 0.5e0;
    T[3][3] = OneSixth;

    T[4][0] = OneSixth;
    T[4][1] = OneSixth;
    T[4][2] = OneSixth;
    T[4][3] = 0.5e0;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
		Q[k][2] = T[k][2]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TetQ4() const
{
	Set::VectorSpace::Vector Q0(3);
	vector<Set::VectorSpace::Vector> Q(16,Q0);
	Set::VectorSpace::Vector T0(4);
	vector<Set::VectorSpace::Vector> T(16,T0);

    T[0][0] = 0.7716429020672371e0;
    T[0][1] = 0.7611903264425430e-1;
    T[0][2] = 0.7611903264425430e-1;
    T[0][3] = 0.7611903264425430e-1;

    T[1][0] = 0.7611903264425430e-1;
    T[1][1] = 0.7716429020672371e0;
    T[1][2] = 0.7611903264425430e-1;
    T[1][3] = 0.7611903264425430e-1;

    T[2][0] = 0.7611903264425430e-1;
    T[2][1] = 0.7611903264425430e-1;
    T[2][2] = 0.7716429020672371e0;
    T[2][3] = 0.7611903264425430e-1;

    T[3][0] = 0.7611903264425430e-1;
    T[3][1] = 0.7611903264425430e-1;
    T[3][2] = 0.7611903264425430e-1;
    T[3][3] = 0.7716429020672371e0;

    T[4][0] = 0.1197005277978019e0;
    T[4][1] = 0.7183164526766925e-1;
    T[4][2] = 0.4042339134672644e0;
    T[4][3] = 0.4042339134672644e0;

    T[5][0] = 0.7183164526766925e-1;
    T[5][1] = 0.1197005277978019e0;
    T[5][2] = 0.4042339134672644e0;
    T[5][3] = 0.4042339134672644e0;

    T[6][0] = 0.1197005277978019e0;
    T[6][1] = 0.4042339134672644e0;
    T[6][2] = 0.7183164526766925e-1;
    T[6][3] = 0.4042339134672644e0;

    T[7][0] = 0.7183164526766925e-1;
    T[7][1] = 0.4042339134672644e0;
    T[7][2] = 0.1197005277978019e0;
    T[7][3] = 0.4042339134672644e0;

    T[8][0] = 0.4042339134672644e0;
    T[8][1] = 0.1197005277978019e0;
    T[8][2] = 0.7183164526766925e-1;
    T[8][3] = 0.4042339134672644e0;

    T[9][0] = 0.4042339134672644e0;
    T[9][1] = 0.7183164526766925e-1;
    T[9][2] = 0.1197005277978019e0;
    T[9][3] = 0.4042339134672644e0;

    T[10][0] = 0.1197005277978019e0;
    T[10][1] = 0.4042339134672644e0;
    T[10][2] = 0.4042339134672644e0;
    T[10][3] = 0.7183164526766925e-1;

    T[11][0] = 0.7183164526766925e-1;
    T[11][1] = 0.4042339134672644e0;
    T[11][2] = 0.4042339134672644e0;
    T[11][3] = 0.1197005277978019e0;

    T[12][0] = 0.4042339134672644e0;
    T[12][1] = 0.1197005277978019e0;
    T[12][2] = 0.4042339134672644e0;
    T[12][3] = 0.7183164526766925e-1;

    T[13][0] = 0.4042339134672644e0;
    T[13][1] = 0.7183164526766925e-1;
    T[13][2] = 0.4042339134672644e0;
    T[13][3] = 0.1197005277978019e0;

    T[14][0] = 0.4042339134672644e0;
    T[14][1] = 0.4042339134672644e0;
    T[14][2] = 0.1197005277978019e0;
    T[14][3] = 0.7183164526766925e-1;

    T[15][0] = 0.4042339134672644e0;
    T[15][1] = 0.4042339134672644e0;
    T[15][2] = 0.7183164526766925e-1;
    T[15][3] = 0.1197005277978019e0;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
		Q[k][2] = T[k][2]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TetQ5() const
{
	Set::VectorSpace::Vector Q0(3);
	vector<Set::VectorSpace::Vector> Q(17,Q0);
	Set::VectorSpace::Vector T0(4);
	vector<Set::VectorSpace::Vector> T(17,T0);

    T[0][0] = 0.25e0;
    T[0][1] = 0.25e0;
    T[0][2] = 0.25e0;
    T[0][3] = 0.25e0;

    T[1][0] = 0.7316369079576180e0;
    T[1][1] = 0.8945436401412733e-1;
    T[1][2] = 0.8945436401412733e-1;
    T[1][3] = 0.8945436401412733e-1;

    T[2][0] = 0.8945436401412733e-1;
    T[2][1] = 0.7316369079576180e0;
    T[2][2] = 0.8945436401412733e-1;
    T[2][3] = 0.8945436401412733e-1;

    T[3][0] = 0.8945436401412733e-1;
    T[3][1] = 0.8945436401412733e-1;
    T[3][2] = 0.7316369079576180e0;
    T[3][3] = 0.8945436401412733e-1;

    T[4][0] = 0.8945436401412733e-1;
    T[4][1] = 0.8945436401412733e-1;
    T[4][2] = 0.8945436401412733e-1;
    T[4][3] = 0.7316369079576180e0;

    T[5][0] = 0.1325810999384657e0;
    T[5][1] = 0.2454003792903000e-1;
    T[5][2] = 0.4214394310662522e0;
    T[5][3] = 0.4214394310662522e0;

    T[6][0] = 0.2454003792903000e-1;
    T[6][1] = 0.1325810999384657e0;
    T[6][2] = 0.4214394310662522e0;
    T[6][3] = 0.4214394310662522e0;

    T[7][0] = 0.1325810999384657e0;
    T[7][1] = 0.4214394310662522e0;
    T[7][2] = 0.2454003792903000e-1;
    T[7][3] = 0.4214394310662522e0;

    T[8][0] = 0.2454003792903000e-1;
    T[8][1] = 0.4214394310662522e0;
    T[8][2] = 0.1325810999384657e0;
    T[8][3] = 0.4214394310662522e0;

    T[9][0] = 0.4214394310662522e0;
    T[9][1] = 0.1325810999384657e0;
    T[9][2] = 0.2454003792903000e-1;
    T[9][3] = 0.4214394310662522e0;

    T[10][0] = 0.4214394310662522e0;
    T[10][1] = 0.2454003792903000e-1;
    T[10][2] = 0.1325810999384657e0;
    T[10][3] = 0.4214394310662522e0;

    T[11][0] = 0.1325810999384657e0;
    T[11][1] = 0.4214394310662522e0;
    T[11][2] = 0.4214394310662522e0;
    T[11][3] = 0.2454003792903000e-1;

    T[12][0] = 0.2454003792903000e-1;
    T[12][1] = 0.4214394310662522e0;
    T[12][2] = 0.4214394310662522e0;
    T[12][3] = 0.1325810999384657e0;

    T[13][0] = 0.4214394310662522e0;
    T[13][1] = 0.1325810999384657e0;
    T[13][2] = 0.4214394310662522e0;
    T[13][3] = 0.2454003792903000e-1;

    T[14][0] = 0.4214394310662522e0;
    T[14][1] = 0.2454003792903000e-1;
    T[14][2] = 0.4214394310662522e0;
    T[14][3] = 0.1325810999384657e0;

    T[15][0] = 0.4214394310662522e0;
    T[15][1] = 0.4214394310662522e0;
    T[15][2] = 0.1325810999384657e0;
    T[15][3] = 0.2454003792903000e-1;

    T[16][0] = 0.4214394310662522e0;
    T[16][1] = 0.4214394310662522e0;
    T[16][2] = 0.2454003792903000e-1;
    T[16][3] = 0.1325810999384657e0;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
		Q[k][2] = T[k][2]; 
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::TetQ6() const
{
	Set::VectorSpace::Vector Q0(3);
	vector<Set::VectorSpace::Vector> Q(29,Q0);
	Set::VectorSpace::Vector T0(4);
	vector<Set::VectorSpace::Vector> T(29,T0);

    T[0][0] = 0.25e0;
    T[0][1] = 0.25e0;
    T[0][2] = 0.25e0;
    T[0][3] = 0.25e0;

    T[1][0] = 0.8277192480479295e0;
    T[1][1] = 0.5742691731735683e-1;
    T[1][2] = 0.5742691731735683e-1;
    T[1][3] = 0.5742691731735683e-1;

    T[2][0] = 0.5742691731735683e-1;
    T[2][1] = 0.8277192480479295e0;
    T[2][2] = 0.5742691731735683e-1;
    T[2][3] = 0.5742691731735683e-1;

    T[3][0] = 0.5742691731735683e-1;
    T[3][1] = 0.5742691731735683e-1;
    T[3][2] = 0.8277192480479295e0;
    T[3][3] = 0.5742691731735683e-1;

    T[4][0] = 0.5742691731735683e-1;
    T[4][1] = 0.5742691731735683e-1;
    T[4][2] = 0.5742691731735683e-1;
    T[4][3] = 0.8277192480479295e0;

    T[5][0] = 0.5135188412556341e-1;
    T[5][1] = 0.4860510285706072e0;
    T[5][2] = 0.2312985436519147e0;
    T[5][3] = 0.2312985436519147e0;

    T[6][0] = 0.4860510285706072e0;
    T[6][1] = 0.5135188412556341e-1;
    T[6][2] = 0.2312985436519147e0;
    T[6][3] = 0.2312985436519147e0;

    T[7][0] = 0.5135188412556341e-1;
    T[7][1] = 0.2312985436519147e0;
    T[7][2] = 0.4860510285706072e0;
    T[7][3] = 0.2312985436519147e0;

    T[8][0] = 0.4860510285706072e0;
    T[8][1] = 0.2312985436519147e0;
    T[8][2] = 0.5135188412556341e-1;
    T[8][3] = 0.2312985436519147e0;

    T[9][0] = 0.2312985436519147e0;
    T[9][1] = 0.5135188412556341e-1;
    T[9][2] = 0.4860510285706072e0;
    T[9][3] = 0.2312985436519147e0;

    T[10][0] = 0.2312985436519147e0;
    T[10][1] = 0.4860510285706072e0;
    T[10][2] = 0.5135188412556341e-1;
    T[10][3] = 0.2312985436519147e0;

    T[11][0] = 0.5135188412556341e-1;
    T[11][1] = 0.2312985436519147e0;
    T[11][2] = 0.2312985436519147e0;
    T[11][3] = 0.4860510285706072e0;

    T[12][0] = 0.4860510285706072e0;
    T[12][1] = 0.2312985436519147e0;
    T[12][2] = 0.2312985436519147e0;
    T[12][3] = 0.5135188412556341e-1;

    T[13][0] = 0.2312985436519147e0;
    T[13][1] = 0.5135188412556341e-1;
    T[13][2] = 0.2312985436519147e0;
    T[13][3] = 0.4860510285706072e0;

    T[14][0] = 0.2312985436519147e0;
    T[14][1] = 0.4860510285706072e0;
    T[14][2] = 0.2312985436519147e0;
    T[14][3] = 0.5135188412556341e-1;

    T[15][0] = 0.2312985436519147e0;
    T[15][1] = 0.2312985436519147e0;
    T[15][2] = 0.5135188412556341e-1;
    T[15][3] = 0.4860510285706072e0;

    T[16][0] = 0.2312985436519147e0;
    T[16][1] = 0.2312985436519147e0;
    T[16][2] = 0.4860510285706072e0;
    T[16][3] = 0.5135188412556341e-1;

    T[17][0] = 0.2967538129690260e0;
    T[17][1] = 0.6081079894015281e0;
    T[17][2] = 0.4756909881472290e-1;
    T[17][3] = 0.4756909881472290e-1;

    T[18][0] = 0.6081079894015281e0;
    T[18][1] = 0.2967538129690260e0;
    T[18][2] = 0.4756909881472290e-1;
    T[18][3] = 0.4756909881472290e-1;

    T[19][0] = 0.2967538129690260e0;
    T[19][1] = 0.4756909881472290e-1;
    T[19][2] = 0.6081079894015281e0;
    T[19][3] = 0.4756909881472290e-1;

    T[20][0] = 0.6081079894015281e0;
    T[20][1] = 0.4756909881472290e-1;
    T[20][2] = 0.2967538129690260e0;
    T[20][3] = 0.4756909881472290e-1;

    T[21][0] = 0.4756909881472290e-1;
    T[21][1] = 0.2967538129690260e0;
    T[21][2] = 0.6081079894015281e0;
    T[21][3] = 0.4756909881472290e-1;

    T[22][0] = 0.4756909881472290e-1;
    T[22][1] = 0.6081079894015281e0;
    T[22][2] = 0.2967538129690260e0;
    T[22][3] = 0.4756909881472290e-1;

    T[23][0] = 0.2967538129690260e0;
    T[23][1] = 0.4756909881472290e-1;
    T[23][2] = 0.4756909881472290e-1;
    T[23][3] = 0.6081079894015281e0;

    T[24][0] = 0.6081079894015281e0;
    T[24][1] = 0.4756909881472290e-1;
    T[24][2] = 0.4756909881472290e-1;
    T[24][3] = 0.2967538129690260e0;

    T[25][0] = 0.4756909881472290e-1;
    T[25][1] = 0.2967538129690260e0;
    T[25][2] = 0.4756909881472290e-1;
    T[25][3] = 0.6081079894015281e0;

    T[26][0] = 0.4756909881472290e-1;
    T[26][1] = 0.6081079894015281e0;
    T[26][2] = 0.4756909881472290e-1;
    T[26][3] = 0.2967538129690260e0;

    T[27][0] = 0.4756909881472290e-1;
    T[27][1] = 0.4756909881472290e-1;
    T[27][2] = 0.2967538129690260e0;
    T[27][3] = 0.6081079894015281e0;

    T[28][0] = 0.4756909881472290e-1;
    T[28][1] = 0.4756909881472290e-1;
    T[28][2] = 0.6081079894015281e0;
    T[28][3] = 0.2967538129690260e0;

	for (unsigned int k=0; k<Q.size(); k++) 
	{
		Q[k][0] = T[k][0]; 
		Q[k][1] = T[k][1]; 
		Q[k][2] = T[k][2]; 
	}

	return Q;
}

const vector<double>
Rule::TetW1() const
{
	vector<double> W(1);

    W[0]     = 1.0e0;

	return W;
}

const vector<double>
Rule::TetW2() const
{
	vector<double> W(4);

    W[0]    = 0.25e0;
    W[1]    = 0.25e0;
    W[2]    = 0.25e0;
    W[3]    = 0.25e0;

	return W;
}

const vector<double>
Rule::TetW3() const
{
	vector<double> W(5);

    const double OneSixth=1.0/6.0;

    W[0]    = -0.8e0;
    W[1]    = 0.45e0;
    W[2]    = 0.45e0;
    W[3]    = 0.45e0;
    W[4]    = 0.45e0;

	return W;
}

const vector<double>
Rule::TetW4() const
{
	vector<double> W(16);

    W[0]    = 0.5037379410012282e-1;
    W[1]    = 0.5037379410012282e-1;
    W[2]    = 0.5037379410012282e-1;
    W[3]    = 0.5037379410012282e-1;
    W[4]    = 0.6654206863329239e-1;
    W[5]    = 0.6654206863329239e-1;
    W[6]    = 0.6654206863329239e-1;
    W[7]    = 0.6654206863329239e-1;
    W[8]    = 0.6654206863329239e-1;
    W[9]    = 0.6654206863329239e-1;
    W[10]   = 0.6654206863329239e-1;
    W[11]   = 0.6654206863329239e-1;
    W[12]   = 0.6654206863329239e-1;
    W[13]   = 0.6654206863329239e-1;
    W[14]   = 0.6654206863329239e-1;
    W[15]   = 0.6654206863329239e-1;

	return W;
}

const vector<double>
Rule::TetW5() const
{
	vector<double> W(17);

    W[0]    = 0.1884185567365411e0;
    W[1]    = 0.6703858372604275e-1;
    W[2]    = 0.6703858372604275e-1;
    W[3]    = 0.6703858372604275e-1;
    W[4]    = 0.6703858372604275e-1;
    W[5]    = 0.4528559236327399e-1;
    W[6]    = 0.4528559236327399e-1;
    W[7]    = 0.4528559236327399e-1;
    W[8]    = 0.4528559236327399e-1;
    W[9]    = 0.4528559236327399e-1;
    W[10]   = 0.4528559236327399e-1;
    W[11]   = 0.4528559236327399e-1;
    W[12]   = 0.4528559236327399e-1;
    W[13]   = 0.4528559236327399e-1;
    W[14]   = 0.4528559236327399e-1;
    W[15]   = 0.4528559236327399e-1;
    W[16]   = 0.4528559236327399e-1;

	return W;
}

const vector<double>
Rule::TetW6() const
{
	vector<double> W(29);

    W[0]    = 0.9040129046014750e-1;
    W[1]    = 0.1911983427899124e-1;
    W[2]    = 0.1911983427899124e-1;
    W[3]    = 0.1911983427899124e-1;
    W[4]    = 0.1911983427899124e-1;
    W[5]    = 0.4361493840666568e-1;
    W[6]    = 0.4361493840666568e-1;
    W[7]    = 0.4361493840666568e-1;
    W[8]    = 0.4361493840666568e-1;
    W[9]    = 0.4361493840666568e-1;
    W[10]   = 0.4361493840666568e-1;
    W[11]   = 0.4361493840666568e-1;
    W[12]   = 0.4361493840666568e-1;
    W[13]   = 0.4361493840666568e-1;
    W[14]   = 0.4361493840666568e-1;
    W[15]   = 0.4361493840666568e-1;
    W[16]   = 0.4361493840666568e-1;
    W[17]   = 0.2581167596199161e-1;
    W[18]   = 0.2581167596199161e-1;
    W[19]   = 0.2581167596199161e-1;
    W[20]   = 0.2581167596199161e-1;
    W[21]   = 0.2581167596199161e-1;
    W[22]   = 0.2581167596199161e-1;
    W[23]   = 0.2581167596199161e-1;
    W[24]   = 0.2581167596199161e-1;
    W[25]   = 0.2581167596199161e-1;
    W[26]   = 0.2581167596199161e-1;
    W[27]   = 0.2581167596199161e-1;
    W[28]   = 0.2581167596199161e-1;

	return W;
}

}

}

}