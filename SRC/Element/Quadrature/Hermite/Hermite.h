// Hermite.h: interface for the Rule class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_QUADRATURE_HERMITE__INCLUDED_)
#define ELEMENT_QUADRATURE_HERMITE__INCLUDED_

#pragma once

#include "../Quadrature.h"

using namespace std;

namespace Element
{
namespace Quadrature
{
namespace Hermite
{
class Rule : public Element::Quadrature::Rule
{
public:

	Rule();
	virtual ~Rule();
	Rule(const unsigned int &, const unsigned int &);
	const vector<Set::VectorSpace::Vector> GetQ();
	const vector<double> GetW();

private:
	
	// QN_n means order=N and dim=n, the same for WN_n
	const vector<Set::VectorSpace::Vector> Q3_6() const;
	const vector<Set::VectorSpace::Vector> Q3_9() const;
	const vector<Set::VectorSpace::Vector> Q3_39() const;
	const vector<Set::VectorSpace::Vector> Q5_3() const;
	const vector<Set::VectorSpace::Vector> Q5_6() const;
	
	const vector<double> W3_6() const;
	const vector<double> W3_9() const;
	const vector<double> W3_39() const;
	const vector<double> W5_3() const;
	const vector<double> W5_6() const;
	
private:

	unsigned int n;
	unsigned int N;

private:

	Rule(const Rule &);
	Rule & operator = (const Rule &);
};

}

}

}

#endif // !defined(ELEMENT_QUADRATURE_HERMITE__INCLUDED_)
