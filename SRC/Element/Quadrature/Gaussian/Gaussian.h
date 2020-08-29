// Gaussian.h: interface for the Rule class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_QUADRATURE_GAUSSIAN__INCLUDED_)
#define ELEMENT_QUADRATURE_GAUSSIAN__INCLUDED_

#pragma once

#include "../Quadrature.h"

using namespace std;

namespace Element
{
namespace Quadrature
{
namespace Gaussian
{
class Rule : public Element::Quadrature::Rule
{
public:

	Rule();
	virtual ~Rule();
	Rule(const unsigned int &, const unsigned int &);
	vector<Set::VectorSpace::Vector> GetQ() const;
	vector<double> GetW() const;

private:

	const vector<Set::VectorSpace::Vector> Q1() const;
	const vector<Set::VectorSpace::Vector> Q3() const;
	const vector<Set::VectorSpace::Vector> Q5() const;
	const vector<Set::VectorSpace::Vector> Q7() const;
	const vector<Set::VectorSpace::Vector> Q9() const;
	const vector<Set::VectorSpace::Vector> Q11() const;
        const vector<Set::VectorSpace::Vector> Q13() const;
	const vector<Set::VectorSpace::Vector> Q15() const;
	const vector<Set::VectorSpace::Vector> Q17() const;
	const vector<Set::VectorSpace::Vector> Q19() const;
	const vector<Set::VectorSpace::Vector> Q21() const;
	const vector<Set::VectorSpace::Vector> Q23() const;
        const vector<Set::VectorSpace::Vector> Q25() const;
	const vector<Set::VectorSpace::Vector> Q27() const;
	const vector<Set::VectorSpace::Vector> Q29() const;
	const vector<Set::VectorSpace::Vector> Q31() const;

	const vector<double> W1() const;
	const vector<double> W3() const;
	const vector<double> W5() const;
	const vector<double> W7() const;
	const vector<double> W9() const;
	const vector<double> W11() const;
	const vector<double> W13() const;
	const vector<double> W15() const;
	const vector<double> W17() const;
	const vector<double> W19() const;
	const vector<double> W21() const;
	const vector<double> W23() const;
	const vector<double> W25() const;
	const vector<double> W27() const;
	const vector<double> W29() const;
	const vector<double> W31() const;

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

#endif // !defined(ELEMENT_QUADRATURE_GAUSSIAN__INCLUDED_)


