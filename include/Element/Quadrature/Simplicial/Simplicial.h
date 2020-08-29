// Simplicial.h: interface for the Rule class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_QUADRATURE_SIMPLICIAL__INCLUDED_)
#define ELEMENT_QUADRATURE_SIMPLICIAL__INCLUDED_

#pragma once

#include "../Quadrature.h"

namespace Element
{
namespace Quadrature
{
namespace Simplicial
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

	const vector<Set::VectorSpace::Vector> TriQ1() const;
	const vector<Set::VectorSpace::Vector> TriQ2() const;
	const vector<Set::VectorSpace::Vector> TriQ3() const;
	const vector<Set::VectorSpace::Vector> TriQ4() const;
	const vector<Set::VectorSpace::Vector> TriQ5() const;
	const vector<Set::VectorSpace::Vector> TriQ6() const;
	const vector<Set::VectorSpace::Vector> TriQ7() const;
	const vector<Set::VectorSpace::Vector> TriQ8() const;
	const vector<Set::VectorSpace::Vector> TriQ9() const;

	const vector<double> TriW1() const;
	const vector<double> TriW2() const;
	const vector<double> TriW3() const;
	const vector<double> TriW4() const;
	const vector<double> TriW5() const;
	const vector<double> TriW6() const;
	const vector<double> TriW7() const;
	const vector<double> TriW8() const;
	const vector<double> TriW9() const;

	const vector<Set::VectorSpace::Vector> TetQ1() const;
	const vector<Set::VectorSpace::Vector> TetQ2() const;
	const vector<Set::VectorSpace::Vector> TetQ3() const;
	const vector<Set::VectorSpace::Vector> TetQ4() const;
	const vector<Set::VectorSpace::Vector> TetQ5() const;
	const vector<Set::VectorSpace::Vector> TetQ6() const;

	const vector<double> TetW1() const;
	const vector<double> TetW2() const;
	const vector<double> TetW3() const;
	const vector<double> TetW4() const;
	const vector<double> TetW5() const;
	const vector<double> TetW6() const;

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

#endif // !defined(ELEMENT_QUADRATURE_SIMPLICIAL__INCLUDED_)
