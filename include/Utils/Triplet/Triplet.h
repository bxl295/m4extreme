// Triplet.h: interface for the Triplet class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(TRIPLE_H_INCLUDED_)
#define TRIPLE_H_INCLUDED_

#pragma once

using namespace std;

template<class first_type, class second_type, class third_type> class triplet
{
public: 

	triplet() {}
	triplet(const first_type &first_, const second_type &second_, const third_type &third_) :
		first(first_), second(second_), third(third_) {}
	virtual ~triplet() {}
	triplet(const triplet<first_type,second_type,third_type> &t) :
		first(t.first), second(t.second), third(t.third) {}

	first_type first;
	second_type second;
	third_type third;

private:

	triplet<first_type,second_type,third_type> & 
		operator = (const triplet<first_type,second_type,third_type> &);
};

template<class first_type, class second_type, class third_type> 
triplet<first_type, second_type, third_type>
make_triplet(const first_type &first_, const second_type &second_, const third_type &third_)
{
	return triplet<first_type, second_type, third_type>(first_,second_,third_); 
}

#endif // !defined(TRIPLE_H_INCLUDED_)
