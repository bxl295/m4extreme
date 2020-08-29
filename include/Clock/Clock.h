// Clock.h: interface for the Clock class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(CLOCK__INCLUDED_)
#define CLOCK__INCLUDED_

#pragma once

class Clock  
{
public:

	Clock();
	virtual ~Clock();
	Clock(Clock &);
	Clock *Clone();
	const double & Time() const;
	double & Time();
	const double & TimeOld() const;
	double & TimeOld();
	const double & DTime() const;
	double & DTime();
	void operator++();

private:

	double t;
	double told;
	double dt;
};

#endif // !defined(CLOCK__INCLUDED_)
