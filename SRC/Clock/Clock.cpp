// Clock.cpp: implementation of the Clock class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./Clock.h"

Clock::Clock() 
: t(0.0), told(0.0), dt(0.0) {}

Clock::~Clock() {}

Clock::Clock(Clock &rhs) 
: t(rhs.t), told(rhs.told), dt(rhs.dt) {}

Clock *Clock::Clone()
{
	return new Clock(*this);
}

const double & 
Clock::Time() const
{
	return t;
}

double & 
Clock::Time()
{
	return t;
}

const double & 
Clock::TimeOld() const
{
	return told;
}

double & 
Clock::TimeOld()
{
	return told;
}

const double & 
Clock::DTime() const
{
	return dt;
}

double & 
Clock::DTime()
{
	return dt;
}

void 
Clock::operator++()
{
	told = t; t += dt;
}
