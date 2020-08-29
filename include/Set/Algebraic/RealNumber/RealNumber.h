// RealNumber.h: interface for the RealNumber class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_REALNUMBER_H__INCLUDED_)
#define SET_REALNUMBER_H__INCLUDED_

#pragma once

#include <math.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////
// RealNumber functions
//////////////////////////////////////////////////////////////////////

inline void Random(double &x){x = 2.0*((double)rand()/(double)RAND_MAX) - 1.0;}
inline void Null(double &x){x = 0.0;}
inline void Identity(double &x){x = 1.0;}
inline double Inverse(const double &x){return 1.0/x;}
inline double Norm(const double &x){return fabs(x);}
inline double Involution(const double &x){return x;}
inline double RealPart(const double &x){return x;}
inline double Sqrt(const double &x){return sqrt(x);}
inline double Exp(const double &x){return exp(x);}
inline double Log(const double &x){return log(x);}

#endif // !defined(SET_REALNUMBER_H__INCLUDED_)

