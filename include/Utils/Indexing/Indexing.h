// Indexing.h: interface for the Indexing namespace.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(INDEXING_H__INCLUDED_)
#define INDEXING_H__INCLUDED_

#pragma once

namespace Indexing
{
	double ** New(double *,unsigned int,unsigned int);
	double *** New(double *,unsigned int,unsigned int,unsigned int);
	double **** New(double *,unsigned int,unsigned int,unsigned int,unsigned int);
	double ***** New(double *,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int);
	double ****** New(double *,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int);
	double ** New(double *,unsigned int *,unsigned int);
	double **** New(double *,unsigned int *,unsigned int *,unsigned int,unsigned int);
	void Delete(double **A,unsigned int,unsigned int);
	void Delete(double ***A,unsigned int,unsigned int,unsigned int);
	void Delete(double ****A,unsigned int,unsigned int,unsigned int,unsigned int);
	void Delete(double *****A,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int);
	void Delete(double ******A,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int);
	void Delete(double **A,unsigned int *,unsigned int);
	void Delete(double ****A,unsigned int *,unsigned int *,unsigned int,unsigned int);
}

#endif // !defined(INDEXING_H__INCLUDED_)
