// Indexing.cpp: implementation of the Indexing namespace.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./Indexing.h"

namespace Indexing
{
double ** New(double *head,
		 unsigned int dim1, 
		 unsigned int dim2)
{
	double **A; double *p;
	unsigned i2;
	A = new double * [dim2];
	for (i2=0, p=head; i2<dim2; i2++, p+=dim1)
	{
		A[i2] = p;
	}
	return A;
}

double *** New(double *head,
		  unsigned int dim1, 
		  unsigned int dim2,
		  unsigned int dim3)
{
	double ***A; double *p;
	unsigned int i2, i3;
	A = new double ** [dim3];
	for (i3=0, p=head; i3<dim3; i3++)
	{
		A[i3] = new double * [dim2];
		for (i2=0; i2<dim2; i2++, p+=dim1)
		{
			A[i3][i2] = p;
		}
	}
	return A;
}

double **** New(double *head,
		   unsigned int dim1,
		   unsigned int dim2,
		   unsigned int dim3,
		   unsigned int dim4)
{
	double ****A; double *p;
	unsigned int i2, i3, i4;
	A = new double *** [dim4];
	for (i4=0, p=head; i4<dim4; i4++)
	{
		A[i4] = new double ** [dim3];
		for (i3=0; i3<dim3; i3++)
		{
			A[i4][i3] = new double * [dim2];
			for (i2=0; i2<dim2; i2++, p+=dim1)
			{
				A[i4][i3][i2] = p;
			}
		}
	}
	return A;
}

double ***** New(double *head,
		   unsigned int dim1,
		   unsigned int dim2,
		   unsigned int dim3,
		   unsigned int dim4,
		   unsigned int dim5)
{
	double *****A; double *p;
	unsigned int i2, i3, i4, i5;
	A = new double **** [dim5];
	for (i5=0, p=head; i5<dim5; i5++)
	{
		A[i5] = new double *** [dim4];
		for (i4=0; i4<dim4; i4++)
		{
			A[i5][i4] = new double ** [dim3];
			for (i3=0; i3<dim3; i3++)
			{
				A[i5][i4][i3] = new double * [dim2];
				for (i2=0; i2<dim2; i2++, p+=dim1)
				{
					A[i5][i4][i3][i2] = p;
				}
			}
		}
	}
	return A;
}

double ****** New(double *head,
		   unsigned int dim1,
		   unsigned int dim2,
		   unsigned int dim3,
		   unsigned int dim4,
		   unsigned int dim5,
		   unsigned int dim6)
{
	double ******A; double *p;
	unsigned int i2, i3, i4, i5, i6;
	A = new double ***** [dim6];
	for (i6=0, p=head; i6<dim6; i6++)
	{
		A[i6] = new double **** [dim5];
		for (i5=0; i5<dim5; i5++)
		{
			A[i6][i5] = new double *** [dim4];
			for (i4=0; i4<dim4; i4++)
			{
				A[i6][i5][i4] = new double ** [dim3];
				for (i3=0; i3<dim3; i3++)
				{
					A[i6][i5][i4][i3] = new double * [dim2];
					for (i2=0; i2<dim2; i2++, p+=dim1)
					{
						A[i6][i5][i4][i3][i2] = p;
					}
				}
			}
		}
	}
	return A;
}

double ** New(double *head,
		 unsigned int *dim1,
		 unsigned int dim2)
{
	double **A; double *p;
	unsigned i2;
	A = new double * [dim2];
	for (i2=0, p=head; i2<dim2; p+=dim1[i2], i2++)
	{
		A[i2] = p;
	}
	return A;
}

double **** New(double *head,
		   unsigned int *dim1,
		   unsigned int *dim2,
		   unsigned int dim3,
		   unsigned int dim4)
{
	double ****A; double *p;
	unsigned int i2, i3, i4;
	A = new double *** [dim4];
	for (i4=0, p=head; i4<dim4; i4++)
	{
		A[i4] = new double ** [dim3];
		for (i3=0; i3<dim3; i3++)
		{
			A[i4][i3] = new double * [dim2[i4]];
			for (i2=0; i2<dim2[i4]; i2++, p+=dim1[i3])
			{
				A[i4][i3][i2] = p;
			}
		}
	}
	return A;
}

void  Delete(double **A,
			 unsigned int dim1,
			 unsigned int dim2)
{
	delete [] A;
}

void  Delete(double ***A,
			 unsigned int dim1,
			 unsigned int dim2,
			 unsigned int dim3)
{	
	unsigned int i3;
	for (i3=0; i3<dim3; i3++)
		delete [] A[i3];
	delete [] A;
}

void  Delete(double ****A,
			 unsigned int dim1,
			 unsigned int dim2,
			 unsigned int dim3,
			 unsigned int dim4)
{
	unsigned int i3, i4;
	for (i4=0; i4<dim4; i4++)
		for (i3=0; i3<dim3; i3++)
			delete [] A[i4][i3];
	for (i4=0; i4<dim4; i4++)
		delete [] A[i4];
	delete [] A;
}

void  Delete(double *****A,
			 unsigned int dim1,
			 unsigned int dim2,
			 unsigned int dim3,
			 unsigned int dim4,
			 unsigned int dim5)
{
	unsigned int i3, i4, i5;
	for (i5=0; i5<dim5; i5++)
		for (i4=0; i4<dim4; i4++)
			for (i3=0; i3<dim3; i3++)
				delete [] A[i5][i4][i3];
	for (i5=0; i5<dim5; i5++)
		for (i4=0; i4<dim4; i4++)
			delete [] A[i5][i4];
	for (i5=0; i5<dim5; i5++)
		delete [] A[i5];
	delete [] A;
}

void  Delete(double ******A,
			 unsigned int dim1,
			 unsigned int dim2,
			 unsigned int dim3,
			 unsigned int dim4,
			 unsigned int dim5,
			 unsigned int dim6)
{
	unsigned int i3, i4, i5, i6;
	for (i6=0; i6<dim6; i6++)
		for (i5=0; i5<dim5; i5++)
			for (i4=0; i4<dim4; i4++)
				for (i3=0; i3<dim3; i3++)
					delete [] A[i6][i5][i4][i3];
	for (i6=0; i6<dim6; i6++)
		for (i5=0; i5<dim5; i5++)
			for (i4=0; i4<dim4; i4++)
				delete [] A[i6][i5][i4];
	for (i6=0; i6<dim6; i6++)
		for (i5=0; i5<dim5; i5++)
			delete [] A[i6][i5];
	for (i6=0; i6<dim6; i6++)
		delete [] A[i6];
	delete [] A;
}

void  Delete(double **A,
			 unsigned int *dim1,
			 unsigned int dim2)
{
	delete [] A;
}

void  Delete(double ****A,
			 unsigned int *dim1,
			 unsigned int *dim2,
			 unsigned int dim3,
			 unsigned int dim4)
{
	unsigned int i3, i4;
	for (i4=0; i4<dim4; i4++)
		for (i3=0; i3<dim3; i3++)
			delete [] A[i4][i3];
	for (i4=0; i4<dim4; i4++)
		delete [] A[i4];
	delete [] A;
}

}
