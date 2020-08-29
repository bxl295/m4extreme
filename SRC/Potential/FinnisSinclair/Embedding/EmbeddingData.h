// Material properties for the FinnisSinclair::Embedding class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

//	material[0] -> A [eV]

#define FSEF_DATA_SETS 3

struct FSEFDataType
{
	char designation[80];
	double material[1];
};

static struct FSEFDataType FSEFData[FSEF_DATA_SETS] = 
{
	{"VANADIUM",	// 
		{ 2.010637 }},

	{"NIOBIUM",	// 
		{ 3.013789 }},

	{"MOLYBDENUM",	// 
		{ 1.887117 }},
};
