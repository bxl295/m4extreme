// Material properties for the FinnisSinclair::Density class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

//	material[0] -> d [A]

#define FSCD_DATA_SETS 3

struct FSCDDataType
{
	char designation[80];
	double material[1];
};

static struct FSCDDataType FSCDData[FSCD_DATA_SETS] = 
{
	{"VANADIUM",	// 
		{ 3.692767 }},

	{"NIOBIUM",	// 
		{ 3.915354 }},

	{"MOLYBDENUM",	// 
		{ 4.114825 }},
};
