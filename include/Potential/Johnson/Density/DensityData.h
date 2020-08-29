// Material properties for the Johnson::Density class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

//	material[0] -> fe    [-]
//	material[1] -> beta  [-]
//	material[2] -> ao    [A]

#define JCD_DATA_SETS 3

struct JCDDataType
{
	char designation[80];
	double material[3];
};

static struct JCDDataType JCDData[JCD_DATA_SETS] = 
{
	{"COPPER",	// 
		{ 0.30, 5.85, 3.61 }},
               
        {"NICKEL",	// 
		{ 0.41, 6.41, 3.52 }},
                
        {"PALLADIUM",	// 
		{ 0.27, 5.91, 3.89 }}
};
