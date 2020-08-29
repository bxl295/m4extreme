// Material properties for the Johnson::Embedding class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////
//	material[0] -> Ec    [eV]
//	material[1] -> alpha [-]
//	material[2] -> beta  [-]
//	material[3] -> gamma [-]
//	material[4] -> fe    [-]
//	material[5] -> phie  [eV]

#define JEF_DATA_SETS 3

struct JEFDataType
{
	char designation[80];
	double material[6];
};

static struct JEFDataType JEFData[JEF_DATA_SETS] = 
{
	{"COPPER",	// 
		{ 3.54, 5.09, 5.85, 8.00, 0.30, 0.59 }},
                
        {"NICKEL",	// 
		{ 4.45, 4.98, 6.41, 8.86, 0.41, 0.74 }},
        
        {"PALLADIUM",	// 
		{ 3.91, 6.42, 5.91, 8.23, 0.27, 0.65 }}
};
