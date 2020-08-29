// Material properties for the Johnson::PairPotential class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////
//	material[0] -> phie  [eV]
//	material[1] -> gamma [-]
//	material[2] -> ao    [A]

#define JPP_DATA_SETS 3

struct JPPDataType
{
	char designation[80];
	double material[3];
};

static struct JPPDataType JPPData[JPP_DATA_SETS] = 
{
	{"COPPER",	// 
		{ 0.59, 8.00, 3.61}},
                
        {"NICKEL",	// 
		{ 0.74, 8.86, 3.52}},
         
        {"PALLADIUM",	// 
		{ 0.65, 8.23, 3.89}}        
};
