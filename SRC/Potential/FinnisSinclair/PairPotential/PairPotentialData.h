// Material properties for the FinnisSinclair::PairPotential class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

//	material[0] -> c     [A]
//	material[1] -> c0    [-]
//	material[2] -> c1    [-]
//	material[3] -> c2    [-]
//	material[4] -> fnd   [A]
//	material[5] -> B     [eV/A]
//      material[6] -> alpha [1/A]

#define FSPP_DATA_SETS 3

struct FSPPDataType
{
	char designation[80];
	double material[7];
};

static struct FSPPDataType FSPPData[FSPP_DATA_SETS] = 
{
	{"VANADIUM",	// 
		{ 3.8,  -0.8816318,   1.4907756, -0.397637,  2.6320, 23.0, 0.5 }},
		
	{"NIOBIUM",	// 
		{ 4.2,  -1.5640104,   2.0055779, -0.4663764, 2.8585, 48.0, 0.8 }},
		
	{"MOLYBDENUM",	// 
		{ 3.25, 43.4475218, -31.9332978,  6.0804249, 2.7255, 1223.0, 3.9 }},		
};
