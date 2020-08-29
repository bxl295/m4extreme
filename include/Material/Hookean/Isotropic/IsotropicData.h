// Material properties for the Hookean::Isotropic class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

//	material[0] -> First Lame constant Lambda
//	material[1] -> Second Lame constant Mu (shear modulus)

#define HI_DATA_SETS 19

struct HIDataType
{
	char designation[80];
	double material[2];
};

static struct HIDataType HIData[HI_DATA_SETS] = 
{
	{"ALUMINUM",	// 0
		{  59.3e9, 26.5e9 }},
	{"CARBON",		// 1
		{  85.0e9,536.0e9 }},
	{"COPPER",		// 2
		{ 100.6e9, 54.6e9 }},
	{"CROMIUM",		// 3
		{  77.8e9,121.0e9 }},
	{"GERMANIUM",	// 4
		{  37.6e9, 56.4e9 }},
	{"GOLD",		// 5
		{ 146.0e9, 31.0e9 }},
	{"IRON",		// 6
		{ 121.0e9, 86.0e9 }},
	{"LEAD",		// 7
		{  34.8e9, 10.1e9 }},
	{"MOLYBDENUM",	// 8
		{ 189.0e9,123.0e9 }},
	{"NICKEL",		// 9
		{ 117.0e9, 94.7e9 }},
	{"NIOBIUM",		// 10
		{ 145.0e9, 39.6e9 }},
	{"POTASSIUM",	// 11
		{  2.85e9, 1.74e9 }},
	{"SILICON",		// 12
		{  52.4e9, 68.1e9 }},
	{"SILVER",		// 13
		{  81.1e9, 33.8e9 }},
	{"SODIUM",		// 14
		{  2.35e9, 3.80e9 }},
	{"THORIUM",		// 15
		{  35.1e9, 34.0e9 }},
	{"TANTALUM",	// 16
		{ 149.0e9, 70.7e9 }},
	{"TUNSTEN",		// 17
		{ 201.0e9,160.0e9 }},
	{"VANADIUM",	// 18
		{ 124.0e9, 47.3e9 }},
};
