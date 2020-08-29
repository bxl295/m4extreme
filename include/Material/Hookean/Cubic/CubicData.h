// Material properties for the Hookean::Cubic class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

//	material[0] -> C11
//	material[1] -> C12
//	material[2] -> C44

#define HC_DATA_SETS 20

struct HCDataType
{
	char designation[80];
	double material[3];
};

static struct HCDataType HCData[HC_DATA_SETS] = 
{
	{"ALUMINUM",	// 0
		{  108.2e9,  61.3e9,  28.5e9 }},
	{"CARBON",		// 1
		{ 1076.0e9, 125.0e9, 576.0e9 }},
	{"COPPER",		// 2
		{  168.4e9, 121.4e9,  75.4e9 }},
	{"CROMIUM",		// 3
		{  350.0e9,  57.8e9, 101.0e9 }},
	{"GERMANIUM",	// 4
		{  128.9e9,  48.3e9,  67.1e9 }},
	{"GOLD",		// 5
		{  186.0e9, 157.0e9,  42.0e9 }},
	{"IRON",		// 6
		{  242.0e9, 146.5e9, 112.0e9 }},
	{"LEAD",		// 7
		{  46.6e9,  39.2e9,  14.4e9 }},
	{"LITHIUM",		// 8
		{  13.42e9, 11.30e9, 8.89e9 }},
	{"MOLYBDENUM",	// 9
		{ 460.0e9, 176.0e9, 110.0e9 }},
	{"NICKEL",		// 10
		{ 246.5e9, 147.3e9, 124.7e9 }},
	{"NIOBIUM",		// 11
		{ 246.5e9, 134.5e9, 28.73e9 }},
	{"POTASSIUM",	// 12
		{  4.57e9,  3.74e9,  2.63e9 }},
	{"SILICON",		// 13
		{ 165.7e9,  63.9e9,  79.6e9 }},
	{"SILVER",		// 14
		{ 124.0e9,  93.4e9,  46.1e9 }},
	{"SODIUM",		// 15
		{  6.03e9,  4.59e9,  5.86e9 }},
	{"THORIUM",		// 16
		{  75.3e9,  48.9e9,  47.8e9 }},
	{"TANTALUM",	// 17
		{ 266.8e9, 161.1e9, 82.49e9 }},
	{"TUNSTEN",		// 18
		{ 523.3e9, 204.5e9, 160.7e9 }},
	{"VANADIUM",	// 19
		{ 231.0e9, 120.2e9, 43.76e9 }},
};
