// MatLib.h: interface for the Material library.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_MATLIB_H__INCLUDED_)
#define MATERIAL_MATLIB_H__INCLUDED_

#pragma once

#include "./Material.h"
#include "./Factory.h"
#include "./Hookean/Hookean.h"
#include "./NeoHookean/Factory.h"
#include "./Gas/Factory.h"
#include "./Gas/EoS/EoSLib.h"
#include "./Uniaxial/UniLib.h"
#include "./PlaneStrain/PSLib.h"
#include "./Source/SouLib.h"
#include "./Incremental/Factory.h"
#include "./UniaxialStrain/USLib.h"
#include "./Hadamard/Factory.h"
#include "./Shear/SheLib.h"
#include "./Taylor/Factory.h"
#include "./Parallel/Factory.h"
#include "./Failure/Factory.h"
#include "./MooneyRivlin/Factory.h"
#include "./Deviatoric/FiniteKinematics/Factory.h"

#endif // !defined(MATERIAL_MATLIB_H__INCLUDED_)
