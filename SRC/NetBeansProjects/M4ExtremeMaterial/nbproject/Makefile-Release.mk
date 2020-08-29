#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1900861571/BlatzKo.o \
	${OBJECTDIR}/_ext/1239377311/Exponential.o \
	${OBJECTDIR}/_ext/1496518520/FK2LK.o \
	${OBJECTDIR}/_ext/602149135/FKQuadratic.o \
	${OBJECTDIR}/_ext/1492953487/Ideal.o \
	${OBJECTDIR}/_ext/1490977580/LK2FK.o \
	${OBJECTDIR}/_ext/1303900055/Polytropic.o \
	${OBJECTDIR}/_ext/373281996/Quadratic.o \
	${OBJECTDIR}/_ext/540442924/SolidPolytropic.o \
	${OBJECTDIR}/_ext/1266288529/Gas.o \
	${OBJECTDIR}/_ext/234059076/Hadamard.o \
	${OBJECTDIR}/_ext/1561980706/Cubic.o \
	${OBJECTDIR}/_ext/1077958954/Isotropic.o \
	${OBJECTDIR}/_ext/962705554/Incremental.o \
	${OBJECTDIR}/_ext/586777391/MooneyRivlin.o \
	${OBJECTDIR}/_ext/1173186849/NeoHookean.o \
	${OBJECTDIR}/_ext/1692428367/Parallel.o \
	${OBJECTDIR}/_ext/4989155/FiniteKinematics.o \
	${OBJECTDIR}/_ext/965253097/LinearizedKinematics.o \
	${OBJECTDIR}/_ext/644047494/Constant.o \
	${OBJECTDIR}/_ext/884222033/Exponential.o \
	${OBJECTDIR}/_ext/682901418/FK2LK.o \
	${OBJECTDIR}/_ext/2087146480/Constant.o \
	${OBJECTDIR}/_ext/255001977/Exponential.o \
	${OBJECTDIR}/_ext/1222444386/Polymerization.o \
	${OBJECTDIR}/_ext/2101588343/LinearThermalSoftening.o \
	${OBJECTDIR}/_ext/908839656/LogLaw.o \
	${OBJECTDIR}/_ext/1944682412/MeltingTemperature.o \
	${OBJECTDIR}/_ext/1641348407/PowerLaw.o \
	${OBJECTDIR}/_ext/1070833312/ShearViscosity.o \
	${OBJECTDIR}/_ext/1233628238/ThermalSoftening.o \
	${OBJECTDIR}/_ext/1669469797/TwoStepPowerLaw.o \
	${OBJECTDIR}/_ext/1742962370/FiniteKinematics.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -Wno-deprecated
CXXFLAGS=-m64 -Wno-deprecated

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=--64

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Release/GNU-Linux-x86/libM4ExtremeMaterial.a

dist/Release/GNU-Linux-x86/libM4ExtremeMaterial.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremeMaterial.a
	${AR} -rv dist/Release/GNU-Linux-x86/libM4ExtremeMaterial.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/libM4ExtremeMaterial.a

${OBJECTDIR}/_ext/1900861571/BlatzKo.o: ../../Material/Gas/EoS/BlatzKo/BlatzKo.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1900861571
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1900861571/BlatzKo.o ../../Material/Gas/EoS/BlatzKo/BlatzKo.cpp

${OBJECTDIR}/_ext/1239377311/Exponential.o: ../../Material/Gas/EoS/Exponential/Exponential.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1239377311
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1239377311/Exponential.o ../../Material/Gas/EoS/Exponential/Exponential.cpp

${OBJECTDIR}/_ext/1496518520/FK2LK.o: ../../Material/Gas/EoS/FK2LK/FK2LK.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1496518520
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1496518520/FK2LK.o ../../Material/Gas/EoS/FK2LK/FK2LK.cpp

${OBJECTDIR}/_ext/602149135/FKQuadratic.o: ../../Material/Gas/EoS/FKQuadratic/FKQuadratic.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/602149135
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/602149135/FKQuadratic.o ../../Material/Gas/EoS/FKQuadratic/FKQuadratic.cpp

${OBJECTDIR}/_ext/1492953487/Ideal.o: ../../Material/Gas/EoS/Ideal/Ideal.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1492953487
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1492953487/Ideal.o ../../Material/Gas/EoS/Ideal/Ideal.cpp

${OBJECTDIR}/_ext/1490977580/LK2FK.o: ../../Material/Gas/EoS/LK2FK/LK2FK.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1490977580
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1490977580/LK2FK.o ../../Material/Gas/EoS/LK2FK/LK2FK.cpp

${OBJECTDIR}/_ext/1303900055/Polytropic.o: ../../Material/Gas/EoS/Polytropic/Polytropic.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1303900055
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1303900055/Polytropic.o ../../Material/Gas/EoS/Polytropic/Polytropic.cpp

${OBJECTDIR}/_ext/373281996/Quadratic.o: ../../Material/Gas/EoS/Quadratic/Quadratic.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/373281996
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/373281996/Quadratic.o ../../Material/Gas/EoS/Quadratic/Quadratic.cpp

${OBJECTDIR}/_ext/540442924/SolidPolytropic.o: ../../Material/Gas/EoS/SolidPolytropic/SolidPolytropic.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/540442924
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/540442924/SolidPolytropic.o ../../Material/Gas/EoS/SolidPolytropic/SolidPolytropic.cpp

${OBJECTDIR}/_ext/1266288529/Gas.o: ../../Material/Gas/Gas.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1266288529
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1266288529/Gas.o ../../Material/Gas/Gas.cpp

${OBJECTDIR}/_ext/234059076/Hadamard.o: ../../Material/Hadamard/Hadamard.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/234059076
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/234059076/Hadamard.o ../../Material/Hadamard/Hadamard.cpp

${OBJECTDIR}/_ext/1561980706/Cubic.o: ../../Material/Hookean/Cubic/Cubic.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1561980706
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1561980706/Cubic.o ../../Material/Hookean/Cubic/Cubic.cpp

${OBJECTDIR}/_ext/1077958954/Isotropic.o: ../../Material/Hookean/Isotropic/Isotropic.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1077958954
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1077958954/Isotropic.o ../../Material/Hookean/Isotropic/Isotropic.cpp

${OBJECTDIR}/_ext/962705554/Incremental.o: ../../Material/Incremental/Incremental.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/962705554
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/962705554/Incremental.o ../../Material/Incremental/Incremental.cpp

${OBJECTDIR}/_ext/586777391/MooneyRivlin.o: ../../Material/MooneyRivlin/MooneyRivlin.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/586777391
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/586777391/MooneyRivlin.o ../../Material/MooneyRivlin/MooneyRivlin.cpp

${OBJECTDIR}/_ext/1173186849/NeoHookean.o: ../../Material/NeoHookean/NeoHookean.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1173186849
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1173186849/NeoHookean.o ../../Material/NeoHookean/NeoHookean.cpp

${OBJECTDIR}/_ext/1692428367/Parallel.o: ../../Material/Parallel/Parallel.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1692428367
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1692428367/Parallel.o ../../Material/Parallel/Parallel.cpp

${OBJECTDIR}/_ext/4989155/FiniteKinematics.o: ../../Material/PlaneStrain/FiniteKinematics/FiniteKinematics.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/4989155
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/4989155/FiniteKinematics.o ../../Material/PlaneStrain/FiniteKinematics/FiniteKinematics.cpp

${OBJECTDIR}/_ext/965253097/LinearizedKinematics.o: ../../Material/PlaneStrain/LinearizedKinematics/LinearizedKinematics.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/965253097
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/965253097/LinearizedKinematics.o ../../Material/PlaneStrain/LinearizedKinematics/LinearizedKinematics.cpp

${OBJECTDIR}/_ext/644047494/Constant.o: ../../Material/Shear/Constant/Constant.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/644047494
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/644047494/Constant.o ../../Material/Shear/Constant/Constant.cpp

${OBJECTDIR}/_ext/884222033/Exponential.o: ../../Material/Shear/Exponential/Exponential.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/884222033
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/884222033/Exponential.o ../../Material/Shear/Exponential/Exponential.cpp

${OBJECTDIR}/_ext/682901418/FK2LK.o: ../../Material/Shear/FK2LK/FK2LK.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/682901418
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/682901418/FK2LK.o ../../Material/Shear/FK2LK/FK2LK.cpp

${OBJECTDIR}/_ext/2087146480/Constant.o: ../../Material/Source/Constant/Constant.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2087146480
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2087146480/Constant.o ../../Material/Source/Constant/Constant.cpp

${OBJECTDIR}/_ext/255001977/Exponential.o: ../../Material/Source/Exponential/Exponential.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/255001977
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/255001977/Exponential.o ../../Material/Source/Exponential/Exponential.cpp

${OBJECTDIR}/_ext/1222444386/Polymerization.o: ../../Material/Source/Polymerization/Polymerization.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1222444386
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1222444386/Polymerization.o ../../Material/Source/Polymerization/Polymerization.cpp

${OBJECTDIR}/_ext/2101588343/LinearThermalSoftening.o: ../../Material/Uniaxial/LinearThermalSoftening/LinearThermalSoftening.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2101588343
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2101588343/LinearThermalSoftening.o ../../Material/Uniaxial/LinearThermalSoftening/LinearThermalSoftening.cpp

${OBJECTDIR}/_ext/908839656/LogLaw.o: ../../Material/Uniaxial/LogLaw/LogLaw.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/908839656
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/908839656/LogLaw.o ../../Material/Uniaxial/LogLaw/LogLaw.cpp

${OBJECTDIR}/_ext/1944682412/MeltingTemperature.o: ../../Material/Uniaxial/MeltingTemperature/MeltingTemperature.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1944682412
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1944682412/MeltingTemperature.o ../../Material/Uniaxial/MeltingTemperature/MeltingTemperature.cpp

${OBJECTDIR}/_ext/1641348407/PowerLaw.o: ../../Material/Uniaxial/PowerLaw/PowerLaw.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1641348407
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1641348407/PowerLaw.o ../../Material/Uniaxial/PowerLaw/PowerLaw.cpp

${OBJECTDIR}/_ext/1070833312/ShearViscosity.o: ../../Material/Uniaxial/ShearViscosity/ShearViscosity.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1070833312
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1070833312/ShearViscosity.o ../../Material/Uniaxial/ShearViscosity/ShearViscosity.cpp

${OBJECTDIR}/_ext/1233628238/ThermalSoftening.o: ../../Material/Uniaxial/ThermalSoftening/ThermalSoftening.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1233628238
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1233628238/ThermalSoftening.o ../../Material/Uniaxial/ThermalSoftening/ThermalSoftening.cpp

${OBJECTDIR}/_ext/1669469797/TwoStepPowerLaw.o: ../../Material/Uniaxial/TwoStepPowerLaw/TwoStepPowerLaw.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1669469797
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1669469797/TwoStepPowerLaw.o ../../Material/Uniaxial/TwoStepPowerLaw/TwoStepPowerLaw.cpp

${OBJECTDIR}/_ext/1742962370/FiniteKinematics.o: ../../Material/UniaxialStrain/FiniteKinematics/FiniteKinematics.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1742962370
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1742962370/FiniteKinematics.o ../../Material/UniaxialStrain/FiniteKinematics/FiniteKinematics.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremeMaterial.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
