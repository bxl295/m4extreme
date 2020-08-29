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
	${OBJECTDIR}/_ext/1461301167/Conforming.o \
	${OBJECTDIR}/_ext/1461301167/LumpedMass.o \
	${OBJECTDIR}/_ext/1461301167/Source.o \
	${OBJECTDIR}/_ext/1776091633/Interpolation.o \
	${OBJECTDIR}/_ext/1322199414/MLS.o \
	${OBJECTDIR}/_ext/518236901/MaxEnt.o \
	${OBJECTDIR}/_ext/911651282/Polynomial.o \
	${OBJECTDIR}/_ext/686212420/Legendre.o \
	${OBJECTDIR}/_ext/987580459/MultiIndex.o \
	${OBJECTDIR}/_ext/998436752/Polynomial.o \
	${OBJECTDIR}/_ext/999469090/Isoparametric.o \
	${OBJECTDIR}/_ext/999469090/LumpedMass.o \
	${OBJECTDIR}/_ext/999469090/Source.o \
	${OBJECTDIR}/_ext/53635318/LumpedMass.o \
	${OBJECTDIR}/_ext/53635318/MaterialPoint.o \
	${OBJECTDIR}/_ext/45780461/Gaussian.o \
	${OBJECTDIR}/_ext/446349168/Hermite.o \
	${OBJECTDIR}/_ext/2136010243/Quadrature.o \
	${OBJECTDIR}/_ext/1781950777/Simplicial.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Release/GNU-Linux-x86/libM4ExtremeElement.a

dist/Release/GNU-Linux-x86/libM4ExtremeElement.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremeElement.a
	${AR} -rv dist/Release/GNU-Linux-x86/libM4ExtremeElement.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/libM4ExtremeElement.a

${OBJECTDIR}/_ext/1461301167/Conforming.o: ../../Element/Conforming/Conforming.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1461301167
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1461301167/Conforming.o ../../Element/Conforming/Conforming.cpp

${OBJECTDIR}/_ext/1461301167/LumpedMass.o: ../../Element/Conforming/LumpedMass.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1461301167
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1461301167/LumpedMass.o ../../Element/Conforming/LumpedMass.cpp

${OBJECTDIR}/_ext/1461301167/Source.o: ../../Element/Conforming/Source.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1461301167
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1461301167/Source.o ../../Element/Conforming/Source.cpp

${OBJECTDIR}/_ext/1776091633/Interpolation.o: ../../Element/Interpolation/Interpolation.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1776091633
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1776091633/Interpolation.o ../../Element/Interpolation/Interpolation.cpp

${OBJECTDIR}/_ext/1322199414/MLS.o: ../../Element/Interpolation/MLS/MLS.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1322199414
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1322199414/MLS.o ../../Element/Interpolation/MLS/MLS.cpp

${OBJECTDIR}/_ext/518236901/MaxEnt.o: ../../Element/Interpolation/MaxEnt/MaxEnt.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/518236901
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/518236901/MaxEnt.o ../../Element/Interpolation/MaxEnt/MaxEnt.cpp

${OBJECTDIR}/_ext/911651282/Polynomial.o: ../../Element/Interpolation/Polynomial/Polynomial.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/911651282
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/911651282/Polynomial.o ../../Element/Interpolation/Polynomial/Polynomial.cpp

${OBJECTDIR}/_ext/686212420/Legendre.o: ../../Element/Interpolation/Utils/Legendre/Legendre.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/686212420
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/686212420/Legendre.o ../../Element/Interpolation/Utils/Legendre/Legendre.cpp

${OBJECTDIR}/_ext/987580459/MultiIndex.o: ../../Element/Interpolation/Utils/MultiIndex/MultiIndex.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/987580459
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/987580459/MultiIndex.o ../../Element/Interpolation/Utils/MultiIndex/MultiIndex.cpp

${OBJECTDIR}/_ext/998436752/Polynomial.o: ../../Element/Interpolation/Utils/Polynomial/Polynomial.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/998436752
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/998436752/Polynomial.o ../../Element/Interpolation/Utils/Polynomial/Polynomial.cpp

${OBJECTDIR}/_ext/999469090/Isoparametric.o: ../../Element/Isoparametric/Isoparametric.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/999469090
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/999469090/Isoparametric.o ../../Element/Isoparametric/Isoparametric.cpp

${OBJECTDIR}/_ext/999469090/LumpedMass.o: ../../Element/Isoparametric/LumpedMass.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/999469090
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/999469090/LumpedMass.o ../../Element/Isoparametric/LumpedMass.cpp

${OBJECTDIR}/_ext/999469090/Source.o: ../../Element/Isoparametric/Source.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/999469090
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/999469090/Source.o ../../Element/Isoparametric/Source.cpp

${OBJECTDIR}/_ext/53635318/LumpedMass.o: ../../Element/MaterialPoint/LumpedMass.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/53635318
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/53635318/LumpedMass.o ../../Element/MaterialPoint/LumpedMass.cpp

${OBJECTDIR}/_ext/53635318/MaterialPoint.o: ../../Element/MaterialPoint/MaterialPoint.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/53635318
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/53635318/MaterialPoint.o ../../Element/MaterialPoint/MaterialPoint.cpp

${OBJECTDIR}/_ext/45780461/Gaussian.o: ../../Element/Quadrature/Gaussian/Gaussian.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/45780461
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/45780461/Gaussian.o ../../Element/Quadrature/Gaussian/Gaussian.cpp

${OBJECTDIR}/_ext/446349168/Hermite.o: ../../Element/Quadrature/Hermite/Hermite.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/446349168
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/446349168/Hermite.o ../../Element/Quadrature/Hermite/Hermite.cpp

${OBJECTDIR}/_ext/2136010243/Quadrature.o: ../../Element/Quadrature/Quadrature.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2136010243
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2136010243/Quadrature.o ../../Element/Quadrature/Quadrature.cpp

${OBJECTDIR}/_ext/1781950777/Simplicial.o: ../../Element/Quadrature/Simplicial/Simplicial.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1781950777
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1781950777/Simplicial.o ../../Element/Quadrature/Simplicial/Simplicial.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremeElement.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
