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
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/989297486/Indexing.o \
	${OBJECTDIR}/_ext/657789960/Cholesky.o \
	${OBJECTDIR}/_ext/1146615637/Crout.o \
	${OBJECTDIR}/_ext/500686949/EigenSym.o \
	${OBJECTDIR}/_ext/1695411707/LinearPolynomial2d.o \
	${OBJECTDIR}/_ext/751750679/Polynomial1d.o


# C Compiler Flags
CFLAGS=-m64

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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Debug/GNU-Linux-x86/libM4ExtremeUtils.a

dist/Debug/GNU-Linux-x86/libM4ExtremeUtils.a: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeUtils.a
	${AR} -rv dist/Debug/GNU-Linux-x86/libM4ExtremeUtils.a ${OBJECTFILES} 
	$(RANLIB) dist/Debug/GNU-Linux-x86/libM4ExtremeUtils.a

${OBJECTDIR}/_ext/989297486/Indexing.o: ../../Utils/Indexing/Indexing.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/989297486
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/989297486/Indexing.o ../../Utils/Indexing/Indexing.cpp

${OBJECTDIR}/_ext/657789960/Cholesky.o: ../../Utils/LinearAlgebra/Cholesky/Cholesky.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/657789960
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/657789960/Cholesky.o ../../Utils/LinearAlgebra/Cholesky/Cholesky.cpp

${OBJECTDIR}/_ext/1146615637/Crout.o: ../../Utils/LinearAlgebra/Crout/Crout.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1146615637
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1146615637/Crout.o ../../Utils/LinearAlgebra/Crout/Crout.cpp

${OBJECTDIR}/_ext/500686949/EigenSym.o: ../../Utils/LinearAlgebra/EigenSym/EigenSym.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/500686949
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/500686949/EigenSym.o ../../Utils/LinearAlgebra/EigenSym/EigenSym.cpp

${OBJECTDIR}/_ext/1695411707/LinearPolynomial2d.o: ../../Utils/Regression/LinearPolynomial2d/LinearPolynomial2d.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1695411707
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1695411707/LinearPolynomial2d.o ../../Utils/Regression/LinearPolynomial2d/LinearPolynomial2d.cpp

${OBJECTDIR}/_ext/751750679/Polynomial1d.o: ../../Utils/Regression/Polynomial1d/Polynomial1d.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/751750679
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/751750679/Polynomial1d.o ../../Utils/Regression/Polynomial1d/Polynomial1d.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeUtils.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
