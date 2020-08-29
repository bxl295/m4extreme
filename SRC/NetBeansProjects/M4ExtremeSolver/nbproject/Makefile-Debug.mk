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
	${OBJECTDIR}/_ext/1479195460/ExplicitDynamics.o \
	${OBJECTDIR}/_ext/533610305/NewtonRaphson.o \
	${OBJECTDIR}/_ext/1153682067/Secant.o \
	${OBJECTDIR}/_ext/769325372/NewtonRaphson.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Debug/GNU-Linux-x86/libM4ExtremeSolver.a

dist/Debug/GNU-Linux-x86/libM4ExtremeSolver.a: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeSolver.a
	${AR} -rv dist/Debug/GNU-Linux-x86/libM4ExtremeSolver.a ${OBJECTFILES} 
	$(RANLIB) dist/Debug/GNU-Linux-x86/libM4ExtremeSolver.a

${OBJECTDIR}/_ext/1479195460/ExplicitDynamics.o: ../../Solver/ExplicitDynamics/ExplicitDynamics.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1479195460
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1479195460/ExplicitDynamics.o ../../Solver/ExplicitDynamics/ExplicitDynamics.cpp

${OBJECTDIR}/_ext/533610305/NewtonRaphson.o: ../../Solver/LineSearch/NewtonRaphson/NewtonRaphson.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/533610305
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/533610305/NewtonRaphson.o ../../Solver/LineSearch/NewtonRaphson/NewtonRaphson.cpp

${OBJECTDIR}/_ext/1153682067/Secant.o: ../../Solver/LineSearch/Secant/Secant.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1153682067
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1153682067/Secant.o ../../Solver/LineSearch/Secant/Secant.cpp

${OBJECTDIR}/_ext/769325372/NewtonRaphson.o: ../../Solver/NewtonRaphson/NewtonRaphson.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/769325372
	${RM} $@.d
	$(COMPILE.cc) -g -I../.. -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/769325372/NewtonRaphson.o ../../Solver/NewtonRaphson/NewtonRaphson.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeSolver.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
