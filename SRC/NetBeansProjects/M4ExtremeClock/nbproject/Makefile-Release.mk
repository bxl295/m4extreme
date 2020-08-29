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
	${OBJECTDIR}/_ext/1590049358/Clock.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Release/GNU-Linux-x86/libM4ExtremeClock.a

dist/Release/GNU-Linux-x86/libM4ExtremeClock.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremeClock.a
	${AR} -rv dist/Release/GNU-Linux-x86/libM4ExtremeClock.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/libM4ExtremeClock.a

${OBJECTDIR}/_ext/1590049358/Clock.o: nbproject/Makefile-${CND_CONF}.mk ../../Clock/Clock.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1590049358
	${RM} $@.d
	$(COMPILE.cc) -O3 -w -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1590049358/Clock.o ../../Clock/Clock.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremeClock.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
