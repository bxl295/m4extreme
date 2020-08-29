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
CND_CONF=Debug_mt
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/523314341/LumpedMass.o \
	${OBJECTDIR}/_ext/1466398220/Static.o \
	${OBJECTDIR}/_ext/1664855157/Operators.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Debug/GNU-Linux-x86/libM4ExtremeModel_mt.a

dist/Debug/GNU-Linux-x86/libM4ExtremeModel_mt.a: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeModel_mt.a
	${AR} -rv dist/Debug/GNU-Linux-x86/libM4ExtremeModel_mt.a ${OBJECTFILES} 
	$(RANLIB) dist/Debug/GNU-Linux-x86/libM4ExtremeModel_mt.a

${OBJECTDIR}/_ext/523314341/LumpedMass.o: ../../Model/LumpedMass/LumpedMass.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/523314341
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -I../.. -I../../External/stlib -I../../External/CommonCPP -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/523314341/LumpedMass.o ../../Model/LumpedMass/LumpedMass.cpp

${OBJECTDIR}/_ext/1466398220/Static.o: ../../Model/Static/Static.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1466398220
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -I../.. -I../../External/stlib -I../../External/CommonCPP -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1466398220/Static.o ../../Model/Static/Static.cpp

${OBJECTDIR}/_ext/1664855157/Operators.o: ../../Model/Utils/Operators/Operators.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1664855157
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../../External/SuperLU_4.0 -I../.. -I../../External/stlib -I../../External/CommonCPP -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1664855157/Operators.o ../../Model/Utils/Operators/Operators.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeModel_mt.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
