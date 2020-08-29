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
	${OBJECTDIR}/_ext/294816988/Category.o \
	${OBJECTDIR}/_ext/2104748320/Diagonal.o \
	${OBJECTDIR}/_ext/1602231934/SkewSymmetric.o \
	${OBJECTDIR}/_ext/2121662996/Symmetric.o \
	${OBJECTDIR}/_ext/904390889/Vector.o \
	${OBJECTDIR}/_ext/639431374/Array.o \
	${OBJECTDIR}/_ext/656456771/Table.o \
	${OBJECTDIR}/_ext/853990608/Category.o \
	${OBJECTDIR}/_ext/1875076136/MapId.o \
	${OBJECTDIR}/_ext/1018913861/Cartesian.o \
	${OBJECTDIR}/_ext/2116455460/Orthonormal.o \
	${OBJECTDIR}/_ext/1626970377/LinearMapping.o \
	${OBJECTDIR}/_ext/2094882623/SymmetricSpace.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Debug/GNU-Linux-x86/libM4ExtremeSet.a

dist/Debug/GNU-Linux-x86/libM4ExtremeSet.a: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeSet.a
	${AR} -rv dist/Debug/GNU-Linux-x86/libM4ExtremeSet.a ${OBJECTFILES} 
	$(RANLIB) dist/Debug/GNU-Linux-x86/libM4ExtremeSet.a

${OBJECTDIR}/_ext/294816988/Category.o: ../../Set/Algebraic/VectorSpace/Category/Category.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/294816988
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/294816988/Category.o ../../Set/Algebraic/VectorSpace/Category/Category.cpp

${OBJECTDIR}/_ext/2104748320/Diagonal.o: ../../Set/Algebraic/VectorSpace/Category/Diagonal/Diagonal.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2104748320
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2104748320/Diagonal.o ../../Set/Algebraic/VectorSpace/Category/Diagonal/Diagonal.cpp

${OBJECTDIR}/_ext/1602231934/SkewSymmetric.o: ../../Set/Algebraic/VectorSpace/Category/SkewSymmetric/SkewSymmetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1602231934
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1602231934/SkewSymmetric.o ../../Set/Algebraic/VectorSpace/Category/SkewSymmetric/SkewSymmetric.cpp

${OBJECTDIR}/_ext/2121662996/Symmetric.o: ../../Set/Algebraic/VectorSpace/Category/Symmetric/Symmetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2121662996
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2121662996/Symmetric.o ../../Set/Algebraic/VectorSpace/Category/Symmetric/Symmetric.cpp

${OBJECTDIR}/_ext/904390889/Vector.o: ../../Set/Algebraic/VectorSpace/Vector/Vector.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/904390889
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/904390889/Vector.o ../../Set/Algebraic/VectorSpace/Vector/Vector.cpp

${OBJECTDIR}/_ext/639431374/Array.o: ../../Set/Indexed/Array/Array.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/639431374
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/639431374/Array.o ../../Set/Indexed/Array/Array.cpp

${OBJECTDIR}/_ext/656456771/Table.o: ../../Set/Indexed/Table/Table.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/656456771
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/656456771/Table.o ../../Set/Indexed/Table/Table.cpp

${OBJECTDIR}/_ext/853990608/Category.o: ../../Set/Manifold/Category/Category.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/853990608
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/853990608/Category.o ../../Set/Manifold/Category/Category.cpp

${OBJECTDIR}/_ext/1875076136/MapId.o: ../../Set/Manifold/Category/MapId/MapId.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1875076136
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1875076136/MapId.o ../../Set/Manifold/Category/MapId/MapId.cpp

${OBJECTDIR}/_ext/1018913861/Cartesian.o: ../../Set/Manifold/Euclidean/Cartesian/Cartesian.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1018913861
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1018913861/Cartesian.o ../../Set/Manifold/Euclidean/Cartesian/Cartesian.cpp

${OBJECTDIR}/_ext/2116455460/Orthonormal.o: ../../Set/Manifold/Euclidean/Orthonormal/Orthonormal.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2116455460
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2116455460/Orthonormal.o ../../Set/Manifold/Euclidean/Orthonormal/Orthonormal.cpp

${OBJECTDIR}/_ext/1626970377/LinearMapping.o: ../../Set/Manifold/LinearMapping/LinearMapping.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1626970377
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1626970377/LinearMapping.o ../../Set/Manifold/LinearMapping/LinearMapping.cpp

${OBJECTDIR}/_ext/2094882623/SymmetricSpace.o: ../../Set/Manifold/SymmetricSpace/SymmetricSpace.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2094882623
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2094882623/SymmetricSpace.o ../../Set/Manifold/SymmetricSpace/SymmetricSpace.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeSet.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
