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
	${OBJECTDIR}/_ext/1133419759/Angular.o \
	${OBJECTDIR}/_ext/687601941/Spring.o \
	${OBJECTDIR}/_ext/1044210476/EmbeddedAtom.o \
	${OBJECTDIR}/_ext/170428714/Constant.o \
	${OBJECTDIR}/_ext/1397222683/Coulomb.o \
	${OBJECTDIR}/_ext/265784371/HalfSpace.o \
	${OBJECTDIR}/_ext/480320374/Hausdorff.o \
	${OBJECTDIR}/_ext/2036083058/TimeDependent.o \
	${OBJECTDIR}/_ext/354914912/Density.o \
	${OBJECTDIR}/_ext/1445083441/Embedding.o \
	${OBJECTDIR}/_ext/626059812/PairPotential.o \
	${OBJECTDIR}/_ext/1596633513/OneBody.o \
	${OBJECTDIR}/_ext/1446192127/Contact.o \
	${OBJECTDIR}/_ext/21131172/LennardJones.o \
	${OBJECTDIR}/_ext/853668624/Radial.o \
	${OBJECTDIR}/_ext/563632591/ThreeBody.o \
	${OBJECTDIR}/_ext/1187184835/TwoBody.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Release/GNU-Linux-x86/libM4ExtremePotential.a

dist/Release/GNU-Linux-x86/libM4ExtremePotential.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremePotential.a
	${AR} -rv dist/Release/GNU-Linux-x86/libM4ExtremePotential.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/libM4ExtremePotential.a

${OBJECTDIR}/_ext/1133419759/Angular.o: ../../Potential/Angular/Angular.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1133419759
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1133419759/Angular.o ../../Potential/Angular/Angular.cpp

${OBJECTDIR}/_ext/687601941/Spring.o: ../../Potential/Angular/Spring/Spring.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/687601941
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/687601941/Spring.o ../../Potential/Angular/Spring/Spring.cpp

${OBJECTDIR}/_ext/1044210476/EmbeddedAtom.o: ../../Potential/EmbeddedAtom/EmbeddedAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1044210476
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1044210476/EmbeddedAtom.o ../../Potential/EmbeddedAtom/EmbeddedAtom.cpp

${OBJECTDIR}/_ext/170428714/Constant.o: ../../Potential/Field/Constant/Constant.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/170428714
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/170428714/Constant.o ../../Potential/Field/Constant/Constant.cpp

${OBJECTDIR}/_ext/1397222683/Coulomb.o: ../../Potential/Field/Coulomb/Coulomb.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1397222683
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1397222683/Coulomb.o ../../Potential/Field/Coulomb/Coulomb.cpp

${OBJECTDIR}/_ext/265784371/HalfSpace.o: ../../Potential/Field/HalfSpace/HalfSpace.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/265784371
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/265784371/HalfSpace.o ../../Potential/Field/HalfSpace/HalfSpace.cpp

${OBJECTDIR}/_ext/480320374/Hausdorff.o: ../../Potential/Field/Hausdorff/Hausdorff.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/480320374
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/480320374/Hausdorff.o ../../Potential/Field/Hausdorff/Hausdorff.cpp

${OBJECTDIR}/_ext/2036083058/TimeDependent.o: ../../Potential/Field/Time_Dependent_Field/TimeDependent.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2036083058
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2036083058/TimeDependent.o ../../Potential/Field/Time_Dependent_Field/TimeDependent.cpp

${OBJECTDIR}/_ext/354914912/Density.o: ../../Potential/FinnisSinclair/Density/Density.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/354914912
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/354914912/Density.o ../../Potential/FinnisSinclair/Density/Density.cpp

${OBJECTDIR}/_ext/1445083441/Embedding.o: ../../Potential/FinnisSinclair/Embedding/Embedding.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1445083441
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1445083441/Embedding.o ../../Potential/FinnisSinclair/Embedding/Embedding.cpp

${OBJECTDIR}/_ext/626059812/PairPotential.o: ../../Potential/FinnisSinclair/PairPotential/PairPotential.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/626059812
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/626059812/PairPotential.o ../../Potential/FinnisSinclair/PairPotential/PairPotential.cpp

${OBJECTDIR}/_ext/1596633513/OneBody.o: ../../Potential/OneBody/OneBody.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1596633513
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1596633513/OneBody.o ../../Potential/OneBody/OneBody.cpp

${OBJECTDIR}/_ext/1446192127/Contact.o: ../../Potential/Radial/Contact/Contact.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1446192127
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1446192127/Contact.o ../../Potential/Radial/Contact/Contact.cpp

${OBJECTDIR}/_ext/21131172/LennardJones.o: ../../Potential/Radial/LennardJones/LennardJones.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/21131172
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/21131172/LennardJones.o ../../Potential/Radial/LennardJones/LennardJones.cpp

${OBJECTDIR}/_ext/853668624/Radial.o: ../../Potential/Radial/Radial.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/853668624
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/853668624/Radial.o ../../Potential/Radial/Radial.cpp

${OBJECTDIR}/_ext/563632591/ThreeBody.o: ../../Potential/ThreeBody/ThreeBody.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/563632591
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/563632591/ThreeBody.o ../../Potential/ThreeBody/ThreeBody.cpp

${OBJECTDIR}/_ext/1187184835/TwoBody.o: ../../Potential/TwoBody/TwoBody.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1187184835
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../.. -I../../External/jama -I../../External/tnt -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1187184835/TwoBody.o ../../Potential/TwoBody/TwoBody.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Release/GNU-Linux-x86/libM4ExtremePotential.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
