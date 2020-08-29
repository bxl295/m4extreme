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
	${OBJECTDIR}/_ext/1504064444/CellComplex.o \
	${OBJECTDIR}/_ext/102390297/ChainComplex.o \
	${OBJECTDIR}/_ext/1413817713/Array.o \
	${OBJECTDIR}/_ext/60986063/SparseArray.o \
	${OBJECTDIR}/_ext/43960666/SparseTable.o \
	${OBJECTDIR}/_ext/1430843110/Table.o \
	${OBJECTDIR}/_ext/2076363841/SimplicialComplex.o \
	${OBJECTDIR}/_ext/2025373939/NormalForm.o \
	${OBJECTDIR}/_ext/576602419/SparseNormalForm.o \
	${OBJECTDIR}/_ext/1801724035/GeneralGeometry.o \
	${OBJECTDIR}/_ext/1291050933/Simplex.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/Debug/GNU-Linux-x86/libM4ExtremeGeometry.a

dist/Debug/GNU-Linux-x86/libM4ExtremeGeometry.a: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeGeometry.a
	${AR} -rv dist/Debug/GNU-Linux-x86/libM4ExtremeGeometry.a ${OBJECTFILES} 
	$(RANLIB) dist/Debug/GNU-Linux-x86/libM4ExtremeGeometry.a

${OBJECTDIR}/_ext/1504064444/CellComplex.o: ../../Geometry/Algebraic/CellComplex/CellComplex.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1504064444
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1504064444/CellComplex.o ../../Geometry/Algebraic/CellComplex/CellComplex.cpp

${OBJECTDIR}/_ext/102390297/ChainComplex.o: ../../Geometry/Algebraic/ChainComplex/ChainComplex.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/102390297
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/102390297/ChainComplex.o ../../Geometry/Algebraic/ChainComplex/ChainComplex.cpp

${OBJECTDIR}/_ext/1413817713/Array.o: ../../Geometry/Algebraic/Indexed/Array/Array.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1413817713
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1413817713/Array.o ../../Geometry/Algebraic/Indexed/Array/Array.cpp

${OBJECTDIR}/_ext/60986063/SparseArray.o: ../../Geometry/Algebraic/Indexed/SparseArray/SparseArray.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/60986063
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/60986063/SparseArray.o ../../Geometry/Algebraic/Indexed/SparseArray/SparseArray.cpp

${OBJECTDIR}/_ext/43960666/SparseTable.o: ../../Geometry/Algebraic/Indexed/SparseTable/SparseTable.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/43960666
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/43960666/SparseTable.o ../../Geometry/Algebraic/Indexed/SparseTable/SparseTable.cpp

${OBJECTDIR}/_ext/1430843110/Table.o: ../../Geometry/Algebraic/Indexed/Table/Table.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1430843110
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1430843110/Table.o ../../Geometry/Algebraic/Indexed/Table/Table.cpp

${OBJECTDIR}/_ext/2076363841/SimplicialComplex.o: ../../Geometry/Algebraic/SimplicialComplex/SimplicialComplex.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2076363841
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2076363841/SimplicialComplex.o ../../Geometry/Algebraic/SimplicialComplex/SimplicialComplex.cpp

${OBJECTDIR}/_ext/2025373939/NormalForm.o: ../../Geometry/Algebraic/Utils/NormalForm/NormalForm.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2025373939
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2025373939/NormalForm.o ../../Geometry/Algebraic/Utils/NormalForm/NormalForm.cpp

${OBJECTDIR}/_ext/576602419/SparseNormalForm.o: ../../Geometry/Algebraic/Utils/SparseNormalForm/SparseNormalForm.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/576602419
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/576602419/SparseNormalForm.o ../../Geometry/Algebraic/Utils/SparseNormalForm/SparseNormalForm.cpp

${OBJECTDIR}/_ext/1801724035/GeneralGeometry.o: ../../Geometry/GeneralGeometry/GeneralGeometry.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1801724035
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1801724035/GeneralGeometry.o ../../Geometry/GeneralGeometry/GeneralGeometry.cpp

${OBJECTDIR}/_ext/1291050933/Simplex.o: ../../Geometry/Utils/Simplex/Simplex.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1291050933
	${RM} $@.d
	$(COMPILE.cc) -g -I../../External/jama -I../../External/tnt -I../.. -I../../External/stlib -Wno-deprecated -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1291050933/Simplex.o ../../Geometry/Utils/Simplex/Simplex.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/Debug/GNU-Linux-x86/libM4ExtremeGeometry.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
