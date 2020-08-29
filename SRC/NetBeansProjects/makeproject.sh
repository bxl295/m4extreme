#!/bin/sh 

#
#
# Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
# All rights reserved
# see file License.txt for license details
#

cur=$(pwd)
src=$cur/..
clean=0
libraries=("Clock" "Element" "Potential" "Model" "Material" "Geometry" "Set" "Solver" "Utils")

if [ $# -lt 2 ]; then
	echo "missing arguments"
	echo "./makeproject <configuration> <project name> [clean: 0 or 1]"
	echo "configuration: Debug or Release"
	echo "project name : All ${libraries[@]}"
	echo
	exit 1
fi

conf=$1
project=$2

if [ $# -gt 2 ]; then
    clean=$3
fi

makeproject() {    
    if [ ! -d $projectdir ]; then
	echo "Couldnot find project folder $projectdir"
	return 1
    else
	echo
	echo "-- Start building the library $projectdir ........."
	echo
    fi

    cd $projectdir

    if [ "$clean" != "0" ]; then
	make CONF=$conf clean
    fi

    # compile and build the library
    make CONF=$conf
 
    # check if the library has been built
    count=$(find $projectdir/dist -name "*.a"|wc -l)
    if [ "$count" -lt "1"  ]; then
        echo
	echo "Failed to compile $projectdir"
        echo
	exit
    else
	echo
	echo "---------- SUCCEED ----------"
	echo
	echo
    fi

    return 0
}

if [ "$project" == "All" ]; then
    for libname in ${libraries[@]}; do
	projectdir=$cur/M4Extreme$libname
	makeproject
    done
else    
    projectdir=$cur/M4Extreme$project
    makeproject
fi
