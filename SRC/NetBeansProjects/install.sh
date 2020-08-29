#!/bin/sh 

#
#
# Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
# All rights reserved
# see file License.txt for license details
#

cur=$(pwd)
src=$cur/..
installdir=$cur/../..

if [ $# -lt 1 ]; then
	echo "missing arguments"
	echo "./install.sh <configuration: Debug or Release> [installation path]"
	exit 1
fi

conf=$1

if [ $# -gt 1 ]; then
    installdir=$4
fi

echo 
echo "Start installation in $installdir ........."
echo

cd $installdir
libdir=$installdir/lib/$conf
if [ ! -d $libdir ]; then
    echo "cannot find the target library folder $libdir"
    echo "creating the target library folder..."
    mkdir -p $libdir
else
    echo "clean up the target library folder $libdir"
    rm -rf $libdir/*
fi

find $cur -name "*.a" | xargs -i cp {} "$libdir/"

includedir=$installdir/include
if [ ! -d $includedir ]; then
    echo "cannot find the target include folder $includedir"
    echo "creating the target include folder..."
    mkdir -p $includedir
else
    echo "clean up the target include folder $includedir"
    rm -rf $includedir/*
fi

cd $src
find . -name '*.h' -print | cpio -pmudL "$includedir/"
find . -name '*.ipp' -print | cpio -pmudL "$includedir/"
