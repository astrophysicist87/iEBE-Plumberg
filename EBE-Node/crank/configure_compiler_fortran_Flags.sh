#! /usr/bin/env bash

FC=`which ifort`;
FFLAGS=" -O3 -heap-arrays -cpp"
if [ "$FC" == "" ]; then
   FC=`which gfortran`;
   FFLAGS=" -O3 -cpp"
fi
echo $FFLAGS
