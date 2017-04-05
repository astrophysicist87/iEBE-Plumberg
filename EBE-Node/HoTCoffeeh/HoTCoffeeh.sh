#! /usr/bin/env bash

RunSVWR=$1
RunCFWR=$2
AllArgs="${@:3}"

#echo $AllArgs

if [ "$RunSVWR" = true ]
then
	./svwr.e $AllArgs
fi

if [ "$RunCFWR" = true ]
then
	#valgrind --tool=massif ./cfwr.e $AllArgs
	./cfwr.e $AllArgs
fi
