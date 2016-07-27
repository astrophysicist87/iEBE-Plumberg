#! /usr/bin/env bash

RunSVWR=$1
RunCFWR=$2
AllArgs="${@:3}"

if [ "$RunSVWR" = true ]
then
	./svwr.e $AllArgs
fi

if [ "$RunCFWR" = true ]
then
	./cfwr.e $AllArgs
fi
