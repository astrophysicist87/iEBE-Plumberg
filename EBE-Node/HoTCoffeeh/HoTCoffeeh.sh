#! /usr/bin/env bash

RunSVWR=$1
RunCFWR=$2

if [ "$RunSVWR" = true ]
then
	./svwr.e
fi

if [ "$RunCFWR" = true ]
then
	./cfwr.e
fi
