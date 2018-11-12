#! /usr/bin/env bash

RunSVWR=$1
RunCFWR=$2
AllArgs="${@:3}"

#echo $AllArgs

if [ "$RunSVWR" = true ]
then
	#valgrind --error-limit=no --track-origins=yes --leak-check=full ./svwr.e $AllArgs
	./svwr.e $AllArgs
fi

if [ "$RunCFWR" = true ]
then
	#valgrind --tool=massif ./cfwr.e $AllArgs
	#valgrind --error-limit=no --track-origins=yes --leak-check=full ./cfwr.e $AllArgs
	#valgrind ./cfwr.e $AllArgs
	#./cfwr.e $AllArgs
        cfwrJobID=`qsub run_cfwr.pbs`
	qsub -W after:${cfwrJobID} ./done.sh
fi
