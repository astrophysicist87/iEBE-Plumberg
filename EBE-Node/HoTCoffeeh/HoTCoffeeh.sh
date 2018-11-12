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

cfwrName=cfwr-`echo $$`

if [ "$RunCFWR" = true ]
then
	#valgrind --tool=massif ./cfwr.e $AllArgs
	#valgrind --error-limit=no --track-origins=yes --leak-check=full ./cfwr.e $AllArgs
	#valgrind ./cfwr.e $AllArgs
	#./cfwr.e $AllArgs
        qsub -N $cfwrName run_cfwr.pbs
	qsub -hold_jid $cfwrName -N dummy ./done.sh
fi

# forces HoTCoffeeh.sh and done.sh to wait for run_cfwr.pbs to finish
#qsub -hold_jid $cfwrName -N dummy ./done.sh
