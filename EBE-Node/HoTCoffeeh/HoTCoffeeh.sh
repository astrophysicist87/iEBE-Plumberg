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
	echo 'Submitting cfwr.e now...'
        cfwrJobID=`qsub run_cfwr.pbs`
	echo 'cfwrJobID=', $cfwrJobID
	until [ `qstat -f $cfwrJobID | grep -c 'job_state = C'` -eq 1 ]
	do
		sleep 60
	done
	echo 'Finished everything in HoTCoffeeh.sh!'
fi
