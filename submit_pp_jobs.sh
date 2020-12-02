#! /usr/bin/env bash

for cen in '0-10%' '10-20%' '20-30%' '30-40%' '40-50%' '50-60%' '60-70%' '70-80%' '80-90%' '90-100%'
do
	python updateParameterDict.py ParameterDict_EA_pp_C${cen}_run2.py initial_condition_control:'centrality'=${cen}
	./generateJobs_local.py 1 10000 PlayGround_EA_pp_C${cen}_run2 RESULTS_EA_pp_C${cen}_run2 03:00:00 no ParameterDict_EA_pp_C${cen}_run2.py
	./submitJobs_local.py
	sleep 5
done

echo 'Finished all in' `pwd`
