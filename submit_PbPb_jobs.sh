#! /usr/bin/env bash

# Tune "finalfactor" first
#for ff in 35 40 45 50 53 54 55 56 57 60 65 70 75 80
#do
#	python updateParameterDict.py ParameterDict_EA_PbPb_FF${ff}.py superMCParameters:finalFactor=${ff}
#	./generateJobs_local.py 1 5000 PlayGround_EA_PbPb_FF${ff} RESULTS_EA_PbPb_FF${ff} 03:00:00 no ParameterDict_EA_PbPb_FF${ff}.py
#	./submitJobs_local.py
#	sleep 5
#done
#finalFactor obtained: 55.3037

for cen in '0-100%' '0-10%' '10-20%' '20-30%' '30-40%' '40-50%' '50-60%' '60-70%' '70-80%' '80-90%' '90-100%'
do
        python updateParameterDict.py ParameterDict_EA_PbPb_C${cen}.py initial_condition_control:'centrality'=${cen}
        ./generateJobs_local.py 1 10000 PlayGround_EA_PbPb_C${cen} RESULTS_EA_PbPb_C${cen} 24:00:00 no ParameterDict_EA_PbPb_C${cen}.py
        ./submitJobs_local.py
	sleep 5
done

echo 'Finished all in' `pwd`
