#! /usr/bin/env bash
#-------------------

sys=$1

case "$sys" in
	pp)
		for cen in '0-10%' '10-20%' '20-30%' '30-40%' '40-50%' '50-60%' '60-70%' '70-80%' '80-90%' '90-100%' '0-100%'
		do
			count=10000
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}.py \
				initial_condition_control:'centrality'=${cen} \
				superMCParameters:'Aproj'=1 \
				superMCParameters:'Atarg'=1 \
				superMCParameters:'ecm'=7000 \
				superMCParameters:'finalFactor'=80.377 \
				superMCParameters:'maxx'='10.0' \
				superMCParameters:'maxy'='10.0' \
				superMCParameters:'dx'='0.1' \
				superMCParameters:'dy'='0.1' \
				hydroParameters:'iLS'=100 \
				hydroParameters:'dx'='0.1' \
				hydroParameters:'dy'='0.1'

			./generateJobs_cluster.py 1 $count qgp \
				PlayGround_EA_${sys}_C${cen} \
				RESULTS_EA_${sys}_C${cen} \
				18:00:00 "no" \
				ParameterDict_EA_${sys}_C${cen}.py \

			./submitJobs_cluster.py

			sleep 5
		done
		;;
	pPb)
		for cen in '0-10%' '10-20%' '20-30%' '30-40%' '40-50%' '50-60%' '60-70%' '70-80%' '80-90%' '90-100%' '0-100%'
		do
			count=10000
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}.py \
				initial_condition_control:'centrality'=${cen} \
				superMCParameters:'Aproj'=1 \
				superMCParameters:'Atarg'=208 \
				superMCParameters:'ecm'=5020 \
				superMCParameters:'finalFactor'='54.3968' \
				superMCParameters:'maxx'='15.0' \
				superMCParameters:'maxy'='15.0' \
				superMCParameters:'dx'='0.1' \
				superMCParameters:'dy'='0.1' \
				hydroParameters:'iLS'=150 \
				hydroParameters:'dx'='0.1' \
				hydroParameters:'dy'='0.1'

			./generateJobs_cluster.py 1 $count qgp \
				PlayGround_EA_${sys}_C${cen} \
				RESULTS_EA_${sys}_C${cen} \
				18:00:00 "no" \
				ParameterDict_EA_${sys}_C${cen}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;
	CuCu)
		for cen in '0-10%' '10-20%' '20-30%' '30-40%' '40-50%' '50-60%' '60-70%' '70-80%' '80-90%' '90-100%' '0-100%'
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}.py \
				initial_condition_control:'centrality'=${cen} \
				superMCParameters:'Aproj'=64 \
				superMCParameters:'Atarg'=64 \
				superMCParameters:'ecm'=200 \
				superMCParameters:'finalFactor'=28.656 \
				superMCParameters:'maxx'='20.0' \
				superMCParameters:'maxy'='20.0' \
				superMCParameters:'dx'='0.1' \
				superMCParameters:'dy'='0.1' \
				hydroParameters:'iLS'=200 \
				hydroParameters:'dx'='0.1' \
				hydroParameters:'dy'='0.1'

			./generateJobs_cluster.py 1 10000 qgp \
				PlayGround_EA_${sys}_C${cen} \
				RESULTS_EA_${sys}_C${cen} \
				24:00:00 "no" \
				ParameterDict_EA_${sys}_C${cen}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;
	AuAu)
		for cen in '0-10%' '10-20%' '20-30%' '30-40%' '40-50%' '50-60%' '60-70%' '70-80%' '80-90%' '90-100%' '0-100%'
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}.py \
				initial_condition_control:'centrality'=${cen} \
				superMCParameters:'Aproj'=197 \
				superMCParameters:'Atarg'=197 \
				superMCParameters:'ecm'=200 \
				superMCParameters:'finalFactor'=28.656 \
				superMCParameters:'maxx'='20.0' \
				superMCParameters:'maxy'='20.0' \
				superMCParameters:'dx'='0.1' \
				superMCParameters:'dy'='0.1' \
				hydroParameters:'iLS'=200 \
				hydroParameters:'dx'='0.1' \
				hydroParameters:'dy'='0.1'

			./generateJobs_cluster.py 1 10000 qgp \
				PlayGround_EA_${sys}_C${cen} \
				RESULTS_EA_${sys}_C${cen} \
				24:00:00 "no" \
				ParameterDict_EA_${sys}_C${cen}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;
	PbPb)
		for cen in '0-10%' '10-20%' '20-30%' '30-40%' '40-50%' '50-60%' '60-70%' '70-80%' '80-90%' '90-100%' '0-100%'
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}.py \
				initial_condition_control:'centrality'=${cen} \
				superMCParameters:'Aproj'=208 \
				superMCParameters:'Atarg'=208 \
				superMCParameters:'ecm'=2760 \
				superMCParameters:'finalFactor'=55.1918 \
				superMCParameters:'maxx'='20.0' \
				superMCParameters:'maxy'='20.0' \
				superMCParameters:'dx'='0.1' \
				superMCParameters:'dy'='0.1' \
				hydroParameters:'iLS'=200 \
				hydroParameters:'dx'='0.1' \
				hydroParameters:'dy'='0.1'

			./generateJobs_cluster.py 1 10000 qgp \
				PlayGround_EA_${sys}_C${cen} \
				RESULTS_EA_${sys}_C${cen} \
				24:00:00 "no" \
				ParameterDict_EA_${sys}_C${cen}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;
	*)
		echo $"Usage: $0 {pp|pPb|CuCu|AuAu|PbPb}"
		exit 1
 
esac
