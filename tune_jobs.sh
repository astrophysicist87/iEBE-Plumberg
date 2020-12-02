#! /usr/bin/env bash
#-------------------

sys=$1

case "$sys" in
	pp)
		cen='0-100%'
		#Target dNchdeta: 6.0
		for ff in 78 79 80 81 82
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}_ff${ff}.py \
					initial_condition_control:'centrality'=${cen} \
					superMCParameters:'Aproj'=1 \
					superMCParameters:'Atarg'=1 \
					superMCParameters:'ecm'=7000 \
					superMCParameters:'finalFactor'=${ff} \
					superMCParameters:'maxx'='10.0' \
					superMCParameters:'maxy'='10.0' \
					superMCParameters:'dx'='0.1' \
					superMCParameters:'dy'='0.1' \
					hydroParameters:'iLS'=100 \
					hydroParameters:'dx'='0.1' \
					hydroParameters:'dy'='0.1' \
					HoTCoffeehControl:'runHoTCoffeeh'=False

			./generateJobs_cluster.py 1 50000 qgp \
			  PlayGround_EA_${sys}_C${cen}_ff${ff} \
			  RESULTS_EA_${sys}_C${cen}_ff${ff} \
			  03:00:00 \
			  ParameterDict_EA_${sys}_C${cen}_ff${ff}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;
	pPb)
		cen='0-100%'
		#Target dNchdeta: 17.5
		for ff in 35 40 45 50 53 54 55 56 57 60 65 70 75 80
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}_ff${ff}.py \
					initial_condition_control:'centrality'=${cen} \
					superMCParameters:'Aproj'=1 \
					superMCParameters:'Atarg'=208 \
					superMCParameters:'ecm'=5020 \
					superMCParameters:'finalFactor'=${ff} \
					superMCParameters:'maxx'='10.0' \
					superMCParameters:'maxy'='10.0' \
					superMCParameters:'dx'='0.1' \
					superMCParameters:'dy'='0.1' \
					hydroParameters:'iLS'=100 \
					hydroParameters:'dx'='0.1' \
					hydroParameters:'dy'='0.1' \
					HoTCoffeehControl:'runHoTCoffeeh'=False

			./generateJobs_cluster.py 1 50000 qgp \
			  PlayGround_EA_${sys}_C${cen}_ff${ff} \
			  RESULTS_EA_${sys}_C${cen}_ff${ff} \
			  03:00:00 \
			  ParameterDict_EA_${sys}_C${cen}_ff${ff}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;

	CuCu)
		cen='0-3%'
		#Target dNchdeta: 200ish
		for ff in 5 10 11 12 13 14 15 16 17 18 19 20 25
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}_ff${ff}.py \
					initial_condition_control:'centrality'=${cen} \
					superMCParameters:'Aproj'=63 \
					superMCParameters:'Atarg'=63 \
					superMCParameters:'ecm'=200 \
					superMCParameters:'finalFactor'=${ff} \
                	superMCParameters:'dx'='0.1' \
                	superMCParameters:'dy'='0.1' \
					superMCParameters:'maxx'='15.0' \
					superMCParameters:'maxy'='15.0' \
					hydroParameters:'iLS'=150  \
                	hydroParameters:'dx'='0.1' \
                	hydroParameters:'dy'='0.1' \
					HoTCoffeehControl:'runHoTCoffeeh'=False

			./generateJobs_cluster.py 1 10000 qgp \
			  PlayGround_EA_${sys}_C${cen}_ff${ff} \
			  RESULTS_EA_${sys}_C${cen}_ff${ff} \
			  03:00:00 \
			  ParameterDict_EA_${sys}_C${cen}_ff${ff}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;

	AuAu)
		cen='0-5%'
		#Target dNchdeta: 690(?)
		for ff in 15 20 25 26 27 28 29 30 31 32 33 34 35 40 45 50
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}_ff${ff}.py \
					initial_condition_control:'centrality'=${cen} \
					superMCParameters:'Aproj'=197 \
					superMCParameters:'Atarg'=197 \
					superMCParameters:'ecm'=200 \
					superMCParameters:'finalFactor'=${ff} \
                	superMCParameters:'dx'='0.1' \
                	superMCParameters:'dy'='0.1' \
					superMCParameters:'maxx'='20.0' \
					superMCParameters:'maxy'='20.0' \
					hydroParameters:'iLS'=200 \
                	hydroParameters:'dx'='0.1' \
                	hydroParameters:'dy'='0.1' \
					HoTCoffeehControl:'runHoTCoffeeh'=False

			./generateJobs_cluster.py 1 10000 qgp \
			  PlayGround_EA_${sys}_C${cen}_ff${ff} \
			  RESULTS_EA_${sys}_C${cen}_ff${ff} \
			  03:00:00 \
			  ParameterDict_EA_${sys}_C${cen}_ff${ff}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;

	PbPb)
		cen='0-5%'
		#Target dNchdeta: 1601
		for ff in 35 40 45 50 53 54 55 56 57 60 65 70 75 80
		do
			python updateParameterDict.py ParameterDict_EA_${sys}_C${cen}_ff${ff}.py \
					initial_condition_control:'centrality'=${cen} \
					superMCParameters:'Aproj'=208 \
					superMCParameters:'Atarg'=208 \
					superMCParameters:'ecm'=2760 \
					superMCParameters:'finalFactor'=${ff} \
					superMCParameters:'maxx'='20.0' \
					superMCParameters:'maxy'='20.0' \
                	superMCParameters:'dx'='0.1' \
                	superMCParameters:'dy'='0.1' \
					hydroParameters:'iLS'=200 \
                	hydroParameters:'dx'='0.1' \
                	hydroParameters:'dy'='0.1' \
					HoTCoffeehControl:'runHoTCoffeeh'=False

			./generateJobs_cluster.py 1 10000 qgp \
			  PlayGround_EA_${sys}_C${cen}_ff${ff} \
			  RESULTS_EA_${sys}_C${cen}_ff${ff} \
			  03:00:00 \
			  ParameterDict_EA_${sys}_C${cen}_ff${ff}.py

			./submitJobs_cluster.py

			sleep 5
		done
		;;

	*)
		echo $"Usage: $0 {pp|pPb|CuCu|AuAu|PbPb}"
		exit 1
 
esac
