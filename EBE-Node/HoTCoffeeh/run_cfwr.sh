#! /usr/bin/env bash

rm -rf ./results

#rf=0.60
rf=0.00
#for rf in 0.00 0.10
#do
	#direc=results_`echo $rf`
	for ((i=1; i<=1; i++))
	do
		mkdir results
		workingDirectory='../../RESULTS_Edec300/results/results-'`echo $i`
		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat ./results
		
		#nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=48 CF_npT=15 CF_npY=15 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=13 qxnpts=7 qynpts=7 qznpts=7 &> ./results/all.out
		nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=4 CF_npT=15 CF_npY=3 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=3 qxnpts=3 qynpts=3 qznpts=3 &> ./results/all.out
		#nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=48 CF_npT=15 CF_npY=15 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=5 qxnpts=1 qynpts=1 qznpts=5 use_extrapolation=0 &> ./results/all.out
		#nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=4 CF_npT=3 CF_npY=3 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=3 qxnpts=3 qynpts=3 qznpts=3 use_extrapolation=0 &> ./results/all.out
	
		#mv results $direc
	done
#done
