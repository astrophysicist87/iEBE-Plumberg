#! /usr/bin/env bash

rm -rf ./results

for rf in 0.00
do
	direc=CHECK_simple_thermal_`echo $rf`
	rm -rf $direc
	for ((i=1; i<=1; i++))
	do
		mkdir results
		workingDirectory='../../RESULTS_Edec300/results/results-'`echo $i`
		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat ./results
		
		nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=12 CF_npT=5 CF_npY=5 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=13 qxnpts=7 qynpts=7 qznpts=7 delta_qx=0.02 delta_qy=0.02 delta_qz=0.02 &> ./results/all.out
		#nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=4 CF_npT=3 CF_npY=3 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` flesh_out_cf=0 qtnpts=13 qxnpts=3 qynpts=3 qznpts=3 delta_qx=0.05 delta_qy=0.05 delta_qz=0.075  &> ./results/all.out
		#nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=4 CF_npT=15 CF_npY=15 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=1 qxnpts=1 qynpts=1 qznpts=1 &> ./results/all.out
	
		#FYI: need CF_npY>=~21(ish) to get (barely) acceptable accuracy; subject to change upon calculation of full CF...

		mv results $direc
	done
done

#for rf in 0.00
#do
#	direc=CHECK_results_`echo $rf`_v2
#	rm -rf $direc
#	for ((i=1; i<=1; i++))
#	do
#		mkdir results
#		workingDirectory='../../RESULTS_Edec300/results/results-'`echo $i`
#		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat ./results
#		
#		nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=12 CF_npT=15 CF_npY=15 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=1 qxnpts=1 qynpts=1 qznpts=3 delta_qz=0.001 &> ./results/all.out
#	
#		#FYI: need CF_npY>=~21(ish) to get (barely) acceptable accuracy; subject to change upon calculation of full CF...
#
#		mv results $direc
#	done
#done
