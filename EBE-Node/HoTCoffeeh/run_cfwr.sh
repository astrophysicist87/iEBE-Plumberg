#! /usr/bin/env bash

rm -rf ./results

for rf in 0.00
do
	direc=results_`echo $rf`_qz_only
	rm -rf $direc
	for ((i=1; i<=1; i++))
	do
		mkdir results
		workingDirectory='../../RESULTS_Edec300/results/results-'`echo $i`
		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat ./results
		
		nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=36 CF_npT=15 CF_npY=21 qtnpts=1 qxnpts=1 qynpts=1 qznpts=101 delta_qz=0.002 \
														fit_with_projected_cfvals=1 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf`  &> ./results/all.out
		#nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=36 CF_npT=15 CF_npY=21 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=13 qxnpts=11 qynpts=11 qznpts=1 delta_qx=0.02 delta_qy=0.02 delta_qz=0.02 &> ./results/all.out
	
		#FYI: need CF_npY>=~21(ish) to get (barely) acceptable accuracy; subject to change upon calculation of full CF...

		mv results $direc
	done
done




for rf in 0.10
do
	direc=results_`echo $rf`_FULLGRID
	rm -rf $direc
	for ((i=1; i<=1; i++))
	do
		mkdir results
		workingDirectory='../../RESULTS_Edec300/results/results-'`echo $i`
		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat ./results
		
		nice -n 0 time bash ./HoTCoffeeh.sh false true CF_npphi=36 CF_npT=15 CF_npY=21 chosenParticlesMode=0 CF_resonanceThreshold=`echo $rf` qtnpts=13 qxnpts=11 qynpts=11 qznpts=11 delta_qx=0.02 delta_qy=0.02 delta_qz=0.02 &> ./results/all.out
	
		#FYI: need CF_npY>=~21(ish) to get (barely) acceptable accuracy; subject to change upon calculation of full CF...

		mv results $direc
	done
done
