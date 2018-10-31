#! /usr/bin/env bash

rm -rf ./results

for ((i=2; i<=1000; i++))
do
	mkdir results
	workingDirectory='../../RESULTS_Edec300/results/results-'`echo $i`
	cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat $workingDirectory/correlfunct3D_Pion_+.dat ./results
	
	nice -n 0 bash ./HoTCoffeeh.sh false true SV_npphi=48 chosenParticlesMode=0 qtnpts=13 delta_qz=0.015 calculate_CF_mode=2 flag_negative_S=1 CF_resonanceThreshold=0.6 particle_diff_tolerance=0.0 nKT=101 KTmin=0.01 KTmax=1.01 use_plane_psi_order=0 CF_npT=15 fit_with_projected_cfvals=1 use_log_fit=1 CF_npphi=48 flesh_out_cf=1 use_lambda=1 tolerance=0.0 qynpts=7 ignore_long_lived_resonances=1 delta_qx=0.0125 delta_qy=0.0125 include_bulk_pi=1 grouping_particles=0 SV_resonanceThreshold=0.6 max_lifetime=100.0 include_delta_f=1 SV_npT=15 use_extrapolation=1 n_order=4 qxnpts=7 qznpts=7 nKphi=48 &> ./results/tmp.out
	
	mv ./results/lambdas.dat $workingDirectory
	rm -rf ./results
done
