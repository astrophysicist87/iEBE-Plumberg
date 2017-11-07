#!/bin/bash

baseDirectory=/home/plumberg.1/Plumberg_iEBE/iEBE-stable/EBE-Node
homeDirectory=$baseDirectory/code_checks
outfilename=$homeDirectory/"submit_code_checks_pbs_jobs_record_`date +%F`.out"
outfile=`get_filename $outfilename`
srcDirec=$baseDirectory/HoTCoffeeh_code_checks

i=1
workingDirectory='/home/plumberg.1/Plumberg_iEBE/iEBE-stable/RESULTS_Edec300/results/results-'`echo $i`

npt0=21
npphi0=36
npy0=15
nqt0=17

#create directories
for axis in X Y Z
do
	mkdir $homeDirectory/AXIS_`echo $axis`
	for resfrac in 0.00 0.10 0.20 0.25 0.50 0.60 0.75 1.00
	do
		mkdir $homeDirectory/AXIS_`echo $axis`/RESFRAC_`echo $resfrac`
	done
done

#submit jobs
for axis in X
do
	nqx0=11
	nqy0=1
	nqz0=1
	for resfrac in 0.00 0.10 0.20 0.25 0.50 0.60 0.75 1.00
	do
		lwd=$homeDirectory/AXIS_`echo $axis`/RESFRAC_`echo $resfrac`
		mkdir $lwd/results

		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat $lwd/results

		newPBSscriptName=run_cfwr.pbs
		cp run_cfwr.pbs $lwd/$newPBSscriptName
		(
			cd $lwd
	
			\cp $srcDirec/cfwr.e .
			\cp $srcDirec/parameters.dat .
			\cp -r $srcDirec/EOS .
	
			echo 'Results directory and submission ID:' $i `qsub -v workingDirectory=$lwd,NPT=$npt0,NPPHI=$npphi0,NPY=$npy0,NQT=$nqt0,NQX=$nqx0,NQY=$nqy0,NQZ=$nqz0,RESFRAC=$resfrac $newPBSscriptName` >> $outfile
			cd ..;
		)
	done
done

for axis in Y
do
	nqx0=1
	nqy0=11
	nqz0=1
	for resfrac in 0.00 0.10 0.20 0.25 0.50 0.60 0.75 1.00
	do
		lwd=$homeDirectory/AXIS_`echo $axis`/RESFRAC_`echo $resfrac`
		mkdir $lwd/results

		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat $lwd/results

		newPBSscriptName=run_cfwr.pbs
		cp run_cfwr.pbs $lwd/$newPBSscriptName
		(
			cd $lwd
	
			\cp $srcDirec/cfwr.e .
			\cp $srcDirec/parameters.dat .
			\cp -r $srcDirec/EOS .
	
			echo 'Results directory and submission ID:' $i `qsub -v workingDirectory=$lwd,NPT=$npt0,NPPHI=$npphi0,NPY=$npy0,NQT=$nqt0,NQX=$nqx0,NQY=$nqy0,NQZ=$nqz0,RESFRAC=$resfrac $newPBSscriptName` >> $outfile
			cd ..;
		)
	done
done

for axis in Z
do
	nqx0=1
	nqy0=1
	nqz0=11
	for resfrac in 0.00 0.10 0.20 0.25 0.50 0.60 0.75 1.00
	do
		lwd=$homeDirectory/AXIS_`echo $axis`/RESFRAC_`echo $resfrac`
		mkdir $lwd/results

		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat $lwd/results

		newPBSscriptName=run_cfwr.pbs
		cp run_cfwr.pbs $lwd/$newPBSscriptName
		(
			cd $lwd
	
			\cp $srcDirec/cfwr.e .
			\cp $srcDirec/parameters.dat .
			\cp -r $srcDirec/EOS .
	
			echo 'Results directory and submission ID:' $i `qsub -v workingDirectory=$lwd,NPT=$npt0,NPPHI=$npphi0,NPY=$npy0,NQT=$nqt0,NQX=$nqx0,NQY=$nqy0,NQZ=$nqz0,RESFRAC=$resfrac $newPBSscriptName` >> $outfile
			cd ..;
		)
	done
done

qstat -u plumberg.1 >> $outfile
