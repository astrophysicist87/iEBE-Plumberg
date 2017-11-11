#!/bin/bash

baseDirectory=/home/plumberg.1/Plumberg_iEBE/iEBE-stable/EBE-Node
homeDirectory=$baseDirectory/NEW_results
outfilename=$homeDirectory/"submit_all_pbs_jobs_record_`date +%F`.out"
jobIDsfilename=$homeDirectory/"jobIDs_`date +%F`.out"
outfile=`get_filename $outfilename`
jobIDsfile=`get_filename $jobIDsfilename`
srcDirec=$baseDirectory/HoTCoffeeh

i=1
workingDirectory='/home/plumberg.1/Plumberg_iEBE/iEBE-stable/RESULTS_Edec300/results/results-'`echo $i`

#create directories
for axis in X
do
	###########
	for nqt0 in 21
	do
		###########
		direcName0=$homeDirectory/results
		if [ ! -d "$direcName0" ]
		then
			mkdir $direcName0
		fi
		###########
		for resfrac in 0.00 0.10 0.20
		do
			direcName=$direcName0/RESFRAC_`echo $resfrac`
			if [ -d "$direcName" ]
			then
				rm -rf $direcName
			fi
			mkdir $direcName
		done
		###########
	done
	###########
done

#submit jobs
for ((i=1; i<=1000; i++))
do
		npt0=15
		npphi0=36
		npy0=15
		nqt0=17
		nqx0=7
		nqy0=7
		nqz0=7
		resfrac=0.60

		lwd=$homeDirectory/AXIS_`echo $axis`_pT`echo $npt0`_pY`echo $npy0`_qt`echo $nqt0`_CHEBINTERP/RESFRAC_`echo $resfrac`
		mkdir $lwd/results

		cp $workingDirectory/decdat2.dat $workingDirectory/decdat_mu.dat $workingDirectory/surface.dat $lwd/results

		newPBSscriptName=run_cfwr.pbs
		cp run_cfwr.pbs $lwd/$newPBSscriptName
		(
			cd $lwd

			\cp $srcDirec/cfwr.e .
			\cp $srcDirec/parameters.dat .
			\cp -r $srcDirec/EOS .

			echo 'Results directory and submission ID:' $i \
					`qsub -v workingDirectory=$lwd,NPT=$npt0,NPPHI=$npphi0,NPY=$npy0,NQT=$nqt0,NQX=$nqx0,NQY=$nqy0,NQZ=$nqz0,RESFRAC=$resfrac $newPBSscriptName` >> $outfile
			cd ..;
		)
done


qstat -u plumberg.1 >> $outfile
qstat -u plumberg.1 | grep plumberg | awk -F. '{print $1}' >> $jobIDsfile

# End of file
