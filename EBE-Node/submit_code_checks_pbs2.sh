#!/bin/bash

baseDirectory=/home/plumberg.1/Plumberg_iEBE/iEBE-stable/EBE-Node
homeDirectory=$baseDirectory/code_checks2
outfilename=$homeDirectory/"submit_code_checks_pbs2_jobs_record_`date +%F`.out"
jobIDsfilename=$homeDirectory/"jobIDs_`date +%F`.out"
outfile=`get_filename $outfilename`
srcDirec=$baseDirectory/HoTCoffeeh

i=1
workingDirectory='/home/plumberg.1/Plumberg_iEBE/iEBE-stable/RESULTS_Edec300/results/results-'`echo $i`

npt0=15
npphi0=36
npy0=15

#create directories
for axis in XYZ
do
	###########
	for nqt0 in 7
	do
		###########
		direcName0=$homeDirectory/AXIS_`echo $axis`_pT`echo $npt0`_pY`echo $npy0`_qt`echo $nqt0`
		if [ ! -d "$direcName0" ]
		then
			mkdir $direcName0
		fi
		###########
		for resfrac in 0.00
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
for axis in XYZ
do
	for nqt0 in 7
	do
		nqx0=3
		nqy0=3
		nqz0=3
		for resfrac in 0.00
		do
			lwd=$homeDirectory/AXIS_`echo $axis`_pT`echo $npt0`_pY`echo $npy0`_qt`echo $nqt0`/RESFRAC_`echo $resfrac`
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
	done
done


qstat -u plumberg.1 >> $outfile
qstat -u plumberg.1 | grep plumberg | awk -F. '{print $1}' >> $jobIDsfilename

# End of file
