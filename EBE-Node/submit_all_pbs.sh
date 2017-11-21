#!/bin/bash

#set up directory structure
baseDirectory=/home/plumberg.1/Plumberg_iEBE/iEBE-stable/EBE-Node
srcDirec=$baseDirectory/HoTCoffeeh

homeDirectory=$baseDirectory/NEW_results
rm -rf $homeDirectory
mkdir $homeDirectory

outfilename=$homeDirectory/"submit_all_pbs_jobs_record_`date +%F`.out"
outfile=`get_filename $outfilename`

jobIDsfilename=$homeDirectory/"jobIDs_`date +%F`.out"
jobIDsfile=`get_filename $jobIDsfilename`

#submit jobs
for ((i=1; i<=1; i++))
do
		npt0=5
		npphi0=6
		npy0=5
		nqt0=1
		nqx0=1
		nqy0=1
		nqz0=1
		resfrac=0.60

		workingDirectory='/home/plumberg.1/Plumberg_iEBE/iEBE-stable/all_hydro_results/results-'`echo $i`

		lwd=$homeDirectory/results-`echo $i`
		mkdir $lwd
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
