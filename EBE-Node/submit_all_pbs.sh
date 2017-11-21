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

#run at most 12 jobs at a time
nMaxProcessesRunning=12

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
			echo 'Submitted' $i 'at' `date` >> $outfile
			#update list of currently running jobs in case I need to kill all currently running jobs
			qstat -u plumberg.1 | grep plumberg.1 | awk '$(NF-1)=="R"' | awk -F. '{print $1}' > $jobIDsfile
			sleep 3
		)

		nProcessesRunning=`qstat -u plumberg.1 | grep plumberg.1 | awk '$(NF-1)=="R"' | wc -l`
		until [ "$nProcessesRunning" -lt "$nMaxProcessesRunning" ]
		do
			echo $nProcessesRunning '==' $nMaxProcessesRunning "processes currently running at" `date` >> $outfile
			sleep 10
			nProcessesRunning=`qstat -u plumberg.1 | grep plumberg.1 | awk '$(NF-1)=="R"' | wc -l`
		done
		echo $nProcessesRunning '<' $nMaxProcessesRunning "processes currently running at" `date` >> $outfile
done

# End of file
