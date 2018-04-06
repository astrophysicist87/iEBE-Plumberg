#! /usr/bin/env bash

workingDirectory="../EBE-Node/CHUNCHECKS/AXIS_XYZ_rho/RESFRAC_0.00/results"
qsub -v workingDirectory=$workingDirectory,currentDirectory=`pwd` run.pbs

workingDirectory="../EBE-Node/CHUNCHECKS/AXIS_XYZ_omega/RESFRAC_0.00/results"
qsub -v workingDirectory=$workingDirectory,currentDirectory=`pwd` run.pbs
