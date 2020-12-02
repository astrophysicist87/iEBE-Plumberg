#! /usr/bin/env python
"""
    This script duplicates the EBE-Node folder and generate a collection of sbatch
    files to be batch-submitted. For efficiency all codes inside EBE-Node should
    be compiled.
"""

from sys import argv, exit
from os import makedirs, path, unlink
from shutil import copytree, copy, rmtree
from subprocess import call

from check_prerequisites import checkEnvironment, checkExecutables, greetings

# check argv
try:
    # set parameters
    numberOfJobs = int(argv[1])
    numberOfEventsPerJob = int(argv[2])
    # cluster name
    cluster_name = str(argv[3])

    # set optional parameters
    argId = 3

    argId += 1
    if len(argv)>=argId+1: # set working folder
        workingFolder = path.abspath(argv[argId])
    else:
        workingFolder = path.abspath("./PlayGround")

    argId += 1
    if len(argv)>=argId+1: # folder to store results
        resultsFolder = path.abspath(argv[argId])
    else:
        resultsFolder = path.abspath("./RESULTS")

    argId += 1
    if len(argv)>=argId+1: # set wall time
        walltime = argv[argId]
    else:
        walltime = "%d:00:00" % (3*numberOfEventsPerJob) # 3 hours per job

    argId += 1
    if len(argv)>=argId+1: # whether to compress final results folder
        parameterDictFilename = argv[argId]
    else:
        parameterDictFilename = "ParameterDict.py"
except:
    print('Usage: ./generateJobs_cluster.py number_of_jobs number_of_events_per_job cluster_name [working_folder="./PlayGround"] [results_folder="./RESULTS"] [walltime="03:00:00" (per event)] [parameterDictFilename="ParameterDict.py"]')
    exit()

# save config files
open("saved_configs.py", "w").writelines("""
iEbeConfigs = {
    "number_of_jobs"            :   %d,
    "number_of_events_per_job"  :   %d,
    "working_folder"            :   "%s",
    "results_folder"            :   "%s",
    "walltime"                  :   "%s",
    "parameter_dict_filename"   :   "%s",
}
""" % (numberOfJobs, numberOfEventsPerJob, workingFolder, resultsFolder, walltime, parameterDictFilename)
)

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

# print welcome message
print(yellow)
greetings(1)
print(purple + "\n" + "-"*80 + "\n>>>>> Welcome to the v-USPhydro event generator! <<<<<\n" + "-"*80 + normal)

# check prerequisites
print(green + "\n>>>>> Checking for required libraries <<<<<\n" + normal)
if not checkEnvironment():
    print("Prerequisites not met. Install the required library first please. Aborting.")
    exit()

# check existence of executables
print(green + "\n>>>>> Checking for existence of executables <<<<<\n" + normal)
if not checkExecutables():
    print("Not all executables can be generated. Aborting.")
    exit()

# clean up check_prerequisites.pyc
if path.exists("check_prerequisites.pyc"): unlink("check_prerequisites.pyc")

# generate events
print(green + "\n>>>>> Generating events <<<<<\n" + normal)

# prepare directories
if not path.exists(resultsFolder): makedirs(resultsFolder)
if path.exists(workingFolder): rmtree(workingFolder)
makedirs(workingFolder)

ebeNodeFolder = "EBE-Node"
crankFolderName = "crank"
crankFolder = path.join(ebeNodeFolder, crankFolderName)

# copy parameter file into the crank folder
#copy("ParameterDict.py", crankFolder)
copy(parameterDictFilename, 
     path.join(crankFolder, "ParameterDict.py"))

# backup parameter files to the result folder
copy(path.join(crankFolder, "SequentialEventDriver.py"), resultsFolder)
copy(path.join(crankFolder, "ParameterDict.py"), resultsFolder)

# duplicate EBE-Node folder to working directory, write .sbatch file
for i in range(1, numberOfJobs+1):
    targetWorkingFolder = path.join(workingFolder, "job-%d" % i)
    # copy folder
    copytree(ebeNodeFolder, targetWorkingFolder)
    open(path.join(targetWorkingFolder, "job-%d.sbatch" % i), "w").write(
"""#!/bin/bash
#SBATCH -N 1
#SBATCH -J iEBE-%d
#SBATCH -t %s
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -o iEBE.o%s
#SBATCH -e iEBE.e%s
cd %s
(cd %s
    ulimit -n 1000
    python ./SequentialEventDriver_shell.py %d 1> RunRecord.txt 2> ErrorRecord.txt
    cp RunRecord.txt ErrorRecord.txt ../finalResults/
)
mv ./finalResults %s/job-%d
""" % (i, walltime, '%j', '%j', targetWorkingFolder, crankFolderName, numberOfEventsPerJob, resultsFolder, i)
    )

from importlib import import_module
#import ParameterDict
parameterDictModule = import_module(path.splitext(parameterDictFilename)[0])
initial_condition_type = (
    parameterDictModule.initial_condition_control['initial_condition_type'])
if initial_condition_type == 'pre-generated':
    initial_file_path = (parameterDictModule.initial_condition_control[
                             'pre-generated_initial_file_path'])
    call("./copy_pre_generated_initial_conditions.sh %d %d %s %s" 
         % (numberOfJobs, numberOfEventsPerJob, initial_file_path, 
            workingFolder), shell=True)

print("Jobs generated. Submit them using submitJobs scripts.")



###########################################################################
