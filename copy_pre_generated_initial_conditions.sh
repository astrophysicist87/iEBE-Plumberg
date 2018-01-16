#!/usr/bin/env bash

number_of_jobs=$1
number_of_events_per_job=$2
initial_condition_path=$3
running_folder=$4

usage="./copy_pre_generated_initial_conditions number_of_jobs number_of_events_per_job initial_condition_path running_folder"

if [ -z "$number_of_jobs" ]; then
    echo $usage
    exit 1
fi
if [ -z "$number_of_events_per_job" ]; then
    echo $usage
    exit 1
fi
if [ -z "$initial_condition_path" ]; then
    echo $usage
    exit 1
fi
if [ -z "$running_folder" ]; then
    echo $usage
    exit 1
fi

event_id=1
for ijob in $(seq 1 $number_of_jobs)
do
    mkdir $running_folder/job-$ijob/initial_conditions
    for ievent in $(seq 1 $number_of_events_per_job)
    do
        #cp $initial_condition_path/sdAvg_order_$event_id\_block.dat $running_folder/job-$ijob/initial_conditions/
        cp $initial_condition_path/job-$ijob/event-$ievent/sdAvg_order_$ievent\_block.dat $running_folder/job-$ijob/initial_conditions/
		echo '***' COPYING '***' $initial_condition_path/job-$ijob/event-$ievent/sdAvg_order_$ievent\_block.dat '***' TO '***' $running_folder/job-$ijob/initial_conditions/
		pwd
        ((event_id++))
    done
done
