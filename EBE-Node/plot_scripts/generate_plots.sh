#!/bin/bash

directoryToProcess=../NEW_results
numberOfEvents=1000

python collect_results.py $directoryToProcess $numberOfEvents

python plot_correlation_function.py $directoryToProcess

python generate_histograms_GF.py $directoryToProcess

echo 'Finished generating plots.'
