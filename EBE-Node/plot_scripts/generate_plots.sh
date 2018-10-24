#!/bin/bash

directoryToProcess=../NEW_results
numberOfEvents=1000

python collect_results.py $directoryToProcess $numberOfEvents

python plot_correlation_function.py $directoryToProcess

python generate_histograms_GF.py $directoryToProcess $numberOfEvents

python generate_histograms.py ../SVWR_results/

python plot_CF_compare_flesh_out.py ../HoTCoffeeh/fit_correlation_function_src ../thesis_plots

python plot_CF_compare_reso_extrap.py

echo 'Finished generating plots.'
