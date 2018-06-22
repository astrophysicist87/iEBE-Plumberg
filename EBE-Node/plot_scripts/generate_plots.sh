#!/bin/bash

directoryToProcess=NEW_results
numberOfEvents=100

###this should average correlator, etc.
python collect_results.py $directoryToProcess $numberOfEvents

python compare_CF_resonance_vs_thermal.py $directoryToProcess

python generate_histograms_GF.py $directoryToProcess/complete_FOsurface_properties_GF_`echo $numberOfEvents`evs_COSneq0.dat
