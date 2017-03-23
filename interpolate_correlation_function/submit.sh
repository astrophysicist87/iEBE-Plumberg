#! /usr/bin/env bash

nohup ./interpolate_CF ../EBE-Node/HoTCoffeeh_thermal/LARGE_results_THERMAL 1> averaged_thermalCF_slices.out 2> interpolate_thermalCF.err &
