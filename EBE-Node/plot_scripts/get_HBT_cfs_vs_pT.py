#!/usr/bin/env python

# first command-line argument: results directory to process and store output in
# second command-line argument: number of events in results directory to process
# third command-line argument: whether to use thermal or full resonance results

from numpy import *
from subprocess import call
from scipy import interpolate
import sys, gauss

### command-line arguments
direc = sys.argv[1]
nevs = int(sys.argv[2])
useThermalRadii = ( sys.argv[3] == "True" )

### initial parameters
maxorder = 0
npT, npphi = 15, 36
nqx, nqy, nqz = 7, 7, 7
VEC_df_stem = ['']
ntrig = 1	#cos only
#grid0Stem = ['', '_grid0']

methodStem = ['GF']
#methodStem = ['SVWR']

HBTcolsC = [2, 4, 8]
#HBTcolsS = [4, 6, 10]
#useThermalRadii = False
if useThermalRadii:
	HBTcolsC = [3, 5, 9]

nRadii = len(HBTcolsC)
nCols = 7
columnsToUse = [0,2,3,4,5,6,7]

nKT = 101
KTmin, KTmax = 0.01, 1.01
KT = linspace(KTmin, KTmax, nKT)
pphiPts, pphiWts = gauss.gausspts(npphi, 0, 0.0, 2.0*pi)

###############################
###############################
## do interpolation
def generate_pphiAvgd_radii_vs_dense_pT(localdata, useSmoothing):
	datax = localdata[:, 0]		# grid of pT points
	#datay = localdata[:, 1:].T	# data to be interpolated
	
	densedatax = linspace(KTmin, KTmax, nKT)
	densedatay = zeros([nKT, nCols-1])
	#f = interpolate.interp1d(datax, datay, kind='linear')
	#f = interpolate.interp1d(datax, datay, kind='cubic')
	#densedatay = f(densedatax)
	for rad in xrange(nCols-1):
		datay = localdata[:,rad+1]
		kToUse,sToUse=3,0
		if useSmoothing and rad==3:	# for R^2_{l,0} only
			kToUse,sToUse=2,1
		tck = interpolate.splrep(datax, datay, s=sToUse,k=kToUse)
		densedatay[:,rad] = interpolate.splev(densedatax, tck, der=0)
		#print densedatax.shape, densedatay.shape

	return vstack((densedatax,densedatay.T)).T

###############################
###############################
## HBT coefficients
def collect_HBT_coefficients(method, df_stem):
	subdirec="results"
	path=direc + '/' + subdirec + '-'
	for ev in xrange(1,nevs+1):
		# define filenames
		thermalStem = ''
		if useThermalRadii:
			thermalStem = '_ev%(ev)d_thermal' % {'ev': ev}
		infile = path + '%(ev)d/' % {'ev': ev} + 'results/HBTradii_%(ms)s' % {'ms': method} + thermalStem + df_stem + '_grid0.dat'
		
		# load data
		data = loadtxt(infile, usecols=tuple(columnsToUse)).reshape([npT, npphi, nCols])
		data = swapaxes(data, 1, 2)
		
		# average over pphi
		coeffs = data[:,:,0]
		coeffs[:,1:nCols] = dot(data[:,1:nCols],pphiWts) / (2.0*pi)
		
		# save to same directory
		outfile = path + '%(ev)d/' % {'ev': ev} + 'results/HBTradii_%(ms)s' % {'ms': method} + thermalStem + df_stem + '_pphiAvgd.dat'
		savetxt(outfile, coeffs, fmt = '%0.8f   %0.6f   %0.6f   %0.6f   %0.6f   %0.6f   %0.6f')
		print 'Saved to', outfile

		usingSmoothing = False
		densepTcoeffs = generate_pphiAvgd_radii_vs_dense_pT(coeffs, usingSmoothing)
		
		# save to same directory
		outfile = path + '%(ev)d/' % {'ev': ev} + 'results/HBTradii_%(ms)s_cfs' % {'ms': method} + thermalStem + df_stem + '_pphiAvgd.dat'
		savetxt(outfile, densepTcoeffs, fmt = '%0.8f   %0.6f   %0.6f   %0.6f   %0.6f   %0.6f   %0.6f')
		print 'Saved to', outfile
		
		usingSmoothing = True
		densepTcoeffs = generate_pphiAvgd_radii_vs_dense_pT(coeffs, usingSmoothing)
		
		# save to same directory
		outfile = path + '%(ev)d/' % {'ev': ev} + 'results/HBTradii_%(ms)s_smoothed_cfs' % {'ms': method} + thermalStem + df_stem + '_pphiAvgd.dat'
		savetxt(outfile, densepTcoeffs, fmt = '%0.8f   %0.6f   %0.6f   %0.6f   %0.6f   %0.6f   %0.6f')
		print 'Saved to', outfile

	return
		

###############################
###############################
if __name__ == "__main__":
	for iMETHODSTEM in methodStem:
		for iDFSTEM in VEC_df_stem:
			collect_HBT_coefficients(iMETHODSTEM, iDFSTEM)
	


# End of file
