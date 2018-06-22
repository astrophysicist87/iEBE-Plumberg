#!/usr/bin/env python

# first command-line argument: results directory to process and store output in

from numpy import *
import sys
from subprocess import call

direc = sys.argv[1]

#nevs = 100
nevs = int(sys.argv[2])
maxorder = 3
nKT = 101
npT, npphi = 15, 36
nqx, nqy, nqz = 7, 7, 7
eps = 0.000001
VEC_df_stem = ['']
ntrig = 1	#cos only
HBTcolsC = [2, 4, 8]
nRadiiAndSVs = len(HBTcolsC)
nCFcols = 10
#HBTcolsS = [4, 6, 10]
#grid0Stem = ['', '_grid0']
methodStem = ['GF']
#methodStem = ['SVWR']
KT = linspace(0.01, 1.01, nKT)

###############################
###############################
###############################
###############################
###############################
## HBT coefficients
def collect_HBT_coefficients_fromDirec(path, allCOScoefficients, method, df_stem):
	for ev in xrange(1,nevs+1):
		# define filenames
		infile = path + '%(ev)d/' % {'ev': ev} + 'results/HBTradii_%(ms)s_cfs' % {'ms':method} + df_stem + '.dat'
		
		# load data
		data = loadtxt(infile, usecols=tuple(HBTcolsC)).reshape([nKT,maxorder+1,nRadiiAndSVs])
		
		for iOrder in xrange(maxorder+1):
		    for iKT in xrange(nKT):
			allCOScoefficients[ev-1, iOrder, iKT] = data[iKT, iOrder]


def collect_HBT_coefficients(method, df_stem):
	subdirec="results"
	path=direc + '/' + subdirec + '-'
	allCOScoefficients = zeros([nevs, maxorder+1, nKT, nRadiiAndSVs])
	collect_HBT_coefficients_fromDirec(path, allCOScoefficients, method, df_stem)
	return allCOScoefficients

def generate_complete_FOsurface_properties(COSdata, method, df_stem):
	# outputs magnitudes of given harmonics in format:
	# K_T, R^2_s, R^2_o, R^2_l
	subdirec="results"
	path=direc + '/' + subdirec + '-'	
	#set_KT_pts(path, 1)
	allCOSresults = zeros([nKT, 2*nRadiiAndSVs+1])
	for iOrder in xrange(maxorder+1):
		for iKT in xrange(nKT):
			allCOSresults[iKT] = hstack((KT[iKT], mean(COSdata[:,iOrder, iKT],0), var(COSdata[:,iOrder, iKT],0)))
		outCOSfile = direc + '/' + 'complete_FOsurface_properties_%(ms)s_%(nevs)devs_COSneq%(ord)d' % {'ms': method, "nevs": nevs, "ord": iOrder} + df_stem + '.dat'
		print 'Saving to', outCOSfile
		savetxt(outCOSfile, allCOSresults, fmt='%f')

###############################
###############################
###############################
###############################
###############################
## HBT correlation functions
def collect_HBT_CFs_fromDirec(path):
	# define filenames
	infile = path + '%(ev)d/' % {'ev': 1} + 'results/correlfunct3D_Pion_+.dat'
	print 'Doing', 1
	# load data
	data = loadtxt(infile)
	
	CFnums = (data[:,5]**2) * ( data[:,6] + data[:,7] + data[:,8] )
	CFdens = (data[:,5]**2)

	for ev in xrange(2,nevs+1):
		# define filenames
		infile = path + '%(ev)d/' % {'ev': ev} + 'results/correlfunct3D_Pion_+.dat'
		print 'Doing', ev
		# load data
		data = loadtxt(infile)
		
		CFnums += (data[:,5]**2) * ( data[:,6] + data[:,7] + data[:,8] )
		CFdens += (data[:,5]**2)
	
	averageCF = data[:,0:6]
	averageCF[:,5] = 1.0 + CFnums / CFdens
	return averageCF



def generate_average_HBTCF(method, df_stem):
	subdirec="results"
	path=direc + '/' + subdirec + '-'
	averageCF = collect_HBT_CFs_fromDirec(path)
	
	outaverageCFfile = direc + '/' + 'avgCorrelationFunction_Pion_+_%(nevs)devs' % {'ms': method, "nevs": nevs} + df_stem + '.dat'
	print 'Saving to', outaverageCFfile
	savetxt(outaverageCFfile, averageCF, fmt='%14.8e   %14.8e   %9.3e   %9.3e   %9.3e   %12.8e')




if __name__ == "__main__":
	for iMETHODSTEM in methodStem:
		for iDFSTEM in VEC_df_stem:
			coeffs = collect_HBT_coefficients(iMETHODSTEM, iDFSTEM)
			generate_complete_FOsurface_properties(coeffs, iMETHODSTEM, iDFSTEM)
			generate_average_HBTCF(iMETHODSTEM, iDFSTEM)
	


# End of file
