#!/usr/bin/env python

# first command-line argument: value of shear viscosity (for choosing results folder to process)

from numpy import *
#import sys
from subprocess import call

nevs = 1000
ebs = '0.08'
maxorder = 3
nKT = 101
eps = 0.000001
VEC_df_stem = ['']
ntrig = 1	#cos only
HBTcolsC = [3, 5, 9]
nRadiiAndSVs = len(HBTcolsC)
#HBTcolsS = [4, 6, 10]
#grid0Stem = ['', '_grid0']
#methodStem = ['GF']
methodStem = ['SVWR']
TPOstem = 'TPO_'
#KT = zeros(nKT)
KT = linspace(0.01, 1.01, nKT)

#def set_KT_pts(path, ev):
#    #infile = path + '%(ev)d/pT_pts.dat' % {"ev": ev}
#    #KT = loadtxt(infile)
#	KT = linspace(0.05, 0.7, nKT)

def collect_HBT_coefficients_fromDirec(path, allCOScoefficients, method, df_stem):
	for ev in xrange(1,nevs+1):
		# define filenames
		infile = path + '%(ev)d/' % {'ev': ev} + TPOstem + 'HBTradii_%(ms)s_cfs_ev%(ev)d' % {'ms':method, 'ev': ev} + df_stem + '.dat'
		
		# load data
		data = loadtxt(infile, usecols=tuple(HBTcolsC)).reshape([nKT,maxorder+1,nRadiiAndSVs])
		
		for iOrder in xrange(maxorder+1):
		    for iKT in xrange(nKT):
			allCOScoefficients[ev-1, iOrder, iKT] = data[iKT, iOrder]


def collect_HBT_coefficients(method, df_stem):
	direc='RESULTS_etaBYs_%(ebs)s' % {"ebs": ebs}
	subdirec="results"
	path=direc + '/' + subdirec + '-'
	allCOScoefficients = zeros([nevs, maxorder+1, nKT, nRadiiAndSVs])
	collect_HBT_coefficients_fromDirec(path, allCOScoefficients, method, df_stem)
        return allCOScoefficients


def generate_complete_FOsurface_properties(COSdata, method, df_stem):
	# outputs magnitudes of given harmonics in format:
	# K_T, R^2_s, R^2_o, R^2_l
	direc='RESULTS_etaBYs_%(ebs)s' % {"ebs": ebs}
	subdirec="results"
	path=direc + '/' + subdirec + '-'	
	#set_KT_pts(path, 1)
	allCOSresults = zeros([nKT, 2*nRadiiAndSVs+1])
	for iOrder in xrange(maxorder+1):
		for iKT in xrange(nKT):
			allCOSresults[iKT] = hstack((KT[iKT], mean(COSdata[:,iOrder, iKT],0), var(COSdata[:,iOrder, iKT],0)))
		outCOSfile = direc + '/' + TPOstem + 'complete_FOsurface_properties_%(ms)s_etaBYs_%(ebs)s_%(nevs)devs_COSneq%(ord)d' % {'ms': method, "ebs": ebs, "nevs": nevs, "ord": iOrder} + df_stem + '.dat'
		print 'Saving to', outCOSfile
		savetxt(outCOSfile, allCOSresults, fmt='%f')


if __name__ == "__main__":
		
    	for iMETHODSTEM in methodStem:
           for iDFSTEM in VEC_df_stem:
  	    	coeffs = collect_HBT_coefficients(iMETHODSTEM, iDFSTEM)
  		generate_complete_FOsurface_properties(coeffs, iMETHODSTEM, iDFSTEM)
	


# End of file
