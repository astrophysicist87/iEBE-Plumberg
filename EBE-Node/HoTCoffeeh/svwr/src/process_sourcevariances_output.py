####################################################################
# Author: Christopher Plumberg
# Date: August 24, 2015
####################################################################
# script intended for plotting/analysis of output files
# (containing interpolation grid points) from sourcevariances codes
####################################################################

from numpy import *
import sys
from os import path
#from glob import glob
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

n_resonances = 320
# contains particle_info array index of chosen resonance for analysis
chosen_resonance_to_plot = 1

####################################################################
def relSymDiff(a,b):
	return abs(200.*(a-b)/(a+b))

####################################################################
def relRtDiff(a,b):
	return abs(100.*(a-b)/(b+1.e-100))

####################################################################
def plotSVdata(givenSVdata, pxgrid, pygrid, filename):
	#fig, axs = plt.subplots(1, 1)
	plotfontsize = 12
	fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	ax = fig.gca(projection='3d')
	localcm = plt.cm.get_cmap('rainbow')
	frac = 0.9
	signs = 0.5 * frac * (sign(givenSVdata) + 1.0) + 0.5 * frac
	ax.scatter(pxgrid, pygrid, log(abs(givenSVdata)), c=signs, edgecolor='', s=5, cmap=localcm)
	ax.view_init(elev=90., azim=0.0)
	# sc = ...; plt.colorbar(sc)
	fig.suptitle(filename, fontsize=plotfontsize)
	#plt.show()

	return plt

####################################################################
def plotHBTcfsdata(HBTcfsFilepath, chosenR2ij, chosenHarmonic):
	# chosenHarmonic = 'cos' or 'sin'
	HBTcfsData = loadtxt(HBTcfsFilepath)[:,1:]	# first column is dummy index
	columnEntries = ['pT', 'order',
			'R2s(cos)', 'R2s(sin)', 'R2o(cos)', 'R2o(sin)', 'R2os(cos)', 'R2os(sin)',
			'R2l(cos)', 'R2l(sin)', 'R2sl(cos)', 'R2sl(sin)', 'R2ol(cos)', 'R2ol(sin)']
	allCols = dict([(x,y) for (x,y) in zip(columnEntries, range(len(columnEntries)))])
	chosenColumn = chosenR2ij + '(' + chosenHarmonic + ')'
	chosenR2ijIndex = allCols[chosenColumn]
	pTvals = HBTcfsData[:,0]
	chosenR2ijvals = HBTcfsData[:,chosenR2ijIndex]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(pTvals, chosenR2ijvals)
	plt.show()

	return plt


####################################################################
def oldMain():
	# set the names of all input files from the command line
	#thermal_sourcevariance_grid_filenames = map(lambda x: path.basename(x), sys.argv[1:])
	thermal_sourcevariance_grid_filenames = sys.argv[1:]
	
	# get the number of source variances to analyze
	number_of_sourcevariances = len(thermal_sourcevariance_grid_filenames)

	# current working directory
	currentWorkingDirectory = path.dirname(thermal_sourcevariance_grid_filenames[0])

	# load grid points in pT - pphi space
	# (this part will have to be modified for more general pT and pphi files)
	pTpts, pTwts = loadtxt(currentWorkingDirectory + '/pT_pts.dat').T
	pphipts, pphiwts = loadtxt(currentWorkingDirectory + '/pphi_pts.dat').T
	pTgrid, pphigrid = meshgrid(pTpts, pphipts)
	npT = len(pTpts)
	npphi = len(pphipts)
	pxgrid = pTgrid * cos(pphigrid)
	pygrid = pTgrid * sin(pphigrid)
	
	useSVdata = True

	# for each source variance, read in appropriate file and plot results
	if useSVdata:
		for iSV in range(number_of_sourcevariances):
			inFile = thermal_sourcevariance_grid_filenames[iSV]
			temp = loadtxt(inFile)
			#print temp.shape
			SVdata = loadtxt(inFile).reshape([n_resonances, npT, npphi])
			for ires in xrange(n_resonances):
				outfile = path.splitext(inFile)[0] + '_res' + str(ires) +'_reformatted.dat'
				print 'Saving to $CWD/' + path.basename(outfile)
				# this format is useful for gnuplot, which is a little easier to use for right now...
				finalResult = vstack((pxgrid.flatten(), pygrid.flatten(), SVdata[ires].flatten())).T
				savetxt(outfile, finalResult, fmt='%10.8f   '*3)
	else:
		for iSV in range(number_of_sourcevariances):
			inFile = thermal_sourcevariance_grid_filenames[iSV]
			data = loadtxt(inFile)
			outfile = path.splitext(inFile)[0] + '_reformatted.dat'
			print 'Saving to $CWD/' + path.basename(outfile)
			# this format is useful for gnuplot, which is a little easier to use for right now...
			finalResult = vstack((pxgrid.flatten(), pygrid.flatten(), data.flatten())).T
			savetxt(outfile, finalResult, fmt='%10.8f   '*3)

	print 'Finished all.'


####################################################################
def plotTestInterpolationData(testInterpolationPath):
	testInterpolationData = loadtxt(testInterpolationPath)[:,1:]	# first column is dummy index
	columnEntries = ['pT', 'pphi', 'exact', 'NEWinterp', 'OLDinterp']
	allCols = dict([(x,y) for (x,y) in zip(columnEntries, range(len(columnEntries)))])
	npt = 50
	npphi = 50
	testInterpolationData = testInterpolationData.reshape([npt,npphi,len(columnEntries)])
	
	pTvals = testInterpolationData[:,0,allCols['pT']]
	pphivals = testInterpolationData[0,:,allCols['pphi']]
	exactData = mean(testInterpolationData[:,:,allCols['exact']],1)
	NEWinterpData = mean(testInterpolationData[:,:,allCols['NEWinterp']],1)
	OLDinterpData = mean(testInterpolationData[:,:,allCols['OLDinterp']],1)
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(pTvals, exactData)
	ax.plot(pTvals, NEWinterpData)
	ax.plot(pTvals, OLDinterpData)
	#ax.plot(pTvals, relRtDiff(OLDinterpData, exactData))
	#ax.plot(pTvals, relRtDiff(NEWinterpData, exactData))
	#print relRtDiff(OLDinterpData, exactData)
	#print relRtDiff(NEWinterpData, exactData)
	ax.set_yscale('log')
	plt.show()

	return plt


####################################################################
if __name__ == "__main__":
	workingParentDirectory = '/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.08'
	for R2ij in ['R2s','R2o','R2l']:
		workingSubDirectory = 'COPY_resfrac_1.00_results-avg-1'
		workingPath = workingParentDirectory + '/' + workingSubDirectory
		plotHBTcfsdata(workingPath + '/HBTradii_cfs_ev1.dat', R2ij, 'cos')
	#plotTestInterpolationData(workingPath + '/interpolation_comparison_monval_9060225.dat')
	#oldMain()
















# End of file
