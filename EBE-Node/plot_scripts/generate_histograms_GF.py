#!/usr/bin/env python
from numpy import *
from pylab import *
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

################################################################################
## Initialize stuff here
################################################################################
ebsvals = ['0.08']
dfstems = ['', '_no_df']

chosenOrder = 0
#chosenKTinds = [1, 3, 5, 7, 9, 11]
chosenKTinds = [0, 20, 40, 60, 80, 100]

methodStem='GF'

numberOfEvents = int(sys.argv[2])
nKT = 101
nOrder = 4    #for now
nRadii = 3  #R2o, R2s, R2l (remaining three are zero by boost-invariance or azimuthal symmetry)
nTrig = 1 # cos only
directions = ['s','o','l']

panelLabels = ['(a)', '(b)', '(c)']
panelCounter = 0

workingDirectory = sys.argv[1]
subDirectoryStem = 'results-'

#########################################################################################
## End of initializations
#########################################################################################
#########################################################################################
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")
#########################################################################################
#########################################################################################

def plot_EbE_R2ij(radiiForEvents, direction, ymin, ymax, chosenFigSize):
	# Set-up
	global panelCounter
	numplots = 6
	plotfontsize = 24
	fig, ax = plt.subplots(1, 1, figsize=chosenFigSize)
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	plt.axhline(0.0, color='black', linewidth=1.5)

	plt.tick_params(axis='both', which='major', labelsize=plotfontsize-8)
	#plt.tick_params(axis='both', which='minor', labelsize=8)
	
	R2ijKT01n = radiiForEvents[:,chosenKTinds[0]] / mean(radiiForEvents[:,chosenKTinds[0]])
	#print radiiForEvents
	#print radiiForEvents.shape
	R2ijKT02n = radiiForEvents[:,chosenKTinds[1]] / mean(radiiForEvents[:,chosenKTinds[1]])
	R2ijKT03n = radiiForEvents[:,chosenKTinds[2]] / mean(radiiForEvents[:,chosenKTinds[2]])
	R2ijKT04n = radiiForEvents[:,chosenKTinds[3]] / mean(radiiForEvents[:,chosenKTinds[3]])
	R2ijKT05n = radiiForEvents[:,chosenKTinds[4]] / mean(radiiForEvents[:,chosenKTinds[4]])
	R2ijKT06n = radiiForEvents[:,chosenKTinds[5]] / mean(radiiForEvents[:,chosenKTinds[5]])

	minKT01n, maxKT01n = min(R2ijKT01n), max(R2ijKT01n)
	minKT02n, maxKT02n = min(R2ijKT02n), max(R2ijKT02n)
	minKT03n, maxKT03n = min(R2ijKT03n), max(R2ijKT03n)
	minKT04n, maxKT04n = min(R2ijKT04n), max(R2ijKT04n)
	minKT05n, maxKT05n = min(R2ijKT05n), max(R2ijKT05n)
	minKT06n, maxKT06n = min(R2ijKT06n), max(R2ijKT06n)

	plotmin = min(array([minKT01n, minKT02n, minKT03n, minKT04n, minKT05n, minKT06n]))
	plotmax = max(array([maxKT01n, maxKT02n, maxKT03n, maxKT04n, maxKT05n, maxKT06n]))

	#print plotmin, plotmax

	nbins = 20
	dKT = (plotmax - plotmin)/nbins

	plotbins = linspace(plotmin, plotmax, num = nbins)

	#bin_centers = 0.5*(bin_edge[1:] + bin_edges[:-1])
	colormap = plt.cm.gist_ncar
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, numplots)])

	# the histogram of the data
	nKT01n, binsKT01n, patchesKT01n = ax.hist(R2ijKT01n, plotbins, histtype='step', normed=1)
	nKT02n, binsKT02n, patchesKT02n = ax.hist(R2ijKT02n, plotbins, histtype='step', normed=1)
	nKT03n, binsKT03n, patchesKT03n = ax.hist(R2ijKT03n, plotbins, histtype='step', normed=1)
	nKT04n, binsKT04n, patchesKT04n = ax.hist(R2ijKT04n, plotbins, histtype='step', normed=1)
	nKT05n, binsKT05n, patchesKT05n = ax.hist(R2ijKT05n, plotbins, histtype='step', normed=1)
	nKT06n, binsKT06n, patchesKT06n = ax.hist(R2ijKT06n, plotbins, histtype='step', normed=1)

	l01, = ax.plot(delete(binsKT01n,-1) + 0.5*dKT, nKT01n, '-o', markersize=7, label='$K_T = 0.01$ GeV')
	l02, = ax.plot(delete(binsKT02n,-1) + 0.5*dKT, nKT02n, '-s', markersize=7, label='$K_T = 0.21$ GeV')
	l03, = ax.plot(delete(binsKT03n,-1) + 0.5*dKT, nKT03n, '-p', markersize=10, label='$K_T = 0.41$ GeV')
	l04, = ax.plot(delete(binsKT04n,-1) + 0.5*dKT, nKT04n, '-D', markersize=7, label='$K_T = 0.61$ GeV')
	l05, = ax.plot(delete(binsKT05n,-1) + 0.5*dKT, nKT05n, '-*', markersize=10, label='$K_T = 0.81$ GeV')
	l06, = ax.plot(delete(binsKT06n,-1) + 0.5*dKT, nKT06n, '-^', markersize=8, label='$K_T = 1.01$ GeV')

	#adds some nice bars in the background of the histogram to give plot some structure
	dummy1, dummy2, dummy3 = ax.hist(R2ijKT01n, plotbins, histtype='bar', normed=1, edgecolor='black', color='white', alpha=0.5)

	ax.axis([plotmin, plotmax, ymin, ymax])
	#print array([minKT01n, minKT02n, minKT03n, minKT04n, minKT05n, minKT06n]), array([maxKT01n, maxKT02n, maxKT03n, maxKT04n, maxKT05n, maxKT06n])
	ax.set_xlabel(r'$R^2_{%(ij)s}\!/\left<R^2_{%(ij)s}\!\right>_{\mathrm{ev}}$' % {"ij": direction}, {'fontsize': plotfontsize + 5})
	#plt.xticks(color='k', size=plotfontsize-3)
	#plt.yticks(color='k', size=plotfontsize)
	minorLocator=MultipleLocator(0.05)
	#ax.xaxis.set_minor_locator(minorLocator)
	#ax.xaxis.set_tick_params(length=7.5, width=1.5)
	#ax.xaxis.set_tick_params(which='minor',length=5.0, width=1.0)
	if direction=='o':
		ax.legend(loc=1)
	
	ax.set_ylabel(r'P($R^2_{%(ij)s}\!/\left<R^2_{%(ij)s}\!\right>_{\mathrm{ev}}$)' % {"ij": direction}, {'fontsize': plotfontsize + 5})
	#ax.set_xticklabels(fontsize=plotfontsize-8)

	text(0.075, 0.925, panelLabels[panelCounter], ha='center', va='center', transform=ax.transAxes, fontsize=plotfontsize)
	panelCounter += 1

	#axs[0].set_ylabel(r'P($R^2_{ij}\!/\left<R^2_{ij}\!\right>_{\mathrm{ev}}$)', {'fontsize': plotfontsize + 5})
	#axs[1].set_yticklabels([])
	#axs[2].set_yticklabels([])

	#handles = [l01, l02, l03, l04, l05, l06]
	#handles = [l01, l02]
	#labels = [r'$K_T = 0.1$ GeV', r'$K_T = 0.2$ GeV', r'$K_T = 0.3$ GeV', r'$K_T = 0.4$ GeV', r'$K_T = 0.5$ GeV', r'$K_T = 0.6$ GeV']
	#fig.legend(handles, labels, bbox_to_anchor=(0.5,0.9))
	#plt.show(block=False)

	filename = workingDirectory + '/EbE_%(ms)s_R2%(ij)s%(n)i_vs_KT_%(nev)ievs.pdf'\
									% {"ij": direction, "ms": methodStem, "n": chosenOrder, "nev": numberOfEvents}
	print 'Saving to', filename
	plt.savefig(filename, format='pdf', bbox_inches='tight')

	#plt.close()




if __name__ == "__main__":
	data = zeros([numberOfEvents, nOrder, nKT, nRadii * nTrig + 2])
	currentDirectory = workingDirectory + '/' + subDirectoryStem
	for event in xrange(1,numberOfEvents+1):
		filename = currentDirectory + str(event) + '/results/HBTradii_' + methodStem + '_cfs.dat'
		data[event-1] = reshape(loadtxt(filename, usecols=(1, 2, 4, 6, 8)), [nOrder, nKT, nRadii * nTrig + 2])

	plot_EbE_R2ij(data[:,chosenOrder,:,1], 's', -0.25, 6.75, (8, 7))
	plot_EbE_R2ij(data[:,chosenOrder,:,2], 'o', -0.25, 6.75, (8, 7))
	plot_EbE_R2ij(data[:,chosenOrder,:,4], 'l', -0.25, 6.75, (8, 7))
	#pause()



# End of file
