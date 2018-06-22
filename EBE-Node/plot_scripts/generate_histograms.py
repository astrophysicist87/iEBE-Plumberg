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
chosenKTinds = [0, 20, 40, 60, 80, 100]

methodStem='_SVWR'
#methodStem='_CFWR'

numberOfEvents = 1000
nKT = 101
nOrder = 4    #for now
nKphi = 48
nRadii = 3  #R2o, R2s, R2l (remaining two are zero by boost-invariance or azimuthal symmetry)
nTrig = 1 # cos only

panelLabels = ['(a)', '(b)']
panelCounter = 0

workingDirectory = 'RESULTS_etaBYs_%(ebsstring)s/' % {"ebsstring": ebsvals[0]}
subDirectoryStem = 'results-'

################################################################################
## End of initializations
################################################################################

def plot_EbE_R2o(ebs, dfstem, radiiForEvents):
	numplots = 6
	plotfontsize = 18
	ax=plt.axes([0.1,0.15,0.85,0.8])
	plt.axhline(0.0, color='black', linewidth=1.5)
	
	direc = 'RESULTS_etaBYs_%(ebsstring)s/' % {"ebsstring": ebs}
	
	R2oKT00n = radiiForEvents[:,0] / mean(radiiForEvents[:,0])
	R2oKT02n = radiiForEvents[:,20] / mean(radiiForEvents[:,20])
	R2oKT04n = radiiForEvents[:,40] / mean(radiiForEvents[:,40])
	R2oKT06n = radiiForEvents[:,60] / mean(radiiForEvents[:,60])
	R2oKT08n = radiiForEvents[:,80] / mean(radiiForEvents[:,80])
	R2oKT10n = radiiForEvents[:,100] / mean(radiiForEvents[:,100])

	minKT00n, maxKT00n = min(R2oKT00n), max(R2oKT00n)
	minKT02n, maxKT02n = min(R2oKT02n), max(R2oKT02n)
	minKT04n, maxKT04n = min(R2oKT04n), max(R2oKT04n)
	minKT06n, maxKT06n = min(R2oKT06n), max(R2oKT06n)
	minKT08n, maxKT08n = min(R2oKT08n), max(R2oKT08n)
	minKT10n, maxKT10n = min(R2oKT10n), max(R2oKT10n)

	plotmin = min(array([minKT00n, minKT02n, minKT04n, minKT06n, minKT08n, minKT10n]))
	plotmax = max(array([maxKT00n, maxKT02n, maxKT04n, maxKT06n, maxKT08n, maxKT10n]))

	#print plotmin, plotmax

	nbins = 20
	dKT = (plotmax - plotmin)/nbins

	plotbins = linspace(plotmin, plotmax, num = nbins)

	#bin_centers = 0.5*(bin_edge[1:] + bin_edges[:-1])
	colormap = plt.cm.gist_ncar
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, numplots)])

	# the histogram of the data
	nKT00n, binsKT00n, patchesKT00n = plt.hist(R2oKT00n, plotbins, histtype='step', normed=1)
	nKT02n, binsKT02n, patchesKT02n = plt.hist(R2oKT02n, plotbins, histtype='step', normed=1)
	nKT04n, binsKT04n, patchesKT04n = plt.hist(R2oKT04n, plotbins, histtype='step', normed=1)
	nKT06n, binsKT06n, patchesKT06n = plt.hist(R2oKT06n, plotbins, histtype='step', normed=1)
	nKT08n, binsKT08n, patchesKT08n = plt.hist(R2oKT08n, plotbins, histtype='step', normed=1)
	nKT10n, binsKT10n, patchesKT10n = plt.hist(R2oKT10n, plotbins, histtype='step', normed=1)

	plt.plot(delete(binsKT00n,-1) + 0.5*dKT, nKT00n, '-o', label='$K_T = 0.0$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT02n,-1) + 0.5*dKT, nKT02n, '-o', label='$K_T = 0.2$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT04n,-1) + 0.5*dKT, nKT04n, '-o', label='$K_T = 0.4$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT06n,-1) + 0.5*dKT, nKT06n, '-o', label='$K_T = 0.6$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT08n,-1) + 0.5*dKT, nKT08n, '-o', label='$K_T = 0.8$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT10n,-1) + 0.5*dKT, nKT10n, '-o', label='$K_T = 1.0$ GeV' % {"df": dfstem})

	#adds some nice bars in the background of the histogram to give plot some structure
	dummy1, dummy2, dummy3 = plt.hist(R2oKT00n, plotbins, histtype='bar', normed=1, edgecolor='black', color='white', alpha=0.5)

	plt.xlabel(r'$R^2_{o}/\left<R^2_{o}\!\right>_{\mathrm{ev}}$', {'fontsize': plotfontsize + 5})
	plt.ylabel(r'P($R^2_{o}/\left<R^2_{o}\!\right>_{\mathrm{ev}}$)', {'fontsize': plotfontsize + 5})
	plt.xticks(color='k', size=plotfontsize)
	plt.yticks(color='k', size=plotfontsize)
	minorLocator=MultipleLocator(0.05)
	ax.xaxis.set_minor_locator(minorLocator)
	ax.xaxis.set_tick_params(length=7.5, width=1.5)
	ax.xaxis.set_tick_params(which='minor',length=5.0, width=1.0)
	plt.xlim([plotmin, plotmax])
	plt.ylim([-0.25, 12.5])
	#plt.text(0.65, 4.5, 'Normalized event-by-event\n\t\tdistribution for $R^2_{o}$', fontsize=plotfontsize)
	#hl = plt.legend(loc=2, fontsize=plotfontsize)
	#hl.draw_frame(False)

	plt.legend()

	#plt.show()

	filename = '/home/plumberg.1/EBE-results/database/EbE%(ms)s_R2o%(n)i_vs_KT_1000evs%(df)s.pdf' % {"ms": methodStem, "n": chosenOrder, "df": dfstem}
	print 'Saving to', filename
	plt.savefig(filename, format='pdf')
	#plt.savefig('/home/plumberg.1/EBE-results/database/EbE_R2o_vs_KT_1000evs%(df)s.eps', format='eps')

	plt.close()


def plot_EbE_R2s(ebs, dfstem, radiiForEvents):
	numplots = 6
	plotfontsize = 18
	ax=plt.axes([0.1,0.15,0.85,0.8])
	plt.axhline(0.0, color='black', linewidth=1.5)
	
	direc = 'RESULTS_etaBYs_%(ebsstring)s/' % {"ebsstring": ebs}

	R2sKT00n = radiiForEvents[:,0] / mean(radiiForEvents[:,0])
	R2sKT02n = radiiForEvents[:,20] / mean(radiiForEvents[:,20])
	R2sKT04n = radiiForEvents[:,40] / mean(radiiForEvents[:,40])
	R2sKT06n = radiiForEvents[:,60] / mean(radiiForEvents[:,60])
	R2sKT08n = radiiForEvents[:,80] / mean(radiiForEvents[:,80])
	R2sKT10n = radiiForEvents[:,100] / mean(radiiForEvents[:,100])

	minKT00n, maxKT00n = min(R2sKT00n), max(R2sKT00n)
	minKT02n, maxKT02n = min(R2sKT02n), max(R2sKT02n)
	minKT04n, maxKT04n = min(R2sKT04n), max(R2sKT04n)
	minKT06n, maxKT06n = min(R2sKT06n), max(R2sKT06n)
	minKT08n, maxKT08n = min(R2sKT08n), max(R2sKT08n)
	minKT10n, maxKT10n = min(R2sKT10n), max(R2sKT10n)

	plotmin = min(array([minKT00n, minKT02n, minKT04n, minKT06n, minKT08n, minKT10n]))
	plotmax = max(array([maxKT00n, maxKT02n, maxKT04n, maxKT06n, maxKT08n, maxKT10n]))

	nbins = 20
	dKT = (plotmax - plotmin)/nbins

	plotbins = linspace(plotmin, plotmax, num=nbins)

	#bin_centers = 0.5*(bin_edge[1:] + bin_edges[:-1])
	colormap = plt.cm.gist_ncar
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, numplots)])

	# the histogram of the data
	nKT00n, binsKT00n, patchesKT00n = plt.hist(R2sKT00n, plotbins, histtype='step', normed=1)
	nKT02n, binsKT02n, patchesKT02n = plt.hist(R2sKT02n, plotbins, histtype='step', normed=1)
	nKT04n, binsKT04n, patchesKT04n = plt.hist(R2sKT04n, plotbins, histtype='step', normed=1)
	nKT06n, binsKT06n, patchesKT06n = plt.hist(R2sKT06n, plotbins, histtype='step', normed=1)
	nKT08n, binsKT08n, patchesKT08n = plt.hist(R2sKT08n, plotbins, histtype='step', normed=1)
	nKT10n, binsKT10n, patchesKT10n = plt.hist(R2sKT10n, plotbins, histtype='step', normed=1)

	plt.plot(delete(binsKT00n,-1) + 0.5*dKT, nKT00n, '-o', label='$K_T = 0.0$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT02n,-1) + 0.5*dKT, nKT02n, '-o', label='$K_T = 0.2$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT04n,-1) + 0.5*dKT, nKT04n, '-o', label='$K_T = 0.4$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT06n,-1) + 0.5*dKT, nKT06n, '-o', label='$K_T = 0.6$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT08n,-1) + 0.5*dKT, nKT08n, '-o', label='$K_T = 0.8$ GeV' % {"df": dfstem})
	plt.plot(delete(binsKT10n,-1) + 0.5*dKT, nKT10n, '-o', label='$K_T = 1.0$ GeV' % {"df": dfstem})

	#adds some nice bars in the background of the histogram to give plot some structure
	dummy1, dummy2, dummy3 = plt.hist(R2sKT00n, plotbins, histtype='bar', normed=1, edgecolor='black', color='white', alpha=0.5)

	plt.xlabel(r'$R^2_{s}/\left<R^2_{s}\!\right>_{\mathrm{ev}}$', {'fontsize': plotfontsize + 5})
	plt.ylabel(r'P($R^2_{s}/\left<R^2_{s}\!\right>_{\mathrm{ev}}$)', {'fontsize': plotfontsize + 5})
	plt.xticks(color='k', size=plotfontsize)
	plt.yticks(color='k', size=plotfontsize)
	minorLocator=MultipleLocator(0.05)
	ax.xaxis.set_minor_locator(minorLocator)
	ax.xaxis.set_tick_params(length=7.5, width=1.5)
	ax.xaxis.set_tick_params(which='minor',length=5.0, width=1.0)
	plt.xlim([plotmin, plotmax])
	#plt.xlim([0.85, 1.25])
	plt.ylim([-0.25, 10.5])
	#plt.text(0.79, 5.35, 'Normalized event-by-event\n\t\tdistribution for $R^2_{s}$', fontsize=plotfontsize)
	#hl = plt.legend(loc=2, fontsize=plotfontsize)
	#hl.draw_frame(False)

	plt.legend()

	#plt.show()

	filename = '/home/plumberg.1/EBE-results/database/EbE%(ms)s_R2s%(n)i_vs_KT_1000evs%(df)s.pdf' % {"ms": methodStem, "n": chosenOrder, "df": dfstem}
	print 'Saving to', filename
	plt.savefig(filename, format='pdf')
	#plt.savefig('/home/plumberg.1/EBE-results/database/EbE_R2s_vs_KT_1000evs%(df)s.eps', format='eps')

	plt.close()




def plot_EbE_R2l(ebs, dfstem, radiiForEvents):
	numplots = 6
	plotfontsize = 18
	ax=plt.axes([0.1,0.15,0.85,0.8])
	plt.axhline(0.0, color='black', linewidth=1.5)
	
	direc = 'RESULTS_etaBYs_%(ebsstring)s/' % {"ebsstring": ebs}

	R2lKT00n = radiiForEvents[:,0] / mean(radiiForEvents[:,0])
	R2lKT02n = radiiForEvents[:,20] / mean(radiiForEvents[:,20])
	R2lKT04n = radiiForEvents[:,40] / mean(radiiForEvents[:,40])
	R2lKT06n = radiiForEvents[:,60] / mean(radiiForEvents[:,60])
	R2lKT08n = radiiForEvents[:,80] / mean(radiiForEvents[:,80])
	R2lKT10n = radiiForEvents[:,100] / mean(radiiForEvents[:,100])

	minKT00n, maxKT00n = min(R2lKT00n), max(R2lKT00n)
	minKT02n, maxKT02n = min(R2lKT02n), max(R2lKT02n)
	minKT04n, maxKT04n = min(R2lKT04n), max(R2lKT04n)
	minKT06n, maxKT06n = min(R2lKT06n), max(R2lKT06n)
	minKT08n, maxKT08n = min(R2lKT08n), max(R2lKT08n)
	minKT10n, maxKT10n = min(R2lKT10n), max(R2lKT10n)

	plotmin, plotmax = min(array([minKT00n, minKT02n, minKT04n, minKT06n, minKT08n, minKT10n])), max(array([maxKT00n, maxKT02n, maxKT04n, maxKT06n, maxKT08n, maxKT10n]))

	#print plotmin, plotmax

	nbins = 20
	dKT = (plotmax - plotmin)/nbins

	plotbins = linspace(plotmin, plotmax, num=nbins)

	#bin_centers = 0.5*(bin_edge[1:] + bin_edges[:-1])
	colormap = plt.cm.gist_ncar
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, numplots)])

	# the histogram of the data
	nKT00n, binsKT00n, patchesKT00n = plt.hist(R2lKT00n, plotbins, histtype='step', normed=1)
	nKT02n, binsKT02n, patchesKT02n = plt.hist(R2lKT02n, plotbins, histtype='step', normed=1)
	nKT04n, binsKT04n, patchesKT04n = plt.hist(R2lKT04n, plotbins, histtype='step', normed=1)
	nKT06n, binsKT06n, patchesKT06n = plt.hist(R2lKT06n, plotbins, histtype='step', normed=1)
	nKT08n, binsKT08n, patchesKT08n = plt.hist(R2lKT08n, plotbins, histtype='step', normed=1)
	nKT10n, binsKT10n, patchesKT10n = plt.hist(R2lKT10n, plotbins, histtype='step', normed=1)

	plt.plot(delete(binsKT00n,-1) + 0.5*dKT, nKT00n, '-o', label='$K_T = 0.0$ GeV')
	plt.plot(delete(binsKT02n,-1) + 0.5*dKT, nKT02n, '-o', label='$K_T = 0.2$ GeV')
	plt.plot(delete(binsKT04n,-1) + 0.5*dKT, nKT04n, '-o', label='$K_T = 0.4$ GeV')
	plt.plot(delete(binsKT06n,-1) + 0.5*dKT, nKT06n, '-o', label='$K_T = 0.6$ GeV')
	plt.plot(delete(binsKT08n,-1) + 0.5*dKT, nKT08n, '-o', label='$K_T = 0.8$ GeV')
	plt.plot(delete(binsKT10n,-1) + 0.5*dKT, nKT10n, '-o', label='$K_T = 1.0$ GeV')

	#adds some nice bars in the background of the histogram to give plot some structure
	dummy1, dummy2, dummy3 = plt.hist(R2lKT04n, plotbins, histtype='bar', normed=1, edgecolor='black', color='white', alpha=0.5)

	plt.xlabel(r'$R^2_{l}/\left<R^2_{l} \!\right>_{\mathrm{ev}}$', {'fontsize': plotfontsize + 5})
	plt.ylabel(r'P($R^2_{l}/\left<R^2_{l}\!\right>_{\mathrm{ev}}$)', {'fontsize': plotfontsize + 5})
	plt.xticks(color='k', size=plotfontsize)
	plt.yticks(color='k', size=plotfontsize)
	minorLocator=MultipleLocator(0.05)
	ax.xaxis.set_minor_locator(minorLocator)
	ax.xaxis.set_tick_params(length=7.5, width=1.5)
	ax.xaxis.set_tick_params(which='minor',length=5.0, width=1.0)
	plt.xlim([plotmin, plotmax])
	plt.ylim([-0.25, 12.5])
	#plt.text(0.85, 4.75, 'Normalized event-by-event\n\t\tdistribution for $R^2_{l}$', fontsize=plotfontsize)
	#hl = plt.legend(loc=2, fontsize=plotfontsize)
	#hl.draw_frame(False)

	plt.legend()

	#plt.show()

	filename = '/home/plumberg.1/EBE-results/database/EbE%(ms)s_R2l%(n)i_vs_KT_1000evs%(df)s.pdf' % {"ms": methodStem, "n": chosenOrder, "df": dfstem}
	print 'Saving to', filename
	plt.savefig(filename, format='pdf')
	#plt.savefig('/home/plumberg.1/EBE-results/database/EbE_R2l_vs_KT_1000evs%(df)s.eps', format='eps')

	plt.close()


#########################################################################################
#########################################################################################
#########################################################################################

def plot_EbE_R2ij_v2(ebs, radiiForEvents, radiiForEventsNODF, direction, ymin, ymax, chosenFigSize):
	#Set-up
	global panelCounter
	numplots = 6
	plotfontsize = 18
	fig = plt.figure(figsize=chosenFigSize)
	#fig = plt.figure()
	fig.subplots_adjust(wspace=0.0, hspace=0.0, bottom=0.15, top=0.95)
	ax = fig.add_subplot(121)
	plt.axhline(0.0, color='black', linewidth=1.5)
	
	direc = 'RESULTS_etaBYs_%(ebsstring)s/' % {"ebsstring": ebs}
	
	R2ijKT00n = radiiForEvents[:,0] / mean(radiiForEvents[:,0])
	R2ijKT02n = radiiForEvents[:,20] / mean(radiiForEvents[:,20])
	R2ijKT04n = radiiForEvents[:,40] / mean(radiiForEvents[:,40])
	R2ijKT06n = radiiForEvents[:,60] / mean(radiiForEvents[:,60])
	R2ijKT08n = radiiForEvents[:,80] / mean(radiiForEvents[:,80])
	R2ijKT10n = radiiForEvents[:,100] / mean(radiiForEvents[:,100])

	minKT00n, maxKT00n = min(R2ijKT00n), max(R2ijKT00n)
	minKT02n, maxKT02n = min(R2ijKT02n), max(R2ijKT02n)
	minKT04n, maxKT04n = min(R2ijKT04n), max(R2ijKT04n)
	minKT06n, maxKT06n = min(R2ijKT06n), max(R2ijKT06n)
	minKT08n, maxKT08n = min(R2ijKT08n), max(R2ijKT08n)
	minKT10n, maxKT10n = min(R2ijKT10n), max(R2ijKT10n)

	plotmin = min(array([minKT00n, minKT02n, minKT04n, minKT06n, minKT08n, minKT10n]))
	plotmax = 0.995*max(array([maxKT00n, maxKT02n, maxKT04n, maxKT06n, maxKT08n, maxKT10n]))

	#print plotmin, plotmax

	nbins = 20
	dKT = (plotmax - plotmin)/nbins

	plotbins = linspace(plotmin, plotmax, num = nbins)

	#bin_centers = 0.5*(bin_edge[1:] + bin_edges[:-1])
	colormap = plt.cm.gist_ncar
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, numplots)])

	# the histogram of the data
	nKT00n, binsKT00n, patchesKT00n = ax.hist(R2ijKT00n, plotbins, histtype='step', normed=1)
	nKT02n, binsKT02n, patchesKT02n = ax.hist(R2ijKT02n, plotbins, histtype='step', normed=1)
	nKT04n, binsKT04n, patchesKT04n = ax.hist(R2ijKT04n, plotbins, histtype='step', normed=1)
	nKT06n, binsKT06n, patchesKT06n = ax.hist(R2ijKT06n, plotbins, histtype='step', normed=1)
	nKT08n, binsKT08n, patchesKT08n = ax.hist(R2ijKT08n, plotbins, histtype='step', normed=1)
	nKT10n, binsKT10n, patchesKT10n = ax.hist(R2ijKT10n, plotbins, histtype='step', normed=1)

	l00, = ax.plot(delete(binsKT00n,-1) + 0.5*dKT, nKT00n, '-o', label='$K_T = 0.0$ GeV')
	l02, = ax.plot(delete(binsKT02n,-1) + 0.5*dKT, nKT02n, '-o', label='$K_T = 0.2$ GeV')
	l04, = ax.plot(delete(binsKT04n,-1) + 0.5*dKT, nKT04n, '-o', label='$K_T = 0.4$ GeV')
	l06, = ax.plot(delete(binsKT06n,-1) + 0.5*dKT, nKT06n, '-o', label='$K_T = 0.6$ GeV')
	l08, = ax.plot(delete(binsKT08n,-1) + 0.5*dKT, nKT08n, '-o', label='$K_T = 0.8$ GeV')
	l10, = ax.plot(delete(binsKT10n,-1) + 0.5*dKT, nKT10n, '-o', label='$K_T = 1.0$ GeV')

	#adds some nice bars in the background of the histogram to give plot some structure
	dummy1, dummy2, dummy3 = ax.hist(R2ijKT00n, plotbins, histtype='bar', normed=1, edgecolor='black', color='white', alpha=0.5)
	
	text(0.075, 0.9, panelLabels[0], ha='center', va='center', transform=ax.transAxes, fontsize=plotfontsize)
	
	ax.axis([plotmin, plotmax, ymin, ymax])
	ax.set_xlabel(r'$R^2_{%(ij)s}\!/\left<R^2_{%(ij)s}\!\right>_{\mathrm{ev}}$' % {"ij": direction}, {'fontsize': plotfontsize + 5})
	ax.set_ylabel(r'P($R^2_{%(ij)s}\!/\left<R^2_{%(ij)s}\!\right>_{\mathrm{ev}}$)' % {"ij": direction}, {'fontsize': plotfontsize + 5})
	plt.xticks(color='k', size=plotfontsize-3)
	plt.yticks(color='k', size=plotfontsize)
	minorLocator=MultipleLocator(0.05)
	ax.xaxis.set_minor_locator(minorLocator)
	ax.xaxis.set_tick_params(length=7.5, width=1.5)
	ax.xaxis.set_tick_params(which='minor',length=5.0, width=1.0)

	#second half
	ax2 = fig.add_subplot(122)
	plt.axhline(0.0, color='black', linewidth=1.5)
	
	direc = 'RESULTS_etaBYs_%(ebsstring)s/' % {"ebsstring": ebs}
	
	R2ijKT00n = radiiForEventsNODF[:,0] / mean(radiiForEventsNODF[:,0])
	R2ijKT02n = radiiForEventsNODF[:,20] / mean(radiiForEventsNODF[:,20])
	R2ijKT04n = radiiForEventsNODF[:,40] / mean(radiiForEventsNODF[:,40])
	R2ijKT06n = radiiForEventsNODF[:,60] / mean(radiiForEventsNODF[:,60])
	R2ijKT08n = radiiForEventsNODF[:,80] / mean(radiiForEventsNODF[:,80])
	R2ijKT10n = radiiForEventsNODF[:,100] / mean(radiiForEventsNODF[:,100])

	minKT00n, maxKT00n = min(R2ijKT00n), max(R2ijKT00n)
	minKT02n, maxKT02n = min(R2ijKT02n), max(R2ijKT02n)
	minKT04n, maxKT04n = min(R2ijKT04n), max(R2ijKT04n)
	minKT06n, maxKT06n = min(R2ijKT06n), max(R2ijKT06n)
	minKT08n, maxKT08n = min(R2ijKT08n), max(R2ijKT08n)
	minKT10n, maxKT10n = min(R2ijKT10n), max(R2ijKT10n)

	plotmin = min(array([minKT00n, minKT02n, minKT04n, minKT06n, minKT08n, minKT10n]))
	plotmax = 0.995*max(array([maxKT00n, maxKT02n, maxKT04n, maxKT06n, maxKT08n, maxKT10n]))

	#print plotmin, plotmax

	nbins = 20
	dKT = (plotmax - plotmin)/nbins

	plotbins = linspace(plotmin, plotmax, num = nbins)

	#bin_centers = 0.5*(bin_edge[1:] + bin_edges[:-1])
	colormap = plt.cm.gist_ncar
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, numplots)])

	# the histogram of the data
	nKT00n, binsKT00n, patchesKT00n = ax2.hist(R2ijKT00n, plotbins, histtype='step', normed=1)
	nKT02n, binsKT02n, patchesKT02n = ax2.hist(R2ijKT02n, plotbins, histtype='step', normed=1)
	nKT04n, binsKT04n, patchesKT04n = ax2.hist(R2ijKT04n, plotbins, histtype='step', normed=1)
	nKT06n, binsKT06n, patchesKT06n = ax2.hist(R2ijKT06n, plotbins, histtype='step', normed=1)
	nKT08n, binsKT08n, patchesKT08n = ax2.hist(R2ijKT08n, plotbins, histtype='step', normed=1)
	nKT10n, binsKT10n, patchesKT10n = ax2.hist(R2ijKT10n, plotbins, histtype='step', normed=1)

	ax2.plot(delete(binsKT00n,-1) + 0.5*dKT, nKT00n, '-o', label='$K_T = 0.01$ GeV')
	ax2.plot(delete(binsKT02n,-1) + 0.5*dKT, nKT02n, '-o', label='$K_T = 0.21$ GeV')
	ax2.plot(delete(binsKT04n,-1) + 0.5*dKT, nKT04n, '-o', label='$K_T = 0.41$ GeV')
	ax2.plot(delete(binsKT06n,-1) + 0.5*dKT, nKT06n, '-o', label='$K_T = 0.61$ GeV')
	ax2.plot(delete(binsKT08n,-1) + 0.5*dKT, nKT08n, '-o', label='$K_T = 0.81$ GeV')
	ax2.plot(delete(binsKT10n,-1) + 0.5*dKT, nKT10n, '-o', label='$K_T = 1.01$ GeV')

	#adds some nice bars in the background of the histogram to give plot some structure
	dummy1, dummy2, dummy3 = ax2.hist(R2ijKT00n, plotbins, histtype='bar', normed=1, edgecolor='black', color='white', alpha=0.5)
	
	text(0.9, 0.9, panelLabels[1], ha='center', va='center', transform=ax2.transAxes, fontsize=plotfontsize)
	
	ax2.axis([plotmin, plotmax, ymin, ymax])
	ax2.set_xlabel(r'$R^2_{%(ij)s}\!/\left<R^2_{%(ij)s}\!\right>_{\mathrm{ev}}$' % {"ij": direction}, {'fontsize': plotfontsize + 5})
	minorLocator=MultipleLocator(0.05)
	ax2.xaxis.set_minor_locator(minorLocator)
	ax2.xaxis.set_tick_params(length=7.5, width=1.5)
	ax2.xaxis.set_tick_params(which='minor',length=5.0, width=1.0)
	ax2.yaxis.set_ticks_position('right')
	plt.xticks(color='k', size=plotfontsize-3)
	plt.yticks(color='k', size=plotfontsize)

	handles = [l00, l02, l04, l06, l08, l10]
	labels = [r'$K_T = 0.01$ GeV', r'$K_T = 0.21$ GeV', r'$K_T = 0.41$ GeV', r'$K_T = 0.61$ GeV', r'$K_T = 0.81$ GeV', r'$K_T = 1.01$ GeV']
	fig.legend(handles, labels, bbox_to_anchor=(0.5,0.9))
	#plt.show()
	#plt.tight_layout()

	filename = '/home/plumberg.1/EBE-results/database/EbE%(ms)s_R2%(ij)s%(n)i_vs_KT_1000evs_both.pdf' % {"ij": direction, "ms": methodStem, "n": chosenOrder}
	print 'Saving to', filename
	plt.savefig(filename, format='pdf', bbox_inches='tight')

	plt.close()




if __name__ == "__main__":
	#for dfs in dfstems:
	#	data = zeros([numberOfEvents, nKT, nOrder, nRadii * nTrig + 2])
	#	currentDirectory = workingDirectory + subDirectoryStem
	#	for event in xrange(1,numberOfEvents+1):
	#		filename = currentDirectory + str(event) + '/HBTradii_SVWR_cfs_ev' + str(event) + dfs + '.dat'
	#		data[event-1] = reshape(loadtxt(filename, usecols=(1, 2, 3, 5, 9)), [nKT, nOrder, nRadii * nTrig + 2])
	#	for etaBYs in ebsvals:
	#		plot_EbE_R2s(etaBYs, dfs, data[:,:,chosenOrder,2])
	#		plot_EbE_R2o(etaBYs, dfs, data[:,:,chosenOrder,3])
	#		plot_EbE_R2l(etaBYs, dfs, data[:,:,chosenOrder,4])
	data = zeros([numberOfEvents, nKT, nOrder, nRadii * nTrig + 2])
	dataNODF = zeros([numberOfEvents, nKT, nOrder, nRadii * nTrig + 2])
	currentDirectory = workingDirectory + subDirectoryStem
	for event in xrange(1,numberOfEvents+1):
		filename = currentDirectory + str(event) + '/HBTradii_SVWR_cfs_ev' + str(event)  + '.dat'
		data[event-1] = reshape(loadtxt(filename, usecols=(1, 2, 3, 5, 9)), [nKT, nOrder, nRadii * nTrig + 2])
		filename = currentDirectory + str(event) + '/HBTradii_SVWR_cfs_ev' + str(event)  + '_no_df.dat'
		dataNODF[event-1] = reshape(loadtxt(filename, usecols=(1, 2, 3, 5, 9)), [nKT, nOrder, nRadii * nTrig + 2])

	for etaBYs in ebsvals:
		plot_EbE_R2ij_v2(etaBYs, data[:,:,chosenOrder,2], dataNODF[:,:,chosenOrder,2], 's', -0.25, 11.5, (14,6))
		plot_EbE_R2ij_v2(etaBYs, data[:,:,chosenOrder,3], dataNODF[:,:,chosenOrder,3], 'o', -0.25, 12.5, (14,6))
		plot_EbE_R2ij_v2(etaBYs, data[:,:,chosenOrder,4], dataNODF[:,:,chosenOrder,4], 'l', -0.25, 11.5, (14,6))



# End of file
