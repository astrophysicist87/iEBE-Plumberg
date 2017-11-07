#!/usr/bin/env python

from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.patches import Ellipse
import sys, os

mpl.rcParams['pdf.fonttype'] = 42

#################################################################
# Parameters to run script with
#################################################################
#file information
#filename = sys.argv[1]

#grid information
npT, npphi = 15, 36
#npT, npphi = 1, 1
nqx, nqy, nqz = 51, 51, 51
#nqx, nqy, nqz = 1, 51, 1
#nqx, nqy, nqz = 1, 1, 51
nCols = 6
chosenCols = [0, 1, 2, 3, 4, 9]
dims = [npT, npphi, nqx, nqy, nqz, nCols]

HBTinds = [[1,2,5],[2,0,4],[5,4,3]]
dfStems = ['', '_no_df']
qAxes = ['X', 'Y', 'Z']
qAxesLC = ['x', 'y', 'z']
FOthresholds = ['0.90', '1.00']
resExtrapOpts = ['with', 'without']
resFracs = ['0.00', '0.60', '0.80', '0.90', '1.00']
gridSizes = ['small', 'large']

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

panelLabels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
panelCounter = 0


hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

#data = loadtxt(filename, usecols=tuple(chosenCols)).reshape(dims)

def generate_plotdata(localdata, iAxis):
	# shape data arrays appropriately
	dims = localdata.shape[0:3]
	#print localdata.shape, dims
	
	datax = localdata[:, (dims[1]-1)/2, (dims[2]-1)/2, iAxis]
	datay = localdata[:, (dims[1]-1)/2, (dims[2]-1)/2, 3]
	
	if iAxis==1:
		datax = localdata[(dims[0]-1)/2, :, (dims[2]-1)/2, iAxis]
		datay = localdata[(dims[0]-1)/2, :, (dims[2]-1)/2, 3]
	elif iAxis==2:
		datax = localdata[(dims[0]-1)/2, (dims[1]-1)/2, :, iAxis]
		datay = localdata[(dims[0]-1)/2, (dims[1]-1)/2, :, 3]
	
	xlower, xupper = min(datax), max(datax)
	ylower, yupper = min(datay), max(datay)
	
	f = interpolate.interp1d(datax, datay, kind='linear')
	plotdatax = linspace(xlower, xupper, 1001)
	plotdatay = f(plotdatax)

	return plotdatax, plotdatay, [xlower, xupper, 1.0, 2.0]

#################################################################
#################################################################

def generate_CF_plot(ipT, ipphi, includeFitCurve):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)        

	fitData = loadtxt('HBTradii_GF_evavg_grid0.dat').reshape([npT, npphi, 8])
	lambdas = loadtxt('lambdas.dat').reshape([npT, npphi])
	
	#print 'Loading data file...'
	#loadedData = loadtxt('correlfunct3D_Pion_+_fleshed_out.dat', usecols=(2,3,4,11)).reshape([npT, npphi, nqx, nqy, nqz, 4])
	loadedData = loadtxt('correlfunct3D_Pion_+_fleshed_out_%(ipT)d_%(ipphi)d.dat' % {'ipT': ipT, 'ipphi': ipphi}, usecols=(2,3,4,11)).reshape([1, 1, nqx, nqy, nqz, 4])
	#print 'Data file loaded.'
	#print loadedData.shape, loadedData[0,0].shape
	
	dims = [npT, npphi, nqx, nqy, nqz]
	
	for iAxis in xrange(3):
		#print 'Doing iAxis =', iAxis
		#print loadedData.shape, loadedData[ipT,ipphi].shape
		#cols = [2+iAxis, 11]	# chooses q-axes, CF value from CF files
		plotdatax, plotdatay, lims = generate_plotdata(loadedData[0,0], iAxis)
		#print iAxis, lambdas[ipT,ipphi], fitData[ipT, ipphi]
		if includeFitCurve:
			plotfity = 1.0 + lambdas[ipT,ipphi]*(getFitToCF(plotdatax, fitData[ipT, ipphi], iAxis) - 1.0)
		axs[iAxis].plot(plotdatax, plotdatay, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].plot(plotdatax, plotfity, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})
		axs[iAxis].axis(lims)
		ptpphiString=r'$K_T = %(pt)0.2f$ GeV, $\Phi_K = %(pphi)0.1f$' % {'pt': fitData[ipT, ipphi, 0], 'pphi': fitData[ipT, ipphi, 1]}
		if iAxis==1:
			text(0.5, 0.15, ptpphiString, ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=24)
			#text(0.5, 0.25, r'$p_T = %(pt)0.2f$ GeV,' % {'pt': fitData[ipT, ipphi, 0]}, ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=24)
			#text(0.5, 0.15, r'$p_{\phi} = %(pphi)0.2f$' % {'pphi': fitData[ipT, ipphi, 1]}, ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=24)
		
	axs[0].set_ylabel(r'$C_{\mathrm{avg}}$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_fit_and_evavg_%(ipT)d_%(ipphi)d.pdf' % {'ipT': ipT, 'ipphi': ipphi}, format='pdf', bbox_inches='tight')
	print 'Saved to', 'CF_fit_and_evavg_%(ipT)d_%(ipphi)d.pdf' % {'ipT': ipT, 'ipphi': ipphi}


#################################################################
#################################################################

def generate_all_plots():
	generate_CF_plot(0, 0, True)

#################################################################
#################################################################

if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
