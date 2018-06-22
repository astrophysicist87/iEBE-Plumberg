#!/usr/bin/env python

from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy import interpolate
import sys, os

mpl.rcParams['pdf.fonttype'] = 42
labelsize = 15
mpl.rcParams['xtick.labelsize'] = labelsize
mpl.rcParams['ytick.labelsize'] = labelsize

#################################################################
# Parameters to run script with
#################################################################
#file information
#filename = sys.argv[1]

qAxisOpts = ['XYZ']
resFracsOpts = ['0.60']
qtPtsOpts = ['23', '51']
pyPtsOpts = ['15', '21']
PSOpts = ['12']

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

GeVToMeV = 1000.0
#GeVToMeV = 1.0	#uses GeV instead of MeV
panelLabels = ['(a)', '(b)', '(c)']
panelCounter = 0

comparisonPath = './' \
					+ 'AXIS_%(opt1)s_qt%(opt3)s_pT15_pY%(opt4)s_PS%(opt5)s/' \
					+ 'RESFRAC_%(opt2)s/' \
					+ 'results/HBTradii_GF_cfs.dat'

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def generate_Rij_plot():
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(3, 1, figsize=(15,15))
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0
	lw2 = 1.0
	maxOrder = 3
	nKT = 101
	chosenCols = [0,1,2,4,8]
	chosenOrder = 0
	nCols = len(chosenCols)
	dims = [nKT*(maxOrder+1), nCols]
	
	handles = []
	
	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=float(len(qtPtsOpts)*len(pyPtsOpts)*len(resFracsOpts)*len(PSOpts)))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	idx = 0
	for iAxis in xrange(len(qAxisOpts)):
		for inqt in xrange(len(qtPtsOpts)):
			for inpy in xrange(len(pyPtsOpts)):
				for iRF in xrange(len(resFracsOpts)):
					for iPS in xrange(len(PSOpts)):
						cmpPath = comparisonPath % {'opt1': qAxisOpts[iAxis],\
													  'opt2': resFracsOpts[iRF],\
													  'opt3': qtPtsOpts[inqt],\
													  'opt4': pyPtsOpts[inpy],\
													  'opt5': PSOpts[iPS]}
						inHBTdata = loadtxt(cmpPath, usecols=tuple(chosenCols)).reshape(dims)
						HBTdata = inHBTdata[where(inHBTdata[:,1] == chosenOrder)]

						#print cmpPath
						#print inHBTdata
						colorVal = scalarMap.to_rgba(float(idx))
						for iLocAxis in xrange(3):
							axs[iLocAxis].plot(HBTdata[:,0], HBTdata[:,2+iLocAxis], linestyle='-', color=colorVal, linewidth=1.5, label=qtPtsOpts[inqt]+','+pyPtsOpts[inpy])
							lims = [min(HBTdata[:,0]), max(HBTdata[:,0]), 0.0, 1.1*max(HBTdata[:,2+iLocAxis])]
							axs[iLocAxis].axis(lims)

						axs[2].set_xlabel(r'$K_T$ (GeV)', {'fontsize': plotfontsize + 5})
						axs[0].set_ylabel(r'$R^2_s$', {'fontsize': plotfontsize + 5})
						axs[1].set_ylabel(r'$R^2_o$', {'fontsize': plotfontsize + 5})
						axs[2].set_ylabel(r'$R^2_l$', {'fontsize': plotfontsize + 5})
					
						idx += 1
	
	#h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	#handles += [h]
	#h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	#handles += [h]
	#labels = [r'60% of resonance decay $\pi^+$s', r'100% of resonance decay $\pi^+$s']
	#fig.legend(handles, labels, bbox_to_anchor=(0.7,0.9), ncol=2)
	axs[2].legend(loc=0, prop={'size': plotfontsize+5})
	
	plt.show()


def generate_all_plots():
	generate_Rij_plot()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
