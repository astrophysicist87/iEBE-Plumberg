#!/usr/bin/env python

from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.patches import Ellipse
import sys, os

mpl.rcParams['pdf.fonttype'] = 42

### command-line arguments
thermalFileToPlot = sys.argv[1]
resonanceFileToPlot = sys.argv[2]

# other initial parameters
usingSmoothing = True

#########################################################################################
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")
#########################################################################################

def makeDEAPlot():
	# set-up
	plotfontsize = 12
	ms = 8
	fig, axs = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0)        
	
	thermalData = loadtxt(thermalFileToPlot)
	resonanceData = loadtxt(resonanceFileToPlot)
	
	axs.plot(thermalData[:,0], thermalData[:,1],\
				linestyle='--', marker='s', markevery=4, markersize=ms, color='red',\
				linewidth=1.5)
	axs.plot(thermalData[:,0], thermalData[:,2],\
				linestyle='--', marker='o', markevery=4, markersize=ms, color='green',\
				linewidth=1.5)
	axs.plot(thermalData[:,0], thermalData[:,3],\
				linestyle='--', marker='v', markevery=4, markersize=ms, color='blue',\
				linewidth=1.5)
	axs.plot(resonanceData[:,0], resonanceData[:,1],\
				linestyle='-', marker='s', markevery=4, markersize=ms, color='red',\
				linewidth=1.5, label=r'$\left< R^2_s \!\right>$')
	axs.plot(resonanceData[:,0], resonanceData[:,2],\
				linestyle='-', marker='o', markevery=4, markersize=ms, color='green',\
				linewidth=1.5, label=r'$\left< R^2_o \!\right>$')
	axs.plot(resonanceData[:,0], resonanceData[:,3],\
				linestyle='-', marker='v', markevery=4, markersize=ms, color='blue',\
				linewidth=1.5, label=r'$\left< R^2_l \!\right>$')
	
	axs.set_xlim(0.001, 1.0)
	axs.set_ylim(bottom=0.0)
	axs.set_xlabel(r'$K_T$ (GeV)' % {'fontsize': plotfontsize + 5})
	axs.set_ylabel(r'$\left< R^2_{ij}\right>$ (' + 'fm' + r'$^2\!)$', {'fontsize': plotfontsize + 5})
	axs.legend(loc=0, prop={'size': plotfontsize+5})
	axs.text(0.925, 0.15, '(a)', ha='center', va='center', transform=axs.transAxes, fontsize=plotfontsize + 10)
	
	#plt.show(block=False)
	smoothingStem = ''
	if usingSmoothing:
		smoothingStem = '_smoothed'
	outputFileName = './HBT_GF_w_v_wo_res_DEA' + smoothingStem + '.pdf'
	plt.savefig(outputFileName, format='pdf', bbox_inches='tight')
	print 'Saved to', outputFileName


def makerelsigPlot():
	# set-up
	plotfontsize = 12
	ms = 8
	fig, axs = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0)        
	
	thermalData = loadtxt(thermalFileToPlot)
	resonanceData = loadtxt(resonanceFileToPlot)
	
	#axs.plot(thermalData[:,0], sqrt(thermalData[:,4])/thermalData[:,1],\
	#			linestyle='--', marker='s', markevery=4, markersize=ms, color='red', linewidth=1.5)
	#axs.plot(thermalData[:,0], sqrt(thermalData[:,5])/thermalData[:,2],\
	#			linestyle='--', marker='o', markevery=4, markersize=ms, color='green', linewidth=1.5)
	#axs.plot(thermalData[:,0], sqrt(thermalData[:,6])/thermalData[:,3],\
	#			linestyle='--', marker='v', markevery=4, markersize=ms, color='blue', linewidth=1.5)
	axs.plot(resonanceData[:,0], sqrt(resonanceData[:,4])/resonanceData[:,1],\
				linestyle='-', marker='s', markevery=4, markersize=ms, color='red',\
				linewidth=1.5, label=r'$\sigma_s/\left< R^2_s \!\right>$')
	axs.plot(resonanceData[:,0], sqrt(resonanceData[:,5])/resonanceData[:,2],\
				linestyle='-', marker='o', markevery=4, markersize=ms, color='green',\
				linewidth=1.5, label=r'$\sigma_o/\left< R^2_o \!\right>$')
	axs.plot(resonanceData[:,0], sqrt(resonanceData[:,6])/resonanceData[:,3],\
				linestyle='-', marker='v', markevery=4, markersize=ms, color='blue',\
				linewidth=1.5, label=r'$\sigma_l/\left< R^2_l \!\right>$')
	
	axs.set_xlim(0.001, 1.0)
	axs.set_ylim(bottom=0.0)
	axs.set_xlabel(r'$K_T$ (GeV)' % {'fontsize': plotfontsize + 5})
	axs.set_ylabel(r'$\sigma_{ij}/\left< R^2_{ij} \right>$', {'fontsize': plotfontsize + 5})
	axs.legend(loc=0, prop={'size': plotfontsize+5})
	axs.text(0.925, 0.075, '(b)', ha='center', va='center', transform=axs.transAxes, fontsize=plotfontsize + 10)
	
	#plt.show(block=False)
	smoothingStem = ''
	if usingSmoothing:
		smoothingStem = '_smoothed'
	outputFileName = './HBT_GF_w_res_relsig' + smoothingStem + '.pdf'
	plt.savefig(outputFileName, format='pdf', bbox_inches='tight')
	print 'Saved to', outputFileName



if __name__ == "__main__":
	makeDEAPlot()
	makerelsigPlot()
	
	#pause()
	


# End of file
