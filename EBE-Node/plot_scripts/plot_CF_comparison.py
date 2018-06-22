#!/usr/bin/env python
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import mpl_toolkits.axes_grid.inset_locator as inloc
from scipy import interpolate
from scipy.interpolate import griddata

mpl.rcParams['pdf.fonttype'] = 42

npT, npphi = 15, 36
nqpts = 51

qAxisOpts = ['X', 'Y', 'Z']
qAxesLC = ['x', 'y', 'z']
#projectionOpts = ['_unprojected', '']
#projectionChoice = '_unprojected'
#projectionChoice = ''
#resFracsOpts = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80']
#resFracsOpts = ['0.00', '0.25', '0.50', '0.60', '0.70', '0.75']
resFracsOpts = ['0.60']
qtPtsOpts = ['23']
pyPtsOpts = ['15', '21']
PSOpts = ['12']

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

GeVToMeV = 1000.0
#GeVToMeV = 1.0	#uses GeV instead of MeV
panelLabels = ['(a)', '(b)', '(c)']
panelCounter = 0

comparisonPath = './' \
					+ 'AXIS_%(opt1)s_qt%(opt4)s_pT15_pY%(opt5)s_PS%(opt6)s/' \
					+ 'RESFRAC_%(opt2)s/' \
					+ 'results/correlfunct3D%(opt3)s_Pion_+.dat'

#############################################################################
#############################################################################
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

#############################################################################
#############################################################################

def generate_plotdata(path, ipT, ipphi, cols):
	nCols = len(cols)
	#print 'cols = ', cols
	# shape data arrays appropriately
	if cols[0]==2:
		nqx, nqy, nqz = nqpts, 1, 1
	elif cols[0]==3:
		nqx, nqy, nqz = 1, nqpts, 1
	elif cols[0]==4:
		nqx, nqy, nqz = 1, 1, nqpts
	
	#try:
	data = loadtxt(path, usecols=tuple(cols)).reshape([npT, npphi, nqx, nqy, nqz, nCols])

	if nCols==2:
		datax = data[ipT, ipphi, :, 0, 0, 0]
		datay = data[ipT, ipphi, :, 0, 0, 1]

		if cols[0]==3:
			#print 'Made it here'
			datax = data[ipT, ipphi, 0, :, 0, 0]
			datay = data[ipT, ipphi, 0, :, 0, 1]
			#print datax, datay
		elif cols[0]==4:
			datax = data[ipT, ipphi, 0, 0, :, 0]
			datay = data[ipT, ipphi, 0, 0, :, 1]
	else:
		datax = data[ipT, ipphi, :, 0, 0, 0]
		datay = 0.0*data[ipT, ipphi, :, 0, 0, 1]   \
	            + 0.0*data[ipT, ipphi, :, 0, 0, 2] \
	            + data[ipT, ipphi, :, 0, 0, 3]

		if cols[0]==3:
			datax = data[ipT, ipphi, 0, :, 0, 0]
			datay = 0.0*data[ipT, ipphi, 0, :, 0, 1]   \
	                + 0.0*data[ipT, ipphi, 0, :, 0, 2] \
	                + data[ipT, ipphi, 0, :, 0, 3]
		elif cols[0]==4:
			datax = data[ipT, ipphi, 0, 0, :, 0]
			datay = 0.0*data[ipT, ipphi, 0, 0, :, 1]   \
	                + 0.0*data[ipT, ipphi, 0, 0, :, 2] \
	                + data[ipT, ipphi, 0, 0, :, 3]

	xlower, xupper = min(datax), max(datax)
	ylower, yupper = min(datay), max(datay)

	f = interpolate.interp1d(datax, datay, kind='linear')
	plotdatax = linspace(xlower, xupper, 1001)
	plotdatay = f(plotdatax)

	#plotdatax = datax
	#plotdatay = datay
	#except:
	#	xlower, xupper = -0.1, 0.1
	#	plotdatax = linspace(xlower, xupper, 1001)
	#	plotdatay = 0.0*plotdatax

	#return plotdatax, plotdatay, [xlower, xupper, ylower, yupper]
	#return plotdatax, plotdatay, [xlower, xupper, min(plotdatay), max(plotdatay)]
	#return plotdatax, plotdatay, [xlower, xupper, 0.0, 0.5]
	#return plotdatax, plotdatay, [xlower, xupper, 1.8, 2.05]
	return plotdatax, plotdatay, [xlower, xupper, 1.0, 2.0]
	#return plotdatax, plotdatay, [GeVToMeV * xlower, GeVToMeV * xupper, 1.0, 2.0]


#############################################################################
#############################################################################

def generate_resfrac_comparison(projectionChoice, ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []
	
	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=float(len(qtPtsOpts)*len(pyPtsOpts)*len(resFracsOpts)*len(PSOpts)))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	idx = 0
	for inqt in xrange(len(qtPtsOpts)):
		for inpy in xrange(len(pyPtsOpts)):
			for iRF in xrange(len(resFracsOpts)):
				for iPS in xrange(len(PSOpts)):
					cmpPaths = [comparisonPath % {'opt1': qAxisOpts[iAxis],\
												  'opt2': resFracsOpts[iRF],\
												  'opt3': projectionChoice,\
												  'opt4': qtPtsOpts[inqt],\
												  'opt5': pyPtsOpts[inpy],\
												  'opt6': PSOpts[iPS]} for iAxis in xrange(len(qAxisOpts))]
					#print cmpPaths
					for iAxis in xrange(len(qAxisOpts)):
						cols = [2+iAxis, 9]		# chooses q-axes, CF value from CF files
						#cols = [2+iAxis, 6, 7, 8]	# chooses q-axes, CF value from CF files
						colorVal = scalarMap.to_rgba(float(idx))
						#colorVal = qAxisColors[inqt]
						plotdatax, plotdatay, lims = generate_plotdata(cmpPaths[iAxis], ipT, ipphi, cols)
						axs[iAxis].plot(plotdatax, plotdatay, linestyle='-', color=colorVal, linewidth=1.5, label=resFracsOpts[iRF])
						axs[iAxis].axis(lims)
						axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})
					
					idx += 1
	
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	#labels = [r'60% of resonance decay $\pi^+$s', r'100% of resonance decay $\pi^+$s']
	#fig.legend(handles, labels, bbox_to_anchor=(0.7,0.9), ncol=2)
	axs[0].legend(loc=0, prop={'size': plotfontsize+5})
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	
	plt.show(block=False)
	#plt.savefig('CF_resfracs_comp_%(opt1)s_vs_%(opt2)s_0_0.pdf' % {'opt1': resFracs[ind1], 'opt2': resFracs[ind2]}, format='pdf', bbox_inches='tight')


#############################################################################
#############################################################################

def generate_all_plots():
	#projectionChoice = '_unprojected'
	#for ipT in xrange(6,8,1):
	#	generate_resfrac_comparison(projectionChoice, ipT, 0)
	projectionChoice = ''
	for ipT in xrange(6,8,1):
		generate_resfrac_comparison(projectionChoice, ipT, 0)
	pause()

#############################################################################
#############################################################################

if __name__ == "__main__":
	generate_all_plots()

	print 'Finished all.'

# End of file
