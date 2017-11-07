#!/usr/bin/env python
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid.inset_locator as inloc
from scipy import interpolate
from scipy.interpolate import griddata

mpl.rcParams['pdf.fonttype'] = 42

npT, npphi, npY = 15, 36, 21
nqpts = 11

qAxisOpts = ['X', 'Y', 'Z']
qAxesLC = ['x', 'y', 'z']
projectionOpts = ['', '_unprojected']
#resFracsOpts = ['0.00', '0.10', '0.20', '0.60', '1.00']
resFracsOpts = ['0.00', '0.10', '0.20', '0.60']

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

GeVToMeV = 1000.0
#GeVToMeV = 1.0	#uses GeV instead of MeV
panelLabels = ['(a)', '(b)', '(c)']
panelCounter = 0

nCols = 2
comparisonPath = './' \
					+ 'AXIS_%(opt1)s/' \
					+ 'RESFRAC_%(opt2)s/' \
					+ 'results/correlfunct3D%(opt3)s_Pion_+.dat'

#############################################################################
#############################################################################
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

#############################################################################
#############################################################################

def generate_plotdata(path, ipT, ipphi, cols):
	# shape data arrays appropriately
	if cols[0]==2:
		nqx, nqy, nqz = nqpts, 1, 1
	elif cols[0]==3:
		nqx, nqy, nqz = 1, nqpts, 1
	elif cols[0]==4:
		nqx, nqy, nqz = 1, 1, nqpts
	
	data = loadtxt(path, usecols=tuple(cols)).reshape([npT, npphi, nqx, nqy, nqz, nCols])

	datax = data[ipT, ipphi, :, 0, 0, 0]
	datay = data[ipT, ipphi, :, 0, 0, 1]

	if cols[0]==3:
		datax = data[ipT, ipphi, 0, :, 0, 0]
		datay = data[ipT, ipphi, 0, :, 0, 1]
	elif cols[0]==4:
		datax = data[ipT, ipphi, 0, 0, :, 0]
		datay = data[ipT, ipphi, 0, 0, :, 1]
	
	xlower, xupper = min(datax), max(datax)
	ylower, yupper = min(datay), max(datay)
	
	f = interpolate.interp1d(datax, datay, kind='cubic')
	plotdatax = linspace(xlower, xupper, 1001)
	plotdatay = f(plotdatax)

	#plotdatax = datax
	#plotdatay = datay
	
	#return plotdatax, plotdatay, [xlower, xupper, ylower, yupper]
	return plotdatax, plotdatay, [xlower, xupper, 1.0, 2.0]
	#return plotdatax, plotdatay, [GeVToMeV * xlower, GeVToMeV * xupper, 1.0, 2.0]

#############################################################################
#############################################################################

def generate_resfrac_comparison(ind1, ind2, ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = [comparisonPath % {'opt1': qAxisOpts[iAxis], 'opt2': resFracsOpts[ind1], 'opt3': projectionOpts[1]} for iAxis in xrange(3)]
	cmpPaths2 = [comparisonPath % {'opt1': qAxisOpts[iAxis], 'opt2': resFracsOpts[ind2], 'opt3': projectionOpts[1]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 9]	# chooses q-axes, CF value from CF files
		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxisOpts[iAxis])
		axs[iAxis].axis(lims1)
		
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxisOpts[iAxis])
		axs[iAxis].axis(lims2)
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})
	
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'60% of resonance decay $\pi^+$s', r'100% of resonance decay $\pi^+$s']
	#fig.legend(handles, labels, bbox_to_anchor=(0.7,0.9), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	
	plt.show(block=False)
	#plt.savefig('CF_resfracs_comp_%(opt1)s_vs_%(opt2)s_0_0.pdf' % {'opt1': resFracs[ind1], 'opt2': resFracs[ind2]}, format='pdf', bbox_inches='tight')


#############################################################################
#############################################################################



'''def generate_resExtrapOpts_comparison(ind1, ind2, ipT, ipphi):
	# set-up
	global panelCounter
	panelCounter = 0
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = [comparisonPath % {'opt1': resExtrapOpts[ind1], 'opt2': FOthresholds[0], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[1]} for iAxis in xrange(3)]
	cmpPaths2 = [comparisonPath % {'opt1': resExtrapOpts[ind2], 'opt2': FOthresholds[0], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[1]} for iAxis in xrange(3)]
	cmpPaths3 = [comparisonPath % {'opt1': resExtrapOpts[ind1], 'opt2': FOthresholds[0], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[4]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 9]	# chooses q-axes, CF value from CF files
		#print 'Comparing', cmpPaths1[iAxis], 'and', cmpPaths2[iAxis]
		plotdata3x, plotdata3y, lims3 = generate_plotdata(cmpPaths3[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(GeVToMeV * plotdata3x, plotdata3y, linestyle='-', color='black', linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims3)

		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(GeVToMeV * plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=2.0, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(GeVToMeV * plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=2.0, label=qAxes[iAxis])
		axs[iAxis].axis(lims2)

		text(0.1, 0.5, panelLabels[panelCounter], ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=plotfontsize+5)
		panelCounter += 1
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (MeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})

	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'With resonance extrapolation', r'Without resonance extrapolation$\,$']	#need some extra whitespace
	fig.legend(handles, labels, bbox_to_anchor=(0.7,0.925), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.set_ylabel(r'$\left< R^2_s \!\right>$ (fm$^2\!$)', {'fontsize': plotfontsize + 5})
	#ax1.text(0.85, 0.85,'(a)', transform=ax1.transAxes, fontsize=plotfontsize + 5)
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_resExtrapOpts_comp_%(opt1)s_vs_%(opt2)s_0_0.pdf' % {'opt1': resExtrapOpts[ind1], 'opt2': resExtrapOpts[ind2]}, format='pdf', bbox_inches='tight')
'''


'''
def generate_Cev_slices(ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	paths = ['code_checks/qt_spacing/TESTqtSpacing0_%(ax)s_results-avg-1/correlfunct3D_Pion_+.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]
'''

'''
def evaluate_correlation_function(CFdata, iqx, iqy, iqz, pT, pphi):
	pTpts = CFdata[:,:,iqx, iqy, iqz, 0]
	pphipts = CFdata[:,:,iqx, iqy, iqz, 1]
	pTpts = reshape(pTpts, pTpts.size)
	pphipts = reshape(pphipts, pphipts.size)
	points = vstack((pTpts, pphipts)).T
	CFdataSlice = CFdata[:,:,iqx, iqy, iqz, 5]
	values = reshape(CFdataSlice, CFdataSlice.size)
	print points.shape
	print values.shape
	return griddata(points, values, array([[pT,pphi]]), method='nearest')
'''

#############################################################################
#############################################################################

def generate_all_plots():
	#chosenpT, chosenpphi = 0, 0
	generate_resfrac_comparison(0, 1, 0, 1)
	generate_resfrac_comparison(0, 1, 4, 1)
	generate_resfrac_comparison(0, 1, 8, 1)
	#generate_resExtrapOpts_comparison(0, 1, chosenpT, chosenpphi)
	pause()

#############################################################################
#############################################################################

if __name__ == "__main__":
	generate_all_plots()

	print 'Finished all.'

# End of file
