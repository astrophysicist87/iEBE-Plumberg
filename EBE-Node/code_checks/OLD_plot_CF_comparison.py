#!/usr/bin/env python
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid.inset_locator as inloc
from scipy import interpolate
from scipy.interpolate import griddata

mpl.rcParams['pdf.fonttype'] = 42

npT, npphi = 7, 11

dfStems = ['', '_no_df']
qAxes = ['X', 'Y', 'Z']
qAxesLC = ['x', 'y', 'z']
FOthresholds = ['0.90', '1.00']
resExtrapOpts = ['with', 'without']
resFracs = ['0.00', '0.60', '0.80', '0.90', '1.00']
gridSizes = ['small', 'large']

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

GeVToMeV = 1000.0
#GeVToMeV = 1.0	#uses GeV instead of MeV
panelLabels = ['(a)', '(b)', '(c)']
panelCounter = 0

nCols = 2
comparisonPath = 'OLD_code_checks/' \
					+ '%(opt1)s_resonance_extrapolation/' \
					+ 'FOintegral_extrap_%(opt2)s/' \
					+ '%(opt3)s_grid/' \
					+ 'TESTPBS_%(opt4)s_%(opt5)s_results-avg-1/correlfunct3D_Pion_+.dat'

def generate_plotdata(path, ipT, ipphi, cols):
	# shape data arrays appropriately
	if cols[0]==2:
		nqx, nqy, nqz = 51, 1, 1
	elif cols[0]==3:
		nqx, nqy, nqz = 1, 51, 1
	elif cols[0]==4:
		nqx, nqy, nqz = 1, 1, 51
	
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
	return plotdatax, plotdatay, [GeVToMeV * xlower, GeVToMeV * xupper, 1.0, 2.0]


def generate_resfrac_comparison(ind1, ind2, ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = [comparisonPath % {'opt1': resExtrapOpts[0], 'opt2': FOthresholds[0], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[ind1]} for iAxis in xrange(3)]
	cmpPaths2 = [comparisonPath % {'opt1': resExtrapOpts[0], 'opt2': FOthresholds[0], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[ind2]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 9]	# chooses q-axes, CF value from CF files
		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims2)
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})
	
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'60% of resonance decay $\pi^+$s', r'100% of resonance decay $\pi^+$s']
	fig.legend(handles, labels, bbox_to_anchor=(0.7,0.9), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.set_ylabel(r'$\left< R^2_s \!\right>$ (fm$^2\!$)', {'fontsize': plotfontsize + 5})
	#ax1.text(0.85, 0.85,'(a)', transform=ax1.transAxes, fontsize=plotfontsize + 5)
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_resfracs_comp_%(opt1)s_vs_%(opt2)s_0_0.pdf' % {'opt1': resFracs[ind1], 'opt2': resFracs[ind2]}, format='pdf', bbox_inches='tight')



def generate_FOthreshhold_comparison(ind1, ind2, ipT, ipphi):
	# set-up
	global panelCounter
	plotfontsize = 14
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = [comparisonPath % {'opt1': resExtrapOpts[0], 'opt2': FOthresholds[ind1], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[1]} for iAxis in xrange(3)]
	cmpPaths2 = [comparisonPath % {'opt1': resExtrapOpts[0], 'opt2': FOthresholds[ind2], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[1]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 9]	# chooses q-axes, CF value from CF files
		#print 'Comparing', cmpPaths1[iAxis], 'and', cmpPaths2[iAxis]
		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(GeVToMeV * plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		text(0.1, 0.5, panelLabels[panelCounter], ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=plotfontsize+5)
		panelCounter += 1
		
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(GeVToMeV * plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		# update axes limits and axes labels
		axs[iAxis].axis(lims2)
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (MeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})
		axs[iAxis].set_xticks(arange(-60, 61, 20))
		axs[iAxis].set_xticklabels(['-60','-40','-20','0','20','40','60'],fontsize=plotfontsize-2)
		axs[iAxis].set_yticks(arange(1.0, 2.0+1.e-10, 0.2))
		axs[iAxis].set_yticklabels(ndarray.tolist(linspace(1.0,2.0,6)),fontsize=plotfontsize-2)

	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'$f_{\mathrm{thresh}}=0.9$', r'$f_{\mathrm{thresh}}=1.0$']
	fig.legend(handles, labels, bbox_to_anchor=(0.55,0.9), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	axs[2].set_yticks(arange(1.0, 2.0+1.e-10, 0.2))
	axs[2].set_yticklabels(ndarray.tolist(linspace(1.0,2.0,6)),fontsize=plotfontsize-2)
	#ax1.set_xticklabels([])
	#ax1.set_ylabel(r'$\left< R^2_s \!\right>$ (fm$^2\!$)', {'fontsize': plotfontsize + 5})
	#ax1.text(0.85, 0.85,'(a)', transform=ax1.transAxes, fontsize=plotfontsize + 5)
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_FOthreshholds_comp_%(opt1)s_vs_%(opt2)s_0_0.pdf' % {'opt1': FOthresholds[ind1], 'opt2': FOthresholds[ind2]}, format='pdf', bbox_inches='tight')


def generate_resExtrapOpts_comparison(ind1, ind2, ipT, ipphi):
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

def generate_qt_spacing_comparison(ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = ['code_checks/qt_spacing/TESTqtSpacing0_%(ax)s_results-avg-1/correlfunct3D_Pion_+.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]
	cmpPaths2 = ['code_checks/qt_spacing/TESTqtSpacing1_%(ax)s_results-avg-1/correlfunct3D_Pion_+.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 9]	# chooses q-axes, CF value from CF files
		#print 'Comparing', cmpPaths1[iAxis], 'and', cmpPaths2[iAxis]
		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		# update axes limits and axes labels
		axs[iAxis].axis(lims2)
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})

	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'$N_{q_t}=43$', r'$N_{q_t}=75$']
	fig.legend(handles, labels, bbox_to_anchor=(0.5,0.3), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_qtspacing_Nqt_43_vs_75_3_2.pdf', format='pdf', bbox_inches='tight')


def generate_Cev_slices(ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	paths = ['code_checks/qt_spacing/TESTqtSpacing0_%(ax)s_results-avg-1/correlfunct3D_Pion_+.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]


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


def generate_flesh_out_comparison(ipT, ipphi):
	# set-up
	global panelCounter
	panelCounter = 0
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = ['OLD_code_checks/flesh_out_tests/sparse_grid_XYZ_results-avg-1/correlfunct1D_%(ax)s_Pion_+_fleshed_out.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]
	cmpPaths2 = ['OLD_code_checks/flesh_out_tests/dense_grid_%(ax)s_results-avg-1/correlfunct3D_Pion_+.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 11]	# chooses q-axes, CF value from CF files
		#print 'Comparing', cmpPaths1[iAxis], 'and', cmpPaths2[iAxis]
		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(GeVToMeV * plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		cols2 = [2+iAxis, 9]
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols2)
		axs[iAxis].plot(GeVToMeV * plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		# update axes limits and axes labels
		axs[iAxis].axis(lims2)
		text(0.1, 0.5, panelLabels[panelCounter], ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=plotfontsize+5)
		panelCounter += 1
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (MeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})

	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'Sparse grid, fleshed out', r'Dense grid']
	fig.legend(handles, labels, bbox_to_anchor=(0.5,0.9), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_flesh_out_0_0.pdf', format='pdf', bbox_inches='tight')


def generate_PEA_plot(ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)        
	path = 'correlfunct3D_Pion_+.dat'
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 9]	# chooses q-axes, CF value from CF files
		#print 'Comparing', cmpPaths1[iAxis], 'and', cmpPaths2[iAxis]
		plotdatax, plotdatay, lims = generate_plotdata(path, ipT, ipphi, cols)
		axs[iAxis].plot(plotdatax, plotdatay, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims)
		
	axs[0].set_ylabel(r'$C_{\mathrm{ev}}$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	plt.show()
	#plt.savefig('CF_PEA_%(ipT)d_%(ipphi)d.pdf' % {'ipT': ipT, 'ipphi': ipphi}, format='pdf', bbox_inches='tight')



def generate_all_plots():
	chosenpT, chosenpphi = 0, 0
	#generate_resfrac_comparison(1, 4, chosenpT, chosenpphi)
	#generate_FOthreshhold_comparison(0, 1, chosenpT, chosenpphi)
	generate_resExtrapOpts_comparison(0, 1, chosenpT, chosenpphi)
	#generate_qt_spacing_comparison(3, 2)
	#generate_flesh_out_comparison(chosenpT, chosenpphi)
	#generate_PEA_plot(0,0)

if __name__ == "__main__":
	generate_all_plots()
	#chosenKT, chosenKphi = 0.1, 0.0
	#smallNPT, smallNPPHI = 7, 11
	#largeNPT, largeNPPHI = 15, 31
	#nqx, nqy, nqz = 25, 1, 1
	#path1 = 'code_checks/with_resonance_extrapolation/FOintegral_extrap_0.90/small_grid/TESTPBS_X_0.60_results-avg-1/correlfunct3D_Pion_+.dat'
	#path1 = 'small_grid_results-avg-1/correlfunct3D_Pion_+.dat'
	#path2 = 'code_checks/with_resonance_extrapolation/FOintegral_extrap_0.90/large_grid/TESTPBS_X_0.60_results-avg-1/correlfunct3D_Pion_+.dat'
	#CFdata1 = loadtxt(path1, usecols=(0,1,2,3,4,9)).reshape([smallNPT, smallNPPHI, nqx, nqy, nqz, 6])
	#CFdata2 = loadtxt(path2, usecols=(0,1,2,3,4,9)).reshape([largeNPT, largeNPPHI, nqx, nqy, nqz, 6])
	#CFatK1 = zeros([nqx, nqy, nqz])
	#CFatK2 = zeros([nqx, nqy, nqz])
	#for iqx in xrange(1):
	#	for iqy in xrange(1):
	#		for iqz in xrange(1):
	#			CFatK1[iqx,iqy,iqz] = evaluate_correlation_function(CFdata1, iqx, iqy, iqz, chosenKT, chosenKphi)
	#			CFatK2[iqx,iqy,iqz] = evaluate_correlation_function(CFdata2, iqx, iqy, iqz, chosenKT, chosenKphi)

	#savetxt('CFvSMALL.dat', CFatK1)
	#savetxt('CFvLARGE.dat', CFatK2)
	print 'Finished all.'

# End of file
