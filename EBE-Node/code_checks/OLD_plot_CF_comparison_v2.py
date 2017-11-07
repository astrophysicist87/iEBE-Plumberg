#!/usr/bin/env python
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import griddata

mpl.rcParams['pdf.fonttype'] = 42

hbarC = 0.197327053
loadSeparately = True

#[[R2o, R2os, R2ol],[R2os, R2s, R2sl],[R2ol, R2sl, R2l]]
HBTinds = [[1,2,5],[2,0,4],[5,4,3]]

#npT, npphi = 7, 11
npT, npphi = 15, 36

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

nCols = 2
comparisonPath = 'code_checks/' \
					+ '%(opt1)s_resonance_extrapolation/' \
					+ 'FOintegral_extrap_%(opt2)s/' \
					+ '%(opt3)s_grid/' \
					+ 'TESTPBS_%(opt4)s_%(opt5)s_results-avg-1/correlfunct3D_Pion_+.dat'

loadedData = zeros([npT, npphi, 51, 51, 51, nCols])

def generate_plotdata_v1(path, ipT, ipphi, cols):
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
	return plotdatax, plotdatay, [xlower, xupper, 1.0, 2.0]


def generate_plotdata(path, ipT, ipphi, cols, dims):
	# shape data arrays appropriately
	
	data = loadtxt(path, usecols=tuple(cols)).reshape(dims+[nCols])

	datax = data[ipT, ipphi, :, (dims[3]-1)/2, (dims[4]-1)/2, 0]
	datay = data[ipT, ipphi, :, (dims[3]-1)/2, (dims[4]-1)/2, 1]

	if cols[0]==3:
		datax = data[ipT, ipphi, (dims[2]-1)/2, :, (dims[4]-1)/2, 0]
		datay = data[ipT, ipphi, (dims[2]-1)/2, :, (dims[4]-1)/2, 1]
	elif cols[0]==4:
		datax = data[ipT, ipphi, (dims[2]-1)/2, (dims[3]-1)/2, :, 0]
		datay = data[ipT, ipphi, (dims[2]-1)/2, (dims[3]-1)/2, :, 1]
	
	xlower, xupper = min(datax), max(datax)
	ylower, yupper = min(datay), max(datay)
	
	f = interpolate.interp1d(datax, datay, kind='linear')
	plotdatax = linspace(xlower, xupper, 1001)
	plotdatay = f(plotdatax)

	return plotdatax, plotdatay, [xlower, xupper, 1.0, 2.0]


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
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = [comparisonPath % {'opt1': resExtrapOpts[0], 'opt2': FOthresholds[ind1], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[1]} for iAxis in xrange(3)]
	cmpPaths2 = [comparisonPath % {'opt1': resExtrapOpts[0], 'opt2': FOthresholds[ind2], 'opt3': gridSizes[0], 'opt4': qAxes[iAxis], 'opt5': resFracs[1]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 9]	# chooses q-axes, CF value from CF files
		#print 'Comparing', cmpPaths1[iAxis], 'and', cmpPaths2[iAxis]
		plotdata1x, plotdata1y, lims1 = generate_plotdata_v1(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		plotdata2x, plotdata2y, lims2 = generate_plotdata_v1(cmpPaths2[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		# update axes limits and axes labels
		axs[iAxis].axis(lims2)
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})

	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'$f_{\mathrm{thresh}}=0.9$', r'$f_{\mathrm{thresh}}=1.0$']
	fig.legend(handles, labels, bbox_to_anchor=(0.55,0.9), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.set_ylabel(r'$\left< R^2_s \!\right>$ (fm$^2\!$)', {'fontsize': plotfontsize + 5})
	#ax1.text(0.85, 0.85,'(a)', transform=ax1.transAxes, fontsize=plotfontsize + 5)
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	plt.show()
	#plt.savefig('CF_FOthreshholds_comp_%(opt1)s_vs_%(opt2)s_0_0.pdf' % {'opt1': FOthresholds[ind1], 'opt2': FOthresholds[ind2]}, format='pdf', bbox_inches='tight')


def generate_resExtrapOpts_comparison(ind1, ind2, ipT, ipphi):
	# set-up
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
		axs[iAxis].plot(plotdata3x, plotdata3y, linestyle='-', color='black', linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims3)

		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=2.0, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=2.0, label=qAxes[iAxis])
		axs[iAxis].axis(lims2)

		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})

	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'With resonance extrapolation', r'Without resonance extrapolation']
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
	plotfontsize = 12
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	handles = []

	cmpPaths1 = ['code_checks/flesh_out_tests/sparse_grid_XYZ_results-avg-1/correlfunct1D_%(ax)s_Pion_+_fleshed_out.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]
	cmpPaths2 = ['code_checks/flesh_out_tests/dense_grid_%(ax)s_results-avg-1/correlfunct3D_Pion_+.dat' % {'ax': qAxes[iAxis]} for iAxis in xrange(3)]
	
	for iAxis in xrange(3):
		cols = [2+iAxis, 11]	# chooses q-axes, CF value from CF files
		#print 'Comparing', cmpPaths1[iAxis], 'and', cmpPaths2[iAxis]
		plotdata1x, plotdata1y, lims1 = generate_plotdata(cmpPaths1[iAxis], ipT, ipphi, cols)
		axs[iAxis].plot(plotdata1x, plotdata1y, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].axis(lims1)
		
		cols2 = [2+iAxis, 9]
		plotdata2x, plotdata2y, lims2 = generate_plotdata(cmpPaths2[iAxis], ipT, ipphi, cols2)
		axs[iAxis].plot(plotdata2x, plotdata2y, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		# update axes limits and axes labels
		axs[iAxis].axis(lims2)
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize + 5})

	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='-', color='black', linewidth=1.5)
	handles += [h]
	h, = axs[0].plot(linspace(-1,1,11), linspace(-100,-100,11), linestyle='--', color='black', linewidth=1.5)
	handles += [h]
	labels = [r'Sparse grid - fleshed out', r'Dense grid']
	fig.legend(handles, labels, bbox_to_anchor=(0.5,0.9), ncol=2)
	axs[0].set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	#ax1.set_xticklabels([])
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_flesh_out_0_0.pdf', format='pdf', bbox_inches='tight')


def getFitToCF(qpoints, pANDradii, axis):
	args = zeros(len(qpoints))
	localpT, localpphi = pANDradii[0:2]
	qipts = zeros([3,len(qpoints)])	#osl
	
	if axis==0:	#x-axis
		qipts[0] = cos(localpphi)*qpoints
		qipts[1] = sin(localpphi)*qpoints
	elif axis==1:	#y-axis
		qipts[0] = -sin(localpphi)*qpoints
		qipts[1] = cos(localpphi)*qpoints
	elif axis==2:	#z-axis
		qipts[2] = qpoints
	
	for d1 in xrange(3):
		for d2 in xrange(3):
			args += qipts[d1] * qipts[d2] * pANDradii[HBTinds[d1][d2]+2]
			# "+2" is shift due to pT and pphi being in pANDradii array
	
	args *= -1./(hbarC**2)	#get the units and sign right
	return 1.+exp(args)


def generate_PEA_plot(ipT, ipphi, includeFitCurve):
	# set-up
	global panelCounter
	plotfontsize = 24
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)        

	GeVToMeV = 1000.0

	fitData = loadtxt('HBTradii_GF_evavg_grid0.dat').reshape([npT, npphi, 8])
	lambdas = loadtxt('lambdas.dat').reshape([npT, npphi])
	
	print 'Loading data file...'
	#loadedData = loadtxt('correlfunct3D_Pion_+_fleshed_out.dat', usecols=(2,3,4,11)).reshape([npT, npphi, 51, 51, 51, 4])
	loadedData = loadtxt('correlfunct3D_Pion_+_fleshed_out_%(ipt)s_%(ipphi)s.dat' % {'ipt': ipT, 'ipphi':ipphi}, usecols=(2,3,4,11)).reshape([1, 1, 51, 51, 51, 4])
	print 'Data file loaded.'
	#print loadedData.shape, loadedData[0,0].shape
	
	#dims = [npT, npphi, nqx, nqy, nqz]
	#dims = [7, 11, 51, 51, 51]
	
	for iAxis in xrange(3):
		#print 'Doing iAxis =', iAxis
		#print loadedData.shape, loadedData[ipT,ipphi].shape
		#cols = [2+iAxis, 11]	# chooses q-axes, CF value from CF files
		#plotdatax, plotdatay, lims = generate_plotdata(loadedData[ipT,ipphi], iAxis)
		plotdatax, plotdatay, lims = generate_plotdata(loadedData[0,0], iAxis)
		# force x, y, and z axes to have same scale
		#print lims
		lims[0] = -37.5
		lims[1] = 37.5	
		if includeFitCurve:
			plotfity = 1.0 + lambdas[ipT,ipphi]*(getFitToCF(plotdatax, fitData[ipT, ipphi], iAxis) - 1.0)
		axs[iAxis].plot(GeVToMeV * plotdatax, plotdatay, linestyle='-', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		axs[iAxis].plot(GeVToMeV * plotdatax, plotfity, linestyle='--', color=qAxisColors[iAxis], linewidth=1.5, label=qAxes[iAxis])
		#axs[iAxis].plot(GeVToMeV * plotdatax, plotdatay, linestyle='-', color=qAxisColors[iAxis], linewidth=2.5, label=qAxes[iAxis])
		#axs[iAxis].plot(GeVToMeV * plotdatax, 1.0 + plotTdatay / Tspectra, linestyle='--', color=qAxisColors[iAxis], linewidth=2.5, label=qAxes[iAxis])
		axs[iAxis].set_xlabel(r'$q_%(opt1)s$ (MeV)' % {'opt1': qAxesLC[iAxis]}, {'fontsize': plotfontsize})
		axs[iAxis].set_xticks(arange(-40, 40, 20))
		axs[iAxis].set_xticklabels(['','-20','0','20',''],fontsize=plotfontsize-8)
		axs[iAxis].set_yticks(arange(1.0, 2.0+1.e-10, 0.2))
		axs[iAxis].set_yticklabels(ndarray.tolist(linspace(1.0,2.0,6)),fontsize=plotfontsize-8)
		axs[iAxis].axis(lims)
		text(0.075, 0.925, panelLabels[panelCounter], ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=plotfontsize)
		panelCounter += 1
		ptpphiString=r'$K_T = %(pt)0.1f$ GeV, $\Phi_K = %(pphi)0.0f$' % {'pt': fitData[ipT, ipphi, 0], 'pphi': fitData[ipT, ipphi, 1]}
		if iAxis==1:
			text(0.5, 0.15, ptpphiString, ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=plotfontsize)

		
	axs[0].set_ylabel(r'$C_{\mathrm{avg}}$', {'fontsize': plotfontsize + 5})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	axs[2].set_yticklabels(ndarray.tolist(linspace(1.0,2.0,6)),fontsize=plotfontsize-8)
	#ax1.set_xticklabels([])
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CF_fit_and_PEA_%(ipT)d_%(ipphi)d_EDIT.pdf' % {'ipT': ipT, 'ipphi': ipphi}, format='pdf', bbox_inches='tight')
	print 'Saved to', 'CF_fit_and_PEA_%(ipT)d_%(ipphi)d_EDIT.pdf' % {'ipT': ipT, 'ipphi': ipphi}



def generate_all_plots():
	chosenpT, chosenpphi = 0, 0
	#generate_resfrac_comparison(1, 4, chosenpT, chosenpphi)
	generate_FOthreshhold_comparison(0, 1, chosenpT, chosenpphi)
	#generate_resExtrapOpts_comparison(0, 1, chosenpT, chosenpphi)
	#generate_qt_spacing_comparison(3, 2)
	#generate_flesh_out_comparison(0, 0)
	#generate_PEA_plot(0, 0, True)
	#generate_PEA_plot(5, 0, True)
	#generate_PEA_plot(8, 0, True)

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
