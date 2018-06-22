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


def eval_corr_fit_func(xpts,lam,R):
	return 1.0 + lam*exp(-xpts**2*R**2/hbarC**2)

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


def generate_lambdas_plot(whichIPTtoCircle):
	# set-up
	plotfontsize = 12
	#fig, ax = plt.subplots(1, 1, aspect='equal')
	fig = plt.figure()
	ax = fig.add_subplot(111, aspect='equal')
	#fig.subplots_adjust(wspace=0.0, hspace=0.0)        
	
	pTData = loadtxt('HBTradii_GF_evavg_grid0.dat', usecols=(0,1)).reshape([npT, npphi, 2])
	lambdas = loadtxt('lambdas.dat').reshape([npT, npphi])
	
	chosenIPPHI = 0
	
	xdata = pTData[:,0,0]
	ydata = lambdas[:,0]
	
	#xlower, xupper = min(xdata), max(xdata)
	#ylower, yupper = min(ydata), max(ydata)
	lims = [0.01, 1.01, 0.0, 1.0]
	
	f = interpolate.interp1d(xdata, ydata, kind='cubic')
	plotdatax = linspace(0.01, 1.01, 1001)
	plotdatay = f(plotdatax)
	ax.plot(plotdatax, plotdatay, linestyle='-', color='purple', linewidth=2.5)
	#ax.axhline(0.0, color='black', linewidth=1)
	ax.axis(lims)
	
	cxc, cyc = 0.5, 0.7
	if whichIPTtoCircle==0:
		cxc, cyc = 0.15, 0.7
	elif whichIPTtoCircle==5:
		cxc, cyc = 0.4, 0.8
	elif whichIPTtoCircle==8:
		cxc, cyc = 0.85, 0.85
	
	ax.set_xlabel(r'$K_T$ (GeV)', {'fontsize': plotfontsize + 5})
	ax.set_ylabel(r'$\lambda(K_T)$', {'fontsize': plotfontsize + 5})
	text(0.5, 0.15, r'$\Phi_K = %(pphi)0.0f$' % {'pphi': 0.0}, ha='center', va='center', transform=ax.transAxes, fontsize=24)
	
	#circle = plt.Circle((cxc, cyc), 0.15, color='r', fill=False, linewidth=3.0)
	#ax.add_artist(circle)
	
	#plt.show()
	plt.savefig('lambdas_evavg_%(ipT)d_%(ipphi)d_EDIT.pdf' % {'ipT': whichIPTtoCircle, 'ipphi': chosenIPPHI}, format='pdf', bbox_inches='tight')
	print 'Saved to', 'lambdas_evavg_%(ipT)d_%(ipphi)d_EDIT.pdf' % {'ipT': whichIPTtoCircle, 'ipphi': chosenIPPHI}



def generate_R2ij_plot(whichIPTtoCircle):
	# set-up
	plotfontsize = 12
	#fig, ax = plt.subplots(1, 1, aspect='equal')
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#fig.subplots_adjust(wspace=0.0, hspace=0.0)        
	
	data = loadtxt('HBTradii_GF_evavg_grid0.dat').reshape([npT, npphi, 8])
	pTData = data[:,:,[0,1]]
	R2sdata = data[:,:,2]
	R2odata = data[:,:,3]
	R2ldata = data[:,:,5]
	
	chosenIPPHI = 0
	
	xdata = pTData[:,0,0]
	ysdata = R2sdata[:,0]
	yodata = R2odata[:,0]
	yldata = R2ldata[:,0]
	
	xlower, xupper = min(xdata), max(xdata)
	ylower = min(hstack((ysdata, yodata, yldata)))
	yupper = max(hstack((ysdata, yodata, yldata)))
	lims = [0.01, 1.01, 0.0, yupper+1.0]
	
	fs = interpolate.interp1d(xdata, ysdata, kind='cubic')
	fo = interpolate.interp1d(xdata, yodata, kind='cubic')
	fl = interpolate.interp1d(xdata, yldata, kind='cubic')
	plotdatax = linspace(0.01, 1.01, 1001)
	plotdatays = fs(plotdatax)
	plotdatayo = fo(plotdatax)
	plotdatayl = fl(plotdatax)
	ax.plot(plotdatax, plotdatays, linestyle='-', color='red', linewidth=2.5, label=r'$R^2_s$')
	ax.plot(plotdatax, plotdatayo, linestyle='-', color='green', linewidth=2.5, label=r'$R^2_o$')
	ax.plot(plotdatax, plotdatayl, linestyle='-', color='blue', linewidth=2.5, label=r'$R^2_l$')
	ax.axis(lims)
	
	cxc, cyc = 0.5, 25.0
	if whichIPTtoCircle==0:
		cxc, cyc = 0.06, 35.0
	elif whichIPTtoCircle==5:
		cxc, cyc = 0.4, 17.5
	elif whichIPTtoCircle==8:
		cxc, cyc = 0.95, 15.0
	
	ax.set_xlabel(r'$K_T$ (GeV)', {'fontsize': plotfontsize + 5})
	ax.set_ylabel(r'$R^2_{ij}(K_T, p_{\phi} = 0.0)$ (fm$^2$)', {'fontsize': plotfontsize + 5})
	#text(0.5, 0.15, r'$p_{\phi} = %(pphi)0.1f$' % {'pphi': 0.0}, ha='center', va='center', transform=ax.transAxes, fontsize=24)
	
	#circle = plt.Circle((cxc, cyc), 0.15, color='r', fill=False, linewidth=3.0)
	#ax.add_artist(circle)
	#ellipse = Ellipse(xy=(cxc, cyc), width=0.1, height=30.0, angle=0.0, color='r', fill=False, linewidth=3.0)
	#ax.add_artist(ellipse)
	
	ax.legend(loc=0, prop={'size': plotfontsize+5})
	
	#plt.show()
	plt.savefig('R2ij_evavg_%(ipT)d_%(ipphi)d.pdf' % {'ipT': whichIPTtoCircle, 'ipphi': chosenIPPHI}, format='pdf', bbox_inches='tight')
	print 'Saved to', 'R2ij_evavg_%(ipT)d_%(ipphi)d.pdf' % {'ipT': whichIPTtoCircle, 'ipphi': chosenIPPHI}


def compare_thermal_and_resonance_CF(ipT, ipphi):
	# set-up
	global panelCounter
	plotfontsize = 24
	fig, axs = plt.subplots(1, 3, figsize=(15,5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0)       
	
	GeVToMeV = 1000.0

	fitData = loadtxt('HBTradii_GF_evavg_grid0.dat').reshape([npT, npphi, 8])
	lambdas = loadtxt('lambdas.dat').reshape([npT, npphi])
	
	#print 'Loading data file...'
	#loadedData = loadtxt('correlfunct3D_Pion_+_fleshed_out.dat', usecols=(2,3,4,11)).reshape([npT, npphi, nqx, nqy, nqz, 4])
	loadedData = loadtxt('correlfunct3D_Pion_+_fleshed_out_%(ipT)d_%(ipphi)d.dat' % {'ipT': ipT, 'ipphi': ipphi}, usecols=(2,3,4,11)).reshape([1, 1, nqx, nqy, nqz, 4])
	loadedThermalData = loadtxt('correlfunct3D_Pion_+_fleshed_out_%(ipT)d_%(ipphi)d.dat' % {'ipT': ipT, 'ipphi': ipphi}, usecols=(2,3,4,8)).reshape([1, 1, nqx, nqy, nqz, 4])
	#print 'Data file loaded.'
	#print loadedData.shape, loadedData[0,0].shape
	
	dims = [npT, npphi, nqx, nqy, nqz]
	
	for iAxis in xrange(3):
		#print 'Doing iAxis =', iAxis
		#print loadedData.shape, loadedData[ipT,ipphi].shape
		plotdatax, plotdatay, lims = generate_plotdata(loadedData[0,0], iAxis)
		plotdatax, plotTdatay, lims = generate_plotdata(loadedThermalData[0,0], iAxis)
		
		# force x, y, and z axes to have same scale
		#print lims
		lims[0] = -37.5
		lims[1] = 37.5	
	
		#print iAxis, lambdas[ipT,ipphi], fitData[ipT, ipphi]
		Tspectra = loadedThermalData[0,0,nqx/2,nqy/2,nqz/2,3]
		#print Tspectra
		#if includeFitCurve:
		#	plotfity = 1.0 + lambdas[ipT,ipphi]*(getFitToCF(plotdatax, fitData[ipT, ipphi], iAxis) - 1.0)
		axs[iAxis].plot(GeVToMeV * plotdatax, plotdatay, linestyle='-', color=qAxisColors[iAxis], linewidth=2.5, label=qAxes[iAxis])
		axs[iAxis].plot(GeVToMeV * plotdatax, 1.0 + plotTdatay / Tspectra, linestyle='--', color=qAxisColors[iAxis], linewidth=2.5, label=qAxes[iAxis])
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
			#text(0.5, 0.25, r'$p_T = %(pt)0.2f$ GeV,' % {'pt': fitData[ipT, ipphi, 0]}, ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=24)
			#text(0.5, 0.15, r'$p_{\phi} = %(pphi)0.2f$' % {'pphi': fitData[ipT, ipphi, 1]}, ha='center', va='center', transform=axs[iAxis].transAxes, fontsize=24)
		
	axs[0].set_ylabel(r'$C_{\mathrm{avg}}$', {'fontsize': plotfontsize})
	axs[1].set_yticklabels([])
	axs[2].yaxis.set_ticks_position('right')
	axs[2].set_yticklabels(ndarray.tolist(linspace(1.0,2.0,6)),fontsize=plotfontsize-8)
	#ax1.set_xticklabels([])
	#ax1.legend(loc=3, prop={'size': plotfontsize})
	
	#plt.show()
	plt.savefig('CFevavg_w_v_wo_res_%(ipT)d_%(ipphi)d_EDIT.pdf' % {'ipT': ipT, 'ipphi': ipphi}, format='pdf', bbox_inches='tight')
	print 'Saved to', 'CFevavg_w_v_wo_res_%(ipT)d_%(ipphi)d_EDIT.pdf' % {'ipT': ipT, 'ipphi': ipphi}



def extraPlot(ipT, ipphi):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0)        

	fitData = loadtxt('HBTradii_GF_evavg_grid0.dat').reshape([npT, npphi, 8])

	fitLargeQlambda, fitLargeQR = 0.583404254910389, 5.33478579266871
	fitSmallQlambda, fitSmallQR = 0.784140991401117, 8.05580378843792
	fitAllQlambda, fitAllQR = 0.737582423165663, 6.46478233610704
	
	#print 'Loading data file...'
	loadedData = loadtxt('correlfunct3D_Pion_+_fleshed_out_%(ipT)d_%(ipphi)d.dat' % {'ipT': ipT, 'ipphi': ipphi}, usecols=(2,3,4,11)).reshape([1, 1, nqx, nqy, nqz, 4])
	#print 'Data file loaded.'
	#print loadedData.shape, loadedData[0,0].shape
	
	dims = [npT, npphi, nqx, nqy, nqz]
	
	plotdatax, plotdatay, lims = generate_plotdata(loadedData[0,0], 0)
	axs.plot(plotdatax, plotdatay, linestyle='-', color='black', linewidth=1.5, label=r'$C_{\mathrm{avg}}$')
	axs.plot(plotdatax, eval_corr_fit_func(plotdatax,fitSmallQlambda,fitSmallQR), linestyle='-', color='red', linewidth=1.5, label=r'$|q_x| \leq$ 20 MeV')
	axs.plot(plotdatax, eval_corr_fit_func(plotdatax,fitLargeQlambda,fitLargeQR), linestyle='-', color='green', linewidth=1.5, label=r'$|q_x| \geq$ 20 MeV')
	axs.plot(plotdatax, eval_corr_fit_func(plotdatax,fitAllQlambda,fitAllQR), linestyle='-', color='blue', linewidth=1.5, label=r'All $q_x$')
	axs.set_xlabel(r'$q_%(opt1)s$ (GeV)' % {'opt1': qAxesLC[0]}, {'fontsize': plotfontsize + 5})
	axs.axis(lims)
	ptpphiString=r'$K_T = %(pt)0.2f$ GeV, $\Phi_K = %(pphi)0.1f$' % {'pt': fitData[ipT, ipphi, 0], 'pphi': fitData[ipT, ipphi, 1]}
		
	#axs.set_ylabel(r'$C_{\mathrm{avg}}$', {'fontsize': plotfontsize + 5})
	axs.legend(loc=0, prop={'size': plotfontsize+5})
	
	#plt.show()
	plt.savefig('CFevavg_fitrange_X_%(ipT)d_%(ipphi)d.pdf' % {'ipT': ipT, 'ipphi': ipphi}, format='pdf', bbox_inches='tight')
	print 'Saved to', 'CFevavg_fitrange_X_%(ipT)d_%(ipphi)d.pdf' % {'ipT': ipT, 'ipphi': ipphi}


def generate_all_plots():
	#generate_CF_plot(0, 0, True)
	#generate_CF_plot(5, 0, True)
	#generate_CF_plot(8, 0, True)
	#generate_lambdas_plot(0)
	#generate_lambdas_plot(5)
	#generate_lambdas_plot(8)
	#generate_R2ij_plot(0)
	#generate_R2ij_plot(5)
	#generate_R2ij_plot(8)
	#compare_thermal_and_resonance_CF(0, 0)
	#compare_thermal_and_resonance_CF(5, 0)
	#compare_thermal_and_resonance_CF(8, 0)
	extraPlot(0, 0)


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
