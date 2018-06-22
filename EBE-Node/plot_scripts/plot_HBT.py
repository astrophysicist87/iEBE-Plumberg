#!/usr/bin/env python

from numpy import *
from pylab import *
import matplotlib.pyplot as plt
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

#grid information
nKT = 4
nqpts = 13
nqaxes = 3
nCols = 5
dims = [nqaxes, nqpts, nKT, nCols]

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

panelLabels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
panelCounter = 0

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

qaxes = ['o','s','l']
	

def Cpts(qpts, fitRadius2, fitlambda):
	return 1.0 + fitlambda*exp(-qpts**2 *fitRadius2 / hbarC**2)


def generate_CF_plot():
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(2, 2, figsize=(11,8.5))
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0
	lw2 = 1.0
	nKT = 101
	chosenCols = [1,3,5,9]
	qaxisCols = [1,0,2]	#order of out-side-long columns in data file
	CS_qaxisCols = [0,2,1]	#order of out-side-long columns in data file
	
	#CFdata = loadtxt('../../interpolate_correlation_function/tmp.out').reshape(dims)
	#CFdata = loadtxt('../../interpolate_correlation_function/tmp_leadingBCs.out').reshape(dims)
	CFdata = loadtxt('../../interpolate_correlation_function/averaged_CF_slices.out').reshape(dims)
	pTpts, pTwts = loadtxt('pT_gauss_table.dat').T
	pphipts, pphiwts = loadtxt('phi_gauss_table.dat').T
	fitRadiidata = loadtxt('HBTradii_GF_cfs_evavg.dat', usecols=tuple(chosenCols)).reshape([nKT, len(chosenCols)])
	fitLambdadata = loadtxt('lambdas.dat').reshape([len(pTpts),len(pphipts)])

	CHUNSHEN_CFdata_KT100 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_correlation_function_KT_0_0.2.dat').reshape([31,31,31,8])
	CHUNSHEN_CFdata_KT300 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_correlation_function_KT_0.2_0.4.dat').reshape([31,31,31,8])
	CHUNSHEN_CFdata_KT500 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_correlation_function_KT_0.4_0.6.dat').reshape([31,31,31,8])
	CHUNSHEN_CFdata_KT700 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_correlation_function_KT_0.6_0.8.dat').reshape([31,31,31,8])
	
	CS_qaxes = [(slice(15,16),slice(0,31),slice(15,16),[0,6]), (slice(0,31), slice(15,16), slice(15,16), [2, 6]), (slice(15,16),slice(15,16),slice(0,31),[1,6])]
	
	# need to get lambdas on same grid as HBTdata, so use 2D interpolation
	lambdaF = interpolate.interp2d(pTpts, pphipts, fitLambdadata, kind='linear')
	lambdaFalt = interpolate.interp1d(pTpts, fitLambdadata.dot(pphiwts) / (2.0*pi), kind='cubic')
	KTpts = linspace(0.1, 0.7, 4)
	#KTpts = linspace(0.01, 0.61, 4)
	chosenFitLambda = lambdaFalt(KTpts)
	chosenFitRadii = asarray([ (fitRadiidata[ where( abs( fitRadiidata[:,0] - KT ) < 1.e-3 ) ].flatten())[1:] for KT in KTpts])
		
	iqaxis = 0
	offset = 0.0
	for i in [0,1]:
		for j in [0,1]:
			if (i==1 and j==1):
				offset = 1000.0
				iqaxis -= 1
			CS_data_KT100 = CHUNSHEN_CFdata_KT100[CS_qaxes[CS_qaxisCols[iqaxis]]].reshape([31,2])
			CS_data_KT300 = CHUNSHEN_CFdata_KT300[CS_qaxes[CS_qaxisCols[iqaxis]]].reshape([31,2])
			CS_data_KT500 = CHUNSHEN_CFdata_KT500[CS_qaxes[CS_qaxisCols[iqaxis]]].reshape([31,2])
			CS_data_KT700 = CHUNSHEN_CFdata_KT700[CS_qaxes[CS_qaxisCols[iqaxis]]].reshape([31,2])
			#axs[i,j].semilogy(CFdata[iqaxis,0,:,iqaxis+1], Cpts(CFdata[iqaxis,0,:,iqaxis+1], chosenFitRadii[0,qaxisCols[iqaxis]], chosenFitLambda[0]) - 1.0 + offset, linestyle='-', color='black', linewidth=lw, label=r'$K_T = 100$ MeV')
			#axs[i,j].semilogy(CFdata[iqaxis,1,:,iqaxis+1], Cpts(CFdata[iqaxis,1,:,iqaxis+1], chosenFitRadii[1,qaxisCols[iqaxis]], chosenFitLambda[1]) - 1.0 + offset, linestyle='--', color='red', linewidth=lw, label=r'$K_T = 300$ MeV')
			#axs[i,j].semilogy(CFdata[iqaxis,2,:,iqaxis+1], Cpts(CFdata[iqaxis,2,:,iqaxis+1], chosenFitRadii[2,qaxisCols[iqaxis]], chosenFitLambda[2]) - 1.0 + offset, linestyle='-.', color='green', linewidth=lw, label=r'$K_T = 500$ MeV')
			#axs[i,j].semilogy(CFdata[iqaxis,3,:,iqaxis+1], Cpts(CFdata[iqaxis,3,:,iqaxis+1], chosenFitRadii[3,qaxisCols[iqaxis]], chosenFitLambda[3]) - 1.0 + offset, linestyle=':', color='blue', linewidth=lw, label=r'$K_T = 700$ MeV')
			axs[i,j].semilogy(CS_data_KT100[:,0], abs(CS_data_KT100[:,1]) + offset, '-', color='black', linewidth=lw, label=r'$K_T = 100$ MeV')
			axs[i,j].semilogy(CS_data_KT300[:,0], abs(CS_data_KT300[:,1]) + offset, '--', color='red', linewidth=lw, label=r'$K_T = 300$ MeV')
			axs[i,j].semilogy(CS_data_KT500[:,0], abs(CS_data_KT500[:,1]) + offset, '-.', color='green', linewidth=lw, label=r'$K_T = 500$ MeV')
			axs[i,j].semilogy(CS_data_KT700[:,0], abs(CS_data_KT700[:,1]) + offset, ':', color='blue', linewidth=lw, label=r'$K_T = 700$ MeV')
			axs[i,j].semilogy(CFdata[iqaxis,:,0,iqaxis+1], CFdata[iqaxis,:,0,4] - 1.0 + offset, 'o-', linewidth=0.5, markevery=2, color='black')
			axs[i,j].semilogy(CFdata[iqaxis,:,1,iqaxis+1], CFdata[iqaxis,:,1,4] - 1.0 + offset, 's-', linewidth=0.5, markevery=2, color='red')
			axs[i,j].semilogy(CFdata[iqaxis,:,2,iqaxis+1], CFdata[iqaxis,:,2,4] - 1.0 + offset, '^-', linewidth=0.5, markevery=2, color='green')
			axs[i,j].semilogy(CFdata[iqaxis,:,3,iqaxis+1], CFdata[iqaxis,:,3,4] - 1.0 + offset, 'v-', linewidth=0.5, markevery=2, color='blue')

			if j==1:
				axs[i,j].yaxis.tick_right()
				axs[i,j].yaxis.set_label_position("right")
			if i==0:
				axs[i,j].set_ylim([0.05,1.01])
			else:
				axs[i,j].set_ylim([0.01,1.01])
			
			axs[i,j].set_xlim([-0.0599,0.0599])
			if (i==1 and j==1):
				continue
			#axs[i,j].legend(loc=0, ncol=2, prop={'size': plotfontsize+5})
			axs[i,j].set_xlabel(r'$q_%(qax)s$ (GeV)' % {'qax': qaxes[iqaxis]}, fontsize = labelsize)
			axs[i,j].set_ylabel(r'$C(q_%(qax)s) - 1$' % {'qax': qaxes[iqaxis]}, fontsize = labelsize)
			iqaxis += 1
	
	axs[1,1].legend(loc=10, frameon=False, prop={'size': plotfontsize+5})
	
	plt.show()
	#outfilename = 'correlation_function_OSL_slices_with_averaged_comparison.pdf'
	#plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	#print 'Saved to', outfilename

	

def generate_Rij_plot():
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(2, 2)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0
	lw2 = 1.0
	nKT = 101
	chosenCols = [1,3,5,9]
	nCols = len(chosenCols)
	dims = [nKT, nCols]
	ylims = [[0.0, 8.0],[0.0,8.0],[0.0,9.9],[0.5,1.5]]
	qaxisCols = [2,1,3]	#order of out-side-long columns in data file
	
	HBTdata = loadtxt('HBTradii_GF_cfs_evavg.dat', usecols=tuple(chosenCols)).reshape(dims)
	CHUNSHEN_HBTdata_KT100 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_radii_kt_0_0.2.dat', skiprows=1)
	CHUNSHEN_HBTdata_KT300 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_radii_kt_0.2_0.4.dat', skiprows=1)
	CHUNSHEN_HBTdata_KT500 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_radii_kt_0.4_0.6.dat', skiprows=1)
	CHUNSHEN_HBTdata_KT700 = loadtxt('cross_check_with_Chris/HBT/take_1/HBT_radii_kt_0.6_0.8.dat', skiprows=1)
	#choose slice corresponding to chi^2 minimum, keep only R_o, R_s, R_l (in that order)
	CHUNSHEN_HBTdata_KT100 = CHUNSHEN_HBTdata_KT100[argmin(CHUNSHEN_HBTdata_KT100[:,6])][2:5]
	CHUNSHEN_HBTdata_KT300 = CHUNSHEN_HBTdata_KT300[argmin(CHUNSHEN_HBTdata_KT300[:,6])][2:5]
	CHUNSHEN_HBTdata_KT500 = CHUNSHEN_HBTdata_KT500[argmin(CHUNSHEN_HBTdata_KT500[:,6])][2:5]
	CHUNSHEN_HBTdata_KT700 = CHUNSHEN_HBTdata_KT700[argmin(CHUNSHEN_HBTdata_KT700[:,6])][2:5]
	
	KTpts = linspace(0.1, 0.7, 4)
	CHUNSHEN_HBTdata = vstack((CHUNSHEN_HBTdata_KT100, CHUNSHEN_HBTdata_KT300, CHUNSHEN_HBTdata_KT500, CHUNSHEN_HBTdata_KT700))
	
	iqaxis = 0
	for i in [0,1]:
		for j in [0,1]:
			if (i==1 and j==1):
				axs[i,j].plot(sqrt(0.13957**2+HBTdata[:,0]**2), sqrt(HBTdata[:,2])/sqrt(HBTdata[:,1]), linestyle='-', color='black', linewidth=lw)
				axs[i,j].plot(sqrt(0.13957**2+KTpts**2), sqrt(CHUNSHEN_HBTdata[:,0])/sqrt(CHUNSHEN_HBTdata[:,1]), 's--', color='black', linewidth=lw)
				axs[i,j].set_xlabel(r'$m_T$ (GeV)', fontsize = labelsize)
				axs[i,j].set_ylabel(r'$R_o/R_s$', fontsize = labelsize)
			else:
				axs[i,j].plot(sqrt(0.13957**2+HBTdata[:,0]**2), sqrt(HBTdata[:,qaxisCols[iqaxis]]), linestyle='-', color='black', linewidth=lw)
				axs[i,j].plot(sqrt(0.13957**2+KTpts**2), CHUNSHEN_HBTdata[:,iqaxis], 's--', color='black', linewidth=lw)
				axs[i,j].set_xlabel(r'$m_T$ (GeV)', fontsize = labelsize)
				axs[i,j].set_ylabel(r'$R_{%(qax)s}$ (fm)' % {'qax': qaxes[iqaxis]}, fontsize = labelsize)
			if j==1:
				axs[i,j].yaxis.tick_right()
				axs[i,j].yaxis.set_label_position("right")
			
			axs[i,j].set_xlim([0.0,0.99])
			axs[i,j].set_ylim(ylims[iqaxis])
			iqaxis += 1
	
	#axs[1,1].legend(loc=10, frameon=False, prop={'size': plotfontsize+5})
	
	plt.show()
	#outfilename = 'HBTradii_OSL_slices_with_comparison.pdf'
	#plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	#print 'Saved to', outfilename



def generate_all_plots():
	#generate_CF_plot(0)	#out
	#generate_CF_plot(1)	#side
	#generate_CF_plot(2)	#long
	#generate_CF_plot()
	generate_Rij_plot()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
