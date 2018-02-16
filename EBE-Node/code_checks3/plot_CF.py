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

npy = 15
npphi = 36
nqpts = 51

qAxisOpts = ['X', 'Y', 'Z']
qAxesLC = ['x', 'y', 'z']
projectionOpts = ['', '_unprojected']
FFOpts = ['1.0']
resFracsOpts = ['0.10']
PTOpts = ['15']
PYOpts = ['15']
QTOpts = ['17']

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

GeVToMeV = 1000.0
#GeVToMeV = 1.0	#uses GeV instead of MeV
panelLabels = ['(a)', '(b)', '(c)']
panelCounter = 0

#nCols = 2
comparisonPath = './' \
					+ 'AXIS_X_CHECK3/' \
					+ 'RESFRAC_%(opt4)s/' \
					+ 'results/correlfunct3D%(opt5)s_Pion_+.dat'

#############################################################################
#############################################################################
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

#############################################################################
#############################################################################

def generate_plotdata(path, ipT, ipphi, npT, npphi, cols):
	# shape data arrays appropriately
	if cols[0]==2:
		nqx, nqy, nqz = nqpts, 1, 1
	elif cols[0]==3:
		nqx, nqy, nqz = 1, nqpts, 1
	elif cols[0]==4:
		nqx, nqy, nqz = 1, 1, nqpts
	
	nCols = len(cols)
	
	data = loadtxt(path, usecols=tuple(cols)).reshape([npT, npphi, nqx, nqy, nqz, nCols])

	if nCols==2:
		datax = data[ipT, ipphi, :, 0, 0, 0]
		datay = data[ipT, ipphi, :, 0, 0, 1]

		if cols[0]==3:
			datax = data[ipT, ipphi, 0, :, 0, 0]
			datay = data[ipT, ipphi, 0, :, 0, 1]
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

	#datay += 1.0
	
	xlower, xupper = min(datax), max(datax)
	ylower, yupper = min(datay), max(datay)
	
	#f = interpolate.interp1d(datax, datay, kind='cubic')
	f = interpolate.interp1d(datax, datay, kind='linear')
	plotdatax = linspace(xlower, xupper, 1001)
	plotdatay = f(plotdatax)

	#plotdatax = datax
	#plotdatay = datay
	
	#return plotdatax, plotdatay-1.0, [xlower, xupper, ylower-1.0, yupper-1.0]
	return plotdatax, plotdatay, [xlower, xupper, ylower, yupper]
	#return plotdatax, plotdatay, [xlower, xupper, 1.0, 1.05]
	#return plotdatax, plotdatay, [GeVToMeV * xlower, GeVToMeV * xupper, 1.0, 2.0]

#############################################################################
#############################################################################

#############################################################################
#############################################################################

def generate_comparison_plot(ipT, ipphi, iprojection):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	#ax.set_yscale("log")
	fig.subplots_adjust(wspace=0.0, hspace=0.0)

	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=float(len(PTOpts)*len(PYOpts)*len(QTOpts)*len(resFracsOpts)*len(FFOpts)))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	
	idx = 0
	for inpt in xrange(len(PTOpts)):
		for inpy in xrange(len(PYOpts)):
			for inqt in xrange(len(QTOpts)):
				for iRF in xrange(len(resFracsOpts)):
					for iFF in xrange(len(FFOpts)):
						plotPath = comparisonPath % {'opt1': PTOpts[inpt], 'opt2': PYOpts[inpy], 'opt3': QTOpts[inqt], 'opt4': resFracsOpts[iRF], 'opt5': projectionOpts[iprojection], 'opt6': FFOpts[iFF]}
						colorVal = scalarMap.to_rgba(float(idx))
						cols = [2+0, 6, 7, 8, 9]	# chooses q-axes, CF value from CF files
						#cols = [2+0, 9]	# chooses q-axes, CF value from CF files
						plotdatax, plotdatay, lims = generate_plotdata(plotPath, ipT, ipphi, int(PTOpts[inpt]), npphi, cols)
						ax.plot(plotdatax, plotdatay, linestyle='-', color=colorVal, linewidth=1.5, label=r'$p_T = %(opt1)s$' % {'opt1': PTOpts[inpt]})
	#                               label=r'$p_T = %(opt1)s, p_Y = %(opt2)s, q_t = %(opt3)s, RF = %(opt4)s$' \
	#                               % {'opt1': PTOpts[inpt], 'opt2': PYOpts[inpy], 'opt3': QTOpts[inqt], 'opt4': resFracsOpts[iRF]})
						#ax.plot(tplotdatax, 1.0+tplotdatay, linestyle='--', color='b', linewidth=1.5, label='Thermal part')
						ax.axis(lims)
						idx += 1

	ax.set_xlabel(r'$q_x$ (GeV)' % {'opt1': qAxesLC[0]}, {'fontsize': plotfontsize + 5})	
	ax.set_ylabel(r'$C$', {'fontsize': plotfontsize})
	#ax.legend(loc=0, prop={'size': plotfontsize+5})
	#plt.title(r'%(opt5)s' % {'opt4': resFracsOpts[ires], 'opt5': projectionOpts[iprojection]})
	#filename = 'cmp_thermal_kt%(ikt)i.out' % {'ikt': ipT}
	#filename = 'avgSparse_rho_decay_pions_pt%(ipt)i.dat' % {'ipt': ipT}
	#savetxt(filename, vstack((plotdatax, plotdatay)).T)
	
	#plt.show(block=False)
	filename = 'avgSparse_rho_decay_pions_NO_OSC_pt%(ipt)i.pdf' % {'ipt': ipT}
	#filename = 'thermal_rho_CF_pt%(ipt)i.pdf' % {'ipt': ipT}
	plt.savefig(filename, format='pdf', bbox_inches='tight')
	print 'Saved to', filename


#############################################################################
#############################################################################


def generate_all_plots():
	chosenpphi = 0
	#for chosenpT in xrange(0,15,3):
	#for chosenpT in [0, 9, 18]:
	#for chosenpT in [0, 6, 12]:
	for chosenpT in [0, 4, 8]:
		#print chosenpT
		generate_comparison_plot(chosenpT, chosenpphi, 1)
	#pause()

#############################################################################
#############################################################################

if __name__ == "__main__":
	generate_all_plots()

	print 'Finished all.'

# End of file
