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
nqpts = 21

qAxisOpts = ['X', 'Y', 'Z']
qAxesLC = ['x', 'y', 'z']
projectionOpts = ['', '_unprojected']
resFracsOpts = ['0.00', '0.10', '0.20', '0.60', '1.00']
PTOpts = ['15']
#PTOpts = ['31']
PYOpts = ['15']
#QTOpts = ['13', '15', '17', '19', '21', '23', '25', '27', '29', '31']
QTOpts = ['7', '9', '11', '13', '15', '17', '19', '21']

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

GeVToMeV = 1000.0
#GeVToMeV = 1.0	#uses GeV instead of MeV
panelLabels = ['(a)', '(b)', '(c)']
panelCounter = 0

nCols = 2
comparisonPath = './' \
					+ 'AXIS_X_pT%(opt1)s_pY%(opt2)s_qt%(opt3)s/' \
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
	
	#f = interpolate.interp1d(datax, datay, kind='cubic')
	f = interpolate.interp1d(datax, datay, kind='linear')
	plotdatax = linspace(xlower, xupper, 1001)
	plotdatay = f(plotdatax)

	#plotdatax = datax
	#plotdatay = datay
	
	#return plotdatax, plotdatay-1.0, [xlower, xupper, ylower-1.0, yupper-1.0]
	return plotdatax, plotdatay, [xlower, xupper, ylower, yupper]
	#return plotdatax, plotdatay, [xlower, xupper, 1.0, 2.0]
	#return plotdatax, plotdatay, [GeVToMeV * xlower, GeVToMeV * xupper, 1.0, 2.0]

#############################################################################
#############################################################################

def generate_plot(ipT, ipphi, inpt, inqt, ires, iprojection):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	#ax.set_yscale("log")
	fig.subplots_adjust(wspace=0.0, hspace=0.0)
	
	plotPath = comparisonPath % {'opt1': PTOpts[inpt], 'opt2': QTOpts[inqt], 'opt3': resFracsOpts[ires], 'opt4': projectionOpts[iprojection]}
	
	cols = [2+0, 9]	# chooses q-axes, CF value from CF files
	plotdatax, plotdatay, lims = generate_plotdata(plotPath, ipT, ipphi, int(PTOpts[inpt]), npphi, cols)
	cols = [2+0, 6]	# chooses q-axes, CF value from CF files
	tplotdatax, tplotdatay, tlims = generate_plotdata(plotPath, ipT, ipphi, int(PTOpts[inpt]), npphi, cols)
	ax.plot(plotdatax, plotdatay, linestyle='-', color='r', linewidth=1.5, label=r'$p_T = %(opt1)s, q_t = %(opt2)s$, RF = %(opt3)s' % {'opt1': PTOpts[inpt], 'opt2': QTOpts[inqt], 'opt3': resFracsOpts[ires], 'opt4': projectionOpts[iprojection]})
	ax.plot(tplotdatax, 1.0+tplotdatay, linestyle='--', color='b', linewidth=1.5, label='Thermal part')
	ax.axis(lims)
	ax.set_xlabel(r'$q_x$ (GeV)' % {'opt1': qAxesLC[0]}, {'fontsize': plotfontsize + 5})	
	ax.set_ylabel(r'$C$', {'fontsize': plotfontsize + 5})
	ax.legend(loc=0, prop={'size': plotfontsize})
	
	plt.show(block=False)


#############################################################################
#############################################################################

def generate_comparison_plot(ipT, ipphi, ires, iprojection):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	#ax.set_yscale("log")
	fig.subplots_adjust(wspace=0.0, hspace=0.0)

	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=float(len(PTOpts)*len(PYOpts)*len(QTOpts)))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	
	idx = 0
	for inpt in xrange(len(PTOpts)):
		for inpy in xrange(len(PYOpts)):
			for inqt in xrange(len(QTOpts)):
				plotPath = comparisonPath % {'opt1': PTOpts[inpt], 'opt2': PYOpts[inpy], 'opt3': QTOpts[inqt], 'opt4': resFracsOpts[ires], 'opt5': projectionOpts[iprojection]}
				colorVal = scalarMap.to_rgba(float(idx))
				cols = [2+0, 9]	# chooses q-axes, CF value from CF files
				plotdatax, plotdatay, lims = generate_plotdata(plotPath, ipT, ipphi, int(PTOpts[inpt]), npphi, cols)
				ax.plot(plotdatax, plotdatay, linestyle='-', color=colorVal, linewidth=1.5, label=r'$p_T = %(opt1)s, p_Y = %(opt2)s, q_t = %(opt3)s$' % {'opt1': PTOpts[inpt], 'opt2': PYOpts[inpy], 'opt3': QTOpts[inqt]})
				#ax.plot(tplotdatax, 1.0+tplotdatay, linestyle='--', color='b', linewidth=1.5, label='Thermal part')
				ax.axis(lims)
				idx += 1

	ax.set_xlabel(r'$q_x$ (GeV)' % {'opt1': qAxesLC[0]}, {'fontsize': plotfontsize + 5})	
	ax.set_ylabel(r'$C$', {'fontsize': plotfontsize})
	ax.legend(loc=0, prop={'size': plotfontsize+5})
	plt.title(r'RF$=%(opt4)s$, %(opt5)s' % {'opt4': resFracsOpts[ires], 'opt5': projectionOpts[iprojection]})
	
	plt.show(block=False)


#############################################################################
#############################################################################


def generate_all_plots():
	#chosenpT, chosenpphi = 0, 0
	#for i in xrange(7):
	#	generate_plot(chosenpT, chosenpphi, 3, i, 1, 1)
	#
	#generate_comparison_plot(chosenpT, chosenpphi, 0, 1, 1)
	#chosenpT = 5
	#generate_comparison_plot(chosenpT, chosenpphi, 1, 1, 1)
	#chosenpT = 6
	#generate_comparison_plot(chosenpT, chosenpphi, 3, 1, 1)
	#generate_plot(4, 0, 0, 0, 1, 1)
	#generate_plot(6, 0, 1, 0, 1, 1)
	chosenpphi = 0
	#for chosenpT in xrange(0,15,3):
	for chosenpT in xrange(9):
		#print chosenpT
		generate_comparison_plot(chosenpT, chosenpphi, 3, 0)
	pause()

#############################################################################
#############################################################################

if __name__ == "__main__":
	generate_all_plots()

	print 'Finished all.'

# End of file
