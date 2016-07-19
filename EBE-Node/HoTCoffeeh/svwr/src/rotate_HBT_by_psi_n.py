# first command-line argument is directory containing files to be rotated (directory)
# second command-line argument is the order of the angle to rotate by (rotation_angle_order)
# 	--> if rotation_angle_order == 0, each harmonic is rotated by the corresponding psi_n
# third command-line argument is either empty string '' or '_no_df'

from numpy import *
import sys
from os import path
from glob import glob

#df_stem = ''
#df_stem = '_no_df'
df_stem = sys.argv[3]

def rotmat(theta):
	return array([[cos(theta),sin(theta)],[-sin(theta),cos(theta)]])

globList = [
	'plane_psis_ev*[0-9]' + df_stem + '.dat',
	'HBTradii_cfs_ev*[0-9]' + df_stem + '.dat',
	'Sourcefunction_variances_cfs_COS' + df_stem + '.dat',
	'Sourcefunction_variances_cfs_SIN' + df_stem + '.dat'
]

# use this to force to look for specific files
#globList = [
#	'avgplane_psis1to1000_no_df.dat',
#	'CavgHBTradii_cfs_evs1to1000_no_df.dat',
#	'PLUMBERGdummyfile1.dat',
#	'PLUMBERGdummyfile2.dat'
#]




def rotate_HBT_cfs(HBTcfsPATH, plane_psis, rotation_angle_order):
	HBT_cfs = loadtxt(HBTcfsPATH)
	plane_psi = plane_psis[rotation_angle_order,1]
	linenum = 0
	for line in HBT_cfs:
		# gives order harmonic of HBT coeffs. to rotate
		forder = float(line[2])
		if (int(rotation_angle_order)==0):
			# i.e., if n == 0
			plane_psi = plane_psis[line[2],1]
		for iCol in range(3,15,2):
			# picks out each R2ij set of cosine and sine coefficients in file
			cols=[iCol,iCol+1]
			# reset line[cols] to rotated line[cols]
			line[cols] = dot(rotmat(forder*plane_psi),line[cols])
		
		# once entire line is rotated, reset HBT_cfs[linenum] to line
		HBT_cfs[linenum] = line
		linenum+=1
	
	if (int(rotation_angle_order)==0):
		rotation_angle_order_string = 'n'
	else:
		rotation_angle_order_string = str(rotation_angle_order)
	HBTcfs_outfile = HBTcfsPATH + '_neq' + rotation_angle_order_string
	savetxt(HBTcfs_outfile, HBT_cfs, fmt = '%i %0.1f %i ' + '%f '*12)





def rotate_SV_cfs(SVcfsCOSPATH, SVcfsSINPATH, plane_psis, rotation_angle_order):
	SV_cfs_COS = loadtxt(SVcfsCOSPATH)
	SV_cfs_SIN = loadtxt(SVcfsSINPATH)
	plane_psi = plane_psis[rotation_angle_order,1]
	for linenum in xrange(len(SV_cfs_COS)):		# same length as SV_cfs_SIN
		# gives order harmonic of SV coeffs. to rotate
		forder = float(SV_cfs_COS[linenum, 2])
		if (int(rotation_angle_order)==0):
			# i.e., if n == 0
			plane_psi = plane_psis[SV_cfs_COS[linenum, 2],1]
		for iCol in range(3,13):
			# picks out each SV set of cosine and sine coefficients in files
			temp=array([SV_cfs_COS[linenum, iCol], SV_cfs_SIN[linenum, iCol]])
			SV_cfs_COS[linenum, iCol], SV_cfs_SIN[linenum, iCol] = dot(rotmat(forder*plane_psi),temp)
		
	
	if (int(rotation_angle_order)==0):
		rotation_angle_order_string = 'n'
	else:
		rotation_angle_order_string = str(rotation_angle_order)
	SV_cfs_COS_outfile = SVcfsCOSPATH + '_neq' + rotation_angle_order_string
	SV_cfs_SIN_outfile = SVcfsSINPATH + '_neq' + rotation_angle_order_string
	savetxt(SV_cfs_COS_outfile, SV_cfs_COS, fmt = '%i %0.1f %i ' + '%f '*10)
	savetxt(SV_cfs_SIN_outfile, SV_cfs_SIN, fmt = '%i %0.1f %i ' + '%f '*10)






if __name__ == "__main__":
	directory = sys.argv[1]
	rotation_angle_order = sys.argv[2]
	
	planepsisPATH = glob(path.join(directory, globList[0]))[0]
	HBTcfsPATH = glob(path.join(directory, globList[1]))[0]
	SVcfsCOSPATH = ''
	SVcfsSINPATH = ''
	if (path.exists(path.join(directory, globList[2])) and path.exists(path.join(directory, globList[3]))):
		SVcfsCOSPATH = glob(path.join(directory, globList[2]))[0]
		SVcfsSINPATH = glob(path.join(directory, globList[3]))[0]
	
	# check that files exist before rotating them
	if path.exists(planepsisPATH):
		#print 'plane psis file exists'
		plane_psis = loadtxt(planepsisPATH)
		#plane_psi = plane_psis[rotation_angle_order,1]
	else:
		print 'Couldn\'t find file containing plane psis!  Aborting...'
		sys.exit(0)
	
	if path.exists(HBTcfsPATH):
		#print 'HBT cfs file exists'
		rotate_HBT_cfs(HBTcfsPATH, plane_psis, rotation_angle_order)
	
	if (path.exists(SVcfsCOSPATH) and path.exists(SVcfsSINPATH)):
		#print 'SV cfs COS and SV cfs SIN files exist'
		rotate_SV_cfs(SVcfsCOSPATH, SVcfsSINPATH, plane_psis, rotation_angle_order)

# End of file
