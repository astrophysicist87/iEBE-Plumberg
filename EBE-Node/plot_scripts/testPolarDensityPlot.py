import numpy as np
import pylab as plt
import scipy.interpolate as interp
import sys

'''
rPts = np.loadtxt(rptsFile)
thetaPts = np.loadtxt(thetaptsFile)
NrPts, rmin, rmax = len(rPts), min(rPts), max(rPts)
NthetaPts, thetamin, thetamax = len(thetaPts), min(thetaPts), max(thetaPts)
'''
'''
def read_in_data(filename):
    data = np.loadtxt(filename)
    rpts = data[:,0]
    thetapts = data[:,1]
    Rgrid,THETAgrid = np.meshgrid(rpts, thetapts)
    return Rgrid, THETAgrid, data
'''

filename = sys.argv[1]
#data = np.loadtxt('C:\Users\Christopher Plumberg\Desktop\\results-1\\results\lambdas.dat')
data = np.loadtxt(filename, usecols=(0,1,3))

r = np.linspace(min(data[:,0]), 1.00, 100)
theta = np.linspace(0.00, 2.0*np.pi, 1000)
R,THETA = np.meshgrid(r,theta)
#Z = np.exp(-R**2*2.0) * np.cos(3.0*THETA)

npt, npphi = 15, 36
shift = np.asarray([0.0, 2.0*np.pi, 0.0])
dataT = data.reshape([npt, npphi, 3])
dataT = np.concatenate((dataT[:,-1].reshape([npt,1,3])-shift, dataT, dataT[:,0].reshape([npt,1,3]) + shift), axis=1).reshape([npt*(npphi+2), 3])
points = np.vstack((dataT[:,0],dataT[:,1])).T
#points = points.reshape([npt, npphi+2, 2])
values = dataT[:,2]

Z = interp.griddata(points, values, (R, THETA), method='cubic', fill_value=0.0)

'''
print points
print values
print R
print THETA
print Z
'''

plt.figure()
ax = plt.subplot(111, polar=True)
cm = plt.cm.gist_rainbow
#ax.pcolormesh(THETA,R,Z, cmap=plt.cm.gist_rainbow)
plt.pcolormesh(THETA,R,Z, cmap=cm)
#ax.pcolormesh(Y,X,Z, cmap=plt.cm.plasma)
#ax.axis('off')
plt.colorbar()

plt.show()
