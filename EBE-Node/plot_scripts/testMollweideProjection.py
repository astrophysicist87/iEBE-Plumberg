import numpy as np
import pylab as plt

ra = np.linspace(-np.pi, np.pi, 40)
dec= np.linspace(-np.pi/2, np.pi/2, 20)
X,Y = np.meshgrid(ra,dec)
Z = np.sin(X) * np.cos(X) * np.sin(Y) * np.cos(Y)

plt.figure()
ax = plt.subplot(111, projection = 'mollweide')
ax.contourf(X,Y,Z,100,cmap=plt.cm.gist_rainbow)
#ax.contour(X,Y,Z,10,colors='k')

plt.show()