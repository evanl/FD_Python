#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from time import clock, time
import scipy

import fd2_functions as fd2

xw = 0.0
xe = 2250./3.281

ys = 0.0
yn = 2500./3.281

n = 100
nx = n
ny = n

bc = fd2.BoundaryCondition("flux", 0.0, "flux", 0.0, "flux", 0.0, "flux", 0.0)

dx = (xe - xw)/(nx-1)
dy = (yn - ys)/(ny-1)

Y,X = np.mgrid[ ys:yn:ny*1j, xw:xe:nx*1j]

#S = 1000.0 * 9.81 * compressibility * thickness
S = 0.2

recharge = -20. # in/yr 
N = recharge * 2.54 / (365 * 100)

wellvals = []
wellvals.append(270 * 5.451)
wellvals.append(700 * 5.451)
wellvals.append(400 * 5.451)

wellLoc = []
wellLoc.append([729 /3.281, 250/3.281])
wellLoc.append([950/3.281, 1017/3.281])
wellLoc.append([1150/3.281, 1389/3.281])

print wellvals
print wellLoc

H = fd2.SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, t = 0, S = 0, dt = 1.0,
    r= np.zeros((1,1)), h0 = 44 / 3.281, K_B = 10/3.281, infil = True, infilval = N, 
    Qwell1 = wellvals, wellLoc1 = wellLoc, wells = True)

for i in range(nx):
  for j in range(ny):
    X[i,j] = X[i,j] * 3.281
    Y[i,j] = Y[i,j] * 3.281
    H[i,j] = H[i,j] * 3.281

fig = plt.figure()
ax = fig.add_subplot(111)
surf = ax.contourf(X,Y,H,55,cmap = cm.Greys,
    linewidth = (3,))
CB = plt.colorbar(surf ) 
ax.set_title('Hydraulic Head Contours [ft]')
ax.set_xlabel('x-direction [ft]')
ax.set_ylabel('y-direction [ft]')
ax.xaxis.grid(True)
ax.yaxis.grid(True)
plt.show()

