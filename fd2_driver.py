#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from time import clock, time
import scipy

import fd2_functions as fd2

xw = 0.0
xe = 70

ys = 0.0
yn = 70

nx = 15
ny = 15

bc = fd2.BoundaryCondition("head", 10.0, "head", 10.0, "head", 10.0, "head", 10.0)

dx = (xe - xw)/(nx-1)
dy = (yn - ys)/(ny-1)

Y,X = np.mgrid[ ys:yn:ny*1j, xw:xe:nx*1j]

compressibility = 5.0 * pow(10,-10)
thickness = 10.0;
S = 1000.0 * 9.81 * compressibility * thickness

tmax = 0.5
nt = 5

#H1 = fd2.TimeSolve( xw , xe, ys, yn, nx, ny, dx, dy, bc, tmax, nt, S)

HL = fd2.SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, t = 0, S = 0, dt = 1.0, r= np.zeros((1,1)), h0 = 11.0, K_B = 0.00005/5.0)

fig1 = plt.figure()
fig1.suptitle("Center-line in Domain")
ax = fig1.add_subplot(111)
for i in range(nt):
	ax.plot( X[7,:], H1[7,:,i], label=("t = "+str(i * (tmax/(nt-1))) + " days") )
ax.plot(X[7,:] , HL[7,:], label = "leaky case")
ax.set_xlabel('x -direction [m]')
ax.set_ylabel('average head [m]')
ax.legend()


fig = plt.figure()
ax = fig.add_subplot(111)
p1 = ax.plot(X[7,:] , HL[7,:])
ax.set_xlabel('x-direction[m]')
ax.set_ylabel('average hydraulic head [m]')
fig.suptitle("Center-line in Domain")


plt.show()
