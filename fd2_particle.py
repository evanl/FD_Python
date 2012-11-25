#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from time import clock, time
import sys

import fd2_functions as fd2

xw = 0.0; xe = 70; ys = 0.0; yn = 70;

if len(sys.argv) == 1 :
  print "please specify nx"
  sys.exit()

ns = int(sys.argv[1])

nx = ns
ny = ns
print nx , ny 

compressibility = 5.0 * pow(10,-10)
thickness = 10.0;
S = 1000.0 * 9.81 * compressibility * thickness

tmax = 0.5
nt = 5
                          #west,        east        south        north
bc = fd2.BoundaryCondition("head",5.0,"head", 10.0,"flux", 0.0,"flux", 0.0)

dx = (xe - xw)/(nx-1)
dy = (yn - ys)/(ny-1)

tInit = time()

Y,X = np.mgrid[ ys:yn:ny*1j, xw:xe:nx*1j]

H = fd2.SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, t = 0, S = 0, dt = 1.0, r= np.zeros((1,1)), h0 = 0.0 , K_B = 0. )

tFinal = time()


print "Nx = " + str(nx)
print "Solve time=  " + str( tFinal - tInit)

plottype = 0
if plottype !=0:
  if plottype == "line":
    fig = plt.figure()
    ax = fig.add_subplot(111)
    p1 = ax.plot(X[1,:] , H[1,:])
    ax.set_xlabel('x-direction[m]')
    ax.set_ylabel('average hydraulic head [m]')
    fig.suptitle("Center-line in Domain")
  else:
    fig = plt.figure()
    ax = Axes3D(fig)
    surf = ax.plot_surface(X,Y,H, rstride=1, cstride =1, linewidth =0, cmap = cm.jet, antialiased=False)
    CB = plt.colorbar(surf, shrink = 0.8, extend = 'both')
  plt.show()

# particle tracking stuff goes here: 

