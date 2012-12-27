#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from time import clock, time
import sys

import fd2_functions as fd2
import fd2_part_functions as fd2part

k = 50
phi = 0.25

xw = 0.0; xe = 5000; ys = 0.0; yn = 5000;

if len(sys.argv) == 1 :
  print "please specify nx"
  sys.exit()

ns = int(sys.argv[1])

nx = ns
ny = ns
                          #west,        east        south        north
bc = fd2.BoundaryCondition("head",50.0,"head", 25.0,"flux", 0.0,"flux", 0.0)

dx = (xe - xw)/(nx-1)
dy = (yn - ys)/(ny-1)

tInit = time()

Y,X = np.mgrid[ ys:yn:ny*1j, xw:xe:nx*1j]

H = fd2.SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, \
    t = 0, S = 0, dt = 1.0, r= np.zeros((1,1)), h0 = 0.0 , K_B = 0.,\
    well1=True, Qwell1=0.3, wellLoc1=[2500,2500])

tFinal = time()

print "Nx = " + str(nx)
print "Solve time=  " + str( tFinal - tInit)

plottype = "surf"
matplotlib.rcParams.update({'font.size': 14})

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
    ax = fig.add_subplot(111)
    surf = ax.contour(X,Y,H,25,cmap = cm.bone)
#    ax = Axes3D(fig)
#    surf = ax.plot_surface(X,Y,H,25,cmap = cm.jet )
    CB = plt.colorbar(surf, orientation='horizontal', ticks=[25,45])
    CB.ax.set_xticklabels(['25' ,'45'])
    ax.set_xlabel('x-direction')
    ax.set_ylabel('y-direction')
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

# particle tracking 
xp = 100.
yp = 3000.
xpart, ypart, tpart = fd2part.ParticleTrack(xp, yp, H, X, Y, dx, dy, nx, ny, k, phi)


# plotting tools
X1 = np.asarray(xpart)
Y1 = np.asarray(ypart)
fig.add_subplot(111)
p2 = ax.plot(X1,Y1,label = "particle tracks",color = 'r', marker = 'x', linestyle = '--')
l = plt.legend(bbox_to_anchor=(1.05,1), loc=1, borderaxespad=0)
plt.show()
fig.savefig('track.png')
