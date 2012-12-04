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
# particle tracking stuff goes here: 

tpart = 4900
t = 0.

xpart = []
ypart = []
xpart.append( 100.)
ypart.append( 3000.)
outOfDomain = False
# initialize location and velocity in step terms. 
xp = []
vp = []
dxi = []
xp.append(xpart[0])
xp.append(ypart[0])
vp.append(1.)
vp.append(1.)
dxi.append( dx)
dxi.append( dy)

stepcount = 0
i=0
j=0
while X[1,i] < xpart[stepcount]:
  while Y[j,1] < ypart[stepcount]:
    j +=1
  i+=1

v0[0] = vel[0]
v0[1] = vel[2]
v1[0] = vel[1]
v1[1] = vel[3]
x0[0] = X[1,i] - dxi[0]/2
x0[1] = Y[j,1] - dxi[1]/2
x1[0] = X[1,i] + dxi[0]/2
x1[1] = Y[j,1] + dxi[1]/2
  
while t < tpart and outOfDomain == False and 1 < i < nx-1 and 1 < j < ny-1:
  
  while X[1,i] + dx/2 < xpart[stepcount]:
    while Y[j,1] +dy/2 < ypart[stepcount]:
      j +=1
    i+=1
  vel = fd2part.CalcFaceVelocities( H[j,i], H[j,i-1], H[j,i+1], H[j-1,i], H[j+1,i], dx, dy, k, phi)
  #creates vectors to pass into step function
  v0,v1,x0,x1=[0]*2, [0]*2, [0]*2, [0]*2  

  step = fd2part.ParticleStep(xp, vp, v0, v1, x0, x1, dxi)

  print "dtmin \/"
  print step.dtMin
  print "\n\n"
  if step.dtMin != False:
    # increment necessary values to get to next timestep. 
    t += step.dtMin    
    xpart.append(step.xnew[0])
    ypart.append(step.xnew[1])
    xp[0] = step.xnew[0]
    xp[1] = step.xnew[1]
    if step.exitdir == 0:
      if step.exitface == 1 :
        i +=1
        vp[0] = v1[0]
        vp[1] = v1[1]
      else:
        i -=1
        vp[0] = v0[0]
        vp[1] = v0[0]
    elif step.exitdir ==1:
      if step.exitface ==1: 
        j+=1
        vp[1] = v1[1]
        vp[0] = v1[0]
      else:
        j-=1
        vp[1] = v0[1]
        vp[0] = v0[0]
  else:
    outOfDomain = True

  stepcount +=1


X1 = np.asarray(xpart)
Y1 = np.asarray(ypart)
fig.add_subplot(111)
p2 = ax.plot(X1,Y1,label = "particle tracks",color = 'r', marker = 'x', linestyle = '--')
l = plt.legend(bbox_to_anchor=(1.05,1), loc=1, borderaxespad=0)
plt.show()
fig.savefig('track.png')
