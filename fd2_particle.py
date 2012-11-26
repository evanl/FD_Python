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
bc = fd2.BoundaryCondition("head",25.0,"head", 0.0,"flux", 0.0,"flux", 0.0)

dx = (xe - xw)/(nx-1)
dy = (yn - ys)/(ny-1)

tInit = time()

Y,X = np.mgrid[ ys:yn:ny*1j, xw:xe:nx*1j]

H = fd2.SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, t = 0, S = 0, dt = 1.0, r= np.zeros((1,1)), h0 = 0.0 , K_B = 0. )

tFinal = time()


print "Nx = " + str(nx)
print "Solve time=  " + str( tFinal - tInit)

plottype = "surf"
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
  #plt.show()
plt.savefig('head.png')
# particle tracking stuff goes here: 

tpart = 1
t = 0

xpart = []
ypart = []
xpart.append( 5.)
ypart.append( 5.)


i=0
j=0
while X[1,i] < xpart[0]:
  while Y[j,1] < ypart[0]:
    j +=1
  i+=1

print i, j 
print X[1,i], Y[j,1]

t = 0
outOfDomain = False
# initialize location and velocity in step terms. 
xp = []
vp = []
dxi = []
xp.append(xpart[0])
xp.append(ypart[0])
vp.append(0.000001)
vp.append(0.000001)
dxi.append( dx)
dxi.append( dx)

while t < tpart and outOfDomain == False :
  
  vel = fd2.CalcFaceVelocities( H[j,i], H[j,i-1], H[j,i+1], H[j-1,i], H[j+1,i], 0.001, 0.25)
  #creates vectors to pass into step function
  v0,v1,x0,x1=[0]*2, [0]*2, [0]*2, [0]*2  
  
  v0[0] = vel[0]
  v0[1] = vel[2]
  v1[0] = vel[1]
  v1[1] = vel[3]
  
  x0[0] = X[1,i] - dxi[0]
  x0[1] = Y[j,1] - dxi[1]
  x1[0] = X[1,i] + dxi[0]
  x1[1] = Y[j,1] + dxi[1]
  

  step = fd2.ParticleStep(xp, vp, v0, v1, x0, x1, dxi)

  if step.dTmin != False:
    # increment necessary values to get to next timestep. 
    t += step.dtMin    
    xpart.append(step.xnew[0])
    ypart.append(step.xnew[1])
    xp[0] = step.xnew[0]
    xp[1] = step.xnew[1]
    if step.exitdir == 0:
      if step.exitface == 1 :
        i +=1
      else:
        i -=1
    elif step.exitdir ==1:
      if step.exitface ==1: 
        j+=1
      else:
        j-=1
  else:
    outOfDomain = True

  
