import numpy as np


def CalcFaceVelocities( hC, hW, hE, hS, hN, dx, dy, k, phi):
  vel = []
  vel.append(k * ( hW - hC) / (dx * phi))
  vel.append( k * ( hC - hE) / (dx * phi))
  vel.append( k * ( hS - hC) / (dy * phi))
  vel.append( k * ( hC - hN) / (dy * phi))
  #returns w e s n
  return vel 
    
def CheckExit(vx0, vx1, vxp):
  # returns False if no exit is possible in this direction. 
  # returns 0 if an exit across face 0 is possible. 
  # returns 1 if an exit across face 1 is possible. 
  
  # positive velocities
  print " \n CheckExit"
  if (vx0 > 1.e-6 and vx1 > 1.e-6):
    return 1

  #negative velocities
  elif(vx0 < -1.e-6 and vx1 < -1.e-6):
    return 0

  # sink condition
  elif( vx0 > 1.e-6 and vx1 < -1.e-6):
    return "flag"

  # flow divide condition
  elif(vx0 < -1.e-6 and vx1 > 1.e-6):
    if (vxp > 0 ):
      return 1
    else:
      return 0

  # no flow on left  edge. 
  elif ( vx0 * vx0 < 1.e-12):
    if vx1 > 0:
      return 1
    else:
      return "flag"

  # no flow on right edge
  elif (vx1 * vx1 < 1.e-12):
    if vx0 < 0:
      return 0
    else:
      return "flag"
  else:
    print "Invalid velocity condition or input error"
    exit(1)

def CalcTravelTime(Ax, x0, x1, xp, vx0, vx1, vxp):
  # assumes exit face conditions have been checked
  print "\n CalcTravelTime"
  if (vx0 !=0 and vx1 !=0):
    print "AX  vx1  vx0  vxp"
    print Ax, vx1, vx0, vxp
    if (abs(vx1 - vx0) > 1.):
      if (vx0 > 0 and vx1 >0):
        print " --------------------------1"
        return np.log(vx1/vxp) / Ax
      else:
        print " --------------------------2"
        return np.log(vx0/vxp) / Ax
    else:
      if vx0 > 0:
        print "===========================3"
        print x1, xp, vx0
        return (x1 - xp) / vx0
      if vx0 < 0:
        print "============================4"
        print x0, xp, vx0
        return (x0 - xp) / vx0
  else:
    return float('Inf')                 

def CalcPosition (x0, Ax, vx0, vxp, dt):
  if Ax > 1e-2:
    return  x0 + (vxp * np.exp( Ax * dt) - vx0) / Ax
  else:
    return x0 + dt * vx0

class StepReturn(object):

  def __init__ (self, xnew, dtMin, exitdir, exitface):
    self.xnew = xnew
    self.dtMin = dtMin
    self.exitdir = exitdir
    self.exitface = exitface

def ParticleStep(xp, vp, v0, v1, x0, x1, dx):
  # note, each of the inputs is a two element list.
  # eg. xp[0] = xp, xp[1] = yp
  A = []
  checkexit = []
  dt = [10000]*2

  print "\n particleStep"
  print "location"
  print xp[0], xp[1]

  for i in range(2):
   
    checkexit.append( CheckExit(v0[i],v1[i],vp[i]))

    print "velocities, exit status"
    print v0[i], v1[i], checkexit[i]
    
    A.append((v1[i] - v0[i])/dx[i])
    if checkexit[i] != "flag":
      dt[i] = CalcTravelTime(A[i], x0[i], x1[i], xp[i], v0[i], v1[i], vp[i])
      if dt[i] < 0:
        dt[i] = dt[i] * -1
  # All velocities are into the current cell. 
  if (checkexit[0] == "flag" and checkexit[1] == "flag"):
    return StepReturn(False,False,False,False)
  print "dt \/"
  print dt
  dtMin = min(dt)

  # some clever stuff here. exitdir is binary, 0 is x direction, 1 is y direction
  # exitface is 0 for lower, 1 for upper. 
  exitdir = dt.index(min(dt))
  exitface = checkexit[exitdir]
  print "exitdir"
  print exitdir
  print "exitface"
  print exitface
  xnew = []
  for i in range(2):
    if (v0[i] !=0.0 and v1[i] !=0 ):
      xnew.append(CalcPosition(xp[i], A[i], v0[i], vp[i], dtMin))
    else:
      xnew.append(xp[i])

  if abs(xnew[0]-xp[0]) > dx[0]:
    print "x stepped out of bounds"
    #exit(1)
  if abs(xnew[1] - xp[1]) > dx[1]:
    print "y stepped out of bounds"
    #exit(1)
       
  step = StepReturn( xnew, dtMin, exitdir, exitface)
  return step


def ParticleTrack(xpin, ypin, H, X, Y, dx, dy, nx, ny, k, phi):
  """ takes in an initial position along with a grid of head values 
  and grid cells. Starts at t=0.
  Returns a list of x particle locations.
  """
  tmax = 4900
  t = 0.

  xpart = []
  ypart = []
  tpart = []
  xpart.append(xpin)
  ypart.append(ypin) 
  tpart.append(t)
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
  while X[1,i] < xpin:
    i+=1
  while Y[j,1] < ypin:
    j +=1

    
  while (t < tmax and outOfDomain == False and 
        1 < i < nx-1 and 1 < j < ny-1):
    
    vel = CalcFaceVelocities( H[j,i], H[j,i-1], H[j,i+1], H[j-1,i],
                                      H[j+1,i], dx, dy, k, phi)
    #creates vectors to pass into step function

    v0,v1,x0,x1=[0]*2, [0]*2, [0]*2, [0]*2  

    v0[0] = vel[0]
    v0[1] = vel[2]
    v1[0] = vel[1]
    v1[1] = vel[3]
    x0[0] = X[1,i] - dxi[0]/2
    x0[1] = Y[j,1] - dxi[1]/2
    x1[0] = X[1,i] + dxi[0]/2
    x1[1] = Y[j,1] + dxi[1]/2
    print i, j 
    print "edges x0, x1, y0, y1"
    print x0[0], x0[1], x1[0], x1[1]

    print "xp", xp
    step = ParticleStep(xp, vp, v0, v1, x0, x1, dxi)

    print "dtmin \/"
    print step.dtMin
    if step.dtMin != False:
      # increment necessary values to get to next timestep. 
      t += step.dtMin    

      print "xpinit, ypinit"
      print xpart[-1], ypart[-1]

      xpart.append(step.xnew[0])
      ypart.append(step.xnew[1])
      print "xpend, ypend"
      print xpart[-1], ypart[-1]

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

    tpart.append(t)
    vp[0] = (xpart[-1] - xpart[-2])/(tpart[-1] - tpart[-2])
    vp[1] = (ypart[-1] - ypart[-1])/(tpart[-1] - tpart[-2])
    stepcount +=1
    print "\n\n"

   
  return xpart, ypart, tpart
