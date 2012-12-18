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
  if (vx0 !=0 and vx1 !=0):
    print "AX  vx1  vx0  vxp"
    print Ax, vx1, vx0, vxp
    if (abs(vx1-vx0) > 1.e-2):
      if (vx0 > 0 and vx1 >0):
        return np.log(vx1/vxp) / Ax
      else:
        return np.log(vx0/vxp) / Ax
    else:
      if vx0 > 0:
        return (x1 - xp) / vx0
      if vx0 < 0:
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

  
  step = StepReturn( xnew, dtMin, exitdir, exitface)
  return step
