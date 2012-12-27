import numpy as np


class ParticleField(object):

  def __init__(xp_in, yp_in, headMatrix, xLocations, yLocations,\
               dx, dy, nx, ny, conductivity, porosity):
    """
    xp_in and yp_in are initial coordinates of particle
    headMatrix is a numpy array of head values
    where xLocations and yLocations correspond to the grid centers
    dx and dy are grid spacing.
    nx and ny are stored for convenience
    conductivity and porosity are constant across the domain in this version.
    """
    self.xp = [xp_in]
    sely.yp = [yp_in]
    self.head = headMatrix
    self.xgrid = xLocations
    self.ygrid = yLocations
    self.dx = dx
    self.dy = dy
    self.nx = nx
    self.ny = ny
    self.k = conductivity
    self.phi = porosity
    self.t = [0.]
    self.vp = [0. 0.]
    self.outofdomain = False


  def ParticleTrack(tMax):
    """ 
    Runs the particle tracking for the parameters that were initialized. 
    Follows the pollock algorithm prett closely this time.
    """

    # Assign particle to cell 
    # assigns the index of the cell being used in the arrays. 
    i, j = AssignInitialCell()

    while outOfDomain == False and self.t < tMax:

      cel = GridCell(self.xgrid[i,j], self.ygrid[i,j], self.dx, self.dy\
          self.k, self.phi)

      
      # compute cell face velocities
      cel.CalcFaceVelocities(self.head[i,j], self.head[i-1,j]\
          self.head[i+1,j], self.head[i,j-1], self.head[i,j+1])

      # determine exit faces
      cel.AssignExitFaces(vp)
      #need some way to make a discharge point possible.


      # Compute transit time, determine actual exit face
      travel = CalcTravelTime(cel)


      # compute exit point coordinates
      #xpout is 2d list
      xpout = CalcPosition(cel, min(travel))

      
      # Determine new cell
      i_add, j_add= AssignNewCell(cel, travel, xpout)
      i = i + i_add
      j = j + j_add


      # write particle coordinates and times to list
      # some domain check goes here: 
      self.xp.append(xp_out)
      self.yp.append(yp_out)
      

class GridCell(xCenter, yCenter, dx, dy, conductivity, porosity):
  """ Grid cell contains the edge coordinates, face velocities and 
  if the particle is a capture point. 
  Takes in the center and cell widths in the following form:
  xCenter, yCenter, dx, dy
  and stores the edge faces as a 2d list.
  x0[0] = xwest
  x0[1] = ysouth
  x1[0] = xeast
  x1[1] = ynorth
  Also stores the porosity and conductivity for transmissivity equations.
  """
  self.x0 = [xCenter - dx/2, yCenter - dy/2]
  self.x1 = [xCenter + dx/2, yCenter + dy/2]
  self.exit0 = [True, True]
  self.exit1 = [True, True]
  self.dx = [dx, dy]
  self.v0 = [0., 0.]
  self.v1 = [0., 0.]
  self.A = [0., 0.]
  self.k = conductivity
  self.phi = porosity

  def CalcFaceVelocities(hC, hW, hE, hS, hN):
    """ 
    Calculates the velocities on each face by linearly interpolating between
    the heads of the 5-point finite difference
    """
    
    self.v0[0] = k * (hW - hC) / (self.dx[0] * self.phi)
    self.v0[1] = k * (hC - hE) / (self.dx[1] * self.phi)
    self.v1[0] = k * (hS - hC) / (self.dx[0] * self.phi)
    self.v1[1] = k * (hC - hN) / (self.dx[1] * self.phi)
    
    for i in range(2):
      A[i] = (v1[i] - v0[i])/ dx[i]
    
  def AssignExitFaces(vp):
    """
    removes the possiblity of exiting through faces due to certain velocity conditions. 
    vp must be a two-dimensional list where
    vp[0] = vx
    vp[1] = vy
    """
    v0 = self.v0
    v1 = self.v1
    tol = 1.e-6
    for i in range(2):

      # positive velocities
      if (v0[i] > tol and v1[i] > tol):
        exit0[i] = False
      
      #negative velocities
      elif (v0[i] < -tol and v1[i] < -tol):
        exit1[i] = False
      
      #sink condition
      elif (v0[i] > tol and v1[i] < -tol):
        exit0[i] = False
        exit1[i] = False
      
      # flow divide condition
      elif (v0[i] < - tol and v1[i] > tol ):
        if vp[i] > 0. :
          exit0[i] = False
        else:
          exit1[i] = False

      # no flow on left edge
      elif (v0[i] * v0[i] < tol * tol ):
        exit0[i] = False
        if v1[i] < 0. : 
          exit1[i] = False 

      #no flow on right edge
      elif (v1[i] * v1[i] < tol * tol ):
        exit1[i] = False
        if v0[i] > 0. :
          exit0[i] = False

      else:
        print "invalid velocity condition or input error"
        print "vx0, vx1, vp"
        for i in range(2):
          print vx0[i], vx1[i], vp[i]
        exit(1)

